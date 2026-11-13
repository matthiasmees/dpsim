/* SPDX-License-Identifier: MPL-2.0 */
#include <dpsim/HipLinearSSNFusedKernels.h>

#include <hip/hip_runtime.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace DPsim {
namespace {
struct DeviceLinearSSNModel {
  int kind = 0;
  int grounded[2]{};
  int nodes[6]{};
  int referenceTerminal = 0;
  int hasInductance = 0;
  int hasCapacitance = 0;
  double ratio = 1.0;
  double inductorAdmittanceScalar = 0.0;
  double capacitorAdmittanceScalar = 0.0;
  double matrices[8][9]{};
  double state[18]{};
  double history[15]{};
};

struct DeviceLinearSSNBatch {
  DeviceLinearSSNModel *models = nullptr;
  double *histories = nullptr;
  double *solution = nullptr;
  std::size_t modelCount = 0;
  std::size_t solutionCapacity = 0;
};

constexpr unsigned int blockSize = 128;

void checkHip(hipError_t status, const char *operation) {
  if (status != hipSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             hipGetErrorString(status));
}

void copyMatrix(const CPS::Matrix &source, double *destination) {
  std::fill(destination, destination + 9, 0.0);
  if (source.rows() > 3 || source.cols() > 3)
    throw std::invalid_argument("HIP fused SSN matrix exceeds 3x3");
  for (Eigen::Index row = 0; row < source.rows(); ++row)
    for (Eigen::Index col = 0; col < source.cols(); ++col)
      destination[col * 3 + row] = source(row, col);
}

void packModel(const CPS::EMT::LinearSSNModel &source,
               DeviceLinearSSNModel &destination) {
  using Kind = CPS::EMT::LinearSSNModel::Kind;
  destination.kind = static_cast<int>(source.kind);
  for (int index = 0; index < 2; ++index)
    destination.grounded[index] = source.grounded[index] ? 1 : 0;
  for (int index = 0; index < 6; ++index)
    destination.nodes[index] = static_cast<int>(source.nodeIndices[index]);
  destination.referenceTerminal = static_cast<int>(source.referenceTerminal);
  destination.ratio = source.ratioFromReference;
  destination.hasInductance = source.hasLoadInductance ? 1 : 0;
  destination.hasCapacitance = source.hasLoadCapacitance ? 1 : 0;
  destination.inductorAdmittanceScalar = source.inductorAdmittanceScalar;
  destination.capacitorAdmittanceScalar = source.capacitorAdmittanceScalar;
  if (source.state.size() > 18)
    throw std::invalid_argument("HIP fused SSN state exceeds 18 values");
  for (Eigen::Index index = 0; index < source.state.size(); ++index)
    destination.state[index] = source.state(index);
  switch (source.kind) {
  case Kind::PiLine:
    copyMatrix(source.seriesResistance, destination.matrices[0]);
    copyMatrix(source.inductorAdmittance, destination.matrices[1]);
    copyMatrix(source.seriesSystemInverse, destination.matrices[2]);
    copyMatrix(source.seriesAdmittance, destination.matrices[3]);
    copyMatrix(source.capAdmittance, destination.matrices[4]);
    break;
  case Kind::ImpedanceLoad:
    break;
  case Kind::Transformer:
    copyMatrix(source.seriesResistanceHalf, destination.matrices[0]);
    copyMatrix(source.inductorAdmittance, destination.matrices[1]);
    copyMatrix(source.seriesSystemInverse, destination.matrices[2]);
    copyMatrix(source.seriesAdmittance, destination.matrices[3]);
    copyMatrix(source.magnetizingInductorAdmittance,
               destination.matrices[4]);
    copyMatrix(source.midpointInverse, destination.matrices[5]);
    break;
  }
}

__device__ void matrixVector(const double *matrix, const double *vector,
                             double *result) {
  for (int row = 0; row < 3; ++row) {
    result[row] = 0.0;
    for (int col = 0; col < 3; ++col)
      result[row] += matrix[col * 3 + row] * vector[col];
  }
}

__device__ void matrixSumProduct(const double *matrix, const double *first,
                                 const double *second, double *result) {
  double sum[3];
  matrixVector(matrix, first, sum);
  for (int i = 0; i < 3; ++i)
    sum[i] += second[i];
  for (int i = 0; i < 3; ++i)
    result[i] = sum[i];
}

__global__ void updateHistoryKernel(DeviceLinearSSNModel *models,
                                    std::size_t count) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= count)
    return;
  auto &model = models[index];
  if (model.kind == 0) {
    double inductorHistory[3];
    matrixSumProduct(model.matrices[1], &model.state[3], &model.state[0],
                     inductorHistory);
    matrixVector(model.matrices[2], inductorHistory, &model.history[0]);
    double cap0[3], cap1[3];
    matrixVector(model.matrices[4], &model.state[6], cap0);
    matrixVector(model.matrices[4], &model.state[12], cap1);
    for (int phase = 0; phase < 3; ++phase) {
      model.history[3 + phase] = -model.state[9 + phase] - cap0[phase];
      model.history[6 + phase] = -model.state[15 + phase] - cap1[phase];
    }
  } else if (model.kind == 1) {
    for (int phase = 0; phase < 3; ++phase) {
      model.history[phase] = model.hasInductance
                                 ? model.inductorAdmittanceScalar *
                                           model.state[phase] +
                                       model.state[3 + phase]
                                 : 0.0;
      model.history[3 + phase] = model.hasCapacitance
                                     ? -model.state[9 + phase] -
                                           model.capacitorAdmittanceScalar *
                                               model.state[6 + phase]
                                     : 0.0;
    }
  } else {
    double h0inner[3], h1inner[3], h0[3], h1[3], hmag[3];
    matrixSumProduct(model.matrices[1], &model.state[0], &model.state[3],
                     h0inner);
    matrixSumProduct(model.matrices[1], &model.state[6], &model.state[9],
                     h1inner);
    matrixVector(model.matrices[2], h0inner, h0);
    matrixVector(model.matrices[2], h1inner, h1);
    matrixSumProduct(model.matrices[4], &model.state[12], &model.state[15],
                     hmag);
    double midpointHistory[3], midpointSource[3], seriesMidpoint[3];
    for (int phase = 0; phase < 3; ++phase)
      midpointHistory[phase] = -h0[phase] + h1[phase] + hmag[phase];
    matrixVector(model.matrices[5], midpointHistory, midpointSource);
    matrixVector(model.matrices[3], midpointSource, seriesMidpoint);
    for (int phase = 0; phase < 3; ++phase) {
      model.history[phase] = h0[phase];
      model.history[3 + phase] = h1[phase];
      model.history[6 + phase] = hmag[phase];
      model.history[9 + phase] = h0[phase] + seriesMidpoint[phase];
      model.history[12 + phase] =
          model.ratio * (seriesMidpoint[phase] - h1[phase]);
    }
  }
}

__global__ void gatherHistoryKernel(const DeviceLinearSSNModel *models,
                                    double *histories,
                                    std::size_t modelCount) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= modelCount * 15)
    return;
  const std::size_t model = index / 15;
  const std::size_t item = index % 15;
  histories[index] = models[model].history[item];
}

__global__ void accumulateHistoryRhsKernel(
    const DeviceLinearSSNModel *models, std::size_t modelCount,
    double *rightHandSide) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= modelCount)
    return;
  const auto &model = models[index];
  if (model.kind == 0) {
    for (int phase = 0; phase < 3; ++phase) {
      if (!model.grounded[0])
        atomicAdd(&rightHandSide[model.nodes[phase]],
                  model.history[phase] - model.history[3 + phase]);
      if (!model.grounded[1])
        atomicAdd(&rightHandSide[model.nodes[3 + phase]],
                  -model.history[phase] - model.history[6 + phase]);
    }
  } else if (model.kind == 1) {
    if (!model.grounded[0])
      for (int phase = 0; phase < 3; ++phase)
        atomicAdd(&rightHandSide[model.nodes[phase]],
                  -(model.history[phase] + model.history[3 + phase]));
  } else {
    const int reference = model.referenceTerminal;
    const int other = 1 - reference;
    if (!model.grounded[reference])
      for (int phase = 0; phase < 3; ++phase)
        atomicAdd(&rightHandSide[model.nodes[reference * 3 + phase]],
                  model.history[9 + phase]);
    if (!model.grounded[other])
      for (int phase = 0; phase < 3; ++phase)
        atomicAdd(&rightHandSide[model.nodes[other * 3 + phase]],
                  model.history[12 + phase]);
  }
}

__device__ double nodeValue(const DeviceLinearSSNModel &model,
                            const double *solution, std::size_t size,
                            int localIndex, bool grounded) {
  if (grounded)
    return 0.0;
  const int node = model.nodes[localIndex];
  return node >= 0 && static_cast<std::size_t>(node) < size ? solution[node]
                                                             : 0.0;
}

__global__ void updateStateKernel(DeviceLinearSSNModel *models,
                                  std::size_t count, const double *solution,
                                  std::size_t solutionSize) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= count)
    return;
  auto &model = models[index];
  if (model.kind == 0) {
    double v0[3], v1[3], interfaceVoltage[3], current[3], voltage[3];
    for (int phase = 0; phase < 3; ++phase) {
      v0[phase] = nodeValue(model, solution, solutionSize, phase,
                            model.grounded[0] != 0);
      v1[phase] = nodeValue(model, solution, solutionSize, 3 + phase,
                            model.grounded[1] != 0);
      interfaceVoltage[phase] = v1[phase] - v0[phase];
    }
    double seriesCurrent[3];
    matrixVector(model.matrices[3], interfaceVoltage, seriesCurrent);
    for (int phase = 0; phase < 3; ++phase)
      current[phase] = seriesCurrent[phase] + model.history[phase];
    matrixVector(model.matrices[0], current, voltage);
    for (int phase = 0; phase < 3; ++phase) {
      model.state[phase] = current[phase];
      model.state[3 + phase] = interfaceVoltage[phase] - voltage[phase];
    }
    double cap0[3], cap1[3];
    matrixVector(model.matrices[4], v0, cap0);
    matrixVector(model.matrices[4], v1, cap1);
    for (int phase = 0; phase < 3; ++phase) {
      model.state[6 + phase] = v0[phase];
      model.state[9 + phase] = cap0[phase] + model.history[3 + phase];
      model.state[12 + phase] = v1[phase];
      model.state[15 + phase] = cap1[phase] + model.history[6 + phase];
    }
  } else if (model.kind == 1) {
    for (int phase = 0; phase < 3; ++phase) {
      const double voltage = nodeValue(model, solution, solutionSize, phase,
                                        model.grounded[0] != 0);
      if (model.hasInductance) {
        model.state[phase] = voltage;
        model.state[3 + phase] = model.inductorAdmittanceScalar * voltage +
                                 model.history[phase];
      }
      if (model.hasCapacitance) {
        model.state[6 + phase] = voltage;
        model.state[9 + phase] = model.capacitorAdmittanceScalar * voltage +
                                 model.history[3 + phase];
      }
    }
  } else {
    double v0[3], v1[3], referredOther[3], reference[3];
    for (int phase = 0; phase < 3; ++phase) {
      v0[phase] = nodeValue(model, solution, solutionSize, phase,
                            model.grounded[0] != 0);
      v1[phase] = nodeValue(model, solution, solutionSize, 3 + phase,
                            model.grounded[1] != 0);
      reference[phase] = model.referenceTerminal == 0 ? v0[phase] : v1[phase];
      referredOther[phase] =
          model.ratio * (model.referenceTerminal == 0 ? v1[phase] : v0[phase]);
    }
    double midpointHistory[3], sumVoltage[3], midpointVoltage[3];
    double sRef[3], sOther[3];
    matrixVector(model.matrices[3], reference, sRef);
    matrixVector(model.matrices[3], referredOther, sOther);
    for (int phase = 0; phase < 3; ++phase) {
      midpointHistory[phase] = -model.history[phase] +
                               model.history[3 + phase] +
                               model.history[6 + phase];
      sumVoltage[phase] = sRef[phase] + sOther[phase] + midpointHistory[phase];
    }
    matrixVector(model.matrices[5], sumVoltage, midpointVoltage);
    double total0[3], total1[3], current0Base[3], current1Base[3];
    for (int phase = 0; phase < 3; ++phase) {
      total0[phase] = midpointVoltage[phase] - reference[phase];
      total1[phase] = referredOther[phase] - midpointVoltage[phase];
    }
    matrixVector(model.matrices[3], total0, current0Base);
    matrixVector(model.matrices[3], total1, current1Base);
    double drop0[3], drop1[3], magCurrent[3];
    for (int phase = 0; phase < 3; ++phase) {
      current0Base[phase] += model.history[phase];
      current1Base[phase] += model.history[3 + phase];
    }
    matrixVector(model.matrices[0], current0Base, drop0);
    matrixVector(model.matrices[0], current1Base, drop1);
    double magVoltage[3];
    for (int phase = 0; phase < 3; ++phase)
      magVoltage[phase] = -midpointVoltage[phase];
    matrixVector(model.matrices[4], magVoltage, magCurrent);
    for (int phase = 0; phase < 3; ++phase) {
      model.state[0 + phase] = total0[phase] - drop0[phase];
      model.state[3 + phase] = current0Base[phase];
      model.state[6 + phase] = total1[phase] - drop1[phase];
      model.state[9 + phase] = current1Base[phase];
      model.state[12 + phase] = magVoltage[phase];
      model.state[15 + phase] = magCurrent[phase] + model.history[6 + phase];
    }
  }
}
} // namespace

void *hipLinearSSNCreate(const CPS::EMT::LinearSSNModel *models,
                         std::size_t modelCount) {
  if (!models || modelCount == 0)
    throw std::invalid_argument("HIP fused SSN model list is empty");
  std::vector<DeviceLinearSSNModel> packed(modelCount);
  for (std::size_t index = 0; index < modelCount; ++index)
    packModel(models[index], packed[index]);
  auto *batch = new DeviceLinearSSNBatch;
  batch->modelCount = modelCount;
  try {
    checkHip(hipMalloc(&batch->models,
                       modelCount * sizeof(DeviceLinearSSNModel)),
             "hipMalloc(fused SSN models)");
    checkHip(hipMalloc(&batch->histories,
                       modelCount * 15 * sizeof(double)),
             "hipMalloc(fused SSN history buffer)");
    checkHip(hipMemcpy(batch->models, packed.data(),
                       modelCount * sizeof(DeviceLinearSSNModel),
                       hipMemcpyHostToDevice),
             "hipMemcpy(fused SSN models)");
  } catch (...) {
    if (batch->models)
      (void)hipFree(batch->models);
    if (batch->histories)
      (void)hipFree(batch->histories);
    delete batch;
    throw;
  }
  return batch;
}

void hipLinearSSNUpdateHistories(void *deviceModels, std::size_t modelCount,
                                 double *hostHistoryRows) {
  if (!deviceModels || !hostHistoryRows)
    throw std::invalid_argument("HIP fused SSN history buffers are required");
  auto *batch = static_cast<DeviceLinearSSNBatch *>(deviceModels);
  if (modelCount != batch->modelCount)
    throw std::invalid_argument("HIP fused SSN model count changed");
  const auto blocks = static_cast<unsigned int>(
      (modelCount + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(updateHistoryKernel, dim3(blocks), dim3(blockSize), 0,
                     nullptr, batch->models, modelCount);
  checkHip(hipGetLastError(), "HIP fused SSN history kernel launch");
  const auto historyBlocks = static_cast<unsigned int>(
      (modelCount * 15 + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(gatherHistoryKernel, dim3(historyBlocks), dim3(blockSize),
                     0, nullptr, batch->models, batch->histories, modelCount);
  checkHip(hipGetLastError(), "HIP fused SSN history gather launch");
  checkHip(hipMemcpy(hostHistoryRows, batch->histories,
                     modelCount * 15 * sizeof(double), hipMemcpyDeviceToHost),
           "hipMemcpy(fused SSN histories)");
}

void hipLinearSSNUpdateHistoriesDevice(void *deviceModels,
                                       std::size_t modelCount) {
  if (!deviceModels)
    throw std::invalid_argument("HIP fused SSN device batch is required");
  auto *batch = static_cast<DeviceLinearSSNBatch *>(deviceModels);
  if (modelCount != batch->modelCount)
    throw std::invalid_argument("HIP fused SSN model count changed");
  const auto blocks = static_cast<unsigned int>(
      (modelCount + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(updateHistoryKernel, dim3(blocks), dim3(blockSize), 0,
                     nullptr, batch->models, modelCount);
  checkHip(hipGetLastError(), "HIP fused SSN history kernel launch");
}

void hipLinearSSNAccumulateRhs(void *deviceModels,
                               double *deviceRightHandSide) {
  if (!deviceModels || !deviceRightHandSide)
    throw std::invalid_argument("HIP fused SSN RHS buffers are required");
  auto *batch = static_cast<DeviceLinearSSNBatch *>(deviceModels);
  const auto blocks = static_cast<unsigned int>(
      (batch->modelCount + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(accumulateHistoryRhsKernel, dim3(blocks), dim3(blockSize),
                     0, nullptr, batch->models, batch->modelCount,
                     deviceRightHandSide);
  checkHip(hipGetLastError(), "HIP fused SSN RHS accumulation launch");
}

void hipLinearSSNUpdateStates(void *deviceModels, std::size_t modelCount,
                              const double *hostSolution,
                              std::size_t solutionSize) {
  if (!deviceModels || !hostSolution || solutionSize == 0)
    throw std::invalid_argument("HIP fused SSN state update inputs are empty");
  auto *batch = static_cast<DeviceLinearSSNBatch *>(deviceModels);
  if (modelCount != batch->modelCount)
    throw std::invalid_argument("HIP fused SSN model count changed");
  if (solutionSize > batch->solutionCapacity) {
    if (batch->solution)
      checkHip(hipFree(batch->solution), "hipFree(old fused SSN solution)");
    batch->solution = nullptr;
    checkHip(hipMalloc(&batch->solution, solutionSize * sizeof(double)),
             "hipMalloc(fused SSN solution)");
    batch->solutionCapacity = solutionSize;
  }
  try {
    checkHip(hipMemcpy(batch->solution, hostSolution,
                       solutionSize * sizeof(double), hipMemcpyHostToDevice),
             "hipMemcpy(fused SSN solution)");
    const auto blocks = static_cast<unsigned int>(
        (modelCount + blockSize - 1) / blockSize);
    hipLaunchKernelGGL(updateStateKernel, dim3(blocks), dim3(blockSize), 0,
                       nullptr,
                       batch->models, modelCount, batch->solution,
                       solutionSize);
    checkHip(hipGetLastError(), "HIP fused SSN state kernel launch");
    checkHip(hipDeviceSynchronize(), "HIP fused SSN state update");
  } catch (...) {
    throw;
  }
}

void hipLinearSSNUpdateStatesDevice(void *deviceModels,
                                    std::size_t modelCount,
                                    const double *deviceSolution,
                                    std::size_t solutionSize) {
  if (!deviceModels || !deviceSolution || solutionSize == 0)
    throw std::invalid_argument("HIP fused SSN device solution is required");
  auto *batch = static_cast<DeviceLinearSSNBatch *>(deviceModels);
  if (modelCount != batch->modelCount)
    throw std::invalid_argument("HIP fused SSN model count changed");
  const auto blocks = static_cast<unsigned int>(
      (modelCount + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(updateStateKernel, dim3(blocks), dim3(blockSize), 0,
                     nullptr, batch->models, modelCount, deviceSolution,
                     solutionSize);
  checkHip(hipGetLastError(), "HIP fused SSN device state kernel launch");
}

void hipLinearSSNDestroy(void *deviceModels) noexcept {
  if (!deviceModels)
    return;
  auto *batch = static_cast<DeviceLinearSSNBatch *>(deviceModels);
  if (batch->models)
    (void)hipFree(batch->models);
  if (batch->histories)
    (void)hipFree(batch->histories);
  if (batch->solution)
    (void)hipFree(batch->solution);
  delete batch;
}
} // namespace DPsim
