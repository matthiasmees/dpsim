/* SPDX-License-Identifier: MPL-2.0 */
#include <dpsim/HipSG4BatchKernels.h>

#include <hip/hip_runtime.h>

#include <stdexcept>
#include <string>

namespace DPsim {
namespace {
using State = CPS::EMT::Ph3::HipSG4State;
constexpr unsigned int blockSize = 128;

struct DeviceBatch {
  State *states = nullptr;
  double *solution = nullptr;
  int *matrixPositions = nullptr;
  int *restorePositions = nullptr;
  double *restoreValues = nullptr;
  double *matrixValues = nullptr;
  double *rightHandSide = nullptr;
  std::size_t count = 0;
  std::size_t solutionCapacity = 0;
  std::size_t restoreCount = 0;
};

void checkHip(hipError_t status, const char *operation) {
  if (status != hipSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             hipGetErrorString(status));
}

__device__ void multiply3(const double *a, const double *b, double *c) {
  for (int col = 0; col < 3; ++col)
    for (int row = 0; row < 3; ++row) {
      c[col * 3 + row] = 0.0;
      for (int k = 0; k < 3; ++k)
        c[col * 3 + row] += a[k * 3 + row] * b[col * 3 + k];
    }
}

__device__ void multiplyVector3(const double *a, const double *x,
                                double *y) {
  for (int row = 0; row < 3; ++row) {
    y[row] = 0.0;
    for (int col = 0; col < 3; ++col)
      y[row] += a[col * 3 + row] * x[col];
  }
}

__device__ void makeParkMatrices(double theta, double *abcToDq0,
                                 double *dq0ToAbc) {
  const double thetaB = theta - 2.0 * 3.141592653589793238462643383279502884 / 3.0;
  const double thetaC = theta + 2.0 * 3.141592653589793238462643383279502884 / 3.0;
  const double angles[3] = {theta, thetaB, thetaC};
  for (int phase = 0; phase < 3; ++phase) {
    const double c = cos(angles[phase]);
    const double s = sin(angles[phase]);
    abcToDq0[phase * 3 + 0] = (2.0 / 3.0) * c;
    abcToDq0[phase * 3 + 1] = -(2.0 / 3.0) * s;
    abcToDq0[phase * 3 + 2] = 1.0 / 3.0;
    dq0ToAbc[0 * 3 + phase] = c;
    dq0ToAbc[1 * 3 + phase] = -s;
    dq0ToAbc[2 * 3 + phase] = 1.0;
  }
}

__global__ void preStepKernel(State *states, std::size_t count, double time) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= count)
    return;
  State &s = states[index];
  s.electricalTorque =
      s.vdq0[0] * s.idq0[0] + s.vdq0[1] * s.idq0[1];
  s.omega = s.omega + s.timestep *
                           (1.0 / (2.0 * s.inertia) *
                            (s.previousMechanicalTorque - s.electricalTorque -
                             s.damping * (s.omega - 1.0)));
  s.theta = s.theta + s.timestep * (s.omega * s.baseOmegaMech);
  s.delta = s.delta + s.timestep * (s.omega - 1.0) * s.baseOmegaMech;

  if (time > 0.0) {
    s.edq0Transient[0] = -s.idq0[1] * s.lqTransient + s.vdq0[0];
    s.edq0Transient[1] = s.idq0[0] * s.ldTransient + s.vdq0[1];
    s.edq0Transient[2] = 0.0;
  }

  double abcToDq0[9], dq0ToAbc[9], intermediate[9];
  makeParkMatrices(s.theta, abcToDq0, dq0ToAbc);
  multiply3(dq0ToAbc, s.conductanceDq0, intermediate);
  multiply3(intermediate, abcToDq0, s.conductanceAbc);

  s.historyDq0[0] =
      s.adTransient * s.idq0[1] + s.bdTransient * s.edq0Transient[0];
  s.historyDq0[1] = s.aqTransient * s.idq0[0] +
                    s.bqTransient * s.edq0Transient[1] +
                    s.dqTransient * s.previousFieldVoltage +
                    s.dqTransient * s.fieldVoltage;
  s.historyDq0[2] = 0.0;
  multiplyVector3(dq0ToAbc, s.historyDq0, s.emfAbc);
  for (int phase = 0; phase < 3; ++phase)
    s.emfAbc[phase] *= s.baseVoltage;
  multiplyVector3(s.conductanceAbc, s.emfAbc, s.currentSourceAbc);
}

__global__ void scatterMatrixKernel(const State *states,
                                    const int *matrixPositions,
                                    double *matrixValues,
                                    std::size_t count) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= count)
    return;
  const State &s = states[index];
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column)
      atomicAdd(&matrixValues[matrixPositions[index * 9 + row * 3 + column]],
                s.conductanceAbc[column * 3 + row]);
}

__global__ void restoreMatrixKernel(const int *positions,
                                    const double *baseline,
                                    double *matrixValues,
                                    std::size_t count) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index < count)
    matrixValues[positions[index]] = baseline[index];
}

__global__ void accumulateSg4RhsKernel(const State *states,
                                       double *rightHandSide,
                                       std::size_t count) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= count)
    return;
  const State &s = states[index];
  for (int phase = 0; phase < 3; ++phase)
    atomicAdd(&rightHandSide[s.busRows[phase]], s.currentSourceAbc[phase]);
}

__global__ void postStepKernel(State *states, std::size_t count,
                               const double *solution,
                               std::size_t solutionSize) {
  const std::size_t index = static_cast<std::size_t>(blockIdx.x) * blockDim.x +
                            threadIdx.x;
  if (index >= count)
    return;
  State &s = states[index];
  double abcToDq0[9], dq0ToAbc[9];
  makeParkMatrices(s.theta, abcToDq0, dq0ToAbc);
  for (int phase = 0; phase < 3; ++phase) {
    const int row = s.busRows[phase];
    s.terminalVoltageAbc[phase] =
        !s.grounded && row >= 0 && static_cast<std::size_t>(row) < solutionSize
            ? solution[row]
            : 0.0;
  }
  multiplyVector3(abcToDq0, s.terminalVoltageAbc, s.voltageDq0);
  for (int axis = 0; axis < 3; ++axis)
    s.voltageDq0[axis] /= s.baseVoltage;

  double conductanceCurrent[3];
  multiplyVector3(s.conductanceAbc, s.terminalVoltageAbc, conductanceCurrent);
  for (int phase = 0; phase < 3; ++phase)
    s.terminalCurrentAbc[phase] =
        s.currentSourceAbc[phase] - conductanceCurrent[phase];
  multiplyVector3(abcToDq0, s.terminalCurrentAbc, s.currentDq0);
  for (int axis = 0; axis < 3; ++axis) {
    s.currentDq0[axis] /= s.baseCurrent;
    s.vdq0[axis] = s.voltageDq0[axis];
    s.idq0[axis] = s.currentDq0[axis];
  }
}
} // namespace

void *hipSG4BatchCreate(const State *hostStates, std::size_t count) {
  if (!hostStates || count == 0)
    throw std::invalid_argument("HIP SG4 batch requires at least one state");
  auto *batch = new DeviceBatch;
  batch->count = count;
  try {
    checkHip(hipMalloc(&batch->states, count * sizeof(State)),
             "hipMalloc(SG4 states)");
    checkHip(hipMemcpy(batch->states, hostStates, count * sizeof(State),
                       hipMemcpyHostToDevice),
             "hipMemcpy(SG4 initial state)");
  } catch (...) {
    if (batch->states)
      (void)hipFree(batch->states);
    delete batch;
    throw;
  }
  return batch;
}

void hipSG4BatchPreStep(void *opaqueBatch, double time, State *hostStates) {
  auto *batch = static_cast<DeviceBatch *>(opaqueBatch);
  if (!batch || (batch->matrixValues && !batch->matrixPositions) ||
      (!hostStates && !batch->matrixValues))
    throw std::invalid_argument("HIP SG4 pre-step buffers are required");
  const auto blocks =
      static_cast<unsigned int>((batch->count + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(preStepKernel, dim3(blocks), dim3(blockSize), 0, nullptr,
                     batch->states, batch->count, time);
  checkHip(hipGetLastError(), "HIP SG4 pre-step kernel launch");
  if (hostStates)
    checkHip(hipMemcpy(hostStates, batch->states, batch->count * sizeof(State),
                       hipMemcpyDeviceToHost),
             "hipMemcpy(SG4 pre-step outputs)");
}

void hipSG4BatchScatterMatrix(void *opaqueBatch) {
  auto *batch = static_cast<DeviceBatch *>(opaqueBatch);
  if (!batch || !batch->matrixPositions || !batch->matrixValues)
    throw std::invalid_argument("HIP SG4 matrix buffers are not configured");
  if (batch->restoreCount != 0) {
    const auto restoreBlocks = static_cast<unsigned int>(
        (batch->restoreCount + blockSize - 1) / blockSize);
    hipLaunchKernelGGL(restoreMatrixKernel, dim3(restoreBlocks),
                       dim3(blockSize), 0, nullptr, batch->restorePositions,
                       batch->restoreValues, batch->matrixValues,
                       batch->restoreCount);
    checkHip(hipGetLastError(), "HIP SG4 matrix restore launch");
  }
  const auto blocks = static_cast<unsigned int>(
      (batch->count + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(scatterMatrixKernel, dim3(blocks), dim3(blockSize), 0,
                     nullptr, batch->states, batch->matrixPositions,
                     batch->matrixValues, batch->count);
  checkHip(hipGetLastError(), "HIP SG4 matrix scatter launch");
}

void hipSG4BatchConfigureAssembly(void *opaqueBatch,
                                 const int *hostMatrixPositions,
                                 double *deviceMatrixValues,
                                 double *deviceRightHandSide,
                                 const int *hostRestorePositions,
                                 const double *hostRestoreValues,
                                 std::size_t restoreCount) {
  auto *batch = static_cast<DeviceBatch *>(opaqueBatch);
  if (!batch || !hostMatrixPositions || !deviceMatrixValues ||
      !deviceRightHandSide ||
      (restoreCount != 0 && (!hostRestorePositions || !hostRestoreValues)))
    throw std::invalid_argument("HIP SG4 assembly buffers are required");
  if (batch->matrixValues)
    throw std::logic_error("HIP SG4 assembly is already configured");
  checkHip(hipMalloc(&batch->matrixPositions,
                     batch->count * 9 * sizeof(int)),
           "hipMalloc(SG4 matrix positions)");
  try {
    checkHip(hipMemcpy(batch->matrixPositions, hostMatrixPositions,
                       batch->count * 9 * sizeof(int),
                       hipMemcpyHostToDevice),
             "hipMemcpy(SG4 matrix positions)");
  } catch (...) {
    (void)hipFree(batch->matrixPositions);
    batch->matrixPositions = nullptr;
    throw;
  }
  if (restoreCount != 0) {
    checkHip(hipMalloc(&batch->restorePositions,
                       restoreCount * sizeof(int)),
             "hipMalloc(SG4 restore positions)");
    try {
      checkHip(hipMalloc(&batch->restoreValues,
                         restoreCount * sizeof(double)),
               "hipMalloc(SG4 restore values)");
      checkHip(hipMemcpy(batch->restorePositions, hostRestorePositions,
                         restoreCount * sizeof(int), hipMemcpyHostToDevice),
               "hipMemcpy(SG4 restore positions)");
      checkHip(hipMemcpy(batch->restoreValues, hostRestoreValues,
                         restoreCount * sizeof(double), hipMemcpyHostToDevice),
               "hipMemcpy(SG4 restore values)");
    } catch (...) {
      if (batch->restorePositions) {
        (void)hipFree(batch->restorePositions);
        batch->restorePositions = nullptr;
      }
      if (batch->restoreValues) {
        (void)hipFree(batch->restoreValues);
        batch->restoreValues = nullptr;
      }
      throw;
    }
  }
  batch->restoreCount = restoreCount;
  batch->matrixValues = deviceMatrixValues;
  batch->rightHandSide = deviceRightHandSide;
}

void hipSG4BatchAccumulateRhs(void *opaqueBatch) {
  auto *batch = static_cast<DeviceBatch *>(opaqueBatch);
  if (!batch || !batch->rightHandSide)
    throw std::invalid_argument("HIP SG4 RHS buffer is not configured");
  const auto blocks = static_cast<unsigned int>(
      (batch->count + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(accumulateSg4RhsKernel, dim3(blocks), dim3(blockSize), 0,
                     nullptr, batch->states, batch->rightHandSide,
                     batch->count);
  checkHip(hipGetLastError(), "HIP SG4 RHS accumulation launch");
}

void hipSG4BatchPostStep(void *opaqueBatch, const double *hostSolution,
                         std::size_t solutionSize, State *hostStates) {
  auto *batch = static_cast<DeviceBatch *>(opaqueBatch);
  if (!batch || !hostSolution || !hostStates || solutionSize == 0)
    throw std::invalid_argument("HIP SG4 post-step inputs are required");
  if (solutionSize > batch->solutionCapacity) {
    if (batch->solution)
      checkHip(hipFree(batch->solution), "hipFree(old SG4 solution)");
    batch->solution = nullptr;
    checkHip(hipMalloc(&batch->solution, solutionSize * sizeof(double)),
             "hipMalloc(SG4 solution)");
    batch->solutionCapacity = solutionSize;
  }
  checkHip(hipMemcpy(batch->solution, hostSolution,
                     solutionSize * sizeof(double), hipMemcpyHostToDevice),
           "hipMemcpy(SG4 solution)");
  const auto blocks =
      static_cast<unsigned int>((batch->count + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(postStepKernel, dim3(blocks), dim3(blockSize), 0, nullptr,
                     batch->states, batch->count, batch->solution,
                     solutionSize);
  checkHip(hipGetLastError(), "HIP SG4 post-step kernel launch");
  checkHip(hipMemcpy(hostStates, batch->states, batch->count * sizeof(State),
                     hipMemcpyDeviceToHost),
           "hipMemcpy(SG4 post-step outputs)");
}

void hipSG4BatchPostStepDevice(void *opaqueBatch, const double *deviceSolution,
                               std::size_t solutionSize, State *hostStates) {
  auto *batch = static_cast<DeviceBatch *>(opaqueBatch);
  if (!batch || !deviceSolution || solutionSize == 0)
    throw std::invalid_argument("HIP SG4 post-step inputs are required");
  const auto blocks = static_cast<unsigned int>(
      (batch->count + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(postStepKernel, dim3(blocks), dim3(blockSize), 0, nullptr,
                     batch->states, batch->count, deviceSolution,
                     solutionSize);
  checkHip(hipGetLastError(), "HIP SG4 post-step kernel launch");
  if (hostStates)
    checkHip(hipMemcpy(hostStates, batch->states,
                       batch->count * sizeof(State), hipMemcpyDeviceToHost),
             "hipMemcpy(SG4 post-step outputs)");
}

void hipSG4BatchDestroy(void *opaqueBatch) noexcept {
  auto *batch = static_cast<DeviceBatch *>(opaqueBatch);
  if (!batch)
    return;
  if (batch->states)
    (void)hipFree(batch->states);
  if (batch->solution)
    (void)hipFree(batch->solution);
  if (batch->matrixPositions)
    (void)hipFree(batch->matrixPositions);
  if (batch->restorePositions)
    (void)hipFree(batch->restorePositions);
  if (batch->restoreValues)
    (void)hipFree(batch->restoreValues);
  delete batch;
}
} // namespace DPsim
