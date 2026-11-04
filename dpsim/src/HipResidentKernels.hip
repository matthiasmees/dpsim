/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/HipResidentKernels.h>

#include <hip/hip_runtime.h>

#include <stdexcept>
#include <string>

namespace DPsim {
namespace {
constexpr unsigned int blockSize = 256;

__global__ void accumulateIndexedRhsKernel(double *rightHandSide,
                                           const double *contributions,
                                           const int *rows,
                                           std::size_t count) {
  const std::size_t index =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (index < count)
    atomicAdd(&rightHandSide[rows[index]], contributions[index]);
}

__global__ void multiplyVectorKernel(double *values, const double *scale,
                                     std::size_t count) {
  const std::size_t index =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (index < count)
    values[index] *= scale[index];
}

__global__ void scatterSparseValuesKernel(double *values,
                                          const int *positions,
                                          const double *updates,
                                          std::size_t count) {
  const std::size_t index =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (index < count)
    values[positions[index]] = updates[index];
}

__device__ double blockDotProduct(const double *left, const double *right,
                                  int count, double *partial) {
  const unsigned int lane = threadIdx.x;
  double sum = 0.0;
  for (int index = static_cast<int>(lane); index < count;
       index += static_cast<int>(blockDim.x))
    sum += left[index] * right[index];
  partial[lane] = sum;
  __syncthreads();
  for (unsigned int stride = blockDim.x / 2; stride != 0; stride >>= 1) {
    if (lane < stride)
      partial[lane] += partial[lane + stride];
    __syncthreads();
  }
  const double result = partial[0];
  __syncthreads();
  return result;
}

__device__ void residentSpmv(const int *rowPointers,
                             const int *columnIndices, const double *values,
                             const double *rowScale, const double *x,
                             double *y, int count) {
  const unsigned int lane = threadIdx.x;
  for (int row = static_cast<int>(lane); row < count;
       row += static_cast<int>(blockDim.x)) {
    double sum = 0.0;
    for (int entry = rowPointers[row]; entry < rowPointers[row + 1]; ++entry)
      sum += values[entry] * x[columnIndices[entry]];
    y[row] = sum * rowScale[row];
  }
  __syncthreads();
}

__global__ void residentBiCGStabKernel(
    const int *rowPointers, const int *columnIndices, const double *values,
    const double *rowScale, const double *rightHandSide, double *solution,
    double *residual, double *shadowResidual, double *direction,
    double *product, double *intermediate, double *intermediateProduct,
    int count, int maximumIterations, double relativeTolerance,
    int *status) {
  __shared__ double partial[blockSize];
  __shared__ double alphaShared;
  __shared__ double omegaShared;
  __shared__ double rhoPreviousShared;
  __shared__ double rhoShared;
  __shared__ double betaShared;
  __shared__ double targetResidualShared;
  __shared__ int doneShared;
  __shared__ int resultShared;
  const unsigned int lane = threadIdx.x;

  for (int index = static_cast<int>(lane); index < count;
       index += static_cast<int>(blockDim.x)) {
    solution[index] = 0.0;
    residual[index] = rightHandSide[index] * rowScale[index];
    shadowResidual[index] = residual[index];
    direction[index] = 0.0;
    product[index] = 0.0;
    intermediate[index] = 0.0;
    intermediateProduct[index] = 0.0;
  }
  __syncthreads();

  const double initialNormSquared =
      blockDotProduct(residual, residual, count, partial);
  if (lane == 0) {
    const double initialNorm = sqrt(initialNormSquared);
    targetResidualShared = relativeTolerance * fmax(1.0, initialNorm);
    doneShared = initialNorm <= targetResidualShared;
    resultShared = doneShared ? 1 : 0;
    rhoPreviousShared = 1.0;
    alphaShared = 1.0;
    omegaShared = 1.0;
  }
  __syncthreads();

  for (int iteration = 0; iteration < maximumIterations; ++iteration) {
    if (doneShared)
      break;
    const double rho = blockDotProduct(shadowResidual, residual, count, partial);
    if (lane == 0) {
      rhoShared = rho;
      if (!isfinite(rho) || fabs(rho) <= 2.2204460492503131e-16) {
        doneShared = 1;
        resultShared = -iteration - 1;
      } else if (iteration == 0) {
        betaShared = 0.0;
      } else if (!isfinite(omegaShared) ||
                 fabs(omegaShared) <= 2.2204460492503131e-16) {
        doneShared = 1;
        resultShared = -iteration - 1;
      } else {
        betaShared = (rho / rhoPreviousShared) *
                     (alphaShared / omegaShared);
      }
    }
    __syncthreads();
    if (doneShared)
      break;

    if (iteration == 0) {
      for (int index = static_cast<int>(lane); index < count;
           index += static_cast<int>(blockDim.x))
        direction[index] = residual[index];
    } else {
      for (int index = static_cast<int>(lane); index < count;
           index += static_cast<int>(blockDim.x)) {
        const double previous = direction[index] -
                                omegaShared * product[index];
        direction[index] = residual[index] + betaShared * previous;
      }
    }
    __syncthreads();
    residentSpmv(rowPointers, columnIndices, values, rowScale, direction,
                 product, count);
    const double shadowProduct =
        blockDotProduct(shadowResidual, product, count, partial);
    if (lane == 0) {
      if (!isfinite(shadowProduct) ||
          fabs(shadowProduct) <= 2.2204460492503131e-16) {
        doneShared = 1;
        resultShared = -iteration - 1;
      } else {
        alphaShared = rhoShared / shadowProduct;
      }
    }
    __syncthreads();
    if (doneShared)
      break;

    for (int index = static_cast<int>(lane); index < count;
         index += static_cast<int>(blockDim.x))
      intermediate[index] =
          residual[index] - alphaShared * product[index];
    __syncthreads();
    const double intermediateNormSquared =
        blockDotProduct(intermediate, intermediate, count, partial);
    if (lane == 0 && sqrt(intermediateNormSquared) <= targetResidualShared) {
      doneShared = 1;
      resultShared = iteration + 1;
    }
    __syncthreads();
    if (doneShared) {
      for (int index = static_cast<int>(lane); index < count;
           index += static_cast<int>(blockDim.x))
        solution[index] += alphaShared * direction[index];
      __syncthreads();
      break;
    }

    residentSpmv(rowPointers, columnIndices, values, rowScale, intermediate,
                 intermediateProduct, count);
    const double tDotT =
        blockDotProduct(intermediateProduct, intermediateProduct, count,
                        partial);
    const double tDotS =
        blockDotProduct(intermediateProduct, intermediate, count, partial);
    if (lane == 0) {
      if (!isfinite(tDotT) || tDotT <= 2.2204460492503131e-16 ||
          !isfinite(tDotS)) {
        doneShared = 1;
        resultShared = -iteration - 1;
      } else {
        omegaShared = tDotS / tDotT;
      }
    }
    __syncthreads();
    if (doneShared)
      break;

    for (int index = static_cast<int>(lane); index < count;
         index += static_cast<int>(blockDim.x)) {
      solution[index] += alphaShared * direction[index] +
                         omegaShared * intermediate[index];
      residual[index] = intermediate[index] -
                        omegaShared * intermediateProduct[index];
    }
    __syncthreads();
    const double residualNormSquared =
        blockDotProduct(residual, residual, count, partial);
    if (lane == 0) {
      if (sqrt(residualNormSquared) <= targetResidualShared) {
        doneShared = 1;
        resultShared = iteration + 1;
      }
      rhoPreviousShared = rhoShared;
    }
    __syncthreads();
  }

  if (lane == 0) {
    if (!doneShared)
      resultShared = 0;
    *status = resultShared;
  }
}

void checkHip(hipError_t status, const char *operation) {
  if (status != hipSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             hipGetErrorString(status));
}
} // namespace

void hipClearAndAccumulateIndexedRhs(double *deviceRightHandSide,
                                     std::size_t rowCount,
                                     const double *deviceContributions,
                                     const int *deviceRows,
                                     std::size_t contributionCount,
                                     hipStream_t stream) {
  if (!deviceRightHandSide || rowCount == 0)
    throw std::invalid_argument("HIP RHS buffer must have at least one row");
  if (contributionCount != 0 && (!deviceContributions || !deviceRows))
    throw std::invalid_argument("HIP RHS contributions and rows are required");

  checkHip(hipMemsetAsync(deviceRightHandSide, 0,
                          rowCount * sizeof(double), stream),
           "hipMemsetAsync(right-hand side)");
  if (contributionCount == 0)
    return;

  const auto blocks = static_cast<unsigned int>(
      (contributionCount + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(accumulateIndexedRhsKernel, dim3(blocks),
                     dim3(blockSize), 0, stream, deviceRightHandSide,
                     deviceContributions, deviceRows, contributionCount);
  checkHip(hipGetLastError(), "HIP RHS accumulation kernel launch");
}

void hipMultiplyVectorInPlace(double *deviceVector, const double *deviceScale,
                              std::size_t elementCount, hipStream_t stream) {
  if (elementCount == 0)
    return;
  if (!deviceVector || !deviceScale)
    throw std::invalid_argument("HIP vector and scale buffers are required");
  const auto blocks = static_cast<unsigned int>(
      (elementCount + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(multiplyVectorKernel, dim3(blocks), dim3(blockSize), 0,
                     stream, deviceVector, deviceScale, elementCount);
  checkHip(hipGetLastError(), "HIP vector scaling kernel launch");
}

void hipScatterSparseValues(double *deviceValues, const int *devicePositions,
                            const double *deviceUpdates,
                            std::size_t updateCount, hipStream_t stream) {
  if (updateCount == 0)
    return;
  if (!deviceValues || !devicePositions || !deviceUpdates)
    throw std::invalid_argument("HIP sparse scatter buffers are required");
  const auto blocks = static_cast<unsigned int>(
      (updateCount + blockSize - 1) / blockSize);
  hipLaunchKernelGGL(scatterSparseValuesKernel, dim3(blocks), dim3(blockSize),
                     0, stream, deviceValues, devicePositions, deviceUpdates,
                     updateCount);
  checkHip(hipGetLastError(), "HIP sparse scatter kernel launch");
}

void hipSolveBiCGStabResident(
    const int *deviceRowPointers, const int *deviceColumnIndices,
    const double *deviceValues, const double *deviceRowScale,
    const double *deviceRightHandSide, double *deviceSolution,
    double *deviceResidual, double *deviceShadowResidual,
    double *deviceDirection, double *deviceProduct, double *deviceIntermediate,
    double *deviceIntermediateProduct, int dimension, int maximumIterations,
    double relativeTolerance, int *deviceStatus, hipStream_t stream) {
  if (dimension <= 0 || maximumIterations <= 0 || !deviceRowPointers ||
      !deviceColumnIndices || !deviceValues || !deviceRowScale ||
      !deviceRightHandSide || !deviceSolution || !deviceResidual ||
      !deviceShadowResidual || !deviceDirection || !deviceProduct ||
      !deviceIntermediate || !deviceIntermediateProduct || !deviceStatus)
    throw std::invalid_argument("HIP resident BiCGSTAB buffers are required");
  hipLaunchKernelGGL(
      residentBiCGStabKernel, dim3(1), dim3(blockSize), 0, stream,
      deviceRowPointers, deviceColumnIndices, deviceValues, deviceRowScale,
      deviceRightHandSide, deviceSolution, deviceResidual, deviceShadowResidual,
      deviceDirection, deviceProduct, deviceIntermediate,
      deviceIntermediateProduct, dimension, maximumIterations, relativeTolerance,
      deviceStatus);
  checkHip(hipGetLastError(), "HIP resident BiCGSTAB kernel launch");
}
} // namespace DPsim
