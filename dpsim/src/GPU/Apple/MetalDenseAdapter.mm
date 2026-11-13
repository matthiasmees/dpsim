/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>

#include <dpsim/GPU/Apple/MetalDenseAdapter.h>

#include <stdexcept>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <vector>

namespace DPsim {

namespace {
struct ScopedAutoreleasePool {
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  ~ScopedAutoreleasePool() { [pool drain]; }
};
} // namespace

struct MetalDenseAdapter::Impl {
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLBuffer> factorBuffer = nil;
  id<MTLBuffer> pivotBuffer = nil;
  id<MTLBuffer> statusBuffer = nil;
  id<MTLBuffer> rhsBuffer = nil;
  id<MTLBuffer> solutionBuffer = nil;
  MPSMatrix *factors = nil;
  MPSMatrix *pivots = nil;
  MPSMatrix *rhsMatrix = nil;
  MPSMatrix *solutionMatrix = nil;
  MPSMatrixDecompositionLU *decomposition = nil;
  MPSMatrixSolveLU *solve = nil;
  std::size_t order = 0;
  std::vector<float> denseHost;
  std::vector<float> rhsHost;
  bool factorizationPending = false;

  ~Impl() {
    [factors release];
    [pivots release];
    [rhsMatrix release];
    [solutionMatrix release];
    [decomposition release];
    [solve release];
    [factorBuffer release];
    [pivotBuffer release];
    [statusBuffer release];
    [rhsBuffer release];
    [solutionBuffer release];
    [queue release];
    [device release];
  }

  void ensureWorkspace(id<MTLDevice> device, std::size_t n) {
    if (order == n && factors && pivots && rhsMatrix && solutionMatrix &&
        decomposition && solve)
      return;
    [factors release];
    [pivots release];
    [rhsMatrix release];
    [solutionMatrix release];
    [decomposition release];
    [solve release];
    [factorBuffer release];
    [pivotBuffer release];
    [statusBuffer release];
    [rhsBuffer release];
    [solutionBuffer release];

    order = n;
    denseHost.resize(n * n);
    rhsHost.resize(n);
    const NSUInteger rowBytes = static_cast<NSUInteger>(n * sizeof(float));
    const NSUInteger denseBytes = static_cast<NSUInteger>(n * n * sizeof(float));
    const NSUInteger vectorBytes = static_cast<NSUInteger>(n * sizeof(float));
    factorBuffer = [device newBufferWithLength:denseBytes
                                        options:MTLResourceStorageModeShared];
    pivotBuffer = [device newBufferWithLength:n * sizeof(uint32_t)
                                       options:MTLResourceStorageModeShared];
    statusBuffer = [device newBufferWithLength:sizeof(int32_t)
                                        options:MTLResourceStorageModeShared];
    rhsBuffer = [device newBufferWithLength:vectorBytes
                                     options:MTLResourceStorageModeShared];
    solutionBuffer = [device newBufferWithLength:vectorBytes
                                          options:MTLResourceStorageModeShared];
    if (!factorBuffer || !pivotBuffer || !statusBuffer || !rhsBuffer ||
        !solutionBuffer)
      throw std::runtime_error("MetalDenseAdapter: failed to allocate Metal buffers");

    MPSMatrixDescriptor *factorDescriptor =
        [MPSMatrixDescriptor matrixDescriptorWithRows:n
                                             columns:n
                                            rowBytes:rowBytes
                                           dataType:MPSDataTypeFloat32];
    factors = [[MPSMatrix alloc] initWithBuffer:factorBuffer
                                     descriptor:factorDescriptor];
    MPSMatrixDescriptor *pivotDescriptor =
        [MPSMatrixDescriptor matrixDescriptorWithRows:n
                                             columns:1
                                            rowBytes:sizeof(uint32_t)
                                           dataType:MPSDataTypeUInt32];
    pivots = [[MPSMatrix alloc] initWithBuffer:pivotBuffer
                                    descriptor:pivotDescriptor];
    MPSMatrixDescriptor *vectorDescriptor =
        [MPSMatrixDescriptor matrixDescriptorWithRows:n
                                             columns:1
                                            rowBytes:sizeof(float)
                                           dataType:MPSDataTypeFloat32];
    rhsMatrix = [[MPSMatrix alloc] initWithBuffer:rhsBuffer
                                       descriptor:vectorDescriptor];
    solutionMatrix = [[MPSMatrix alloc] initWithBuffer:solutionBuffer
                                            descriptor:vectorDescriptor];
    decomposition = [[MPSMatrixDecompositionLU alloc]
        initWithDevice:device rows:n columns:n];
    solve = [[MPSMatrixSolveLU alloc] initWithDevice:device
                                           transpose:NO
                                                order:n
                              numberOfRightHandSides:1];
    if (!factors || !pivots || !rhsMatrix || !solutionMatrix ||
        !decomposition || !solve)
      throw std::runtime_error("MetalDenseAdapter: failed to create MPS LU objects");
  }
};

MetalDenseAdapter::MetalDenseAdapter(CPS::Logger::Log log)
    : DirectLinearSolver(log), mImpl(std::make_unique<Impl>()) {
  mImpl->device = [MTLCreateSystemDefaultDevice() retain];
  if (!mImpl->device)
    throw std::runtime_error("MetalDenseAdapter: no Metal device available");
  mImpl->queue = [mImpl->device newCommandQueue];
  if (!mImpl->queue)
    throw std::runtime_error("MetalDenseAdapter: unable to create command queue");
}

MetalDenseAdapter::~MetalDenseAdapter() = default;

void MetalDenseAdapter::preprocessing(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> & /*listVariableSystemMatrixEntries*/) {
  if (systemMatrix.rows() != systemMatrix.cols())
    throw std::invalid_argument("MetalDenseAdapter requires a square matrix");
  mImpl->order = 0;
  mImpl->factorizationPending = false;
}

void MetalDenseAdapter::factorize(SparseMatrix &systemMatrix) {
  ScopedAutoreleasePool autoreleasePool;
  const std::size_t n = static_cast<std::size_t>(systemMatrix.rows());
  if (n == 0 || systemMatrix.cols() != systemMatrix.rows())
    throw std::invalid_argument("MetalDenseAdapter requires a non-empty square matrix");
  mImpl->ensureWorkspace(mImpl->device, n);

  auto &dense = mImpl->denseHost;
  std::fill(dense.begin(), dense.end(), 0.0f);
  for (Int row = 0; row < systemMatrix.outerSize(); ++row) {
    for (SparseMatrix::InnerIterator it(systemMatrix, row); it; ++it)
      dense[static_cast<std::size_t>(it.row()) * n +
            static_cast<std::size_t>(it.col())] = static_cast<float>(it.value());
  }

  std::memcpy(mImpl->factorBuffer.contents, dense.data(),
              n * n * sizeof(float));
  // Defer encoding until solveInto so the numeric factorization and solve can
  // share one command buffer and one CPU/GPU synchronization per time step.
  mImpl->factorizationPending = true;

}

void MetalDenseAdapter::refactorize(SparseMatrix &systemMatrix) {
  factorize(systemMatrix);
}

void MetalDenseAdapter::partialRefactorize(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> & /*listVariableSystemMatrixEntries*/) {
  // MPS does not expose KLU-style numeric partial refactorization. This
  // experimental adapter therefore performs a complete dense LU each update.
  factorize(systemMatrix);
}

Matrix MetalDenseAdapter::solve(Matrix &rightSideVector) {
  Matrix result;
  solveInto(rightSideVector, result);
  return result;
}

void MetalDenseAdapter::solveInto(Matrix &rightSideVector, Matrix &solution) {
  ScopedAutoreleasePool autoreleasePool;
  const std::size_t n = mImpl->order;
  if (!mImpl->factors || rightSideVector.size() != static_cast<Eigen::Index>(n))
    throw std::runtime_error("MetalDenseAdapter: solver not factorized or RHS size mismatch");

  auto &rhs = mImpl->rhsHost;
  for (std::size_t i = 0; i < n; ++i)
    rhs[i] = static_cast<float>(rightSideVector.data()[i]);
  std::memcpy(mImpl->rhsBuffer.contents, rhs.data(), n * sizeof(float));

  id<MTLCommandBuffer> command = [mImpl->queue commandBuffer];
  if (!command)
    throw std::runtime_error("MetalDenseAdapter: failed to create command buffer");
  if (mImpl->factorizationPending)
    [mImpl->decomposition encodeToCommandBuffer:command
                                  sourceMatrix:mImpl->factors
                                  resultMatrix:mImpl->factors
                                  pivotIndices:mImpl->pivots
                                        status:mImpl->statusBuffer];
  [mImpl->solve encodeToCommandBuffer:command
                        sourceMatrix:mImpl->factors
                rightHandSideMatrix:mImpl->rhsMatrix
                       pivotIndices:mImpl->pivots
                     solutionMatrix:mImpl->solutionMatrix];
  [command commit];
  [command waitUntilCompleted];
  if (command.error)
    throw std::runtime_error("MetalDenseAdapter: Metal solve command failed");
  if (mImpl->factorizationPending) {
    const auto status =
        *static_cast<const int32_t *>(mImpl->statusBuffer.contents);
    if (status != MPSMatrixDecompositionStatusSuccess)
      throw std::runtime_error("MetalDenseAdapter: MPS LU factorization failed");
    mImpl->factorizationPending = false;
  }

  solution.resize(static_cast<Eigen::Index>(n), 1);
  const auto *values = static_cast<const float *>(mImpl->solutionBuffer.contents);
  for (std::size_t i = 0; i < n; ++i)
    solution.data()[i] = static_cast<Real>(values[i]);

}

} // namespace DPsim
