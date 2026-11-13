/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <dpsim/GPU/Apple/MetalSparseBiCGSTABAdapter.h>

#include <Eigen/SparseCholesky>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace DPsim {
namespace {
constexpr char kBiCGSTABKernel[] = R"metal(
#include <metal_stdlib>
using namespace metal;

struct Params {
  uint n; uint maxIterations; float relativeTolerance; float absoluteTolerance;
  uint denseStorage;
};

inline void reduceSum(threadgroup float *scratch, uint tid, uint threads) {
  threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
  for (uint stride = threads >> 1; stride > 0; stride >>= 1) {
    if (tid < stride) scratch[tid] += scratch[tid + stride];
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
  }
}

inline void reduceSumPair(threadgroup float *scratch, uint tid, uint threads) {
  threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
  for (uint stride = threads >> 1; stride > 0; stride >>= 1) {
    if (tid < stride) {
      scratch[tid] += scratch[tid + stride];
      scratch[threads + tid] += scratch[threads + tid + stride];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
  }
}

inline float csrAt(device const uint *rowOffsets,
                   device const uint *columnIndices,
                   device const float *values, uint row, uint column) {
  for (uint k = rowOffsets[row]; k < rowOffsets[row + 1]; ++k)
    if (columnIndices[k] == column) return values[k];
  return 0.0f;
}

inline float multiplyRow(device const uint *rowOffsets,
                         device const uint *columnIndices,
                         device const float *values, device const float *input,
                         uint row, uint n, bool denseStorage) {
  float sum = 0.0f;
  if (denseStorage) {
    const uint base = row * n;
    for (uint col = 0; col < n; ++col)
      sum += values[base + col] * input[col];
  } else {
    for (uint k = rowOffsets[row]; k < rowOffsets[row + 1]; ++k)
      sum += values[k] * input[columnIndices[k]];
  }
  return sum;
}

inline void buildBlockJacobi(device const uint *rowOffsets,
                             device const uint *columnIndices,
                             device const float *values,
                             device float *inverseBlocks, uint tid, uint n,
                             bool denseStorage) {
  if (tid >= n || tid % 3 != 0) return;
  const uint i = tid;
  const float a00 = denseStorage ? values[i * n + i]
                                  : csrAt(rowOffsets, columnIndices, values, i, i);
  const float a01 = denseStorage ? values[i * n + i + 1]
                                  : csrAt(rowOffsets, columnIndices, values, i, i + 1);
  const float a02 = denseStorage ? values[i * n + i + 2]
                                  : csrAt(rowOffsets, columnIndices, values, i, i + 2);
  const float a10 = denseStorage ? values[(i + 1) * n + i]
                                  : csrAt(rowOffsets, columnIndices, values, i + 1, i);
  const float a11 = denseStorage ? values[(i + 1) * n + i + 1]
                                  : csrAt(rowOffsets, columnIndices, values, i + 1, i + 1);
  const float a12 = denseStorage ? values[(i + 1) * n + i + 2]
                                  : csrAt(rowOffsets, columnIndices, values, i + 1, i + 2);
  const float a20 = denseStorage ? values[(i + 2) * n + i]
                                  : csrAt(rowOffsets, columnIndices, values, i + 2, i);
  const float a21 = denseStorage ? values[(i + 2) * n + i + 1]
                                  : csrAt(rowOffsets, columnIndices, values, i + 2, i + 1);
  const float a22 = denseStorage ? values[(i + 2) * n + i + 2]
                                  : csrAt(rowOffsets, columnIndices, values, i + 2, i + 2);
  const float c00 = a11 * a22 - a12 * a21;
  const float c01 = a02 * a21 - a01 * a22;
  const float c02 = a01 * a12 - a02 * a11;
  const float c10 = a12 * a20 - a10 * a22;
  const float c11 = a00 * a22 - a02 * a20;
  const float c12 = a02 * a10 - a00 * a12;
  const float c20 = a10 * a21 - a11 * a20;
  const float c21 = a01 * a20 - a00 * a21;
  const float c22 = a00 * a11 - a01 * a10;
  const float det = a00 * c00 + a01 * c10 + a02 * c20;
  const uint base = (i / 3) * 9;
  if (abs(det) > 1e-30f && isfinite(det)) {
    inverseBlocks[base] = c00 / det;
    inverseBlocks[base + 1] = c01 / det;
    inverseBlocks[base + 2] = c02 / det;
    inverseBlocks[base + 3] = c10 / det;
    inverseBlocks[base + 4] = c11 / det;
    inverseBlocks[base + 5] = c12 / det;
    inverseBlocks[base + 6] = c20 / det;
    inverseBlocks[base + 7] = c21 / det;
    inverseBlocks[base + 8] = c22 / det;
  } else {
    inverseBlocks[base] = 1.0f;
    inverseBlocks[base + 1] = 0.0f;
    inverseBlocks[base + 2] = 0.0f;
    inverseBlocks[base + 3] = 0.0f;
    inverseBlocks[base + 4] = 1.0f;
    inverseBlocks[base + 5] = 0.0f;
    inverseBlocks[base + 6] = 0.0f;
    inverseBlocks[base + 7] = 0.0f;
    inverseBlocks[base + 8] = 1.0f;
  }
}

inline void applyBlockJacobi(device const float *inverseBlocks,
                             device const float *input,
                             device float *output, uint tid, uint n) {
  if (tid >= n || tid % 3 != 0) return;
  const uint i = tid;
  const uint base = (i / 3) * 9;
  const float r0 = input[i], r1 = input[i + 1], r2 = input[i + 2];
  output[i] = inverseBlocks[base] * r0 + inverseBlocks[base + 1] * r1 +
              inverseBlocks[base + 2] * r2;
  output[i + 1] = inverseBlocks[base + 3] * r0 + inverseBlocks[base + 4] * r1 +
                  inverseBlocks[base + 5] * r2;
  output[i + 2] = inverseBlocks[base + 6] * r0 + inverseBlocks[base + 7] * r1 +
                  inverseBlocks[base + 8] * r2;
}

kernel void bicgstab(device const uint *rowOffsets [[buffer(0)]],
                device const uint *columnIndices [[buffer(1)]],
                device const float *values [[buffer(2)]],
                device const float *rhs [[buffer(3)]],
                device float *x [[buffer(4)]],
                device float *work [[buffer(5)]],
                device uint *status [[buffer(6)]],
                constant Params &params [[buffer(7)]],
                threadgroup float *scratch [[threadgroup(0)]],
                uint tid [[thread_index_in_threadgroup]],
                uint threads [[threads_per_threadgroup]]) {
  const uint n = params.n;
  device float *r = work;
  device float *rhat = work + n;
  device float *p = work + 2 * n;
  device float *v = work + 3 * n;
  device float *s = work + 4 * n;
  device float *t = work + 5 * n;
  device float *phat = work + 6 * n;
  device float *shat = work + 7 * n;
  device float *inverseBlocks = work + 8 * n;
  const uint scalarBase = 2 * threads;
  const uint rhoOldIndex = scalarBase;
  const uint alphaIndex = scalarBase + 1;
  const uint omegaIndex = scalarBase + 2;
  const uint rhoIndex = scalarBase + 3;
  const uint betaIndex = scalarBase + 4;
  const uint tempIndex = scalarBase + 5;
  const uint flagIndex = scalarBase + 6;
  const uint normIndex = scalarBase + 7;

  if (tid == 0) status[1] = 0;
  if (tid < n) {
    const float ax = multiplyRow(rowOffsets, columnIndices, values, x, tid, n,
                                 params.denseStorage != 0);
    r[tid] = rhs[tid] - ax;
    rhat[tid] = r[tid];
    p[tid] = 0.0f;
    v[tid] = 0.0f;
    scratch[tid] = rhs[tid] * rhs[tid];
  } else {
    scratch[tid] = 0.0f;
  }
  buildBlockJacobi(rowOffsets, columnIndices, values, inverseBlocks, tid, n,
                   params.denseStorage != 0);
  reduceSum(scratch, tid, threads);
  if (tid == 0) {
    scratch[normIndex] = sqrt(max(scratch[0], 0.0f));
    scratch[rhoOldIndex] = 1.0f;
    scratch[alphaIndex] = 1.0f;
    scratch[omegaIndex] = 1.0f;
  }
  threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
  if (tid < n) scratch[tid] = r[tid] * r[tid];
  else scratch[tid] = 0.0f;
  reduceSum(scratch, tid, threads);
  if (tid == 0) {
    const float residualNorm = sqrt(max(scratch[0], 0.0f));
    status[2] = as_type<uint>(residualNorm);
    status[3] = as_type<uint>(scratch[normIndex]);
    // Initially rhat == r, so the first BiCGSTAB rho is already the squared
    // residual norm from the reduction immediately above.
    scratch[rhoIndex] = scratch[0];
    scratch[flagIndex] = residualNorm <= max(params.absoluteTolerance,
        params.relativeTolerance * scratch[normIndex]) ? 1.0f : 0.0f;
  }
  threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);

  for (uint iteration = 0; iteration < params.maxIterations; ++iteration) {
    if (scratch[flagIndex] != 0.0f) break;
    if (iteration == 0) {
      // p and v are both zero on the first iteration, making beta immaterial;
      // avoid recomputing rhat dot r and its threadgroup reduction.
      if (tid == 0) {
        if (abs(scratch[rhoIndex]) < 1e-30f)
          scratch[flagIndex] = -1.0f;
        scratch[betaIndex] = 0.0f;
      }
      threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    } else {
      if (tid < n) scratch[tid] = rhat[tid] * r[tid];
      else scratch[tid] = 0.0f;
      reduceSum(scratch, tid, threads);
      if (tid == 0) {
        const float rho = scratch[0];
        if (abs(rho) < 1e-30f || abs(scratch[omegaIndex]) < 1e-30f) {
          scratch[flagIndex] = -1.0f;
        } else {
          scratch[rhoIndex] = rho;
          scratch[betaIndex] = (rho / scratch[rhoOldIndex]) *
                               (scratch[alphaIndex] / scratch[omegaIndex]);
        }
      }
      threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    }
    if (scratch[flagIndex] < 0.0f) break;

    if (tid < n)
      p[tid] = r[tid] + scratch[betaIndex] * (p[tid] - scratch[omegaIndex] * v[tid]);
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    applyBlockJacobi(inverseBlocks, p, phat, tid, n);
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (tid < n) {
      const float sum = multiplyRow(rowOffsets, columnIndices, values, phat,
                                    tid, n, params.denseStorage != 0);
      v[tid] = sum;
      scratch[tid] = rhat[tid] * sum;
    } else {
      scratch[tid] = 0.0f;
    }
    reduceSum(scratch, tid, threads);
    if (tid == 0) {
      const float denominator = scratch[0];
      if (abs(denominator) < 1e-30f || !isfinite(denominator)) {
        scratch[flagIndex] = -1.0f;
      } else {
        scratch[alphaIndex] = scratch[rhoIndex] / denominator;
      }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (scratch[flagIndex] < 0.0f) break;

    if (tid < n) {
      s[tid] = r[tid] - scratch[alphaIndex] * v[tid];
      scratch[tid] = s[tid] * s[tid];
    } else {
      scratch[tid] = 0.0f;
    }
    reduceSum(scratch, tid, threads);
    if (tid == 0) {
      const float sNorm = sqrt(max(scratch[0], 0.0f));
      status[2] = as_type<uint>(sNorm);
      if (sNorm <= max(params.absoluteTolerance,
                       params.relativeTolerance * scratch[normIndex]))
        scratch[flagIndex] = 1.0f;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (scratch[flagIndex] > 0.0f) {
      if (tid < n) x[tid] += scratch[alphaIndex] * phat[tid];
      if (tid == 0) status[1] = iteration + 1;
      threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
      break;
    }
    applyBlockJacobi(inverseBlocks, s, shat, tid, n);
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (tid < n) {
      const float sum = multiplyRow(rowOffsets, columnIndices, values, shat,
                                    tid, n, params.denseStorage != 0);
      t[tid] = sum;
      scratch[tid] = sum * s[tid];
      scratch[threads + tid] = sum * sum;
    } else {
      scratch[tid] = 0.0f;
      scratch[threads + tid] = 0.0f;
    }
    reduceSumPair(scratch, tid, threads);
    if (tid == 0) {
      scratch[tempIndex] = scratch[0];
      const float denominator = scratch[threads];
      if (denominator < 1e-30f || !isfinite(denominator))
        scratch[flagIndex] = -1.0f;
      else
        scratch[omegaIndex] = scratch[tempIndex] / denominator;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
    if (scratch[flagIndex] < 0.0f) break;
    if (tid < n) {
      x[tid] += scratch[alphaIndex] * phat[tid] + scratch[omegaIndex] * shat[tid];
      r[tid] = s[tid] - scratch[omegaIndex] * t[tid];
      scratch[tid] = r[tid] * r[tid];
    } else {
      scratch[tid] = 0.0f;
    }
    reduceSum(scratch, tid, threads);
    if (tid == 0) {
      const float residualNorm = sqrt(max(scratch[0], 0.0f));
      status[2] = as_type<uint>(residualNorm);
      if (residualNorm <= max(params.absoluteTolerance,
                              params.relativeTolerance * scratch[normIndex]))
        scratch[flagIndex] = 1.0f;
      scratch[rhoOldIndex] = scratch[rhoIndex];
      status[1] = iteration + 1;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup | mem_flags::mem_device);
  }
  if (tid == 0 && scratch[flagIndex] < 0.0f) status[0] = 2;
  else if (tid == 0 && scratch[flagIndex] == 0.0f && status[0] == 0)
    status[0] = 1;
}
)metal";

struct KernelParams {
  uint32_t n;
  uint32_t maxIterations;
  float relativeTolerance;
  float absoluteTolerance;
  uint32_t denseStorage;
};

struct ScopedAutoreleasePool {
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  ~ScopedAutoreleasePool() { [pool drain]; }
};

NSUInteger nextPowerOfTwo(NSUInteger x) {
  NSUInteger result = 1;
  while (result < x)
    result <<= 1;
  return result;
}
} // namespace

struct MetalSparseBiCGSTABAdapter::Impl {
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLComputePipelineState> pipeline = nil;
  id<MTLBuffer> rowOffsets = nil;
  id<MTLBuffer> columnIndices = nil;
  id<MTLBuffer> values = nil;
  id<MTLBuffer> rhs = nil;
  id<MTLBuffer> solution = nil;
  id<MTLBuffer> work = nil;
  id<MTLBuffer> status = nil;
  id<MTLBuffer> paramsBuffer = nil;
  std::size_t order = 0;
  std::size_t nonZeros = 0;
  NSUInteger threads = 0;
  KernelParams params{0, 4096, 1e-6f, 1e-7f, 0};

  ~Impl() {
    [status release];
    [paramsBuffer release];
    [work release];
    [solution release];
    [rhs release];
    [values release];
    [columnIndices release];
    [rowOffsets release];
    [pipeline release];
    [queue release];
    [device release];
  }

  void ensureBuffers(std::size_t n, std::size_t nnz) {
    if (order == n && nonZeros == nnz && pipeline)
      return;
    [status release]; status = nil;
    [paramsBuffer release]; paramsBuffer = nil;
    [work release]; work = nil;
    [solution release]; solution = nil;
    [rhs release]; rhs = nil;
    [values release]; values = nil;
    [columnIndices release]; columnIndices = nil;
    [rowOffsets release]; rowOffsets = nil;
    [pipeline release]; pipeline = nil;

    order = n;
    nonZeros = nnz;
    threads = nextPowerOfTwo(static_cast<NSUInteger>(n));
    MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
    options.mathMode = MTLMathModeSafe;
    NSError *error = nil;
    NSString *source = [NSString stringWithUTF8String:kBiCGSTABKernel];
    id<MTLLibrary> library = [device newLibraryWithSource:source options:options error:&error];
    [options release];
    if (!library)
      throw std::runtime_error(std::string("MetalSparseBiCGSTABAdapter shader compile failed: ") +
                               (error.localizedDescription.UTF8String ?: "unknown error"));
    id<MTLFunction> function = [library newFunctionWithName:@"bicgstab"];
    MTLComputePipelineDescriptor *descriptor =
        [[MTLComputePipelineDescriptor alloc] init];
    descriptor.computeFunction = function;
    descriptor.supportIndirectCommandBuffers = YES;
    pipeline = [device newComputePipelineStateWithDescriptor:descriptor
                                                     options:0
                                                  reflection:nil
                                                       error:&error];
    [descriptor release];
    [function release];
    [library release];
    if (!pipeline)
      throw std::runtime_error(std::string("MetalSparseBiCGSTABAdapter pipeline creation failed: ") +
                               (error.localizedDescription.UTF8String ?: "unknown error"));
    if (threads > pipeline.maxTotalThreadsPerThreadgroup)
      throw std::invalid_argument("MetalSparseBiCGSTABAdapter matrix order exceeds one-threadgroup solver limit");

    rowOffsets = [device newBufferWithLength:(n + 1) * sizeof(uint32_t)
                                      options:MTLResourceStorageModeShared];
    columnIndices = [device newBufferWithLength:nnz * sizeof(uint32_t)
                                         options:MTLResourceStorageModeShared];
    values = [device newBufferWithLength:nnz * sizeof(float)
                                  options:MTLResourceStorageModeShared];
    rhs = [device newBufferWithLength:n * sizeof(float)
                               options:MTLResourceStorageModeShared];
    solution = [device newBufferWithLength:n * sizeof(float)
                                    options:MTLResourceStorageModeShared];
    work = [device newBufferWithLength:11 * n * sizeof(float)
                                options:MTLResourceStorageModeShared];
    status = [device newBufferWithLength:4 * sizeof(uint32_t)
                                  options:MTLResourceStorageModeShared];
    params.n = static_cast<uint32_t>(n);
    paramsBuffer = [device newBufferWithLength:sizeof(params)
                                         options:MTLResourceStorageModeShared];
    if (!rowOffsets || !columnIndices || !values || !rhs || !solution || !work ||
        !status || !paramsBuffer)
      throw std::runtime_error("MetalSparseBiCGSTABAdapter failed to allocate buffers");
    std::memcpy(paramsBuffer.contents, &params, sizeof(params));
  }

  void encode(id<MTLCommandBuffer> command, id<MTLBuffer> matrixValues,
              id<MTLBuffer> rightHandSide, id<MTLBuffer> output) {
    if (!pipeline || !command || !matrixValues || !rightHandSide || !output)
      throw std::runtime_error("MetalSparseBiCGSTABAdapter encode called before preprocessing");
    if (matrixValues.length < nonZeros * sizeof(float) ||
        rightHandSide.length < order * sizeof(float) ||
        output.length < order * sizeof(float))
      throw std::invalid_argument("MetalSparseBiCGSTABAdapter device buffer is too small");
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    if (!encoder)
      throw std::runtime_error("MetalSparseBiCGSTABAdapter unable to create compute encoder");
    [encoder setComputePipelineState:pipeline];
    [encoder setBuffer:rowOffsets offset:0 atIndex:0];
    [encoder setBuffer:columnIndices offset:0 atIndex:1];
    [encoder setBuffer:matrixValues offset:0 atIndex:2];
    [encoder setBuffer:rightHandSide offset:0 atIndex:3];
    [encoder setBuffer:output offset:0 atIndex:4];
    [encoder setBuffer:work offset:0 atIndex:5];
    [encoder setBuffer:status offset:0 atIndex:6];
    [encoder setBuffer:paramsBuffer offset:0 atIndex:7];
    [encoder setThreadgroupMemoryLength:(2 * threads + 8) * sizeof(float) atIndex:0];
    [encoder dispatchThreadgroups:MTLSizeMake(1, 1, 1)
             threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
    [encoder endEncoding];
  }
};

MetalSparseBiCGSTABAdapter::MetalSparseBiCGSTABAdapter(CPS::Logger::Log log)
    : DirectLinearSolver(log), mImpl(std::make_unique<Impl>()) {
  mImpl->device = [MTLCreateSystemDefaultDevice() retain];
  if (!mImpl->device)
    throw std::runtime_error("MetalSparseBiCGSTABAdapter: no Metal device available");
  mImpl->queue = [mImpl->device newCommandQueue];
  if (!mImpl->queue)
    throw std::runtime_error("MetalSparseBiCGSTABAdapter: unable to create command queue");
}

MetalSparseBiCGSTABAdapter::~MetalSparseBiCGSTABAdapter() = default;

void MetalSparseBiCGSTABAdapter::setConvergenceCriteria(
    uint32_t maxIterations, float relativeTolerance, float absoluteTolerance) {
  if (mImpl->pipeline || maxIterations == 0 || relativeTolerance <= 0.0f ||
      absoluteTolerance < 0.0f)
    throw std::invalid_argument(
        "MetalSparseBiCGSTABAdapter convergence criteria must be set before preprocessing");
  mImpl->params.maxIterations = maxIterations;
  mImpl->params.relativeTolerance = relativeTolerance;
  mImpl->params.absoluteTolerance = absoluteTolerance;
}

void MetalSparseBiCGSTABAdapter::preprocessing(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &) {
  if (systemMatrix.rows() <= 0 || systemMatrix.rows() != systemMatrix.cols())
    throw std::invalid_argument("MetalSparseBiCGSTABAdapter requires a non-empty square matrix");
  systemMatrix.makeCompressed();
  const auto n = static_cast<std::size_t>(systemMatrix.rows());
  mImpl->params.denseStorage =
      static_cast<std::size_t>(systemMatrix.nonZeros()) == n * n ? 1u : 0u;
  const auto nnz = static_cast<std::size_t>(systemMatrix.nonZeros());
  mImpl->ensureBuffers(n, nnz);

  auto *rowStart = static_cast<uint32_t *>(mImpl->rowOffsets.contents);
  auto *columns = static_cast<uint32_t *>(mImpl->columnIndices.contents);
  std::fill_n(static_cast<uint32_t *>(mImpl->status.contents), 4, 0u);
  std::fill_n(static_cast<float *>(mImpl->solution.contents), n, 0.0f);
  std::fill_n(rowStart, n + 1, 0u);
  const auto *outer = systemMatrix.outerIndexPtr();
  const auto *inner = systemMatrix.innerIndexPtr();
  for (std::size_t i = 0; i <= n; ++i)
    rowStart[i] = static_cast<uint32_t>(outer[i]);
  for (std::size_t k = 0; k < nnz; ++k)
    columns[k] = static_cast<uint32_t>(inner[k]);

  // The live MNA matrix can be nonsymmetric and indefinite (notably with the
  // conventional VBR generator stamps), so use BiCGSTAB rather than assuming
  // the SPD conditions required by CG.
}

void MetalSparseBiCGSTABAdapter::factorize(SparseMatrix &systemMatrix) {
  if (static_cast<std::size_t>(systemMatrix.rows()) != mImpl->order)
    throw std::runtime_error("MetalSparseBiCGSTABAdapter factorize called before preprocessing");
  auto *values = static_cast<float *>(mImpl->values.contents);
  const auto *source = systemMatrix.valuePtr();
  for (std::size_t k = 0; k < mImpl->nonZeros; ++k)
    values[k] = static_cast<float>(source[k]);
}

void MetalSparseBiCGSTABAdapter::refactorize(SparseMatrix &systemMatrix) {
  factorize(systemMatrix);
}

void MetalSparseBiCGSTABAdapter::partialRefactorize(
    SparseMatrix &systemMatrix, std::vector<std::pair<UInt, UInt>> &) {
  factorize(systemMatrix);
}

Matrix MetalSparseBiCGSTABAdapter::solve(Matrix &rightSideVector) {
  Matrix solution;
  solveInto(rightSideVector, solution);
  return solution;
}

void MetalSparseBiCGSTABAdapter::solveInto(Matrix &rightSideVector, Matrix &solution) {
  if (rightSideVector.size() != static_cast<Eigen::Index>(mImpl->order))
    throw std::invalid_argument("MetalSparseBiCGSTABAdapter RHS size mismatch");
  auto *rhs = static_cast<float *>(mImpl->rhs.contents);
  auto *x = static_cast<float *>(mImpl->solution.contents);
  for (std::size_t i = 0; i < mImpl->order; ++i) {
    rhs[i] = static_cast<float>(rightSideVector.data()[i]);
  }
  id<MTLCommandBuffer> command = [mImpl->queue commandBuffer];
  mImpl->encode(command, mImpl->values, mImpl->rhs, mImpl->solution);
  [command commit];
  [command waitUntilCompleted];
  if (command.error)
    throw std::runtime_error("MetalSparseBiCGSTABAdapter command failed");
  const auto *status = static_cast<const uint32_t *>(mImpl->status.contents);
  if (status[0] != 0u) {
    float residual = 0.0f, rhsNorm = 0.0f;
    std::memcpy(&residual, &status[2], sizeof(float));
    std::memcpy(&rhsNorm, &status[3], sizeof(float));
    const auto *rowStart = static_cast<const uint32_t *>(mImpl->rowOffsets.contents);
    const auto *columns = static_cast<const uint32_t *>(mImpl->columnIndices.contents);
    const auto *values = static_cast<const float *>(mImpl->values.contents);
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(mImpl->nonZeros);
    for (std::size_t row = 0; row < mImpl->order; ++row)
      for (uint32_t k = rowStart[row]; k < rowStart[row + 1]; ++k)
        triplets.emplace_back(static_cast<Int>(row), static_cast<Int>(columns[k]),
                              static_cast<Real>(values[k]));
    SparseMatrix debugMatrix(static_cast<Int>(mImpl->order),
                             static_cast<Int>(mImpl->order));
    debugMatrix.setFromTriplets(triplets.begin(), triplets.end());
    Real maxAsymmetry = 0.0;
    for (Int row = 0; row < debugMatrix.rows(); ++row)
      for (SparseMatrix::InnerIterator it(debugMatrix, row); it; ++it)
        maxAsymmetry = std::max<Real>(
            maxAsymmetry,
            std::abs(it.value() - debugMatrix.coeff(it.col(), it.row())));
    Eigen::SimplicialLDLT<SparseMatrix> reference;
    reference.compute(debugMatrix);
    std::string referenceStatus = "LDLT-failed";
    if (reference.info() == Eigen::Success) {
      const Real minPivot = reference.vectorD().minCoeff();
      const Matrix referenceSolution = reference.solve(rightSideVector);
      const Real referenceResidual =
          (rightSideVector - debugMatrix * referenceSolution).norm() /
          std::max<Real>(rightSideVector.norm(), 1e-30);
      referenceStatus = "LDLT-min-pivot=" + std::to_string(minPivot) +
                       ", asymmetry=" + std::to_string(maxAsymmetry) +
                       ", reference-relative-residual=" +
                       std::to_string(referenceResidual);
    }
    throw std::runtime_error(
        "MetalSparseBiCGSTABAdapter failed; status=" +
        std::to_string(status[0]) + ", iterations=" +
        std::to_string(status[1]) + ", relative residual=" +
        std::to_string(residual / std::max(rhsNorm, 1e-30f)) + ", " +
        referenceStatus);
  }
  solution.resize(static_cast<Eigen::Index>(mImpl->order), 1);
  for (std::size_t i = 0; i < mImpl->order; ++i)
    solution.data()[i] = static_cast<Real>(x[i]);
}

void MetalSparseBiCGSTABAdapter::encodeSolve(void *commandBuffer,
                                       void *matrixValuesBuffer,
                                       void *rightSideBuffer,
                                       void *solutionBuffer) {
  mImpl->encode((id<MTLCommandBuffer>)commandBuffer,
                (id<MTLBuffer>)matrixValuesBuffer,
                (id<MTLBuffer>)rightSideBuffer,
                (id<MTLBuffer>)solutionBuffer);
}

void MetalSparseBiCGSTABAdapter::encodeZeroInitialGuess(void *commandBuffer) {
  if (!commandBuffer || !mImpl->solution)
    throw std::invalid_argument(
        "MetalSparseBiCGSTABAdapter zero-guess requested before preprocessing");
  id<MTLBlitCommandEncoder> encoder =
      [(id<MTLCommandBuffer>)commandBuffer blitCommandEncoder];
  if (!encoder)
    throw std::runtime_error(
        "MetalSparseBiCGSTABAdapter could not create zero-guess encoder");
  [encoder fillBuffer:mImpl->solution
                range:NSMakeRange(0, mImpl->solution.length)
                value:0];
  [encoder endEncoding];
}

void *MetalSparseBiCGSTABAdapter::rowOffsetsBuffer() const { return mImpl->rowOffsets; }
void *MetalSparseBiCGSTABAdapter::columnIndicesBuffer() const { return mImpl->columnIndices; }
void *MetalSparseBiCGSTABAdapter::matrixValuesBuffer() const { return mImpl->values; }
void *MetalSparseBiCGSTABAdapter::rightSideBuffer() const { return mImpl->rhs; }
void *MetalSparseBiCGSTABAdapter::solutionBuffer() const { return mImpl->solution; }
void *MetalSparseBiCGSTABAdapter::statusBuffer() const { return mImpl->status; }
void *MetalSparseBiCGSTABAdapter::commandQueue() const { return mImpl->queue; }
MetalSparseBiCGSTABAdapter::IndirectResources
MetalSparseBiCGSTABAdapter::indirectResources() const {
  return {mImpl->pipeline, mImpl->work, mImpl->status, mImpl->paramsBuffer,
          static_cast<uint32_t>(mImpl->order),
          static_cast<uint32_t>(mImpl->threads)};
}

} // namespace DPsim
