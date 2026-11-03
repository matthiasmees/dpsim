/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>

#include <dpsim/GPU/Apple/MetalSparseLUAdapter.h>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace DPsim {
namespace {
constexpr char kExpandCSR[] = R"metal(
#include <metal_stdlib>
using namespace metal;
kernel void expandCSR(device const uint *rowOffsets [[buffer(0)]],
                     device const uint *columns [[buffer(1)]],
                     device const float *values [[buffer(2)]],
                     device float *dense [[buffer(3)]],
                     constant uint &order [[buffer(4)]],
                     constant uint &rowStride [[buffer(5)]],
                     uint row [[thread_position_in_grid]]) {
  if (row >= order) return;
  device float *denseRow = dense + row * rowStride;
  for (uint col = 0; col < order; ++col) denseRow[col] = 0.0f;
  for (uint k = rowOffsets[row]; k < rowOffsets[row + 1]; ++k)
    denseRow[columns[k]] = values[k];
}
)metal";
}

struct MetalSparseLUAdapter::Impl {
  id<MTLDevice> device = nil;
  id<MTLBuffer> rowOffsets = nil;
  id<MTLBuffer> columns = nil;
  id<MTLBuffer> factorBuffer = nil;
  id<MTLBuffer> pivotBuffer = nil;
  id<MTLBuffer> statusBuffer = nil;
  id<MTLComputePipelineState> expandPipeline = nil;
  MPSMatrix *factors = nil;
  MPSMatrix *pivots = nil;
  MPSMatrixDecompositionLU *decomposition = nil;
  MPSMatrixSolveLU *solve = nil;
  uint32_t order = 0;
  uint32_t rowStrideFloats = 0;

  ~Impl() {
    [factors release];
    [pivots release];
    [decomposition release];
    [solve release];
    [expandPipeline release];
    [rowOffsets release];
    [columns release];
    [factorBuffer release];
    [pivotBuffer release];
    [statusBuffer release];
    [device release];
  }
};

MetalSparseLUAdapter::MetalSparseLUAdapter(
    void *devicePointer, uint32_t order,
    const std::vector<uint32_t> &rowOffsets,
    const std::vector<uint32_t> &columnIndices)
    : mImpl(std::make_unique<Impl>()) {
  if (!devicePointer || order == 0 || rowOffsets.size() != order + 1 ||
      rowOffsets.back() != columnIndices.size())
    throw std::invalid_argument("MetalSparseLUAdapter received invalid CSR structure");
  mImpl->device = [(id<MTLDevice>)devicePointer retain];
  mImpl->order = order;
  mImpl->rowStrideFloats = (order + 15u) & ~15u;
  mImpl->rowOffsets = [mImpl->device newBufferWithBytes:rowOffsets.data()
                                                 length:rowOffsets.size() * sizeof(uint32_t)
                                                options:MTLResourceStorageModeShared];
  mImpl->columns = [mImpl->device newBufferWithBytes:columnIndices.data()
                                              length:std::max<std::size_t>(1, columnIndices.size()) * sizeof(uint32_t)
                                             options:MTLResourceStorageModeShared];
  const std::size_t denseBytes = static_cast<std::size_t>(order) *
                                 mImpl->rowStrideFloats * sizeof(float);
  mImpl->factorBuffer = [mImpl->device newBufferWithLength:denseBytes
                                                   options:MTLResourceStorageModeShared];
  mImpl->pivotBuffer = [mImpl->device newBufferWithLength:order * sizeof(uint32_t)
                                                  options:MTLResourceStorageModeShared];
  mImpl->statusBuffer = [mImpl->device newBufferWithLength:sizeof(int32_t)
                                                   options:MTLResourceStorageModeShared];
  if (!mImpl->rowOffsets || !mImpl->columns || !mImpl->factorBuffer ||
      !mImpl->pivotBuffer || !mImpl->statusBuffer)
    throw std::runtime_error("MetalSparseLUAdapter failed to allocate GPU buffers");

  NSError *error = nil;
  MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
  id<MTLLibrary> library = [mImpl->device
      newLibraryWithSource:[NSString stringWithUTF8String:kExpandCSR]
                   options:options error:&error];
  [options release];
  if (!library)
    throw std::runtime_error(std::string("MetalSparseLUAdapter shader compilation failed: ") +
                             (error.localizedDescription.UTF8String ?: "unknown error"));
  id<MTLFunction> function = [library newFunctionWithName:@"expandCSR"];
  mImpl->expandPipeline = [mImpl->device newComputePipelineStateWithFunction:function
                                                                     error:&error];
  [function release];
  [library release];
  if (!mImpl->expandPipeline)
    throw std::runtime_error(std::string("MetalSparseLUAdapter pipeline creation failed: ") +
                             (error.localizedDescription.UTF8String ?: "unknown error"));

  MPSMatrixDescriptor *factorDescriptor =
      [MPSMatrixDescriptor matrixDescriptorWithRows:order
                                           columns:order
                                          rowBytes:static_cast<NSUInteger>(mImpl->rowStrideFloats) * sizeof(float)
                                         dataType:MPSDataTypeFloat32];
  mImpl->factors = [[MPSMatrix alloc] initWithBuffer:mImpl->factorBuffer
                                          descriptor:factorDescriptor];
  MPSMatrixDescriptor *pivotDescriptor =
      [MPSMatrixDescriptor matrixDescriptorWithRows:order columns:1
                                           rowBytes:sizeof(uint32_t)
                                          dataType:MPSDataTypeUInt32];
  mImpl->pivots = [[MPSMatrix alloc] initWithBuffer:mImpl->pivotBuffer
                                        descriptor:pivotDescriptor];
  mImpl->decomposition = [[MPSMatrixDecompositionLU alloc]
      initWithDevice:mImpl->device rows:order columns:order];
  mImpl->solve = [[MPSMatrixSolveLU alloc] initWithDevice:mImpl->device
                                                transpose:NO order:order
                               numberOfRightHandSides:1];
  if (!mImpl->factors || !mImpl->pivots || !mImpl->decomposition || !mImpl->solve)
    throw std::runtime_error("MetalSparseLUAdapter failed to create MPS LU objects");
}

MetalSparseLUAdapter::~MetalSparseLUAdapter() = default;

void MetalSparseLUAdapter::encodeSolve(void *commandBufferPointer,
                                       void *sparseValuesPointer,
                                       void *rhsPointer,
                                       void *solutionPointer) {
  id<MTLCommandBuffer> command = (id<MTLCommandBuffer>)commandBufferPointer;
  id<MTLBuffer> sparseValues = (id<MTLBuffer>)sparseValuesPointer;
  id<MTLBuffer> rhsBuffer = (id<MTLBuffer>)rhsPointer;
  id<MTLBuffer> solutionBuffer = (id<MTLBuffer>)solutionPointer;
  if (!command || !sparseValues || !rhsBuffer || !solutionBuffer)
    throw std::invalid_argument("MetalSparseLUAdapter encodeSolve received a null resource");
  id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
  if (!encoder) throw std::runtime_error("MetalSparseLUAdapter could not create expansion encoder");
  const uint32_t order = mImpl->order;
  const uint32_t rowStride = mImpl->rowStrideFloats;
  [encoder setComputePipelineState:mImpl->expandPipeline];
  [encoder setBuffer:mImpl->rowOffsets offset:0 atIndex:0];
  [encoder setBuffer:mImpl->columns offset:0 atIndex:1];
  [encoder setBuffer:sparseValues offset:0 atIndex:2];
  [encoder setBuffer:mImpl->factorBuffer offset:0 atIndex:3];
  [encoder setBytes:&order length:sizeof(order) atIndex:4];
  [encoder setBytes:&rowStride length:sizeof(rowStride) atIndex:5];
  const NSUInteger threads = std::min<NSUInteger>(
      std::max<NSUInteger>(1, mImpl->expandPipeline.maxTotalThreadsPerThreadgroup), 128);
  [encoder dispatchThreads:MTLSizeMake(order, 1, 1)
      threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
  [encoder endEncoding];

  MPSMatrixDescriptor *vectorDescriptor =
      [MPSMatrixDescriptor matrixDescriptorWithRows:mImpl->order columns:1
                                           rowBytes:sizeof(float)
                                          dataType:MPSDataTypeFloat32];
  MPSMatrix *rhs = [[MPSMatrix alloc] initWithBuffer:rhsBuffer descriptor:vectorDescriptor];
  MPSMatrix *solution = [[MPSMatrix alloc] initWithBuffer:solutionBuffer descriptor:vectorDescriptor];
  if (!rhs || !solution) {
    [rhs release]; [solution release];
    throw std::runtime_error("MetalSparseLUAdapter could not wrap solve vectors");
  }
  [mImpl->decomposition encodeToCommandBuffer:command
                                sourceMatrix:mImpl->factors
                                resultMatrix:mImpl->factors
                                  pivotIndices:mImpl->pivots
                                        status:mImpl->statusBuffer];
  [mImpl->solve encodeToCommandBuffer:command sourceMatrix:mImpl->factors
                rightHandSideMatrix:rhs pivotIndices:mImpl->pivots
                      solutionMatrix:solution];
  [rhs release];
  [solution release];
}

void MetalSparseLUAdapter::checkStatus() const {
  const auto status = *static_cast<const int32_t *>(mImpl->statusBuffer.contents);
  if (status != MPSMatrixDecompositionStatusSuccess)
    throw std::runtime_error("MetalSparseLUAdapter MPS LU factorization failed; status=" +
                             std::to_string(status));
}

} // namespace DPsim
