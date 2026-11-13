/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace {
constexpr char kShader[] = R"metal(
#include <metal_stdlib>
#include <metal_command_buffer>
using namespace metal;

kernel void incrementByOne(device atomic_uint *value [[buffer(0)]]) {
  atomic_fetch_add_explicit(value, 1u, memory_order_relaxed);
}
kernel void incrementByTwo(device atomic_uint *value [[buffer(0)]]) {
  atomic_fetch_add_explicit(value, 2u, memory_order_relaxed);
}
struct CommandArguments {
  compute_pipeline_state pipelineOne [[id(0)]];
  compute_pipeline_state pipelineTwo [[id(1)]];
  device uint *value [[id(2)]];
  command_buffer commands [[id(3)]];
};
kernel void buildCommand(device CommandArguments &arguments [[buffer(0)]],
                         uint index [[thread_position_in_grid]]) {
  if (index >= 2) return;
  compute_command command(arguments.commands, index);
  command.reset();
  command.set_compute_pipeline_state(index == 0 ? arguments.pipelineOne
                                                  : arguments.pipelineTwo);
  if (index == 1) command.set_barrier();
  command.set_kernel_buffer(arguments.value, 0);
  command.concurrent_dispatch_threadgroups(uint3(1, 1, 1), uint3(1, 1, 1));
}
)metal";
} // namespace

int main() {
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  try {
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device)
      throw std::runtime_error("No Metal device available");
    NSError *error = nil;
    MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
    id<MTLLibrary> library = [device newLibraryWithSource:
        [NSString stringWithUTF8String:kShader] options:options error:&error];
    [options release];
    if (!library)
      throw std::runtime_error(error.localizedDescription.UTF8String);

    id<MTLFunction> markOneFunction = [library newFunctionWithName:@"incrementByOne"];
    id<MTLFunction> markTwoFunction = [library newFunctionWithName:@"incrementByTwo"];
    id<MTLFunction> buildFunction = [library newFunctionWithName:@"buildCommand"];
    MTLComputePipelineDescriptor *descriptor =
        [[MTLComputePipelineDescriptor alloc] init];
    descriptor.computeFunction = markOneFunction;
    descriptor.supportIndirectCommandBuffers = YES;
    id<MTLComputePipelineState> markOnePipeline =
        [device newComputePipelineStateWithDescriptor:descriptor
                                               options:0 reflection:nil error:&error];
    descriptor.computeFunction = markTwoFunction;
    id<MTLComputePipelineState> markTwoPipeline =
        [device newComputePipelineStateWithDescriptor:descriptor
                                               options:0 reflection:nil error:&error];
    [descriptor release];
    id<MTLComputePipelineState> buildPipeline =
        [device newComputePipelineStateWithFunction:buildFunction error:&error];
    if (!markOnePipeline || !markTwoPipeline || !buildPipeline)
      throw std::runtime_error(error.localizedDescription.UTF8String);

    MTLIndirectCommandBufferDescriptor *icbDescriptor =
        [[MTLIndirectCommandBufferDescriptor alloc] init];
    icbDescriptor.commandTypes = MTLIndirectCommandTypeConcurrentDispatch;
    icbDescriptor.inheritPipelineState = NO;
    icbDescriptor.inheritBuffers = NO;
    icbDescriptor.maxKernelBufferBindCount = 1;
    id<MTLIndirectCommandBuffer> indirect =
        [device newIndirectCommandBufferWithDescriptor:icbDescriptor
                                        maxCommandCount:2
                                               options:MTLResourceStorageModeShared];
    [icbDescriptor release];
    if (!indirect)
      throw std::runtime_error("Could not allocate indirect command buffer");

    uint32_t initial = 0;
    id<MTLBuffer> result = [device newBufferWithLength:sizeof(initial)
                                               options:MTLResourceStorageModeShared];
    std::memcpy(result.contents, &initial, sizeof(initial));
    id<MTLCommandQueue> queue = [device newCommandQueue];
    id<MTLCommandBuffer> commandBuffer = [queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
    id<MTLArgumentEncoder> argumentEncoder =
        [buildFunction newArgumentEncoderWithBufferIndex:0];
    id<MTLBuffer> arguments = [device newBufferWithLength:argumentEncoder.encodedLength
                                                  options:MTLResourceStorageModeShared];
    [argumentEncoder setArgumentBuffer:arguments offset:0];
    [argumentEncoder setComputePipelineState:markOnePipeline atIndex:0];
    [argumentEncoder setComputePipelineState:markTwoPipeline atIndex:1];
    [argumentEncoder setBuffer:result offset:0 atIndex:2];
    [argumentEncoder setIndirectCommandBuffer:indirect atIndex:3];
    [encoder setComputePipelineState:buildPipeline];
    [encoder setBuffer:arguments offset:0 atIndex:0];
    [encoder useResource:indirect usage:MTLResourceUsageWrite];
    [encoder useResource:result usage:MTLResourceUsageRead];
    [encoder dispatchThreadgroups:MTLSizeMake(1, 1, 1)
             threadsPerThreadgroup:MTLSizeMake(2, 1, 1)];
    [encoder endEncoding];

    encoder = [commandBuffer computeCommandEncoder];
    [encoder useResource:indirect usage:MTLResourceUsageRead];
    [encoder useResource:result usage:MTLResourceUsageWrite];
    [encoder executeCommandsInBuffer:indirect withRange:NSMakeRange(0, 2)];
    [encoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    if (commandBuffer.error)
      throw std::runtime_error(commandBuffer.error.localizedDescription.UTF8String);
    const uint32_t output = *static_cast<const uint32_t *>(result.contents);
    if (output != 3)
      throw std::runtime_error("GPU-generated indirect dispatch output was " +
                               std::to_string(output));
    std::cout << "GPU-generated Metal indirect dispatches with selected pipelines passed; output=" << output
              << '\n';

    [arguments release];
    [argumentEncoder release];
    [result release];
    [queue release];
    [indirect release];
    [buildPipeline release];
    [markTwoPipeline release];
    [markOnePipeline release];
    [buildFunction release];
    [markTwoFunction release];
    [markOneFunction release];
    [library release];
    [device release];
    [pool drain];
    return 0;
  } catch (const std::exception &exception) {
    std::cerr << exception.what() << '\n';
    [pool drain];
    return 1;
  }
}
