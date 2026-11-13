/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <dpsim/GPU/Apple/MetalLinearSSNGroupAdapter.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace DPsim {
using namespace CPS;
namespace {
constexpr char kLinearGroupKernels[] = R"metal(
#include <metal_stdlib>
using namespace metal;

struct Params { uint stateCount; uint terminalCount; };
kernel void groupHistory(device const uint *offsets [[buffer(0)]],
                         device const uint *columns [[buffer(1)]],
                         device const float *valuesHi [[buffer(2)]],
                         device const float *valuesLo [[buffer(3)]],
                         device const float *stateHi [[buffer(4)]],
                         device const float *stateLo [[buffer(5)]],
                         device float *history [[buffer(6)]],
                         constant Params &p [[buffer(7)]],
                         device const uint *terminalIndices [[buffer(8)]],
                         device float *rhs [[buffer(9)]],
                         uint globalIndex [[thread_position_in_grid]],
                         uint lane [[thread_index_in_simdgroup]]) {
  const uint row = globalIndex / 32;
  if (row >= p.terminalCount) return;
  float sum = 0.0f;
  for (uint k = offsets[row] + lane; k < offsets[row + 1]; k += 32) {
    sum = fma(valuesHi[k], stateHi[columns[k]], sum);
  }
  sum = simd_sum(sum);
  if (lane == 0) {
    history[row] = sum;
    // Each retained matrix index is unique, so rows assemble into disjoint
    // RHS slots without atomics or a separate serial scatter dispatch.
    rhs[terminalIndices[row]] = sum;
  }
}

kernel void groupStateUpdate(device const uint *stateOffsets [[buffer(0)]],
                             device const uint *stateColumns [[buffer(1)]],
                             device const float *stateValuesHi [[buffer(2)]],
                             device const float *stateValuesLo [[buffer(3)]],
                             device const uint *inputOffsets [[buffer(4)]],
                             device const uint *inputColumns [[buffer(5)]],
                             device const float *inputValuesHi [[buffer(6)]],
                             device const float *inputValuesLo [[buffer(7)]],
                             device const uint *terminalIndices [[buffer(8)]],
                             device const float *solution [[buffer(9)]],
                             device const float *stateHi [[buffer(10)]],
                             device const float *stateLo [[buffer(11)]],
                             constant Params &p [[buffer(12)]],
                             device float *nextStateHi [[buffer(13)]],
                             device float *nextStateLo [[buffer(14)]],
                             uint globalIndex [[thread_position_in_grid]],
                             uint lane [[thread_index_in_simdgroup]]) {
  const uint row = globalIndex / 32;
  if (row >= p.stateCount) return;
  float next = 0.0f;
  for (uint k = stateOffsets[row] + lane; k < stateOffsets[row + 1]; k += 32) {
    next = fma(stateValuesHi[k], stateHi[stateColumns[k]], next);
  }
  for (uint k = inputOffsets[row] + lane; k < inputOffsets[row + 1]; k += 32) {
    next = fma(inputValuesHi[k],
               solution[terminalIndices[inputColumns[k]]], next);
  }
  next = simd_sum(next);
  if (lane == 0) {
    nextStateHi[row] = next;
    nextStateLo[row] = 0.0f;
  }
}
)metal";

struct Params {
  uint32_t stateCount;
  uint32_t terminalCount;
};

struct CsrFloat {
  std::vector<uint32_t> offsets;
  std::vector<uint32_t> columns;
  std::vector<float> values;
  std::vector<float> low;
};

CsrFloat toCsrFloat(const SparseMatrixRow &matrix) {
  SparseMatrixRow compressed = matrix;
  compressed.makeCompressed();
  CsrFloat result;
  result.offsets.resize(static_cast<std::size_t>(compressed.rows()) + 1);
  result.columns.resize(static_cast<std::size_t>(compressed.nonZeros()));
  result.values.resize(static_cast<std::size_t>(compressed.nonZeros()));
  result.low.resize(static_cast<std::size_t>(compressed.nonZeros()));
  for (Eigen::Index i = 0; i <= compressed.rows(); ++i)
    result.offsets[static_cast<std::size_t>(i)] =
        static_cast<uint32_t>(compressed.outerIndexPtr()[i]);
  for (Eigen::Index k = 0; k < compressed.nonZeros(); ++k) {
    result.columns[static_cast<std::size_t>(k)] =
        static_cast<uint32_t>(compressed.innerIndexPtr()[k]);
    const float high = static_cast<float>(compressed.valuePtr()[k]);
    result.values[static_cast<std::size_t>(k)] = high;
    result.low[static_cast<std::size_t>(k)] = static_cast<float>(
        compressed.valuePtr()[k] - static_cast<Real>(high));
  }
  return result;
}

id<MTLBuffer> makeBuffer(id<MTLDevice> device, const void *data,
                         std::size_t bytes) {
  if (bytes == 0)
    bytes = sizeof(float);
  id<MTLBuffer> buffer = [device newBufferWithLength:bytes
                                              options:MTLResourceStorageModeShared];
  if (!buffer)
    throw std::runtime_error("MetalLinearSSNGroupAdapter buffer allocation failed");
  if (data && bytes > 0)
    std::memcpy(buffer.contents, data, bytes);
  return buffer;
}

NSUInteger threadgroupWidth(id<MTLComputePipelineState> pipeline) {
  return std::min<NSUInteger>(128, pipeline.maxTotalThreadsPerThreadgroup);
}

void dispatchRows(id<MTLComputeCommandEncoder> encoder,
                  id<MTLComputePipelineState> pipeline, NSUInteger rows) {
  const NSUInteger width = threadgroupWidth(pipeline);
  [encoder setComputePipelineState:pipeline];
  [encoder dispatchThreads:MTLSizeMake(rows * 32, 1, 1)
      threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
}
} // namespace

struct MetalLinearSSNGroupAdapter::Impl {
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLComputePipelineState> historyPipeline = nil;
  id<MTLComputePipelineState> updatePipeline = nil;
  id<MTLBuffer> historyOffsets = nil;
  id<MTLBuffer> historyColumns = nil;
  id<MTLBuffer> historyValues = nil;
  id<MTLBuffer> historyValuesLow = nil;
  id<MTLBuffer> stateOffsets = nil;
  id<MTLBuffer> stateColumns = nil;
  id<MTLBuffer> stateValues = nil;
  id<MTLBuffer> stateValuesLow = nil;
  id<MTLBuffer> inputOffsets = nil;
  id<MTLBuffer> inputColumns = nil;
  id<MTLBuffer> inputValues = nil;
  id<MTLBuffer> inputValuesLow = nil;
  id<MTLBuffer> terminalIndices = nil;
  id<MTLBuffer> state = nil;
  id<MTLBuffer> stateLow = nil;
  id<MTLBuffer> nextState = nil;
  id<MTLBuffer> nextStateLow = nil;
  id<MTLBuffer> history = nil;
  id<MTLBuffer> terminalAdmittance = nil;
  id<MTLBuffer> paramsBuffer = nil;
  Params params{};

  ~Impl() {
    [paramsBuffer release];
    [terminalAdmittance release];
    [history release];
    [nextState release];
    [nextStateLow release];
    [state release];
    [stateLow release];
    [terminalIndices release];
    [inputValuesLow release]; [inputValues release]; [inputColumns release]; [inputOffsets release];
    [stateValuesLow release]; [stateValues release]; [stateColumns release]; [stateOffsets release];
    [historyValuesLow release]; [historyValues release]; [historyColumns release]; [historyOffsets release];
    [updatePipeline release]; [historyPipeline release];
    [queue release]; [device release];
  }

  id<MTLComputePipelineState> makePipeline(id<MTLLibrary> library,
                                           NSString *name) {
    id<MTLFunction> function = [library newFunctionWithName:name];
    if (!function)
      throw std::runtime_error(std::string("Missing Metal kernel ") +
                               name.UTF8String);
    NSError *error = nil;
    MTLComputePipelineDescriptor *descriptor =
        [[MTLComputePipelineDescriptor alloc] init];
    descriptor.computeFunction = function;
    descriptor.supportIndirectCommandBuffers = YES;
    id<MTLComputePipelineState> pipeline =
        [device newComputePipelineStateWithDescriptor:descriptor
                                              options:0
                                           reflection:nil
                                                error:&error];
    [descriptor release];
    [function release];
    if (!pipeline)
      throw std::runtime_error(std::string("Metal kernel pipeline failed: ") +
                               (error.localizedDescription.UTF8String ?: "unknown error"));
    return pipeline;
  }
};

MetalLinearSSNGroupAdapter::MetalLinearSSNGroupAdapter(
    const CPS::EMT::LinearSSNGroup::PrecomputedStateMaps &model)
    : mImpl(std::make_unique<Impl>()) {
  if (model.historyMap.rows() <= 0 || model.stateMap.rows() <= 0 ||
      model.inputMap.rows() != model.stateMap.rows() ||
      model.historyMap.cols() != model.stateMap.rows() ||
      model.inputMap.cols() != model.terminalAdmittance.rows() ||
      model.terminalAdmittance.rows() != model.terminalAdmittance.cols() ||
      model.initialState.rows() != model.stateMap.rows() ||
      model.terminalMatrixIndices.size() !=
          static_cast<std::size_t>(model.historyMap.rows()))
    throw std::invalid_argument("MetalLinearSSNGroupAdapter invalid realization dimensions");

  mImpl->device = [MTLCreateSystemDefaultDevice() retain];
  if (!mImpl->device)
    throw std::runtime_error("MetalLinearSSNGroupAdapter: no Metal device available");
  mImpl->queue = [mImpl->device newCommandQueue];
  if (!mImpl->queue)
    throw std::runtime_error("MetalLinearSSNGroupAdapter: command queue creation failed");

  NSError *error = nil;
  MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
  options.mathMode = MTLMathModeSafe;
  id<MTLLibrary> library = [mImpl->device
      newLibraryWithSource:[NSString stringWithUTF8String:kLinearGroupKernels]
                   options:options
                     error:&error];
  [options release];
  if (!library)
    throw std::runtime_error(std::string("Metal linear-group shader compile failed: ") +
                             (error.localizedDescription.UTF8String ?: "unknown error"));
  mImpl->historyPipeline = mImpl->makePipeline(library, @"groupHistory");
  mImpl->updatePipeline = mImpl->makePipeline(library, @"groupStateUpdate");
  [library release];

  const CsrFloat history = toCsrFloat(model.historyMap);
  const CsrFloat state = toCsrFloat(model.stateMap);
  const CsrFloat input = toCsrFloat(model.inputMap);
  const auto reportMap = [](const char *name, const SparseMatrixRow &matrix,
                            const CsrFloat &csr) {
    const double entries = static_cast<double>(matrix.rows()) * matrix.cols();
    std::cout << "Metal SSN group " << name << " map: " << matrix.rows()
              << 'x' << matrix.cols() << ", nnz=" << matrix.nonZeros()
              << ", density=" << (entries > 0.0 ? matrix.nonZeros() / entries : 0.0)
              << ", avg row nnz="
              << (matrix.rows() > 0 ? static_cast<double>(csr.values.size()) / matrix.rows() : 0.0)
              << '\n';
  };
  reportMap("history", model.historyMap, history);
  reportMap("state", model.stateMap, state);
  reportMap("input", model.inputMap, input);
  auto &impl = *mImpl;
  impl.historyOffsets = makeBuffer(impl.device, history.offsets.data(),
                                   history.offsets.size() * sizeof(uint32_t));
  impl.historyColumns = makeBuffer(impl.device, history.columns.data(),
                                   history.columns.size() * sizeof(uint32_t));
  impl.historyValues = makeBuffer(impl.device, history.values.data(),
                                  history.values.size() * sizeof(float));
  impl.historyValuesLow = makeBuffer(impl.device, history.low.data(),
                                      history.low.size() * sizeof(float));
  impl.stateOffsets = makeBuffer(impl.device, state.offsets.data(),
                                 state.offsets.size() * sizeof(uint32_t));
  impl.stateColumns = makeBuffer(impl.device, state.columns.data(),
                                 state.columns.size() * sizeof(uint32_t));
  impl.stateValues = makeBuffer(impl.device, state.values.data(),
                                state.values.size() * sizeof(float));
  impl.stateValuesLow = makeBuffer(impl.device, state.low.data(),
                                    state.low.size() * sizeof(float));
  impl.inputOffsets = makeBuffer(impl.device, input.offsets.data(),
                                 input.offsets.size() * sizeof(uint32_t));
  impl.inputColumns = makeBuffer(impl.device, input.columns.data(),
                                 input.columns.size() * sizeof(uint32_t));
  impl.inputValues = makeBuffer(impl.device, input.values.data(),
                                input.values.size() * sizeof(float));
  impl.inputValuesLow = makeBuffer(impl.device, input.low.data(),
                                    input.low.size() * sizeof(float));
  std::vector<uint32_t> terminalIndices(model.terminalMatrixIndices.size());
  std::transform(model.terminalMatrixIndices.begin(),
                 model.terminalMatrixIndices.end(), terminalIndices.begin(),
                 [](UInt index) { return static_cast<uint32_t>(index); });
  std::vector<bool> terminalIndexSeen(terminalIndices.size(), false);
  for (uint32_t index : terminalIndices) {
    if (index >= terminalIndexSeen.size() || terminalIndexSeen[index])
      throw std::invalid_argument(
          "MetalLinearSSNGroupAdapter requires terminal indices to cover the compact RHS exactly once");
    terminalIndexSeen[index] = true;
  }
  impl.terminalIndices = makeBuffer(impl.device, terminalIndices.data(),
                                   terminalIndices.size() * sizeof(uint32_t));
  std::vector<float> initialState(static_cast<std::size_t>(model.initialState.rows()));
  std::vector<float> initialStateLow(initialState.size());
  for (Eigen::Index i = 0; i < model.initialState.rows(); ++i)
  {
    const float high = static_cast<float>(model.initialState(i, 0));
    initialState[static_cast<std::size_t>(i)] = high;
    initialStateLow[static_cast<std::size_t>(i)] = static_cast<float>(
        model.initialState(i, 0) - static_cast<Real>(high));
  }
  impl.state = makeBuffer(impl.device, initialState.data(),
                          initialState.size() * sizeof(float));
  impl.stateLow = makeBuffer(impl.device, initialStateLow.data(),
                             initialStateLow.size() * sizeof(float));
  impl.nextState = makeBuffer(impl.device, nullptr,
                              initialState.size() * sizeof(float));
  impl.nextStateLow = makeBuffer(impl.device, nullptr,
                                 initialState.size() * sizeof(float));
  impl.history = makeBuffer(impl.device, nullptr,
                            terminalIndices.size() * sizeof(float));
  std::vector<float> admittance(static_cast<std::size_t>(
      model.terminalAdmittance.rows() * model.terminalAdmittance.cols()));
  for (Eigen::Index row = 0; row < model.terminalAdmittance.rows(); ++row)
    for (Eigen::Index col = 0; col < model.terminalAdmittance.cols(); ++col)
      admittance[static_cast<std::size_t>(row * model.terminalAdmittance.cols() + col)] =
          static_cast<float>(model.terminalAdmittance(row, col));
  impl.terminalAdmittance = makeBuffer(impl.device, admittance.data(),
                                       admittance.size() * sizeof(float));
  impl.params.stateCount = static_cast<uint32_t>(model.stateMap.rows());
  impl.params.terminalCount = static_cast<uint32_t>(model.historyMap.rows());
  impl.paramsBuffer = makeBuffer(impl.device, &impl.params, sizeof(impl.params));
}

MetalLinearSSNGroupAdapter::~MetalLinearSSNGroupAdapter() = default;

void MetalLinearSSNGroupAdapter::encodeHistoryAndRhs(void *commandBuffer,
                                                      void *rhsBuffer) {
  if (!commandBuffer || !rhsBuffer)
    throw std::invalid_argument("MetalLinearSSNGroupAdapter requires command and RHS buffers");
  id<MTLCommandBuffer> command = (id<MTLCommandBuffer>)commandBuffer;
  id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
  if (!encoder)
    throw std::runtime_error("MetalLinearSSNGroupAdapter history encoder creation failed");
  [encoder setBuffer:mImpl->historyOffsets offset:0 atIndex:0];
  [encoder setBuffer:mImpl->historyColumns offset:0 atIndex:1];
  [encoder setBuffer:mImpl->historyValues offset:0 atIndex:2];
  [encoder setBuffer:mImpl->historyValuesLow offset:0 atIndex:3];
  [encoder setBuffer:mImpl->state offset:0 atIndex:4];
  [encoder setBuffer:mImpl->stateLow offset:0 atIndex:5];
  [encoder setBuffer:mImpl->history offset:0 atIndex:6];
  [encoder setBuffer:mImpl->paramsBuffer offset:0 atIndex:7];
  [encoder setBuffer:mImpl->terminalIndices offset:0 atIndex:8];
  [encoder setBuffer:(id<MTLBuffer>)rhsBuffer offset:0 atIndex:9];
  dispatchRows(encoder, mImpl->historyPipeline, mImpl->params.terminalCount);
  [encoder endEncoding];
}

void MetalLinearSSNGroupAdapter::encodeStateUpdate(void *commandBuffer,
                                                   void *solutionBuffer) {
  if (!commandBuffer || !solutionBuffer)
    throw std::invalid_argument("MetalLinearSSNGroupAdapter requires command and solution buffers");
  id<MTLComputeCommandEncoder> encoder =
      [(id<MTLCommandBuffer>)commandBuffer computeCommandEncoder];
  if (!encoder)
    throw std::runtime_error("MetalLinearSSNGroupAdapter state encoder creation failed");
  [encoder setBuffer:mImpl->stateOffsets offset:0 atIndex:0];
  [encoder setBuffer:mImpl->stateColumns offset:0 atIndex:1];
  [encoder setBuffer:mImpl->stateValues offset:0 atIndex:2];
  [encoder setBuffer:mImpl->stateValuesLow offset:0 atIndex:3];
  [encoder setBuffer:mImpl->inputOffsets offset:0 atIndex:4];
  [encoder setBuffer:mImpl->inputColumns offset:0 atIndex:5];
  [encoder setBuffer:mImpl->inputValues offset:0 atIndex:6];
  [encoder setBuffer:mImpl->inputValuesLow offset:0 atIndex:7];
  [encoder setBuffer:mImpl->terminalIndices offset:0 atIndex:8];
  [encoder setBuffer:(id<MTLBuffer>)solutionBuffer offset:0 atIndex:9];
  [encoder setBuffer:mImpl->state offset:0 atIndex:10];
  [encoder setBuffer:mImpl->stateLow offset:0 atIndex:11];
  [encoder setBuffer:mImpl->paramsBuffer offset:0 atIndex:12];
  [encoder setBuffer:mImpl->nextState offset:0 atIndex:13];
  [encoder setBuffer:mImpl->nextStateLow offset:0 atIndex:14];
  dispatchRows(encoder, mImpl->updatePipeline, mImpl->params.stateCount);
  [encoder endEncoding];
  // Future timestep kernels must read the just-encoded output. Metal preserves
  // command ordering, so rotate the handles immediately without waiting.
  std::swap(mImpl->state, mImpl->nextState);
  std::swap(mImpl->stateLow, mImpl->nextStateLow);
}

void *MetalLinearSSNGroupAdapter::stateBuffer() const { return mImpl->state; }
void *MetalLinearSSNGroupAdapter::historyBuffer() const { return mImpl->history; }
void *MetalLinearSSNGroupAdapter::terminalIndicesBuffer() const { return mImpl->terminalIndices; }
void *MetalLinearSSNGroupAdapter::terminalAdmittanceBuffer() const { return mImpl->terminalAdmittance; }
void *MetalLinearSSNGroupAdapter::commandQueue() const { return mImpl->queue; }
std::size_t MetalLinearSSNGroupAdapter::stateCount() const { return mImpl->params.stateCount; }
std::size_t MetalLinearSSNGroupAdapter::terminalCount() const { return mImpl->params.terminalCount; }
MetalLinearSSNGroupAdapter::IndirectResources
MetalLinearSSNGroupAdapter::indirectResources() const {
  return {mImpl->historyPipeline, mImpl->updatePipeline,
          mImpl->historyOffsets, mImpl->historyColumns, mImpl->historyValues,
          mImpl->historyValuesLow, mImpl->stateOffsets, mImpl->stateColumns,
          mImpl->stateValues, mImpl->stateValuesLow, mImpl->inputOffsets,
          mImpl->inputColumns, mImpl->inputValues, mImpl->inputValuesLow,
          mImpl->terminalIndices, mImpl->state, mImpl->stateLow,
          mImpl->nextState, mImpl->nextStateLow, mImpl->history,
          mImpl->paramsBuffer, mImpl->params.stateCount,
          mImpl->params.terminalCount};
}
void MetalLinearSSNGroupAdapter::finishIndirectStateUpdates(
    std::size_t updateCount) {
  if (updateCount & 1u) {
    std::swap(mImpl->state, mImpl->nextState);
    std::swap(mImpl->stateLow, mImpl->nextStateLow);
  }
}
std::vector<float> MetalLinearSSNGroupAdapter::copyStateToHost() const {
  const auto *state = static_cast<const float *>(mImpl->state.contents);
  const auto *low = static_cast<const float *>(mImpl->stateLow.contents);
  std::vector<float> result(mImpl->params.stateCount);
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = state[i] + low[i];
  return result;
}

} // namespace DPsim
