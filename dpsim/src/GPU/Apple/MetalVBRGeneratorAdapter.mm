/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <dpsim/GPU/Apple/MetalVBRGeneratorAdapter.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace DPsim {
using CPS::Real;
namespace {
constexpr char kVBRKernels[] = R"metal(
#include <metal_stdlib>
using namespace metal;

struct StepParams {
  uint count; float time; uint order; uint traceOffset; uint traceWidth;
};

struct Generator {
  float voltageDq0[3];
  float currentDq0[3];
  float transientVoltageDq0[3];
  float conductanceDq0[9];
  uint nodeIndices[3];
  uint stampPositions[9];
  float fixedStampValues[9];
  float conductanceAbc[9];
  float sourceAbc[3];
  float timeStep;
  float baseVoltage;
  float baseCurrent;
  float baseImpedance;
  float baseOmegaElec;
  float baseOmegaMech;
  float inertia;
  float damping;
  float ldTransient;
  float lqTransient;
  float adTransient;
  float bdTransient;
  float aqTransient;
  float bqTransient;
  float dqTransient;
  float fieldVoltage;
  float previousFieldVoltage;
  float mechanicalTorque;
  float previousMechanicalTorque;
  float electricalTorque;
  float rotorSpeed;
  float rotorAngle;
  float loadAngle;
};

inline void makeTransforms(float theta, thread float *abcToDq,
                           thread float *dqToAbc) {
  const float angle[3] = {theta, theta - 2.0f * M_PI_F / 3.0f,
                          theta + 2.0f * M_PI_F / 3.0f};
  for (uint phase = 0; phase < 3; ++phase) {
    const float c = cos(angle[phase]);
    const float s = sin(angle[phase]);
    abcToDq[phase] = (2.0f / 3.0f) * c;
    abcToDq[3 + phase] = -(2.0f / 3.0f) * s;
    abcToDq[6 + phase] = 1.0f / 3.0f;
    dqToAbc[3 * phase] = c;
    dqToAbc[3 * phase + 1] = -s;
    dqToAbc[3 * phase + 2] = 1.0f;
  }
}

inline void multiply3(const thread float *a, const thread float *b,
                      thread float *out) {
  for (uint row = 0; row < 3; ++row)
    for (uint col = 0; col < 3; ++col) {
      float value = 0.0f;
      for (uint k = 0; k < 3; ++k)
        value += a[3 * row + k] * b[3 * k + col];
      out[3 * row + col] = value;
    }
}

inline void multiplyVector3(const thread float *a, const thread float *x,
                            thread float *out) {
  for (uint row = 0; row < 3; ++row) {
    out[row] = 0.0f;
    for (uint col = 0; col < 3; ++col)
      out[row] += a[3 * row + col] * x[col];
  }
}

inline void updateElectricalPreStep(device Generator &g,
                                    device float *matrixValues,
                                    device float *rhs, float time) {
  if (time > 0.0f) {
    g.transientVoltageDq0[0] =
        -g.currentDq0[1] * g.lqTransient + g.voltageDq0[0];
    g.transientVoltageDq0[1] =
        g.currentDq0[0] * g.ldTransient + g.voltageDq0[1];
    g.transientVoltageDq0[2] = 0.0f;
  }

  float abcToDq[9], dqToAbc[9];
  makeTransforms(g.rotorAngle, abcToDq, dqToAbc);
  float gDq0[9];
  for (uint i = 0; i < 9; ++i) gDq0[i] = g.conductanceDq0[i];
  float temp[9], conductanceAbc[9];
  multiply3(dqToAbc, gDq0, temp);
  multiply3(temp, abcToDq, conductanceAbc);

  const float eDq[3] = {
      g.adTransient * g.currentDq0[1] +
          g.bdTransient * g.transientVoltageDq0[0],
      g.aqTransient * g.currentDq0[0] +
          g.bqTransient * g.transientVoltageDq0[1] +
          g.dqTransient * (g.previousFieldVoltage + g.fieldVoltage),
      0.0f};
  float eAbc[3];
  multiplyVector3(dqToAbc, eDq, eAbc);
  for (uint phase = 0; phase < 3; ++phase)
    eAbc[phase] *= g.baseVoltage;
  float sourceAbc[3];
  multiplyVector3(conductanceAbc, eAbc, sourceAbc);

  for (uint row = 0; row < 3; ++row) {
    g.sourceAbc[row] = sourceAbc[row];
    rhs[g.nodeIndices[row]] += sourceAbc[row];
    for (uint col = 0; col < 3; ++col) {
      g.conductanceAbc[3 * row + col] = conductanceAbc[3 * row + col];
      matrixValues[g.stampPositions[3 * row + col]] =
          g.fixedStampValues[3 * row + col] + conductanceAbc[3 * row + col];
    }
  }
}

kernel void generatorPreStep(device Generator *generators [[buffer(0)]],
                             device float *matrixValues [[buffer(1)]],
                             device float *rhs [[buffer(2)]],
                             constant StepParams &step [[buffer(3)]],
                             uint gid [[thread_position_in_grid]]) {
  if (gid >= step.count) return;
  device Generator &g = generators[gid];
  // Match Base::ReducedOrderSynchronGenerator<Real>::mnaCompPreStep exactly:
  // electrical torque from prior dq state, forward-Euler rotor update, then
  // SG4 electrical history and the angle-dependent abc stamp.
  g.electricalTorque = g.voltageDq0[0] * g.currentDq0[0] +
                       g.voltageDq0[1] * g.currentDq0[1];
  g.rotorSpeed += g.timeStep * (1.0f / (2.0f * g.inertia)) *
      (g.previousMechanicalTorque - g.electricalTorque -
       g.damping * (g.rotorSpeed - 1.0f));
  g.rotorAngle += g.timeStep * (g.rotorSpeed * g.baseOmegaMech);
  g.loadAngle += g.timeStep * (g.rotorSpeed - 1.0f) * g.baseOmegaMech;
  updateElectricalPreStep(g, matrixValues, rhs, step.time);
}

kernel void generatorPostStep(device Generator *generators [[buffer(0)]],
                              device const float *solution [[buffer(1)]],
                              constant StepParams &step [[buffer(2)]],
                              device float *traceValues [[buffer(3)]],
                              device const uint *solverStatus [[buffer(4)]],
                              device uint *iterationValues [[buffer(5)]],
                              uint gid [[thread_position_in_grid]]) {
  if (gid >= step.count) return;
  device Generator &g = generators[gid];
  float abcToDq[9], dqToAbc[9];
  makeTransforms(g.rotorAngle, abcToDq, dqToAbc);
  float conductanceAbc[9];
  for (uint i = 0; i < 9; ++i) conductanceAbc[i] = g.conductanceAbc[i];
  float terminalVoltage[3];
  for (uint phase = 0; phase < 3; ++phase) {
    terminalVoltage[phase] = solution[g.nodeIndices[phase]];
    if (step.traceWidth != 0)
      traceValues[step.traceOffset + gid * 3 + phase] = terminalVoltage[phase];
  }
  if (gid == 0 && step.traceWidth != 0)
    iterationValues[step.traceOffset / step.traceWidth] = solverStatus[1];
  float voltageDq0[3];
  multiplyVector3(abcToDq, terminalVoltage, voltageDq0);
  for (uint i = 0; i < 3; ++i)
    g.voltageDq0[i] = voltageDq0[i] / g.baseVoltage;

  float conductanceCurrent[3];
  multiplyVector3(conductanceAbc, terminalVoltage, conductanceCurrent);
  float terminalCurrent[3];
  for (uint phase = 0; phase < 3; ++phase)
    terminalCurrent[phase] = g.sourceAbc[phase] - conductanceCurrent[phase];
  float currentDq0[3];
  multiplyVector3(abcToDq, terminalCurrent, currentDq0);
  for (uint i = 0; i < 3; ++i)
    g.currentDq0[i] = currentDq0[i] / g.baseCurrent;
}
)metal";

struct StepParams {
  uint32_t count;
  float time;
  uint32_t order;
  uint32_t traceOffset;
  uint32_t traceWidth;
};

struct DeviceGenerator {
  float voltageDq0[3];
  float currentDq0[3];
  float transientVoltageDq0[3];
  float conductanceDq0[9];
  uint32_t nodeIndices[3];
  uint32_t stampPositions[9];
  float fixedStampValues[9];
  float conductanceAbc[9];
  float sourceAbc[3];
  float timeStep;
  float baseVoltage;
  float baseCurrent;
  float baseImpedance;
  float baseOmegaElec;
  float baseOmegaMech;
  float inertia;
  float damping;
  float ldTransient;
  float lqTransient;
  float adTransient;
  float bdTransient;
  float aqTransient;
  float bqTransient;
  float dqTransient;
  float fieldVoltage;
  float previousFieldVoltage;
  float mechanicalTorque;
  float previousMechanicalTorque;
  float electricalTorque;
  float rotorSpeed;
  float rotorAngle;
  float loadAngle;
};

id<MTLBuffer> makeBuffer(id<MTLDevice> device, const void *data,
                         std::size_t bytes) {
  if (bytes == 0) bytes = sizeof(float);
  id<MTLBuffer> buffer = [device newBufferWithLength:bytes
                                              options:MTLResourceStorageModeShared];
  if (!buffer)
    throw std::runtime_error("MetalVBRGeneratorAdapter buffer allocation failed");
  if (data) std::memcpy(buffer.contents, data, bytes);
  return buffer;
}

void dispatchRows(id<MTLComputeCommandEncoder> encoder,
                  id<MTLComputePipelineState> pipeline, NSUInteger rows) {
  const NSUInteger width = std::min<NSUInteger>(128,
                                                 pipeline.maxTotalThreadsPerThreadgroup);
  [encoder setComputePipelineState:pipeline];
  [encoder dispatchThreads:MTLSizeMake(rows, 1, 1)
      threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
}

bool finiteData(const MetalVBRGeneratorAdapter::Model::MetalDeviceData &model) {
  const Real values[] = {
      model.timeStep, model.baseVoltage, model.baseCurrent, model.baseImpedance,
      model.baseOmegaElec, model.baseOmegaMech, model.inertia, model.damping,
      model.ldTransient, model.lqTransient, model.adTransient, model.bdTransient,
      model.aqTransient, model.bqTransient, model.dqTransient, model.fieldVoltage,
      model.previousFieldVoltage, model.mechanicalTorque,
      model.previousMechanicalTorque, model.electricalTorque, model.rotorSpeed,
      model.rotorAngle, model.loadAngle};
  return std::all_of(std::begin(values), std::end(values),
                     [](Real value) { return std::isfinite(value); }) &&
         std::all_of(model.voltageDq0.begin(), model.voltageDq0.end(),
                     [](Real value) { return std::isfinite(value); }) &&
         std::all_of(model.currentDq0.begin(), model.currentDq0.end(),
                     [](Real value) { return std::isfinite(value); }) &&
         std::all_of(model.transientVoltageDq0.begin(),
                     model.transientVoltageDq0.end(),
                     [](Real value) { return std::isfinite(value); }) &&
         std::all_of(model.conductanceDq0.begin(),
                     model.conductanceDq0.end(),
                     [](Real value) { return std::isfinite(value); }) &&
         std::all_of(model.conductanceAbc.begin(),
                     model.conductanceAbc.end(),
                     [](Real value) { return std::isfinite(value); });
}
} // namespace

struct MetalVBRGeneratorAdapter::Impl {
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLComputePipelineState> prePipeline = nil;
  id<MTLComputePipelineState> postPipeline = nil;
  id<MTLBuffer> generators = nil;
  id<MTLBuffer> dummyTrace = nil;
  id<MTLBuffer> dummyStatus = nil;
  id<MTLBuffer> dummyIterations = nil;
  NSUInteger count = 0;

  ~Impl() {
    [dummyIterations release];
    [dummyStatus release];
    [dummyTrace release];
    [generators release];
    [postPipeline release];
    [prePipeline release];
    [queue release];
    [device release];
  }

  id<MTLComputePipelineState> makePipeline(id<MTLLibrary> library,
                                           NSString *name) {
    id<MTLFunction> function = [library newFunctionWithName:name];
    if (!function)
      throw std::runtime_error(std::string("Missing Metal generator kernel ") +
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
      throw std::runtime_error(std::string("Metal generator pipeline failed: ") +
                               (error.localizedDescription.UTF8String ?: "unknown error"));
    return pipeline;
  }
};

MetalVBRGeneratorAdapter::MetalVBRGeneratorAdapter(
    const std::vector<Model::MetalDeviceData> &models,
    const std::vector<StampPositions> &stampPositions)
    : mImpl(std::make_unique<Impl>()) {
  if (models.empty() || models.size() != stampPositions.size())
    throw std::invalid_argument("MetalVBRGeneratorAdapter requires matching non-empty generator data and stamp slots");
  std::set<UInt> occupiedNodes;
  std::vector<DeviceGenerator> packed(models.size());
  for (std::size_t i = 0; i < models.size(); ++i) {
    const auto &model = models[i];
    if (!finiteData(model) || model.timeStep <= 0. || model.baseVoltage == 0. ||
        model.baseCurrent == 0. || model.inertia <= 0.)
      throw std::invalid_argument("MetalVBRGeneratorAdapter received invalid SG4 data");
    for (UInt phase = 0; phase < 3; ++phase)
      if (!occupiedNodes.insert(model.matrixNodeIndices[phase]).second)
        throw std::invalid_argument("Metal VBR path requires distinct generator buses");

    auto &out = packed[i];
    for (UInt j = 0; j < 3; ++j) {
      out.voltageDq0[j] = static_cast<float>(model.voltageDq0[j]);
      out.currentDq0[j] = static_cast<float>(model.currentDq0[j]);
      out.transientVoltageDq0[j] =
          static_cast<float>(model.transientVoltageDq0[j]);
      out.nodeIndices[j] = static_cast<uint32_t>(model.matrixNodeIndices[j]);
    }
    for (UInt j = 0; j < 9; ++j) {
      out.conductanceDq0[j] = static_cast<float>(model.conductanceDq0[j]);
      out.conductanceAbc[j] = static_cast<float>(model.conductanceAbc[j]);
      out.stampPositions[j] = stampPositions[i].positions[j];
      out.fixedStampValues[j] = stampPositions[i].fixedValues[j];
    }
    out.timeStep = static_cast<float>(model.timeStep);
    out.baseVoltage = static_cast<float>(model.baseVoltage);
    out.baseCurrent = static_cast<float>(model.baseCurrent);
    out.baseImpedance = static_cast<float>(model.baseImpedance);
    out.baseOmegaElec = static_cast<float>(model.baseOmegaElec);
    out.baseOmegaMech = static_cast<float>(model.baseOmegaMech);
    out.inertia = static_cast<float>(model.inertia);
    out.damping = static_cast<float>(model.damping);
    out.ldTransient = static_cast<float>(model.ldTransient);
    out.lqTransient = static_cast<float>(model.lqTransient);
    out.adTransient = static_cast<float>(model.adTransient);
    out.bdTransient = static_cast<float>(model.bdTransient);
    out.aqTransient = static_cast<float>(model.aqTransient);
    out.bqTransient = static_cast<float>(model.bqTransient);
    out.dqTransient = static_cast<float>(model.dqTransient);
    out.fieldVoltage = static_cast<float>(model.fieldVoltage);
    out.previousFieldVoltage = static_cast<float>(model.previousFieldVoltage);
    out.mechanicalTorque = static_cast<float>(model.mechanicalTorque);
    out.previousMechanicalTorque =
        static_cast<float>(model.previousMechanicalTorque);
    out.electricalTorque = static_cast<float>(model.electricalTorque);
    out.rotorSpeed = static_cast<float>(model.rotorSpeed);
    out.rotorAngle = static_cast<float>(model.rotorAngle);
    out.loadAngle = static_cast<float>(model.loadAngle);
  }

  mImpl->device = [MTLCreateSystemDefaultDevice() retain];
  if (!mImpl->device)
    throw std::runtime_error("MetalVBRGeneratorAdapter: no Metal device available");
  mImpl->queue = [mImpl->device newCommandQueue];
  if (!mImpl->queue)
    throw std::runtime_error("MetalVBRGeneratorAdapter: command queue creation failed");
  NSError *error = nil;
  MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
  options.mathMode = MTLMathModeSafe;
  id<MTLLibrary> library = [mImpl->device
      newLibraryWithSource:[NSString stringWithUTF8String:kVBRKernels]
                   options:options
                     error:&error];
  [options release];
  if (!library)
    throw std::runtime_error(std::string("Metal VBR shader compile failed: ") +
                             (error.localizedDescription.UTF8String ?: "unknown error"));
  mImpl->prePipeline = mImpl->makePipeline(library, @"generatorPreStep");
  mImpl->postPipeline = mImpl->makePipeline(library, @"generatorPostStep");
  [library release];
  mImpl->generators = makeBuffer(mImpl->device, packed.data(),
                                 packed.size() * sizeof(DeviceGenerator));
  const uint32_t zero = 0;
  mImpl->dummyTrace = makeBuffer(mImpl->device, &zero, sizeof(zero));
  mImpl->dummyStatus = makeBuffer(mImpl->device, &zero, sizeof(zero));
  mImpl->dummyIterations = makeBuffer(mImpl->device, &zero, sizeof(zero));
  mImpl->count = static_cast<NSUInteger>(packed.size());
}

MetalVBRGeneratorAdapter::~MetalVBRGeneratorAdapter() = default;

void MetalVBRGeneratorAdapter::encodePreStep(void *commandBuffer,
                                              void *matrixValuesBuffer,
                                              void *rightSideBuffer,
                                              double time) {
  if (!commandBuffer || !matrixValuesBuffer || !rightSideBuffer)
    throw std::invalid_argument("MetalVBRGeneratorAdapter requires command, matrix, and RHS buffers");
  id<MTLComputeCommandEncoder> encoder =
      [(id<MTLCommandBuffer>)commandBuffer computeCommandEncoder];
  if (!encoder)
    throw std::runtime_error("Metal generator pre-step encoder creation failed");
  const StepParams params{static_cast<uint32_t>(mImpl->count),
                          static_cast<float>(time), 0, 0, 0};
  [encoder setBuffer:mImpl->generators offset:0 atIndex:0];
  [encoder setBuffer:(id<MTLBuffer>)matrixValuesBuffer offset:0 atIndex:1];
  [encoder setBuffer:(id<MTLBuffer>)rightSideBuffer offset:0 atIndex:2];
  [encoder setBytes:&params length:sizeof(params) atIndex:3];
  dispatchRows(encoder, mImpl->prePipeline, mImpl->count);
  [encoder endEncoding];
}

void MetalVBRGeneratorAdapter::encodePostStep(void *commandBuffer,
                                               void *solutionBuffer) {
  if (!commandBuffer || !solutionBuffer)
    throw std::invalid_argument("MetalVBRGeneratorAdapter requires command and solution buffers");
  id<MTLComputeCommandEncoder> encoder =
      [(id<MTLCommandBuffer>)commandBuffer computeCommandEncoder];
  if (!encoder)
    throw std::runtime_error("Metal generator post-step encoder creation failed");
  const StepParams params{static_cast<uint32_t>(mImpl->count), 0.f, 0, 0, 0};
  [encoder setBuffer:mImpl->generators offset:0 atIndex:0];
  [encoder setBuffer:(id<MTLBuffer>)solutionBuffer offset:0 atIndex:1];
  [encoder setBytes:&params length:sizeof(params) atIndex:2];
  [encoder setBuffer:mImpl->dummyTrace offset:0 atIndex:3];
  [encoder setBuffer:mImpl->dummyStatus offset:0 atIndex:4];
  [encoder setBuffer:mImpl->dummyIterations offset:0 atIndex:5];
  dispatchRows(encoder, mImpl->postPipeline, mImpl->count);
  [encoder endEncoding];
}

void MetalVBRGeneratorAdapter::encodePostStepWithTrace(
    void *commandBuffer, void *solutionBuffer, void *traceValuesBuffer,
    void *solverStatusBuffer, void *iterationValuesBuffer,
    uint32_t stepIndex, uint32_t traceWidth) {
  if (!commandBuffer || !solutionBuffer || !traceValuesBuffer ||
      !solverStatusBuffer || !iterationValuesBuffer || traceWidth == 0)
    throw std::invalid_argument("MetalVBRGeneratorAdapter trace post-step requires valid buffers and width");
  id<MTLComputeCommandEncoder> encoder =
      [(id<MTLCommandBuffer>)commandBuffer computeCommandEncoder];
  if (!encoder)
    throw std::runtime_error("Metal generator trace post-step encoder creation failed");
  const StepParams params{static_cast<uint32_t>(mImpl->count), 0.f, 0,
                          stepIndex * traceWidth, traceWidth};
  [encoder setBuffer:mImpl->generators offset:0 atIndex:0];
  [encoder setBuffer:(id<MTLBuffer>)solutionBuffer offset:0 atIndex:1];
  [encoder setBytes:&params length:sizeof(params) atIndex:2];
  [encoder setBuffer:(id<MTLBuffer>)traceValuesBuffer offset:0 atIndex:3];
  [encoder setBuffer:(id<MTLBuffer>)solverStatusBuffer offset:0 atIndex:4];
  [encoder setBuffer:(id<MTLBuffer>)iterationValuesBuffer offset:0 atIndex:5];
  dispatchRows(encoder, mImpl->postPipeline, mImpl->count);
  [encoder endEncoding];
}

std::vector<MetalVBRGeneratorAdapter::Model::MetalDeviceData>
MetalVBRGeneratorAdapter::copyStateToHost() const {
  const auto *packed =
      static_cast<const DeviceGenerator *>(mImpl->generators.contents);
  std::vector<Model::MetalDeviceData> result(mImpl->count);
  for (std::size_t i = 0; i < result.size(); ++i) {
    const auto &in = packed[i];
    auto &out = result[i];
    for (UInt j = 0; j < 3; ++j) {
      out.voltageDq0[j] = in.voltageDq0[j];
      out.currentDq0[j] = in.currentDq0[j];
      out.transientVoltageDq0[j] = in.transientVoltageDq0[j];
      out.matrixNodeIndices[j] = in.nodeIndices[j];
    }
    for (UInt j = 0; j < 9; ++j)
      out.conductanceDq0[j] = in.conductanceDq0[j];
    for (UInt j = 0; j < 9; ++j)
      out.conductanceAbc[j] = in.conductanceAbc[j];
    for (UInt j = 0; j < 3; ++j)
      out.sourceAbc[j] = in.sourceAbc[j];
    out.timeStep = in.timeStep;
    out.baseVoltage = in.baseVoltage;
    out.baseCurrent = in.baseCurrent;
    out.baseImpedance = in.baseImpedance;
    out.baseOmegaElec = in.baseOmegaElec;
    out.baseOmegaMech = in.baseOmegaMech;
    out.inertia = in.inertia;
    out.damping = in.damping;
    out.ldTransient = in.ldTransient;
    out.lqTransient = in.lqTransient;
    out.adTransient = in.adTransient;
    out.bdTransient = in.bdTransient;
    out.aqTransient = in.aqTransient;
    out.bqTransient = in.bqTransient;
    out.dqTransient = in.dqTransient;
    out.fieldVoltage = in.fieldVoltage;
    out.previousFieldVoltage = in.previousFieldVoltage;
    out.mechanicalTorque = in.mechanicalTorque;
    out.previousMechanicalTorque = in.previousMechanicalTorque;
    out.electricalTorque = in.electricalTorque;
    out.rotorSpeed = in.rotorSpeed;
    out.rotorAngle = in.rotorAngle;
    out.loadAngle = in.loadAngle;
  }
  return result;
}

void *MetalVBRGeneratorAdapter::stateBuffer() const { return mImpl->generators; }
void *MetalVBRGeneratorAdapter::commandQueue() const { return mImpl->queue; }
std::size_t MetalVBRGeneratorAdapter::generatorCount() const { return mImpl->count; }
MetalVBRGeneratorAdapter::IndirectResources
MetalVBRGeneratorAdapter::indirectResources() const {
  return {mImpl->prePipeline, mImpl->postPipeline, mImpl->generators,
          static_cast<uint32_t>(mImpl->count)};
}

} // namespace DPsim
