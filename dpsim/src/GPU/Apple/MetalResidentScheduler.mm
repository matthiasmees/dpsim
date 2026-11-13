/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <dpsim/GPU/Apple/MetalLinearSSNGroupAdapter.h>
#include <dpsim/GPU/Apple/MetalResidentScheduler.h>
#include <dpsim/GPU/Apple/MetalSparseLUAdapter.h>
#include <dpsim/GPU/Apple/MetalSparseBiCGSTABAdapter.h>
#include <dpsim/GPU/Apple/MetalVBRGeneratorAdapter.h>

#include <Eigen/LU>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace DPsim {
namespace {
constexpr char kResidentKernels[] = R"metal(
#include <metal_stdlib>
using namespace metal;
kernel void clearRhs(device float *rhs [[buffer(0)]],
                     constant uint &n [[buffer(1)]],
                     uint i [[thread_position_in_grid]]) {
  if (i < n) rhs[i] = 0.0f;
}
kernel void captureTrace(device const float *solution [[buffer(0)]],
                         device const uint *indices [[buffer(1)]],
                         device float *trace [[buffer(2)]],
                         constant uint &offset [[buffer(3)]],
                         device const uint *solverStatus [[buffer(4)]],
                         device uint *iterationTrace [[buffer(5)]],
                         constant uint &step [[buffer(6)]],
                         uint i [[thread_position_in_grid]]) {
  if (i == 0) iterationTrace[step] = solverStatus[1];
  trace[offset + i] = solution[indices[i]];
}

)metal";

constexpr char kResidentIcbKernels[] = R"metal(
#include <metal_stdlib>
#include <metal_command_buffer>
using namespace metal;

struct StepInfo {
  uint count; float time; uint order; uint traceOffset; uint traceWidth;
};
struct BuildInfo {
  uint commandCount; uint firstStep; uint generatorCount; uint order;
  uint traceWidth; uint clearGroups; uint historyGroups; uint traceGroups;
  uint stateGroups; uint generatorGroups; uint solverThreads;
  float startTime; float timeStep;
};

struct IcbArguments {
  compute_pipeline_state clearPipeline [[id(0)]];
  compute_pipeline_state historyPipeline [[id(1)]];
  compute_pipeline_state generatorPrePipeline [[id(3)]];
  compute_pipeline_state solverPipeline [[id(5)]];
  compute_pipeline_state statePipeline [[id(7)]];
  compute_pipeline_state generatorPostPipeline [[id(8)]];
  command_buffer commands [[id(9)]];
  device float *rhs [[id(10)]];
  device const uint *csrOffsets [[id(11)]];
  device const uint *csrColumns [[id(12)]];
  device float *matrixValues [[id(13)]];
  device float *solution [[id(14)]];
  device float *solverWork [[id(15)]];
  device uint *solverStatus [[id(16)]];
  constant uint *solverParams [[id(17)]];
  device const uint *historyOffsets [[id(18)]];
  device const uint *historyColumns [[id(19)]];
  device const float *historyValues [[id(20)]];
  device const float *historyValuesLow [[id(21)]];
  device const uint *stateOffsets [[id(22)]];
  device const uint *stateColumns [[id(23)]];
  device const float *stateValues [[id(24)]];
  device const float *stateValuesLow [[id(25)]];
  device const uint *inputOffsets [[id(26)]];
  device const uint *inputColumns [[id(27)]];
  device const float *inputValues [[id(28)]];
  device const float *inputValuesLow [[id(29)]];
  device const uint *terminalIndices [[id(30)]];
  device float *state [[id(31)]];
  device float *stateLow [[id(32)]];
  device float *nextState [[id(33)]];
  device float *nextStateLow [[id(34)]];
  device float *alternateState [[id(35)]];
  device float *alternateStateLow [[id(36)]];
  device float *alternateNextState [[id(37)]];
  device float *alternateNextStateLow [[id(38)]];
  device float *history [[id(39)]];
  constant uint *groupParams [[id(40)]];
  device uchar *generators [[id(41)]];
  device StepInfo *stepInfosWrite [[id(42)]];
  constant StepInfo *stepInfosRead [[id(43)]];
  device const uint *traceIndices [[id(44)]];
  device float *traceValues [[id(45)]];
  device uint *iterationValues [[id(46)]];
};

kernel void clearRhsForIcb(device float *rhs [[buffer(0)]],
                           constant StepInfo &step [[buffer(1)]],
                           uint i [[thread_position_in_grid]]) {
  if (i < step.order) rhs[i] = 0.0f;
}
kernel void buildTimestepCommands(device IcbArguments &a [[buffer(0)]],
                                  constant BuildInfo &b [[buffer(1)]],
                                  uint commandIndex [[thread_position_in_grid]]) {
  if (commandIndex >= b.commandCount) return;
  const uint op = commandIndex % 5;
  const uint localStep = commandIndex / 5;
  const uint stepNumber = b.firstStep + localStep;
  if (op == 0) {
    a.stepInfosWrite[localStep] = {
      b.generatorCount, b.startTime + float(stepNumber) * b.timeStep,
      b.order, stepNumber * b.traceWidth, b.traceWidth};
  }

  compute_command c(a.commands, commandIndex);
  c.reset();
  // Serialize at actual data dependencies. The post-solve trace, group state,
  // and generator state updates are independent and may overlap.
  if (commandIndex != 0 &&
      (op == 0 || op == 1 || op == 2 || op == 3))
    c.set_barrier();
  const constant StepInfo *step = a.stepInfosRead + localStep;
  const uint solverThreads = b.solverThreads;
  if (op == 0) {
    const bool alternate = (localStep & 1) != 0;
    c.set_compute_pipeline_state(a.historyPipeline);
    c.set_kernel_buffer(a.historyOffsets, 0);
    c.set_kernel_buffer(a.historyColumns, 1);
    c.set_kernel_buffer(a.historyValues, 2);
    c.set_kernel_buffer(a.historyValuesLow, 3);
    c.set_kernel_buffer(alternate ? a.alternateState : a.state, 4);
    c.set_kernel_buffer(alternate ? a.alternateStateLow : a.stateLow, 5);
    c.set_kernel_buffer(a.history, 6);
    c.set_kernel_buffer(a.groupParams, 7);
    c.set_kernel_buffer(a.terminalIndices, 8);
    c.set_kernel_buffer(a.rhs, 9);
    c.concurrent_dispatch_threadgroups(uint3(b.historyGroups, 1, 1), uint3(128, 1, 1));
  } else if (op == 1) {
    c.set_compute_pipeline_state(a.generatorPrePipeline);
    c.set_kernel_buffer(a.generators, 0);
    c.set_kernel_buffer(a.matrixValues, 1);
    c.set_kernel_buffer(a.rhs, 2);
    c.set_kernel_buffer(step, 3);
    c.concurrent_dispatch_threadgroups(uint3(b.generatorGroups, 1, 1), uint3(128, 1, 1));
  } else if (op == 2) {
    c.set_compute_pipeline_state(a.solverPipeline);
    c.set_kernel_buffer(a.csrOffsets, 0);
    c.set_kernel_buffer(a.csrColumns, 1);
    c.set_kernel_buffer(a.matrixValues, 2);
    c.set_kernel_buffer(a.rhs, 3);
    c.set_kernel_buffer(a.solution, 4);
    c.set_kernel_buffer(a.solverWork, 5);
    c.set_kernel_buffer(a.solverStatus, 6);
    c.set_kernel_buffer(a.solverParams, 7);
    c.set_threadgroup_memory_length((2 * solverThreads + 8) * sizeof(float), 0);
    c.concurrent_dispatch_threadgroups(uint3(1, 1, 1), uint3(solverThreads, 1, 1));
  } else if (op == 3) {
    const bool alternate = (localStep & 1) != 0;
    c.set_compute_pipeline_state(a.statePipeline);
    c.set_kernel_buffer(a.stateOffsets, 0);
    c.set_kernel_buffer(a.stateColumns, 1);
    c.set_kernel_buffer(a.stateValues, 2);
    c.set_kernel_buffer(a.stateValuesLow, 3);
    c.set_kernel_buffer(a.inputOffsets, 4);
    c.set_kernel_buffer(a.inputColumns, 5);
    c.set_kernel_buffer(a.inputValues, 6);
    c.set_kernel_buffer(a.inputValuesLow, 7);
    c.set_kernel_buffer(a.terminalIndices, 8);
    c.set_kernel_buffer(a.solution, 9);
    c.set_kernel_buffer(alternate ? a.alternateState : a.state, 10);
    c.set_kernel_buffer(alternate ? a.alternateStateLow : a.stateLow, 11);
    c.set_kernel_buffer(a.groupParams, 12);
    c.set_kernel_buffer(alternate ? a.alternateNextState : a.nextState, 13);
    c.set_kernel_buffer(alternate ? a.alternateNextStateLow : a.nextStateLow, 14);
    c.concurrent_dispatch_threadgroups(uint3(b.stateGroups, 1, 1), uint3(128, 1, 1));
  } else {
    c.set_compute_pipeline_state(a.generatorPostPipeline);
    c.set_kernel_buffer(a.generators, 0);
    c.set_kernel_buffer(a.solution, 1);
    c.set_kernel_buffer(step, 2);
    c.set_kernel_buffer(a.traceValues, 3);
    c.set_kernel_buffer(a.solverStatus, 4);
    c.set_kernel_buffer(a.iterationValues, 5);
    c.concurrent_dispatch_threadgroups(uint3(b.generatorGroups, 1, 1), uint3(128, 1, 1));
  }
}
)metal";

struct CsrFloat {
  std::vector<uint32_t> offsets;
  std::vector<uint32_t> columns;
  std::vector<float> values;
};

id<MTLBuffer> makeBuffer(id<MTLDevice> device, const void *source,
                         std::size_t bytes) {
  if (bytes == 0)
    bytes = sizeof(uint32_t);
  id<MTLBuffer> result = [device newBufferWithLength:bytes
                                             options:MTLResourceStorageModeShared];
  if (!result)
    throw std::runtime_error("MetalResidentScheduler buffer allocation failed");
  if (source && bytes > 0)
    std::memcpy(result.contents, source, bytes);
  return result;
}

uint32_t findPosition(const CPS::SparseMatrixRow &matrix, UInt row, UInt column) {
  for (Eigen::Index k = matrix.outerIndexPtr()[row];
       k < matrix.outerIndexPtr()[row + 1]; ++k) {
    if (static_cast<UInt>(matrix.innerIndexPtr()[k]) == column)
      return static_cast<uint32_t>(k);
  }
  throw std::runtime_error("Metal resident matrix is missing a generator stamp entry");
}
} // namespace

struct MetalResidentScheduler::Impl {
  CPS::EMT::LinearSSNGroup::Ptr group;
  std::vector<std::shared_ptr<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>> generators;
  CPS::EMT::SimNode::List nodes;
  Real dt;
  UInt expectedSteps;
  CPS::String tracePath;
  Bool wholeRunFastPath;
  Bool useSparseInputGpuLu;
  Bool useGpuGeneratedIcb;
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLComputePipelineState> clearPipeline = nil;
  id<MTLComputePipelineState> tracePipeline = nil;
  id<MTLComputePipelineState> clearIcbPipeline = nil;
  id<MTLComputePipelineState> icbBuilderPipeline = nil;
  id<MTLArgumentEncoder> icbArgumentEncoder = nil;
  id<MTLBuffer> icbArguments = nil;
  id<MTLBuffer> icbStepInfos = nil;
  id<MTLIndirectCommandBuffer> generatedIcb = nil;
  id<MTLBuffer> traceIndices = nil;
  id<MTLBuffer> traceValues = nil;
  id<MTLBuffer> iterationValues = nil;
  std::array<double, 8> phaseGpuSeconds{};
  std::array<UInt, 8> phaseTimingSamples{};
  static constexpr UInt phaseFirstSampleStep = 10000;
  static constexpr UInt phaseSampleSteps = 8;
  static constexpr UInt phaseCount = 8;
  static constexpr const char *phaseNames[phaseCount] = {
      "clear_rhs", "group_history_rhs", "generator_pre",
      "zero_initial_guess", "bicgstab_solve", "trace_capture",
      "group_state_update", "generator_post"};
  std::unique_ptr<MetalSparseBiCGSTABAdapter> solver;
  std::unique_ptr<MetalSparseLUAdapter> sparseInputGpuLu;
  std::unique_ptr<MetalLinearSSNGroupAdapter> groupAdapter;
  std::unique_ptr<MetalVBRGeneratorAdapter> generatorAdapter;
  id<MTLCommandBuffer> command = nil;
  id<MTLCommandBuffer> lastCommitted = nil;
  std::vector<id<MTLCommandBuffer>> pendingCommands;
  UInt encodedSteps = 0;
  UInt committedSteps = 0;
  UInt lastWaitedSteps = 0;
  double completedGpuSeconds = 0.0;
  UInt completedCommandBuffers = 0;
  Int nextScheduledStep = 0;
  UInt order = 0;
  UInt fullOrder = 0;
  UInt traceWidth = 0;
  std::vector<uint32_t> csrOffsets;
  std::vector<uint32_t> csrColumns;
  std::vector<UInt> activeOriginalIndices;
  std::vector<UInt> originalToCompact;
  Matrix firstStepExpectedRhs;
  Matrix firstStepExpectedGroupRhs;
  Matrix firstStepExpectedGeneratorRhs;
  std::vector<UInt> groupTerminalIndices;
  Bool checkedFirstChunk = false;

  UInt compactIndex(UInt original) const {
    if (original >= originalToCompact.size() ||
        originalToCompact[original] == std::numeric_limits<UInt>::max())
      throw std::runtime_error("Metal resident matrix contains an index outside the retained terminal set");
    return originalToCompact[original];
  }

  CPS::SparseMatrixRow compactSparseMatrix(
      const CPS::SparseMatrixRow &fullMatrix) const {
    CPS::SparseMatrixRow compact(static_cast<Int>(order), static_cast<Int>(order));
    for (Eigen::Index row = 0; row < fullMatrix.outerSize(); ++row) {
      for (CPS::SparseMatrixRow::InnerIterator entry(fullMatrix, row); entry; ++entry) {
        if (entry.value() == 0.0)
          continue;
        const UInt compactRow = compactIndex(static_cast<UInt>(row));
        const UInt compactCol = compactIndex(static_cast<UInt>(entry.col()));
        compact.coeffRef(compactRow, compactCol) = entry.value();
      }
    }
    compact.makeCompressed();
    return compact;
  }

  Matrix compactVector(const Matrix &fullVector) const {
    if (fullVector.rows() != static_cast<Eigen::Index>(fullOrder))
      throw std::invalid_argument("Metal resident RHS has an unexpected full MNA order");
    Matrix compact(static_cast<Eigen::Index>(order), fullVector.cols());
    for (UInt i = 0; i < order; ++i)
      for (Eigen::Index col = 0; col < fullVector.cols(); ++col)
        compact(i, col) = fullVector(activeOriginalIndices[i], col);
    return compact;
  }

  Matrix expandVector(const Matrix &compactValue) const {
    if (compactValue.rows() != static_cast<Eigen::Index>(order))
      throw std::invalid_argument("Metal resident solution has an unexpected compact order");
    Matrix full = Matrix::Zero(fullOrder, compactValue.cols());
    for (UInt i = 0; i < order; ++i)
      for (Eigen::Index col = 0; col < compactValue.cols(); ++col)
        full(activeOriginalIndices[i], col) = compactValue(i, col);
    return full;
  }

  Impl(CPS::EMT::LinearSSNGroup::Ptr linearGroup,
       std::vector<std::shared_ptr<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>> machines,
       CPS::EMT::SimNode::List simNodes, Real step, UInt steps, CPS::String path,
       Bool runWholeInterval, Bool useGpuLu, Bool useIcb)
      : group(std::move(linearGroup)), generators(std::move(machines)),
        nodes(std::move(simNodes)), dt(step), expectedSteps(steps),
        tracePath(std::move(path)), wholeRunFastPath(runWholeInterval),
        useSparseInputGpuLu(useGpuLu), useGpuGeneratedIcb(useIcb) {}

  void initialize(const CPS::Task::List &tasks) {
    if (!group || generators.empty() || nodes.empty() || dt <= 0.0 ||
        expectedSteps == 0)
      throw std::invalid_argument("Metal resident scheduler received invalid model settings");
    (void)tasks;

    UInt maxIndex = 0;
    for (const auto &node : nodes)
      for (UInt index : node->matrixNodeIndices())
        maxIndex = std::max(maxIndex, index);
    for (const auto &gen : generators)
      for (UInt index : gen->metalDeviceData().matrixNodeIndices)
        maxIndex = std::max(maxIndex, index);
    fullOrder = maxIndex + 1;
    const auto stateMaps = group->precomputedStateMaps();
    if (stateMaps.terminalMatrixIndices.empty())
      throw std::runtime_error("Metal resident path requires a fully initialized linear SSN group");
    activeOriginalIndices = stateMaps.terminalMatrixIndices;
    order = static_cast<UInt>(activeOriginalIndices.size());
    if (fullOrder == 0 || order == 0 || order > 512)
      throw std::runtime_error("Metal resident IEEE118 solver order is unsupported");
    originalToCompact.assign(fullOrder, std::numeric_limits<UInt>::max());
    for (UInt local = 0; local < order; ++local) {
      const UInt original = activeOriginalIndices[local];
      if (original >= fullOrder ||
          originalToCompact[original] != std::numeric_limits<UInt>::max())
        throw std::runtime_error("Metal resident retained terminal indices are invalid or duplicated");
      originalToCompact[original] = local;
    }

    CPS::SparseMatrixRow matrix(static_cast<Int>(fullOrder), static_cast<Int>(fullOrder));
    std::vector<std::array<UInt, 3>> genIndices;
    genIndices.reserve(generators.size());
    for (const auto &gen : generators) {
      auto data = gen->metalDeviceData();
      genIndices.push_back(data.matrixNodeIndices);
      for (UInt index : data.matrixNodeIndices)
        (void)compactIndex(index);
    }
    CPS::SparseMatrixRow groupStamp(static_cast<Int>(fullOrder), static_cast<Int>(fullOrder));
    group->mnaApplySystemMatrixStamp(groupStamp);
    matrix = std::move(groupStamp);

    // SG4 stamps are updated in place by the device kernel. Create all nine
    // positions explicitly in the fixed matrix before compression so cached
    // CSR offsets address these slots and nothing from a CPU-side initial
    // generator stamp can leak into the time-varying coefficients.
    for (const auto &indices : genIndices)
      for (UInt row = 0; row < 3; ++row)
        for (UInt col = 0; col < 3; ++col)
          matrix.coeffRef(indices[row], indices[col]) += 0.0;

    // Position caching is valid only after the CSR/row-major structure is
    // finalized; makeCompressed() may reorder entries otherwise.
    matrix.makeCompressed();
    matrix = compactSparseMatrix(matrix);
    for (auto &indices : genIndices)
      for (UInt &index : indices)
        index = compactIndex(index);
    for (const auto &indices : genIndices)
      for (UInt row = 0; row < 3; ++row)
        for (UInt col = 0; col < 3; ++col)
          matrix.coeffRef(indices[row], indices[col]) += 0.0;
    matrix.makeCompressed();
    const double matrixDensity = static_cast<double>(matrix.nonZeros()) /
        (static_cast<double>(order) * static_cast<double>(order));
    const Bool useDenseMatrixStorage = matrixDensity >= 0.70;
    if (useDenseMatrixStorage) {
      for (UInt row = 0; row < order; ++row)
        for (UInt col = 0; col < order; ++col)
          matrix.coeffRef(row, col) += 0.0;
      matrix.makeCompressed();
    }
    std::cout << "Metal resident reduced MNA matrix: " << order << " x "
              << order << ", initial sparsity=" << matrixDensity
              << ", storage=" << (useDenseMatrixStorage ? "dense" : "CSR")
              << '\n';
    csrOffsets.resize(static_cast<std::size_t>(order) + 1);
    csrColumns.resize(static_cast<std::size_t>(matrix.nonZeros()));
    for (UInt row = 0; row <= order; ++row)
      csrOffsets[row] = static_cast<uint32_t>(matrix.outerIndexPtr()[row]);
    for (Eigen::Index k = 0; k < matrix.nonZeros(); ++k)
      csrColumns[static_cast<std::size_t>(k)] =
          static_cast<uint32_t>(matrix.innerIndexPtr()[k]);
    std::vector<MetalVBRGeneratorAdapter::StampPositions> positions;
    positions.reserve(generators.size());
    for (std::size_t g = 0; g < generators.size(); ++g) {
      MetalVBRGeneratorAdapter::StampPositions one{};
      for (UInt row = 0; row < 3; ++row)
        for (UInt col = 0; col < 3; ++col) {
          const UInt slot = row * 3 + col;
          one.positions[slot] = findPosition(
              matrix, genIndices[g][row], genIndices[g][col]);
          one.fixedValues[slot] = static_cast<float>(
              matrix.valuePtr()[one.positions[slot]]);
        }
      positions.push_back(one);
    }
    solver = std::make_unique<MetalSparseBiCGSTABAdapter>();
    solver->setConvergenceCriteria(512, 1e-7f, 1e-8f);
    std::vector<std::pair<UInt, UInt>> variableEntries;
    solver->preprocessing(matrix, variableEntries);
    solver->factorize(matrix);

    auto compactStateMaps = stateMaps;
    for (UInt &index : compactStateMaps.terminalMatrixIndices)
      index = compactIndex(index);
    groupAdapter = std::make_unique<MetalLinearSSNGroupAdapter>(compactStateMaps);
    groupTerminalIndices = compactStateMaps.terminalMatrixIndices;
    std::vector<CPS::EMT::Ph3::SynchronGenerator4OrderVBR::MetalDeviceData> data;
    data.reserve(generators.size());
    for (const auto &gen : generators) {
      auto deviceData = gen->metalDeviceData();
      for (UInt &index : deviceData.matrixNodeIndices)
        index = compactIndex(index);
      data.push_back(deviceData);
    }
    generatorAdapter = std::make_unique<MetalVBRGeneratorAdapter>(data, positions);

    // One-time validation oracle from the model's own pre-step and compact
    // RHS stamp implementations. Device snapshots above retain the exact
    // pre-step state; these CPU components are not scheduled by this executor.
    Matrix fullExpectedRhs = Matrix::Zero(fullOrder, 1);
    Matrix fullExpectedGroupRhs = Matrix::Zero(fullOrder, 1);
    Matrix fullExpectedGeneratorRhs = Matrix::Zero(fullOrder, 1);
    group->mnaPreStep(dt, 0);
    group->mnaAddCompactRightSideVectorStamp(fullExpectedGroupRhs);
    fullExpectedRhs += fullExpectedGroupRhs;
    for (const auto &gen : generators) {
      gen->mnaPreStep(dt, 0);
      gen->mnaApplyRightSideVectorStamp(fullExpectedGeneratorRhs);
    }
    fullExpectedRhs += fullExpectedGeneratorRhs;
    firstStepExpectedRhs = compactVector(fullExpectedRhs);
    firstStepExpectedGroupRhs = compactVector(fullExpectedGroupRhs);
    firstStepExpectedGeneratorRhs = compactVector(fullExpectedGeneratorRhs);

    device = MTLCreateSystemDefaultDevice();
    if (!device)
      throw std::runtime_error("Metal resident path requested without an Apple Metal device");
    queue = [device newCommandQueue];
    if (!queue)
      throw std::runtime_error("Metal resident scheduler could not create a command queue");
    if (useSparseInputGpuLu)
      sparseInputGpuLu = std::make_unique<MetalSparseLUAdapter>(
          (void *)device, order, csrOffsets, csrColumns);
    MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
    NSError *error = nil;
    id<MTLLibrary> library = [device newLibraryWithSource:
        [NSString stringWithUTF8String:kResidentKernels] options:options error:&error];
    [options release];
    if (!library)
      throw std::runtime_error(std::string("Metal resident shader compilation failed: ") +
                               (error.localizedDescription.UTF8String ?: "unknown error"));
    id<MTLFunction> function = [library newFunctionWithName:@"clearRhs"];
    clearPipeline = [device newComputePipelineStateWithFunction:function error:&error];
    [function release];
    function = [library newFunctionWithName:@"captureTrace"];
    tracePipeline = [device newComputePipelineStateWithFunction:function error:&error];
    [function release];
    [library release];
    if (!clearPipeline || !tracePipeline)
      throw std::runtime_error("Metal resident scheduler failed to create its trace/clear kernels");

    traceWidth = static_cast<UInt>(nodes.size() * 3);
    std::vector<uint32_t> indices(traceWidth);
    for (std::size_t b = 0; b < nodes.size(); ++b) {
      const auto &node = nodes[b];
      const auto nodeIndices = node->matrixNodeIndices();
      if (nodeIndices.size() < 3)
        throw std::runtime_error("Metal trace node is not a three-phase node");
      for (std::size_t phase = 0; phase < 3; ++phase)
        indices[b * 3 + phase] = compactIndex(nodeIndices[phase]);
    }
    traceIndices = makeBuffer(device, indices.data(),
                              indices.size() * sizeof(uint32_t));
    traceValues = makeBuffer(device, nullptr,
                             static_cast<std::size_t>(expectedSteps) *
                                 indices.size() * sizeof(float));
    iterationValues = makeBuffer(device, nullptr,
                                 static_cast<std::size_t>(expectedSteps) *
                                     sizeof(uint32_t));
    if (useGpuGeneratedIcb) {
      if (useSparseInputGpuLu)
        throw std::invalid_argument(
            "GPU-generated timestep ICB currently requires the Metal BiCGSTAB path");
      initializeGpuGeneratedIcb();
    }
  }

  template <typename EncodePhase>
  void runTimedPhase(UInt phaseIndex, EncodePhase &&encodePhase) {
    id<MTLCommandBuffer> phaseCommand = [queue commandBuffer];
    if (!phaseCommand)
      throw std::runtime_error("Could not allocate phase-timing command buffer");
    encodePhase(phaseCommand);
    [phaseCommand commit];
    [phaseCommand waitUntilCompleted];
    if (phaseCommand.error)
      throw std::runtime_error("GPU phase-timing command failed: " +
          std::string(phaseCommand.error.localizedDescription.UTF8String ?: "unknown error"));
    if (phaseCommand.GPUEndTime >= phaseCommand.GPUStartTime &&
        phaseCommand.GPUStartTime > 0.0) {
      phaseGpuSeconds[phaseIndex] +=
          phaseCommand.GPUEndTime - phaseCommand.GPUStartTime;
      ++phaseTimingSamples[phaseIndex];
      completedGpuSeconds += phaseCommand.GPUEndTime - phaseCommand.GPUStartTime;
      ++completedCommandBuffers;
    }
  }

  void reportPhaseTimings() const {
    if (std::accumulate(phaseTimingSamples.begin(), phaseTimingSamples.end(), 0u) == 0)
      return;
    std::cout << "GPU per-phase times (separate command buffers, "
              << phaseSampleSteps << " mid-run timesteps; diagnostic only):\n";
    for (UInt phase = 0; phase < phaseCount; ++phase) {
      std::cout << "  " << phaseNames[phase] << ": ";
      if (phaseTimingSamples[phase])
        std::cout << phaseGpuSeconds[phase] * 1.0e6 /
                        phaseTimingSamples[phase]
                  << " us/step (" << phaseTimingSamples[phase] << " samples)\n";
      else
        std::cout << "one-time initialization / not sampled\n";
    }
  }

  id<MTLComputePipelineState> makeIcbPipeline(id<MTLLibrary> library,
                                               NSString *functionName) {
    id<MTLFunction> function = [library newFunctionWithName:functionName];
    if (!function)
      throw std::runtime_error(std::string("Missing GPU ICB kernel ") +
                               functionName.UTF8String);
    MTLComputePipelineDescriptor *descriptor =
        [[MTLComputePipelineDescriptor alloc] init];
    descriptor.computeFunction = function;
    descriptor.supportIndirectCommandBuffers = YES;
    NSError *error = nil;
    id<MTLComputePipelineState> result =
        [device newComputePipelineStateWithDescriptor:descriptor
                                              options:0
                                           reflection:nil
                                                error:&error];
    [descriptor release];
    [function release];
    if (!result)
      throw std::runtime_error(std::string("Failed to create ICB pipeline ") +
                               functionName.UTF8String + ": " +
                               (error.localizedDescription.UTF8String ?: "unknown error"));
    return result;
  }

  void initializeGpuGeneratedIcb() {
    NSError *error = nil;
    MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
    id<MTLLibrary> library = [device
        newLibraryWithSource:[NSString stringWithUTF8String:kResidentIcbKernels]
                     options:options error:&error];
    [options release];
    if (!library)
      throw std::runtime_error(std::string("GPU timestep ICB shader compilation failed: ") +
                               (error.localizedDescription.UTF8String ?: "unknown error"));
    clearIcbPipeline = makeIcbPipeline(library, @"clearRhsForIcb");
    icbBuilderPipeline = makeIcbPipeline(library, @"buildTimestepCommands");
    id<MTLFunction> builder = [library newFunctionWithName:@"buildTimestepCommands"];
    icbArgumentEncoder = [builder newArgumentEncoderWithBufferIndex:0];
    [builder release];
    [library release];
    if (!icbArgumentEncoder)
      throw std::runtime_error("GPU timestep ICB argument encoder creation failed");
    icbArguments = [device newBufferWithLength:icbArgumentEncoder.encodedLength
                                       options:MTLResourceStorageModeShared];
    if (!icbArguments)
      throw std::runtime_error("GPU timestep ICB argument buffer allocation failed");

    auto group = groupAdapter->indirectResources();
    auto gen = generatorAdapter->indirectResources();
    auto solve = solver->indirectResources();
    [icbArgumentEncoder setArgumentBuffer:icbArguments offset:0];
    [icbArgumentEncoder setComputePipelineState:clearIcbPipeline atIndex:0];
    [icbArgumentEncoder setComputePipelineState:(id<MTLComputePipelineState>)group.historyPipeline atIndex:1];
    [icbArgumentEncoder setComputePipelineState:(id<MTLComputePipelineState>)gen.prePipeline atIndex:3];
    [icbArgumentEncoder setComputePipelineState:(id<MTLComputePipelineState>)solve.pipeline atIndex:5];
    [icbArgumentEncoder setComputePipelineState:(id<MTLComputePipelineState>)group.updatePipeline atIndex:7];
    [icbArgumentEncoder setComputePipelineState:(id<MTLComputePipelineState>)gen.postPipeline atIndex:8];

    // IDs 10..45 mirror IcbArguments in kResidentIcbKernels.
    id<MTLBuffer> args[] = {
        (id<MTLBuffer>)solver->rightSideBuffer(),
        (id<MTLBuffer>)solver->rowOffsetsBuffer(),
        (id<MTLBuffer>)solver->columnIndicesBuffer(),
        (id<MTLBuffer>)solver->matrixValuesBuffer(),
        (id<MTLBuffer>)solver->solutionBuffer(),
        (id<MTLBuffer>)solve.work,
        (id<MTLBuffer>)solve.status,
        (id<MTLBuffer>)solve.params,
        (id<MTLBuffer>)group.historyOffsets,
        (id<MTLBuffer>)group.historyColumns,
        (id<MTLBuffer>)group.historyValues,
        (id<MTLBuffer>)group.historyValuesLow,
        (id<MTLBuffer>)group.stateOffsets,
        (id<MTLBuffer>)group.stateColumns,
        (id<MTLBuffer>)group.stateValues,
        (id<MTLBuffer>)group.stateValuesLow,
        (id<MTLBuffer>)group.inputOffsets,
        (id<MTLBuffer>)group.inputColumns,
        (id<MTLBuffer>)group.inputValues,
        (id<MTLBuffer>)group.inputValuesLow,
        (id<MTLBuffer>)group.terminalIndices,
        (id<MTLBuffer>)group.state,
        (id<MTLBuffer>)group.stateLow,
        (id<MTLBuffer>)group.nextState,
        (id<MTLBuffer>)group.nextStateLow,
        (id<MTLBuffer>)group.nextState,
        (id<MTLBuffer>)group.nextStateLow,
        (id<MTLBuffer>)group.state,
        (id<MTLBuffer>)group.stateLow,
        (id<MTLBuffer>)group.history,
        (id<MTLBuffer>)group.params,
        (id<MTLBuffer>)gen.generators,
        nil, // stepInfosWrite, allocated for each run
        nil, // stepInfosRead, same buffer with read-only shader access
        traceIndices,
        traceValues,
        iterationValues};
    for (NSUInteger i = 0; i < sizeof(args) / sizeof(args[0]); ++i)
      if (args[i])
        [icbArgumentEncoder setBuffer:args[i] offset:0 atIndex:i + 10];
  }

  id<MTLCommandBuffer> beginChunkIfNeeded() {
    if (command)
      return command;
    command = [queue commandBuffer];
    if (!command)
      throw std::runtime_error("Metal resident scheduler failed to allocate command buffer");
    encodedSteps = 0;
    return command;
  }

  void waitForPendingCommands() {
    std::string commandError;
    for (auto pending : pendingCommands) {
      [pending waitUntilCompleted];
      if (pending.GPUEndTime >= pending.GPUStartTime &&
          pending.GPUStartTime > 0.0) {
        completedGpuSeconds += pending.GPUEndTime - pending.GPUStartTime;
        ++completedCommandBuffers;
      }
      if (commandError.empty() && pending.error)
        commandError = pending.error.localizedDescription.UTF8String
            ? pending.error.localizedDescription.UTF8String
            : "unknown Metal command-buffer error";
      [pending release];
    }
    pendingCommands.clear();
    lastWaitedSteps = committedSteps;
    if (!commandError.empty())
      throw std::runtime_error("Metal resident command chunk failed: " + commandError);
  }

  void encodeSingleStepWithPhaseTiming(Real time, Int count) {
    // Isolate representative mid-run steps into separate command buffers.
    // This works on devices where dispatch-boundary Metal counters are absent.
    if (command) {
      [command commit];
      pendingCommands.push_back([command retain]);
      committedSteps += encodedSteps;
      command = nil;
    }
    if (!pendingCommands.empty())
      waitForPendingCommands();

    runTimedPhase(1, [&](id<MTLCommandBuffer> cmd) {
      groupAdapter->encodeHistoryAndRhs((void *)cmd, solver->rightSideBuffer());
    });
    runTimedPhase(2, [&](id<MTLCommandBuffer> cmd) {
      generatorAdapter->encodePreStep((void *)cmd, solver->matrixValuesBuffer(),
                                      solver->rightSideBuffer(), time);
    });
    runTimedPhase(4, [&](id<MTLCommandBuffer> cmd) {
      if (sparseInputGpuLu)
        sparseInputGpuLu->encodeSolve((void *)cmd, solver->matrixValuesBuffer(),
                                      solver->rightSideBuffer(), solver->solutionBuffer());
      else
        solver->encodeSolve((void *)cmd, solver->matrixValuesBuffer(),
                            solver->rightSideBuffer(), solver->solutionBuffer());
    });
    runTimedPhase(5, [&](id<MTLCommandBuffer> cmd) {
      id<MTLComputeCommandEncoder> encoder = [cmd computeCommandEncoder];
      const uint32_t offset = static_cast<uint32_t>(count) * traceWidth;
      [encoder setComputePipelineState:tracePipeline];
      [encoder setBuffer:(id<MTLBuffer>)solver->solutionBuffer() offset:0 atIndex:0];
      [encoder setBuffer:traceIndices offset:0 atIndex:1];
      [encoder setBuffer:traceValues offset:0 atIndex:2];
      [encoder setBytes:&offset length:sizeof(offset) atIndex:3];
      [encoder setBuffer:(id<MTLBuffer>)solver->statusBuffer() offset:0 atIndex:4];
      [encoder setBuffer:iterationValues offset:0 atIndex:5];
      const uint32_t step = static_cast<uint32_t>(count);
      [encoder setBytes:&step length:sizeof(step) atIndex:6];
      [encoder dispatchThreads:MTLSizeMake(traceWidth, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min<NSUInteger>(traceWidth, tracePipeline.maxTotalThreadsPerThreadgroup), 1, 1)];
      [encoder endEncoding];
    });
    runTimedPhase(6, [&](id<MTLCommandBuffer> cmd) {
      groupAdapter->encodeStateUpdate((void *)cmd, solver->solutionBuffer());
    });
    runTimedPhase(7, [&](id<MTLCommandBuffer> cmd) {
      generatorAdapter->encodePostStep((void *)cmd, solver->solutionBuffer());
    });
    committedSteps = static_cast<UInt>(count + 1);
    lastWaitedSteps = committedSteps;
  }

  void encodeSingleStep(Real time, Int count) {
    if (static_cast<UInt>(count) >= expectedSteps)
      throw std::runtime_error("Metal resident trace buffer is smaller than the requested run");
    if (!wholeRunFastPath && !useSparseInputGpuLu &&
        count >= static_cast<Int>(phaseFirstSampleStep) &&
        count < static_cast<Int>(phaseFirstSampleStep + phaseSampleSteps)) {
      encodeSingleStepWithPhaseTiming(time, count);
      return;
    }
    auto cmd = beginChunkIfNeeded();
    groupAdapter->encodeHistoryAndRhs((void *)cmd,
                                      solver->rightSideBuffer());
    generatorAdapter->encodePreStep((void *)cmd, solver->matrixValuesBuffer(),
                                    solver->rightSideBuffer(), time);
    if (sparseInputGpuLu)
      sparseInputGpuLu->encodeSolve((void *)cmd, solver->matrixValuesBuffer(),
                                    solver->rightSideBuffer(), solver->solutionBuffer());
    else
      solver->encodeSolve((void *)cmd, solver->matrixValuesBuffer(),
                          solver->rightSideBuffer(), solver->solutionBuffer());
    id<MTLComputeCommandEncoder> encoder = [cmd computeCommandEncoder];
    const uint32_t offset = static_cast<uint32_t>(count) * traceWidth;
    [encoder setComputePipelineState:tracePipeline];
    [encoder setBuffer:(id<MTLBuffer>)solver->solutionBuffer() offset:0 atIndex:0];
    [encoder setBuffer:traceIndices offset:0 atIndex:1];
    [encoder setBuffer:traceValues offset:0 atIndex:2];
    [encoder setBytes:&offset length:sizeof(offset) atIndex:3];
    [encoder setBuffer:(id<MTLBuffer>)solver->statusBuffer() offset:0 atIndex:4];
    [encoder setBuffer:iterationValues offset:0 atIndex:5];
    const uint32_t step = static_cast<uint32_t>(count);
    [encoder setBytes:&step length:sizeof(step) atIndex:6];
    [encoder dispatchThreads:MTLSizeMake(traceWidth, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min<NSUInteger>(traceWidth, tracePipeline.maxTotalThreadsPerThreadgroup), 1, 1)];
    [encoder endEncoding];
    groupAdapter->encodeStateUpdate((void *)cmd, solver->solutionBuffer());
    generatorAdapter->encodePostStep((void *)cmd, solver->solutionBuffer());
    ++encodedSteps;
    if (encodedSteps == 64 || count <= 1 ||
        static_cast<UInt>(count + 1) == expectedSteps) {
      [cmd commit];
      [lastCommitted release];
      lastCommitted = [cmd retain];
      pendingCommands.push_back([cmd retain]);
      command = nil;
      committedSteps += encodedSteps;
      // Keep a bounded number of timesteps in flight. All batch command
      // buffers share one ordered queue, so intermediate barriers/readbacks
      // are unnecessary; status and state remain device-resident.
      const Bool shouldWait = count <= 1 ||
          committedSteps - lastWaitedSteps >= 1024 ||
          committedSteps == expectedSteps;
      if (shouldWait)
        waitForPendingCommands();
      if (!shouldWait)
        return;
      const auto *status = static_cast<const uint32_t *>(
          [(id<MTLBuffer>)solver->statusBuffer() contents]);
      if (sparseInputGpuLu) {
        sparseInputGpuLu->checkStatus();
      } else if (status[0] != 0) {
          float residual = 0.0f, rhsNorm = 0.0f;
          std::memcpy(&residual, &status[2], sizeof(float));
          std::memcpy(&rhsNorm, &status[3], sizeof(float));
          throw std::runtime_error("Metal resident BiCGSTAB failed at timestep " +
              std::to_string(count) + " (status=" + std::to_string(status[0]) +
              ", iterations=" + std::to_string(status[1]) +
              ", relative residual=" +
              std::to_string(residual / std::max(rhsNorm, 1e-30f)) + ")");
      }
      if (!checkedFirstChunk) {
        const auto *values = static_cast<const float *>(
            [(id<MTLBuffer>)solver->matrixValuesBuffer() contents]);
        const auto *rhs = static_cast<const float *>(
            [(id<MTLBuffer>)solver->rightSideBuffer() contents]);
        const auto *gpuSolution = static_cast<const float *>(
            [(id<MTLBuffer>)solver->solutionBuffer() contents]);
        Matrix dense = Matrix::Zero(order, order);
        Matrix hostRhs(order, 1), hostGpuSolution(order, 1);
        for (UInt row = 0; row < order; ++row) {
          hostRhs(row) = rhs[row];
          hostGpuSolution(row) = gpuSolution[row];
          for (UInt k = csrOffsets[row]; k < csrOffsets[row + 1]; ++k)
            dense(row, csrColumns[k]) = values[k];
        }
        const Matrix hostReference = dense.partialPivLu().solve(hostRhs);
        const Real relativeDifference =
            (hostGpuSolution - hostReference).norm() /
            std::max<Real>(hostReference.norm(), 1e-30);
        const Real relativeResidual =
            (dense * hostGpuSolution - hostRhs).norm() /
            std::max<Real>(hostRhs.norm(), 1e-30);
        const Real rhsDifference =
            (hostRhs - firstStepExpectedRhs).norm() /
            std::max<Real>(firstStepExpectedRhs.norm(), 1e-30);
        Matrix gpuGroupRhs = Matrix::Zero(order, 1);
        const auto *gpuHistory = static_cast<const float *>(
            [(id<MTLBuffer>)groupAdapter->historyBuffer() contents]);
        for (UInt terminal = 0; terminal < groupTerminalIndices.size(); ++terminal)
          gpuGroupRhs(groupTerminalIndices[terminal]) += gpuHistory[terminal];
        const Matrix gpuGeneratorRhs = hostRhs - gpuGroupRhs;
        const Real groupRhsDifference =
            (gpuGroupRhs - firstStepExpectedGroupRhs).norm() /
            std::max<Real>(firstStepExpectedGroupRhs.norm(), 1e-30);
        const Real generatorRhsDifference =
            (gpuGeneratorRhs - firstStepExpectedGeneratorRhs).norm() /
            std::max<Real>(firstStepExpectedGeneratorRhs.norm(), 1e-30);
        float kernelResidual = 0.0f, kernelRhsNorm = 0.0f;
        std::memcpy(&kernelResidual, &status[2], sizeof(float));
        std::memcpy(&kernelRhsNorm, &status[3], sizeof(float));
        std::cout << "Metal resident first-chunk CPU cross-check: relative "
                     "solution difference=" << relativeDifference
                  << ", CPU residual=" << relativeResidual
                  << ", CPU RHS difference=" << rhsDifference
                  << ", group RHS difference=" << groupRhsDifference
                  << ", generator RHS difference=" << generatorRhsDifference
                  << (sparseInputGpuLu ? ", direct GPU LU solve"
                                       : ", kernel residual=" +
                                             std::to_string(kernelResidual /
                                                 std::max(kernelRhsNorm, 1e-30f)) +
                                             ", iterations=" +
                                                 std::to_string(status[1]))
                  << ", status words=" << status[0] << ':' << status[1]
                  << ':' << status[2] << ':' << status[3]
                  << ", host RHS norm=" << hostRhs.norm()
                  << ", GPU solution norm=" << hostGpuSolution.norm()
                  << '\n';
        CPS::Attribute<Matrix>::Ptr referenceLeft =
            CPS::AttributeStatic<Matrix>::make();
        **referenceLeft = expandVector(hostGpuSolution);
        group->mnaPostStep(dt, 0, referenceLeft);
        for (const auto &gen : generators)
          gen->mnaPostStep(dt, 0, referenceLeft);
        const Matrix cpuGroupState = group->precomputedStateMaps().initialState;
        const auto gpuGroupState = groupAdapter->copyStateToHost();
        Matrix gpuGroupStateMatrix(cpuGroupState.rows(), 1);
        for (Eigen::Index i = 0; i < cpuGroupState.rows(); ++i)
          gpuGroupStateMatrix(i) = gpuGroupState[static_cast<std::size_t>(i)];
        const Real groupStateDifference =
            (gpuGroupStateMatrix - cpuGroupState).norm() /
            std::max<Real>(cpuGroupState.norm(), 1e-30);
        const auto gpuGeneratorState = generatorAdapter->copyStateToHost();
        Real generatorStateDifference = 0.0;
        Real generatorStateReferenceNorm = 0.0;
        for (std::size_t i = 0; i < generators.size(); ++i) {
          const auto cpu = generators[i]->metalDeviceData();
          const auto &gpu = gpuGeneratorState[i];
          for (UInt phase = 0; phase < 3; ++phase) {
            const Real dv = gpu.voltageDq0[phase] - cpu.voltageDq0[phase];
            const Real di = gpu.currentDq0[phase] - cpu.currentDq0[phase];
            generatorStateDifference += dv * dv + di * di;
            generatorStateReferenceNorm += cpu.voltageDq0[phase] * cpu.voltageDq0[phase] +
                                           cpu.currentDq0[phase] * cpu.currentDq0[phase];
          }
          const Real ds = gpu.rotorSpeed - cpu.rotorSpeed;
          const Real da = gpu.rotorAngle - cpu.rotorAngle;
          generatorStateDifference += ds * ds + da * da;
          generatorStateReferenceNorm += cpu.rotorSpeed * cpu.rotorSpeed +
                                         cpu.rotorAngle * cpu.rotorAngle;
        }
        generatorStateDifference = std::sqrt(generatorStateDifference /
            std::max<Real>(generatorStateReferenceNorm, 1e-30));
        std::cout << "Metal resident first-step component-state cross-check: "
                     "group relative difference=" << groupStateDifference
                  << ", generator relative difference=" << generatorStateDifference
                  << '\n';
        checkedFirstChunk = true;
        if (!std::isfinite(relativeDifference) ||
            !std::isfinite(relativeResidual) || relativeDifference > 1e-3 ||
            relativeResidual > 1e-3 || !std::isfinite(rhsDifference) ||
            rhsDifference > 1e-3 || groupRhsDifference > 1e-3 ||
            generatorRhsDifference > 1e-3)
          throw std::runtime_error(
              "Metal resident first-step solve failed the independent CPU check");
      } else if (count == 1) {
        Matrix fullCpuRhs = Matrix::Zero(fullOrder, 1);
        group->mnaPreStep(time, count);
        group->mnaAddCompactRightSideVectorStamp(fullCpuRhs);
        for (const auto &gen : generators) {
          gen->mnaPreStep(time, count);
          const auto data = gen->metalDeviceData();
          const auto componentRhs = gen->getRightVector();
          for (UInt phase = 0; phase < 3; ++phase)
            fullCpuRhs(data.matrixNodeIndices[phase]) +=
                (**componentRhs)(data.matrixNodeIndices[phase], 0);
        }
        CPS::SparseMatrixRow fullCpuMatrix(static_cast<Int>(fullOrder),
                                           static_cast<Int>(fullOrder));
        CPS::SparseMatrixRow fullCpuGroupMatrix(static_cast<Int>(fullOrder),
                                                static_cast<Int>(fullOrder));
        group->mnaApplySystemMatrixStamp(fullCpuMatrix);
        group->mnaApplySystemMatrixStamp(fullCpuGroupMatrix);
        for (const auto &gen : generators)
          gen->mnaApplySystemMatrixStamp(fullCpuMatrix);
        fullCpuMatrix.makeCompressed();
        fullCpuGroupMatrix.makeCompressed();
        CPS::SparseMatrixRow cpuMatrix = compactSparseMatrix(fullCpuMatrix);
        CPS::SparseMatrixRow cpuGroupMatrix = compactSparseMatrix(fullCpuGroupMatrix);
        Matrix cpuRhs = compactVector(fullCpuRhs);
        const auto *values = static_cast<const float *>(
            [(id<MTLBuffer>)solver->matrixValuesBuffer() contents]);
        const auto *gpuRhs = static_cast<const float *>(
            [(id<MTLBuffer>)solver->rightSideBuffer() contents]);
        const auto *gpuSolution = static_cast<const float *>(
            [(id<MTLBuffer>)solver->solutionBuffer() contents]);
        Matrix dense = Matrix::Zero(order, order);
        Matrix hostRhs(order, 1), hostSolution(order, 1);
        Real matrixError2 = 0.0, matrixNorm2 = 0.0;
        Real fixedMatrixError2 = 0.0, fixedMatrixNorm2 = 0.0;
        Real maxGeneratorMatrixDelta = 0.0;
        UInt maxMatrixRow = 0, maxMatrixCol = 0;
        std::vector<uint8_t> generatorPosition(values == nullptr ? 0 : csrColumns.size(), 0);
        for (const auto &gen : generators) {
          const auto originalIndices = gen->metalDeviceData().matrixNodeIndices;
          std::array<UInt, 3> indices{};
          for (UInt phase = 0; phase < 3; ++phase)
            indices[phase] = compactIndex(originalIndices[phase]);
          for (UInt localRow = 0; localRow < 3; ++localRow)
            for (UInt localCol = 0; localCol < 3; ++localCol) {
              const UInt row = indices[localRow];
              const UInt col = indices[localCol];
              for (UInt k = csrOffsets[row]; k < csrOffsets[row + 1]; ++k)
                if (csrColumns[k] == col)
                  generatorPosition[k] = 1;
            }
        }
        for (UInt row = 0; row < order; ++row) {
          hostRhs(row) = gpuRhs[row];
          hostSolution(row) = gpuSolution[row];
          for (UInt k = csrOffsets[row]; k < csrOffsets[row + 1]; ++k) {
            const UInt col = csrColumns[k];
            dense(row, col) = values[k];
            const Real expected = cpuMatrix.coeff(row, col);
            const Real delta = values[k] - expected;
            matrixError2 += delta * delta;
            matrixNorm2 += expected * expected;
            if (std::abs(delta) > maxGeneratorMatrixDelta) {
              maxGeneratorMatrixDelta = std::abs(delta);
              maxMatrixRow = row;
              maxMatrixCol = col;
            }
            if (!generatorPosition[k]) {
              const Real fixedExpected = cpuGroupMatrix.coeff(row, col);
              const Real fixedDelta = values[k] - fixedExpected;
              fixedMatrixError2 += fixedDelta * fixedDelta;
              fixedMatrixNorm2 += fixedExpected * fixedExpected;
            }
          }
        }
        const Real rhsError = (hostRhs - cpuRhs).norm() /
                              std::max<Real>(cpuRhs.norm(), 1e-30);
        const Real matrixError = std::sqrt(matrixError2 /
                                    std::max<Real>(matrixNorm2, 1e-30));
        const Real fixedMatrixError = std::sqrt(fixedMatrixError2 /
                                    std::max<Real>(fixedMatrixNorm2, 1e-30));
        Matrix denseCpuMatrix = Matrix::Zero(order, order);
        for (Eigen::Index row = 0; row < cpuMatrix.outerSize(); ++row)
          for (CPS::SparseMatrixRow::InnerIterator it(cpuMatrix, row); it; ++it)
            denseCpuMatrix(it.row(), it.col()) = it.value();
        const Matrix trueCpuSolution = denseCpuMatrix.partialPivLu().solve(cpuRhs);
        const Real trueSolutionError =
            (hostSolution - trueCpuSolution).norm() /
            std::max<Real>(trueCpuSolution.norm(), 1e-30);
        const auto gpuGeneratorState = generatorAdapter->copyStateToHost();
        Real generatorMatrixError2 = 0.0, generatorMatrixNorm2 = 0.0;
        Real generatorRhsError2 = 0.0, generatorRhsNorm2 = 0.0;
        Real globalGeneratorStampError2 = 0.0;
        Real globalGeneratorStampNorm2 = 0.0;
        for (std::size_t g = 0; g < generators.size(); ++g) {
          const auto cpuData = generators[g]->metalDeviceData();
          const auto &gpuData = gpuGeneratorState[g];
          const auto cpuComponentRhs = generators[g]->getRightVector();
          for (UInt row = 0; row < 3; ++row) {
            const Real sourceDelta = gpuData.sourceAbc[row] -
                (**cpuComponentRhs)(cpuData.matrixNodeIndices[row], 0);
            generatorRhsError2 += sourceDelta * sourceDelta;
            generatorRhsNorm2 += gpuData.sourceAbc[row] * gpuData.sourceAbc[row];
            for (UInt col = 0; col < 3; ++col) {
              const Real expected = cpuData.conductanceAbc[row * 3 + col];
              const Real delta = gpuData.conductanceAbc[row * 3 + col] - expected;
              generatorMatrixError2 += delta * delta;
              generatorMatrixNorm2 += expected * expected;
              const UInt matrixRow = compactIndex(cpuData.matrixNodeIndices[row]);
              const UInt matrixCol = compactIndex(cpuData.matrixNodeIndices[col]);
              for (UInt k = csrOffsets[matrixRow]; k < csrOffsets[matrixRow + 1]; ++k)
                if (csrColumns[k] == matrixCol) {
                  const Real fixed = cpuGroupMatrix.coeff(matrixRow, matrixCol);
                  const Real stampDelta = values[k] - fixed -
                      gpuData.conductanceAbc[row * 3 + col];
                  globalGeneratorStampError2 += stampDelta * stampDelta;
                  globalGeneratorStampNorm2 +=
                      gpuData.conductanceAbc[row * 3 + col] *
                      gpuData.conductanceAbc[row * 3 + col];
                }
            }
          }
        }
        const Real generatorMatrixError = std::sqrt(generatorMatrixError2 /
            std::max<Real>(generatorMatrixNorm2, 1e-30));
        const Real globalGeneratorStampError = std::sqrt(
            globalGeneratorStampError2 /
            std::max<Real>(globalGeneratorStampNorm2, 1e-30));
        const Real generatorRhsError = std::sqrt(generatorRhsError2 /
            std::max<Real>(generatorRhsNorm2, 1e-30));
        Matrix cpuGroupOnlyRhs = cpuRhs;
        for (const auto &gen : generators) {
          const auto genData = gen->metalDeviceData();
          const auto genRhs = gen->getRightVector();
          for (UInt phase = 0; phase < 3; ++phase) {
            const UInt index = compactIndex(genData.matrixNodeIndices[phase]);
            cpuGroupOnlyRhs(index) -=
                (**genRhs)(genData.matrixNodeIndices[phase], 0);
          }
        }
        Matrix gpuGroupOnlyRhs = Matrix::Zero(order, 1);
        Matrix cpuGroupHistory(groupTerminalIndices.size(), 1);
        for (UInt terminal = 0; terminal < groupTerminalIndices.size(); ++terminal)
          cpuGroupHistory(terminal) =
              cpuGroupOnlyRhs(groupTerminalIndices[terminal]);
        const auto *groupHistory = static_cast<const float *>(
            [(id<MTLBuffer>)groupAdapter->historyBuffer() contents]);
        Matrix gpuGroupHistory(cpuGroupHistory.rows(), 1);
        for (UInt terminal = 0; terminal < groupTerminalIndices.size(); ++terminal)
        {
          gpuGroupHistory(terminal) = groupHistory[terminal];
          gpuGroupOnlyRhs(groupTerminalIndices[terminal]) += groupHistory[terminal];
        }
        const Real groupOnlyRhsError =
            (gpuGroupOnlyRhs - cpuGroupOnlyRhs).norm() /
            std::max<Real>(cpuGroupOnlyRhs.norm(), 1e-30);
        const Real groupHistoryError =
            (gpuGroupHistory - cpuGroupHistory).norm() /
            std::max<Real>(cpuGroupHistory.norm(), 1e-30);
        Real maxGpuMatrixValue = 0.0, maxCpuMatrixValue = 0.0;
        Real maxFixedMatrixValue = 0.0;
        for (UInt k = csrOffsets[maxMatrixRow]; k < csrOffsets[maxMatrixRow + 1]; ++k)
          if (csrColumns[k] == maxMatrixCol) {
            maxGpuMatrixValue = values[k];
            maxCpuMatrixValue = cpuMatrix.coeff(maxMatrixRow, maxMatrixCol);
            maxFixedMatrixValue = cpuGroupMatrix.coeff(maxMatrixRow, maxMatrixCol);
          }
        Real maxDeviceGeneratorValue = 0.0;
        for (std::size_t g = 0; g < generators.size(); ++g) {
          const auto originalIndices = generators[g]->metalDeviceData().matrixNodeIndices;
          std::array<UInt, 3> indices{};
          for (UInt phase = 0; phase < 3; ++phase)
            indices[phase] = compactIndex(originalIndices[phase]);
          for (UInt row = 0; row < 3; ++row)
            for (UInt col = 0; col < 3; ++col)
              if (indices[row] == maxMatrixRow && indices[col] == maxMatrixCol)
                maxDeviceGeneratorValue +=
                    gpuGeneratorState[g].conductanceAbc[row * 3 + col];
        }
        std::cout << "Metal resident timestep-2 oracle: relative matrix error="
                  << matrixError << " (fixed-network=" << fixedMatrixError
                  << ", generator-block=" << generatorMatrixError
                  << ", global-generator-stamp=" << globalGeneratorStampError
                  << "), RHS error=" << rhsError << " (generator="
                  << generatorRhsError << ", group=" << groupOnlyRhsError
                  << ", direct group history=" << groupHistoryError
                  << " (GPU norm=" << gpuGroupHistory.norm()
                  << ", CPU norm=" << cpuGroupHistory.norm() << ')'
                  << "), max matrix delta=" << maxGeneratorMatrixDelta
                  << " at (" << maxMatrixRow << ',' << maxMatrixCol << ')'
                  << " [GPU=" << maxGpuMatrixValue
                  << ", CPU=" << maxCpuMatrixValue
                  << ", fixed=" << maxFixedMatrixValue << ']'
                  << " [device gen=" << maxDeviceGeneratorValue << ']'
                  << ", CPU solution error=" << trueSolutionError << '\n';
        if (!std::isfinite(matrixError) || !std::isfinite(rhsError) ||
            !std::isfinite(trueSolutionError) || matrixError > 1e-4 ||
            rhsError > 1e-4 || trueSolutionError > 1e-3)
          throw std::runtime_error(
              "Metal resident timestep-2 CPU oracle failed");
      }
    }
  }

  void encodeAndRunGpuGeneratedIcb(Real startTime, Int firstStep,
                                   UInt stepCount) {
    if (!useGpuGeneratedIcb || !icbBuilderPipeline || !icbArgumentEncoder ||
        firstStep < 0 || stepCount == 0)
      throw std::invalid_argument("Invalid GPU-generated timestep ICB request");
    constexpr UInt commandsPerStep = 5;
    if (stepCount > std::numeric_limits<UInt>::max() / commandsPerStep)
      throw std::overflow_error("GPU timestep ICB command count overflow");
    const UInt commandCount = stepCount * commandsPerStep;
    const UInt infoBytes = stepCount * sizeof(uint32_t) * 5u;
    [generatedIcb release];
    generatedIcb = nil;
    [icbStepInfos release];
    icbStepInfos = nil;
    icbStepInfos = [device newBufferWithLength:infoBytes
                                        options:MTLResourceStorageModeShared];
    if (!icbStepInfos)
      throw std::runtime_error("GPU timestep ICB step-info buffer allocation failed");

    MTLIndirectCommandBufferDescriptor *descriptor =
        [[MTLIndirectCommandBufferDescriptor alloc] init];
    descriptor.commandTypes = MTLIndirectCommandTypeConcurrentDispatch;
    descriptor.inheritPipelineState = NO;
    descriptor.inheritBuffers = NO;
    descriptor.maxKernelBufferBindCount = 15;
    generatedIcb = [device newIndirectCommandBufferWithDescriptor:descriptor
                                                  maxCommandCount:commandCount
                                                         options:MTLResourceStorageModeShared];
    [descriptor release];
    if (!generatedIcb)
      throw std::runtime_error("GPU timestep ICB allocation failed");

    [icbArgumentEncoder setArgumentBuffer:icbArguments offset:0];
    [icbArgumentEncoder setIndirectCommandBuffer:generatedIcb atIndex:9];
    [icbArgumentEncoder setBuffer:icbStepInfos offset:0 atIndex:42];
    [icbArgumentEncoder setBuffer:icbStepInfos offset:0 atIndex:43];
    const auto group = groupAdapter->indirectResources();
    const auto solve = solver->indirectResources();
    const UInt groupCount = group.stateCount;
    struct BuildInfo {
      uint32_t commandCount, firstStep, generatorCount, order;
      uint32_t traceWidth, clearGroups, historyGroups, traceGroups;
      uint32_t stateGroups, generatorGroups, solverThreads;
      float startTime, timeStep;
    } buildInfo{
        commandCount, static_cast<uint32_t>(firstStep),
        static_cast<uint32_t>(generatorAdapter->generatorCount()), order,
        traceWidth, (order + 127u) / 128u,
        (static_cast<uint32_t>(group.terminalCount) + 3u) / 4u,
        (traceWidth + 127u) / 128u, (groupCount + 3u) / 4u,
        (static_cast<uint32_t>(generatorAdapter->generatorCount()) + 127u) / 128u,
        solve.threadCount, static_cast<float>(startTime), static_cast<float>(dt)};

    id<MTLCommandBuffer> chunk = [queue commandBuffer];
    if (!chunk)
      throw std::runtime_error("Could not allocate GPU timestep ICB command buffer");
    id<MTLComputeCommandEncoder> encoder = [chunk computeCommandEncoder];
    if (!encoder)
      throw std::runtime_error("Could not create GPU timestep ICB builder encoder");
    [encoder setComputePipelineState:icbBuilderPipeline];
    [encoder setBuffer:icbArguments offset:0 atIndex:0];
    [encoder setBytes:&buildInfo length:sizeof(buildInfo) atIndex:1];
    [encoder useResource:generatedIcb usage:MTLResourceUsageWrite];
    [encoder useResource:icbStepInfos usage:MTLResourceUsageWrite];
    const NSUInteger builderWidth = std::min<NSUInteger>(
        128, icbBuilderPipeline.maxTotalThreadsPerThreadgroup);
    [encoder dispatchThreads:MTLSizeMake(commandCount, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(builderWidth, 1, 1)];
    [encoder endEncoding];

    encoder = [chunk computeCommandEncoder];
    if (!encoder)
      throw std::runtime_error("Could not create GPU timestep ICB execution encoder");
    [encoder useResource:generatedIcb usage:MTLResourceUsageRead];
    [encoder useResource:icbStepInfos usage:MTLResourceUsageRead];
    const auto gen = generatorAdapter->indirectResources();
    id<MTLResource> resources[] = {
        (id<MTLResource>)solver->rightSideBuffer(),
        (id<MTLResource>)solver->rowOffsetsBuffer(),
        (id<MTLResource>)solver->columnIndicesBuffer(),
        (id<MTLResource>)solver->matrixValuesBuffer(),
        (id<MTLResource>)solver->solutionBuffer(),
        (id<MTLResource>)solve.work, (id<MTLResource>)solve.status,
        (id<MTLResource>)solve.params,
        (id<MTLResource>)group.historyOffsets,
        (id<MTLResource>)group.historyColumns,
        (id<MTLResource>)group.historyValues,
        (id<MTLResource>)group.historyValuesLow,
        (id<MTLResource>)group.stateOffsets,
        (id<MTLResource>)group.stateColumns,
        (id<MTLResource>)group.stateValues,
        (id<MTLResource>)group.stateValuesLow,
        (id<MTLResource>)group.inputOffsets,
        (id<MTLResource>)group.inputColumns,
        (id<MTLResource>)group.inputValues,
        (id<MTLResource>)group.inputValuesLow,
        (id<MTLResource>)group.terminalIndices,
        (id<MTLResource>)group.state,
        (id<MTLResource>)group.stateLow,
        (id<MTLResource>)group.nextState,
        (id<MTLResource>)group.nextStateLow,
        (id<MTLResource>)group.history,
        (id<MTLResource>)group.params,
        (id<MTLResource>)gen.generators,
        (id<MTLResource>)traceIndices,
        (id<MTLResource>)traceValues};
    for (id<MTLResource> resource : resources)
      [encoder useResource:resource usage:MTLResourceUsageRead | MTLResourceUsageWrite];
    [encoder executeCommandsInBuffer:generatedIcb
                           withRange:NSMakeRange(0, commandCount)];
    [encoder endEncoding];
    [chunk commit];
    [chunk waitUntilCompleted];
    if (chunk.GPUEndTime >= chunk.GPUStartTime && chunk.GPUStartTime > 0.0) {
      completedGpuSeconds += chunk.GPUEndTime - chunk.GPUStartTime;
      ++completedCommandBuffers;
    }
    if (chunk.error)
      throw std::runtime_error(std::string("GPU timestep ICB execution failed: ") +
                               (chunk.error.localizedDescription.UTF8String ?: "unknown error"));

    const auto *status = static_cast<const uint32_t *>(
        [(id<MTLBuffer>)solver->statusBuffer() contents]);
    if (status[0] != 0) {
      float residual = 0.0f, rhsNorm = 0.0f;
      std::memcpy(&residual, &status[2], sizeof(float));
      std::memcpy(&rhsNorm, &status[3], sizeof(float));
      throw std::runtime_error("GPU timestep ICB BiCGSTAB failed (status=" +
          std::to_string(status[0]) + ", relative residual=" +
          std::to_string(residual / std::max(rhsNorm, 1e-30f)) + ")");
    }
    groupAdapter->finishIndirectStateUpdates(stepCount);
    committedSteps += stepCount;
    lastWaitedSteps = committedSteps;
    std::cout << "Metal GPU-generated ICB executed " << stepCount
              << " timesteps using " << commandCount
              << " indirect compute commands in one command buffer\n";
  }

  void encodeStepBatch(Real time, Int firstCount, UInt count) {
    for (UInt offset = 0; offset < count; ++offset)
      encodeSingleStep(time + static_cast<Real>(offset) * dt,
                       firstCount + static_cast<Int>(offset));
  }

  void finish() {
    if (command) {
      [command commit];
      [lastCommitted release];
      lastCommitted = [command retain];
      pendingCommands.push_back([command retain]);
      command = nil;
      committedSteps += encodedSteps;
    }
    if (!pendingCommands.empty())
      waitForPendingCommands();
    reportPhaseTimings();
    if (!sparseInputGpuLu) {
      const auto *iterations =
          static_cast<const uint32_t *>(iterationValues.contents);
      uint64_t iterationTotal = 0;
      uint32_t maxIterations = 0;
      for (UInt step = 0; step < expectedSteps; ++step) {
        iterationTotal += iterations[step];
        maxIterations = std::max(maxIterations, iterations[step]);
      }
      std::cout << "Metal BiCGSTAB iteration counts: mean="
                << static_cast<double>(iterationTotal) / expectedSteps
                << ", max=" << maxIterations << " over " << expectedSteps
                << " steps\n";
    }
    std::cout << "Metal GPU execution timing: " << completedGpuSeconds
              << " s across " << completedCommandBuffers
              << " command buffer(s)\n";
    if (committedSteps != expectedSteps)
      throw std::runtime_error("Metal resident scheduler ended before all expected timesteps");
    const auto *status = static_cast<const uint32_t *>(
        [(id<MTLBuffer>)solver->statusBuffer() contents]);
    if (sparseInputGpuLu)
      sparseInputGpuLu->checkStatus();
    else if (status[0] != 0)
      throw std::runtime_error("Metal resident BiCGSTAB failed on the final timestep");
    const auto traceWriteStart = std::chrono::steady_clock::now();
    std::ofstream output(tracePath);
    if (!output)
      throw std::runtime_error("Could not open Metal resident trace output: " + tracePath);
    output << "step,time_s";
    for (const auto &node : nodes)
      for (int phase = 0; phase < 3; ++phase)
        output << "," << node->name() << "_phase" << phase;
    output << '\n';
    const auto *trace = static_cast<const float *>(traceValues.contents);
    std::string writeBuffer;
    writeBuffer.reserve(1u << 20);
    auto flushBuffer = [&]() {
      output.write(writeBuffer.data(),
                   static_cast<std::streamsize>(writeBuffer.size()));
      writeBuffer.clear();
    };
    auto appendFloat = [&writeBuffer](auto value) {
      char chars[48];
      const auto result = std::to_chars(chars, chars + sizeof(chars), value,
                                        std::chars_format::general, 6);
      if (result.ec != std::errc{})
        throw std::runtime_error("Metal resident trace float formatting failed");
      writeBuffer.append(chars, result.ptr);
    };
    auto appendInteger = [&writeBuffer](UInt value) {
      char chars[24];
      const auto result = std::to_chars(chars, chars + sizeof(chars), value);
      if (result.ec != std::errc{})
        throw std::runtime_error("Metal resident trace integer formatting failed");
      writeBuffer.append(chars, result.ptr);
    };
    for (UInt step = 0; step < expectedSteps; ++step) {
      appendInteger(step);
      writeBuffer.push_back(',');
      appendFloat(static_cast<Real>(step + 1) * dt);
      for (UInt col = 0; col < traceWidth; ++col) {
        writeBuffer.push_back(',');
        appendFloat(trace[static_cast<std::size_t>(step) * traceWidth + col]);
      }
      writeBuffer.push_back('\n');
      if (writeBuffer.size() >= (1u << 20))
        flushBuffer();
    }
    if (!writeBuffer.empty())
      flushBuffer();
    output.close();
    const auto traceWriteSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - traceWriteStart).count();
    std::cout << "Metal resident trace export wall time: "
              << traceWriteSeconds << " s\n";
  }
};

MetalResidentScheduler::MetalResidentScheduler(
    CPS::EMT::LinearSSNGroup::Ptr group,
    std::vector<std::shared_ptr<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>> generators,
    CPS::EMT::SimNode::List nodes, Real timeStep, UInt expectedSteps,
    CPS::String tracePath, Bool wholeRunFastPath, Bool useSparseInputGpuLu,
    Bool useGpuGeneratedIcb)
    : Scheduler(), mImpl(std::make_unique<Impl>(
                          std::move(group), std::move(generators),
                          std::move(nodes), timeStep, expectedSteps,
                          std::move(tracePath), wholeRunFastPath,
                          useSparseInputGpuLu, useGpuGeneratedIcb)) {}

MetalResidentScheduler::~MetalResidentScheduler() = default;

void MetalResidentScheduler::createSchedule(const CPS::Task::List &tasks,
                                             const Edges &, const Edges &) {
  mImpl->initialize(tasks);
}

void MetalResidentScheduler::step(Real time, Int timeStepCount) {
  // The simulation framework invokes Scheduler::step once per timestep. The
  // device-resident path can prepare a bounded block of dependent GPU work at
  // the first callback in that block, then leave the intervening callbacks
  // empty. Keep the first two steps isolated because they run the independent
  // CPU numerical oracles; subsequent blocks are encoded at once.
  if (timeStepCount < mImpl->nextScheduledStep)
    return;
  if (timeStepCount != mImpl->nextScheduledStep)
    throw std::runtime_error("Metal resident scheduler received a discontinuous timestep sequence");
  const UInt batchSize = timeStepCount <= 1
                             ? 1u
                             : std::min<UInt>(64u, mImpl->expectedSteps -
                                                       static_cast<UInt>(timeStepCount));
  if (batchSize == 0)
    return;
  mImpl->encodeStepBatch(time, timeStepCount, batchSize);
  mImpl->nextScheduledStep += static_cast<Int>(batchSize);
}

Bool MetalResidentScheduler::runSimulation(Real startTime, Int firstStep,
                                           UInt stepCount) {
  if (!mImpl->wholeRunFastPath || firstStep != 0 ||
      mImpl->nextScheduledStep != 0 ||
      stepCount != mImpl->expectedSteps)
    return false;

  if (mImpl->useGpuGeneratedIcb) {
    if (stepCount <= 2) {
      if (stepCount > 0)
        mImpl->encodeStepBatch(startTime, 0, stepCount);
      mImpl->nextScheduledStep = static_cast<Int>(stepCount);
      return true;
    }
    // Keep the first two steps on the thoroughly checked path: CPU comparisons
    // verify assembly, solve, component state, and the dynamic generator stamp.
    mImpl->encodeStepBatch(startTime, 0, 2);
    mImpl->nextScheduledStep = 2;
    mImpl->encodeAndRunGpuGeneratedIcb(startTime, 2, stepCount - 2);
    mImpl->nextScheduledStep = static_cast<Int>(stepCount);
    return true;
  }

  while (static_cast<UInt>(mImpl->nextScheduledStep) < stepCount) {
    const Int count = mImpl->nextScheduledStep;
    const UInt remaining = stepCount - static_cast<UInt>(count);
    const UInt batchSize = count <= 1 ? 1u : std::min<UInt>(64u, remaining);
    mImpl->encodeStepBatch(startTime + static_cast<Real>(count) * mImpl->dt,
                           count, batchSize);
    mImpl->nextScheduledStep += static_cast<Int>(batchSize);
  }
  return true;
}

void MetalResidentScheduler::stop() { mImpl->finish(); }

} // namespace DPsim
