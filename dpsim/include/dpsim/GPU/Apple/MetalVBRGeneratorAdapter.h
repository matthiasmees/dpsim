/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator4OrderVBR.h>

namespace DPsim {

/// Device-resident pre/post-step kernels for the controller-free IEEE SG4 VBR
/// benchmark generators. All state is uploaded once; each step updates the
/// VBR conductance/RHS before MNA and consumes solved terminal voltages after.
class MetalVBRGeneratorAdapter {
public:
  using Model = CPS::EMT::Ph3::SynchronGenerator4OrderVBR;
  struct IndirectResources {
    void *prePipeline;
    void *postPipeline;
    void *generators;
    uint32_t count;
  };
  struct StampPositions {
    std::array<uint32_t, 9> positions;
    std::array<float, 9> fixedValues;
  };

  MetalVBRGeneratorAdapter(const std::vector<Model::MetalDeviceData> &models,
                           const std::vector<StampPositions> &stampPositions);
  ~MetalVBRGeneratorAdapter();

  MetalVBRGeneratorAdapter(const MetalVBRGeneratorAdapter &) = delete;
  MetalVBRGeneratorAdapter &
  operator=(const MetalVBRGeneratorAdapter &) = delete;

  /// Encode mechanical/electrical pre-step, dynamic 3x3 MNA stamps, and Norton
  /// source RHS accumulation. Generator terminal buses must be disjoint.
  void encodePreStep(void *commandBuffer, void *matrixValuesBuffer,
                     void *rightSideBuffer, double time);

  /// Encode terminal-voltage/current and dq state updates after the MNA solve.
  void encodePostStep(void *commandBuffer, void *solutionBuffer);
  /// Resident-run variant that also records terminal voltages and the solve's
  /// iteration count in the existing trace buffers.
  void encodePostStepWithTrace(void *commandBuffer, void *solutionBuffer,
                               void *traceValuesBuffer,
                               void *solverStatusBuffer,
                               void *iterationValuesBuffer, uint32_t stepIndex,
                               uint32_t traceWidth);

  /// Copy final device state to the host. Call only after the producing Metal
  /// command buffer has completed (e.g. for end-of-run logging/validation).
  std::vector<Model::MetalDeviceData> copyStateToHost() const;

  /// Borrowed Metal objects, valid for this adapter's lifetime.
  void *stateBuffer() const;
  void *commandQueue() const;
  std::size_t generatorCount() const;
  IndirectResources indirectResources() const;

private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;
};

} // namespace DPsim
