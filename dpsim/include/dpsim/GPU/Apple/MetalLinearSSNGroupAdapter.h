/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#pragma once

#include <cstdint>
#include <memory>

#include <dpsim-models/EMT/EMT_LinearSSNGroup.h>

namespace DPsim {

/// Device-resident timestep kernels for a fully reduced linear SSN group.
/// The group's discrete maps and state are uploaded once; history/RHS assembly
/// and state updates can then be encoded on the same command buffer as the MNA
/// solve without copying state or voltages back to the CPU.
class MetalLinearSSNGroupAdapter {
public:
  struct IndirectResources {
    void *historyPipeline;
    void *updatePipeline;
    void *historyOffsets;
    void *historyColumns;
    void *historyValues;
    void *historyValuesLow;
    void *stateOffsets;
    void *stateColumns;
    void *stateValues;
    void *stateValuesLow;
    void *inputOffsets;
    void *inputColumns;
    void *inputValues;
    void *inputValuesLow;
    void *terminalIndices;
    void *state;
    void *stateLow;
    void *nextState;
    void *nextStateLow;
    void *history;
    void *params;
    uint32_t stateCount;
    uint32_t terminalCount;
  };

  explicit MetalLinearSSNGroupAdapter(
      const CPS::EMT::LinearSSNGroup::PrecomputedStateMaps &model);
  ~MetalLinearSSNGroupAdapter();

  MetalLinearSSNGroupAdapter(const MetalLinearSSNGroupAdapter &) = delete;
  MetalLinearSSNGroupAdapter &
  operator=(const MetalLinearSSNGroupAdapter &) = delete;

  /// Encode h = historyMap*x and add h at its global MNA rows in rhsBuffer.
  /// rhsBuffer must be a Metal float buffer of the full MNA order. The caller
  /// is responsible for zeroing/assembling the rest of the RHS beforehand.
  void encodeHistoryAndRhs(void *commandBuffer, void *rhsBuffer);

  /// Encode x+ = stateMap*x + inputMap*u, gathering u from the solved full
  /// MNA vector. No command-buffer wait or host readback is performed.
  void encodeStateUpdate(void *commandBuffer, void *solutionBuffer);

  /// Borrowed Metal objects, valid for this adapter's lifetime.
  void *stateBuffer() const;
  void *historyBuffer() const;
  void *terminalIndicesBuffer() const;
  void *terminalAdmittanceBuffer() const;
  void *commandQueue() const;
  std::size_t stateCount() const;
  std::size_t terminalCount() const;
  IndirectResources indirectResources() const;
  void finishIndirectStateUpdates(std::size_t updateCount);
  std::vector<float> copyStateToHost() const;

private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;
};

} // namespace DPsim
