/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#pragma once

#include <cstdint>
#include <memory>

#include <dpsim/DirectLinearSolver.h>

namespace DPsim {

/// Experimental Metal preconditioned BiCGSTAB solver for sparse MNA matrices.
/// Its device-side encode API allows a
/// Metal timestep executor to solve using buffers that remain GPU-resident.
/// The DirectLinearSolver API is retained as a CPU-facing validation adapter.
class MetalSparseBiCGSTABAdapter : public DirectLinearSolver {
public:
  struct IndirectResources {
    void *pipeline;
    void *work;
    void *status;
    void *params;
    uint32_t order;
    uint32_t threadCount;
  };
  explicit MetalSparseBiCGSTABAdapter(CPS::Logger::Log log = nullptr);
  ~MetalSparseBiCGSTABAdapter() override;

  /// Configure the device Krylov iteration cap and relative/absolute residual
  /// tolerances. Must be called before preprocessing; defaults target the
  /// accuracy used by the CPU-facing validation adapter.
  void setConvergenceCriteria(uint32_t maxIterations, float relativeTolerance,
                              float absoluteTolerance);

  void preprocessing(SparseMatrix &systemMatrix,
                     std::vector<std::pair<UInt, UInt>>
                         &listVariableSystemMatrixEntries) override;
  void factorize(SparseMatrix &systemMatrix) override;
  void refactorize(SparseMatrix &systemMatrix) override;
  void partialRefactorize(SparseMatrix &systemMatrix,
                          std::vector<std::pair<UInt, UInt>>
                              &listVariableSystemMatrixEntries) override;
  Matrix solve(Matrix &rightSideVector) override;
  void solveInto(Matrix &rightSideVector, Matrix &solution) override;

  /// Encode BiCGSTAB into an existing Metal command buffer. Matrix values, RHS,
  /// solution and status are Metal buffers; no CPU synchronization/readback is
  /// performed by this method. Static CSR row/column buffers are owned by this
  /// adapter and can be passed to companion Metal component kernels.
  void encodeSolve(void *commandBuffer, void *matrixValuesBuffer,
                   void *rightSideBuffer, void *solutionBuffer);
  /// Clear the device solution vector in-command-buffer, e.g. when the
  /// caller wants a zero initial Krylov guess for a timestep.
  void encodeZeroInitialGuess(void *commandBuffer);
  void *rowOffsetsBuffer() const;
  void *columnIndicesBuffer() const;
  void *matrixValuesBuffer() const;
  void *rightSideBuffer() const;
  void *solutionBuffer() const;
  void *statusBuffer() const;
  void *commandQueue() const;
  IndirectResources indirectResources() const;

private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;
};

} // namespace DPsim
