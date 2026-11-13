/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#pragma once

#include <memory>

#include <dpsim/DirectLinearSolver.h>

namespace DPsim {

/// Experimental Apple Metal dense LU adapter. The current MPS matrix API uses
/// single-precision values; use only for benchmarking until application-level
/// accuracy has been established.
class MetalDenseAdapter : public DirectLinearSolver {
public:
  explicit MetalDenseAdapter(CPS::Logger::Log log = nullptr);
  ~MetalDenseAdapter() override;

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

private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;
};

} // namespace DPsim
