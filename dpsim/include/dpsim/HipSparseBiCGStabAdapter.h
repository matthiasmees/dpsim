/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim/DirectLinearSolver.h>
#include <dpsim/HipSparseDeviceSolver.h>

#include <hip/hip_runtime_api.h>

#include <type_traits>

namespace DPsim {
/// Sparse HIP BiCGSTAB solver. CSR structure and all iteration vectors remain
/// on the device; the host-facing solve() is retained for existing components.
class HipSparseBiCGStabAdapter : public DirectLinearSolver,
                                 public HipSparseDeviceSolver {
  int mSize = 0;
  int mNonZeros = 0;
  int *mRowPointers = nullptr;
  int *mColumnIndices = nullptr;
  double *mValues = nullptr;
  double *mRowScale = nullptr;
  double *mRightHandSide = nullptr;
  double *mSolution = nullptr;
  double *mResidual = nullptr;
  double *mShadowResidual = nullptr;
  double *mDirection = nullptr;
  double *mProduct = nullptr;
  double *mIntermediate = nullptr;
  double *mIntermediateProduct = nullptr;
  int *mStatus = nullptr;
  int *mUpdatePositions = nullptr;
  double *mUpdateValues = nullptr;
  std::size_t mUpdateCapacity = 0;
  bool mMatrixReady = false;

  void allocate(Int size, Int nonZeros);
  void uploadMatrix(const SparseMatrix &systemMatrix, bool uploadPattern);
  void checkHip(hipError_t status, const char *operation) const;

public:
  using DirectLinearSolver::DirectLinearSolver;
  HipSparseBiCGStabAdapter();
  explicit HipSparseBiCGStabAdapter(CPS::Logger::Log log);
  ~HipSparseBiCGStabAdapter() override;

  double *deviceRightHandSide() noexcept override { return mRightHandSide; }
  double *deviceSolution() noexcept override { return mSolution; }
  double *deviceValues() noexcept override { return mValues; }
  void clearDeviceRightHandSide() override;
  void uploadRightHandSide(const Matrix &rightSideVector) override;
  void copyDeviceSolutionToHost(Matrix &solution) override;
  void refactorizeDevice() override {}

  void preprocessing(SparseMatrix &systemMatrix,
                     std::vector<std::pair<UInt, UInt>>
                         &listVariableSystemMatrixEntries) override;
  void factorize(SparseMatrix &systemMatrix) override;
  void refactorize(SparseMatrix &systemMatrix) override;
  void partialRefactorize(SparseMatrix &systemMatrix,
                          std::vector<std::pair<UInt, UInt>>
                              &listVariableSystemMatrixEntries) override;
  void solveDevice(bool verifyConvergence = true) override;
  Matrix solve(Matrix &rightSideVector) override;
};
} // namespace DPsim
