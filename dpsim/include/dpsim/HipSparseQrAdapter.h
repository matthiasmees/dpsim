/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim/DirectLinearSolver.h>

#include <hip/hip_runtime_api.h>
#include <hipsolver/hipsolver.h>
#include <hipsparse/hipsparse.h>

namespace DPsim {
/// HIP sparse QR reference adapter. The matrix and vectors are stored on the
/// device; hipSOLVER performs its sparse QR factorization and solve together.
/// Unlike KLU this API does not expose reusable numeric factors.
class HipSparseQrAdapter : public DirectLinearSolver {
  hipsolverSpHandle_t mSolver = nullptr;
  hipsparseMatDescr_t mMatrixDescriptor = nullptr;
  int mSize = 0;
  int mNonZeros = 0;
  int *mRowPointers = nullptr;
  int *mColumnIndices = nullptr;
  double *mValues = nullptr;
  double *mRightHandSide = nullptr;
  double *mSolution = nullptr;
  bool mMatrixReady = false;

  void allocate(Int size, Int nonZeros);
  void uploadMatrix(const SparseMatrix &systemMatrix, bool uploadPattern);
  void checkHip(hipError_t status, const char *operation) const;
  void checkSolver(hipsolverStatus_t status, const char *operation) const;
  void checkSparse(hipsparseStatus_t status, const char *operation) const;

public:
  using DirectLinearSolver::DirectLinearSolver;
  HipSparseQrAdapter();
  explicit HipSparseQrAdapter(CPS::Logger::Log log);
  ~HipSparseQrAdapter() override;

  double *deviceRightHandSide() noexcept { return mRightHandSide; }
  double *deviceSolution() noexcept { return mSolution; }

  void preprocessing(SparseMatrix &systemMatrix,
                     std::vector<std::pair<UInt, UInt>>
                         &listVariableSystemMatrixEntries) override;
  void factorize(SparseMatrix &systemMatrix) override;
  void refactorize(SparseMatrix &systemMatrix) override;
  void partialRefactorize(SparseMatrix &systemMatrix,
                          std::vector<std::pair<UInt, UInt>>
                              &listVariableSystemMatrixEntries) override;
  void solveDevice();
  Matrix solve(Matrix &rightSideVector) override;
};
} // namespace DPsim
