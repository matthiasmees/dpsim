/* SPDX-License-Identifier: MPL-2.0 */
#pragma once

#include <dpsim/DirectLinearSolver.h>
#include <dpsim/HipSparseDeviceSolver.h>

#include <hip/hip_runtime_api.h>
#include <rocsolver/rocsolver.h>

#include <vector>

namespace DPsim {
/// Sparse LU solver that seeds rocSOLVER CSRRF from an Eigen SparseLU
/// factorization, then refactors the fixed CSR pattern on the GPU.
class HipSparseLURefactorAdapter : public DirectLinearSolver,
                                   public HipSparseDeviceSolver {
  int mSize = 0;
  int mNonZeros = 0;
  int mFactorNonZeros = 0;
  int *mRowPointers = nullptr;
  int *mColumnIndices = nullptr;
  double *mValues = nullptr;
  double *mRightHandSide = nullptr;
  double *mSolution = nullptr;
  int *mFactorRowPointers = nullptr;
  int *mFactorColumnIndices = nullptr;
  double *mFactorValues = nullptr;
  int *mRowPivots = nullptr;
  int *mColumnPivots = nullptr;
  std::vector<int> mHostRowPointers;
  std::vector<int> mHostColumnIndices;
  rocblas_handle mHandle = nullptr;
  rocsolver_rfinfo mRfInfo = nullptr;
  bool mMatrixReady = false;
  std::size_t mVerifiedSolveCount = 0;

  void checkHip(hipError_t status, const char *operation) const;
  void checkSolver(rocblas_status status, const char *operation) const;
  void allocate(Int size, Int nonZeros);
  void uploadMatrix(const SparseMatrix &systemMatrix, bool uploadPattern);
  void initializeRefactorization(const SparseMatrix &systemMatrix);

public:
  HipSparseLURefactorAdapter();
  explicit HipSparseLURefactorAdapter(CPS::Logger::Log log);
  ~HipSparseLURefactorAdapter() override;

  double *deviceRightHandSide() noexcept override { return mRightHandSide; }
  double *deviceSolution() noexcept override { return mSolution; }
  double *deviceValues() noexcept override { return mValues; }
  void clearDeviceRightHandSide() override;
  void uploadRightHandSide(const Matrix &rightSideVector) override;
  void copyDeviceSolutionToHost(Matrix &solution) override;
  void refactorizeDevice() override;
  void solveDevice(bool verifyConvergence = true) override;

  void preprocessing(SparseMatrix &systemMatrix,
                     std::vector<std::pair<UInt, UInt>>
                         &listVariableSystemMatrixEntries) override;
  void factorize(SparseMatrix &systemMatrix) override;
  void refactorize(SparseMatrix &systemMatrix) override;
  void partialRefactorize(SparseMatrix &systemMatrix,
                          std::vector<std::pair<UInt, UInt>>
                              &listVariableSystemMatrixEntries) override;
  Matrix solve(Matrix &rightSideVector) override;
};
} // namespace DPsim
