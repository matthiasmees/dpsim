/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/HipSparseQrAdapter.h>

#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace DPsim {
namespace {
constexpr double singularityTolerance = 1e-12;

void release(void *pointer) {
  if (pointer) {
    const auto status = hipFree(pointer);
    (void)status;
  }
}
} // namespace

HipSparseQrAdapter::HipSparseQrAdapter()
    : DirectLinearSolver(
          CPS::Logger::get("HipSparseQrAdapter", CPS::Logger::Level::off)) {
  checkSolver(hipsolverSpCreate(&mSolver), "hipsolverSpCreate");
  checkSparse(hipsparseCreateMatDescr(&mMatrixDescriptor),
              "hipsparseCreateMatDescr");
  checkSparse(
      hipsparseSetMatType(mMatrixDescriptor, HIPSPARSE_MATRIX_TYPE_GENERAL),
      "hipsparseSetMatType");
  checkSparse(
      hipsparseSetMatIndexBase(mMatrixDescriptor, HIPSPARSE_INDEX_BASE_ZERO),
      "hipsparseSetMatIndexBase");
}

HipSparseQrAdapter::HipSparseQrAdapter(CPS::Logger::Log log)
    : DirectLinearSolver(log) {
  checkSolver(hipsolverSpCreate(&mSolver), "hipsolverSpCreate");
  checkSparse(hipsparseCreateMatDescr(&mMatrixDescriptor),
              "hipsparseCreateMatDescr");
  checkSparse(
      hipsparseSetMatType(mMatrixDescriptor, HIPSPARSE_MATRIX_TYPE_GENERAL),
      "hipsparseSetMatType");
  checkSparse(
      hipsparseSetMatIndexBase(mMatrixDescriptor, HIPSPARSE_INDEX_BASE_ZERO),
      "hipsparseSetMatIndexBase");
}

HipSparseQrAdapter::~HipSparseQrAdapter() {
  release(mRowPointers);
  release(mColumnIndices);
  release(mValues);
  release(mRightHandSide);
  release(mSolution);
  if (mMatrixDescriptor)
    hipsparseDestroyMatDescr(mMatrixDescriptor);
  if (mSolver)
    hipsolverSpDestroy(mSolver);
}

void HipSparseQrAdapter::checkHip(hipError_t status,
                                  const char *operation) const {
  if (status != hipSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             hipGetErrorString(status));
}

void HipSparseQrAdapter::checkSolver(hipsolverStatus_t status,
                                     const char *operation) const {
  if (status != HIPSOLVER_STATUS_SUCCESS)
    throw std::runtime_error(std::string(operation) +
                             " failed with hipSOLVER status " +
                             std::to_string(static_cast<int>(status)));
}

void HipSparseQrAdapter::checkSparse(hipsparseStatus_t status,
                                     const char *operation) const {
  if (status != HIPSPARSE_STATUS_SUCCESS)
    throw std::runtime_error(std::string(operation) +
                             " failed with hipSPARSE status " +
                             std::to_string(static_cast<int>(status)));
}

void HipSparseQrAdapter::allocate(Int size, Int nonZeros) {
  if (size <= 0 || nonZeros <= 0 || size > std::numeric_limits<int>::max() ||
      nonZeros > std::numeric_limits<int>::max())
    throw std::runtime_error("HIP sparse QR received invalid dimensions");
  if (mSize == size && mNonZeros == nonZeros && mRowPointers)
    return;

  release(mRowPointers);
  release(mColumnIndices);
  release(mValues);
  release(mRightHandSide);
  release(mSolution);
  mRowPointers = nullptr;
  mColumnIndices = nullptr;
  mValues = nullptr;
  mRightHandSide = nullptr;
  mSolution = nullptr;
  mMatrixReady = false;
  mSize = static_cast<int>(size);
  mNonZeros = static_cast<int>(nonZeros);
  checkHip(hipMalloc(&mRowPointers, (mSize + 1) * sizeof(int)),
           "hipMalloc(CSR row pointers)");
  checkHip(hipMalloc(&mColumnIndices, mNonZeros * sizeof(int)),
           "hipMalloc(CSR column indices)");
  checkHip(hipMalloc(&mValues, mNonZeros * sizeof(double)),
           "hipMalloc(CSR values)");
  checkHip(hipMalloc(&mRightHandSide, mSize * sizeof(double)),
           "hipMalloc(RHS)");
  checkHip(hipMalloc(&mSolution, mSize * sizeof(double)),
           "hipMalloc(solution)");
}

void HipSparseQrAdapter::uploadMatrix(const SparseMatrix &systemMatrix,
                                      bool uploadPattern) {
  SparseMatrix compressed = systemMatrix;
  compressed.makeCompressed();
  if (compressed.rows() != compressed.cols() || compressed.nonZeros() == 0)
    throw std::runtime_error(
        "HIP sparse QR requires a non-empty square matrix");

  const bool dimensionsChanged = mSize != compressed.rows() ||
                                 mNonZeros != compressed.nonZeros() ||
                                 !mRowPointers;
  if (dimensionsChanged)
    allocate(static_cast<Int>(compressed.rows()),
             static_cast<Int>(compressed.nonZeros()));

  if (uploadPattern || dimensionsChanged) {
    using StorageIndex = SparseMatrix::StorageIndex;
    static_assert(std::is_same_v<StorageIndex, int>,
                  "HIP sparse QR expects 32-bit CSR indices");
    checkHip(hipMemcpy(mRowPointers, compressed.outerIndexPtr(),
                       (mSize + 1) * sizeof(int), hipMemcpyHostToDevice),
             "hipMemcpy(CSR row pointers)");
    checkHip(hipMemcpy(mColumnIndices, compressed.innerIndexPtr(),
                       mNonZeros * sizeof(int), hipMemcpyHostToDevice),
             "hipMemcpy(CSR column indices)");
  }
  checkHip(hipMemcpy(mValues, compressed.valuePtr(), mNonZeros * sizeof(double),
                     hipMemcpyHostToDevice),
           "hipMemcpy(CSR values)");
  mMatrixReady = true;
}

void HipSparseQrAdapter::preprocessing(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  uploadMatrix(systemMatrix, true);
}

void HipSparseQrAdapter::factorize(SparseMatrix &systemMatrix) {
  uploadMatrix(systemMatrix, true);
}

void HipSparseQrAdapter::refactorize(SparseMatrix &systemMatrix) {
  uploadMatrix(systemMatrix, false);
}

void HipSparseQrAdapter::partialRefactorize(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  uploadMatrix(systemMatrix, false);
}

void HipSparseQrAdapter::solveDevice() {
  if (!mMatrixReady)
    throw std::runtime_error("HIP sparse QR matrix is not initialized");

  int singularity = -1;
  checkSolver(hipsolverSpDcsrlsvqr(mSolver, mSize, mNonZeros, mMatrixDescriptor,
                                   mValues, mRowPointers, mColumnIndices,
                                   mRightHandSide, singularityTolerance, 0,
                                   mSolution, &singularity),
              "hipsolverSpDcsrlsvqr");
  if (singularity >= 0)
    throw std::runtime_error("HIP sparse QR found a singularity at row " +
                             std::to_string(singularity));
}

Matrix HipSparseQrAdapter::solve(Matrix &rightSideVector) {
  if (!mMatrixReady || rightSideVector.rows() != mSize)
    throw std::runtime_error("HIP sparse QR right-hand side size mismatch");
  checkHip(hipMemcpy(mRightHandSide, rightSideVector.data(),
                     mSize * sizeof(double), hipMemcpyHostToDevice),
           "hipMemcpy(right-hand side)");
  solveDevice();
  Matrix solution(mSize, 1);
  checkHip(hipMemcpy(solution.data(), mSolution, mSize * sizeof(double),
                     hipMemcpyDeviceToHost),
           "hipMemcpy(solution)");
  return solution;
}
} // namespace DPsim
