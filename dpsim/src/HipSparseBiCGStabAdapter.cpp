/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/HipResidentKernels.h>
#include <dpsim/HipSparseBiCGStabAdapter.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace DPsim {
namespace {
constexpr double relativeTolerance = 1e-10;
constexpr int maximumIterations = 1000;

void release(void *pointer) {
  if (pointer) {
    const auto status = hipFree(pointer);
    (void)status;
  }
}
} // namespace

HipSparseBiCGStabAdapter::HipSparseBiCGStabAdapter()
    : DirectLinearSolver(CPS::Logger::get("HipSparseBiCGStabAdapter",
                                          CPS::Logger::Level::off)) {}

HipSparseBiCGStabAdapter::HipSparseBiCGStabAdapter(CPS::Logger::Log log)
    : DirectLinearSolver(log) {}

HipSparseBiCGStabAdapter::~HipSparseBiCGStabAdapter() {
  release(mRowPointers);
  release(mColumnIndices);
  release(mValues);
  release(mRowScale);
  release(mRightHandSide);
  release(mSolution);
  release(mResidual);
  release(mShadowResidual);
  release(mDirection);
  release(mProduct);
  release(mIntermediate);
  release(mIntermediateProduct);
  release(mUpdatePositions);
  release(mUpdateValues);
  release(mStatus);
}

void HipSparseBiCGStabAdapter::checkHip(hipError_t status,
                                        const char *operation) const {
  if (status != hipSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             hipGetErrorString(status));
}

void HipSparseBiCGStabAdapter::allocate(Int size, Int nonZeros) {
  if (size <= 0 || nonZeros <= 0 || size > std::numeric_limits<int>::max() ||
      nonZeros > std::numeric_limits<int>::max())
    throw std::runtime_error("HIP sparse solver received invalid dimensions");
  if (mSize == size && mNonZeros == nonZeros && mRowPointers)
    return;

  release(mRowPointers);
  release(mColumnIndices);
  release(mValues);
  release(mRowScale);
  release(mRightHandSide);
  release(mSolution);
  release(mResidual);
  release(mShadowResidual);
  release(mDirection);
  release(mProduct);
  release(mIntermediate);
  release(mIntermediateProduct);
  release(mUpdatePositions);
  release(mUpdateValues);
  release(mStatus);
  mRowPointers = nullptr;
  mColumnIndices = nullptr;
  mValues = nullptr;
  mRowScale = nullptr;
  mRightHandSide = nullptr;
  mSolution = nullptr;
  mResidual = nullptr;
  mShadowResidual = nullptr;
  mDirection = nullptr;
  mProduct = nullptr;
  mIntermediate = nullptr;
  mIntermediateProduct = nullptr;
  mUpdatePositions = nullptr;
  mUpdateValues = nullptr;
  mStatus = nullptr;
  mUpdateCapacity = 0;
  mMatrixReady = false;

  mSize = static_cast<int>(size);
  mNonZeros = static_cast<int>(nonZeros);
  checkHip(hipMalloc(&mRowPointers, (mSize + 1) * sizeof(int)),
           "hipMalloc(CSR row pointers)");
  checkHip(hipMalloc(&mColumnIndices, mNonZeros * sizeof(int)),
           "hipMalloc(CSR column indices)");
  checkHip(hipMalloc(&mValues, mNonZeros * sizeof(double)),
           "hipMalloc(CSR values)");
  checkHip(hipMalloc(&mRowScale, mSize * sizeof(double)),
           "hipMalloc(row scale)");
  checkHip(hipMalloc(&mRightHandSide, mSize * sizeof(double)),
           "hipMalloc(RHS)");
  checkHip(hipMalloc(&mSolution, mSize * sizeof(double)),
           "hipMalloc(solution)");
  checkHip(hipMalloc(&mResidual, mSize * sizeof(double)),
           "hipMalloc(residual)");
  checkHip(hipMalloc(&mShadowResidual, mSize * sizeof(double)),
           "hipMalloc(shadow residual)");
  checkHip(hipMalloc(&mDirection, mSize * sizeof(double)),
           "hipMalloc(direction)");
  checkHip(hipMalloc(&mProduct, mSize * sizeof(double)), "hipMalloc(product)");
  checkHip(hipMalloc(&mIntermediate, mSize * sizeof(double)),
           "hipMalloc(intermediate)");
  checkHip(hipMalloc(&mIntermediateProduct, mSize * sizeof(double)),
           "hipMalloc(intermediate product)");
  checkHip(hipMalloc(&mStatus, sizeof(int)), "hipMalloc(solver status)");
}

void HipSparseBiCGStabAdapter::uploadMatrix(const SparseMatrix &systemMatrix,
                                            bool uploadPattern) {
  SparseMatrix compressed = systemMatrix;
  compressed.makeCompressed();
  if (compressed.rows() != compressed.cols())
    throw std::runtime_error("HIP sparse solver requires a square matrix");
  if (compressed.nonZeros() == 0)
    throw std::runtime_error("HIP sparse solver received an empty matrix");

  if (mSize != compressed.rows() || mNonZeros != compressed.nonZeros())
    allocate(static_cast<Int>(compressed.rows()),
             static_cast<Int>(compressed.nonZeros()));

  if (uploadPattern) {
    using StorageIndex = SparseMatrix::StorageIndex;
    static_assert(std::is_same_v<StorageIndex, int>,
                  "HIP sparse adapter expects 32-bit CSR indices");
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
  std::vector<double> rowScale(static_cast<std::size_t>(mSize));
  for (int row = 0; row < mSize; ++row) {
    double rowMaximum = 0.0;
    for (int entry = compressed.outerIndexPtr()[row];
         entry < compressed.outerIndexPtr()[row + 1]; ++entry)
      rowMaximum = std::max(rowMaximum, std::abs(compressed.valuePtr()[entry]));
    if (!std::isfinite(rowMaximum) ||
        rowMaximum <= std::numeric_limits<double>::min())
      throw std::runtime_error("HIP sparse row scaling found an empty or "
                               "invalid row " +
                               std::to_string(row));
    rowScale[row] = 1.0 / rowMaximum;
  }
  checkHip(hipMemcpy(mRowScale, rowScale.data(), mSize * sizeof(double),
                     hipMemcpyHostToDevice),
           "hipMemcpy(row scale)");
  mMatrixReady = true;
}

void HipSparseBiCGStabAdapter::preprocessing(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  uploadMatrix(systemMatrix, true);
}

void HipSparseBiCGStabAdapter::factorize(SparseMatrix &systemMatrix) {
  uploadMatrix(systemMatrix, true);
}

void HipSparseBiCGStabAdapter::refactorize(SparseMatrix &systemMatrix) {
  uploadMatrix(systemMatrix, false);
}

void HipSparseBiCGStabAdapter::partialRefactorize(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  SparseMatrix compressed = systemMatrix;
  compressed.makeCompressed();
  if (compressed.rows() != mSize || compressed.nonZeros() != mNonZeros)
    throw std::runtime_error("HIP sparse partial update structure mismatch");
  std::vector<int> positions;
  std::vector<double> values;
  positions.reserve(listVariableSystemMatrixEntries.size());
  values.reserve(listVariableSystemMatrixEntries.size());
  for (const auto &[row, column] : listVariableSystemMatrixEntries) {
    if (row >= static_cast<UInt>(mSize) || column >= static_cast<UInt>(mSize))
      continue;
    for (int position = compressed.outerIndexPtr()[row];
         position < compressed.outerIndexPtr()[row + 1]; ++position)
      if (compressed.innerIndexPtr()[position] == static_cast<int>(column)) {
        positions.push_back(position);
        values.push_back(compressed.valuePtr()[position]);
        break;
      }
  }
  if (positions.empty())
    return;
  if (positions.size() > mUpdateCapacity) {
    release(mUpdatePositions);
    release(mUpdateValues);
    mUpdatePositions = nullptr;
    mUpdateValues = nullptr;
    checkHip(hipMalloc(&mUpdatePositions, positions.size() * sizeof(int)),
             "hipMalloc(variable CSR positions)");
    checkHip(hipMalloc(&mUpdateValues, positions.size() * sizeof(double)),
             "hipMalloc(variable CSR values)");
    mUpdateCapacity = positions.size();
  }
  checkHip(hipMemcpy(mUpdatePositions, positions.data(),
                     positions.size() * sizeof(int), hipMemcpyHostToDevice),
           "hipMemcpy(variable CSR positions)");
  checkHip(hipMemcpy(mUpdateValues, values.data(),
                     values.size() * sizeof(double), hipMemcpyHostToDevice),
           "hipMemcpy(variable CSR values)");
  hipScatterSparseValues(mValues, mUpdatePositions, mUpdateValues,
                         positions.size());
}

void HipSparseBiCGStabAdapter::uploadRightHandSide(
    const Matrix &rightSideVector) {
  if (rightSideVector.size() != mSize)
    throw std::runtime_error("HIP sparse solver RHS size mismatch");
  checkHip(hipMemcpy(mRightHandSide, rightSideVector.data(),
                     mSize * sizeof(double), hipMemcpyHostToDevice),
           "hipMemcpy(right-hand side)");
}

void HipSparseBiCGStabAdapter::clearDeviceRightHandSide() {
  if (!mRightHandSide || mSize <= 0)
    throw std::runtime_error("HIP sparse solver RHS buffer is not initialized");
  checkHip(hipMemsetAsync(mRightHandSide, 0,
                          static_cast<std::size_t>(mSize) * sizeof(double)),
           "hipMemsetAsync(device RHS)");
}

void HipSparseBiCGStabAdapter::copyDeviceSolutionToHost(Matrix &solution) {
  if (solution.size() != mSize)
    solution.resize(mSize, 1);
  checkHip(hipMemcpy(solution.data(), mSolution, mSize * sizeof(double),
                     hipMemcpyDeviceToHost),
           "hipMemcpy(solution)");
}

void HipSparseBiCGStabAdapter::solveDevice(bool verifyConvergence) {
  if (!mMatrixReady)
    throw std::runtime_error("HIP sparse solver matrix is not initialized");
  hipSolveBiCGStabResident(mRowPointers, mColumnIndices, mValues, mRowScale,
                           mRightHandSide, mSolution, mResidual,
                           mShadowResidual, mDirection, mProduct, mIntermediate,
                           mIntermediateProduct, mSize, maximumIterations,
                           relativeTolerance, mStatus);
  if (!verifyConvergence)
    return;
  int status = 0;
  checkHip(hipMemcpy(&status, mStatus, sizeof(int), hipMemcpyDeviceToHost),
           "hipMemcpy(BiCGSTAB status)");
  if (status <= 0)
    throw std::runtime_error(
        "HIP resident BiCGSTAB failed to converge or broke down (status=" +
        std::to_string(status) + ")");
}

Matrix HipSparseBiCGStabAdapter::solve(Matrix &rightSideVector) {
  if (!mMatrixReady || rightSideVector.rows() != mSize)
    throw std::runtime_error("HIP sparse solver right-hand side size mismatch");
  uploadRightHandSide(rightSideVector);
  solveDevice();
  Matrix solution(mSize, 1);
  copyDeviceSolutionToHost(solution);
  return solution;
}
} // namespace DPsim
