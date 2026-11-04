/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/HipDenseAdapter.h>

#include <stdexcept>
#include <string>

namespace DPsim {
namespace {
void checkHip(hipError_t status, const char *operation) {
  if (status != hipSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             hipGetErrorString(status));
}

void checkSolver(hipsolverStatus_t status, const char *operation) {
  if (status != HIPSOLVER_STATUS_SUCCESS)
    throw std::runtime_error(std::string(operation) +
                             " failed with hipSOLVER status " +
                             std::to_string(static_cast<int>(status)));
}

void releaseDevicePointer(void *pointer) {
  if (pointer) {
    const auto status = hipFree(pointer);
    (void)status;
  }
}
} // namespace

HipDenseAdapter::HipDenseAdapter() {
  checkSolver(hipsolverDnCreate(&mSolver), "hipsolverDnCreate");
}

HipDenseAdapter::HipDenseAdapter(CPS::Logger::Log log)
    : DirectLinearSolver(log) {
  checkSolver(hipsolverDnCreate(&mSolver), "hipsolverDnCreate");
}

HipDenseAdapter::~HipDenseAdapter() {
  if (mSolver)
    hipsolverDnDestroy(mSolver);
  releaseDevicePointer(mMatrix);
  releaseDevicePointer(mRightHandSide);
  releaseDevicePointer(mWorkspace);
  releaseDevicePointer(mPivotSequence);
  releaseDevicePointer(mDeviceInfo);
}

void HipDenseAdapter::allocateDeviceMemory(UInt size) {
  if (size == 0 || size > static_cast<UInt>(std::numeric_limits<int>::max()))
    throw std::runtime_error(
        "HIP dense solver received an invalid matrix size");

  if (mSize == size)
    return;

  releaseDevicePointer(mMatrix);
  releaseDevicePointer(mRightHandSide);
  releaseDevicePointer(mWorkspace);
  releaseDevicePointer(mPivotSequence);
  releaseDevicePointer(mDeviceInfo);

  mMatrix = nullptr;
  mRightHandSide = nullptr;
  mWorkspace = nullptr;
  mPivotSequence = nullptr;
  mDeviceInfo = nullptr;

  mSize = size;
  const auto matrixBytes =
      static_cast<std::size_t>(size) * size * sizeof(double);
  checkHip(hipMalloc(&mMatrix, matrixBytes), "hipMalloc(matrix)");
  checkHip(hipMalloc(&mRightHandSide, size * sizeof(double)),
           "hipMalloc(right-hand side)");
  checkHip(hipMalloc(&mPivotSequence, size * sizeof(int)), "hipMalloc(pivots)");
  checkHip(hipMalloc(&mDeviceInfo, sizeof(int)), "hipMalloc(info)");

  checkSolver(hipsolverDnDgetrf_bufferSize(
                  mSolver, static_cast<int>(size), static_cast<int>(size),
                  mMatrix, static_cast<int>(size), &mWorkspaceSize),
              "hipsolverDnDgetrf_bufferSize");
  if (mWorkspaceSize > 0) {
    checkHip(hipMalloc(&mWorkspace, static_cast<std::size_t>(mWorkspaceSize) *
                                        sizeof(double)),
             "hipMalloc(workspace)");
  }
}

void HipDenseAdapter::copySystemMatrixToDevice(
    const SparseMatrix &systemMatrix) {
  if (systemMatrix.rows() != systemMatrix.cols())
    throw std::runtime_error("HIP dense solver requires a square matrix");
  allocateDeviceMemory(static_cast<UInt>(systemMatrix.rows()));

  Matrix denseMatrix(systemMatrix);
  checkHip(hipMemcpy(mMatrix, denseMatrix.data(),
                     static_cast<std::size_t>(mSize) * mSize * sizeof(double),
                     hipMemcpyHostToDevice),
           "hipMemcpy(system matrix)");
}

void HipDenseAdapter::factorizeOnDevice() {
  checkSolver(hipsolverDnDgetrf(mSolver, static_cast<int>(mSize),
                                static_cast<int>(mSize), mMatrix,
                                static_cast<int>(mSize), mWorkspace,
                                mPivotSequence, mDeviceInfo),
              "hipsolverDnDgetrf");
  checkHip(hipDeviceSynchronize(), "hipDeviceSynchronize(factorization)");

  int info = 0;
  checkHip(hipMemcpy(&info, mDeviceInfo, sizeof(info), hipMemcpyDeviceToHost),
           "hipMemcpy(factorization info)");
  if (info != 0)
    throw std::runtime_error("HIP dense LU factorization failed (info=" +
                             std::to_string(info) + ")");
}

void HipDenseAdapter::preprocessing(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  allocateDeviceMemory(static_cast<UInt>(systemMatrix.rows()));
}

void HipDenseAdapter::factorize(SparseMatrix &systemMatrix) {
  copySystemMatrixToDevice(systemMatrix);
  factorizeOnDevice();
}

void HipDenseAdapter::refactorize(SparseMatrix &systemMatrix) {
  factorize(systemMatrix);
}

void HipDenseAdapter::partialRefactorize(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  factorize(systemMatrix);
}

void HipDenseAdapter::solveDevice() {
  if (mSize == 0 || !mRightHandSide)
    throw std::runtime_error("HIP dense solver has no device right-hand side");

  checkSolver(hipsolverDnDgetrs(mSolver, HIPBLAS_OP_N, static_cast<int>(mSize),
                                1, mMatrix, static_cast<int>(mSize),
                                mPivotSequence, mRightHandSide,
                                static_cast<int>(mSize), mDeviceInfo),
              "hipsolverDnDgetrs");
}

void HipDenseAdapter::copyDeviceSolutionToHost(Matrix &solution) {
  if (mSize == 0 || !mRightHandSide ||
      solution.rows() != static_cast<Eigen::Index>(mSize))
    throw std::runtime_error("HIP dense solver solution size mismatch");

  checkHip(hipDeviceSynchronize(), "hipDeviceSynchronize(solve)");

  int info = 0;
  checkHip(hipMemcpy(&info, mDeviceInfo, sizeof(info), hipMemcpyDeviceToHost),
           "hipMemcpy(solve info)");
  if (info != 0)
    throw std::runtime_error(
        "HIP dense solve failed (info=" + std::to_string(info) + ")");

  checkHip(hipMemcpy(solution.data(), mRightHandSide, mSize * sizeof(double),
                     hipMemcpyDeviceToHost),
           "hipMemcpy(solution)");
}

Matrix HipDenseAdapter::solve(Matrix &rightSideVector) {
  if (mSize == 0 || rightSideVector.rows() != static_cast<Eigen::Index>(mSize))
    throw std::runtime_error("HIP dense solver right-hand side size mismatch");

  Matrix solution = rightSideVector;
  checkHip(hipMemcpy(mRightHandSide, rightSideVector.data(),
                     mSize * sizeof(double), hipMemcpyHostToDevice),
           "hipMemcpy(right-hand side)");
  solveDevice();
  copyDeviceSolutionToHost(solution);
  return solution;
}
} // namespace DPsim
