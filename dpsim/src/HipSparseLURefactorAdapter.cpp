/* SPDX-License-Identifier: MPL-2.0 */
#include <dpsim/HipSparseLURefactorAdapter.h>

#include <klu.h>

#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace DPsim {
namespace {
using ColSparseMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

void release(void *pointer) {
  if (pointer) {
    const auto status = hipFree(pointer);
    (void)status;
  }
}
} // namespace

HipSparseLURefactorAdapter::HipSparseLURefactorAdapter()
    : DirectLinearSolver(CPS::Logger::get("HipSparseLURefactorAdapter",
                                          CPS::Logger::Level::off)) {
  checkSolver(rocblas_create_handle(&mHandle), "rocblas_create_handle");
  try {
    checkSolver(rocsolver_create_rfinfo(&mRfInfo, mHandle),
                "rocsolver_create_rfinfo");
  } catch (...) {
    rocblas_destroy_handle(mHandle);
    mHandle = nullptr;
    throw;
  }
}

HipSparseLURefactorAdapter::HipSparseLURefactorAdapter(CPS::Logger::Log log)
    : DirectLinearSolver(log) {
  checkSolver(rocblas_create_handle(&mHandle), "rocblas_create_handle");
  try {
    checkSolver(rocsolver_create_rfinfo(&mRfInfo, mHandle),
                "rocsolver_create_rfinfo");
  } catch (...) {
    rocblas_destroy_handle(mHandle);
    mHandle = nullptr;
    throw;
  }
}

HipSparseLURefactorAdapter::~HipSparseLURefactorAdapter() {
  release(mRowPointers);
  release(mColumnIndices);
  release(mValues);
  release(mRightHandSide);
  release(mSolution);
  release(mFactorRowPointers);
  release(mFactorColumnIndices);
  release(mFactorValues);
  release(mRowPivots);
  release(mColumnPivots);
  if (mRfInfo)
    rocsolver_destroy_rfinfo(mRfInfo);
  if (mHandle)
    rocblas_destroy_handle(mHandle);
}

void HipSparseLURefactorAdapter::checkHip(hipError_t status,
                                          const char *operation) const {
  if (status != hipSuccess)
    throw std::runtime_error(std::string(operation) + ": " +
                             hipGetErrorString(status));
}

void HipSparseLURefactorAdapter::checkSolver(rocblas_status status,
                                             const char *operation) const {
  if (status != rocblas_status_success)
    throw std::runtime_error(std::string(operation) +
                             " failed with rocBLAS/rocSOLVER status " +
                             std::to_string(static_cast<int>(status)));
}

void HipSparseLURefactorAdapter::allocate(Int size, Int nonZeros) {
  if (size <= 0 || nonZeros <= 0 || size > std::numeric_limits<int>::max() ||
      nonZeros > std::numeric_limits<int>::max())
    throw std::runtime_error("rocSOLVER sparse LU received invalid dimensions");
  if (mSize == size && mNonZeros == nonZeros && mRowPointers)
    return;

  release(mRowPointers);
  release(mColumnIndices);
  release(mValues);
  release(mRightHandSide);
  release(mSolution);
  release(mFactorRowPointers);
  release(mFactorColumnIndices);
  release(mFactorValues);
  release(mRowPivots);
  release(mColumnPivots);
  mRowPointers = nullptr;
  mColumnIndices = nullptr;
  mValues = nullptr;
  mRightHandSide = nullptr;
  mSolution = nullptr;
  mFactorRowPointers = nullptr;
  mFactorColumnIndices = nullptr;
  mFactorValues = nullptr;
  mRowPivots = nullptr;
  mColumnPivots = nullptr;
  mFactorNonZeros = 0;
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
  checkHip(hipMalloc(&mRowPivots, mSize * sizeof(int)),
           "hipMalloc(row pivots)");
  checkHip(hipMalloc(&mColumnPivots, mSize * sizeof(int)),
           "hipMalloc(column pivots)");
}

void HipSparseLURefactorAdapter::uploadMatrix(const SparseMatrix &systemMatrix,
                                              bool uploadPattern) {
  SparseMatrix compressed = systemMatrix;
  compressed.makeCompressed();
  if (compressed.rows() != compressed.cols() || compressed.nonZeros() == 0)
    throw std::runtime_error(
        "rocSOLVER sparse LU requires a non-empty square matrix");
  if (mSize != compressed.rows() || mNonZeros != compressed.nonZeros())
    allocate(static_cast<Int>(compressed.rows()),
             static_cast<Int>(compressed.nonZeros()));

  if (uploadPattern) {
    using StorageIndex = SparseMatrix::StorageIndex;
    static_assert(std::is_same_v<StorageIndex, int>,
                  "rocSOLVER sparse LU expects 32-bit CSR indices");
    mHostRowPointers.assign(compressed.outerIndexPtr(),
                            compressed.outerIndexPtr() + mSize + 1);
    mHostColumnIndices.assign(compressed.innerIndexPtr(),
                              compressed.innerIndexPtr() + mNonZeros);
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
}

void HipSparseLURefactorAdapter::preprocessing(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  (void)listVariableSystemMatrixEntries;
  uploadMatrix(systemMatrix, true);
}

void HipSparseLURefactorAdapter::initializeRefactorization(
    const SparseMatrix &systemMatrix) {
  ColSparseMatrix matrix = systemMatrix;
  matrix.makeCompressed();
  klu_common common{};
  if (!klu_defaults(&common))
    throw std::runtime_error("KLU could not initialize sparse LU settings");
  common.btf = 0;
  common.scale = 0;
  common.ordering = 1;
  klu_symbolic *symbolic = klu_analyze(mSize, matrix.outerIndexPtr(),
                                       matrix.innerIndexPtr(), &common);
  if (!symbolic)
    throw std::runtime_error("KLU could not analyze the rocSOLVER seed matrix");
  klu_numeric *numeric =
      klu_factor(matrix.outerIndexPtr(), matrix.innerIndexPtr(),
                 matrix.valuePtr(), symbolic, &common);
  if (!numeric) {
    klu_free_symbolic(&symbolic, &common);
    throw std::runtime_error("KLU could not factor the rocSOLVER seed matrix");
  }

  std::vector<int> lowerColumnPointers(static_cast<std::size_t>(mSize) + 1);
  std::vector<int> lowerRowIndices(static_cast<std::size_t>(numeric->lnz));
  std::vector<double> lowerValues(static_cast<std::size_t>(numeric->lnz));
  std::vector<int> upperColumnPointers(static_cast<std::size_t>(mSize) + 1);
  std::vector<int> upperRowIndices(static_cast<std::size_t>(numeric->unz));
  std::vector<double> upperValues(static_cast<std::size_t>(numeric->unz));
  std::vector<int> rowPivots(static_cast<std::size_t>(mSize));
  std::vector<int> columnPivots(static_cast<std::size_t>(mSize));
  if (!klu_extract(
          numeric, symbolic, lowerColumnPointers.data(), lowerRowIndices.data(),
          lowerValues.data(), upperColumnPointers.data(),
          upperRowIndices.data(), upperValues.data(), nullptr, nullptr, nullptr,
          rowPivots.data(), columnPivots.data(), nullptr, nullptr, &common)) {
    klu_free_numeric(&numeric, &common);
    klu_free_symbolic(&symbolic, &common);
    throw std::runtime_error("KLU could not extract the rocSOLVER seed LU");
  }

  std::vector<std::map<int, double>> factorRows(
      static_cast<std::size_t>(mSize));
  for (int column = 0; column < mSize; ++column)
    for (int entry = lowerColumnPointers[static_cast<std::size_t>(column)];
         entry < lowerColumnPointers[static_cast<std::size_t>(column) + 1];
         ++entry) {
      const int row = lowerRowIndices[static_cast<std::size_t>(entry)];
      if (row != column)
        factorRows[static_cast<std::size_t>(row)][column] =
            lowerValues[static_cast<std::size_t>(entry)];
    }
  for (int column = 0; column < mSize; ++column)
    for (int entry = upperColumnPointers[static_cast<std::size_t>(column)];
         entry < upperColumnPointers[static_cast<std::size_t>(column) + 1];
         ++entry) {
      const int row = upperRowIndices[static_cast<std::size_t>(entry)];
      factorRows[static_cast<std::size_t>(row)][column] =
          upperValues[static_cast<std::size_t>(entry)];
    }

  std::vector<int> factorRowPointers(static_cast<std::size_t>(mSize) + 1, 0);
  for (int row = 0; row < mSize; ++row)
    factorRowPointers[static_cast<std::size_t>(row) + 1] =
        factorRowPointers[static_cast<std::size_t>(row)] +
        static_cast<int>(factorRows[static_cast<std::size_t>(row)].size());
  mFactorNonZeros = factorRowPointers.back();
  if (mFactorNonZeros == 0)
    throw std::runtime_error("rocSOLVER LU seed factor is empty");
  std::vector<int> factorColumnIndices;
  std::vector<double> factorValues;
  factorColumnIndices.reserve(static_cast<std::size_t>(mFactorNonZeros));
  factorValues.reserve(static_cast<std::size_t>(mFactorNonZeros));
  for (const auto &row : factorRows)
    for (const auto &[column, value] : row) {
      factorColumnIndices.push_back(column);
      factorValues.push_back(value);
    }

  std::vector<double> probeRightHandSide(static_cast<std::size_t>(mSize), 1.0);
  std::vector<double> expected(static_cast<std::size_t>(mSize));
  if (!klu_solve(symbolic, numeric, mSize, 1, probeRightHandSide.data(),
                 &common)) {
    klu_free_numeric(&numeric, &common);
    klu_free_symbolic(&symbolic, &common);
    throw std::runtime_error("KLU seed solve failed for rocSOLVER setup");
  }
  expected = probeRightHandSide;
  klu_free_numeric(&numeric, &common);
  klu_free_symbolic(&symbolic, &common);

  checkHip(
      hipMalloc(&mFactorRowPointers, factorRowPointers.size() * sizeof(int)),
      "hipMalloc(initial LU row pointers)");
  checkHip(hipMalloc(&mFactorColumnIndices,
                     factorColumnIndices.size() * sizeof(int)),
           "hipMalloc(initial LU column indices)");
  checkHip(hipMemcpy(mFactorRowPointers, factorRowPointers.data(),
                     factorRowPointers.size() * sizeof(int),
                     hipMemcpyHostToDevice),
           "hipMemcpy(initial LU row pointers)");
  checkHip(hipMemcpy(mFactorColumnIndices, factorColumnIndices.data(),
                     factorColumnIndices.size() * sizeof(int),
                     hipMemcpyHostToDevice),
           "hipMemcpy(initial LU column indices)");
  checkHip(hipMemcpy(mRowPivots, rowPivots.data(), mSize * sizeof(int),
                     hipMemcpyHostToDevice),
           "hipMemcpy(initial LU row pivots)");
  checkHip(hipMemcpy(mColumnPivots, columnPivots.data(), mSize * sizeof(int),
                     hipMemcpyHostToDevice),
           "hipMemcpy(initial LU column pivots)");

  checkHip(hipMalloc(&mFactorValues, factorValues.size() * sizeof(double)),
           "hipMalloc(initial LU values)");
  checkHip(hipMemcpy(mFactorValues, factorValues.data(),
                     factorValues.size() * sizeof(double),
                     hipMemcpyHostToDevice),
           "hipMemcpy(initial LU values)");
  checkSolver(rocsolver_dcsrrf_analysis(
                  mHandle, mSize, 1, mNonZeros, mRowPointers, mColumnIndices,
                  mValues, mFactorNonZeros, mFactorRowPointers,
                  mFactorColumnIndices, mFactorValues, mRowPivots,
                  mColumnPivots, mSolution, mSize, mRfInfo),
              "rocsolver_dcsrrf_analysis");
  mMatrixReady = true;

  std::vector<double> deviceProbeRhs(static_cast<std::size_t>(mSize), 1.0);
  checkHip(hipMemcpy(mRightHandSide, deviceProbeRhs.data(),
                     mSize * sizeof(double), hipMemcpyHostToDevice),
           "hipMemcpy(refactorization probe RHS)");
  solveDevice();
  std::vector<double> deviceProbe(static_cast<std::size_t>(mSize));
  checkHip(hipMemcpy(deviceProbe.data(), mSolution, mSize * sizeof(double),
                     hipMemcpyDeviceToHost),
           "hipMemcpy(refactorization probe solution)");
  double maxDifference = 0.0;
  double maxExpected = 0.0;
  for (int row = 0; row < mSize; ++row) {
    maxDifference = std::max(
        maxDifference, std::abs(deviceProbe[static_cast<std::size_t>(row)] -
                                expected[static_cast<std::size_t>(row)]));
    maxExpected = std::max(maxExpected,
                           std::abs(expected[static_cast<std::size_t>(row)]));
  }
  if (!std::isfinite(maxDifference) ||
      maxDifference > 1e-8 * std::max(1.0, maxExpected))
    throw std::runtime_error(
        "rocSOLVER sparse LU seed solve differs from KLU (max abs diff=" +
        std::to_string(maxDifference) + ")");
}

void HipSparseLURefactorAdapter::factorize(SparseMatrix &systemMatrix) {
  uploadMatrix(systemMatrix, true);
  initializeRefactorization(systemMatrix);
}

void HipSparseLURefactorAdapter::refactorize(SparseMatrix &systemMatrix) {
  uploadMatrix(systemMatrix, false);
  refactorizeDevice();
}

void HipSparseLURefactorAdapter::partialRefactorize(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  (void)listVariableSystemMatrixEntries;
  refactorize(systemMatrix);
}

void HipSparseLURefactorAdapter::clearDeviceRightHandSide() {
  if (!mRightHandSide || mSize <= 0)
    throw std::runtime_error("rocSOLVER sparse LU RHS is not initialized");
  checkHip(hipMemsetAsync(mRightHandSide, 0,
                          static_cast<std::size_t>(mSize) * sizeof(double)),
           "hipMemsetAsync(device RHS)");
}

void HipSparseLURefactorAdapter::uploadRightHandSide(
    const Matrix &rightSideVector) {
  if (rightSideVector.size() != mSize)
    throw std::runtime_error("rocSOLVER sparse LU RHS size mismatch");
  checkHip(hipMemcpy(mRightHandSide, rightSideVector.data(),
                     mSize * sizeof(double), hipMemcpyHostToDevice),
           "hipMemcpy(right-hand side)");
}

void HipSparseLURefactorAdapter::copyDeviceSolutionToHost(Matrix &solution) {
  if (solution.size() != mSize)
    solution.resize(mSize, 1);
  checkHip(hipMemcpy(solution.data(), mSolution, mSize * sizeof(double),
                     hipMemcpyDeviceToHost),
           "hipMemcpy(solution)");
}

void HipSparseLURefactorAdapter::refactorizeDevice() {
  if (!mMatrixReady)
    throw std::runtime_error("rocSOLVER sparse LU matrix is not initialized");
  checkSolver(rocsolver_dcsrrf_refactlu(mHandle, mSize, mNonZeros, mRowPointers,
                                        mColumnIndices, mValues,
                                        mFactorNonZeros, mFactorRowPointers,
                                        mFactorColumnIndices, mFactorValues,
                                        mRowPivots, mColumnPivots, mRfInfo),
              "rocsolver_dcsrrf_refactlu");
}

void HipSparseLURefactorAdapter::solveDevice(bool verifyConvergence) {
  if (!mMatrixReady)
    throw std::runtime_error("rocSOLVER sparse LU matrix is not initialized");
  checkHip(hipMemcpyAsync(mSolution, mRightHandSide,
                          static_cast<std::size_t>(mSize) * sizeof(double),
                          hipMemcpyDeviceToDevice),
           "hipMemcpyAsync(device RHS to solution)");
  checkSolver(rocsolver_dcsrrf_solve(mHandle, mSize, 1, mFactorNonZeros,
                                     mFactorRowPointers, mFactorColumnIndices,
                                     mFactorValues, mRowPivots, mColumnPivots,
                                     mSolution, mSize, mRfInfo),
              "rocsolver_dcsrrf_solve");
  if (verifyConvergence) {
    std::vector<double> values(static_cast<std::size_t>(mNonZeros));
    std::vector<double> rhs(static_cast<std::size_t>(mSize));
    std::vector<double> solution(static_cast<std::size_t>(mSize));
    checkHip(hipMemcpy(values.data(), mValues,
                       static_cast<std::size_t>(mNonZeros) * sizeof(double),
                       hipMemcpyDeviceToHost),
             "hipMemcpy(residual check matrix values)");
    checkHip(hipMemcpy(rhs.data(), mRightHandSide,
                       static_cast<std::size_t>(mSize) * sizeof(double),
                       hipMemcpyDeviceToHost),
             "hipMemcpy(residual check RHS)");
    checkHip(hipMemcpy(solution.data(), mSolution,
                       static_cast<std::size_t>(mSize) * sizeof(double),
                       hipMemcpyDeviceToHost),
             "hipMemcpy(residual check solution)");

    double residualNorm = 0.0;
    double matrixNorm = 0.0;
    double solutionNorm = 0.0;
    double rhsNorm = 0.0;
    for (int row = 0; row < mSize; ++row) {
      double product = 0.0;
      double rowNorm = 0.0;
      for (int entry = mHostRowPointers[static_cast<std::size_t>(row)];
           entry < mHostRowPointers[static_cast<std::size_t>(row) + 1];
           ++entry) {
        const double value = values[static_cast<std::size_t>(entry)];
        product +=
            value * solution[static_cast<std::size_t>(
                        mHostColumnIndices[static_cast<std::size_t>(entry)])];
        rowNorm += std::abs(value);
      }
      residualNorm = std::max(
          residualNorm, std::abs(product - rhs[static_cast<std::size_t>(row)]));
      matrixNorm = std::max(matrixNorm, rowNorm);
      solutionNorm = std::max(
          solutionNorm, std::abs(solution[static_cast<std::size_t>(row)]));
      rhsNorm = std::max(rhsNorm, std::abs(rhs[static_cast<std::size_t>(row)]));
    }
    const double denominator = matrixNorm * solutionNorm + rhsNorm;
    const double relativeResidual =
        residualNorm /
        std::max(denominator, std::numeric_limits<double>::min());
    ++mVerifiedSolveCount;
    if (!std::isfinite(relativeResidual) || relativeResidual > 1e-8)
      throw std::runtime_error(
          "rocSOLVER sparse LU residual check failed at verified solve " +
          std::to_string(mVerifiedSolveCount) +
          " (relative infinity norm=" + std::to_string(relativeResidual) + ")");
  }
}

Matrix HipSparseLURefactorAdapter::solve(Matrix &rightSideVector) {
  if (!mMatrixReady || rightSideVector.size() != mSize)
    throw std::runtime_error("rocSOLVER sparse LU RHS size mismatch");
  uploadRightHandSide(rightSideVector);
  solveDevice();
  Matrix solution(mSize, 1);
  copyDeviceSolutionToHost(solution);
  return solution;
}
} // namespace DPsim
