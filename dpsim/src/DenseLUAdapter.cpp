/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/DenseLUAdapter.h>

#include <algorithm>

using namespace DPsim;

namespace DPsim {

DenseLUAdapter::~DenseLUAdapter() = default;

void DenseLUAdapter::applyConfiguration() {
  mUseStaticPreconditioner = mConfiguration.getStaticPreconditioner();
}

void DenseLUAdapter::preprocessing(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  (void)systemMatrix;
  mVariableEntries = listVariableSystemMatrixEntries;
  std::sort(mVariableEntries.begin(), mVariableEntries.end());
  mVariableEntries.erase(
      std::unique(mVariableEntries.begin(), mVariableEntries.end()),
      mVariableEntries.end());
}

void DenseLUAdapter::factorize(SparseMatrix &systemMatrix) {
  mReferenceMatrix = Matrix(systemMatrix);
  mCurrentMatrix = mReferenceMatrix;
  mResidual.resize(mReferenceMatrix.rows(), 1);
  mCorrection.resize(mReferenceMatrix.rows(), 1);
  LUFactorized.compute(mReferenceMatrix);
  mStaticPreconditionerActive =
      mUseStaticPreconditioner && !mVariableEntries.empty();
}

void DenseLUAdapter::refactorize(SparseMatrix &systemMatrix) {
  if (!mStaticPreconditionerActive) {
    LUFactorized.compute(Matrix(systemMatrix));
    return;
  }
  // Only declared variable entries can differ from the reference matrix.
  // Keep all static dense coefficients in place instead of copying the whole
  // reduced matrix before updating the generator stamps.
  for (const auto &[row, column] : mVariableEntries)
    mCurrentMatrix(row, column) = systemMatrix.coeff(row, column);
}

void DenseLUAdapter::partialRefactorize(
    SparseMatrix &systemMatrix,
    std::vector<std::pair<UInt, UInt>> &listVariableSystemMatrixEntries) {
  if (!mStaticPreconditionerActive) {
    LUFactorized.compute(Matrix(systemMatrix));
    return;
  }
  mVariableEntries = listVariableSystemMatrixEntries;
  std::sort(mVariableEntries.begin(), mVariableEntries.end());
  mVariableEntries.erase(
      std::unique(mVariableEntries.begin(), mVariableEntries.end()),
      mVariableEntries.end());
  for (const auto &[row, column] : mVariableEntries)
    mCurrentMatrix(row, column) = systemMatrix.coeff(row, column);
}

Matrix DenseLUAdapter::solve(Matrix &mRightHandSideVector) {
  Matrix solution;
  solveInto(mRightHandSideVector, solution);
  return solution;
}

void DenseLUAdapter::solveInto(Matrix &rightSideVector, Matrix &solution) {
  if (!mStaticPreconditionerActive) {
    solution = LUFactorized.solve(rightSideVector);
    return;
  }

  solution = LUFactorized.solve(rightSideVector);

  const Real rhsScale =
      std::max<Real>(1., rightSideVector.cwiseAbs().maxCoeff());
  const Real tolerance = 1e-12 * rhsScale;
  mResidual.setZero();
  for (const auto &[row, column] : mVariableEntries) {
    const Real delta =
        mCurrentMatrix(row, column) - mReferenceMatrix(row, column);
    mResidual(row) -= delta * solution(column);
  }

  for (UInt iteration = 0; iteration < 12; ++iteration) {
    if (mResidual.cwiseAbs().maxCoeff() <= tolerance)
      return;
    mCorrection = LUFactorized.solve(mResidual);
    solution += mCorrection;
    mResidual.setZero();
    for (const auto &[row, column] : mVariableEntries) {
      const Real delta =
          mCurrentMatrix(row, column) - mReferenceMatrix(row, column);
      mResidual(row) -= delta * mCorrection(column);
    }
  }

  // Refinement is an optimization only. Fall back to an exact solve if the
  // previous LU is not an adequate preconditioner for the new matrix.
  LUFactorized.compute(mCurrentMatrix);
  mReferenceMatrix = mCurrentMatrix;
  solution = LUFactorized.solve(rightSideVector);
}
} // namespace DPsim
