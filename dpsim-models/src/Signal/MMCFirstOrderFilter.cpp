// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCFirstOrderFilter.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;
using namespace CPS::Signal;

namespace {
void requireNonNegative(Real value, const char *name) {
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument(String(name) +
                                " must be finite and non-negative.");
}
} // namespace

MMCFirstOrderFilter::MMCFirstOrderFilter(UInt channels)
    : mChannels(channels), mTimeConstant(0.0) {
  if (channels == 0)
    throw std::invalid_argument(
        "First-order filter needs at least one channel.");
}

void MMCFirstOrderFilter::setTimeConstant(Real timeConstant) {
  requireNonNegative(timeConstant, "First-order filter time constant");
  mTimeConstant = timeConstant;
}
UInt MMCFirstOrderFilter::stateSize() const { return mChannels; }
UInt MMCFirstOrderFilter::inputSize() const { return mChannels; }
UInt MMCFirstOrderFilter::outputSize() const { return mChannels; }
std::vector<String> MMCFirstOrderFilter::stateNames() const {
  std::vector<String> result;
  for (UInt i = 0; i < mChannels; ++i)
    result.push_back("filter_1_" + std::to_string(i));
  return result;
}
void MMCFirstOrderFilter::evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                                   Matrix &y) const {
  validateDimensions(x, u);
  dx.resize(stateSize(), 1);
  y.resize(outputSize(), 1);
  if (mTimeConstant <= 0.0) {
    dx.setZero();
    y = u;
    return;
  }
  dx = (u - x) / mTimeConstant;
  y = x;
}

void MMCFirstOrderFilter::evaluateStateDerivative(const Matrix &x,
                                                  const Matrix &u,
                                                  Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCFirstOrderFilter::evaluateOutput(const Matrix &x, const Matrix &u,
                                         Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCFirstOrderFilter::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCFirstOrderFilter::buildStateSpaceModel(const Matrix &x, const Matrix &u,
                                               Matrix &A, Matrix &B, Matrix &C,
                                               Matrix &D, Matrix &E, Matrix &F,
                                               Real relativeStep,
                                               Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCFirstOrderFilter::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                        Real relativeStep,
                                        Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCFirstOrderFilter::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
