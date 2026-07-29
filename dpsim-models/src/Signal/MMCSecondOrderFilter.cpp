// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCSecondOrderFilter.h>

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

MMCSecondOrderFilter::MMCSecondOrderFilter(UInt channels)
    : mChannels(channels), mTimeConstant(0.0) {
  if (channels == 0)
    throw std::invalid_argument(
        "Second-order filter needs at least one channel.");
}
void MMCSecondOrderFilter::setTimeConstant(Real timeConstant) {
  requireNonNegative(timeConstant, "Second-order filter time constant");
  mTimeConstant = timeConstant;
}
UInt MMCSecondOrderFilter::stateSize() const { return 2 * mChannels; }
UInt MMCSecondOrderFilter::inputSize() const { return mChannels; }
UInt MMCSecondOrderFilter::outputSize() const { return mChannels; }
std::vector<String> MMCSecondOrderFilter::stateNames() const {
  std::vector<String> result;
  for (UInt i = 0; i < mChannels; ++i) {
    result.push_back("filter_1_" + std::to_string(i));
    result.push_back("filter_2_" + std::to_string(i));
  }
  return result;
}
void MMCSecondOrderFilter::evaluate(const Matrix &x, const Matrix &u,
                                    Matrix &dx, Matrix &y) const {
  validateDimensions(x, u);
  dx = Matrix::Zero(stateSize(), 1);
  y = Matrix::Zero(outputSize(), 1);
  if (mTimeConstant <= 0.0) {
    y = u;
    return;
  }
  for (UInt channel = 0; channel < mChannels; ++channel) {
    const UInt first = 2 * channel;
    const UInt second = first + 1;
    dx(first, 0) = (u(channel, 0) - x(first, 0)) / mTimeConstant;
    dx(second, 0) = (x(first, 0) - x(second, 0)) / mTimeConstant;
    y(channel, 0) = x(second, 0);
  }
}

void MMCSecondOrderFilter::evaluateStateDerivative(const Matrix &x,
                                                   const Matrix &u,
                                                   Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCSecondOrderFilter::evaluateOutput(const Matrix &x, const Matrix &u,
                                          Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCSecondOrderFilter::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCSecondOrderFilter::buildStateSpaceModel(const Matrix &x,
                                                const Matrix &u, Matrix &A,
                                                Matrix &B, Matrix &C, Matrix &D,
                                                Matrix &E, Matrix &F,
                                                Real relativeStep,
                                                Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCSecondOrderFilter::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                         Real relativeStep,
                                         Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCSecondOrderFilter::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
