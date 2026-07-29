// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCPadeDelay.h>

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

MMCPadeDelay::MMCPadeDelay(UInt channels)
    : mChannels(channels), mTimeDelay(0.0), mEnabled(false) {
  if (channels == 0)
    throw std::invalid_argument("Padé delay needs at least one channel.");
}
void MMCPadeDelay::setDelay(Real timeDelay, Bool enabled) {
  requireNonNegative(timeDelay, "Padé delay");
  mTimeDelay = timeDelay;
  mEnabled = enabled && timeDelay > 0.0;
}
UInt MMCPadeDelay::stateSize() const { return 2 * mChannels; }
UInt MMCPadeDelay::inputSize() const { return mChannels; }
UInt MMCPadeDelay::outputSize() const { return mChannels; }
std::vector<String> MMCPadeDelay::stateNames() const {
  std::vector<String> result;
  for (UInt i = 0; i < mChannels; ++i) {
    result.push_back("delay_" + std::to_string(i) + "_1");
    result.push_back("delay_" + std::to_string(i) + "_2");
  }
  return result;
}
void MMCPadeDelay::evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                            Matrix &y) const {
  validateDimensions(x, u);
  dx = Matrix::Zero(stateSize(), 1);
  y = Matrix::Zero(outputSize(), 1);
  if (!mEnabled) {
    y = u;
    return;
  }
  const Real a0 = 12.0 / (mTimeDelay * mTimeDelay);
  const Real a1 = 6.0 / mTimeDelay;
  for (UInt channel = 0; channel < mChannels; ++channel) {
    const UInt x1 = 2 * channel;
    const UInt x2 = x1 + 1;
    dx(x1, 0) = x(x2, 0);
    dx(x2, 0) = -a0 * x(x1, 0) - a1 * x(x2, 0) + u(channel, 0);
    y(channel, 0) = u(channel, 0) - 2.0 * a1 * x(x2, 0);
  }
}

void MMCPadeDelay::evaluateStateDerivative(const Matrix &x, const Matrix &u,
                                           Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCPadeDelay::evaluateOutput(const Matrix &x, const Matrix &u,
                                  Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCPadeDelay::calculateNumericalJacobians(const Matrix &x, const Matrix &u,
                                               Matrix &A, Matrix &B, Matrix &C,
                                               Matrix &D, Real relativeStep,
                                               Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCPadeDelay::buildStateSpaceModel(const Matrix &x, const Matrix &u,
                                        Matrix &A, Matrix &B, Matrix &C,
                                        Matrix &D, Matrix &E, Matrix &F,
                                        Real relativeStep,
                                        Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization MMCPadeDelay::getStateSpaceModel(const Matrix &x,
                                                  const Matrix &u,
                                                  Real relativeStep,
                                                  Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization
MMCPadeDelay::getSparseStateSpaceModel(const Matrix &x, const Matrix &u,
                                       Real relativeStep, Real absoluteStep,
                                       Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
