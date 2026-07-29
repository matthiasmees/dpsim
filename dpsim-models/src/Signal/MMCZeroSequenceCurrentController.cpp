// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCZeroSequenceCurrentController.h>

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

MMCZeroSequenceCurrentController::MMCZeroSequenceCurrentController()
    : mKp(0.0), mKi(0.0) {}
void MMCZeroSequenceCurrentController::setPI(Real kp, Real ki) {
  requireNonNegative(kp, "Zero-sequence current controller Kp");
  requireNonNegative(ki, "Zero-sequence current controller Ki");
  mKp = kp;
  mKi = ki;
}
UInt MMCZeroSequenceCurrentController::stateSize() const { return 1; }
UInt MMCZeroSequenceCurrentController::inputSize() const { return 3; }
UInt MMCZeroSequenceCurrentController::outputSize() const { return 2; }
std::vector<String> MMCZeroSequenceCurrentController::stateNames() const {
  return {"xi_zcc"};
}
void MMCZeroSequenceCurrentController::evaluate(const Matrix &x,
                                                const Matrix &u, Matrix &dx,
                                                Matrix &y) const {
  validateDimensions(x, u);
  const Real error = u(0, 0) - u(1, 0);
  dx = Matrix::Constant(1, 1, error);
  y = Matrix::Zero(2, 1);
  y(0, 0) = -(mKp * error + mKi * x(0, 0) - u(2, 0) / 2.0);
  y(1, 0) = error;
}

void MMCZeroSequenceCurrentController::evaluateStateDerivative(
    const Matrix &x, const Matrix &u, Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCZeroSequenceCurrentController::evaluateOutput(const Matrix &x,
                                                      const Matrix &u,
                                                      Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCZeroSequenceCurrentController::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCZeroSequenceCurrentController::buildStateSpaceModel(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Matrix &E, Matrix &F, Real relativeStep,
    Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization MMCZeroSequenceCurrentController::getStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep,
    Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization
MMCZeroSequenceCurrentController::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
