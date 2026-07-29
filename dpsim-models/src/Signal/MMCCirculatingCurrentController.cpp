// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCCirculatingCurrentController.h>

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

MMCCirculatingCurrentController::MMCCirculatingCurrentController()
    : mKp(0.0), mKi(0.0) {}
void MMCCirculatingCurrentController::setPI(Real kp, Real ki) {
  requireNonNegative(kp, "Circulating-current controller Kp");
  requireNonNegative(ki, "Circulating-current controller Ki");
  mKp = kp;
  mKi = ki;
}
UInt MMCCirculatingCurrentController::stateSize() const { return 2; }
UInt MMCCirculatingCurrentController::inputSize() const { return 6; }
UInt MMCCirculatingCurrentController::outputSize() const { return 4; }
std::vector<String> MMCCirculatingCurrentController::stateNames() const {
  return {"xi_ccc_d", "xi_ccc_q"};
}
void MMCCirculatingCurrentController::evaluate(const Matrix &x, const Matrix &u,
                                               Matrix &dx, Matrix &y) const {
  validateDimensions(x, u);
  const Real idRef = u(0, 0);
  const Real iqRef = u(1, 0);
  const Real id = u(2, 0);
  const Real iq = u(3, 0);
  const Real omega = u(4, 0);
  const Real lArm = u(5, 0);
  const Real ed = idRef - id;
  const Real eq = iqRef - iq;
  dx = Matrix::Zero(2, 1);
  dx(0, 0) = ed;
  dx(1, 0) = eq;
  y = Matrix::Zero(4, 1);
  y(0, 0) = -(mKp * ed + mKi * x(0, 0) + 2.0 * omega * lArm * iq);
  y(1, 0) = -(mKp * eq + mKi * x(1, 0) - 2.0 * omega * lArm * id);
  y(2, 0) = ed;
  y(3, 0) = eq;
}

void MMCCirculatingCurrentController::evaluateStateDerivative(
    const Matrix &x, const Matrix &u, Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCCirculatingCurrentController::evaluateOutput(const Matrix &x,
                                                     const Matrix &u,
                                                     Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCCirculatingCurrentController::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCCirculatingCurrentController::buildStateSpaceModel(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Matrix &E, Matrix &F, Real relativeStep,
    Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization MMCCirculatingCurrentController::getStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep,
    Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization
MMCCirculatingCurrentController::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
