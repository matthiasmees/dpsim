// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCOutputCurrentController.h>

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

MMCOutputCurrentController::MMCOutputCurrentController() : mKp(0.0), mKi(0.0) {}
void MMCOutputCurrentController::setPI(Real kp, Real ki) {
  requireNonNegative(kp, "Output-current controller Kp");
  requireNonNegative(ki, "Output-current controller Ki");
  mKp = kp;
  mKi = ki;
}
UInt MMCOutputCurrentController::stateSize() const { return 2; }
UInt MMCOutputCurrentController::inputSize() const { return 9; }
UInt MMCOutputCurrentController::outputSize() const { return 4; }
std::vector<String> MMCOutputCurrentController::stateNames() const {
  return {"xi_occ_d", "xi_occ_q"};
}
void MMCOutputCurrentController::evaluate(const Matrix &x, const Matrix &u,
                                          Matrix &dx, Matrix &y) const {
  validateDimensions(x, u);
  const Real idRef = u(0, 0);
  const Real iqRef = u(1, 0);
  const Real id = u(2, 0);
  const Real iq = u(3, 0);
  const Real vd = u(4, 0);
  const Real vq = u(5, 0);
  const Real omega = u(6, 0);
  const Real resistance = u(7, 0);
  const Real inductance = u(8, 0);
  const Real ed = idRef - id;
  const Real eq = iqRef - iq;

  dx = Matrix::Zero(2, 1);
  dx(0, 0) = ed;
  dx(1, 0) = eq;
  y = Matrix::Zero(4, 1);
  y(0, 0) = vd + mKp * ed + mKi * x(0, 0) + resistance * idRef -
            omega * inductance * iqRef;
  y(1, 0) = vq + mKp * eq + mKi * x(1, 0) + resistance * iqRef +
            omega * inductance * idRef;
  y(2, 0) = ed;
  y(3, 0) = eq;
}

void MMCOutputCurrentController::evaluateStateDerivative(const Matrix &x,
                                                         const Matrix &u,
                                                         Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCOutputCurrentController::evaluateOutput(const Matrix &x,
                                                const Matrix &u,
                                                Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCOutputCurrentController::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCOutputCurrentController::buildStateSpaceModel(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Matrix &E, Matrix &F, Real relativeStep,
    Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCOutputCurrentController::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                               Real relativeStep,
                                               Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCOutputCurrentController::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
