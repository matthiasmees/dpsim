// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCReactiveOuterController.h>

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

MMCReactiveOuterController::MMCReactiveOuterController()
    : mMode(MMCReactiveMode::OpenLoop), mKp(0.0), mKi(0.0),
      mReactivePowerReference(0.0), mAcVoltageReference(0.0),
      mOpenLoopCurrentReference(0.0) {}
void MMCReactiveOuterController::setPI(Real kp, Real ki) {
  requireNonNegative(kp, "Reactive controller Kp");
  requireNonNegative(ki, "Reactive controller Ki");
  mKp = kp;
  mKi = ki;
}
void MMCReactiveOuterController::setMode(MMCReactiveMode mode) { mMode = mode; }
void MMCReactiveOuterController::setReferences(Real reactivePowerReference,
                                               Real acVoltageReference,
                                               Real openLoopCurrentReference) {
  if (!std::isfinite(reactivePowerReference) ||
      !std::isfinite(acVoltageReference) ||
      !std::isfinite(openLoopCurrentReference))
    throw std::invalid_argument(
        "Reactive-controller references must be finite.");
  mReactivePowerReference = reactivePowerReference;
  mAcVoltageReference = acVoltageReference;
  mOpenLoopCurrentReference = openLoopCurrentReference;
}
UInt MMCReactiveOuterController::stateSize() const { return 1; }
UInt MMCReactiveOuterController::inputSize() const { return 2; }
UInt MMCReactiveOuterController::outputSize() const { return 3; }
std::vector<String> MMCReactiveOuterController::stateNames() const {
  return {"xi_reactive"};
}
void MMCReactiveOuterController::evaluate(const Matrix &x, const Matrix &u,
                                          Matrix &dx, Matrix &y) const {
  validateDimensions(x, u);
  const Real q = u(0, 0);
  const Real vac = u(1, 0);
  Real error = 0.0;
  Real output = mOpenLoopCurrentReference;
  Real integratorSign = 1.0;
  Bool integrate = false;

  switch (mMode) {
  case MMCReactiveMode::ReactivePower:
    error = mReactivePowerReference - q;
    output = -(mKp * error + mKi * x(0, 0));
    integratorSign = -1.0;
    integrate = true;
    break;
  case MMCReactiveMode::AcVoltage:
    error = mAcVoltageReference - vac;
    output = mKp * error + mKi * x(0, 0);
    integrate = true;
    break;
  case MMCReactiveMode::OpenLoop:
    break;
  }

  dx = Matrix::Zero(1, 1);
  if (integrate)
    dx(0, 0) = error;
  y = Matrix::Zero(3, 1);
  y(0, 0) = output;
  y(1, 0) = error;
  y(2, 0) = integratorSign;
}

void MMCReactiveOuterController::evaluateStateDerivative(const Matrix &x,
                                                         const Matrix &u,
                                                         Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCReactiveOuterController::evaluateOutput(const Matrix &x,
                                                const Matrix &u,
                                                Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCReactiveOuterController::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCReactiveOuterController::buildStateSpaceModel(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Matrix &E, Matrix &F, Real relativeStep,
    Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCReactiveOuterController::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                               Real relativeStep,
                                               Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCReactiveOuterController::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
