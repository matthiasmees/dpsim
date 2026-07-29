// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCActiveOuterController.h>

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

MMCActiveOuterController::MMCActiveOuterController()
    : mMode(MMCActiveMode::OpenLoop), mKp(0.0), mKi(0.0),
      mActivePowerReference(0.0), mDcVoltageReference(0.0), mDroopGain(1.0),
      mOpenLoopCurrentReference(0.0) {}
void MMCActiveOuterController::setPI(Real kp, Real ki) {
  requireNonNegative(kp, "Active controller Kp");
  requireNonNegative(ki, "Active controller Ki");
  mKp = kp;
  mKi = ki;
}
void MMCActiveOuterController::setMode(MMCActiveMode mode) { mMode = mode; }
void MMCActiveOuterController::setReferences(Real activePowerReference,
                                             Real dcVoltageReference,
                                             Real droopGain,
                                             Real openLoopCurrentReference) {
  if (!std::isfinite(activePowerReference) ||
      !std::isfinite(dcVoltageReference) || !std::isfinite(droopGain) ||
      !std::isfinite(openLoopCurrentReference))
    throw std::invalid_argument("Active-controller references must be finite.");
  if (mMode == MMCActiveMode::DcDroop && droopGain == 0.0)
    throw std::invalid_argument("DC-droop gain must be non-zero.");
  mActivePowerReference = activePowerReference;
  mDcVoltageReference = dcVoltageReference;
  mDroopGain = droopGain;
  mOpenLoopCurrentReference = openLoopCurrentReference;
}
UInt MMCActiveOuterController::stateSize() const { return 1; }
UInt MMCActiveOuterController::inputSize() const { return 4; }
UInt MMCActiveOuterController::outputSize() const { return 3; }
std::vector<String> MMCActiveOuterController::stateNames() const {
  return {"xi_active"};
}
void MMCActiveOuterController::evaluate(const Matrix &x, const Matrix &u,
                                        Matrix &dx, Matrix &y) const {
  validateDimensions(x, u);
  const Real p = u(0, 0);
  const Real vdc = u(1, 0);
  const Real vd = u(2, 0);
  const Real held = u(3, 0);

  Real error = 0.0;
  Real output = mOpenLoopCurrentReference;
  Real integratorSign = 1.0;
  Bool integrate = false;

  switch (mMode) {
  case MMCActiveMode::ActivePower:
    error = mActivePowerReference - p;
    output = mKp * error + mKi * x(0, 0);
    integrate = true;
    break;
  case MMCActiveMode::DcVoltage:
    error = vdc - mDcVoltageReference;
    output = mKp * error + mKi * x(0, 0);
    integrate = true;
    break;
  case MMCActiveMode::DcDroop:
    output =
        (mActivePowerReference + (mDcVoltageReference - vdc) / mDroopGain) /
        std::max(1.0, std::abs(vd));
    break;
  case MMCActiveMode::SampledPowerFeedforward:
    output = held;
    break;
  case MMCActiveMode::OpenLoop:
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

void MMCActiveOuterController::evaluateStateDerivative(const Matrix &x,
                                                       const Matrix &u,
                                                       Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCActiveOuterController::evaluateOutput(const Matrix &x, const Matrix &u,
                                              Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCActiveOuterController::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCActiveOuterController::buildStateSpaceModel(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Matrix &E, Matrix &F, Real relativeStep,
    Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCActiveOuterController::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                             Real relativeStep,
                                             Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCActiveOuterController::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
