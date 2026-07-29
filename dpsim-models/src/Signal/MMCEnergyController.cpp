// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCEnergyController.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;
using namespace CPS::Signal;

namespace {
Real regularizeSigned(Real value, Real minimumMagnitude) {
  if (std::abs(value) >= minimumMagnitude)
    return value;
  return value >= 0.0 ? minimumMagnitude : -minimumMagnitude;
}

void requireNonNegative(Real value, const char *name) {
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument(String(name) +
                                " must be finite and non-negative.");
}
} // namespace

MMCEnergyController::MMCEnergyController()
    : mKp(0.0), mKi(0.0), mEnabled(false), mEnergyReference(0.0),
      mMinimumDcVoltage(1.0) {}
void MMCEnergyController::setParameters(Real kp, Real ki, Bool enabled,
                                        Real energyReference,
                                        Real minimumDcVoltage) {
  requireNonNegative(kp, "Energy controller Kp");
  requireNonNegative(ki, "Energy controller Ki");
  if (!std::isfinite(energyReference) || energyReference < 0.0)
    throw std::invalid_argument(
        "Energy reference must be finite and non-negative.");
  if (!std::isfinite(minimumDcVoltage) || minimumDcVoltage <= 0.0)
    throw std::invalid_argument(
        "Minimum DC voltage must be finite and positive.");
  mKp = kp;
  mKi = ki;
  mEnabled = enabled;
  mEnergyReference = energyReference;
  mMinimumDcVoltage = minimumDcVoltage;
}
UInt MMCEnergyController::stateSize() const { return 1; }
UInt MMCEnergyController::inputSize() const { return 3; }
UInt MMCEnergyController::outputSize() const { return 2; }
std::vector<String> MMCEnergyController::stateNames() const {
  return {"xi_energy"};
}
void MMCEnergyController::evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                                   Matrix &y) const {
  validateDimensions(x, u);
  dx = Matrix::Zero(1, 1);
  y = Matrix::Zero(2, 1);
  if (!mEnabled)
    return;
  const Real error = mEnergyReference - u(0, 0);
  const Real powerCommand = mKp * error + mKi * x(0, 0) + u(1, 0);
  const Real vdc = regularizeSigned(u(2, 0), mMinimumDcVoltage);
  dx(0, 0) = error;
  y(0, 0) = powerCommand / (3.0 * vdc);
  y(1, 0) = error;
}

void MMCEnergyController::evaluateStateDerivative(const Matrix &x,
                                                  const Matrix &u,
                                                  Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCEnergyController::evaluateOutput(const Matrix &x, const Matrix &u,
                                         Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCEnergyController::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCEnergyController::buildStateSpaceModel(const Matrix &x, const Matrix &u,
                                               Matrix &A, Matrix &B, Matrix &C,
                                               Matrix &D, Matrix &E, Matrix &F,
                                               Real relativeStep,
                                               Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCEnergyController::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                        Real relativeStep,
                                        Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCEnergyController::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
