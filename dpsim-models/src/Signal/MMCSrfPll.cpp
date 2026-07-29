// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCSrfPll.h>

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

MMCSrfPll::MMCSrfPll() : mKp(0.0), mKi(0.0), mEnabled(false) {}
void MMCSrfPll::setParameters(Real kp, Real ki, Bool enabled) {
  requireNonNegative(kp, "PLL Kp");
  requireNonNegative(ki, "PLL Ki");
  mKp = kp;
  mKi = ki;
  mEnabled = enabled;
}
UInt MMCSrfPll::stateSize() const { return 2; }
UInt MMCSrfPll::inputSize() const { return 1; }
UInt MMCSrfPll::outputSize() const { return 2; }
std::vector<String> MMCSrfPll::stateNames() const {
  return {"xi_pll", "pll_angle_deviation"};
}
void MMCSrfPll::evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                         Matrix &y) const {
  validateDimensions(x, u);
  dx = Matrix::Zero(2, 1);
  y = Matrix::Zero(2, 1);
  if (!mEnabled) {
    y(1, 0) = x(1, 0);
    return;
  }
  const Real error = u(0, 0);
  const Real deltaOmega = mKp * error + mKi * x(0, 0);
  dx(0, 0) = error;
  dx(1, 0) = deltaOmega;
  y(0, 0) = deltaOmega;
  y(1, 0) = x(1, 0);
}

void MMCSrfPll::evaluateStateDerivative(const Matrix &x, const Matrix &u,
                                        Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCSrfPll::evaluateOutput(const Matrix &x, const Matrix &u,
                               Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCSrfPll::calculateNumericalJacobians(const Matrix &x, const Matrix &u,
                                            Matrix &A, Matrix &B, Matrix &C,
                                            Matrix &D, Real relativeStep,
                                            Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCSrfPll::buildStateSpaceModel(const Matrix &x, const Matrix &u,
                                     Matrix &A, Matrix &B, Matrix &C, Matrix &D,
                                     Matrix &E, Matrix &F, Real relativeStep,
                                     Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization MMCSrfPll::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                               Real relativeStep,
                                               Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization
MMCSrfPll::getSparseStateSpaceModel(const Matrix &x, const Matrix &u,
                                    Real relativeStep, Real absoluteStep,
                                    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
