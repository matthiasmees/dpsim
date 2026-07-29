// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCSampledPowerFeedforward.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;
using namespace CPS::Signal;

MMCSampledPowerFeedforward::MMCSampledPowerFeedforward()
    : mCutoffFrequency(1000.0), mSampleTime(40e-6), mMinimumDaxisVoltage(1.0) {}
void MMCSampledPowerFeedforward::setParameters(Real cutoffFrequency,
                                               Real sampleTime,
                                               Real minimumDaxisVoltage) {
  if (!std::isfinite(cutoffFrequency) || cutoffFrequency <= 0.0)
    throw std::invalid_argument(
        "Feedforward cutoff frequency must be positive.");
  if (!std::isfinite(sampleTime) || sampleTime <= 0.0)
    throw std::invalid_argument("Feedforward sample time must be positive.");
  if (!std::isfinite(minimumDaxisVoltage) || minimumDaxisVoltage < 0.0)
    throw std::invalid_argument(
        "Feedforward minimum d-axis voltage must be non-negative.");
  mCutoffFrequency = cutoffFrequency;
  mSampleTime = sampleTime;
  mMinimumDaxisVoltage = minimumDaxisVoltage;
}
UInt MMCSampledPowerFeedforward::stateSize() const { return 3; }
UInt MMCSampledPowerFeedforward::inputSize() const { return 2; }
UInt MMCSampledPowerFeedforward::outputSize() const { return 2; }
std::vector<String> MMCSampledPowerFeedforward::stateNames() const {
  return {"feedforward_vd_filtered", "feedforward_vd_derivative",
          "feedforward_id_held"};
}
void MMCSampledPowerFeedforward::evaluate(const Matrix &x, const Matrix &u,
                                          Matrix &dx, Matrix &y) const {
  validateDimensions(x, u);
  // This is a hybrid block. Its states only change at sample instants through
  // sample(); the continuous-time derivative is zero between samples.
  dx = Matrix::Zero(3, 1);
  y = Matrix::Zero(2, 1);
  y(0, 0) = x(0, 0);
  y(1, 0) = x(2, 0);
}
Matrix MMCSampledPowerFeedforward::initialize(Real measuredVd,
                                              Real activePowerReference) const {
  if (!std::isfinite(measuredVd) || measuredVd <= mMinimumDaxisVoltage)
    throw std::invalid_argument(
        "Feedforward initialization needs a valid positive Vd.");
  if (!std::isfinite(activePowerReference))
    throw std::invalid_argument(
        "Feedforward active-power reference must be finite.");
  Matrix x(3, 1);
  x << measuredVd, 0.0, (2.0 / 3.0) * activePowerReference / measuredVd;
  return x;
}
Matrix MMCSampledPowerFeedforward::sample(const Matrix &x, Real measuredVd,
                                          Real activePowerReference) const {
  return sample(x, measuredVd, activePowerReference, mSampleTime);
}

Matrix MMCSampledPowerFeedforward::sample(const Matrix &x, Real measuredVd,
                                          Real activePowerReference,
                                          Real effectiveSampleTime) const {
  if (x.rows() != 3 || x.cols() != 1 || !x.allFinite())
    throw std::invalid_argument("Feedforward state must be finite 3x1.");
  if (!std::isfinite(measuredVd))
    throw std::invalid_argument("Feedforward measured Vd must be finite.");
  if (!std::isfinite(activePowerReference))
    throw std::invalid_argument(
        "Feedforward active-power reference must be finite.");
  if (!std::isfinite(effectiveSampleTime) || effectiveSampleTime <= 0.0)
    throw std::invalid_argument(
        "Feedforward effective sample time must be finite and positive.");

  const Real omega = 2.0 * PI * mCutoffFrequency;
  const Real decay = std::exp(-omega * effectiveSampleTime);
  const Real error = x(0, 0) - measuredVd;
  const Real updatedError =
      decay * ((1.0 + omega * effectiveSampleTime) * error +
               effectiveSampleTime * x(1, 0));
  const Real updatedDerivative =
      decay * ((1.0 - omega * effectiveSampleTime) * x(1, 0) -
               omega * omega * effectiveSampleTime * error);
  const Real filtered = measuredVd + updatedError;
  if (!std::isfinite(filtered) || filtered <= mMinimumDaxisVoltage)
    throw std::runtime_error("Feedforward filter produced an invalid Vd.");

  Matrix next(3, 1);
  next << filtered, updatedDerivative,
      (2.0 / 3.0) * activePowerReference / filtered;
  return next;
}
Real MMCSampledPowerFeedforward::sampleTime() const { return mSampleTime; }

void MMCSampledPowerFeedforward::evaluateStateDerivative(const Matrix &x,
                                                         const Matrix &u,
                                                         Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCSampledPowerFeedforward::evaluateOutput(const Matrix &x,
                                                const Matrix &u,
                                                Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCSampledPowerFeedforward::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCSampledPowerFeedforward::buildStateSpaceModel(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Matrix &E, Matrix &F, Real relativeStep,
    Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCSampledPowerFeedforward::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                               Real relativeStep,
                                               Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCSampledPowerFeedforward::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
