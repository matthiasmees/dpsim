// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/Signal/MMCStateSpaceBlock.h>

namespace CPS {
namespace Signal {

/// Hybrid sampled active-power feedforward controller.
/// Continuous state representation:
/// [vd_filtered, vd_filtered_derivative, id_ref_held].
class MMCSampledPowerFeedforward final : public MMCStateSpaceBlock {
public:
  MMCSampledPowerFeedforward();
  void setParameters(Real cutoffFrequency, Real sampleTime,
                     Real minimumDaxisVoltage);
  UInt stateSize() const override;
  UInt inputSize() const override;
  UInt outputSize() const override;
  std::vector<String> stateNames() const override;
  void evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                Matrix &y) const override;
  void evaluateStateDerivative(const Matrix &x, const Matrix &u,
                               Matrix &dx) const override;
  void evaluateOutput(const Matrix &x, const Matrix &u,
                      Matrix &y) const override;
  void calculateNumericalJacobians(const Matrix &x, const Matrix &u, Matrix &A,
                                   Matrix &B, Matrix &C, Matrix &D,
                                   Real relativeStep = 1e-6,
                                   Real absoluteStep = 1e-8) const override;
  void buildStateSpaceModel(const Matrix &x, const Matrix &u, Matrix &A,
                            Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                            Matrix &F, Real relativeStep = 1e-6,
                            Real absoluteStep = 1e-8) const override;
  MMCLinearization getStateSpaceModel(const Matrix &x, const Matrix &u,
                                      Real relativeStep = 1e-6,
                                      Real absoluteStep = 1e-8) const override;
  MMCSparseLinearization
  getSparseStateSpaceModel(const Matrix &x, const Matrix &u,
                           Real relativeStep = 1e-6, Real absoluteStep = 1e-8,
                           Real sparseTolerance = 1e-14) const override;

  Matrix initialize(Real measuredVd, Real activePowerReference) const;
  Matrix sample(const Matrix &x, Real measuredVd,
                Real activePowerReference) const;
  Matrix sample(const Matrix &x, Real measuredVd, Real activePowerReference,
                Real effectiveSampleTime) const;
  Real sampleTime() const;

private:
  Real mCutoffFrequency;
  Real mSampleTime;
  Real mMinimumDaxisVoltage;
};

} // namespace Signal
} // namespace CPS
