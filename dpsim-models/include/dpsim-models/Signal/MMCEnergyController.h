// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/Signal/MMCStateSpaceBlock.h>

namespace CPS {
namespace Signal {

/// Total-energy controller.
/// Input: [stored_energy, p_ac, vdc]. State: [xi_energy].
/// Output: [i_sigma_z_ref_raw, energy_error].
class MMCEnergyController final : public MMCStateSpaceBlock {
public:
  MMCEnergyController();
  void setParameters(Real kp, Real ki, Bool enabled, Real energyReference,
                     Real minimumDcVoltage);
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

private:
  Real mKp;
  Real mKi;
  Bool mEnabled;
  Real mEnergyReference;
  Real mMinimumDcVoltage;
};

} // namespace Signal
} // namespace CPS
