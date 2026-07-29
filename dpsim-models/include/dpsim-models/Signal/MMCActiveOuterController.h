// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/Signal/MMCStateSpaceBlock.h>

namespace CPS {
namespace Signal {

enum class MMCActiveMode {
  OpenLoop,
  ActivePower,
  DcVoltage,
  DcDroop,
  SampledPowerFeedforward
};

/// Active outer loop.
/// Input: [p, vdc, vd, held_id_ref].
/// State: [xi_active].
/// Output: [id_ref_raw, error, output_per_integrator_sign].
class MMCActiveOuterController final : public MMCStateSpaceBlock {
public:
  MMCActiveOuterController();
  void setPI(Real kp, Real ki);
  void setMode(MMCActiveMode mode);
  void setReferences(Real activePowerReference, Real dcVoltageReference,
                     Real droopGain, Real openLoopCurrentReference);
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
  MMCActiveMode mMode;
  Real mKp;
  Real mKi;
  Real mActivePowerReference;
  Real mDcVoltageReference;
  Real mDroopGain;
  Real mOpenLoopCurrentReference;
};

} // namespace Signal
} // namespace CPS
