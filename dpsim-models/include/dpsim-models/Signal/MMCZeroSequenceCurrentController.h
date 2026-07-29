// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/Signal/MMCStateSpaceBlock.h>

namespace CPS {
namespace Signal {

/// Zero-sequence circulating-current controller.
/// Input: [i_ref, i_measured, vdc]. State: [xi_z].
/// Output: [v_sigma_z, error].
class MMCZeroSequenceCurrentController final : public MMCStateSpaceBlock {
public:
  MMCZeroSequenceCurrentController();
  void setPI(Real kp, Real ki);
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
};

} // namespace Signal
} // namespace CPS
