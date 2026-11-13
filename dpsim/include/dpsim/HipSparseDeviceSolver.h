/* SPDX-License-Identifier: MPL-2.0 */
#pragma once

#include <dpsim/Definitions.h>

namespace DPsim {
/// Common device-buffer interface used by fully resident HIP MNA assembly.
class HipSparseDeviceSolver {
public:
  virtual ~HipSparseDeviceSolver() = default;
  virtual double *deviceRightHandSide() noexcept = 0;
  virtual double *deviceSolution() noexcept = 0;
  virtual double *deviceValues() noexcept = 0;
  virtual void clearDeviceRightHandSide() = 0;
  virtual void uploadRightHandSide(const Matrix &rightSideVector) = 0;
  virtual void copyDeviceSolutionToHost(Matrix &solution) = 0;
  virtual void refactorizeDevice() = 0;
  virtual void solveDevice(bool verifyConvergence = true) = 0;
};
} // namespace DPsim
