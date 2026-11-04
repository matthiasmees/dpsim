/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim/DirectLinearSolver.h>

#include <hip/hip_runtime_api.h>
#include <hipsolver/hipsolver.h>

namespace DPsim {
class HipDenseAdapter : public DirectLinearSolver {
  hipsolverHandle_t mSolver = nullptr;
  double *mMatrix = nullptr;
  double *mRightHandSide = nullptr;
  double *mWorkspace = nullptr;
  int *mPivotSequence = nullptr;
  int *mDeviceInfo = nullptr;
  UInt mSize = 0;
  int mWorkspaceSize = 0;

  void allocateDeviceMemory(UInt size);
  void copySystemMatrixToDevice(const SparseMatrix &systemMatrix);
  void factorizeOnDevice();

public:
  HipDenseAdapter();
  explicit HipDenseAdapter(CPS::Logger::Log log);
  ~HipDenseAdapter() override;

  /// Device-side RHS storage for HIP kernels in a device-resident solve path.
  /// The pointer remains owned by this adapter and is valid until the matrix
  /// size changes or the adapter is destroyed.
  double *deviceRightHandSide() noexcept { return mRightHandSide; }

  /// Solve the RHS already stored in deviceRightHandSide() in place. The
  /// operation is asynchronous with respect to the host and ordered on the
  /// HIP solver's stream.
  void solveDevice();

  /// Synchronize the device solve, check the solver info value, and copy its
  /// result into caller-owned host storage.
  void copyDeviceSolutionToHost(Matrix &solution);

  void preprocessing(SparseMatrix &systemMatrix,
                     std::vector<std::pair<UInt, UInt>>
                         &listVariableSystemMatrixEntries) override;
  void factorize(SparseMatrix &systemMatrix) override;
  void refactorize(SparseMatrix &systemMatrix) override;
  void partialRefactorize(SparseMatrix &systemMatrix,
                          std::vector<std::pair<UInt, UInt>>
                              &listVariableSystemMatrixEntries) override;
  Matrix solve(Matrix &rightSideVector) override;
};
} // namespace DPsim
