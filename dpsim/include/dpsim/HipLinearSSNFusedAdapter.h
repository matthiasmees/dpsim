/* SPDX-License-Identifier: MPL-2.0 */
#pragma once

#include <cstddef>
#include <vector>

#include <dpsim-models/EMT/EMT_LinearSSNFused.h>

namespace DPsim {
/// Executes LinearSSNFused passive component state updates as batched HIP
/// kernels while preserving the component's existing host-side MNA stamps.
class HipLinearSSNFusedAdapter {
  struct DeviceState;
  CPS::EMT::LinearSSNFused *mComponent = nullptr;
  DeviceState *mDevice = nullptr;

public:
  explicit HipLinearSSNFusedAdapter(CPS::EMT::LinearSSNFused &component);
  ~HipLinearSSNFusedAdapter();
  HipLinearSSNFusedAdapter(const HipLinearSSNFusedAdapter &) = delete;
  HipLinearSSNFusedAdapter &
  operator=(const HipLinearSSNFusedAdapter &) = delete;

  void updateHistories();
  void updateStates(const CPS::Matrix &leftVector);
  void configureDeviceRhs(double *deviceRightHandSide);
  void accumulateDeviceRhs();
  void updateStatesDevice(const CPS::Matrix &leftVector,
                          const double *deviceSolution);
};
} // namespace DPsim
