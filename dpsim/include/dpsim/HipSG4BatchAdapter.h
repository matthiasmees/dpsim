/* SPDX-License-Identifier: MPL-2.0 */
#pragma once

#include <memory>
#include <vector>

#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator4OrderVBR.h>
#include <dpsim/Definitions.h>

namespace DPsim {
/// Batched device-resident state executor for unsupported-controller-free SG4
/// Norton generators in the IEEE118 EMT topology.
class HipSG4BatchAdapter {
  struct DeviceState;
  std::vector<std::shared_ptr<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>>
      mGenerators;
  DeviceState *mDevice = nullptr;
  bool mDeviceAssemblyConfigured = false;

public:
  explicit HipSG4BatchAdapter(
      std::vector<std::shared_ptr<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>>
          generators);
  ~HipSG4BatchAdapter();
  HipSG4BatchAdapter(const HipSG4BatchAdapter &) = delete;
  HipSG4BatchAdapter &operator=(const HipSG4BatchAdapter &) = delete;

  void preStep(CPS::Real time);
  void postStep(const CPS::Matrix &leftVector);
  void configureDeviceAssembly(const SparseMatrix &systemMatrix,
                               const SparseMatrix &baseSystemMatrix,
                               double *deviceMatrixValues,
                               double *deviceRightHandSide);
  void accumulateDeviceRhs();
  void scatterDeviceMatrix();
  void postStepDevice(const CPS::Matrix &leftVector,
                      const double *deviceSolution,
                      bool synchronizeStateToHost);
};
} // namespace DPsim
