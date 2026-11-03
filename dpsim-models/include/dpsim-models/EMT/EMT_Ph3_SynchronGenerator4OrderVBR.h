/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <array>

#include <dpsim-models/EMT/EMT_Ph3_ReducedOrderSynchronGeneratorVBR.h>

namespace CPS {
namespace EMT {
namespace Ph3 {
/// @brief Voltage-Behind-Reactance (VBR) implementation
/// of 4th order synchronous generator model
class SynchronGenerator4OrderVBR
    : public ReducedOrderSynchronGeneratorVBR,
      public SharedFactory<SynchronGenerator4OrderVBR> {
public:
  /// Plain data needed by a specialized device executor for the benchmark's
  /// controller-free, Norton-form SG4 VBR model. Values are captured only
  /// after normal DPsim initialization has established per-unit state and MNA
  /// node indices.
  struct MetalDeviceData {
    std::array<Real, 3> voltageDq0;
    std::array<Real, 3> currentDq0;
    std::array<Real, 3> transientVoltageDq0;
    std::array<Real, 9> conductanceDq0;
    std::array<Real, 9> conductanceAbc;
    std::array<Real, 3> sourceAbc;
    std::array<UInt, 3> matrixNodeIndices;
    Real timeStep;
    Real baseVoltage;
    Real baseCurrent;
    Real baseImpedance;
    Real baseOmegaElec;
    Real baseOmegaMech;
    Real inertia;
    Real damping;
    Real ldTransient;
    Real lqTransient;
    Real adTransient;
    Real bdTransient;
    Real aqTransient;
    Real bqTransient;
    Real dqTransient;
    Real fieldVoltage;
    Real previousFieldVoltage;
    Real mechanicalTorque;
    Real previousMechanicalTorque;
    Real electricalTorque;
    Real rotorSpeed;
    Real rotorAngle;
    Real loadAngle;
  };

  // ### Model specific elements ###
  /// transient voltage
  const Attribute<Matrix>::Ptr mEdq0_t;

protected:
  /// history term of VBR
  Matrix mEhs_vbr;

public:
  ///
  SynchronGenerator4OrderVBR(const String &uid, const String &name,
                             Logger::Level logLevel = Logger::Level::off);
  ///
  SynchronGenerator4OrderVBR(const String &name,
                             Logger::Level logLevel = Logger::Level::off);

  // #### General Functions ####
  ///
  void specificInitialization() final;
  ///
  void stepInPerUnit() final;

  MetalDeviceData metalDeviceData();
};
} // namespace Ph3
} // namespace EMT
} // namespace CPS
