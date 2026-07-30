// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
// PCC-filter revision 2026-07-30:
// The former main-path coupling resistor Rc has been removed.
// A passive damping resistor Rd is connected in series with the shunt
// capacitor Cf at the filter/PCC node.
// P/Q use the generator-positive PCC current i_pcc = -i_intf.
#pragma once

#include <dpsim-models/Base/Base_AvVoltageSourceInverterDQ.h>
#include <dpsim-models/CompositePowerComp.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/EMT/EMT_Ph3_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim-models/Solver/MNAInterface.h>

#include <memory>

namespace CPS {
namespace EMT {
namespace Ph3 {

/// Three-phase grid-forming inverter with P-f and Q-V droop control.
///
/// Generator-positive control equations:
///
///   tau_pq * dP_f/dt = P - P_f
///   tau_pq * dQ_f/dt = Q - Q_f
///   f                  = f_ref + k_p (P_ref - P_f)
///   V_0                = V_ref + k_q (Q_ref - Q_f)
///   d xi_v / dt        = k_iv (V_0 - V)
///   V_1                = V_0 + xi_v
///   d theta / dt       = 2 pi f
///
/// The commanded phase voltages are
///
///   v_a = V_1 cos(theta)
///   v_b = V_1 cos(theta - 2 pi / 3)
///   v_c = V_1 cos(theta + 2 pi / 3).
///
/// V_ref, V_0, V_1, and V are phase-to-neutral peak quantities.
/// P and Q are measured at the external GFM terminal (PCC), using the
/// terminal voltage and generator-positive current flowing into the grid.
///
/// The physical output filter is
///
///   source -- Rf -- Lf -- filter/PCC node
///                              |
///                              Rd
///                              |
///                              Cf
///                              |
///                             GND
///
/// There is no main-path coupling resistor Rc. Rd is a passive damping
/// resistor in series with the shunt capacitor branch. The additional
/// first-order P/Q filter is a measurement filter for the droop controller and
/// prevents instantaneous EMT power ripple from appearing directly in the
/// generated frequency and voltage amplitude.
class GFM_Droop : public Base::AvVoltageSourceInverterDQ,
                  public CompositePowerComp<Real>,
                  public SharedFactory<GFM_Droop> {
protected:
  // Simulation/controller state
  Real mTimeStep = 0.0;
  Real mPreviousOmega = 0.0;
  Real mPreviousVoltageError = 0.0;
  Bool mControlStateInitialized = false;

  // P/Q measurement filter. Default: 10 Hz cutoff.
  Real mPowerFilterTimeConstant = 0.015915494309189534;
  Real mPowerFilterAlpha = 1.0;

  // Optional controller limits. Defaults effectively disable limiting.
  Real mMinimumFrequency = 0.0;
  Real mMaximumFrequency = 1.0e6;
  Real mMinimumVoltage = 0.0;
  Real mMaximumVoltage = 1.0e12;

  // Passive damping resistor in series with the capacitor branch [Ohm/phase].
  Real mCapacitorDampingResistance = 0.0;

  // Electrical subcomponents
  std::shared_ptr<EMT::Ph3::VoltageSource> mSubCtrledVoltageSource;
  std::shared_ptr<EMT::Ph3::Resistor> mSubResistorF;
  std::shared_ptr<EMT::Ph3::Resistor> mSubResistorD;
  std::shared_ptr<EMT::Ph3::Capacitor> mSubCapacitorF;
  std::shared_ptr<EMT::Ph3::Inductor> mSubInductorF;
  std::shared_ptr<EMT::Ph3::Transformer> mConnectionTransformer;

  Bool mWithConnectionTransformer = false;
  Bool mWithControl = true;

  void updateMeasurements();
  void updatePowerMeasurementFilter(Int timeStepCount);
  void updateController(Int timeStepCount);
  void updateOpenLoop(Int timeStepCount);
  void writeVoltageReference();

public:
  // References. Public attributes can be stepped at runtime.
  /// Active-power reference [W], generator-positive.
  const Attribute<Real>::Ptr mActivePowerRef;
  /// Reactive-power reference [var], generator-positive.
  const Attribute<Real>::Ptr mReactivePowerRef;
  /// Nominal frequency reference [Hz].
  const Attribute<Real>::Ptr mFrequencyRef;
  /// PCC voltage-magnitude reference [V phase-to-neutral peak].
  const Attribute<Real>::Ptr mVoltageRef;

  // Controller gains
  /// P-f droop gain [Hz/W]. Set to zero for isochronous 50 Hz operation.
  const Attribute<Real>::Ptr mActivePowerDroop;
  /// Q-V droop gain [V/var].
  const Attribute<Real>::Ptr mReactivePowerDroop;
  /// Voltage-error integral gain [1/s].
  const Attribute<Real>::Ptr mVoltageIntegralGain;

  // Raw measurements at the external GFM terminal (PCC).
  // mIntfCurrent keeps DPsim's consumer-positive interface convention.
  // mPccCurrent is the opposite, generator-positive current flowing from
  // the GFM into the connected grid, and is the current used for P/Q.
  const Attribute<Matrix>::Ptr mPccCurrent;
  const Attribute<Real>::Ptr mElecActivePower;
  const Attribute<Real>::Ptr mElecReactivePower;

  // Filtered droop measurements
  const Attribute<Real>::Ptr mFilteredActivePower;
  const Attribute<Real>::Ptr mFilteredReactivePower;

  // Controller states/outputs
  const Attribute<Real>::Ptr mVoltageMagnitude;
  const Attribute<Real>::Ptr mFrequency;
  const Attribute<Real>::Ptr mOmega;
  const Attribute<Real>::Ptr mTheta;
  const Attribute<Real>::Ptr mVoltageDroopOutput;
  const Attribute<Real>::Ptr mVoltageIntegralState;
  const Attribute<Real>::Ptr mVoltageCommand;

  /// Instantaneous phase-voltage command [V phase-to-neutral peak].
  const Attribute<Matrix>::Ptr mVsref;
  /// Actual voltage across the controlled source.
  const Attribute<Matrix>::Ptr mVs;

  GFM_Droop(String name, Logger::Level logLevel = Logger::Level::off)
      : GFM_Droop(name, name, logLevel, false) {}

  GFM_Droop(String uid, String name,
            Logger::Level logLevel = Logger::Level::off,
            Bool withTrafo = false);

  /// Sets f_ref [Hz], V_ref [phase-neutral peak V], P_ref [W], and Q_ref [var].
  void setParameters(Real frequencyReferenceHz, Real voltageReferencePeak,
                     Real activePowerReference, Real reactivePowerReference);

  /// Sets k_p [Hz/W], k_q [V/var], and k_iv [1/s].
  void setDroopParameters(Real activePowerDroop, Real reactivePowerDroop,
                          Real voltageIntegralGain);

  /// Sets the common first-order P/Q measurement-filter time constant [s].
  /// A value of zero bypasses the filter.
  void setPowerFilterTimeConstant(Real timeConstant);

  /// Sets optional hard limits for anti-windup and numerical protection.
  void setControllerLimits(Real minimumFrequencyHz, Real maximumFrequencyHz,
                           Real minimumVoltagePeak, Real maximumVoltagePeak);

  void setTransformerParameters(Real nomVoltageEnd1, Real nomVoltageEnd2,
                                Real ratedPower, Real ratioAbs, Real ratioPhase,
                                Real resistance, Real inductance, Real omega);

  /// Sets the physical output filter:
  /// series Lf/Rf and a shunt Rd-Cf damping branch at the filter/PCC node.
  ///
  /// All resistance, inductance, and capacitance values are per phase.
  void setFilterParameters(Real Lf, Real Cf, Real Rf, Real Rd);

  /// Disabled control gives an open-loop sinusoidal source at f_ref and V_ref.
  void withControl(Bool controlOn) { mWithControl = controlOn; }

  void initializeParentFromNodesAndTerminals(Real frequency) override;

  void mnaParentInitialize(Real omega, Real timeStep,
                           Attribute<Matrix>::Ptr leftVector) override;

  void mnaCompUpdateCurrent(const Matrix &leftVector) override;
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;

  void mnaParentPreStep(Real time, Int timeStepCount) override;
  void mnaParentPostStep(Real time, Int timeStepCount,
                         Attribute<Matrix>::Ptr &leftVector) override;

  void mnaParentAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) override;

  void
  mnaParentAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                   AttributeBase::List &attributeDependencies,
                                   AttributeBase::List &modifiedAttributes,
                                   Attribute<Matrix>::Ptr &leftVector) override;
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
