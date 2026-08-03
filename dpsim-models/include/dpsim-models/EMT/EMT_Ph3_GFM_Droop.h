// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

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

namespace CPS {
namespace EMT {
namespace Ph3 {

/// Grid-forming voltage source with filtered P-f and Q-V droop control.
///
/// The electrical MNA subcomponents continue to operate in SI units. All
/// controller measurements, references, droop equations, limits, controller
/// states, and commands are evaluated internally in per unit.
///
/// Existing SI setter names are retained. They convert their arguments to the
/// configured per-unit base. Explicit ...PerUnit() setters are provided for
/// scale-independent controller design.
///
/// Base definitions:
///   S_base       : total three-phase apparent power [VA]
///   V_base_LL    : line-to-line RMS voltage [V]
///   V_base_peak  : phase-to-neutral peak voltage [V]
///   I_base_peak  : phase-current peak [A]
///   Z_base       : V_base_LL^2 / S_base [Ohm]
///   omega_base   : 2*pi*f_base [rad/s]
class GFM_Droop : public CompositePowerComp<Real>,
                  public Base::AvVoltageSourceInverterDQ,
                  public SharedFactory<GFM_Droop> {
public:
  GFM_Droop(String uid, String name,
            Logger::Level logLevel = Logger::Level::off,
            Bool withTrafo = false);

  GFM_Droop(String name, Logger::Level logLevel = Logger::Level::off,
            Bool withTrafo = false)
      : GFM_Droop(name, name, logLevel, withTrafo) {}

  // -----------------------------------------------------------------------
  // Base quantities
  // -----------------------------------------------------------------------

  /// Configure the controller and electrical conversion base.
  void setBaseParameters(Real ratedApparentPower,
                         Real ratedVoltageLineToLineRms,
                         Real nominalFrequencyHz);

  Real baseApparentPower() const { return mBaseApparentPower; }
  Real baseVoltageLineToLineRms() const { return mBaseVoltageLineToLineRms; }
  Real baseVoltagePhasePeak() const { return mBaseVoltagePhasePeak; }
  Real baseCurrentPhasePeak() const { return mBaseCurrentPhasePeak; }
  Real baseImpedance() const { return mBaseImpedance; }
  Real baseFrequency() const { return mBaseFrequency; }
  Real baseOmega() const { return mBaseOmega; }

  // -----------------------------------------------------------------------
  // References
  // -----------------------------------------------------------------------

  /// Backward-compatible SI interface.
  void setParameters(Real frequencyReferenceHz, Real voltageReferencePeak,
                     Real activePowerReference, Real reactivePowerReference);

  /// Explicit per-unit interface.
  void setParametersPerUnit(Real frequencyReferencePu, Real voltageReferencePu,
                            Real activePowerReferencePu,
                            Real reactivePowerReferencePu);

  // -----------------------------------------------------------------------
  // Controller tuning
  // -----------------------------------------------------------------------

  /// Backward-compatible SI interface: Hz/W, V/var, and 1/s.
  void setDroopParameters(Real activePowerDroopHzPerW,
                          Real reactivePowerDroopVPerVar,
                          Real voltageIntegralGainPerSecond);

  /// Per-unit droop gains.
  ///
  /// frequency_pu = frequency_ref_pu
  ///              + k_p_pu * (P_ref_pu - P_filtered_pu)
  ///
  /// V0_pu = V_ref_pu
  ///       + k_q_pu * (Q_ref_pu - Q_filtered_pu)
  ///
  /// A 5% droop is therefore entered directly as 0.05.
  void setDroopParametersPerUnit(Real activePowerDroopPu,
                                 Real reactivePowerDroopPu,
                                 Real voltageIntegralGainPerSecond);

  void setPowerFilterTimeConstant(Real timeConstantSeconds);

  /// Backward-compatible SI limits.
  void setControllerLimits(Real minimumFrequencyHz, Real maximumFrequencyHz,
                           Real minimumVoltagePeak, Real maximumVoltagePeak);

  /// Explicit per-unit limits.
  void setControllerLimitsPerUnit(Real minimumFrequencyPu,
                                  Real maximumFrequencyPu,
                                  Real minimumVoltagePu, Real maximumVoltagePu);

  // -----------------------------------------------------------------------
  // Electrical parameters
  // -----------------------------------------------------------------------

  void setTransformerParameters(Real nomVoltageEnd1, Real nomVoltageEnd2,
                                Real ratedPower, Real ratioAbs, Real ratioPhase,
                                Real resistance, Real inductance, Real omega);

  /// Backward-compatible physical filter values.
  void setFilterParameters(Real Lf, Real Cf, Real Rf, Real Rd);

  /// Scale-independent filter values.
  ///
  /// xLfPu = omega_base * Lf / Z_base
  /// bCfPu = omega_base * Cf * Z_base
  /// rRfPu = Rf / Z_base
  /// rRdPu = Rd / Z_base
  void setFilterParametersPerUnit(Real xLfPu, Real bCfPu, Real rRfPu,
                                  Real rRdPu);

  void withControl(Bool controlOn) { mWithControl = controlOn; }

  // -----------------------------------------------------------------------
  // MNA interface
  // -----------------------------------------------------------------------

  void initializeParentFromNodesAndTerminals(Real frequency) override;

  void mnaParentInitialize(Real omega, Real timeStep,
                           Attribute<Matrix>::Ptr leftVector) override;

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

  void mnaCompUpdateCurrent(const Matrix &leftVector) override;
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;

  // -----------------------------------------------------------------------
  // Existing SI attributes retained for compatibility and logging
  // -----------------------------------------------------------------------

  const Attribute<Real>::Ptr mActivePowerRef;
  const Attribute<Real>::Ptr mReactivePowerRef;
  const Attribute<Real>::Ptr mFrequencyRef;
  const Attribute<Real>::Ptr mVoltageRef;
  const Attribute<Real>::Ptr mActivePowerDroop;
  const Attribute<Real>::Ptr mReactivePowerDroop;
  const Attribute<Real>::Ptr mVoltageIntegralGain;

  const Attribute<Matrix>::Ptr mPccCurrent;
  const Attribute<Real>::Ptr mElecActivePower;
  const Attribute<Real>::Ptr mElecReactivePower;
  const Attribute<Real>::Ptr mFilteredActivePower;
  const Attribute<Real>::Ptr mFilteredReactivePower;
  const Attribute<Real>::Ptr mVoltageMagnitude;
  const Attribute<Real>::Ptr mFrequency;
  const Attribute<Real>::Ptr mOmega;
  const Attribute<Real>::Ptr mTheta;
  const Attribute<Real>::Ptr mVoltageDroopOutput;
  const Attribute<Real>::Ptr mVoltageIntegralState;
  const Attribute<Real>::Ptr mVoltageCommand;
  const Attribute<Matrix>::Ptr mVsref;
  const Attribute<Matrix>::Ptr mVs;

  // -----------------------------------------------------------------------
  // Canonical per-unit controller attributes
  // -----------------------------------------------------------------------

  const Attribute<Real>::Ptr mActivePowerRefPu;
  const Attribute<Real>::Ptr mReactivePowerRefPu;
  const Attribute<Real>::Ptr mFrequencyRefPu;
  const Attribute<Real>::Ptr mVoltageRefPu;
  const Attribute<Real>::Ptr mActivePowerDroopPu;
  const Attribute<Real>::Ptr mReactivePowerDroopPu;

  const Attribute<Matrix>::Ptr mPccCurrentPu;
  const Attribute<Real>::Ptr mElecActivePowerPu;
  const Attribute<Real>::Ptr mElecReactivePowerPu;
  const Attribute<Real>::Ptr mFilteredActivePowerPu;
  const Attribute<Real>::Ptr mFilteredReactivePowerPu;
  const Attribute<Real>::Ptr mVoltageMagnitudePu;
  const Attribute<Real>::Ptr mFrequencyPu;
  const Attribute<Real>::Ptr mOmegaPu;
  const Attribute<Real>::Ptr mVoltageDroopOutputPu;
  const Attribute<Real>::Ptr mVoltageIntegralStatePu;
  const Attribute<Real>::Ptr mVoltageCommandPu;
  const Attribute<Matrix>::Ptr mVsrefPu;

protected:
  void requireBaseParameters() const;
  void updateMeasurements();
  void updatePowerMeasurementFilter(Int timeStepCount);
  void updateController(Int timeStepCount);
  void updateOpenLoop(Int timeStepCount);
  void writeVoltageReference();
  void updatePhysicalMirrors();

  // Electrical subcomponents
  std::shared_ptr<EMT::Ph3::VoltageSource> mSubCtrledVoltageSource;
  std::shared_ptr<EMT::Ph3::Resistor> mSubResistorF;
  std::shared_ptr<EMT::Ph3::Resistor> mSubResistorD;
  std::shared_ptr<EMT::Ph3::Capacitor> mSubCapacitorF;
  std::shared_ptr<EMT::Ph3::Inductor> mSubInductorF;
  std::shared_ptr<EMT::Ph3::Transformer> mConnectionTransformer;

  Bool mWithConnectionTransformer = false;
  Bool mWithControl = true;
  Bool mControlStateInitialized = false;
  Bool mBaseParametersSet = false;

  // Base quantities
  Real mBaseApparentPower = 0.0;
  Real mBaseVoltageLineToLineRms = 0.0;
  Real mBaseVoltagePhasePeak = 0.0;
  Real mBaseCurrentPhasePeak = 0.0;
  Real mBaseImpedance = 0.0;
  Real mBaseFrequency = 0.0;
  Real mBaseOmega = 0.0;

  // Filter values in per unit for diagnostics
  Real mFilterInductiveReactancePu = 0.0;
  Real mFilterCapacitiveSusceptancePu = 0.0;
  Real mFilterResistancePu = 0.0;
  Real mCapacitorDampingResistancePu = 0.0;
  Real mCapacitorDampingResistance = 0.0;

  // Discrete controller state
  Real mTimeStep = 0.0;
  Real mPowerFilterTimeConstant = 0.0;
  Real mPowerFilterAlpha = 1.0;
  Real mMinimumFrequencyPu = 0.0;
  Real mMaximumFrequencyPu = 2.0;
  Real mMinimumVoltagePu = 0.0;
  Real mMaximumVoltagePu = 2.0;
  Real mPreviousOmegaPu = 1.0;
  Real mPreviousVoltageErrorPu = 0.0;
};

/// Compatibility alias for code that wants an explicitly named PU model.
/// It is the same extended implementation; no duplicate electrical model exists.
using GFM_Droop_PU = GFM_Droop;

} // namespace Ph3
} // namespace EMT
} // namespace CPS
