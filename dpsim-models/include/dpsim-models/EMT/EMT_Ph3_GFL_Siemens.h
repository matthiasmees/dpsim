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
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim-models/Solver/MNAInterface.h>

#include <memory>

namespace CPS {
namespace EMT {
namespace Ph3 {

/// Grid-following converter reconstructed from the MATLAB/Simulink subsystem
/// "CSI Model" in converter_models.slx.
///
/// Per-unit control structure:
///
///   PCC voltage -> SRF PLL -> theta, omega_pll, |v_pcc|
///
///   P_cmd = P_ref - K_fP (f_pll - f_ref)
///   Q_cmd = Q_ref - K_VQ (|v_pcc| - V_ref)
///
///   [P_cmd,Q_cmd], v_pcc,dq -> complex power division -> i_ref,dq
///
///   filtered i_pcc,dq -> SRF current PI controller
///   filtered v_pcc,dq -> PCC-voltage feedforward
///   Lf cross-decoupling -> v_source,dq
///
///   v_source,dq -> inverse Park transform -> controlled three-phase source
///
/// The PCC-voltage feedforward is deliberately first-order filtered. Raw EMT
/// PCC voltage can contain fast topology-switching and filter transients that
/// should not be copied directly into an averaged converter voltage command.
///
/// Electrical topology:
///
///   ideal 3-phase source -- Rf -- Lf -- PCC
///                                      |
///                                      Cf
///                                      |
///                                     GND
class GFL_Siemens : public CompositePowerComp<Real>,
                    public Base::AvVoltageSourceInverterDQ,
                    public SharedFactory<GFL_Siemens> {
public:
  GFL_Siemens(String uid, String name,
              Logger::Level logLevel = Logger::Level::off);

  GFL_Siemens(String name, Logger::Level logLevel = Logger::Level::off)
      : GFL_Siemens(name, name, logLevel) {}

  void setBaseParameters(Real ratedApparentPower,
                         Real ratedVoltageLineToLineRms,
                         Real nominalFrequencyHz);

  void setReferencesPerUnit(Real frequencyReferencePu, Real voltageReferencePu,
                            Real activePowerReferencePu,
                            Real reactivePowerReferencePu);

  void setReferences(Real frequencyReferenceHz, Real voltageReferencePhasePeak,
                     Real activePowerReference, Real reactivePowerReference);

  void setDroopParametersPerUnit(Real frequencyToActivePowerGainPu,
                                 Real voltageToReactivePowerGainPu);

  void setPllParameters(Real proportionalGain, Real integralGainPerSecond);

  void setCurrentControllerParameters(Real proportionalGain,
                                      Real integralGainPerSecond,
                                      Real pccVoltageFeedforwardGain = 1.0);

  /// First-order measurement filters. A zero time constant bypasses the
  /// corresponding filter. Existing three-argument calls remain valid:
  /// tau_Vff defaults to 1 ms.
  void setMeasurementFilterTimeConstants(
      Real activePowerTimeConstant, Real reactivePowerTimeConstant,
      Real currentTimeConstant, Real voltageFeedforwardTimeConstant = 1.0e-3);

  void
  setControllerLimitsPerUnit(Real minimumFrequencyPu, Real maximumFrequencyPu,
                             Real maximumCurrentReferencePu,
                             Real maximumVoltageCommandPu,
                             Real minimumVoltageForCurrentReferencePu = 0.05);

  void withControl(Bool controlOn) { mWithControl = controlOn; }

  void setFilterParameters(Real filterInductance, Real filterCapacitance,
                           Real filterResistance);

  void setFilterParametersPerUnit(Real xLfPu, Real bCfPu, Real rRfPu);

  Real baseApparentPower() const { return mBaseApparentPower; }
  Real baseVoltageLineToLineRms() const { return mBaseVoltageLineToLineRms; }
  Real baseVoltagePhasePeak() const { return mBaseVoltagePhasePeak; }
  Real baseCurrentPhasePeak() const { return mBaseCurrentPhasePeak; }
  Real baseImpedance() const { return mBaseImpedance; }
  Real baseFrequency() const { return mBaseFrequency; }
  Real baseOmega() const { return mBaseOmega; }

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

  // References and outer gains.
  const Attribute<Real>::Ptr mActivePowerRefPu;
  const Attribute<Real>::Ptr mReactivePowerRefPu;
  const Attribute<Real>::Ptr mFrequencyRefPu;
  const Attribute<Real>::Ptr mVoltageRefPu;
  const Attribute<Real>::Ptr mFrequencyToActivePowerGainPu;
  const Attribute<Real>::Ptr mVoltageToReactivePowerGainPu;

  // PLL and current-controller gains.
  const Attribute<Real>::Ptr mPllKp;
  const Attribute<Real>::Ptr mPllKi;
  const Attribute<Real>::Ptr mCurrentControllerKp;
  const Attribute<Real>::Ptr mCurrentControllerKi;
  const Attribute<Real>::Ptr mCurrentControllerFeedforward;

  // ABC measurements. Currents use generator-positive orientation.
  const Attribute<Matrix>::Ptr mPccCurrent;
  const Attribute<Matrix>::Ptr mPccCurrentPu;
  const Attribute<Matrix>::Ptr mFilterCurrent;
  const Attribute<Matrix>::Ptr mFilterCurrentPu;
  const Attribute<Matrix>::Ptr mCapacitorCurrent;
  const Attribute<Matrix>::Ptr mCapacitorCurrentPu;

  // dq measurements and commands, each [d,q]^T in PU.
  const Attribute<Matrix>::Ptr mPccVoltageDqPu;
  const Attribute<Matrix>::Ptr mFilteredPccVoltageDqPu;
  const Attribute<Matrix>::Ptr mPccCurrentDqPu;
  const Attribute<Matrix>::Ptr mFilteredPccCurrentDqPu;
  const Attribute<Matrix>::Ptr mFilterCurrentDqPu;
  const Attribute<Matrix>::Ptr mCurrentReferenceDqPu;
  const Attribute<Matrix>::Ptr mCurrentErrorDqPu;
  const Attribute<Matrix>::Ptr mCurrentControllerIntegralPu;
  const Attribute<Matrix>::Ptr mVoltageCommandDqPu;

  // Power measurements and outer command values.
  const Attribute<Real>::Ptr mElecActivePowerPu;
  const Attribute<Real>::Ptr mElecReactivePowerPu;
  const Attribute<Real>::Ptr mFilteredActivePowerPu;
  const Attribute<Real>::Ptr mFilteredReactivePowerPu;
  const Attribute<Real>::Ptr mActivePowerCommandPu;
  const Attribute<Real>::Ptr mReactivePowerCommandPu;
  const Attribute<Real>::Ptr mVoltageMagnitudePu;

  // PLL states and outputs.
  const Attribute<Real>::Ptr mPllVoltageErrorPu;
  const Attribute<Real>::Ptr mPllIntegral;
  const Attribute<Real>::Ptr mFrequencyPu;
  const Attribute<Real>::Ptr mOmega;
  const Attribute<Real>::Ptr mTheta;

  // Source command and source voltage.
  const Attribute<Matrix>::Ptr mVsrefPu;
  const Attribute<Matrix>::Ptr mVsref;
  const Attribute<Matrix>::Ptr mVs;

  // SI mirrors.
  const Attribute<Real>::Ptr mElecActivePower;
  const Attribute<Real>::Ptr mElecReactivePower;
  const Attribute<Real>::Ptr mFilteredActivePower;
  const Attribute<Real>::Ptr mFilteredReactivePower;
  const Attribute<Real>::Ptr mFrequency;
  const Attribute<Real>::Ptr mVoltageMagnitude;

protected:
  void requireBaseParameters() const;
  void requireFilterParameters() const;

  void updateMeasurements();
  void updateMeasurementFilters();
  void updateController(Int timeStepCount);
  void updateOpenLoop(Int timeStepCount);
  void updatePhysicalMirrors();
  void writeVoltageReference();

  Matrix abcToDq(const Matrix &abc, Real theta) const;
  Matrix dqToAbc(const Matrix &dq, Real theta) const;

  std::shared_ptr<EMT::Ph3::VoltageSource> mSubControlledVoltageSource;
  std::shared_ptr<EMT::Ph3::Resistor> mSubFilterResistor;
  std::shared_ptr<EMT::Ph3::Inductor> mSubFilterInductor;
  std::shared_ptr<EMT::Ph3::Capacitor> mSubFilterCapacitor;

  Bool mWithControl = true;
  Bool mBaseParametersSet = false;
  Bool mFilterParametersSet = false;
  Bool mReferencesSet = false;
  Bool mControlStateInitialized = false;

  // Base quantities.
  Real mBaseApparentPower = 0.0;
  Real mBaseVoltageLineToLineRms = 0.0;
  Real mBaseVoltagePhasePeak = 0.0;
  Real mBaseCurrentPhasePeak = 0.0;
  Real mBaseImpedance = 0.0;
  Real mBaseFrequency = 0.0;
  Real mBaseOmega = 0.0;

  // Physical and PU filter parameters.
  Real mFilterInductance = 0.0;
  Real mFilterCapacitance = 0.0;
  Real mFilterResistance = 0.0;
  Real mFilterInductiveReactancePu = 0.0;
  Real mFilterCapacitiveSusceptancePu = 0.0;
  Real mFilterResistancePu = 0.0;

  // Timing, measurement filters and limits.
  Real mTimeStep = 0.0;
  Real mActivePowerMeasurementTimeConstant = 0.05;
  Real mReactivePowerMeasurementTimeConstant = 0.05;
  Real mCurrentMeasurementTimeConstant = 0.5e-3;
  Real mVoltageFeedforwardMeasurementTimeConstant = 1.0e-3;
  Real mMinimumFrequencyPu = 0.8;
  Real mMaximumFrequencyPu = 1.2;
  Real mMaximumCurrentReferencePu = 0.0;
  Real mMaximumVoltageCommandPu = 0.0;
  Real mMinimumVoltageForCurrentReferencePu = 0.05;

  // Previous-step quantities for trapezoidal integration.
  Real mPreviousPllVoltageErrorPu = 0.0;
  Real mPreviousFrequencyPu = 1.0;
  Matrix mPreviousCurrentErrorPu = Matrix::Zero(2, 1);
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
