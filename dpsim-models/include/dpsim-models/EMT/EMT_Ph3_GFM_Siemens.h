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

/// Grid-forming VSI based on the cascaded control structure of the MATLAB
/// "New VSI Model" contained in converter_models.slx.
///
/// The complete controller is implemented directly in this class. No Signal
/// controller components or controller subclasses are required.
///
/// Control structure, entirely in per unit:
///
///   measured P/Q -> first-order measurement filters
///   filtered P -> P-f droop -> angle integrator
///   filtered Q -> Q-V droop
///              -> SRF capacitor-voltage PI controller
///              -> SRF filter-current PI controller
///              -> optional first-order PWM delay
///              -> dq-to-abc controlled voltage source
///
/// Both SRF controllers include the cross-decoupling/feedforward terms used by
/// the Simulink model. Unlike the source model, the physical L and C values are
/// converted to their dimensionally correct per-unit reactance/susceptance:
///
///   X_Lf_pu = omega_base * Lf / Z_base
///   B_Cf_pu = omega_base * Cf * Z_base
///
/// Electrical topology, stamped in SI units:
///
///   ideal 3-phase source -- Rf -- Lf -- PCC
///                                      |
///                                      Cf
///                                      |
///                                     GND
///
/// Base definitions:
///   S_base          total three-phase apparent power [VA]
///   V_base_LL       line-to-line RMS voltage [V]
///   V_base_ph_peak  phase-to-neutral peak voltage [V]
///   I_base_ph_peak  phase-current peak [A]
///   Z_base          V_base_LL^2 / S_base [Ohm]
///   omega_base      2*pi*f_base [rad/s]
class GFM_Siemens : public CompositePowerComp<Real>,
                    public Base::AvVoltageSourceInverterDQ,
                    public SharedFactory<GFM_Siemens> {
public:
  GFM_Siemens(String uid, String name,
              Logger::Level logLevel = Logger::Level::off);

  GFM_Siemens(String name, Logger::Level logLevel = Logger::Level::off)
      : GFM_Siemens(name, name, logLevel) {}

  // -----------------------------------------------------------------------
  // Base quantities and references
  // -----------------------------------------------------------------------

  void setBaseParameters(Real ratedApparentPower,
                         Real ratedVoltageLineToLineRms,
                         Real nominalFrequencyHz);

  /// Set canonical per-unit references.
  void setReferencesPerUnit(Real frequencyReferencePu, Real voltageReferencePu,
                            Real activePowerReferencePu,
                            Real reactivePowerReferencePu);

  /// SI convenience interface. Voltage is phase-to-neutral peak.
  void setReferences(Real frequencyReferenceHz, Real voltageReferencePhasePeak,
                     Real activePowerReference, Real reactivePowerReference);

  // -----------------------------------------------------------------------
  // Controller tuning
  // -----------------------------------------------------------------------

  /// Configure P-f and Q-V droop gains.
  ///
  /// omega_droop_pu = omega_ref_pu
  ///                + k_p_pu * (P_ref_pu - P_pu)
  ///
  /// vd_droop_pu = V_ref_pu
  ///             + k_q_pu * (Q_ref_pu - Q_pu)
  void setDroopParametersPerUnit(Real activePowerDroopPu,
                                 Real reactivePowerDroopPu);

  /// First-order low-pass filters applied to measured active and reactive
  /// power before the P-f and Q-V droop equations. Set either time constant
  /// to zero to bypass the corresponding measurement filter.
  void setPowerMeasurementFilterTimeConstants(Real activePowerTimeConstant,
                                              Real reactivePowerTimeConstant);

  /// SRF capacitor-voltage controller.
  /// Defaults extracted from the MATLAB New VSI Model:
  /// Kp=0.52, Ki=1.16 1/s, feedforward=1.
  void setVoltageControllerParameters(Real proportionalGain,
                                      Real integralGainPerSecond,
                                      Real outputCurrentFeedforwardGain = 1.0);

  /// SRF filter-current controller.
  /// Defaults extracted from the MATLAB New VSI Model:
  /// Kp=0.74, Ki=1.19 1/s, feedforward=1.
  void setCurrentControllerParameters(Real proportionalGain,
                                      Real integralGainPerSecond,
                                      Real pccVoltageFeedforwardGain = 1.0);

  /// First-order delay between current-controller voltage command and source.
  /// The MATLAB model uses 1 us independently for d and q.
  void setPwmDelayTimeConstant(Real timeConstantSeconds);

  /// Optional vector and frequency limits. Set a vector limit to 0 to disable
  /// it. Back-calculation is applied whenever a vector limit is active.
  void setControllerLimitsPerUnit(Real minimumFrequencyPu,
                                  Real maximumFrequencyPu,
                                  Real maximumCurrentReferencePu,
                                  Real maximumVoltageCommandPu);

  void withControl(Bool controlOn) { mWithControl = controlOn; }

  // -----------------------------------------------------------------------
  // Electrical filter
  // -----------------------------------------------------------------------

  /// Physical single-phase filter parameters.
  void setFilterParameters(Real filterInductance, Real filterCapacitance,
                           Real filterResistance);

  /// Per-unit filter parameters.
  /// xLfPu = omega_base*Lf/Z_base
  /// bCfPu = omega_base*Cf*Z_base
  /// rRfPu = Rf/Z_base
  void setFilterParametersPerUnit(Real xLfPu, Real bCfPu, Real rRfPu);

  Real baseApparentPower() const { return mBaseApparentPower; }
  Real baseVoltageLineToLineRms() const { return mBaseVoltageLineToLineRms; }
  Real baseVoltagePhasePeak() const { return mBaseVoltagePhasePeak; }
  Real baseCurrentPhasePeak() const { return mBaseCurrentPhasePeak; }
  Real baseImpedance() const { return mBaseImpedance; }
  Real baseFrequency() const { return mBaseFrequency; }
  Real baseOmega() const { return mBaseOmega; }

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
  // Public attributes for runtime stepping, diagnostics and logging
  // -----------------------------------------------------------------------

  // References and gains in PU.
  const Attribute<Real>::Ptr mActivePowerRefPu;
  const Attribute<Real>::Ptr mReactivePowerRefPu;
  const Attribute<Real>::Ptr mFrequencyRefPu;
  const Attribute<Real>::Ptr mVoltageRefPu;
  const Attribute<Real>::Ptr mActivePowerDroopPu;
  const Attribute<Real>::Ptr mReactivePowerDroopPu;

  const Attribute<Real>::Ptr mVoltageControllerKp;
  const Attribute<Real>::Ptr mVoltageControllerKi;
  const Attribute<Real>::Ptr mVoltageControllerFeedforward;
  const Attribute<Real>::Ptr mCurrentControllerKp;
  const Attribute<Real>::Ptr mCurrentControllerKi;
  const Attribute<Real>::Ptr mCurrentControllerFeedforward;

  // ABC measurements. Currents are generator-positive.
  const Attribute<Matrix>::Ptr mPccCurrent;
  const Attribute<Matrix>::Ptr mPccCurrentPu;
  const Attribute<Matrix>::Ptr mFilterCurrent;
  const Attribute<Matrix>::Ptr mFilterCurrentPu;
  const Attribute<Matrix>::Ptr mCapacitorCurrent;
  const Attribute<Matrix>::Ptr mCapacitorCurrentPu;

  // SRF measurements and commands in PU; each is [d,q]^T.
  const Attribute<Matrix>::Ptr mPccVoltageDqPu;
  const Attribute<Matrix>::Ptr mPccCurrentDqPu;
  const Attribute<Matrix>::Ptr mFilterCurrentDqPu;
  const Attribute<Matrix>::Ptr mVoltageReferenceDqPu;
  const Attribute<Matrix>::Ptr mCurrentReferenceDqPu;
  const Attribute<Matrix>::Ptr mVoltageControllerIntegralPu;
  const Attribute<Matrix>::Ptr mCurrentControllerIntegralPu;
  const Attribute<Matrix>::Ptr mVoltageCommandPrePwmDqPu;
  const Attribute<Matrix>::Ptr mVoltageCommandDqPu;

  // Outer-loop measurements and states.
  const Attribute<Real>::Ptr mElecActivePowerPu;
  const Attribute<Real>::Ptr mElecReactivePowerPu;
  const Attribute<Real>::Ptr mFilteredActivePowerPu;
  const Attribute<Real>::Ptr mFilteredReactivePowerPu;
  const Attribute<Real>::Ptr mVoltageMagnitudePu;
  const Attribute<Real>::Ptr mVoltageDroopReferencePu;
  const Attribute<Real>::Ptr mFrequencyDroopReferencePu;
  const Attribute<Real>::Ptr mFrequencyPu;
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
  const Attribute<Real>::Ptr mVoltageMagnitude;
  const Attribute<Real>::Ptr mFrequency;
  const Attribute<Real>::Ptr mOmega;

protected:
  void requireBaseParameters() const;
  void requireFilterParameters() const;

  void updateMeasurements();
  void updatePowerMeasurementFilters();
  void updateController(Int timeStepCount);
  void updateOpenLoop(Int timeStepCount);
  void updatePhysicalMirrors();
  void writeVoltageReference();

  Matrix abcToDq(const Matrix &abc, Real theta) const;
  Matrix dqToAbc(const Matrix &dq, Real theta) const;

  // Electrical subcomponents. These are physical MNA primitives, not
  // controller subclasses.
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

  // Physical and per-unit filter parameters.
  Real mFilterInductance = 0.0;
  Real mFilterCapacitance = 0.0;
  Real mFilterResistance = 0.0;
  Real mFilterInductiveReactancePu = 0.0;
  Real mFilterCapacitiveSusceptancePu = 0.0;
  Real mFilterResistancePu = 0.0;

  // Controller timing and limits.
  Real mTimeStep = 0.0;
  Real mActivePowerMeasurementTimeConstant = 0.1;
  Real mReactivePowerMeasurementTimeConstant = 0.1;
  Real mPwmDelayTimeConstant = 0.0;
  Real mMinimumFrequencyPu = 0.5;
  Real mMaximumFrequencyPu = 1.5;
  Real mMaximumCurrentReferencePu = 0.0;
  Real mMaximumVoltageCommandPu = 0.0;

  // Previous-step states needed for trapezoidal integration.
  Real mPreviousFrequencyPu = 1.0;
  Matrix mPreviousVoltageErrorPu = Matrix::Zero(2, 1);
  Matrix mPreviousCurrentErrorPu = Matrix::Zero(2, 1);
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
