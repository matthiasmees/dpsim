// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
//
// Per-unit controller revision 2026-07-31:
// - Electrical MNA subcomponents remain in SI units.
// - Measurements are transformed to per unit at the PCC.
// - P-f droop, Q-V droop, filtering, limits, integral state, frequency,
//   voltage command, and source reference are calculated in per unit.
// - Existing SI setters and attributes are retained as conversion interfaces.
// - Explicit ...PerUnit() setters and PU logging attributes are added.
// - The first scheduled source-angle step is advanced at timeStepCount == 0.

#include <dpsim-models/EMT/EMT_Ph3_GFM_Droop.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;

namespace {
constexpr Real TWO_PI = 2.0 * PI;
constexpr Real TWO_PI_OVER_THREE = 2.0 * PI / 3.0;
constexpr Real SQRT_TWO_OVER_THREE = 0.81649658092772603273;
constexpr Real SQRT_THREE_OVER_TWO = 0.86602540378443864676;
constexpr Real TWO_OVER_THREE = 2.0 / 3.0;

Real clampFinite(Real value, Real lower, Real upper, Real fallback) {
  if (!std::isfinite(value))
    return fallback;
  return std::clamp(value, lower, upper);
}
} // namespace

EMT::Ph3::GFM_Droop::GFM_Droop(String uid, String name, Logger::Level logLevel,
                               Bool withTrafo)
    : CompositePowerComp<Real>(uid, name, true, true, logLevel),
      // Existing SI interface attributes.
      mActivePowerRef(mAttributes->create<Real>("P_ref", 0.0)),
      mReactivePowerRef(mAttributes->create<Real>("Q_ref", 0.0)),
      mFrequencyRef(mAttributes->create<Real>("f_ref", 50.0)),
      mVoltageRef(mAttributes->create<Real>("V_ref", 0.0)),
      mActivePowerDroop(mAttributes->create<Real>("k_p", 0.0)),
      mReactivePowerDroop(mAttributes->create<Real>("k_q", 0.0)),
      mVoltageIntegralGain(mAttributes->create<Real>("k_iv", 0.0)),
      mPccCurrent(mAttributes->create<Matrix>("i_pcc", Matrix::Zero(3, 1))),
      mElecActivePower(mAttributes->create<Real>("P_elec", 0.0)),
      mElecReactivePower(mAttributes->create<Real>("Q_elec", 0.0)),
      mFilteredActivePower(mAttributes->create<Real>("P_filtered", 0.0)),
      mFilteredReactivePower(mAttributes->create<Real>("Q_filtered", 0.0)),
      mVoltageMagnitude(mAttributes->create<Real>("V_magnitude", 0.0)),
      mFrequency(mAttributes->create<Real>("frequency", 50.0)),
      mOmega(mAttributes->create<Real>("omega", TWO_PI * 50.0)),
      mTheta(mAttributes->create<Real>("theta", 0.0)),
      mVoltageDroopOutput(mAttributes->create<Real>("V0", 0.0)),
      mVoltageIntegralState(
          mAttributes->create<Real>("voltage_integrator", 0.0)),
      mVoltageCommand(mAttributes->create<Real>("V1", 0.0)),
      mVsref(mAttributes->create<Matrix>("Vsref", Matrix::Zero(3, 1))),
      mVs(mAttributes->createDynamic<Matrix>("Vs")),
      // Canonical PU controller attributes.
      mActivePowerRefPu(mAttributes->create<Real>("P_ref_pu", 0.0)),
      mReactivePowerRefPu(mAttributes->create<Real>("Q_ref_pu", 0.0)),
      mFrequencyRefPu(mAttributes->create<Real>("f_ref_pu", 1.0)),
      mVoltageRefPu(mAttributes->create<Real>("V_ref_pu", 1.0)),
      mActivePowerDroopPu(mAttributes->create<Real>("k_p_pu", 0.0)),
      mReactivePowerDroopPu(mAttributes->create<Real>("k_q_pu", 0.0)),
      mPccCurrentPu(
          mAttributes->create<Matrix>("i_pcc_pu", Matrix::Zero(3, 1))),
      mElecActivePowerPu(mAttributes->create<Real>("P_elec_pu", 0.0)),
      mElecReactivePowerPu(mAttributes->create<Real>("Q_elec_pu", 0.0)),
      mFilteredActivePowerPu(mAttributes->create<Real>("P_filtered_pu", 0.0)),
      mFilteredReactivePowerPu(mAttributes->create<Real>("Q_filtered_pu", 0.0)),
      mVoltageMagnitudePu(mAttributes->create<Real>("V_magnitude_pu", 0.0)),
      mFrequencyPu(mAttributes->create<Real>("frequency_pu", 1.0)),
      mOmegaPu(mAttributes->create<Real>("omega_pu", 1.0)),
      mVoltageDroopOutputPu(mAttributes->create<Real>("V0_pu", 0.0)),
      mVoltageIntegralStatePu(
          mAttributes->create<Real>("voltage_integrator_pu", 0.0)),
      mVoltageCommandPu(mAttributes->create<Real>("V1_pu", 0.0)),
      mVsrefPu(mAttributes->create<Matrix>("Vsref_pu", Matrix::Zero(3, 1))) {

  mPhaseType = PhaseType::ABC;
  setTerminalNumber(1);

  if (withTrafo) {
    // v0: source terminal
    // v1: between Rf and Lf
    // v2: filter/transformer node
    // v3: between Rd and Cf
    setVirtualNodeNumber(4);
    mConnectionTransformer = EMT::Ph3::Transformer::make(
        **mName + "_trans", **mName + "_trans", mLogLevel, false);
    addMNASubComponent(mConnectionTransformer,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                       MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  } else {
    // v0: source terminal
    // v1: between Rf and Lf
    // v2: between Rd and Cf; the Lf/Rd node is the external PCC
    setVirtualNodeNumber(3);
  }
  mWithConnectionTransformer = withTrafo;

  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(3, 1);

  mSubResistorF = EMT::Ph3::Resistor::make(**mName + "_resF", mLogLevel);
  mSubResistorD = EMT::Ph3::Resistor::make(**mName + "_resD", mLogLevel);
  mSubCapacitorF = EMT::Ph3::Capacitor::make(**mName + "_capF", mLogLevel);
  mSubInductorF = EMT::Ph3::Inductor::make(**mName + "_indF", mLogLevel);
  mSubCtrledVoltageSource =
      EMT::Ph3::VoltageSource::make(**mName + "_src", mLogLevel);

  addMNASubComponent(mSubResistorF, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  addMNASubComponent(mSubResistorD, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  addMNASubComponent(mSubCapacitorF, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  addMNASubComponent(mSubInductorF, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  // Source pre-step is executed explicitly after the PU controller update.
  addMNASubComponent(mSubCtrledVoltageSource, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  mVs->setReference(mSubCtrledVoltageSource->mIntfVoltage);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", type(), name);
  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Droop uses an entirely per-unit P-f/Q-V controller. "
      "The Rf-Lf and series Rd-Cf electrical filter remains stamped in SI.");
}

void EMT::Ph3::GFM_Droop::requireBaseParameters() const {
  if (!mBaseParametersSet) {
    throw std::runtime_error(
        "GFM_Droop base parameters must be configured before references, "
        "controller gains, limits, filter parameters, or initialization");
  }
}

void EMT::Ph3::GFM_Droop::setBaseParameters(Real ratedApparentPower,
                                            Real ratedVoltageLineToLineRms,
                                            Real nominalFrequencyHz) {
  if (!(ratedApparentPower > 0.0) || !(ratedVoltageLineToLineRms > 0.0) ||
      !(nominalFrequencyHz > 0.0) || !std::isfinite(ratedApparentPower) ||
      !std::isfinite(ratedVoltageLineToLineRms) ||
      !std::isfinite(nominalFrequencyHz)) {
    throw std::invalid_argument(
        "GFM_Droop requires positive finite S_base, V_base_LL_RMS, and "
        "f_base");
  }
  if (mControlStateInitialized) {
    throw std::logic_error(
        "GFM_Droop base parameters cannot be changed after initialization");
  }

  mBaseApparentPower = ratedApparentPower;
  mBaseVoltageLineToLineRms = ratedVoltageLineToLineRms;
  mBaseVoltagePhasePeak = RMS3PH_TO_PEAK1PH * ratedVoltageLineToLineRms;
  mBaseCurrentPhasePeak =
      2.0 * ratedApparentPower / (3.0 * mBaseVoltagePhasePeak);
  mBaseImpedance = ratedVoltageLineToLineRms * ratedVoltageLineToLineRms /
                   ratedApparentPower;
  mBaseFrequency = nominalFrequencyHz;
  mBaseOmega = TWO_PI * nominalFrequencyHz;
  mBaseParametersSet = true;

  // Non-restrictive defaults. Production examples should set explicit limits.
  mMinimumFrequencyPu = 0.0;
  mMaximumFrequencyPu = 2.0;
  mMinimumVoltagePu = 0.0;
  mMaximumVoltagePu = 2.0;

  SPDLOG_LOGGER_INFO(mSLog,
                     "GFM_Droop PU base:"
                     "\n  S_base={} VA"
                     "\n  V_base_LL_RMS={} V"
                     "\n  V_base_phase_peak={} V"
                     "\n  I_base_phase_peak={} A"
                     "\n  Z_base={} Ohm"
                     "\n  f_base={} Hz"
                     "\n  omega_base={} rad/s",
                     mBaseApparentPower, mBaseVoltageLineToLineRms,
                     mBaseVoltagePhasePeak, mBaseCurrentPhasePeak,
                     mBaseImpedance, mBaseFrequency, mBaseOmega);
}

void EMT::Ph3::GFM_Droop::setParameters(Real frequencyReferenceHz,
                                        Real voltageReferencePeak,
                                        Real activePowerReference,
                                        Real reactivePowerReference) {
  requireBaseParameters();

  setParametersPerUnit(frequencyReferenceHz / mBaseFrequency,
                       voltageReferencePeak / mBaseVoltagePhasePeak,
                       activePowerReference / mBaseApparentPower,
                       reactivePowerReference / mBaseApparentPower);
}

void EMT::Ph3::GFM_Droop::setParametersPerUnit(Real frequencyReferencePu,
                                               Real voltageReferencePu,
                                               Real activePowerReferencePu,
                                               Real reactivePowerReferencePu) {
  requireBaseParameters();

  if (!(frequencyReferencePu > 0.0) || !(voltageReferencePu >= 0.0) ||
      !std::isfinite(frequencyReferencePu) ||
      !std::isfinite(voltageReferencePu) ||
      !std::isfinite(activePowerReferencePu) ||
      !std::isfinite(reactivePowerReferencePu)) {
    throw std::invalid_argument(
        "GFM_Droop PU references must be finite; f_ref_pu must be positive "
        "and V_ref_pu must be non-negative");
  }

  mParametersSet = true;

  **mFrequencyRefPu = frequencyReferencePu;
  **mVoltageRefPu = voltageReferencePu;
  **mActivePowerRefPu = activePowerReferencePu;
  **mReactivePowerRefPu = reactivePowerReferencePu;

  // SI mirrors retained for compatibility.
  **mFrequencyRef = frequencyReferencePu * mBaseFrequency;
  **mVoltageRef = voltageReferencePu * mBaseVoltagePhasePeak;
  **mActivePowerRef = activePowerReferencePu * mBaseApparentPower;
  **mReactivePowerRef = reactivePowerReferencePu * mBaseApparentPower;

  SPDLOG_LOGGER_INFO(mSLog,
                     "GFM_Droop PU references:"
                     "\n  f_ref={} pu ({} Hz)"
                     "\n  V_ref={} pu ({} V phase peak)"
                     "\n  P_ref={} pu ({} W)"
                     "\n  Q_ref={} pu ({} var)",
                     **mFrequencyRefPu, **mFrequencyRef, **mVoltageRefPu,
                     **mVoltageRef, **mActivePowerRefPu, **mActivePowerRef,
                     **mReactivePowerRefPu, **mReactivePowerRef);
}

void EMT::Ph3::GFM_Droop::setDroopParameters(
    Real activePowerDroopHzPerW, Real reactivePowerDroopVPerVar,
    Real voltageIntegralGainPerSecond) {
  requireBaseParameters();

  setDroopParametersPerUnit(
      activePowerDroopHzPerW * mBaseApparentPower / mBaseFrequency,
      reactivePowerDroopVPerVar * mBaseApparentPower / mBaseVoltagePhasePeak,
      voltageIntegralGainPerSecond);
}

void EMT::Ph3::GFM_Droop::setDroopParametersPerUnit(
    Real activePowerDroopPu, Real reactivePowerDroopPu,
    Real voltageIntegralGainPerSecond) {
  requireBaseParameters();

  if (!std::isfinite(activePowerDroopPu) ||
      !std::isfinite(reactivePowerDroopPu) ||
      !(voltageIntegralGainPerSecond >= 0.0) ||
      !std::isfinite(voltageIntegralGainPerSecond)) {
    throw std::invalid_argument(
        "GFM_Droop PU gains must be finite and k_iv must be non-negative");
  }

  **mActivePowerDroopPu = activePowerDroopPu;
  **mReactivePowerDroopPu = reactivePowerDroopPu;
  **mVoltageIntegralGain = voltageIntegralGainPerSecond;

  // SI mirrors retained for compatibility and diagnostics.
  **mActivePowerDroop =
      activePowerDroopPu * mBaseFrequency / mBaseApparentPower;
  **mReactivePowerDroop =
      reactivePowerDroopPu * mBaseVoltagePhasePeak / mBaseApparentPower;

  SPDLOG_LOGGER_INFO(mSLog,
                     "GFM_Droop PU gains:"
                     "\n  k_p={} pu_f/pu_P ({} Hz/W)"
                     "\n  k_q={} pu_V/pu_Q ({} V/var)"
                     "\n  k_iv={} 1/s",
                     **mActivePowerDroopPu, **mActivePowerDroop,
                     **mReactivePowerDroopPu, **mReactivePowerDroop,
                     **mVoltageIntegralGain);
}

void EMT::Ph3::GFM_Droop::setPowerFilterTimeConstant(Real timeConstantSeconds) {
  if (!(timeConstantSeconds >= 0.0) || !std::isfinite(timeConstantSeconds)) {
    throw std::invalid_argument(
        "GFM_Droop P/Q filter time constant must be finite and non-negative");
  }

  mPowerFilterTimeConstant = timeConstantSeconds;

  if (timeConstantSeconds > 0.0) {
    SPDLOG_LOGGER_INFO(
        mSLog, "GFM_Droop P/Q measurement filter: tau={} s, cutoff={} Hz",
        timeConstantSeconds, 1.0 / (TWO_PI * timeConstantSeconds));
  } else {
    SPDLOG_LOGGER_INFO(mSLog, "GFM_Droop P/Q measurement filter bypassed");
  }
}

void EMT::Ph3::GFM_Droop::setControllerLimits(Real minimumFrequencyHz,
                                              Real maximumFrequencyHz,
                                              Real minimumVoltagePeak,
                                              Real maximumVoltagePeak) {
  requireBaseParameters();

  setControllerLimitsPerUnit(minimumFrequencyHz / mBaseFrequency,
                             maximumFrequencyHz / mBaseFrequency,
                             minimumVoltagePeak / mBaseVoltagePhasePeak,
                             maximumVoltagePeak / mBaseVoltagePhasePeak);
}

void EMT::Ph3::GFM_Droop::setControllerLimitsPerUnit(Real minimumFrequencyPu,
                                                     Real maximumFrequencyPu,
                                                     Real minimumVoltagePu,
                                                     Real maximumVoltagePu) {
  requireBaseParameters();

  if (!(minimumFrequencyPu >= 0.0) ||
      !(maximumFrequencyPu > minimumFrequencyPu) ||
      !(minimumVoltagePu >= 0.0) || !(maximumVoltagePu > minimumVoltagePu) ||
      !std::isfinite(minimumFrequencyPu) ||
      !std::isfinite(maximumFrequencyPu) || !std::isfinite(minimumVoltagePu) ||
      !std::isfinite(maximumVoltagePu)) {
    throw std::invalid_argument("Invalid GFM_Droop PU controller limits");
  }

  mMinimumFrequencyPu = minimumFrequencyPu;
  mMaximumFrequencyPu = maximumFrequencyPu;
  mMinimumVoltagePu = minimumVoltagePu;
  mMaximumVoltagePu = maximumVoltagePu;

  SPDLOG_LOGGER_INFO(mSLog, "GFM_Droop PU limits: f=[{}, {}] pu, V=[{}, {}] pu",
                     mMinimumFrequencyPu, mMaximumFrequencyPu,
                     mMinimumVoltagePu, mMaximumVoltagePu);
}

void EMT::Ph3::GFM_Droop::setTransformerParameters(
    Real nomVoltageEnd1, Real nomVoltageEnd2, Real ratedPower, Real ratioAbs,
    Real ratioPhase, Real resistance, Real inductance, Real omega) {
  Base::AvVoltageSourceInverterDQ::setTransformerParameters(
      nomVoltageEnd1, nomVoltageEnd2, ratedPower, ratioAbs, ratioPhase,
      resistance, inductance);

  if (mWithConnectionTransformer) {
    mConnectionTransformer->setParameters(
        mTransformerNominalVoltageEnd1, mTransformerNominalVoltageEnd2,
        mTransformerRatedPower, mTransformerRatioAbs, mTransformerRatioPhase,
        CPS::Math::singlePhaseParameterToThreePhase(mTransformerResistance),
        CPS::Math::singlePhaseParameterToThreePhase(mTransformerInductance));
  }

  (void)omega;
}

void EMT::Ph3::GFM_Droop::setFilterParameters(Real Lf, Real Cf, Real Rf,
                                              Real Rd) {
  requireBaseParameters();

  if (!(Lf > 0.0) || !(Cf > 0.0) || !(Rf >= 0.0) || !(Rd > 0.0) ||
      !std::isfinite(Lf) || !std::isfinite(Cf) || !std::isfinite(Rf) ||
      !std::isfinite(Rd)) {
    throw std::invalid_argument(
        "GFM_Droop filter requires Lf>0, Cf>0, Rf>=0, and Rd>0");
  }

  Base::AvVoltageSourceInverterDQ::setFilterParameters(Lf, Cf, Rf, 0.0);
  mCapacitorDampingResistance = Rd;

  mFilterInductiveReactancePu = mBaseOmega * mLf / mBaseImpedance;
  mFilterCapacitiveSusceptancePu = mBaseOmega * mCf * mBaseImpedance;
  mFilterResistancePu = mRf / mBaseImpedance;
  mCapacitorDampingResistancePu = mCapacitorDampingResistance / mBaseImpedance;

  mSubResistorF->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mRf));
  mSubResistorD->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mCapacitorDampingResistance));
  mSubInductorF->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mLf));
  mSubCapacitorF->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mCf));

  SPDLOG_LOGGER_INFO(mSLog,
                     "GFM_Droop electrical filter:"
                     "\n  X_Lf={} pu, B_Cf={} pu, Rf={} pu, Rd={} pu"
                     "\n  Lf={} H, Cf={} F, Rf={} Ohm, Rd={} Ohm",
                     mFilterInductiveReactancePu,
                     mFilterCapacitiveSusceptancePu, mFilterResistancePu,
                     mCapacitorDampingResistancePu, mLf, mCf, mRf,
                     mCapacitorDampingResistance);
}

void EMT::Ph3::GFM_Droop::setFilterParametersPerUnit(Real xLfPu, Real bCfPu,
                                                     Real rRfPu, Real rRdPu) {
  requireBaseParameters();

  if (!(xLfPu > 0.0) || !(bCfPu > 0.0) || !(rRfPu >= 0.0) || !(rRdPu > 0.0) ||
      !std::isfinite(xLfPu) || !std::isfinite(bCfPu) || !std::isfinite(rRfPu) ||
      !std::isfinite(rRdPu)) {
    throw std::invalid_argument(
        "GFM_Droop PU filter requires X_Lf>0, B_Cf>0, Rf>=0, and Rd>0");
  }

  const Real Lf = xLfPu * mBaseImpedance / mBaseOmega;
  const Real Cf = bCfPu / (mBaseOmega * mBaseImpedance);
  const Real Rf = rRfPu * mBaseImpedance;
  const Real Rd = rRdPu * mBaseImpedance;

  setFilterParameters(Lf, Cf, Rf, Rd);
}

void EMT::Ph3::GFM_Droop::initializeParentFromNodesAndTerminals(
    Real frequency) {
  requireBaseParameters();

  if (!mParametersSet) {
    throw std::runtime_error(
        "GFM_Droop references must be configured before initialization");
  }

  MatrixComp intfVoltageComplex = MatrixComp::Zero(3, 1);
  MatrixComp intfCurrentComplex = MatrixComp::Zero(3, 1);

  const Real activePower = terminal(0)->singlePower().real();
  const Real reactivePower = terminal(0)->singlePower().imag();

  intfVoltageComplex(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  intfVoltageComplex(1, 0) = intfVoltageComplex(0, 0) * SHIFT_TO_PHASE_B;
  intfVoltageComplex(2, 0) = intfVoltageComplex(0, 0) * SHIFT_TO_PHASE_C;

  intfCurrentComplex(0, 0) =
      -std::conj(TWO_OVER_THREE * Complex(activePower, reactivePower) /
                 intfVoltageComplex(0, 0));
  intfCurrentComplex(1, 0) = intfCurrentComplex(0, 0) * SHIFT_TO_PHASE_B;
  intfCurrentComplex(2, 0) = intfCurrentComplex(0, 0) * SHIFT_TO_PHASE_C;

  MatrixComp filterInterfaceInitialVoltage = MatrixComp::Zero(3, 1);
  MatrixComp filterInterfaceInitialCurrent = MatrixComp::Zero(3, 1);

  if (mWithConnectionTransformer) {
    filterInterfaceInitialVoltage =
        (intfVoltageComplex -
         Complex(mTransformerResistance,
                 mTransformerInductance * TWO_PI * frequency) *
             intfCurrentComplex) /
        Complex(mTransformerRatioAbs, mTransformerRatioPhase);
    filterInterfaceInitialCurrent =
        intfCurrentComplex *
        Complex(mTransformerRatioAbs, mTransformerRatioPhase);

    mVirtualNodes[2]->setInitialVoltage(PEAK1PH_TO_RMS3PH *
                                        filterInterfaceInitialVoltage);
    mConnectionTransformer->connect({mTerminals[0]->node(), mVirtualNodes[2]});
  } else {
    filterInterfaceInitialVoltage = intfVoltageComplex;
    filterInterfaceInitialCurrent = intfCurrentComplex;
  }

  const Real omegaInit = TWO_PI * frequency;

  // filter/PCC -- Rd -- capacitor node -- Cf -- GND
  const Complex capacitorImpedance(0.0, -1.0 / (omegaInit * mCf));
  const Complex branchImpedance(mCapacitorDampingResistance,
                                capacitorImpedance.imag());

  MatrixComp dampingBranchCurrentInit = MatrixComp::Zero(3, 1);
  MatrixComp capacitorVoltageInit = MatrixComp::Zero(3, 1);

  for (UInt phase = 0; phase < 3; ++phase) {
    dampingBranchCurrentInit(phase, 0) =
        filterInterfaceInitialVoltage(phase, 0) / branchImpedance;
    capacitorVoltageInit(phase, 0) =
        dampingBranchCurrentInit(phase, 0) * capacitorImpedance;
  }

  // Current orientation retained from the working SI class.
  const MatrixComp inductorCurrentInit =
      filterInterfaceInitialCurrent - dampingBranchCurrentInit;

  const MatrixComp vfInit = filterInterfaceInitialVoltage -
                            inductorCurrentInit * Complex(0.0, omegaInit * mLf);

  const MatrixComp vsInit = vfInit - inductorCurrentInit * mRf;

  mVirtualNodes[0]->setInitialVoltage(PEAK1PH_TO_RMS3PH * vsInit);
  mVirtualNodes[1]->setInitialVoltage(PEAK1PH_TO_RMS3PH * vfInit);

  if (mWithConnectionTransformer) {
    mVirtualNodes[2]->setInitialVoltage(PEAK1PH_TO_RMS3PH *
                                        filterInterfaceInitialVoltage);
    mVirtualNodes[3]->setInitialVoltage(PEAK1PH_TO_RMS3PH *
                                        capacitorVoltageInit);
  } else {
    mVirtualNodes[2]->setInitialVoltage(PEAK1PH_TO_RMS3PH *
                                        capacitorVoltageInit);
  }

  **mIntfVoltage = intfVoltageComplex.real();
  **mIntfCurrent = intfCurrentComplex.real();
  **mPccCurrent = -intfCurrentComplex.real();
  **mPccCurrentPu = **mPccCurrent / mBaseCurrentPhasePeak;

  const Complex initialPower =
      1.5 * intfVoltageComplex(0, 0) * std::conj(-intfCurrentComplex(0, 0));

  **mElecActivePowerPu = initialPower.real() / mBaseApparentPower;
  **mElecReactivePowerPu = initialPower.imag() / mBaseApparentPower;
  **mFilteredActivePowerPu = **mElecActivePowerPu;
  **mFilteredReactivePowerPu = **mElecReactivePowerPu;

  **mVoltageMagnitudePu =
      std::abs(intfVoltageComplex(0, 0)) / mBaseVoltagePhasePeak;

  mSubCtrledVoltageSource->setParameters(mVirtualNodes[0]->initialVoltage(),
                                         0.0);

  mSubCtrledVoltageSource->connect({SimNode::GND, mVirtualNodes[0]});
  mSubResistorF->connect({mVirtualNodes[0], mVirtualNodes[1]});

  if (mWithConnectionTransformer) {
    mSubInductorF->connect({mVirtualNodes[1], mVirtualNodes[2]});
    mSubResistorD->connect({mVirtualNodes[2], mVirtualNodes[3]});
    mSubCapacitorF->connect({mVirtualNodes[3], SimNode::GND});
  } else {
    mSubInductorF->connect({mVirtualNodes[1], mTerminals[0]->node()});
    mSubResistorD->connect({mTerminals[0]->node(), mVirtualNodes[2]});
    mSubCapacitorF->connect({mVirtualNodes[2], SimNode::GND});
  }

  for (const auto &subcomp : mSubComponents) {
    subcomp->initialize(mFrequencies);
    subcomp->initializeFromNodesAndTerminals(frequency);
  }

  // Bumpless PU controller initialization.
  **mTheta = std::arg(vsInit(0, 0));

  **mFrequencyPu = clampFinite(
      **mFrequencyRefPu + **mActivePowerDroopPu *
                              (**mActivePowerRefPu - **mFilteredActivePowerPu),
      mMinimumFrequencyPu, mMaximumFrequencyPu, **mFrequencyRefPu);

  **mOmegaPu = **mFrequencyPu;
  mPreviousOmegaPu = **mOmegaPu;

  **mVoltageDroopOutputPu =
      **mVoltageRefPu + **mReactivePowerDroopPu * (**mReactivePowerRefPu -
                                                   **mFilteredReactivePowerPu);

  mPreviousVoltageErrorPu = **mVoltageDroopOutputPu - **mVoltageMagnitudePu;

  **mVoltageCommandPu = std::abs(vsInit(0, 0)) / mBaseVoltagePhasePeak;

  **mVoltageIntegralStatePu = **mVoltageCommandPu - **mVoltageDroopOutputPu;

  **mVsrefPu = vsInit.real() / mBaseVoltagePhasePeak;

  updatePhysicalMirrors();
  **mVsref = **mVsrefPu * mBaseVoltagePhasePeak;

  mControlStateInitialized = true;

  SPDLOG_LOGGER_INFO(mSLog,
                     "GFM_Droop PU PF initialization:"
                     "\n  P={} pu ({} W)"
                     "\n  Q={} pu ({} var)"
                     "\n  V={} pu ({} V phase peak)"
                     "\n  f={} pu ({} Hz)"
                     "\n  theta={} rad"
                     "\n  source amplitude={} pu ({} V phase peak)",
                     **mElecActivePowerPu, **mElecActivePower,
                     **mElecReactivePowerPu, **mElecReactivePower,
                     **mVoltageMagnitudePu, **mVoltageMagnitude, **mFrequencyPu,
                     **mFrequency, **mTheta, **mVoltageCommandPu,
                     **mVoltageCommand);
}

void EMT::Ph3::GFM_Droop::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  requireBaseParameters();

  if (!(timeStep > 0.0) || !std::isfinite(timeStep))
    throw std::invalid_argument("GFM_Droop requires a positive time step");

  mTimeStep = timeStep;

  mPowerFilterAlpha =
      mPowerFilterTimeConstant > 0.0
          ? 1.0 - std::exp(-mTimeStep / mPowerFilterTimeConstant)
          : 1.0;

  const Real relativeOmegaDifference =
      std::abs(omega - mBaseOmega) / std::max<Real>(mBaseOmega, 1.0);

  if (relativeOmegaDifference > 1e-6) {
    SPDLOG_LOGGER_WARN(mSLog,
                       "Simulation omega={} rad/s differs from configured PU "
                       "omega_base={} rad/s",
                       omega, mBaseOmega);
  }

  SPDLOG_LOGGER_INFO(
      mSLog, "GFM_Droop discrete PU P/Q filter: dt={} s, tau={} s, alpha={}",
      mTimeStep, mPowerFilterTimeConstant, mPowerFilterAlpha);

  (void)leftVector;
}

void EMT::Ph3::GFM_Droop::updateMeasurements() {
  const Matrix &vPccSi = **mIntfVoltage;
  const Matrix &iInterfaceSi = **mIntfCurrent;
  Matrix &iPccSi = **mPccCurrent;
  Matrix &iPccPu = **mPccCurrentPu;

  for (UInt phase = 0; phase < 3; ++phase) {
    iPccSi(phase, 0) = -iInterfaceSi(phase, 0);
    iPccPu(phase, 0) = iPccSi(phase, 0) / mBaseCurrentPhasePeak;
  }

  const Real vaPu = vPccSi(0, 0) / mBaseVoltagePhasePeak;
  const Real vbPu = vPccSi(1, 0) / mBaseVoltagePhasePeak;
  const Real vcPu = vPccSi(2, 0) / mBaseVoltagePhasePeak;

  const Real iaPu = iPccPu(0, 0);
  const Real ibPu = iPccPu(1, 0);
  const Real icPu = iPccPu(2, 0);

  const Real vAlphaPu = SQRT_TWO_OVER_THREE * (vaPu - 0.5 * vbPu - 0.5 * vcPu);
  const Real vBetaPu =
      SQRT_TWO_OVER_THREE * SQRT_THREE_OVER_TWO * (vbPu - vcPu);

  const Real iAlphaPu = SQRT_TWO_OVER_THREE * (iaPu - 0.5 * ibPu - 0.5 * icPu);
  const Real iBetaPu =
      SQRT_TWO_OVER_THREE * SQRT_THREE_OVER_TWO * (ibPu - icPu);

  // Because V_base_peak * I_base_peak = 2/3*S_base, the factor 2/3
  // converts the power-invariant Clarke product to the three-phase PU base.
  **mElecActivePowerPu =
      TWO_OVER_THREE * (vAlphaPu * iAlphaPu + vBetaPu * iBetaPu);

  **mElecReactivePowerPu =
      TWO_OVER_THREE * (vBetaPu * iAlphaPu - vAlphaPu * iBetaPu);

  const Real voltageSquaredPu =
      TWO_OVER_THREE * (vaPu * vaPu + vbPu * vbPu + vcPu * vcPu);

  **mVoltageMagnitudePu = std::sqrt(std::max<Real>(0.0, voltageSquaredPu));

  updatePhysicalMirrors();
}

void EMT::Ph3::GFM_Droop::updatePowerMeasurementFilter(Int timeStepCount) {
  if (timeStepCount > 0) {
    **mFilteredActivePowerPu +=
        mPowerFilterAlpha * (**mElecActivePowerPu - **mFilteredActivePowerPu);

    **mFilteredReactivePowerPu +=
        mPowerFilterAlpha *
        (**mElecReactivePowerPu - **mFilteredReactivePowerPu);
  }

  **mFilteredActivePower = **mFilteredActivePowerPu * mBaseApparentPower;

  **mFilteredReactivePower = **mFilteredReactivePowerPu * mBaseApparentPower;
}

void EMT::Ph3::GFM_Droop::updateController(Int timeStepCount) {
  updateMeasurements();
  updatePowerMeasurementFilter(timeStepCount);

  const Real frequencyUnclampedPu =
      **mFrequencyRefPu +
      **mActivePowerDroopPu * (**mActivePowerRefPu - **mFilteredActivePowerPu);

  **mFrequencyPu = clampFinite(frequencyUnclampedPu, mMinimumFrequencyPu,
                               mMaximumFrequencyPu, **mFrequencyRefPu);

  **mOmegaPu = **mFrequencyPu;

  // theta_dot = omega_base * omega_pu.
  // The first scheduled solve is at t=dt with timeStepCount==0.
  **mTheta += 0.5 * mTimeStep * mBaseOmega * (mPreviousOmegaPu + **mOmegaPu);

  **mTheta = std::remainder(**mTheta, TWO_PI);

  **mVoltageDroopOutputPu =
      **mVoltageRefPu + **mReactivePowerDroopPu * (**mReactivePowerRefPu -
                                                   **mFilteredReactivePowerPu);

  const Real voltageErrorPu = **mVoltageDroopOutputPu - **mVoltageMagnitudePu;

  if (timeStepCount > 0) {
    **mVoltageIntegralStatePu += 0.5 * mTimeStep * **mVoltageIntegralGain *
                                 (mPreviousVoltageErrorPu + voltageErrorPu);
  }

  const Real voltageUnclampedPu =
      **mVoltageDroopOutputPu + **mVoltageIntegralStatePu;

  **mVoltageCommandPu = clampFinite(voltageUnclampedPu, mMinimumVoltagePu,
                                    mMaximumVoltagePu, **mVoltageRefPu);

  // Back-calculation anti-windup in PU.
  if (**mVoltageCommandPu != voltageUnclampedPu) {
    **mVoltageIntegralStatePu = **mVoltageCommandPu - **mVoltageDroopOutputPu;
  }

  mPreviousOmegaPu = **mOmegaPu;
  mPreviousVoltageErrorPu = voltageErrorPu;

  updatePhysicalMirrors();
  writeVoltageReference();
}

void EMT::Ph3::GFM_Droop::updateOpenLoop(Int timeStepCount) {
  updateMeasurements();
  updatePowerMeasurementFilter(timeStepCount);

  **mFrequencyPu = clampFinite(**mFrequencyRefPu, mMinimumFrequencyPu,
                               mMaximumFrequencyPu, **mFrequencyRefPu);

  **mOmegaPu = **mFrequencyPu;

  **mTheta += 0.5 * mTimeStep * mBaseOmega * (mPreviousOmegaPu + **mOmegaPu);

  **mTheta = std::remainder(**mTheta, TWO_PI);

  // Hold the PF-derived source-side amplitude in PU.
  **mVoltageCommandPu = clampFinite(**mVoltageCommandPu, mMinimumVoltagePu,
                                    mMaximumVoltagePu, **mVoltageRefPu);

  **mVoltageDroopOutputPu = **mVoltageCommandPu;

  **mVoltageIntegralStatePu = 0.0;

  mPreviousOmegaPu = **mOmegaPu;

  updatePhysicalMirrors();
  writeVoltageReference();
}

void EMT::Ph3::GFM_Droop::updatePhysicalMirrors() {
  **mElecActivePower = **mElecActivePowerPu * mBaseApparentPower;
  **mElecReactivePower = **mElecReactivePowerPu * mBaseApparentPower;

  **mFilteredActivePower = **mFilteredActivePowerPu * mBaseApparentPower;
  **mFilteredReactivePower = **mFilteredReactivePowerPu * mBaseApparentPower;

  **mVoltageMagnitude = **mVoltageMagnitudePu * mBaseVoltagePhasePeak;

  **mFrequency = **mFrequencyPu * mBaseFrequency;
  **mOmega = **mOmegaPu * mBaseOmega;

  **mVoltageDroopOutput = **mVoltageDroopOutputPu * mBaseVoltagePhasePeak;
  **mVoltageIntegralState = **mVoltageIntegralStatePu * mBaseVoltagePhasePeak;
  **mVoltageCommand = **mVoltageCommandPu * mBaseVoltagePhasePeak;
}

void EMT::Ph3::GFM_Droop::writeVoltageReference() {
  const Real theta = **mTheta;
  const Real amplitudePu = **mVoltageCommandPu;

  Matrix &referencePu = **mVsrefPu;
  Matrix &referenceSi = **mVsref;

  referencePu(0, 0) = amplitudePu * std::cos(theta);
  referencePu(1, 0) = amplitudePu * std::cos(theta - TWO_PI_OVER_THREE);
  referencePu(2, 0) = amplitudePu * std::cos(theta + TWO_PI_OVER_THREE);

  referenceSi = referencePu * mBaseVoltagePhasePeak;
}

void EMT::Ph3::GFM_Droop::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mIntfVoltage);
  prevStepDependencies.push_back(mIntfCurrent);

  // Canonical controller inputs.
  attributeDependencies.push_back(mActivePowerRefPu);
  attributeDependencies.push_back(mReactivePowerRefPu);
  attributeDependencies.push_back(mFrequencyRefPu);
  attributeDependencies.push_back(mVoltageRefPu);
  attributeDependencies.push_back(mActivePowerDroopPu);
  attributeDependencies.push_back(mReactivePowerDroopPu);
  attributeDependencies.push_back(mVoltageIntegralGain);

  // PU outputs and states.
  modifiedAttributes.push_back(mPccCurrentPu);
  modifiedAttributes.push_back(mElecActivePowerPu);
  modifiedAttributes.push_back(mElecReactivePowerPu);
  modifiedAttributes.push_back(mFilteredActivePowerPu);
  modifiedAttributes.push_back(mFilteredReactivePowerPu);
  modifiedAttributes.push_back(mVoltageMagnitudePu);
  modifiedAttributes.push_back(mFrequencyPu);
  modifiedAttributes.push_back(mOmegaPu);
  modifiedAttributes.push_back(mTheta);
  modifiedAttributes.push_back(mVoltageDroopOutputPu);
  modifiedAttributes.push_back(mVoltageIntegralStatePu);
  modifiedAttributes.push_back(mVoltageCommandPu);
  modifiedAttributes.push_back(mVsrefPu);

  // Existing SI mirrors.
  modifiedAttributes.push_back(mPccCurrent);
  modifiedAttributes.push_back(mElecActivePower);
  modifiedAttributes.push_back(mElecReactivePower);
  modifiedAttributes.push_back(mFilteredActivePower);
  modifiedAttributes.push_back(mFilteredReactivePower);
  modifiedAttributes.push_back(mVoltageMagnitude);
  modifiedAttributes.push_back(mFrequency);
  modifiedAttributes.push_back(mOmega);
  modifiedAttributes.push_back(mVoltageDroopOutput);
  modifiedAttributes.push_back(mVoltageIntegralState);
  modifiedAttributes.push_back(mVoltageCommand);
  modifiedAttributes.push_back(mVsref);

  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::GFM_Droop::mnaParentPreStep(Real time, Int timeStepCount) {
  if (!mControlStateInitialized) {
    throw std::runtime_error(
        "GFM_Droop PU control state was not initialized from power-flow "
        "data");
  }

  if (mWithControl)
    updateController(timeStepCount);
  else
    updateOpenLoop(timeStepCount);

  mSubCtrledVoltageSource->mVoltageRef->set(PEAK1PH_TO_RMS3PH * **mVsref);

  std::dynamic_pointer_cast<MNAInterface>(mSubCtrledVoltageSource)
      ->mnaPreStep(time, timeStepCount);

  mnaApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::GFM_Droop::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);

  (void)prevStepDependencies;
}

void EMT::Ph3::GFM_Droop::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateCurrent(**leftVector);
  mnaCompUpdateVoltage(**leftVector);

  (void)time;
  (void)timeStepCount;
}

void EMT::Ph3::GFM_Droop::mnaCompUpdateCurrent(const Matrix &leftVector) {
  if (mWithConnectionTransformer) {
    **mIntfCurrent = mConnectionTransformer->mIntfCurrent->get();
  } else {
    // Consumer-positive current entering the complete GFM:
    // i_intf = i_L - i_Rd.
    **mIntfCurrent =
        mSubInductorF->mIntfCurrent->get() - mSubResistorD->mIntfCurrent->get();
  }

  (void)leftVector;
}

void EMT::Ph3::GFM_Droop::mnaCompUpdateVoltage(const Matrix &leftVector) {
  for (const auto &virtualNode : mVirtualNodes)
    virtualNode->mnaUpdateVoltage(leftVector);

  (**mIntfVoltage)(0, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
  (**mIntfVoltage)(1, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
  (**mIntfVoltage)(2, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
}
