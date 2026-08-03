// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/EMT/EMT_Ph3_GFM_Siemens.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;

namespace {
constexpr Real TWO_PI = 2.0 * PI;
constexpr Real TWO_PI_OVER_THREE = 2.0 * PI / 3.0;
constexpr Real TWO_OVER_THREE = 2.0 / 3.0;

Real clampFinite(Real value, Real lower, Real upper, Real fallback) {
  if (!std::isfinite(value))
    return fallback;
  return std::clamp(value, lower, upper);
}

Real firstOrderAlpha(Real timeStep, Real timeConstant) {
  if (!(timeConstant > 0.0))
    return 1.0;
  return 1.0 - std::exp(-timeStep / timeConstant);
}

Matrix limitDqMagnitude(const Matrix &input, Real maximumMagnitude,
                        Bool &limited) {
  Matrix output = input;
  limited = false;

  if (!(maximumMagnitude > 0.0) || !std::isfinite(maximumMagnitude))
    return output;

  const Real magnitude = std::hypot(input(0, 0), input(1, 0));
  if (magnitude > maximumMagnitude && magnitude > 0.0) {
    output *= maximumMagnitude / magnitude;
    limited = true;
  }
  return output;
}

Bool finiteMatrix(const Matrix &value) { return value.allFinite(); }
} // namespace

EMT::Ph3::GFM_Siemens::GFM_Siemens(String uid, String name,
                                   Logger::Level logLevel)
    : CompositePowerComp<Real>(uid, name, true, true, logLevel),
      // MATLAB New VSI Model defaults.
      mActivePowerRefPu(mAttributes->create<Real>("P_ref_pu", 0.0)),
      mReactivePowerRefPu(mAttributes->create<Real>("Q_ref_pu", 0.0)),
      mFrequencyRefPu(mAttributes->create<Real>("f_ref_pu", 1.0)),
      mVoltageRefPu(mAttributes->create<Real>("V_ref_pu", 1.0)),
      mActivePowerDroopPu(mAttributes->create<Real>("k_p_pu", 0.02)),
      mReactivePowerDroopPu(mAttributes->create<Real>("k_q_pu", 0.0311)),
      mVoltageControllerKp(mAttributes->create<Real>("Kp_v", 0.52)),
      mVoltageControllerKi(mAttributes->create<Real>("Ki_v", 1.16)),
      mVoltageControllerFeedforward(mAttributes->create<Real>("Kff_i", 1.0)),
      mCurrentControllerKp(mAttributes->create<Real>("Kp_i", 0.74)),
      mCurrentControllerKi(mAttributes->create<Real>("Ki_i", 1.19)),
      mCurrentControllerFeedforward(mAttributes->create<Real>("Kff_v", 1.0)),
      mPccCurrent(mAttributes->create<Matrix>("i_pcc", Matrix::Zero(3, 1))),
      mPccCurrentPu(
          mAttributes->create<Matrix>("i_pcc_pu", Matrix::Zero(3, 1))),
      mFilterCurrent(mAttributes->create<Matrix>("i_lf", Matrix::Zero(3, 1))),
      mFilterCurrentPu(
          mAttributes->create<Matrix>("i_lf_pu", Matrix::Zero(3, 1))),
      mCapacitorCurrent(
          mAttributes->create<Matrix>("i_cf", Matrix::Zero(3, 1))),
      mCapacitorCurrentPu(
          mAttributes->create<Matrix>("i_cf_pu", Matrix::Zero(3, 1))),
      mPccVoltageDqPu(
          mAttributes->create<Matrix>("v_pcc_dq_pu", Matrix::Zero(2, 1))),
      mPccCurrentDqPu(
          mAttributes->create<Matrix>("i_pcc_dq_pu", Matrix::Zero(2, 1))),
      mFilterCurrentDqPu(
          mAttributes->create<Matrix>("i_lf_dq_pu", Matrix::Zero(2, 1))),
      mVoltageReferenceDqPu(
          mAttributes->create<Matrix>("v_ref_dq_pu", Matrix::Zero(2, 1))),
      mCurrentReferenceDqPu(
          mAttributes->create<Matrix>("i_ref_dq_pu", Matrix::Zero(2, 1))),
      mVoltageControllerIntegralPu(
          mAttributes->create<Matrix>("xi_v_dq_pu", Matrix::Zero(2, 1))),
      mCurrentControllerIntegralPu(
          mAttributes->create<Matrix>("xi_i_dq_pu", Matrix::Zero(2, 1))),
      mVoltageCommandPrePwmDqPu(mAttributes->create<Matrix>(
          "v_cmd_pre_pwm_dq_pu", Matrix::Zero(2, 1))),
      mVoltageCommandDqPu(
          mAttributes->create<Matrix>("v_cmd_dq_pu", Matrix::Zero(2, 1))),
      mElecActivePowerPu(mAttributes->create<Real>("P_elec_pu", 0.0)),
      mElecReactivePowerPu(mAttributes->create<Real>("Q_elec_pu", 0.0)),
      mFilteredActivePowerPu(mAttributes->create<Real>("P_filtered_pu", 0.0)),
      mFilteredReactivePowerPu(mAttributes->create<Real>("Q_filtered_pu", 0.0)),
      mVoltageMagnitudePu(mAttributes->create<Real>("V_magnitude_pu", 0.0)),
      mVoltageDroopReferencePu(mAttributes->create<Real>("V_droop_pu", 1.0)),
      mFrequencyDroopReferencePu(mAttributes->create<Real>("f_droop_pu", 1.0)),
      mFrequencyPu(mAttributes->create<Real>("frequency_pu", 1.0)),
      mTheta(mAttributes->create<Real>("theta", 0.0)),
      mVsrefPu(mAttributes->create<Matrix>("Vsref_pu", Matrix::Zero(3, 1))),
      mVsref(mAttributes->create<Matrix>("Vsref", Matrix::Zero(3, 1))),
      mVs(mAttributes->createDynamic<Matrix>("Vs")),
      mElecActivePower(mAttributes->create<Real>("P_elec", 0.0)),
      mElecReactivePower(mAttributes->create<Real>("Q_elec", 0.0)),
      mFilteredActivePower(mAttributes->create<Real>("P_filtered", 0.0)),
      mFilteredReactivePower(mAttributes->create<Real>("Q_filtered", 0.0)),
      mVoltageMagnitude(mAttributes->create<Real>("V_magnitude", 0.0)),
      mFrequency(mAttributes->create<Real>("frequency", 0.0)),
      mOmega(mAttributes->create<Real>("omega", 0.0)) {

  mPhaseType = PhaseType::ABC;
  setTerminalNumber(1);

  // v0: controlled-source terminal
  // v1: node between Rf and Lf
  // external terminal: PCC and shunt-capacitor node
  setVirtualNodeNumber(2);

  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(3, 1);

  mSubControlledVoltageSource =
      EMT::Ph3::VoltageSource::make(**mName + "_src", mLogLevel);
  mSubFilterResistor = EMT::Ph3::Resistor::make(**mName + "_resF", mLogLevel);
  mSubFilterInductor = EMT::Ph3::Inductor::make(**mName + "_indF", mLogLevel);
  mSubFilterCapacitor = EMT::Ph3::Capacitor::make(**mName + "_capF", mLogLevel);

  addMNASubComponent(mSubFilterResistor,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, false);
  addMNASubComponent(mSubFilterInductor,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
  addMNASubComponent(mSubFilterCapacitor,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  // The source pre-step must run after the complete controller update.
  addMNASubComponent(mSubControlledVoltageSource,
                     MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  mVs->setReference(mSubControlledVoltageSource->mIntfVoltage);

  SPDLOG_LOGGER_INFO(
      mSLog,
      "Create {} {}: per-unit cascaded GFM controller with filtered P/Q "
      "measurements, P-f/Q-V droop, SRF voltage/current PI control, "
      "cross-decoupling and optional PWM delay",
      type(), name);
}

void EMT::Ph3::GFM_Siemens::requireBaseParameters() const {
  if (!mBaseParametersSet) {
    throw std::runtime_error(
        "GFM_Siemens base parameters must be configured first");
  }
}

void EMT::Ph3::GFM_Siemens::requireFilterParameters() const {
  if (!mFilterParametersSet) {
    throw std::runtime_error(
        "GFM_Siemens RLC filter parameters must be configured first");
  }
}

void EMT::Ph3::GFM_Siemens::setBaseParameters(Real ratedApparentPower,
                                              Real ratedVoltageLineToLineRms,
                                              Real nominalFrequencyHz) {
  if (!(ratedApparentPower > 0.0) || !(ratedVoltageLineToLineRms > 0.0) ||
      !(nominalFrequencyHz > 0.0) || !std::isfinite(ratedApparentPower) ||
      !std::isfinite(ratedVoltageLineToLineRms) ||
      !std::isfinite(nominalFrequencyHz)) {
    throw std::invalid_argument(
        "GFM_Siemens requires positive finite S_base, V_base_LL_RMS and "
        "f_base");
  }
  if (mControlStateInitialized) {
    throw std::logic_error(
        "GFM_Siemens base parameters cannot change after initialization");
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

  **mFrequency = nominalFrequencyHz;
  **mOmega = mBaseOmega;

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Siemens base: S={} VA, V_LL={} V RMS, V_ph_peak={} V, "
      "I_ph_peak={} A, Z={} Ohm, f={} Hz",
      mBaseApparentPower, mBaseVoltageLineToLineRms, mBaseVoltagePhasePeak,
      mBaseCurrentPhasePeak, mBaseImpedance, mBaseFrequency);
}

void EMT::Ph3::GFM_Siemens::setReferencesPerUnit(
    Real frequencyReferencePu, Real voltageReferencePu,
    Real activePowerReferencePu, Real reactivePowerReferencePu) {
  requireBaseParameters();

  if (!(frequencyReferencePu > 0.0) || !(voltageReferencePu >= 0.0) ||
      !std::isfinite(frequencyReferencePu) ||
      !std::isfinite(voltageReferencePu) ||
      !std::isfinite(activePowerReferencePu) ||
      !std::isfinite(reactivePowerReferencePu)) {
    throw std::invalid_argument("Invalid GFM_Siemens PU references");
  }

  **mFrequencyRefPu = frequencyReferencePu;
  **mVoltageRefPu = voltageReferencePu;
  **mActivePowerRefPu = activePowerReferencePu;
  **mReactivePowerRefPu = reactivePowerReferencePu;
  mReferencesSet = true;
}

void EMT::Ph3::GFM_Siemens::setReferences(Real frequencyReferenceHz,
                                          Real voltageReferencePhasePeak,
                                          Real activePowerReference,
                                          Real reactivePowerReference) {
  requireBaseParameters();
  setReferencesPerUnit(frequencyReferenceHz / mBaseFrequency,
                       voltageReferencePhasePeak / mBaseVoltagePhasePeak,
                       activePowerReference / mBaseApparentPower,
                       reactivePowerReference / mBaseApparentPower);
}

void EMT::Ph3::GFM_Siemens::setDroopParametersPerUnit(
    Real activePowerDroopPu, Real reactivePowerDroopPu) {
  if (!std::isfinite(activePowerDroopPu) ||
      !std::isfinite(reactivePowerDroopPu)) {
    throw std::invalid_argument("GFM_Siemens droop gains must be finite");
  }
  **mActivePowerDroopPu = activePowerDroopPu;
  **mReactivePowerDroopPu = reactivePowerDroopPu;
}

void EMT::Ph3::GFM_Siemens::setPowerMeasurementFilterTimeConstants(
    Real activePowerTimeConstant, Real reactivePowerTimeConstant) {
  if (!(activePowerTimeConstant >= 0.0) ||
      !(reactivePowerTimeConstant >= 0.0) ||
      !std::isfinite(activePowerTimeConstant) ||
      !std::isfinite(reactivePowerTimeConstant)) {
    throw std::invalid_argument(
        "GFM_Siemens P/Q measurement-filter time constants must be finite "
        "and non-negative");
  }

  mActivePowerMeasurementTimeConstant = activePowerTimeConstant;
  mReactivePowerMeasurementTimeConstant = reactivePowerTimeConstant;
}

void EMT::Ph3::GFM_Siemens::setVoltageControllerParameters(
    Real proportionalGain, Real integralGainPerSecond,
    Real outputCurrentFeedforwardGain) {
  if (!(integralGainPerSecond >= 0.0) || !std::isfinite(proportionalGain) ||
      !std::isfinite(integralGainPerSecond) ||
      !std::isfinite(outputCurrentFeedforwardGain)) {
    throw std::invalid_argument(
        "Invalid GFM_Siemens voltage-controller parameters");
  }
  **mVoltageControllerKp = proportionalGain;
  **mVoltageControllerKi = integralGainPerSecond;
  **mVoltageControllerFeedforward = outputCurrentFeedforwardGain;
}

void EMT::Ph3::GFM_Siemens::setCurrentControllerParameters(
    Real proportionalGain, Real integralGainPerSecond,
    Real pccVoltageFeedforwardGain) {
  if (!(integralGainPerSecond >= 0.0) || !std::isfinite(proportionalGain) ||
      !std::isfinite(integralGainPerSecond) ||
      !std::isfinite(pccVoltageFeedforwardGain)) {
    throw std::invalid_argument(
        "Invalid GFM_Siemens current-controller parameters");
  }
  **mCurrentControllerKp = proportionalGain;
  **mCurrentControllerKi = integralGainPerSecond;
  **mCurrentControllerFeedforward = pccVoltageFeedforwardGain;
}

void EMT::Ph3::GFM_Siemens::setPwmDelayTimeConstant(Real timeConstantSeconds) {
  if (!(timeConstantSeconds >= 0.0) || !std::isfinite(timeConstantSeconds)) {
    throw std::invalid_argument(
        "GFM_Siemens PWM time constant must be finite and non-negative");
  }
  mPwmDelayTimeConstant = timeConstantSeconds;
}

void EMT::Ph3::GFM_Siemens::setControllerLimitsPerUnit(
    Real minimumFrequencyPu, Real maximumFrequencyPu,
    Real maximumCurrentReferencePu, Real maximumVoltageCommandPu) {
  if (!(minimumFrequencyPu >= 0.0) ||
      !(maximumFrequencyPu > minimumFrequencyPu) ||
      !(maximumCurrentReferencePu >= 0.0) ||
      !(maximumVoltageCommandPu >= 0.0) || !std::isfinite(minimumFrequencyPu) ||
      !std::isfinite(maximumFrequencyPu) ||
      !std::isfinite(maximumCurrentReferencePu) ||
      !std::isfinite(maximumVoltageCommandPu)) {
    throw std::invalid_argument("Invalid GFM_Siemens PU limits");
  }

  mMinimumFrequencyPu = minimumFrequencyPu;
  mMaximumFrequencyPu = maximumFrequencyPu;
  mMaximumCurrentReferencePu = maximumCurrentReferencePu;
  mMaximumVoltageCommandPu = maximumVoltageCommandPu;
}

void EMT::Ph3::GFM_Siemens::setFilterParameters(Real filterInductance,
                                                Real filterCapacitance,
                                                Real filterResistance) {
  requireBaseParameters();

  if (!(filterInductance > 0.0) || !(filterCapacitance > 0.0) ||
      !(filterResistance >= 0.0) || !std::isfinite(filterInductance) ||
      !std::isfinite(filterCapacitance) || !std::isfinite(filterResistance)) {
    throw std::invalid_argument(
        "GFM_Siemens filter requires Lf>0, Cf>0 and Rf>=0");
  }

  mFilterInductance = filterInductance;
  mFilterCapacitance = filterCapacitance;
  mFilterResistance = filterResistance;

  mFilterInductiveReactancePu = mBaseOmega * mFilterInductance / mBaseImpedance;
  mFilterCapacitiveSusceptancePu =
      mBaseOmega * mFilterCapacitance * mBaseImpedance;
  mFilterResistancePu = mFilterResistance / mBaseImpedance;

  Base::AvVoltageSourceInverterDQ::setFilterParameters(
      mFilterInductance, mFilterCapacitance, mFilterResistance, 0.0);

  mSubFilterResistor->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mFilterResistance));
  mSubFilterInductor->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mFilterInductance));
  mSubFilterCapacitor->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mFilterCapacitance));

  mFilterParametersSet = true;

  SPDLOG_LOGGER_INFO(mSLog,
                     "GFM_Siemens filter: Lf={} H, Cf={} F, Rf={} Ohm; "
                     "X_Lf={} pu, B_Cf={} pu, Rf={} pu",
                     mFilterInductance, mFilterCapacitance, mFilterResistance,
                     mFilterInductiveReactancePu,
                     mFilterCapacitiveSusceptancePu, mFilterResistancePu);
}

void EMT::Ph3::GFM_Siemens::setFilterParametersPerUnit(Real xLfPu, Real bCfPu,
                                                       Real rRfPu) {
  requireBaseParameters();
  if (!(xLfPu > 0.0) || !(bCfPu > 0.0) || !(rRfPu >= 0.0) ||
      !std::isfinite(xLfPu) || !std::isfinite(bCfPu) || !std::isfinite(rRfPu)) {
    throw std::invalid_argument("Invalid GFM_Siemens PU filter values");
  }

  setFilterParameters(xLfPu * mBaseImpedance / mBaseOmega,
                      bCfPu / (mBaseOmega * mBaseImpedance),
                      rRfPu * mBaseImpedance);
}

Matrix EMT::Ph3::GFM_Siemens::abcToDq(const Matrix &abc, Real theta) const {
  Matrix dq = Matrix::Zero(2, 1);
  const Real thetaA = theta;
  const Real thetaB = theta - TWO_PI_OVER_THREE;
  const Real thetaC = theta + TWO_PI_OVER_THREE;

  dq(0, 0) = TWO_OVER_THREE *
             (abc(0, 0) * std::cos(thetaA) + abc(1, 0) * std::cos(thetaB) +
              abc(2, 0) * std::cos(thetaC));

  dq(1, 0) = -TWO_OVER_THREE *
             (abc(0, 0) * std::sin(thetaA) + abc(1, 0) * std::sin(thetaB) +
              abc(2, 0) * std::sin(thetaC));
  return dq;
}

Matrix EMT::Ph3::GFM_Siemens::dqToAbc(const Matrix &dq, Real theta) const {
  Matrix abc = Matrix::Zero(3, 1);
  const Real thetaA = theta;
  const Real thetaB = theta - TWO_PI_OVER_THREE;
  const Real thetaC = theta + TWO_PI_OVER_THREE;

  abc(0, 0) = dq(0, 0) * std::cos(thetaA) - dq(1, 0) * std::sin(thetaA);
  abc(1, 0) = dq(0, 0) * std::cos(thetaB) - dq(1, 0) * std::sin(thetaB);
  abc(2, 0) = dq(0, 0) * std::cos(thetaC) - dq(1, 0) * std::sin(thetaC);
  return abc;
}

void EMT::Ph3::GFM_Siemens::initializeParentFromNodesAndTerminals(
    Real frequency) {
  requireBaseParameters();
  requireFilterParameters();

  if (!mReferencesSet) {
    // The MATLAB model uses 1 pu voltage/frequency and zero P/Q references.
    setReferencesPerUnit(1.0, 1.0, 0.0, 0.0);
  }

  MatrixComp pccVoltageComplex = MatrixComp::Zero(3, 1);
  MatrixComp interfaceCurrentComplex = MatrixComp::Zero(3, 1);

  const Real activePower = terminal(0)->singlePower().real();
  const Real reactivePower = terminal(0)->singlePower().imag();

  pccVoltageComplex(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  pccVoltageComplex(1, 0) = pccVoltageComplex(0, 0) * SHIFT_TO_PHASE_B;
  pccVoltageComplex(2, 0) = pccVoltageComplex(0, 0) * SHIFT_TO_PHASE_C;

  // Composite interface current is consumer-positive. The controller currents
  // are generator-positive.
  interfaceCurrentComplex(0, 0) =
      -std::conj(TWO_OVER_THREE * Complex(activePower, reactivePower) /
                 pccVoltageComplex(0, 0));
  interfaceCurrentComplex(1, 0) =
      interfaceCurrentComplex(0, 0) * SHIFT_TO_PHASE_B;
  interfaceCurrentComplex(2, 0) =
      interfaceCurrentComplex(0, 0) * SHIFT_TO_PHASE_C;

  const MatrixComp pccCurrentComplex = -interfaceCurrentComplex;
  const Real omegaInit = TWO_PI * frequency;

  // Current from PCC into the shunt capacitor and source-side filter current.
  const MatrixComp capacitorCurrentComplex =
      Complex(0.0, omegaInit * mFilterCapacitance) * pccVoltageComplex;
  const MatrixComp filterCurrentComplex =
      pccCurrentComplex + capacitorCurrentComplex;

  const Complex filterImpedance(mFilterResistance,
                                omegaInit * mFilterInductance);
  const MatrixComp sourceVoltageComplex =
      pccVoltageComplex + filterImpedance * filterCurrentComplex;
  const MatrixComp resistorInductorNodeVoltageComplex =
      sourceVoltageComplex - mFilterResistance * filterCurrentComplex;

  mVirtualNodes[0]->setInitialVoltage(PEAK1PH_TO_RMS3PH * sourceVoltageComplex);
  mVirtualNodes[1]->setInitialVoltage(PEAK1PH_TO_RMS3PH *
                                      resistorInductorNodeVoltageComplex);

  mSubControlledVoltageSource->setParameters(mVirtualNodes[0]->initialVoltage(),
                                             0.0);
  mSubControlledVoltageSource->connect(
      {CPS::SimNode<Real>::GND, mVirtualNodes[0]});

  mSubFilterResistor->connect({mVirtualNodes[0], mVirtualNodes[1]});
  mSubFilterInductor->connect({mVirtualNodes[1], mTerminals[0]->node()});

  mSubFilterCapacitor->connect(
      {mTerminals[0]->node(), CPS::SimNode<Real>::GND});
  for (const auto &subcomponent : mSubComponents) {
    subcomponent->initialize(mFrequencies);
    subcomponent->initializeFromNodesAndTerminals(frequency);
  }

  **mIntfVoltage = pccVoltageComplex.real();
  **mIntfCurrent = interfaceCurrentComplex.real();
  **mPccCurrent = pccCurrentComplex.real();
  **mFilterCurrent = filterCurrentComplex.real();
  **mCapacitorCurrent = capacitorCurrentComplex.real();

  **mPccCurrentPu = **mPccCurrent / mBaseCurrentPhasePeak;
  **mFilterCurrentPu = **mFilterCurrent / mBaseCurrentPhasePeak;
  **mCapacitorCurrentPu = **mCapacitorCurrent / mBaseCurrentPhasePeak;

  // Choose the internal angle so the initialized source voltage lies on d.
  **mTheta = std::arg(sourceVoltageComplex(0, 0));

  Matrix pccVoltagePu = pccVoltageComplex.real() / mBaseVoltagePhasePeak;
  **mPccVoltageDqPu = abcToDq(pccVoltagePu, **mTheta);
  **mPccCurrentDqPu = abcToDq(**mPccCurrentPu, **mTheta);
  **mFilterCurrentDqPu = abcToDq(**mFilterCurrentPu, **mTheta);

  **mElecActivePowerPu = activePower / mBaseApparentPower;
  **mElecReactivePowerPu = reactivePower / mBaseApparentPower;

  // Initialize the measurement filters exactly at the power-flow operating
  // point to avoid artificial P/Q filter transients at t=0.
  **mFilteredActivePowerPu = **mElecActivePowerPu;
  **mFilteredReactivePowerPu = **mElecReactivePowerPu;

  **mVoltageMagnitudePu =
      std::hypot((**mPccVoltageDqPu)(0, 0), (**mPccVoltageDqPu)(1, 0));

  **mFrequencyDroopReferencePu =
      **mFrequencyRefPu +
      **mActivePowerDroopPu * (**mActivePowerRefPu - **mFilteredActivePowerPu);
  **mFrequencyPu =
      clampFinite(**mFrequencyDroopReferencePu, mMinimumFrequencyPu,
                  mMaximumFrequencyPu, **mFrequencyRefPu);
  mPreviousFrequencyPu = **mFrequencyPu;

  **mVoltageDroopReferencePu =
      **mVoltageRefPu + **mReactivePowerDroopPu * (**mReactivePowerRefPu -
                                                   **mFilteredReactivePowerPu);

  **mVoltageReferenceDqPu = Matrix::Zero(2, 1);
  (**mVoltageReferenceDqPu)(0, 0) = **mVoltageDroopReferencePu;

  const Matrix voltageError = **mVoltageReferenceDqPu - **mPccVoltageDqPu;
  mPreviousVoltageErrorPu = voltageError;

  // Bumpless voltage-controller initialization: force i_ref=i_Lf.
  Matrix voltageControllerIntegral = Matrix::Zero(2, 1);
  const Real omegaPu = **mFrequencyPu;
  const Real vd = (**mPccVoltageDqPu)(0, 0);
  const Real vq = (**mPccVoltageDqPu)(1, 0);
  const Real iod = (**mPccCurrentDqPu)(0, 0);
  const Real ioq = (**mPccCurrentDqPu)(1, 0);

  voltageControllerIntegral(0, 0) =
      (**mFilterCurrentDqPu)(0, 0) -
      (**mVoltageControllerKp * voltageError(0, 0) +
       **mVoltageControllerFeedforward * iod -
       omegaPu * mFilterCapacitiveSusceptancePu * vq);

  voltageControllerIntegral(1, 0) =
      (**mFilterCurrentDqPu)(1, 0) -
      (**mVoltageControllerKp * voltageError(1, 0) +
       **mVoltageControllerFeedforward * ioq +
       omegaPu * mFilterCapacitiveSusceptancePu * vd);

  **mVoltageControllerIntegralPu = voltageControllerIntegral;
  **mCurrentReferenceDqPu = **mFilterCurrentDqPu;

  const Matrix currentError = **mCurrentReferenceDqPu - **mFilterCurrentDqPu;
  mPreviousCurrentErrorPu = currentError;

  const Matrix sourceVoltagePu =
      sourceVoltageComplex.real() / mBaseVoltagePhasePeak;
  const Matrix sourceVoltageDqPu = abcToDq(sourceVoltagePu, **mTheta);

  // Bumpless current-controller initialization: force v_cmd=v_source.
  Matrix currentControllerIntegral = Matrix::Zero(2, 1);
  const Real ild = (**mFilterCurrentDqPu)(0, 0);
  const Real ilq = (**mFilterCurrentDqPu)(1, 0);

  currentControllerIntegral(0, 0) =
      sourceVoltageDqPu(0, 0) - (**mCurrentControllerKp * currentError(0, 0) +
                                 **mCurrentControllerFeedforward * vd -
                                 omegaPu * mFilterInductiveReactancePu * ilq);

  currentControllerIntegral(1, 0) =
      sourceVoltageDqPu(1, 0) - (**mCurrentControllerKp * currentError(1, 0) +
                                 **mCurrentControllerFeedforward * vq +
                                 omegaPu * mFilterInductiveReactancePu * ild);

  **mCurrentControllerIntegralPu = currentControllerIntegral;
  **mVoltageCommandPrePwmDqPu = sourceVoltageDqPu;
  **mVoltageCommandDqPu = sourceVoltageDqPu;

  updatePhysicalMirrors();
  writeVoltageReference();
  mControlStateInitialized = true;

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Siemens PF initialization: P={} pu, Q={} pu, V={} pu, "
      "f={} pu, theta={} rad, source_dq=[{}, {}] pu",
      **mFilteredActivePowerPu, **mFilteredReactivePowerPu,
      **mVoltageMagnitudePu, **mFrequencyPu, **mTheta, sourceVoltageDqPu(0, 0),
      sourceVoltageDqPu(1, 0));
}

void EMT::Ph3::GFM_Siemens::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  requireBaseParameters();
  requireFilterParameters();

  if (!(timeStep > 0.0) || !std::isfinite(timeStep)) {
    throw std::invalid_argument(
        "GFM_Siemens requires a positive finite time step");
  }
  mTimeStep = timeStep;

  const Real relativeOmegaDifference =
      std::abs(omega - mBaseOmega) / std::max<Real>(mBaseOmega, 1.0);
  if (relativeOmegaDifference > 1.0e-6) {
    SPDLOG_LOGGER_WARN(
        mSLog,
        "Simulation omega={} rad/s differs from configured omega_base={} "
        "rad/s",
        omega, mBaseOmega);
  }

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Siemens P/Q measurement filters: dt={} s, tau_P={} s, "
      "tau_Q={} s, alpha_P={}, alpha_Q={}",
      mTimeStep, mActivePowerMeasurementTimeConstant,
      mReactivePowerMeasurementTimeConstant,
      firstOrderAlpha(mTimeStep, mActivePowerMeasurementTimeConstant),
      firstOrderAlpha(mTimeStep, mReactivePowerMeasurementTimeConstant));

  (void)leftVector;
}

void EMT::Ph3::GFM_Siemens::updateMeasurements() {
  **mPccCurrent = -**mIntfCurrent;

  // EMT branch-current orientation is terminal 1 -> terminal 0. The inductor
  // is connected source-side -> PCC and the capacitor PCC -> ground, hence
  // both signs are inverted for the generator-positive controller convention.
  **mFilterCurrent = -mSubFilterInductor->mIntfCurrent->get();
  **mCapacitorCurrent = -mSubFilterCapacitor->mIntfCurrent->get();

  **mPccCurrentPu = **mPccCurrent / mBaseCurrentPhasePeak;
  **mFilterCurrentPu = **mFilterCurrent / mBaseCurrentPhasePeak;
  **mCapacitorCurrentPu = **mCapacitorCurrent / mBaseCurrentPhasePeak;

  const Matrix pccVoltagePu = **mIntfVoltage / mBaseVoltagePhasePeak;
  **mPccVoltageDqPu = abcToDq(pccVoltagePu, **mTheta);
  **mPccCurrentDqPu = abcToDq(**mPccCurrentPu, **mTheta);
  **mFilterCurrentDqPu = abcToDq(**mFilterCurrentPu, **mTheta);

  const Real vd = (**mPccVoltageDqPu)(0, 0);
  const Real vq = (**mPccVoltageDqPu)(1, 0);
  const Real id = (**mPccCurrentDqPu)(0, 0);
  const Real iq = (**mPccCurrentDqPu)(1, 0);

  **mElecActivePowerPu = vd * id + vq * iq;
  **mElecReactivePowerPu = vq * id - vd * iq;
  **mVoltageMagnitudePu = std::hypot(vd, vq);

  updatePhysicalMirrors();
}

void EMT::Ph3::GFM_Siemens::updatePowerMeasurementFilters() {
  const Real activePowerAlpha =
      firstOrderAlpha(mTimeStep, mActivePowerMeasurementTimeConstant);
  const Real reactivePowerAlpha =
      firstOrderAlpha(mTimeStep, mReactivePowerMeasurementTimeConstant);

  **mFilteredActivePowerPu +=
      activePowerAlpha * (**mElecActivePowerPu - **mFilteredActivePowerPu);

  **mFilteredReactivePowerPu +=
      reactivePowerAlpha *
      (**mElecReactivePowerPu - **mFilteredReactivePowerPu);
}

void EMT::Ph3::GFM_Siemens::updateController(Int timeStepCount) {
  updateMeasurements();
  updatePowerMeasurementFilters();

  // ---------------------------------------------------------------------
  // Outer P-f and Q-V droop driven by filtered power measurements.
  // ---------------------------------------------------------------------
  **mFrequencyDroopReferencePu =
      **mFrequencyRefPu +
      **mActivePowerDroopPu * (**mActivePowerRefPu - **mFilteredActivePowerPu);

  **mFrequencyPu =
      clampFinite(**mFrequencyDroopReferencePu, mMinimumFrequencyPu,
                  mMaximumFrequencyPu, **mFrequencyRefPu);

  **mVoltageDroopReferencePu =
      **mVoltageRefPu + **mReactivePowerDroopPu * (**mReactivePowerRefPu -
                                                   **mFilteredReactivePowerPu);

  **mVoltageReferenceDqPu = Matrix::Zero(2, 1);
  (**mVoltageReferenceDqPu)(0, 0) = **mVoltageDroopReferencePu;

  // ---------------------------------------------------------------------
  // SRF capacitor-voltage controller.
  //
  // iLf_ref,d = PI_v,d + Kff_i*iPCC,d - omega*B_Cf*vPCC,q
  // iLf_ref,q = PI_v,q + Kff_i*iPCC,q + omega*B_Cf*vPCC,d
  // ---------------------------------------------------------------------
  const Matrix voltageError = **mVoltageReferenceDqPu - **mPccVoltageDqPu;

  if (timeStepCount > 0) {
    **mVoltageControllerIntegralPu += 0.5 * mTimeStep * **mVoltageControllerKi *
                                      (mPreviousVoltageErrorPu + voltageError);
  }

  const Real omegaPu = **mFrequencyPu;
  const Real vd = (**mPccVoltageDqPu)(0, 0);
  const Real vq = (**mPccVoltageDqPu)(1, 0);
  const Real iod = (**mPccCurrentDqPu)(0, 0);
  const Real ioq = (**mPccCurrentDqPu)(1, 0);

  Matrix currentReferenceUnclamped = Matrix::Zero(2, 1);
  currentReferenceUnclamped(0, 0) =
      **mVoltageControllerKp * voltageError(0, 0) +
      (**mVoltageControllerIntegralPu)(0, 0) +
      **mVoltageControllerFeedforward * iod -
      omegaPu * mFilterCapacitiveSusceptancePu * vq;

  currentReferenceUnclamped(1, 0) =
      **mVoltageControllerKp * voltageError(1, 0) +
      (**mVoltageControllerIntegralPu)(1, 0) +
      **mVoltageControllerFeedforward * ioq +
      omegaPu * mFilterCapacitiveSusceptancePu * vd;

  Bool currentReferenceLimited = false;
  **mCurrentReferenceDqPu =
      limitDqMagnitude(currentReferenceUnclamped, mMaximumCurrentReferencePu,
                       currentReferenceLimited);
  if (currentReferenceLimited) {
    **mVoltageControllerIntegralPu +=
        **mCurrentReferenceDqPu - currentReferenceUnclamped;
  }

  // ---------------------------------------------------------------------
  // SRF filter-current controller.
  //
  // vSrc_ref,d = PI_i,d + Kff_v*vPCC,d - omega*X_Lf*iLf,q
  // vSrc_ref,q = PI_i,q + Kff_v*vPCC,q + omega*X_Lf*iLf,d
  // ---------------------------------------------------------------------
  const Matrix currentError = **mCurrentReferenceDqPu - **mFilterCurrentDqPu;

  if (timeStepCount > 0) {
    **mCurrentControllerIntegralPu += 0.5 * mTimeStep * **mCurrentControllerKi *
                                      (mPreviousCurrentErrorPu + currentError);
  }

  const Real ild = (**mFilterCurrentDqPu)(0, 0);
  const Real ilq = (**mFilterCurrentDqPu)(1, 0);

  Matrix voltageCommandUnclamped = Matrix::Zero(2, 1);
  voltageCommandUnclamped(0, 0) = **mCurrentControllerKp * currentError(0, 0) +
                                  (**mCurrentControllerIntegralPu)(0, 0) +
                                  **mCurrentControllerFeedforward * vd -
                                  omegaPu * mFilterInductiveReactancePu * ilq;

  voltageCommandUnclamped(1, 0) = **mCurrentControllerKp * currentError(1, 0) +
                                  (**mCurrentControllerIntegralPu)(1, 0) +
                                  **mCurrentControllerFeedforward * vq +
                                  omegaPu * mFilterInductiveReactancePu * ild;

  Bool voltageCommandLimited = false;
  **mVoltageCommandPrePwmDqPu = limitDqMagnitude(
      voltageCommandUnclamped, mMaximumVoltageCommandPu, voltageCommandLimited);
  if (voltageCommandLimited) {
    **mCurrentControllerIntegralPu +=
        **mVoltageCommandPrePwmDqPu - voltageCommandUnclamped;
  }

  // Exact-ZOH discretization of the 1/(1+s*tau_pwm) block. This remains
  // stable when the EMT step is larger than the 1 us MATLAB time constant.
  const Real pwmAlpha = firstOrderAlpha(mTimeStep, mPwmDelayTimeConstant);
  **mVoltageCommandDqPu +=
      pwmAlpha * (**mVoltageCommandPrePwmDqPu - **mVoltageCommandDqPu);

  // Advance the internal oscillator for the next scheduled EMT solve.
  **mTheta +=
      0.5 * mTimeStep * mBaseOmega * (mPreviousFrequencyPu + **mFrequencyPu);
  **mTheta = std::remainder(**mTheta, TWO_PI);

  mPreviousFrequencyPu = **mFrequencyPu;
  mPreviousVoltageErrorPu = voltageError;
  mPreviousCurrentErrorPu = currentError;

  if (!finiteMatrix(**mCurrentReferenceDqPu) ||
      !finiteMatrix(**mVoltageCommandDqPu) || !std::isfinite(**mTheta)) {
    throw std::runtime_error(
        "GFM_Siemens controller produced a non-finite state");
  }

  updatePhysicalMirrors();
  writeVoltageReference();
}

void EMT::Ph3::GFM_Siemens::updateOpenLoop(Int timeStepCount) {
  updateMeasurements();
  updatePowerMeasurementFilters();

  **mFrequencyPu = clampFinite(**mFrequencyRefPu, mMinimumFrequencyPu,
                               mMaximumFrequencyPu, **mFrequencyRefPu);
  **mFrequencyDroopReferencePu = **mFrequencyPu;

  **mTheta +=
      0.5 * mTimeStep * mBaseOmega * (mPreviousFrequencyPu + **mFrequencyPu);
  **mTheta = std::remainder(**mTheta, TWO_PI);
  mPreviousFrequencyPu = **mFrequencyPu;

  updatePhysicalMirrors();
  writeVoltageReference();
  (void)timeStepCount;
}

void EMT::Ph3::GFM_Siemens::updatePhysicalMirrors() {
  **mElecActivePower = **mElecActivePowerPu * mBaseApparentPower;
  **mElecReactivePower = **mElecReactivePowerPu * mBaseApparentPower;
  **mFilteredActivePower = **mFilteredActivePowerPu * mBaseApparentPower;
  **mFilteredReactivePower = **mFilteredReactivePowerPu * mBaseApparentPower;
  **mVoltageMagnitude = **mVoltageMagnitudePu * mBaseVoltagePhasePeak;
  **mFrequency = **mFrequencyPu * mBaseFrequency;
  **mOmega = **mFrequencyPu * mBaseOmega;
}

void EMT::Ph3::GFM_Siemens::writeVoltageReference() {
  **mVsrefPu = dqToAbc(**mVoltageCommandDqPu, **mTheta);
  **mVsref = **mVsrefPu * mBaseVoltagePhasePeak;
}

void EMT::Ph3::GFM_Siemens::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mIntfVoltage);
  prevStepDependencies.push_back(mIntfCurrent);
  prevStepDependencies.push_back(mSubFilterInductor->mIntfCurrent);
  prevStepDependencies.push_back(mSubFilterCapacitor->mIntfCurrent);

  attributeDependencies.push_back(mActivePowerRefPu);
  attributeDependencies.push_back(mReactivePowerRefPu);
  attributeDependencies.push_back(mFrequencyRefPu);
  attributeDependencies.push_back(mVoltageRefPu);
  attributeDependencies.push_back(mActivePowerDroopPu);
  attributeDependencies.push_back(mReactivePowerDroopPu);
  attributeDependencies.push_back(mVoltageControllerKp);
  attributeDependencies.push_back(mVoltageControllerKi);
  attributeDependencies.push_back(mVoltageControllerFeedforward);
  attributeDependencies.push_back(mCurrentControllerKp);
  attributeDependencies.push_back(mCurrentControllerKi);
  attributeDependencies.push_back(mCurrentControllerFeedforward);

  modifiedAttributes.push_back(mPccCurrent);
  modifiedAttributes.push_back(mPccCurrentPu);
  modifiedAttributes.push_back(mFilterCurrent);
  modifiedAttributes.push_back(mFilterCurrentPu);
  modifiedAttributes.push_back(mCapacitorCurrent);
  modifiedAttributes.push_back(mCapacitorCurrentPu);
  modifiedAttributes.push_back(mPccVoltageDqPu);
  modifiedAttributes.push_back(mPccCurrentDqPu);
  modifiedAttributes.push_back(mFilterCurrentDqPu);
  modifiedAttributes.push_back(mVoltageReferenceDqPu);
  modifiedAttributes.push_back(mCurrentReferenceDqPu);
  modifiedAttributes.push_back(mVoltageControllerIntegralPu);
  modifiedAttributes.push_back(mCurrentControllerIntegralPu);
  modifiedAttributes.push_back(mVoltageCommandPrePwmDqPu);
  modifiedAttributes.push_back(mVoltageCommandDqPu);
  modifiedAttributes.push_back(mElecActivePowerPu);
  modifiedAttributes.push_back(mElecReactivePowerPu);
  modifiedAttributes.push_back(mFilteredActivePowerPu);
  modifiedAttributes.push_back(mFilteredReactivePowerPu);
  modifiedAttributes.push_back(mVoltageMagnitudePu);
  modifiedAttributes.push_back(mVoltageDroopReferencePu);
  modifiedAttributes.push_back(mFrequencyDroopReferencePu);
  modifiedAttributes.push_back(mFrequencyPu);
  modifiedAttributes.push_back(mTheta);
  modifiedAttributes.push_back(mVsrefPu);
  modifiedAttributes.push_back(mVsref);
  modifiedAttributes.push_back(mElecActivePower);
  modifiedAttributes.push_back(mElecReactivePower);
  modifiedAttributes.push_back(mFilteredActivePower);
  modifiedAttributes.push_back(mFilteredReactivePower);
  modifiedAttributes.push_back(mVoltageMagnitude);
  modifiedAttributes.push_back(mFrequency);
  modifiedAttributes.push_back(mOmega);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::GFM_Siemens::mnaParentPreStep(Real time, Int timeStepCount) {
  if (!mControlStateInitialized) {
    throw std::runtime_error(
        "GFM_Siemens controller was not initialized from power-flow data");
  }

  if (mWithControl)
    updateController(timeStepCount);
  else
    updateOpenLoop(timeStepCount);

  mSubControlledVoltageSource->mVoltageRef->set(PEAK1PH_TO_RMS3PH * **mVsref);

  std::dynamic_pointer_cast<MNAInterface>(mSubControlledVoltageSource)
      ->mnaPreStep(time, timeStepCount);

  mnaApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::GFM_Siemens::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
  (void)prevStepDependencies;
}

void EMT::Ph3::GFM_Siemens::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateCurrent(**leftVector);
  mnaCompUpdateVoltage(**leftVector);
  (void)time;
  (void)timeStepCount;
}

void EMT::Ph3::GFM_Siemens::mnaCompUpdateCurrent(const Matrix &leftVector) {
  // Branch-current orientation is terminal 1 -> terminal 0:
  //   i_L,intf = -i_Lf(generator-positive)
  //   i_C,intf = -i_C(PCC-to-ground)
  // Therefore i_intf = i_L,intf - i_C,intf = -i_PCC.
  **mIntfCurrent = mSubFilterInductor->mIntfCurrent->get() -
                   mSubFilterCapacitor->mIntfCurrent->get();
  (void)leftVector;
}

void EMT::Ph3::GFM_Siemens::mnaCompUpdateVoltage(const Matrix &leftVector) {
  for (const auto &virtualNode : mVirtualNodes)
    virtualNode->mnaUpdateVoltage(leftVector);

  (**mIntfVoltage)(0, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
  (**mIntfVoltage)(1, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
  (**mIntfVoltage)(2, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
}
