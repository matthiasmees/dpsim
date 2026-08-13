// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/EMT/EMT_Ph3_GFL_Siemens.h>

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
} // namespace

EMT::Ph3::GFL_Siemens::GFL_Siemens(String uid, String name,
                                   Logger::Level logLevel)
    : CompositePowerComp<Real>(uid, name, true, true, logLevel),
      mActivePowerRefPu(mAttributes->create<Real>("P_ref_pu", 0.0)),
      mReactivePowerRefPu(mAttributes->create<Real>("Q_ref_pu", 0.0)),
      mFrequencyRefPu(mAttributes->create<Real>("f_ref_pu", 1.0)),
      mVoltageRefPu(mAttributes->create<Real>("V_ref_pu", 1.0)),
      mFrequencyToActivePowerGainPu(mAttributes->create<Real>("K_fP_pu", 20.0)),
      mVoltageToReactivePowerGainPu(mAttributes->create<Real>("K_VQ_pu", 20.0)),
      mPllKp(mAttributes->create<Real>("Kp_pll", 0.449 / (9.1e-3 + 250.0e-6))),
      mPllKi(mAttributes->create<Real>("Ki_pll",
                                       250.0e-6 / (5.9 * (9.1e-3 + 250.0e-6)))),
      mCurrentControllerKp(mAttributes->create<Real>("Kp_i", 0.74 / 10.0)),
      mCurrentControllerKi(mAttributes->create<Real>("Ki_i", 1.19 / 10.0)),
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
      mFilteredPccVoltageDqPu(mAttributes->create<Matrix>(
          "v_pcc_filtered_dq_pu", Matrix::Zero(2, 1))),
      mPccCurrentDqPu(
          mAttributes->create<Matrix>("i_pcc_dq_pu", Matrix::Zero(2, 1))),
      mFilteredPccCurrentDqPu(mAttributes->create<Matrix>(
          "i_pcc_filtered_dq_pu", Matrix::Zero(2, 1))),
      mFilterCurrentDqPu(
          mAttributes->create<Matrix>("i_lf_dq_pu", Matrix::Zero(2, 1))),
      mCurrentReferenceDqPu(
          mAttributes->create<Matrix>("i_ref_dq_pu", Matrix::Zero(2, 1))),
      mCurrentErrorDqPu(
          mAttributes->create<Matrix>("e_i_dq_pu", Matrix::Zero(2, 1))),
      mCurrentControllerIntegralPu(
          mAttributes->create<Matrix>("xi_i_dq_pu", Matrix::Zero(2, 1))),
      mVoltageCommandDqPu(
          mAttributes->create<Matrix>("v_cmd_dq_pu", Matrix::Zero(2, 1))),
      mElecActivePowerPu(mAttributes->create<Real>("P_elec_pu", 0.0)),
      mElecReactivePowerPu(mAttributes->create<Real>("Q_elec_pu", 0.0)),
      mFilteredActivePowerPu(mAttributes->create<Real>("P_filtered_pu", 0.0)),
      mFilteredReactivePowerPu(mAttributes->create<Real>("Q_filtered_pu", 0.0)),
      mActivePowerCommandPu(mAttributes->create<Real>("P_command_pu", 0.0)),
      mReactivePowerCommandPu(mAttributes->create<Real>("Q_command_pu", 0.0)),
      mVoltageMagnitudePu(mAttributes->create<Real>("V_magnitude_pu", 0.0)),
      mPllVoltageErrorPu(mAttributes->create<Real>("pll_vq_error_pu", 0.0)),
      mPllIntegral(mAttributes->create<Real>("pll_integrator", 0.0)),
      mFrequencyPu(mAttributes->create<Real>("frequency_pu", 1.0)),
      mOmega(mAttributes->create<Real>("omega", TWO_PI * 50.0)),
      mTheta(mAttributes->create<Real>("theta", 0.0)),
      mVsrefPu(mAttributes->create<Matrix>("Vsref_pu", Matrix::Zero(3, 1))),
      mVsref(mAttributes->create<Matrix>("Vsref", Matrix::Zero(3, 1))),
      mVs(mAttributes->createDynamic<Matrix>("Vs")),
      mElecActivePower(mAttributes->create<Real>("P_elec", 0.0)),
      mElecReactivePower(mAttributes->create<Real>("Q_elec", 0.0)),
      mFilteredActivePower(mAttributes->create<Real>("P_filtered", 0.0)),
      mFilteredReactivePower(mAttributes->create<Real>("Q_filtered", 0.0)),
      mFrequency(mAttributes->create<Real>("frequency", 50.0)),
      mVoltageMagnitude(mAttributes->create<Real>("V_magnitude", 0.0)) {

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
      "Create {} {}: per-unit CSI-derived GFL with SRF PLL, inverse "
      "frequency/voltage droop, filtered P/Q and PCC-current measurements, "
      "filtered PCC-voltage feedforward, and SRF current control",
      type(), name);
}

void EMT::Ph3::GFL_Siemens::requireBaseParameters() const {
  if (!mBaseParametersSet) {
    throw std::runtime_error(
        "GFL_Siemens base parameters must be configured first");
  }
}

void EMT::Ph3::GFL_Siemens::requireFilterParameters() const {
  if (!mFilterParametersSet) {
    throw std::runtime_error(
        "GFL_Siemens RLC filter parameters must be configured first");
  }
}

void EMT::Ph3::GFL_Siemens::setBaseParameters(Real ratedApparentPower,
                                              Real ratedVoltageLineToLineRms,
                                              Real nominalFrequencyHz) {
  if (!(ratedApparentPower > 0.0) || !(ratedVoltageLineToLineRms > 0.0) ||
      !(nominalFrequencyHz > 0.0) || !std::isfinite(ratedApparentPower) ||
      !std::isfinite(ratedVoltageLineToLineRms) ||
      !std::isfinite(nominalFrequencyHz)) {
    throw std::invalid_argument(
        "GFL_Siemens requires positive finite S_base, V_base_LL_RMS and "
        "f_base");
  }
  if (mControlStateInitialized) {
    throw std::logic_error(
        "GFL_Siemens base parameters cannot change after initialization");
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
  **mFrequencyPu = 1.0;
  **mOmega = mBaseOmega;

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFL_Siemens base: S={} VA, V_LL={} V RMS, V_ph_peak={} V, "
      "I_ph_peak={} A, Z={} Ohm, f={} Hz",
      mBaseApparentPower, mBaseVoltageLineToLineRms, mBaseVoltagePhasePeak,
      mBaseCurrentPhasePeak, mBaseImpedance, mBaseFrequency);
}

void EMT::Ph3::GFL_Siemens::setReferencesPerUnit(
    Real frequencyReferencePu, Real voltageReferencePu,
    Real activePowerReferencePu, Real reactivePowerReferencePu) {
  requireBaseParameters();

  if (!(frequencyReferencePu > 0.0) || !(voltageReferencePu >= 0.0) ||
      !std::isfinite(frequencyReferencePu) ||
      !std::isfinite(voltageReferencePu) ||
      !std::isfinite(activePowerReferencePu) ||
      !std::isfinite(reactivePowerReferencePu)) {
    throw std::invalid_argument("Invalid GFL_Siemens PU references");
  }

  **mFrequencyRefPu = frequencyReferencePu;
  **mVoltageRefPu = voltageReferencePu;
  **mActivePowerRefPu = activePowerReferencePu;
  **mReactivePowerRefPu = reactivePowerReferencePu;
  mReferencesSet = true;
}

void EMT::Ph3::GFL_Siemens::setReferences(Real frequencyReferenceHz,
                                          Real voltageReferencePhasePeak,
                                          Real activePowerReference,
                                          Real reactivePowerReference) {
  requireBaseParameters();
  setReferencesPerUnit(frequencyReferenceHz / mBaseFrequency,
                       voltageReferencePhasePeak / mBaseVoltagePhasePeak,
                       activePowerReference / mBaseApparentPower,
                       reactivePowerReference / mBaseApparentPower);
}

void EMT::Ph3::GFL_Siemens::setDroopParametersPerUnit(
    Real frequencyToActivePowerGainPu, Real voltageToReactivePowerGainPu) {
  if (!std::isfinite(frequencyToActivePowerGainPu) ||
      !std::isfinite(voltageToReactivePowerGainPu)) {
    throw std::invalid_argument("GFL_Siemens droop gains must be finite");
  }
  **mFrequencyToActivePowerGainPu = frequencyToActivePowerGainPu;
  **mVoltageToReactivePowerGainPu = voltageToReactivePowerGainPu;
}

void EMT::Ph3::GFL_Siemens::setPllParameters(Real proportionalGain,
                                             Real integralGainPerSecond) {
  if (!std::isfinite(proportionalGain) || !(integralGainPerSecond >= 0.0) ||
      !std::isfinite(integralGainPerSecond)) {
    throw std::invalid_argument("Invalid GFL_Siemens PLL gains");
  }
  **mPllKp = proportionalGain;
  **mPllKi = integralGainPerSecond;
}

void EMT::Ph3::GFL_Siemens::setCurrentControllerParameters(
    Real proportionalGain, Real integralGainPerSecond,
    Real pccVoltageFeedforwardGain) {
  if (!std::isfinite(proportionalGain) || !(integralGainPerSecond >= 0.0) ||
      !std::isfinite(integralGainPerSecond) ||
      !std::isfinite(pccVoltageFeedforwardGain)) {
    throw std::invalid_argument(
        "Invalid GFL_Siemens current-controller parameters");
  }

  **mCurrentControllerKp = proportionalGain;
  **mCurrentControllerKi = integralGainPerSecond;
  **mCurrentControllerFeedforward = pccVoltageFeedforwardGain;
}

void EMT::Ph3::GFL_Siemens::setMeasurementFilterTimeConstants(
    Real activePowerTimeConstant, Real reactivePowerTimeConstant,
    Real currentTimeConstant, Real voltageFeedforwardTimeConstant) {
  if (!(activePowerTimeConstant >= 0.0) ||
      !(reactivePowerTimeConstant >= 0.0) || !(currentTimeConstant >= 0.0) ||
      !(voltageFeedforwardTimeConstant >= 0.0) ||
      !std::isfinite(activePowerTimeConstant) ||
      !std::isfinite(reactivePowerTimeConstant) ||
      !std::isfinite(currentTimeConstant) ||
      !std::isfinite(voltageFeedforwardTimeConstant)) {
    throw std::invalid_argument(
        "GFL_Siemens measurement-filter time constants must be finite and "
        "non-negative");
  }

  mActivePowerMeasurementTimeConstant = activePowerTimeConstant;
  mReactivePowerMeasurementTimeConstant = reactivePowerTimeConstant;
  mCurrentMeasurementTimeConstant = currentTimeConstant;
  mVoltageFeedforwardMeasurementTimeConstant = voltageFeedforwardTimeConstant;
}

void EMT::Ph3::GFL_Siemens::setControllerLimitsPerUnit(
    Real minimumFrequencyPu, Real maximumFrequencyPu,
    Real maximumCurrentReferencePu, Real maximumVoltageCommandPu,
    Real minimumVoltageForCurrentReferencePu) {
  if (!(minimumFrequencyPu >= 0.0) ||
      !(maximumFrequencyPu > minimumFrequencyPu) ||
      !(maximumCurrentReferencePu >= 0.0) ||
      !(maximumVoltageCommandPu >= 0.0) ||
      !(minimumVoltageForCurrentReferencePu > 0.0) ||
      !std::isfinite(minimumFrequencyPu) ||
      !std::isfinite(maximumFrequencyPu) ||
      !std::isfinite(maximumCurrentReferencePu) ||
      !std::isfinite(maximumVoltageCommandPu) ||
      !std::isfinite(minimumVoltageForCurrentReferencePu)) {
    throw std::invalid_argument("Invalid GFL_Siemens PU limits");
  }

  mMinimumFrequencyPu = minimumFrequencyPu;
  mMaximumFrequencyPu = maximumFrequencyPu;
  mMaximumCurrentReferencePu = maximumCurrentReferencePu;
  mMaximumVoltageCommandPu = maximumVoltageCommandPu;
  mMinimumVoltageForCurrentReferencePu = minimumVoltageForCurrentReferencePu;
}

void EMT::Ph3::GFL_Siemens::setFilterParameters(Real filterInductance,
                                                Real filterCapacitance,
                                                Real filterResistance) {
  requireBaseParameters();

  if (!(filterInductance > 0.0) || !(filterCapacitance > 0.0) ||
      !(filterResistance >= 0.0) || !std::isfinite(filterInductance) ||
      !std::isfinite(filterCapacitance) || !std::isfinite(filterResistance)) {
    throw std::invalid_argument(
        "GFL_Siemens filter requires Lf>0, Cf>0 and Rf>=0");
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
                     "GFL_Siemens filter: Lf={} H, Cf={} F, Rf={} Ohm; "
                     "X_Lf={} pu, B_Cf={} pu, Rf={} pu",
                     mFilterInductance, mFilterCapacitance, mFilterResistance,
                     mFilterInductiveReactancePu,
                     mFilterCapacitiveSusceptancePu, mFilterResistancePu);
}

void EMT::Ph3::GFL_Siemens::setFilterParametersPerUnit(Real xLfPu, Real bCfPu,
                                                       Real rRfPu) {
  requireBaseParameters();
  if (!(xLfPu > 0.0) || !(bCfPu > 0.0) || !(rRfPu >= 0.0) ||
      !std::isfinite(xLfPu) || !std::isfinite(bCfPu) || !std::isfinite(rRfPu)) {
    throw std::invalid_argument("Invalid GFL_Siemens PU filter values");
  }

  setFilterParameters(xLfPu * mBaseImpedance / mBaseOmega,
                      bCfPu / (mBaseOmega * mBaseImpedance),
                      rRfPu * mBaseImpedance);
}

Matrix EMT::Ph3::GFL_Siemens::abcToDq(const Matrix &abc, Real theta) const {
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

Matrix EMT::Ph3::GFL_Siemens::dqToAbc(const Matrix &dq, Real theta) const {
  Matrix abc = Matrix::Zero(3, 1);
  const Real thetaA = theta;
  const Real thetaB = theta - TWO_PI_OVER_THREE;
  const Real thetaC = theta + TWO_PI_OVER_THREE;

  abc(0, 0) = dq(0, 0) * std::cos(thetaA) - dq(1, 0) * std::sin(thetaA);
  abc(1, 0) = dq(0, 0) * std::cos(thetaB) - dq(1, 0) * std::sin(thetaB);
  abc(2, 0) = dq(0, 0) * std::cos(thetaC) - dq(1, 0) * std::sin(thetaC);
  return abc;
}

void EMT::Ph3::GFL_Siemens::initializeParentFromNodesAndTerminals(
    Real frequency) {
  requireBaseParameters();
  requireFilterParameters();

  MatrixComp pccVoltageComplex = MatrixComp::Zero(3, 1);
  MatrixComp interfaceCurrentComplex = MatrixComp::Zero(3, 1);

  const Real activePower = terminal(0)->singlePower().real();
  const Real reactivePower = terminal(0)->singlePower().imag();

  pccVoltageComplex(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  pccVoltageComplex(1, 0) = pccVoltageComplex(0, 0) * SHIFT_TO_PHASE_B;
  pccVoltageComplex(2, 0) = pccVoltageComplex(0, 0) * SHIFT_TO_PHASE_C;

  // Composite interface current is consumer-positive. Controller currents are
  // generator-positive.
  interfaceCurrentComplex(0, 0) =
      -std::conj(TWO_OVER_THREE * Complex(activePower, reactivePower) /
                 pccVoltageComplex(0, 0));
  interfaceCurrentComplex(1, 0) =
      interfaceCurrentComplex(0, 0) * SHIFT_TO_PHASE_B;
  interfaceCurrentComplex(2, 0) =
      interfaceCurrentComplex(0, 0) * SHIFT_TO_PHASE_C;

  const MatrixComp pccCurrentComplex = -interfaceCurrentComplex;
  const Real omegaInit = TWO_PI * frequency;

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

  const MatrixComp kclResidual =
      filterCurrentComplex - pccCurrentComplex - capacitorCurrentComplex;

  const MatrixComp kvlResidual = sourceVoltageComplex - pccVoltageComplex -
                                 filterImpedance * filterCurrentComplex;

  SPDLOG_INFO(
      "\n========== GFL ANALYTICAL PF -> EMT INITIALIZATION =========="
      "\nTerminal operating point:"
      "\n  P = {} W"
      "\n  Q = {} var"
      "\n  frequency = {} Hz"
      "\nPhase-A phasors in PEAK convention:"
      "\n  Vpcc = [{}, {}] V, |V| = {}, angle = {} deg"
      "\n  Ipcc = [{}, {}] A"
      "\n  Icf  = [{}, {}] A"
      "\n  Ilf  = [{}, {}] A"
      "\n  Vsrc = [{}, {}] V"
      "\n  Vrl   = [{}, {}] V"
      "\nAnalytical residuals:"
      "\n  KCL residual = [{}, {}] A, |res| = {}"
      "\n  KVL residual = [{}, {}] V, |res| = {}"
      "\n============================================================",
      activePower, reactivePower, frequency,

      pccVoltageComplex(0, 0).real(), pccVoltageComplex(0, 0).imag(),
      std::abs(pccVoltageComplex(0, 0)),
      std::arg(pccVoltageComplex(0, 0)) * 180.0 / PI,

      pccCurrentComplex(0, 0).real(), pccCurrentComplex(0, 0).imag(),

      capacitorCurrentComplex(0, 0).real(),
      capacitorCurrentComplex(0, 0).imag(),

      filterCurrentComplex(0, 0).real(), filterCurrentComplex(0, 0).imag(),

      sourceVoltageComplex(0, 0).real(), sourceVoltageComplex(0, 0).imag(),

      resistorInductorNodeVoltageComplex(0, 0).real(),
      resistorInductorNodeVoltageComplex(0, 0).imag(),

      kclResidual(0, 0).real(), kclResidual(0, 0).imag(),
      std::abs(kclResidual(0, 0)),

      kvlResidual(0, 0).real(), kvlResidual(0, 0).imag(),
      std::abs(kvlResidual(0, 0)));

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

  // Controller convention: positive current flows from the converter source
  // toward the PCC. The subcomponent interface orientations are opposite.
  const Matrix actualFilterCurrentAbc =
      -mSubFilterInductor->mIntfCurrent->get();

  const Matrix actualCapacitorCurrentAbc =
      -mSubFilterCapacitor->mIntfCurrent->get();

  const Matrix expectedFilterCurrentAbc = filterCurrentComplex.real();

  const Matrix expectedCapacitorCurrentAbc = capacitorCurrentComplex.real();

  const Matrix filterCurrentErrorAbc =
      actualFilterCurrentAbc - expectedFilterCurrentAbc;

  const Matrix capacitorCurrentErrorAbc =
      actualCapacitorCurrentAbc - expectedCapacitorCurrentAbc;

  **mIntfVoltage = pccVoltageComplex.real();
  **mIntfCurrent = interfaceCurrentComplex.real();
  **mPccCurrent = pccCurrentComplex.real();
  **mFilterCurrent = filterCurrentComplex.real();
  **mCapacitorCurrent = capacitorCurrentComplex.real();

  **mPccCurrentPu = **mPccCurrent / mBaseCurrentPhasePeak;
  **mFilterCurrentPu = **mFilterCurrent / mBaseCurrentPhasePeak;
  **mCapacitorCurrentPu = **mCapacitorCurrent / mBaseCurrentPhasePeak;

  // Use a d-axis-aligned implementation. This is physically equivalent to
  // the source model's "90 degrees behind phase A" convention when the same
  // offset is used consistently in forward and inverse transforms.
  **mTheta = std::arg(pccVoltageComplex(0, 0));

  const Matrix pccVoltagePu = pccVoltageComplex.real() / mBaseVoltagePhasePeak;
  **mPccVoltageDqPu = abcToDq(pccVoltagePu, **mTheta);
  **mFilteredPccVoltageDqPu = **mPccVoltageDqPu;
  **mPccCurrentDqPu = abcToDq(**mPccCurrentPu, **mTheta);
  **mFilteredPccCurrentDqPu = **mPccCurrentDqPu;
  **mFilterCurrentDqPu = abcToDq(**mFilterCurrentPu, **mTheta);

  **mElecActivePowerPu = activePower / mBaseApparentPower;
  **mElecReactivePowerPu = reactivePower / mBaseApparentPower;
  **mFilteredActivePowerPu = **mElecActivePowerPu;
  **mFilteredReactivePowerPu = **mElecReactivePowerPu;

  const Real vd = (**mPccVoltageDqPu)(0, 0);
  const Real vq = (**mPccVoltageDqPu)(1, 0);
  **mVoltageMagnitudePu = std::hypot(vd, vq);

  // When no command was supplied explicitly, adopt the PF operating point for
  // bumpless initialization. Explicit references preserve the CSI defaults or
  // user commands.
  if (!mReferencesSet) {
    setReferencesPerUnit(1.0, **mVoltageMagnitudePu, **mElecActivePowerPu,
                         **mElecReactivePowerPu);
  }

  **mPllVoltageErrorPu = vq;

  // Initialize the PLL exactly at the requested steady-state frequency.
  // Even when the Park transform leaves a small numerical v_q residual,
  // Kp*v_q + xi_pll must not create an artificial frequency step.
  **mFrequencyPu = **mFrequencyRefPu;
  **mOmega = **mFrequencyPu * mBaseOmega;
  **mPllIntegral = **mOmega - mBaseOmega - **mPllKp * **mPllVoltageErrorPu;

  mPreviousPllVoltageErrorPu = **mPllVoltageErrorPu;
  mPreviousFrequencyPu = **mFrequencyPu;

  **mActivePowerCommandPu =
      **mActivePowerRefPu -
      **mFrequencyToActivePowerGainPu * (**mFrequencyPu - **mFrequencyRefPu);
  **mReactivePowerCommandPu =
      **mReactivePowerRefPu - **mVoltageToReactivePowerGainPu *
                                  (**mVoltageMagnitudePu - **mVoltageRefPu);

  const Real denominator =
      std::max(vd * vd + vq * vq, mMinimumVoltageForCurrentReferencePu *
                                      mMinimumVoltageForCurrentReferencePu);

  Matrix currentReference = Matrix::Zero(2, 1);
  currentReference(0, 0) =
      (**mActivePowerCommandPu * vd + **mReactivePowerCommandPu * vq) /
      denominator;
  currentReference(1, 0) =
      (**mActivePowerCommandPu * vq - **mReactivePowerCommandPu * vd) /
      denominator;

  Bool currentReferenceLimited = false;
  **mCurrentReferenceDqPu = limitDqMagnitude(
      currentReference, mMaximumCurrentReferencePu, currentReferenceLimited);
  **mCurrentErrorDqPu = **mCurrentReferenceDqPu - **mFilteredPccCurrentDqPu;
  mPreviousCurrentErrorPu = **mCurrentErrorDqPu;

  const Matrix sourceVoltagePu =
      sourceVoltageComplex.real() / mBaseVoltagePhasePeak;
  const Matrix sourceVoltageDqPu = abcToDq(sourceVoltagePu, **mTheta);

  const Real omegaPu = **mFrequencyPu;
  const Real idRef = (**mCurrentReferenceDqPu)(0, 0);
  const Real iqRef = (**mCurrentReferenceDqPu)(1, 0);
  const Real vdFeedforward = (**mFilteredPccVoltageDqPu)(0, 0);
  const Real vqFeedforward = (**mFilteredPccVoltageDqPu)(1, 0);

  Matrix currentControllerIntegral = Matrix::Zero(2, 1);
  currentControllerIntegral(0, 0) =
      sourceVoltageDqPu(0, 0) -
      (**mCurrentControllerKp * (**mCurrentErrorDqPu)(0, 0) +
       **mCurrentControllerFeedforward * vdFeedforward -
       omegaPu * mFilterInductiveReactancePu * iqRef);
  currentControllerIntegral(1, 0) =
      sourceVoltageDqPu(1, 0) -
      (**mCurrentControllerKp * (**mCurrentErrorDqPu)(1, 0) +
       **mCurrentControllerFeedforward * vqFeedforward +
       omegaPu * mFilterInductiveReactancePu * idRef);

  **mCurrentControllerIntegralPu = currentControllerIntegral;
  **mVoltageCommandDqPu = sourceVoltageDqPu;

  const Real sourceVoltageMagnitudePu =
      std::hypot(sourceVoltageDqPu(0, 0), sourceVoltageDqPu(1, 0));
  const Bool voltageCommandLimitedAtInitialization =
      mMaximumVoltageCommandPu > 0.0 &&
      sourceVoltageMagnitudePu > mMaximumVoltageCommandPu;

  // Initialization residuals. These should all be close to machine precision
  // for a genuinely bumpless PF -> EMT transfer.
  const Matrix capacitorCurrentDqPu = abcToDq(**mCapacitorCurrentPu, **mTheta);
  const Matrix filterKclResidualPu =
      **mFilterCurrentDqPu - **mPccCurrentDqPu - capacitorCurrentDqPu;

  const Real initialMeasuredActivePowerPu =
      vd * (**mPccCurrentDqPu)(0, 0) + vq * (**mPccCurrentDqPu)(1, 0);
  const Real initialMeasuredReactivePowerPu =
      vq * (**mPccCurrentDqPu)(0, 0) - vd * (**mPccCurrentDqPu)(1, 0);

  const Real activePowerResidualPu =
      initialMeasuredActivePowerPu - **mElecActivePowerPu;
  const Real reactivePowerResidualPu =
      initialMeasuredReactivePowerPu - **mElecReactivePowerPu;

  const Real pllOmegaResidual =
      mBaseOmega + **mPllKp * **mPllVoltageErrorPu + **mPllIntegral - **mOmega;
  const Real pllDerivative = **mPllKi * **mPllVoltageErrorPu;
  const Matrix currentIntegralDerivative =
      **mCurrentControllerKi * **mCurrentErrorDqPu;

  updatePhysicalMirrors();
  writeVoltageReference();

  const Real pllFrequencyResidual =
      **mOmega -
      (mBaseOmega + **mPllKp * **mPllVoltageErrorPu + **mPllIntegral);

  const Matrix voltageCommandResidual =
      **mVoltageCommandDqPu - sourceVoltageDqPu;

  const Real currentReferenceMagnitude = std::hypot(
      (**mCurrentReferenceDqPu)(0, 0), (**mCurrentReferenceDqPu)(1, 0));

  const Real voltageCommandMagnitude =
      std::hypot((**mVoltageCommandDqPu)(0, 0), (**mVoltageCommandDqPu)(1, 0));

  SPDLOG_INFO("\n========== GFL CONTROLLER INITIALIZATION =========="
              "\nPLL:"
              "\n  theta = {} rad"
              "\n  vq error = {} pu"
              "\n  PLL integral = {}"
              "\n  frequency = {} pu"
              "\n  omega equilibrium residual = {} rad/s"
              "\nPCC:"
              "\n  Vdq = [{}, {}] pu"
              "\n  |V| = {} pu"
              "\n  Ipcc_dq = [{}, {}] pu"
              "\nFilter:"
              "\n  Ilf_dq = [{}, {}] pu"
              "\nReferences and errors:"
              "\n  Pcmd = {} pu"
              "\n  Qcmd = {} pu"
              "\n  Iref_dq = [{}, {}] pu"
              "\n  current error = [{}, {}] pu, norm = {}"
              "\nCurrent controller:"
              "\n  integral = [{}, {}] pu"
              "\n  Vsource expected = [{}, {}] pu"
              "\n  Vcommand = [{}, {}] pu"
              "\n  Vcommand residual = [{}, {}] pu, norm = {}"
              "\nLimits:"
              "\n  |Iref| = {}, limit = {}"
              "\n  |Vcmd| = {}, limit = {}"
              "\n===================================================",

              **mTheta, **mPllVoltageErrorPu, **mPllIntegral, **mFrequencyPu,
              pllFrequencyResidual,

              (**mPccVoltageDqPu)(0, 0), (**mPccVoltageDqPu)(1, 0),
              **mVoltageMagnitudePu,

              (**mPccCurrentDqPu)(0, 0), (**mPccCurrentDqPu)(1, 0),

              (**mFilterCurrentDqPu)(0, 0), (**mFilterCurrentDqPu)(1, 0),

              **mActivePowerCommandPu, **mReactivePowerCommandPu,

              (**mCurrentReferenceDqPu)(0, 0), (**mCurrentReferenceDqPu)(1, 0),

              (**mCurrentErrorDqPu)(0, 0), (**mCurrentErrorDqPu)(1, 0),
              (**mCurrentErrorDqPu).norm(),

              (**mCurrentControllerIntegralPu)(0, 0),
              (**mCurrentControllerIntegralPu)(1, 0),

              sourceVoltageDqPu(0, 0), sourceVoltageDqPu(1, 0),

              (**mVoltageCommandDqPu)(0, 0), (**mVoltageCommandDqPu)(1, 0),

              voltageCommandResidual(0, 0), voltageCommandResidual(1, 0),
              voltageCommandResidual.norm(),

              currentReferenceMagnitude, mMaximumCurrentReferencePu,

              voltageCommandMagnitude, mMaximumVoltageCommandPu);

  mControlStateInitialized = true;

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFL_Siemens PF initialization:"
      "\n  terminal/PCC: P={} pu, Q={} pu, V={} pu, theta={} rad"
      "\n  PLL: vq_error={} pu, xi_pll={} rad/s, f={} pu, "
      "omega_residual={} rad/s, dxi_pll/dt={}"
      "\n  commands: Pcmd={} pu, Qcmd={} pu"
      "\n  v_pcc_dq=[{}, {}] pu"
      "\n  i_pcc_dq=[{}, {}] pu"
      "\n  i_cf_dq=[{}, {}] pu"
      "\n  i_lf_dq=[{}, {}] pu"
      "\n  i_ref_dq=[{}, {}] pu, e_i=[{}, {}] pu"
      "\n  xi_i=[{}, {}] pu, dxi_i/dt=[{}, {}] pu/s"
      "\n  v_source_dq=v_cmd_dq=[{}, {}] pu, |v_cmd|={} pu"
      "\n  residuals: KCL_dq=[{}, {}] pu, dP={} pu, dQ={} pu"
      "\n  limits active at initialization: current={}, voltage={}",
      **mElecActivePowerPu, **mElecReactivePowerPu, **mVoltageMagnitudePu,
      **mTheta, **mPllVoltageErrorPu, **mPllIntegral, **mFrequencyPu,
      pllOmegaResidual, pllDerivative, **mActivePowerCommandPu,
      **mReactivePowerCommandPu, vd, vq, (**mPccCurrentDqPu)(0, 0),
      (**mPccCurrentDqPu)(1, 0), capacitorCurrentDqPu(0, 0),
      capacitorCurrentDqPu(1, 0), (**mFilterCurrentDqPu)(0, 0),
      (**mFilterCurrentDqPu)(1, 0), (**mCurrentReferenceDqPu)(0, 0),
      (**mCurrentReferenceDqPu)(1, 0), (**mCurrentErrorDqPu)(0, 0),
      (**mCurrentErrorDqPu)(1, 0), (**mCurrentControllerIntegralPu)(0, 0),
      (**mCurrentControllerIntegralPu)(1, 0), currentIntegralDerivative(0, 0),
      currentIntegralDerivative(1, 0), sourceVoltageDqPu(0, 0),
      sourceVoltageDqPu(1, 0), sourceVoltageMagnitudePu,
      filterKclResidualPu(0, 0), filterKclResidualPu(1, 0),
      activePowerResidualPu, reactivePowerResidualPu, currentReferenceLimited,
      voltageCommandLimitedAtInitialization);

  if (currentReferenceLimited || voltageCommandLimitedAtInitialization) {
    SPDLOG_LOGGER_WARN(
        mSLog,
        "GFL_Siemens initial operating point is controller-limited. "
        "A bumpless initialization is impossible with the configured limits.");
  }

  const Real residualTolerance = 1.0e-8;
  if (filterKclResidualPu.norm() > residualTolerance ||
      std::abs(activePowerResidualPu) > residualTolerance ||
      std::abs(reactivePowerResidualPu) > residualTolerance ||
      std::abs(pllOmegaResidual) > residualTolerance) {
    SPDLOG_LOGGER_WARN(
        mSLog,
        "GFL_Siemens initialization residual exceeds {}. Inspect the "
        "printed PF/EMT state before interpreting subsequent oscillations.",
        residualTolerance);
  }
}

void EMT::Ph3::GFL_Siemens::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  requireBaseParameters();
  requireFilterParameters();

  if (!(timeStep > 0.0) || !std::isfinite(timeStep)) {
    throw std::invalid_argument(
        "GFL_Siemens requires a positive finite time step");
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
      "GFL_Siemens measurement filters: dt={} s, tau_P={} s, tau_Q={} s, "
      "tau_I={} s, tau_Vff={} s, alpha_P={}, alpha_Q={}, alpha_I={}, "
      "alpha_Vff={}",
      mTimeStep, mActivePowerMeasurementTimeConstant,
      mReactivePowerMeasurementTimeConstant, mCurrentMeasurementTimeConstant,
      mVoltageFeedforwardMeasurementTimeConstant,
      firstOrderAlpha(mTimeStep, mActivePowerMeasurementTimeConstant),
      firstOrderAlpha(mTimeStep, mReactivePowerMeasurementTimeConstant),
      firstOrderAlpha(mTimeStep, mCurrentMeasurementTimeConstant),
      firstOrderAlpha(mTimeStep, mVoltageFeedforwardMeasurementTimeConstant));

  (void)leftVector;
}

void EMT::Ph3::GFL_Siemens::updateMeasurements() {
  **mPccCurrent = -**mIntfCurrent;

  // Branch-current orientation is terminal 1 -> terminal 0. Invert both
  // internal branch currents for a generator-positive controller convention.
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

  // Same algebra as the CSI "Power calculation" subsystem.
  **mElecActivePowerPu = vd * id + vq * iq;
  **mElecReactivePowerPu = vq * id - vd * iq;
  **mVoltageMagnitudePu = std::hypot(vd, vq);
}

void EMT::Ph3::GFL_Siemens::updateMeasurementFilters() {
  const Real activePowerAlpha =
      firstOrderAlpha(mTimeStep, mActivePowerMeasurementTimeConstant);
  const Real reactivePowerAlpha =
      firstOrderAlpha(mTimeStep, mReactivePowerMeasurementTimeConstant);
  const Real currentAlpha =
      firstOrderAlpha(mTimeStep, mCurrentMeasurementTimeConstant);
  const Real voltageFeedforwardAlpha =
      firstOrderAlpha(mTimeStep, mVoltageFeedforwardMeasurementTimeConstant);

  **mFilteredActivePowerPu +=
      activePowerAlpha * (**mElecActivePowerPu - **mFilteredActivePowerPu);
  **mFilteredReactivePowerPu +=
      reactivePowerAlpha *
      (**mElecReactivePowerPu - **mFilteredReactivePowerPu);
  **mFilteredPccCurrentDqPu +=
      currentAlpha * (**mPccCurrentDqPu - **mFilteredPccCurrentDqPu);
  **mFilteredPccVoltageDqPu +=
      voltageFeedforwardAlpha * (**mPccVoltageDqPu - **mFilteredPccVoltageDqPu);
}

void EMT::Ph3::GFL_Siemens::updateController(Int timeStepCount) {
  updateMeasurements();
  updateMeasurementFilters();

  // ---------------------------------------------------------------------
  // SRF PLL:
  //
  //   xi_dot = Ki * vq
  //   omega  = omega_base + Kp * vq + xi
  //   theta_dot = omega
  //
  // Use the same trapezoidal PLL-state integration as updateOpenLoop().
  // The previous implementation omitted the xi_pll update in closed-loop
  // operation, effectively reducing the PLL to its proportional path.
  // ---------------------------------------------------------------------
  **mPllVoltageErrorPu = (**mPccVoltageDqPu)(1, 0);

  **mPllIntegral += 0.5 * mTimeStep * **mPllKi *
                    (mPreviousPllVoltageErrorPu + **mPllVoltageErrorPu);

  const Real omegaUnclamped =
      mBaseOmega + **mPllKp * **mPllVoltageErrorPu + **mPllIntegral;
  const Real frequencyUnclampedPu = omegaUnclamped / mBaseOmega;

  **mFrequencyPu = clampFinite(frequencyUnclampedPu, mMinimumFrequencyPu,
                               mMaximumFrequencyPu, **mFrequencyRefPu);
  **mOmega = **mFrequencyPu * mBaseOmega;

  // PLL anti-windup: when the frequency is clamped, back-calculate the
  // integral state so that omega = omega_base + Kp*vq + xi remains exact.
  if (**mFrequencyPu != frequencyUnclampedPu) {
    **mPllIntegral = **mOmega - mBaseOmega - **mPllKp * **mPllVoltageErrorPu;
  }

  // mnaParentPreStep(step=0) is called for the first solution at t=dt.
  // The initialized theta represents t=0, so advance it on every pre-step.
  // Use the previous and newly calculated PLL frequencies consistently,
  // matching updateOpenLoop().
  **mTheta +=
      0.5 * mTimeStep * mBaseOmega * (mPreviousFrequencyPu + **mFrequencyPu);

  **mTheta = std::remainder(**mTheta, TWO_PI);

  // ---------------------------------------------------------------------
  // CSI inverse droop laws. No PT command filters are retained.
  // ---------------------------------------------------------------------
  **mActivePowerCommandPu =
      **mActivePowerRefPu -
      **mFrequencyToActivePowerGainPu * (**mFrequencyPu - **mFrequencyRefPu);
  **mReactivePowerCommandPu =
      **mReactivePowerRefPu - **mVoltageToReactivePowerGainPu *
                                  (**mVoltageMagnitudePu - **mVoltageRefPu);

  // ---------------------------------------------------------------------
  // Complex power division used in the CSI block:
  //   i_ref = conj((P_cmd + j Q_cmd)/(v_d + j v_q)).
  // ---------------------------------------------------------------------
  const Real vd = (**mPccVoltageDqPu)(0, 0);
  const Real vq = (**mPccVoltageDqPu)(1, 0);
  const Real denominator =
      std::max(vd * vd + vq * vq, mMinimumVoltageForCurrentReferencePu *
                                      mMinimumVoltageForCurrentReferencePu);

  Matrix currentReferenceUnclamped = Matrix::Zero(2, 1);
  currentReferenceUnclamped(0, 0) =
      (**mActivePowerCommandPu * vd + **mReactivePowerCommandPu * vq) /
      denominator;
  currentReferenceUnclamped(1, 0) =
      (**mActivePowerCommandPu * vq - **mReactivePowerCommandPu * vd) /
      denominator;

  Bool currentReferenceLimited = false;
  **mCurrentReferenceDqPu =
      limitDqMagnitude(currentReferenceUnclamped, mMaximumCurrentReferencePu,
                       currentReferenceLimited);

  // ---------------------------------------------------------------------
  // CSI SRF current controller. The Simulink subsystem uses PCC current as
  // the controlled current and current reference in the decoupling paths.
  // ---------------------------------------------------------------------
  **mCurrentErrorDqPu = **mCurrentReferenceDqPu - **mFilteredPccCurrentDqPu;

  **mCurrentControllerIntegralPu +=
      0.5 * mTimeStep * **mCurrentControllerKi *
      (mPreviousCurrentErrorPu + **mCurrentErrorDqPu);

  const Real omegaPu = **mFrequencyPu;
  const Real idRef = (**mCurrentReferenceDqPu)(0, 0);
  const Real iqRef = (**mCurrentReferenceDqPu)(1, 0);

  // Feed forward a sensor-filtered PCC voltage instead of the raw EMT
  // instantaneous dq voltage. This prevents fast switching/network transients
  // from being copied directly into the controlled-source voltage command.
  const Real vdFeedforward = (**mFilteredPccVoltageDqPu)(0, 0);
  const Real vqFeedforward = (**mFilteredPccVoltageDqPu)(1, 0);

  Matrix voltageCommandUnclamped = Matrix::Zero(2, 1);
  voltageCommandUnclamped(0, 0) =
      **mCurrentControllerKp * (**mCurrentErrorDqPu)(0, 0) +
      (**mCurrentControllerIntegralPu)(0, 0) +
      **mCurrentControllerFeedforward * vdFeedforward -
      omegaPu * mFilterInductiveReactancePu * iqRef;
  voltageCommandUnclamped(1, 0) =
      **mCurrentControllerKp * (**mCurrentErrorDqPu)(1, 0) +
      (**mCurrentControllerIntegralPu)(1, 0) +
      **mCurrentControllerFeedforward * vqFeedforward +
      omegaPu * mFilterInductiveReactancePu * idRef;

  Bool voltageCommandLimited = false;
  **mVoltageCommandDqPu = limitDqMagnitude(
      voltageCommandUnclamped, mMaximumVoltageCommandPu, voltageCommandLimited);

  if (voltageCommandLimited) {
    **mCurrentControllerIntegralPu +=
        **mVoltageCommandDqPu - voltageCommandUnclamped;
  }

  // No 1 us voltage-command delay: apply the current-controller output
  // directly to the inverse Park transform.
  // At timeStepCount==0 the initialized angle already represents t=0.
  // Advancing it here would command the source one full EMT step ahead of
  // the external network and create an artificial startup impulse.

  mPreviousPllVoltageErrorPu = **mPllVoltageErrorPu;
  mPreviousFrequencyPu = **mFrequencyPu;
  mPreviousCurrentErrorPu = **mCurrentErrorDqPu;

  if (!(**mCurrentReferenceDqPu).allFinite() ||
      !(**mVoltageCommandDqPu).allFinite() || !std::isfinite(**mTheta) ||
      !std::isfinite(**mFrequencyPu)) {
    throw std::runtime_error(
        "GFL_Siemens controller produced a non-finite state");
  }

  updatePhysicalMirrors();
  writeVoltageReference();

  (void)currentReferenceLimited;
}

void EMT::Ph3::GFL_Siemens::updateOpenLoop(Int timeStepCount) {
  updateMeasurements();
  updateMeasurementFilters();

  **mPllVoltageErrorPu = (**mPccVoltageDqPu)(1, 0);
  **mPllIntegral += 0.5 * mTimeStep * **mPllKi *
                    (mPreviousPllVoltageErrorPu + **mPllVoltageErrorPu);

  const Real frequencyPu =
      (mBaseOmega + **mPllKp * **mPllVoltageErrorPu + **mPllIntegral) /
      mBaseOmega;
  **mFrequencyPu = clampFinite(frequencyPu, mMinimumFrequencyPu,
                               mMaximumFrequencyPu, **mFrequencyRefPu);
  **mOmega = **mFrequencyPu * mBaseOmega;

  // At timeStepCount==0 the initialized angle already represents t=0.
  // Advancing it here would command the source one full EMT step ahead of
  // the external network and create an artificial startup impulse.
  // mnaParentPreStep(step=0) is called for the first solution at t=dt.
  // The initialized theta represents t=0, so it must be advanced on every
  // pre-step, including timeStepCount == 0.
  **mTheta +=
      0.5 * mTimeStep * mBaseOmega * (mPreviousFrequencyPu + **mFrequencyPu);

  **mTheta = std::remainder(**mTheta, TWO_PI);

  mPreviousPllVoltageErrorPu = **mPllVoltageErrorPu;
  mPreviousFrequencyPu = **mFrequencyPu;

  updatePhysicalMirrors();
  writeVoltageReference();
}

void EMT::Ph3::GFL_Siemens::updatePhysicalMirrors() {
  **mElecActivePower = **mElecActivePowerPu * mBaseApparentPower;
  **mElecReactivePower = **mElecReactivePowerPu * mBaseApparentPower;
  **mFilteredActivePower = **mFilteredActivePowerPu * mBaseApparentPower;
  **mFilteredReactivePower = **mFilteredReactivePowerPu * mBaseApparentPower;
  **mFrequency = **mFrequencyPu * mBaseFrequency;
  **mOmega = **mFrequencyPu * mBaseOmega;
  **mVoltageMagnitude = **mVoltageMagnitudePu * mBaseVoltagePhasePeak;
}

void EMT::Ph3::GFL_Siemens::writeVoltageReference() {
  **mVsrefPu = dqToAbc(**mVoltageCommandDqPu, **mTheta);
  **mVsref = **mVsrefPu * mBaseVoltagePhasePeak;
}

void EMT::Ph3::GFL_Siemens::mnaParentAddPreStepDependencies(
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
  attributeDependencies.push_back(mFrequencyToActivePowerGainPu);
  attributeDependencies.push_back(mVoltageToReactivePowerGainPu);
  attributeDependencies.push_back(mPllKp);
  attributeDependencies.push_back(mPllKi);
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
  modifiedAttributes.push_back(mFilteredPccVoltageDqPu);
  modifiedAttributes.push_back(mPccCurrentDqPu);
  modifiedAttributes.push_back(mFilteredPccCurrentDqPu);
  modifiedAttributes.push_back(mFilterCurrentDqPu);
  modifiedAttributes.push_back(mCurrentReferenceDqPu);
  modifiedAttributes.push_back(mCurrentErrorDqPu);
  modifiedAttributes.push_back(mCurrentControllerIntegralPu);
  modifiedAttributes.push_back(mVoltageCommandDqPu);
  modifiedAttributes.push_back(mElecActivePowerPu);
  modifiedAttributes.push_back(mElecReactivePowerPu);
  modifiedAttributes.push_back(mFilteredActivePowerPu);
  modifiedAttributes.push_back(mFilteredReactivePowerPu);
  modifiedAttributes.push_back(mActivePowerCommandPu);
  modifiedAttributes.push_back(mReactivePowerCommandPu);
  modifiedAttributes.push_back(mVoltageMagnitudePu);
  modifiedAttributes.push_back(mPllVoltageErrorPu);
  modifiedAttributes.push_back(mPllIntegral);
  modifiedAttributes.push_back(mFrequencyPu);
  modifiedAttributes.push_back(mOmega);
  modifiedAttributes.push_back(mTheta);
  modifiedAttributes.push_back(mVsrefPu);
  modifiedAttributes.push_back(mVsref);
  modifiedAttributes.push_back(mElecActivePower);
  modifiedAttributes.push_back(mElecReactivePower);
  modifiedAttributes.push_back(mFilteredActivePower);
  modifiedAttributes.push_back(mFilteredReactivePower);
  modifiedAttributes.push_back(mFrequency);
  modifiedAttributes.push_back(mVoltageMagnitude);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::GFL_Siemens::mnaParentPreStep(Real time, Int timeStepCount) {
  if (!mControlStateInitialized) {
    throw std::runtime_error(
        "GFL_Siemens controller was not initialized from power-flow data");
  }

  if (mWithControl)
    updateController(timeStepCount);
  else
    updateOpenLoop(timeStepCount);

  if (timeStepCount <= 2) {
    SPDLOG_LOGGER_INFO(
        mSLog,
        "GFL_Siemens pre-step {} at t={} s:"
        "\n  Vdq=[{}, {}] pu, Vff_dq=[{}, {}] pu"
        "\n  Idq=[{}, {}] pu, Ifilt_dq=[{}, {}] pu"
        "\n  Iref=[{}, {}] pu, Ierr=[{}, {}] pu"
        "\n  P/Q measured=[{}, {}] pu, filtered=[{}, {}] pu"
        "\n  PLL: vq={}, xi={}, f={} pu, theta={}"
        "\n  P/Q command=[{}, {}] pu"
        "\n  xi_i=[{}, {}] pu, Vcmd=[{}, {}] pu",
        timeStepCount, time, (**mPccVoltageDqPu)(0, 0),
        (**mPccVoltageDqPu)(1, 0), (**mFilteredPccVoltageDqPu)(0, 0),
        (**mFilteredPccVoltageDqPu)(1, 0), (**mPccCurrentDqPu)(0, 0),
        (**mPccCurrentDqPu)(1, 0), (**mFilterCurrentDqPu)(0, 0),
        (**mFilterCurrentDqPu)(1, 0), (**mCurrentReferenceDqPu)(0, 0),
        (**mCurrentReferenceDqPu)(1, 0), (**mCurrentErrorDqPu)(0, 0),
        (**mCurrentErrorDqPu)(1, 0), **mElecActivePowerPu,
        **mElecReactivePowerPu, **mFilteredActivePowerPu,
        **mFilteredReactivePowerPu, **mPllVoltageErrorPu, **mPllIntegral,
        **mFrequencyPu, **mTheta, **mActivePowerCommandPu,
        **mReactivePowerCommandPu, (**mCurrentControllerIntegralPu)(0, 0),
        (**mCurrentControllerIntegralPu)(1, 0), (**mVoltageCommandDqPu)(0, 0),
        (**mVoltageCommandDqPu)(1, 0));
  }

  mSubControlledVoltageSource->mVoltageRef->set(PEAK1PH_TO_RMS3PH * **mVsref);

  std::dynamic_pointer_cast<MNAInterface>(mSubControlledVoltageSource)
      ->mnaPreStep(time, timeStepCount);

  mnaApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::GFL_Siemens::mnaParentAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
  (void)prevStepDependencies;
}

void EMT::Ph3::GFL_Siemens::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateCurrent(**leftVector);
  mnaCompUpdateVoltage(**leftVector);

  if (timeStepCount <= 2) {
    const Matrix pccCurrent = -**mIntfCurrent;
    const Matrix filterCurrent = -mSubFilterInductor->mIntfCurrent->get();
    const Matrix capacitorCurrent = -mSubFilterCapacitor->mIntfCurrent->get();
    const Matrix kclResidual = filterCurrent - pccCurrent - capacitorCurrent;

    SPDLOG_LOGGER_INFO(
        mSLog,
        "GFL_Siemens post-step {} at t={} s:"
        "\n  v_pcc_abc=[{}, {}, {}] V"
        "\n  i_pcc_abc=[{}, {}, {}] A"
        "\n  i_lf_abc=[{}, {}, {}] A"
        "\n  i_cf_abc=[{}, {}, {}] A"
        "\n  KCL residual=[{}, {}, {}] A, norm={}",
        timeStepCount, time, (**mIntfVoltage)(0, 0), (**mIntfVoltage)(1, 0),
        (**mIntfVoltage)(2, 0), pccCurrent(0, 0), pccCurrent(1, 0),
        pccCurrent(2, 0), filterCurrent(0, 0), filterCurrent(1, 0),
        filterCurrent(2, 0), capacitorCurrent(0, 0), capacitorCurrent(1, 0),
        capacitorCurrent(2, 0), kclResidual(0, 0), kclResidual(1, 0),
        kclResidual(2, 0), kclResidual.norm());
  }
}

void EMT::Ph3::GFL_Siemens::mnaCompUpdateCurrent(const Matrix &leftVector) {
  // Composite interface current is consumer-positive:
  // i_intf = i_L,intf - i_C,intf = -i_PCC(generator-positive).
  **mIntfCurrent = mSubFilterInductor->mIntfCurrent->get() -
                   mSubFilterCapacitor->mIntfCurrent->get();
  (void)leftVector;
}

void EMT::Ph3::GFL_Siemens::mnaCompUpdateVoltage(const Matrix &leftVector) {
  for (const auto &virtualNode : mVirtualNodes)
    virtualNode->mnaUpdateVoltage(leftVector);

  (**mIntfVoltage)(0, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
  (**mIntfVoltage)(1, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
  (**mIntfVoltage)(2, 0) =
      Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
}
