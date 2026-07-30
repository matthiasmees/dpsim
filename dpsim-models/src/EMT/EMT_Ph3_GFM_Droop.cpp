// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
// PCC-filter revision 2026-07-30:
// The former main-path coupling resistor Rc has been removed.
// A passive damping resistor Rd is connected in series with the shunt
// capacitor Cf at the filter/PCC node.
// P/Q use the generator-positive PCC current i_pcc = -i_intf.

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

Real clampFinite(Real value, Real lower, Real upper, Real fallback) {
  if (!std::isfinite(value))
    return fallback;
  return std::clamp(value, lower, upper);
}
} // namespace

EMT::Ph3::GFM_Droop::GFM_Droop(String uid, String name, Logger::Level logLevel,
                               Bool withTrafo)
    : CompositePowerComp<Real>(uid, name, true, true, logLevel),
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
      mVs(mAttributes->createDynamic<Matrix>("Vs")) {

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

  // Source pre-step is executed explicitly after the controller update.
  addMNASubComponent(mSubCtrledVoltageSource, MNA_SUBCOMP_TASK_ORDER::NO_TASK,
                     MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

  mVs->setReference(mSubCtrledVoltageSource->mIntfVoltage);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", type(), name);
  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Droop uses scalar P-f/Q-V droop, filtered P/Q measurements, "
      "and an Rf-Lf output filter with a passive series Rd-Cf shunt branch");
}

void EMT::Ph3::GFM_Droop::setParameters(Real frequencyReferenceHz,
                                        Real voltageReferencePeak,
                                        Real activePowerReference,
                                        Real reactivePowerReference) {
  if (!(frequencyReferenceHz > 0.0) || !(voltageReferencePeak >= 0.0) ||
      !std::isfinite(frequencyReferenceHz) ||
      !std::isfinite(voltageReferencePeak) ||
      !std::isfinite(activePowerReference) ||
      !std::isfinite(reactivePowerReference)) {
    throw std::invalid_argument(
        "GFM_Droop references must be finite; f_ref must be positive and "
        "V_ref must be non-negative");
  }

  mParametersSet = true;
  **mFrequencyRef = frequencyReferenceHz;
  **mVoltageRef = voltageReferencePeak;
  **mActivePowerRef = activePowerReference;
  **mReactivePowerRef = reactivePowerReference;

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Droop references: f_ref={} Hz, V_ref={} V_peak, P_ref={} W, "
      "Q_ref={} var",
      frequencyReferenceHz, voltageReferencePeak, activePowerReference,
      reactivePowerReference);
}

void EMT::Ph3::GFM_Droop::setDroopParameters(Real activePowerDroop,
                                             Real reactivePowerDroop,
                                             Real voltageIntegralGain) {
  if (!std::isfinite(activePowerDroop) || !std::isfinite(reactivePowerDroop) ||
      !(voltageIntegralGain >= 0.0) || !std::isfinite(voltageIntegralGain)) {
    throw std::invalid_argument(
        "GFM_Droop gains must be finite and k_iv must be non-negative");
  }

  **mActivePowerDroop = activePowerDroop;
  **mReactivePowerDroop = reactivePowerDroop;
  **mVoltageIntegralGain = voltageIntegralGain;

  SPDLOG_LOGGER_INFO(mSLog,
                     "GFM_Droop gains: k_p={} Hz/W, k_q={} V/var, "
                     "k_iv={} 1/s",
                     activePowerDroop, reactivePowerDroop, voltageIntegralGain);
}

void EMT::Ph3::GFM_Droop::setPowerFilterTimeConstant(Real timeConstant) {
  if (!(timeConstant >= 0.0) || !std::isfinite(timeConstant)) {
    throw std::invalid_argument(
        "GFM_Droop P/Q filter time constant must be finite and non-negative");
  }

  mPowerFilterTimeConstant = timeConstant;

  if (timeConstant > 0.0) {
    SPDLOG_LOGGER_INFO(
        mSLog, "GFM_Droop P/Q measurement filter: tau={} s, cutoff={} Hz",
        timeConstant, 1.0 / (TWO_PI * timeConstant));
  } else {
    SPDLOG_LOGGER_INFO(mSLog, "GFM_Droop P/Q measurement filter bypassed");
  }
}

void EMT::Ph3::GFM_Droop::setControllerLimits(Real minimumFrequencyHz,
                                              Real maximumFrequencyHz,
                                              Real minimumVoltagePeak,
                                              Real maximumVoltagePeak) {
  if (!(minimumFrequencyHz >= 0.0) ||
      !(maximumFrequencyHz > minimumFrequencyHz) ||
      !(minimumVoltagePeak >= 0.0) ||
      !(maximumVoltagePeak > minimumVoltagePeak) ||
      !std::isfinite(minimumFrequencyHz) ||
      !std::isfinite(maximumFrequencyHz) ||
      !std::isfinite(minimumVoltagePeak) ||
      !std::isfinite(maximumVoltagePeak)) {
    throw std::invalid_argument("Invalid GFM_Droop controller limits");
  }

  mMinimumFrequency = minimumFrequencyHz;
  mMaximumFrequency = maximumFrequencyHz;
  mMinimumVoltage = minimumVoltagePeak;
  mMaximumVoltage = maximumVoltagePeak;
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
}

void EMT::Ph3::GFM_Droop::setFilterParameters(Real Lf, Real Cf, Real Rf,
                                              Real Rd) {
  if (!(Lf > 0.0) || !(Cf > 0.0) || !(Rf >= 0.0) || !(Rd > 0.0) ||
      !std::isfinite(Lf) || !std::isfinite(Cf) || !std::isfinite(Rf) ||
      !std::isfinite(Rd)) {
    throw std::invalid_argument(
        "GFM_Droop filter requires Lf>0, Cf>0, Rf>=0, and Rd>0");
  }

  // The base class still stores an Rc field for other VSI implementations.
  // This topology has no resistor in the main path after Lf.
  Base::AvVoltageSourceInverterDQ::setFilterParameters(Lf, Cf, Rf, 0.0);
  mCapacitorDampingResistance = Rd;

  mSubResistorF->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mRf));
  mSubResistorD->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mCapacitorDampingResistance));
  mSubInductorF->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mLf));
  mSubCapacitorF->setParameters(
      CPS::Math::singlePhaseParameterToThreePhase(mCf));

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Droop electrical filter: Lf={} H, Cf={} F, Rf={} Ohm, "
      "Rd={} Ohm; Rd is in series with Cf",
      mLf, mCf, mRf, mCapacitorDampingResistance);
}

void EMT::Ph3::GFM_Droop::initializeParentFromNodesAndTerminals(
    Real frequency) {
  MatrixComp intfVoltageComplex = MatrixComp::Zero(3, 1);
  MatrixComp intfCurrentComplex = MatrixComp::Zero(3, 1);

  // terminal(0)->singlePower() is used consistently with the existing VSI
  // initialization. The computed controller powers below are generator-positive.
  const Real activePower = terminal(0)->singlePower().real();
  const Real reactivePower = terminal(0)->singlePower().imag();

  intfVoltageComplex(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  intfVoltageComplex(1, 0) = intfVoltageComplex(0, 0) * SHIFT_TO_PHASE_B;
  intfVoltageComplex(2, 0) = intfVoltageComplex(0, 0) * SHIFT_TO_PHASE_C;

  intfCurrentComplex(0, 0) =
      -std::conj(2.0 / 3.0 * Complex(activePower, reactivePower) /
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

  // Passive damping branch at the filter interface:
  //
  //   filter/PCC node -- Rd -- capacitor node -- Cf -- GND
  //
  // dampingBranchCurrentInit is defined from the filter/PCC node to ground.
  const Complex capacitorImpedance(0.0, -1.0 / (omegaInit * mCf));
  const Complex dampingBranchImpedance(mCapacitorDampingResistance, 0.0);
  const Complex branchImpedance = dampingBranchImpedance + capacitorImpedance;

  MatrixComp dampingBranchCurrentInit = MatrixComp::Zero(3, 1);
  MatrixComp capacitorVoltageInit = MatrixComp::Zero(3, 1);
  for (UInt phase = 0; phase < 3; ++phase) {
    dampingBranchCurrentInit(phase, 0) =
        filterInterfaceInitialVoltage(phase, 0) / branchImpedance;
    capacitorVoltageInit(phase, 0) =
        dampingBranchCurrentInit(phase, 0) * capacitorImpedance;
  }

  // The inductor interface current is positive from the filter/PCC node back
  // towards the controlled source. KCL therefore gives:
  //
  //   i_L = i_intf - i_RdCf,out
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
  // The DPsim interface current is consumer-positive. The droop control
  // uses generator-positive current flowing from the GFM into the grid.
  **mPccCurrent = -intfCurrentComplex.real();

  **mElecActivePower =
      (1.5 * intfVoltageComplex(0, 0) * std::conj(-intfCurrentComplex(0, 0)))
          .real();
  **mElecReactivePower =
      (1.5 * intfVoltageComplex(0, 0) * std::conj(-intfCurrentComplex(0, 0)))
          .imag();

  // Bumpless initialization of the measurement filters.
  **mFilteredActivePower = **mElecActivePower;
  **mFilteredReactivePower = **mElecReactivePower;

  **mVoltageMagnitude = std::abs(intfVoltageComplex(0, 0));

  mSubCtrledVoltageSource->setParameters(mVirtualNodes[0]->initialVoltage(),
                                         0.0);

  mSubCtrledVoltageSource->connect({SimNode::GND, mVirtualNodes[0]});
  mSubResistorF->connect({mVirtualNodes[0], mVirtualNodes[1]});

  if (mWithConnectionTransformer) {
    mSubInductorF->connect({mVirtualNodes[1], mVirtualNodes[2]});
    mSubResistorD->connect({mVirtualNodes[2], mVirtualNodes[3]});
    mSubCapacitorF->connect({mVirtualNodes[3], SimNode::GND});
  } else {
    // The Lf/Rd node is the external PCC. The capacitor itself is behind Rd.
    mSubInductorF->connect({mVirtualNodes[1], mTerminals[0]->node()});
    mSubResistorD->connect({mTerminals[0]->node(), mVirtualNodes[2]});
    mSubCapacitorF->connect({mVirtualNodes[2], SimNode::GND});
  }

  for (const auto &subcomp : mSubComponents) {
    subcomp->initialize(mFrequencies);
    subcomp->initializeFromNodesAndTerminals(frequency);
  }

  // Bumpless initialization: V1 starts at the source-side PF voltage. The
  // droop frequency starts from the filtered PF power.
  **mTheta = std::arg(vsInit(0, 0));
  **mFrequency = clampFinite(
      **mFrequencyRef +
          **mActivePowerDroop * (**mActivePowerRef - **mFilteredActivePower),
      mMinimumFrequency, mMaximumFrequency, **mFrequencyRef);
  **mOmega = TWO_PI * **mFrequency;
  mPreviousOmega = **mOmega;

  **mVoltageDroopOutput =
      **mVoltageRef +
      **mReactivePowerDroop * (**mReactivePowerRef - **mFilteredReactivePower);
  mPreviousVoltageError = **mVoltageDroopOutput - **mVoltageMagnitude;

  **mVoltageCommand = std::abs(vsInit(0, 0));
  **mVoltageIntegralState = **mVoltageCommand - **mVoltageDroopOutput;

  **mVsref = vsInit.real();
  mControlStateInitialized = true;

  SPDLOG_LOGGER_INFO(
      mSLog,
      "GFM_Droop PF initialization: P={} W, Q={} var, P_f={} W, Q_f={} "
      "var, V={} V_peak, f={} Hz, theta={} rad, source amplitude={} V_peak",
      **mElecActivePower, **mElecReactivePower, **mFilteredActivePower,
      **mFilteredReactivePower, **mVoltageMagnitude, **mFrequency, **mTheta,
      **mVoltageCommand);
}

void EMT::Ph3::GFM_Droop::mnaParentInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  if (!(timeStep > 0.0) || !std::isfinite(timeStep))
    throw std::invalid_argument("GFM_Droop requires a positive time step");

  mTimeStep = timeStep;

  // Exact discretization for a first-order filter with held input.
  mPowerFilterAlpha =
      mPowerFilterTimeConstant > 0.0
          ? 1.0 - std::exp(-mTimeStep / mPowerFilterTimeConstant)
          : 1.0;

  SPDLOG_LOGGER_INFO(
      mSLog, "GFM_Droop discrete P/Q filter: dt={} s, tau={} s, alpha={}",
      mTimeStep, mPowerFilterTimeConstant, mPowerFilterAlpha);
}

void EMT::Ph3::GFM_Droop::updateMeasurements() {
  // Measurement point: external GFM terminal (PCC).
  //
  // mIntfVoltage is the PCC phase voltage. mIntfCurrent follows DPsim's
  // component-consumer convention and is positive from the PCC into the GFM.
  // The droop equations use generator-positive current from the GFM into the
  // connected grid, hence i_pcc = -mIntfCurrent.
  const Matrix &vPcc = **mIntfVoltage;
  Matrix &iPcc = **mPccCurrent;
  const Matrix &iInterface = **mIntfCurrent;

  iPcc(0, 0) = -iInterface(0, 0);
  iPcc(1, 0) = -iInterface(1, 0);
  iPcc(2, 0) = -iInterface(2, 0);

  // Power-invariant Clarke quantities. No temporary dq matrices are allocated.
  const Real vAlpha =
      SQRT_TWO_OVER_THREE * (vPcc(0, 0) - 0.5 * vPcc(1, 0) - 0.5 * vPcc(2, 0));
  const Real vBeta =
      SQRT_TWO_OVER_THREE * SQRT_THREE_OVER_TWO * (vPcc(1, 0) - vPcc(2, 0));
  const Real iAlpha =
      SQRT_TWO_OVER_THREE * (iPcc(0, 0) - 0.5 * iPcc(1, 0) - 0.5 * iPcc(2, 0));
  const Real iBeta =
      SQRT_TWO_OVER_THREE * SQRT_THREE_OVER_TWO * (iPcc(1, 0) - iPcc(2, 0));

  // Generator-positive three-phase instantaneous power at the PCC.
  **mElecActivePower = vAlpha * iAlpha + vBeta * iBeta;
  **mElecReactivePower = vBeta * iAlpha - vAlpha * iBeta;

  // For a balanced sinusoidal three-phase voltage this evaluates directly to
  // the phase-to-neutral peak amplitude.
  const Real voltageSquared =
      (2.0 / 3.0) * (vPcc(0, 0) * vPcc(0, 0) + vPcc(1, 0) * vPcc(1, 0) +
                     vPcc(2, 0) * vPcc(2, 0));
  **mVoltageMagnitude = std::sqrt(std::max<Real>(0.0, voltageSquared));
}

void EMT::Ph3::GFM_Droop::updatePowerMeasurementFilter(Int timeStepCount) {
  if (timeStepCount <= 0)
    return;

  **mFilteredActivePower +=
      mPowerFilterAlpha * (**mElecActivePower - **mFilteredActivePower);
  **mFilteredReactivePower +=
      mPowerFilterAlpha * (**mElecReactivePower - **mFilteredReactivePower);
}

void EMT::Ph3::GFM_Droop::updateController(Int timeStepCount) {
  updateMeasurements();
  updatePowerMeasurementFilter(timeStepCount);

  const Real frequencyUnclamped =
      **mFrequencyRef +
      **mActivePowerDroop * (**mActivePowerRef - **mFilteredActivePower);
  **mFrequency = clampFinite(frequencyUnclamped, mMinimumFrequency,
                             mMaximumFrequency, **mFrequencyRef);
  **mOmega = TWO_PI * **mFrequency;

  if (timeStepCount > 0) {
    **mTheta += 0.5 * mTimeStep * (mPreviousOmega + **mOmega);
    **mTheta = std::remainder(**mTheta, TWO_PI);
  }

  **mVoltageDroopOutput =
      **mVoltageRef +
      **mReactivePowerDroop * (**mReactivePowerRef - **mFilteredReactivePower);

  const Real voltageError = **mVoltageDroopOutput - **mVoltageMagnitude;

  if (timeStepCount > 0) {
    **mVoltageIntegralState += 0.5 * mTimeStep * **mVoltageIntegralGain *
                               (mPreviousVoltageError + voltageError);
  }

  const Real voltageUnclamped = **mVoltageDroopOutput + **mVoltageIntegralState;
  **mVoltageCommand = clampFinite(voltageUnclamped, mMinimumVoltage,
                                  mMaximumVoltage, **mVoltageRef);

  // Back-calculation anti-windup without an additional tuning parameter.
  if (**mVoltageCommand != voltageUnclamped)
    **mVoltageIntegralState = **mVoltageCommand - **mVoltageDroopOutput;

  mPreviousOmega = **mOmega;
  mPreviousVoltageError = voltageError;
  writeVoltageReference();
}

void EMT::Ph3::GFM_Droop::updateOpenLoop(Int timeStepCount) {
  updateMeasurements();
  updatePowerMeasurementFilter(timeStepCount);

  **mFrequency = clampFinite(**mFrequencyRef, mMinimumFrequency,
                             mMaximumFrequency, **mFrequencyRef);
  **mOmega = TWO_PI * **mFrequency;

  if (timeStepCount > 0) {
    **mTheta += 0.5 * mTimeStep * (mPreviousOmega + **mOmega);
    **mTheta = std::remainder(**mTheta, TWO_PI);
  }

  **mVoltageDroopOutput = **mVoltageRef;
  **mVoltageIntegralState = 0.0;
  **mVoltageCommand = clampFinite(**mVoltageRef, mMinimumVoltage,
                                  mMaximumVoltage, **mVoltageRef);
  mPreviousOmega = **mOmega;
  writeVoltageReference();
}

void EMT::Ph3::GFM_Droop::writeVoltageReference() {
  const Real theta = **mTheta;
  const Real amplitude = **mVoltageCommand;
  Matrix &reference = **mVsref;

  reference(0, 0) = amplitude * std::cos(theta);
  reference(1, 0) = amplitude * std::cos(theta - TWO_PI_OVER_THREE);
  reference(2, 0) = amplitude * std::cos(theta + TWO_PI_OVER_THREE);
}

void EMT::Ph3::GFM_Droop::mnaParentAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mIntfVoltage);
  prevStepDependencies.push_back(mIntfCurrent);

  attributeDependencies.push_back(mActivePowerRef);
  attributeDependencies.push_back(mReactivePowerRef);
  attributeDependencies.push_back(mFrequencyRef);
  attributeDependencies.push_back(mVoltageRef);
  attributeDependencies.push_back(mActivePowerDroop);
  attributeDependencies.push_back(mReactivePowerDroop);
  attributeDependencies.push_back(mVoltageIntegralGain);

  modifiedAttributes.push_back(mPccCurrent);
  modifiedAttributes.push_back(mElecActivePower);
  modifiedAttributes.push_back(mElecReactivePower);
  modifiedAttributes.push_back(mFilteredActivePower);
  modifiedAttributes.push_back(mFilteredReactivePower);
  modifiedAttributes.push_back(mVoltageMagnitude);
  modifiedAttributes.push_back(mFrequency);
  modifiedAttributes.push_back(mOmega);
  modifiedAttributes.push_back(mTheta);
  modifiedAttributes.push_back(mVoltageDroopOutput);
  modifiedAttributes.push_back(mVoltageIntegralState);
  modifiedAttributes.push_back(mVoltageCommand);
  modifiedAttributes.push_back(mVsref);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::Ph3::GFM_Droop::mnaParentPreStep(Real time, Int timeStepCount) {
  if (!mControlStateInitialized) {
    throw std::runtime_error(
        "GFM_Droop control state was not initialized from power-flow data");
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
}

void EMT::Ph3::GFM_Droop::mnaParentPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateCurrent(**leftVector);
  mnaCompUpdateVoltage(**leftVector);
}

void EMT::Ph3::GFM_Droop::mnaCompUpdateCurrent(const Matrix &leftVector) {
  // Keep mIntfCurrent in DPsim's consumer-positive convention.
  if (mWithConnectionTransformer) {
    **mIntfCurrent = mConnectionTransformer->mIntfCurrent->get();
  } else {
    // Lf is connected as {filter series node, PCC}, so its interface current
    // is positive from the PCC back into the series filter.
    //
    // Rd is connected as {PCC, capacitor node}; its interface current is
    // positive from the capacitor node towards the PCC. Subtracting it adds
    // the physical current flowing from the PCC into the Rd-Cf branch.
    //
    // Hence the total consumer-positive current entering the complete GFM is
    //
    //   i_intf = i_L - i_Rd.
    **mIntfCurrent =
        mSubInductorF->mIntfCurrent->get() - mSubResistorD->mIntfCurrent->get();
  }
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
