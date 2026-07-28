// SPDX-License-Identifier: MPL-2.0
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMCStation.h>

#include <cmath>
#include <stdexcept>

using namespace CPS;
using namespace CPS::EMT::Ph3;

namespace {
void requireFinite(Real value, const char *name) {
  if (!std::isfinite(value))
    throw std::invalid_argument(String(name) + " must be finite.");
}
} // namespace

SSN_MMCStation::SSN_MMCStation(String name, SSN_MMC::Ptr plant,
                               Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mAngle(mAttributes->createDynamic<Real>("angle")),
      mAngularFrequency(mAttributes->createDynamic<Real>("angular_frequency")),
      mActivePowerReferencePu(
          mAttributes->createDynamic<Real>("active_power_reference_pu")),
      mReactivePowerReferencePu(
          mAttributes->createDynamic<Real>("reactive_power_reference_pu")),
      mDcVoltageReference(
          mAttributes->createDynamic<Real>("dc_voltage_reference")),
      mVdPu(mAttributes->create<Real>("vd_pu", 0.0)),
      mVqPu(mAttributes->create<Real>("vq_pu", 0.0)),
      mIdPu(mAttributes->create<Real>("id_pu", 0.0)),
      mIqPu(mAttributes->create<Real>("iq_pu", 0.0)),
      mActivePowerPu(mAttributes->create<Real>("active_power_pu", 0.0)),
      mReactivePowerPu(mAttributes->create<Real>("reactive_power_pu", 0.0)),
      mFilteredActivePowerPu(
          mAttributes->create<Real>("filtered_active_power_pu", 0.0)),
      mFilteredReactivePowerPu(
          mAttributes->create<Real>("filtered_reactive_power_pu", 0.0)),
      mFilteredDcVoltage(mAttributes->create<Real>("filtered_dc_voltage", 0.0)),
      mIdReferencePu(mAttributes->create<Real>("id_reference_pu", 0.0)),
      mIqReferencePu(mAttributes->create<Real>("iq_reference_pu", 0.0)),
      mVdReferencePu(mAttributes->create<Real>("vd_reference_pu", 0.0)),
      mVqReferencePu(mAttributes->create<Real>("vq_reference_pu", 0.0)),
      mConverterPhaseCommand(mAttributes->create<Matrix>(
          "converter_phase_command", Matrix::Zero(3, 1))),
      mPlantDifferentialVoltageCommand(mAttributes->create<Matrix>(
          "plant_differential_voltage_command", Matrix::Zero(2, 1))),
      mModulationMagnitude(
          mAttributes->create<Real>("modulation_magnitude", 0.0)),
      mModulationD(mAttributes->create<Real>("modulation_d", 0.0)),
      mModulationQ(mAttributes->create<Real>("modulation_q", 0.0)),
      mModulationDUnsaturated(
          mAttributes->create<Real>("modulation_d_unsaturated", 0.0)),
      mModulationQUnsaturated(
          mAttributes->create<Real>("modulation_q_unsaturated", 0.0)),
      mAcPower(mAttributes->create<Real>("ac_power", 0.0)),
      mDcPower(mAttributes->create<Real>("dc_power", 0.0)),
      mPowerBalanceError(mAttributes->create<Real>("power_balance_error", 0.0)),
      mControlMode(mAttributes->create<Int>(
          "control_mode",
          static_cast<Int>(ControlMode::DCVoltageReactivePower))),
      mState(
          mAttributes->create<Int>("state", static_cast<Int>(State::Blocked))),
      mEnableDiagnostic(mAttributes->create<Int>(
          "enable_diagnostic", static_cast<Int>(EnableDiagnostic::None))),
      mControllerEnabled(
          mAttributes->create<Bool>("controller_enabled", false)),
      mOuterLoopsEnabled(
          mAttributes->create<Bool>("outer_loops_enabled", true)),
      mCurrentDSaturated(
          mAttributes->create<Bool>("current_d_saturated", false)),
      mCurrentQSaturated(
          mAttributes->create<Bool>("current_q_saturated", false)),
      mModulationSaturated(
          mAttributes->create<Bool>("modulation_saturated", false)),
      mCurrentDUpperSaturated(
          mAttributes->create<Bool>("current_d_upper_saturated", false)),
      mCurrentDLowerSaturated(
          mAttributes->create<Bool>("current_d_lower_saturated", false)),
      mCurrentQUpperSaturated(
          mAttributes->create<Bool>("current_q_upper_saturated", false)),
      mCurrentQLowerSaturated(
          mAttributes->create<Bool>("current_q_lower_saturated", false)),
      mModulationDUpperSaturated(
          mAttributes->create<Bool>("modulation_d_upper_saturated", false)),
      mModulationDLowerSaturated(
          mAttributes->create<Bool>("modulation_d_lower_saturated", false)),
      mModulationQUpperSaturated(
          mAttributes->create<Bool>("modulation_q_upper_saturated", false)),
      mModulationQLowerSaturated(
          mAttributes->create<Bool>("modulation_q_lower_saturated", false)),
      mActiveOuterUpperSaturated(
          mAttributes->create<Bool>("active_outer_upper_saturated", false)),
      mActiveOuterLowerSaturated(
          mAttributes->create<Bool>("active_outer_lower_saturated", false)),
      mDcOuterUpperSaturated(
          mAttributes->create<Bool>("dc_outer_upper_saturated", false)),
      mDcOuterLowerSaturated(
          mAttributes->create<Bool>("dc_outer_lower_saturated", false)),
      mReactiveOuterUpperSaturated(
          mAttributes->create<Bool>("reactive_outer_upper_saturated", false)),
      mReactiveOuterLowerSaturated(
          mAttributes->create<Bool>("reactive_outer_lower_saturated", false)),
      mDcOuterErrorPu(mAttributes->create<Real>("dc_outer_error_pu", 0.0)),
      mDcOuterProportionalContribution(
          mAttributes->create<Real>("dc_outer_proportional", 0.0)),
      mDcOuterIntegralContribution(
          mAttributes->create<Real>("dc_outer_integral", 0.0)),
      mDcOuterUnsaturatedOutput(
          mAttributes->create<Real>("dc_outer_unsaturated_output", 0.0)),
      mDcOuterOutput(mAttributes->create<Real>("dc_outer_output", 0.0)),
      mPlant(std::move(plant)),
      mTransform(Signal::ExternallyAngledDQAdapter::make(name + ".Transform")),
      mPFilter(Signal::DQSymSecondOrderFilter::make(name + ".PFilter")),
      mQFilter(Signal::DQSymSecondOrderFilter::make(name + ".QFilter")),
      mVdcFilter(Signal::DQSymSecondOrderFilter::make(name + ".VdcFilter")),
      mActiveController(Signal::DQSymOuterController::make(name + ".Active")),
      mDcVoltageController(
          Signal::DQSymOuterController::make(name + ".DcVoltage")),
      mReactiveController(
          Signal::DQSymOuterController::make(name + ".Reactive")),
      mCurrentController(
          Signal::DQSymCurrentController::make(name + ".Current")),
      mModulation(Signal::DQSymModulation::make(name + ".Modulation")) {
  if (!mPlant)
    throw std::invalid_argument("SSN_MMCStation requires one SSN_MMC plant.");
  **mAngle = 0.0;
  **mAngularFrequency = 0.0;
  **mActivePowerReferencePu = 0.0;
  **mReactivePowerReferencePu = 0.0;
  **mDcVoltageReference = 0.0;
}

void SSN_MMCStation::setParameters(const SSN_MMCStationParameters &parameters) {
  if (parameters.nominalPower <= 0.0 ||
      parameters.nominalAcLineLineRms <= 0.0 ||
      parameters.nominalDcVoltage <= 0.0 ||
      parameters.nominalFrequencyHz <= 0.0 ||
      parameters.controllerTimeStep <= 0.0 ||
      parameters.measurementFilterFrequencyHz <= 0.0 ||
      parameters.measurementFilterDamping <= 0.0 ||
      !std::isfinite(parameters.dcVoltageIntegralGain) ||
      parameters.dcVoltageIntegralGain < 0.0)
    throw std::invalid_argument("Invalid SSN_MMCStation parameters.");
  mParameters = parameters;

  mPFilter->setParameters(parameters.measurementFilterFrequencyHz,
                          parameters.measurementFilterDamping);
  mQFilter->setParameters(parameters.measurementFilterFrequencyHz,
                          parameters.measurementFilterDamping);
  mVdcFilter->setParameters(parameters.measurementFilterFrequencyHz,
                            parameters.measurementFilterDamping);
  mActiveController->setParameters(
      Signal::DQSymOuterLoopType::ActivePowerToDcVoltage, 0.5 / 3.0, 1.0, 0.8,
      1.2, 1.0, parameters.nominalDcVoltage, 1.0);
  mDcVoltageController->setParameters(
      Signal::DQSymOuterLoopType::DcVoltageToDCurrent, 4.0,
      parameters.dcVoltageIntegralGain, -2.0, 2.0, parameters.nominalDcVoltage,
      1.0, 0.0);
  mReactiveController->setParameters(
      Signal::DQSymOuterLoopType::ReactivePowerToQCurrent, 0.5 / 3.0, 1.0,
      -0.25, 0.25, 1.0, 1.0, 0.0);
  Signal::DQSymCurrentControllerParameters current;
  current.kp = 0.6;
  current.ki = 6.0;
  current.rFeedforwardPu = parameters.armResistancePu / 2.0;
  current.lFeedforwardPu = parameters.armInductancePu / 2.0;
  current.nominalFrequencyHz = parameters.nominalFrequencyHz;
  current.lowerLimitPu = -2.0;
  current.upperLimitPu = 2.0;
  mCurrentController->setParameters(current);
  mModulation->setParameters(parameters.nominalDcVoltage,
                             parameters.nominalAcLineLineRms, -2.0, 2.0);
  mParametersSet = true;
}

void SSN_MMCStation::setControlMode(ControlMode mode) {
  if (**mState == static_cast<Int>(State::Enabled))
    throw std::logic_error("Cannot change MMCStation mode while enabled.");
  **mControlMode = static_cast<Int>(mode);
}

void SSN_MMCStation::setReferences(Real activePowerPu, Real reactivePowerPu,
                                   Real dcVoltageVolts) {
  requireFinite(activePowerPu, "active-power reference");
  requireFinite(reactivePowerPu, "reactive-power reference");
  requireFinite(dcVoltageVolts, "DC-voltage reference");
  **mActivePowerReferencePu = activePowerPu;
  **mReactivePowerReferencePu = reactivePowerPu;
  **mDcVoltageReference = dcVoltageVolts;
}

void SSN_MMCStation::setOuterLoopsEnabled(Bool enabled) {
  if (**mState == static_cast<Int>(State::Enabled))
    throw std::logic_error(
        "Cannot change MMCStation outer-loop selection while enabled.");
  **mOuterLoopsEnabled = enabled;
}

void SSN_MMCStation::setCurrentReferences(Real idReferencePu,
                                          Real iqReferencePu) {
  requireFinite(idReferencePu, "d-current reference");
  requireFinite(iqReferencePu, "q-current reference");
  if (**mOuterLoopsEnabled)
    throw std::logic_error(
        "Direct current references require disabled outer loops.");
  **mIdReferencePu = idReferencePu;
  **mIqReferencePu = iqReferencePu;
}

void SSN_MMCStation::initializeFromOperatingPoint(Real angle,
                                                  Real angularFrequency,
                                                  Real activePowerPu,
                                                  Real reactivePowerPu,
                                                  Real dcVoltageVolts) {
  if (!mParametersSet)
    throw std::logic_error("Set MMCStation parameters before initialization.");
  setReferences(activePowerPu, reactivePowerPu, dcVoltageVolts);
  **mAngle = angle;
  **mAngularFrequency = angularFrequency;
  measurementStep();

  mPFilter->setInitialValue(**mActivePowerPu);
  mQFilter->setInitialValue(**mReactivePowerPu);
  mVdcFilter->setInitialValue(**mPlant->dcVoltageAttribute());

  const Real vBase = std::sqrt(2.0 / 3.0) * mParameters.nominalAcLineLineRms;
  const Real iBase = (2.0 / 3.0) * mParameters.nominalPower / vBase;
  const Real rff = mParameters.armResistancePu / 2.0;
  const Real lff = mParameters.armInductancePu / 2.0;
  const Real omegaPu =
      angularFrequency / (2.0 * PI * mParameters.nominalFrequencyHz);
  const Real dFeedforward = **mVdPu + rff * **mIdPu - omegaPu * lff * **mIqPu;
  const Real qFeedforward = **mVqPu + rff * **mIqPu + omegaPu * lff * **mIdPu;
  (void)iBase;
  mActiveController->setInitialIntegratorState(dcVoltageVolts /
                                               mParameters.nominalDcVoltage);
  const Bool directActivePowerOrientation =
      **mControlMode == static_cast<Int>(ControlMode::ActivePowerReactivePower);
  mDcVoltageController->setInitialIntegratorState(
      directActivePowerOrientation ? **mIdPu : -**mIdPu);
  mReactiveController->setInitialIntegratorState(-**mIqPu);
  mCurrentController->setInitialIntegratorStates(**mVdPu - dFeedforward,
                                                 **mVqPu - qFeedforward);

  **mIdReferencePu = **mIdPu;
  **mIqReferencePu = **mIqPu;
  **mVdReferencePu = **mVdPu;
  **mVqReferencePu = **mVqPu;
  mHeldPlantCommand << **mVdReferencePu * vBase, **mVqReferencePu * vBase;
  **mPlantDifferentialVoltageCommand = mHeldPlantCommand;
  mPlant->setExternalDifferentialVoltageCommand(mHeldPlantCommand(0, 0),
                                                mHeldPlantCommand(1, 0));
  **mState = static_cast<Int>(State::Ready);
  **mControllerEnabled = false;
  setEnableDiagnostic(EnableDiagnostic::None);
}

Bool SSN_MMCStation::requestEnable(Real errorTolerancePu,
                                   Real commandToleranceVolts) {
  mEnableErrorTolerancePu = errorTolerancePu;
  mEnableCommandToleranceVolts = commandToleranceVolts;
  if (**mState != static_cast<Int>(State::Ready)) {
    setEnableDiagnostic(EnableDiagnostic::NotReady);
    return false;
  }
  const Matrix measurements = (Matrix(7, 1) << **mVdPu, **mVqPu, **mIdPu,
                               **mIqPu, **mFilteredActivePowerPu,
                               **mFilteredReactivePowerPu, **mFilteredDcVoltage)
                                  .finished();
  if (!measurements.allFinite()) {
    setEnableDiagnostic(EnableDiagnostic::NonFiniteMeasurement);
    return false;
  }
  if (!std::isfinite(**mActivePowerReferencePu) ||
      !std::isfinite(**mReactivePowerReferencePu) ||
      !std::isfinite(**mDcVoltageReference)) {
    setEnableDiagnostic(EnableDiagnostic::NonFiniteReference);
    return false;
  }
  const Real activeError =
      !**mOuterLoopsEnabled
          ? std::abs(**mIdReferencePu - **mIdPu)
          : (**mControlMode ==
                     static_cast<Int>(ControlMode::ActivePowerReactivePower)
                 ? std::abs(**mActivePowerReferencePu -
                            **mFilteredActivePowerPu)
                 : std::abs(**mDcVoltageReference - **mFilteredDcVoltage) /
                       mParameters.nominalDcVoltage);
  const Real reactiveError =
      !**mOuterLoopsEnabled
          ? std::abs(**mIqReferencePu - **mIqPu)
          : std::abs(**mReactivePowerReferencePu - **mFilteredReactivePowerPu);
  if (activeError > errorTolerancePu || reactiveError > errorTolerancePu) {
    setEnableDiagnostic(EnableDiagnostic::ControlErrorTooLarge);
    return false;
  }
  if ((**mPlantDifferentialVoltageCommand - mHeldPlantCommand).norm() >
      commandToleranceVolts) {
    setEnableDiagnostic(EnableDiagnostic::CommandMismatch);
    return false;
  }
  mPlant->setControlSource(SSN_MMC::ControlSource::ExternalDifferentialVoltage);
  **mState = static_cast<Int>(State::Enabled);
  **mControllerEnabled = true;
  setEnableDiagnostic(EnableDiagnostic::None);
  return true;
}

void SSN_MMCStation::block() {
  **mState = static_cast<Int>(State::Blocked);
  **mControllerEnabled = false;
  mPlant->setControlSource(SSN_MMC::ControlSource::InternalControllers);
}

void SSN_MMCStation::initialize(Real timeStep) {
  if (!mParametersSet)
    throw std::logic_error("MMCStation parameters are not set.");
  if (std::abs(timeStep - mParameters.controllerTimeStep) >
      1e-12 * std::max(1.0, timeStep))
    throw std::invalid_argument(
        "MMCStation simulation and controller time steps must match.");
  if (**mState == static_cast<Int>(State::Blocked)) {
    measurementStep();
    initializeFromOperatingPoint(**mAngle, **mAngularFrequency,
                                 **mActivePowerPu, **mReactivePowerPu,
                                 **mPlant->dcVoltageAttribute());
  }
  mPFilter->initialize(timeStep);
  mQFilter->initialize(timeStep);
  mVdcFilter->initialize(timeStep);
  mActiveController->initialize(timeStep);
  mDcVoltageController->initialize(timeStep);
  mReactiveController->initialize(timeStep);
  mCurrentController->initialize(timeStep);
  mInitialized = true;
  filterAndOuterStep();
  currentStep();
  commandStep();
}

void SSN_MMCStation::measurementStep() {
  if (!mParametersSet)
    throw std::logic_error("MMCStation parameters are not set.");
  mTransform->mVoltageAbc->set(**mPlant->acTerminalVoltageAttribute());
  mTransform->mCurrentAbc->set(**mPlant->acTerminalCurrentAttribute());
  mTransform->mAngle->set(**mAngle);
  mTransform->mAngularFrequency->set(**mAngularFrequency);
  mTransform->step();
  const Real vBase = std::sqrt(2.0 / 3.0) * mParameters.nominalAcLineLineRms;
  const Real iBase = (2.0 / 3.0) * mParameters.nominalPower / vBase;
  **mVdPu = **mTransform->mVd / vBase;
  **mVqPu = **mTransform->mVq / vBase;
  **mIdPu = **mTransform->mId / iBase;
  **mIqPu = **mTransform->mIq / iBase;
  **mActivePowerPu = **mTransform->mActivePower / mParameters.nominalPower;
  **mReactivePowerPu = **mTransform->mReactivePower / mParameters.nominalPower;
  **mAcPower = **mTransform->mActivePower;
  **mDcPower = **mPlant->dcVoltageAttribute() * **mPlant->dcCurrentAttribute();
  // SSN_MMC uses Pac > 0 for AC-side power transfer and Pdc = Vdc*Idc in
  // the model's DC-current orientation. Pac-Pdc is only a port-power
  // mismatch; it intentionally excludes stored-energy derivatives and
  // represented converter losses.
  **mPowerBalanceError = **mAcPower - **mDcPower;
}

void SSN_MMCStation::filterAndOuterStep() {
  if (!mInitialized)
    throw std::logic_error("MMCStation is not initialized.");
  const Bool enabled =
      **mState == static_cast<Int>(State::Enabled) && **mOuterLoopsEnabled;
  mPFilter->mInput->set(**mActivePowerPu);
  mQFilter->mInput->set(**mReactivePowerPu);
  mVdcFilter->mInput->set(**mPlant->dcVoltageAttribute());
  mPFilter->step();
  mQFilter->step();
  mVdcFilter->step();
  **mFilteredActivePowerPu = **mPFilter->mOutput;
  **mFilteredReactivePowerPu = **mQFilter->mOutput;
  **mFilteredDcVoltage = **mVdcFilter->mOutput;

  mReactiveController->mReference->set(**mReactivePowerReferencePu);
  mReactiveController->mMeasurement->set(**mFilteredReactivePowerPu);
  mReactiveController->mEnable->set(enabled);
  mReactiveController->step();
  **mReactiveOuterUpperSaturated = **mReactiveController->mUpperSaturated;
  **mReactiveOuterLowerSaturated = **mReactiveController->mLowerSaturated;
  if (**mOuterLoopsEnabled)
    // DQsym's reactive PI output is positive for Qref-Q > 0, while this
    // station's documented converter-current convention has positive
    // generated Q for negative iq. The station boundary therefore applies
    // the required current-orientation sign.
    **mIqReferencePu = -**mReactiveController->mOutput;

  Real vdcReference = **mDcVoltageReference;
  if (**mControlMode ==
      static_cast<Int>(ControlMode::ActivePowerReactivePower)) {
    mActiveController->mReference->set(**mActivePowerReferencePu);
    mActiveController->mMeasurement->set(**mFilteredActivePowerPu);
    mActiveController->mEnable->set(enabled);
    mActiveController->step();
    **mActiveOuterUpperSaturated = **mActiveController->mUpperSaturated;
    **mActiveOuterLowerSaturated = **mActiveController->mLowerSaturated;
    vdcReference = **mActiveController->mOutput;
  }
  mDcVoltageController->mReference->set(vdcReference);
  mDcVoltageController->mMeasurement->set(**mFilteredDcVoltage);
  mDcVoltageController->mEnable->set(enabled);
  mDcVoltageController->step();
  **mDcOuterErrorPu = **mDcVoltageController->mError;
  **mDcOuterProportionalContribution = 4.0 * **mDcOuterErrorPu;
  **mDcOuterIntegralContribution = **mDcVoltageController->mIntegratorState;
  **mDcOuterUnsaturatedOutput = **mDcVoltageController->mUnsaturatedOutput;
  **mDcOuterOutput = **mDcVoltageController->mOutput;
  **mDcOuterUpperSaturated = **mDcVoltageController->mUpperSaturated;
  **mDcOuterLowerSaturated = **mDcVoltageController->mLowerSaturated;
  if (**mOuterLoopsEnabled) {
    const Bool directActivePowerOrientation =
        **mControlMode ==
        static_cast<Int>(ControlMode::ActivePowerReactivePower);
    // The SSN_MMC plant's validated DC-voltage convention requires negative
    // iDeltaD for DC overvoltage, whereas the pinned DQsym Vdc PI emits a
    // positive scalar for Vdc_meas-Vdc_ref > 0. Apply that plant-boundary
    // orientation only in direct DC-voltage mode. The active-power cascade
    // retains its independently validated generated-power orientation.
    **mIdReferencePu = directActivePowerOrientation
                           ? **mDcVoltageController->mOutput
                           : -**mDcVoltageController->mOutput;
  }
}

void SSN_MMCStation::currentStep() {
  const Bool enabled = **mState == static_cast<Int>(State::Enabled);
  mCurrentController->mIdReference->set(**mIdReferencePu);
  mCurrentController->mIqReference->set(**mIqReferencePu);
  mCurrentController->mId->set(**mIdPu);
  mCurrentController->mIq->set(**mIqPu);
  mCurrentController->mVd->set(**mVdPu);
  mCurrentController->mVq->set(**mVqPu);
  mCurrentController->mFrequencyHz->set(**mAngularFrequency / (2.0 * PI));
  mCurrentController->mEnable->set(enabled);
  mCurrentController->step();
  **mVdReferencePu = **mCurrentController->mVdReference;
  **mVqReferencePu = **mCurrentController->mVqReference;
  **mCurrentDSaturated = **mCurrentController->mDSaturated;
  **mCurrentQSaturated = **mCurrentController->mQSaturated;
  **mCurrentDUpperSaturated = **mCurrentController->mDUpperSaturated;
  **mCurrentDLowerSaturated = **mCurrentController->mDLowerSaturated;
  **mCurrentQUpperSaturated = **mCurrentController->mQUpperSaturated;
  **mCurrentQLowerSaturated = **mCurrentController->mQLowerSaturated;
}

void SSN_MMCStation::commandStep() {
  mModulation->mVdCommand->set(**mVdReferencePu);
  mModulation->mVqCommand->set(**mVqReferencePu);
  mModulation->mDcVoltage->set(**mPlant->dcVoltageAttribute());
  mModulation->mAngle->set(**mAngle);
  mModulation->step();
  **mConverterPhaseCommand = **mModulation->mAbcCommand;
  **mModulationMagnitude = **mModulation->mModulationMagnitude;
  **mModulationD = **mModulation->mDCommand;
  **mModulationQ = **mModulation->mQCommand;
  **mModulationDUnsaturated = **mModulation->mDUnsaturatedCommand;
  **mModulationQUnsaturated = **mModulation->mQUnsaturatedCommand;
  **mModulationSaturated = **mModulation->mSaturated;
  **mModulationDUpperSaturated = **mModulation->mDUpperSaturated;
  **mModulationDLowerSaturated = **mModulation->mDLowerSaturated;
  **mModulationQUpperSaturated = **mModulation->mQUpperSaturated;
  **mModulationQLowerSaturated = **mModulation->mQLowerSaturated;

  if (**mState == static_cast<Int>(State::Enabled)) {
    const Real vBase = std::sqrt(2.0 / 3.0) * mParameters.nominalAcLineLineRms;

    // Apply the modulation block's saturation to the voltage command that is
    // actually sent to the MMC. Previously the saturated modulation was only
    // logged while the plant still received the unsaturated current-controller
    // reference.
    auto realizedAxisVoltage = [](Real voltageReferencePu,
                                  Real unsaturatedModulation,
                                  Real saturatedModulation, Real baseVoltage) {
      if (std::abs(unsaturatedModulation) <= 1e-12)
        return voltageReferencePu * baseVoltage;
      return voltageReferencePu * baseVoltage *
             (saturatedModulation / unsaturatedModulation);
    };

    (**mPlantDifferentialVoltageCommand)(0, 0) = realizedAxisVoltage(
        **mVdReferencePu, **mModulationDUnsaturated, **mModulationD, vBase);
    (**mPlantDifferentialVoltageCommand)(1, 0) = realizedAxisVoltage(
        **mVqReferencePu, **mModulationQUnsaturated, **mModulationQ, vBase);

    mPlant->setExternalDifferentialVoltageCommand(
        (**mPlantDifferentialVoltageCommand)(0, 0),
        (**mPlantDifferentialVoltageCommand)(1, 0));
  } else {
    **mPlantDifferentialVoltageCommand = mHeldPlantCommand;
  }
}

void SSN_MMCStation::setEnableDiagnostic(EnableDiagnostic diagnostic) {
  **mEnableDiagnostic = static_cast<Int>(diagnostic);
}

SSN_MMCStation::MeasurementTask::MeasurementTask(SSN_MMCStation &station)
    : Task(**station.mName + ".Measurements"), mStation(station) {
  mAttributeDependencies = {station.mPlant->interfaceVoltageAttribute(),
                            station.mPlant->interfaceCurrentAttribute(),
                            station.mAngle, station.mAngularFrequency};
  mModifiedAttributes = {
      station.mVdPu,    station.mVqPu,          station.mIdPu,
      station.mIqPu,    station.mActivePowerPu, station.mReactivePowerPu,
      station.mAcPower, station.mDcPower,       station.mPowerBalanceError};
}

SSN_MMCStation::OuterTask::OuterTask(SSN_MMCStation &station)
    : Task(**station.mName + ".FiltersAndOuter"), mStation(station) {
  mAttributeDependencies = {station.mActivePowerPu,
                            station.mReactivePowerPu,
                            station.mPlant->dcVoltageAttribute(),
                            station.mActivePowerReferencePu,
                            station.mReactivePowerReferencePu,
                            station.mDcVoltageReference,
                            station.mControlMode,
                            station.mState,
                            station.mOuterLoopsEnabled};
  mModifiedAttributes = {station.mFilteredActivePowerPu,
                         station.mFilteredReactivePowerPu,
                         station.mFilteredDcVoltage,
                         station.mIdReferencePu,
                         station.mIqReferencePu,
                         station.mActiveOuterUpperSaturated,
                         station.mActiveOuterLowerSaturated,
                         station.mDcOuterUpperSaturated,
                         station.mDcOuterLowerSaturated,
                         station.mReactiveOuterUpperSaturated,
                         station.mReactiveOuterLowerSaturated,
                         station.mDcOuterErrorPu,
                         station.mDcOuterProportionalContribution,
                         station.mDcOuterIntegralContribution,
                         station.mDcOuterUnsaturatedOutput,
                         station.mDcOuterOutput};
}

SSN_MMCStation::CurrentTask::CurrentTask(SSN_MMCStation &station)
    : Task(**station.mName + ".Current"), mStation(station) {
  mAttributeDependencies = {station.mIdReferencePu,
                            station.mIqReferencePu,
                            station.mIdPu,
                            station.mIqPu,
                            station.mVdPu,
                            station.mVqPu,
                            station.mAngularFrequency,
                            station.mState};
  mModifiedAttributes = {
      station.mVdReferencePu,          station.mVqReferencePu,
      station.mCurrentDSaturated,      station.mCurrentQSaturated,
      station.mCurrentDUpperSaturated, station.mCurrentDLowerSaturated,
      station.mCurrentQUpperSaturated, station.mCurrentQLowerSaturated};
}

SSN_MMCStation::CommandTask::CommandTask(SSN_MMCStation &station)
    : Task(**station.mName + ".Command"), mStation(station) {
  mAttributeDependencies = {station.mVdReferencePu, station.mVqReferencePu,
                            station.mPlant->dcVoltageAttribute(),
                            station.mAngle, station.mState};
  mModifiedAttributes = {
      station.mConverterPhaseCommand,
      station.mPlantDifferentialVoltageCommand,
      station.mModulationMagnitude,
      station.mModulationD,
      station.mModulationQ,
      station.mModulationDUnsaturated,
      station.mModulationQUnsaturated,
      station.mModulationSaturated,
      station.mModulationDUpperSaturated,
      station.mModulationDLowerSaturated,
      station.mModulationQUpperSaturated,
      station.mModulationQLowerSaturated,
      station.mPlant->externalDifferentialVoltageAttribute()};
}

Task::List SSN_MMCStation::getTasks() {
  return {std::make_shared<MeasurementTask>(*this),
          std::make_shared<OuterTask>(*this),
          std::make_shared<CurrentTask>(*this),
          std::make_shared<CommandTask>(*this)};
}
