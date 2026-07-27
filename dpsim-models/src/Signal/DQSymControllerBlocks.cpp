// SPDX-License-Identifier: MPL-2.0
#include <dpsim-models/Signal/DQSymControllerBlocks.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;
using namespace CPS::Signal;

namespace {
void requireFinite(Real value, const String &name) {
  if (!std::isfinite(value))
    throw std::invalid_argument(name + " must be finite.");
}

Real clampValue(Real value, Real lower, Real upper) {
  return std::clamp(value, lower, upper);
}

void advanceConditionalPI(Real kp, Real ki, Real lowerLimit, Real upperLimit,
                          Real timeStep, Bool integrate, Bool enable,
                          Real error, Real feedforward, Real &state,
                          Real &previousIntegratorInput, Real &unsaturated,
                          Real &output, Bool &saturated) {
  if (!enable) {
    previousIntegratorInput = 0.0;
    unsaturated = state + feedforward;
    output = clampValue(unsaturated, lowerLimit, upperLimit);
    saturated = (output != unsaturated);
    return;
  }

  const Real rawBeforeIntegration = kp * error + state + feedforward;
  const Bool hold = (rawBeforeIntegration >= upperLimit && error > 0.0) ||
                    (rawBeforeIntegration <= lowerLimit && error < 0.0);
  const Real integratorInput = hold ? 0.0 : ki * error;
  if (integrate)
    state = clampValue(state + 0.5 * timeStep *
                                   (previousIntegratorInput + integratorInput),
                       lowerLimit, upperLimit);
  previousIntegratorInput = integratorInput;
  unsaturated = kp * error + state + feedforward;
  output = clampValue(unsaturated, lowerLimit, upperLimit);
  saturated = (output != unsaturated);
}
} // namespace

ExternallyAngledDQAdapter::ExternallyAngledDQAdapter(String name,
                                                     Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mVoltageAbc(mAttributes->createDynamic<Matrix>("voltage_abc")),
      mCurrentAbc(mAttributes->createDynamic<Matrix>("current_abc")),
      mAngle(mAttributes->createDynamic<Real>("angle")),
      mAngularFrequency(mAttributes->createDynamic<Real>("angular_frequency")),
      mVd(mAttributes->create<Real>("vd", 0.0)),
      mVq(mAttributes->create<Real>("vq", 0.0)),
      mId(mAttributes->create<Real>("id", 0.0)),
      mIq(mAttributes->create<Real>("iq", 0.0)),
      mActivePower(mAttributes->create<Real>("active_power", 0.0)),
      mReactivePower(mAttributes->create<Real>("reactive_power", 0.0)) {
  **mVoltageAbc = Matrix::Zero(3, 1);
  **mCurrentAbc = Matrix::Zero(3, 1);
  **mAngle = 0.0;
  **mAngularFrequency = 0.0;
}

void ExternallyAngledDQAdapter::step() {
  if ((**mVoltageAbc).rows() != 3 || (**mVoltageAbc).cols() != 1 ||
      !(**mVoltageAbc).allFinite() || (**mCurrentAbc).rows() != 3 ||
      (**mCurrentAbc).cols() != 1 || !(**mCurrentAbc).allFinite() ||
      !std::isfinite(**mAngle) || !std::isfinite(**mAngularFrequency))
    throw std::invalid_argument(
        "Externally angled dq adapter inputs must be finite 3x1 abc vectors "
        "and finite angle/frequency.");

  const Real theta = **mAngle;
  auto transform = [theta](const Matrix &abc, Real &d, Real &q) {
    d = 0.0;
    q = 0.0;
    for (UInt phase = 0; phase < 3; ++phase) {
      const Real phaseAngle = theta - static_cast<Real>(phase) * 2.0 * PI / 3.0;
      d += (2.0 / 3.0) * abc(phase, 0) * std::cos(phaseAngle);
      q -= (2.0 / 3.0) * abc(phase, 0) * std::sin(phaseAngle);
    }
  };
  transform(**mVoltageAbc, **mVd, **mVq);
  transform(**mCurrentAbc, **mId, **mIq);
  **mActivePower = 1.5 * (**mVd * **mId + **mVq * **mIq);
  **mReactivePower = 1.5 * (**mVq * **mId - **mVd * **mIq);
}

Matrix ExternallyAngledDQAdapter::dqToAbc(Real d, Real q) const {
  requireFinite(d, "inverse dq d");
  requireFinite(q, "inverse dq q");
  requireFinite(**mAngle, "inverse dq angle");
  Matrix abc(3, 1);
  for (UInt phase = 0; phase < 3; ++phase) {
    const Real phaseAngle =
        **mAngle - static_cast<Real>(phase) * 2.0 * PI / 3.0;
    abc(phase, 0) = d * std::cos(phaseAngle) - q * std::sin(phaseAngle);
  }
  return abc;
}

ExternallyAngledDQAdapter::StepTask::StepTask(
    ExternallyAngledDQAdapter &adapter)
    : Task(**adapter.mName + ".Step"), mAdapter(adapter) {
  mAttributeDependencies = {adapter.mVoltageAbc, adapter.mCurrentAbc,
                            adapter.mAngle, adapter.mAngularFrequency};
  mModifiedAttributes = {adapter.mVd,          adapter.mVq,
                         adapter.mId,          adapter.mIq,
                         adapter.mActivePower, adapter.mReactivePower};
}

Task::List ExternallyAngledDQAdapter::getTasks() {
  return {std::make_shared<StepTask>(*this)};
}

DQSymPIController::DQSymPIController(String name, Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mError(mAttributes->createDynamic<Real>("error")),
      mFeedforward(mAttributes->createDynamic<Real>("feedforward")),
      mEnable(mAttributes->createDynamic<Bool>("enable")),
      mIntegratorState(mAttributes->create<Real>("integrator_state", 0.0)),
      mUnsaturatedOutput(mAttributes->create<Real>("unsaturated_output", 0.0)),
      mOutput(mAttributes->create<Real>("output", 0.0)),
      mSaturated(mAttributes->create<Bool>("saturated", false)),
      mUpperSaturated(mAttributes->create<Bool>("upper_saturated", false)),
      mLowerSaturated(mAttributes->create<Bool>("lower_saturated", false)) {
  **mError = 0.0;
  **mFeedforward = 0.0;
  **mEnable = true;
}

void DQSymPIController::setParameters(Real kp, Real ki, Real lowerLimit,
                                      Real upperLimit) {
  requireFinite(kp, "Kp");
  requireFinite(ki, "Ki");
  requireFinite(lowerLimit, "lower limit");
  requireFinite(upperLimit, "upper limit");
  if (lowerLimit >= upperLimit)
    throw std::invalid_argument("PI lower limit must be below upper limit.");
  mKp = kp;
  mKi = ki;
  mLowerLimit = lowerLimit;
  mUpperLimit = upperLimit;
  mParametersSet = true;
}

void DQSymOuterController::setInitialIntegratorState(Real state) {
  requireFinite(state, "outer-controller initial integrator state");
  mPI->setInitialState(state);
  **mIntegratorState = state;
}

void DQSymPIController::setInitialState(Real integratorState, Real initialError,
                                        Real initialFeedforward) {
  requireFinite(integratorState, "initial integrator state");
  requireFinite(initialError, "initial error");
  requireFinite(initialFeedforward, "initial feedforward");
  **mIntegratorState = integratorState;
  **mError = initialError;
  **mFeedforward = initialFeedforward;
  mPreviousIntegratorInput = 0.0;
}

void DQSymPIController::initialize(Real timeStep) {
  if (!mParametersSet)
    throw std::logic_error("DQSym PI parameters are not set.");
  requireFinite(timeStep, "PI time step");
  if (timeStep <= 0.0)
    throw std::invalid_argument("PI time step must be positive.");
  mTimeStep = timeStep;
  mInitialized = false;
  step();
  mInitialized = true;
}

void DQSymPIController::step() {
  if (!mParametersSet)
    throw std::logic_error("DQSym PI parameters are not set.");
  const Real error = **mError;
  const Real feedforward = **mFeedforward;
  requireFinite(error, "PI error");
  requireFinite(feedforward, "PI feedforward");
  requireFinite(**mIntegratorState, "PI integrator state");

  advanceConditionalPI(mKp, mKi, mLowerLimit, mUpperLimit, mTimeStep,
                       mInitialized, **mEnable, error, feedforward,
                       **mIntegratorState, mPreviousIntegratorInput,
                       **mUnsaturatedOutput, **mOutput, **mSaturated);
  **mUpperSaturated = **mUnsaturatedOutput > mUpperLimit;
  **mLowerSaturated = **mUnsaturatedOutput < mLowerLimit;
  requireFinite(**mOutput, "PI output");
}

DQSymPIController::StepTask::StepTask(DQSymPIController &controller)
    : Task(**controller.mName + ".Step"), mController(controller) {
  mAttributeDependencies = {controller.mError, controller.mFeedforward,
                            controller.mEnable};
  mModifiedAttributes = {
      controller.mIntegratorState, controller.mUnsaturatedOutput,
      controller.mOutput,          controller.mSaturated,
      controller.mUpperSaturated,  controller.mLowerSaturated};
}

Task::List DQSymPIController::getTasks() {
  return {std::make_shared<StepTask>(*this)};
}

DQSymSecondOrderFilter::DQSymSecondOrderFilter(String name,
                                               Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mInput(mAttributes->createDynamic<Real>("input")),
      mState(mAttributes->create<Matrix>("state", Matrix::Zero(2, 1))),
      mOutput(mAttributes->create<Real>("output", 0.0)) {}

void DQSymSecondOrderFilter::setParameters(Real naturalFrequencyHz,
                                           Real dampingRatio) {
  requireFinite(naturalFrequencyHz, "filter natural frequency");
  requireFinite(dampingRatio, "filter damping ratio");
  if (naturalFrequencyHz <= 0.0 || dampingRatio <= 0.0)
    throw std::invalid_argument(
        "Filter frequency and damping ratio must be positive.");
  mNaturalFrequencyHz = naturalFrequencyHz;
  mDampingRatio = dampingRatio;
  mParametersSet = true;
}

void DQSymSecondOrderFilter::setInitialValue(Real value) {
  requireFinite(value, "filter initial value");
  (**mState)(0, 0) = value;
  (**mState)(1, 0) = 0.0;
  **mOutput = value;
  **mInput = value;
  mPreviousState = **mState;
  mPreviousInput = value;
}

void DQSymSecondOrderFilter::initialize(Real timeStep) {
  if (!mParametersSet)
    throw std::logic_error("DQSym filter parameters are not set.");
  requireFinite(timeStep, "filter time step");
  if (timeStep <= 0.0)
    throw std::invalid_argument("Filter time step must be positive.");
  mTimeStep = timeStep;
  mPreviousState = **mState;
  mPreviousInput = **mInput;
  mInitialized = true;
}

void DQSymSecondOrderFilter::step() {
  if (!mInitialized)
    throw std::logic_error("DQSym filter is not initialized.");
  const Real input = **mInput;
  requireFinite(input, "filter input");
  const Real omega = 2.0 * PI * mNaturalFrequencyHz;
  Matrix a(2, 2);
  a << 0.0, 1.0, -omega * omega, -2.0 * mDampingRatio * omega;
  Matrix b(2, 1);
  b << 0.0, omega * omega;
  **mState = Math::StateSpaceTrapezoidal(
      mPreviousState, a, b, mTimeStep, Matrix::Constant(1, 1, input),
      Matrix::Constant(1, 1, mPreviousInput));
  **mOutput = (**mState)(0, 0);
  if (!(**mState).allFinite())
    throw std::runtime_error("DQSym filter state contains NaN or Inf.");
  mPreviousState = **mState;
  mPreviousInput = input;
}

DQSymSecondOrderFilter::StepTask::StepTask(DQSymSecondOrderFilter &filter)
    : Task(**filter.mName + ".Step"), mFilter(filter) {
  mAttributeDependencies = {filter.mInput};
  mModifiedAttributes = {filter.mState, filter.mOutput};
}

Task::List DQSymSecondOrderFilter::getTasks() {
  return {std::make_shared<StepTask>(*this)};
}

DQSymOuterController::DQSymOuterController(String name, Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mReference(mAttributes->createDynamic<Real>("reference")),
      mMeasurement(mAttributes->createDynamic<Real>("measurement")),
      mEnable(mAttributes->createDynamic<Bool>("enable")),
      mError(mAttributes->create<Real>("error", 0.0)),
      mIntegratorState(mAttributes->create<Real>("integrator_state", 0.0)),
      mUnsaturatedOutput(mAttributes->create<Real>("unsaturated_output", 0.0)),
      mOutput(mAttributes->create<Real>("output", 0.0)),
      mSaturated(mAttributes->create<Bool>("saturated", false)),
      mUpperSaturated(mAttributes->create<Bool>("upper_saturated", false)),
      mLowerSaturated(mAttributes->create<Bool>("lower_saturated", false)),
      mPI(DQSymPIController::make(name + ".PI", logLevel)) {
  **mEnable = true;
}

void DQSymOuterController::setParameters(DQSymOuterLoopType type, Real kp,
                                         Real ki, Real lowerLimit,
                                         Real upperLimit, Real normalization,
                                         Real outputScale,
                                         Real initialIntegratorState) {
  if (!std::isfinite(normalization) || normalization <= 0.0 ||
      !std::isfinite(outputScale))
    throw std::invalid_argument(
        "Invalid DQSym outer-controller normalization or scale.");
  mType = type;
  mNormalization = normalization;
  mOutputScale = outputScale;
  mPI->setParameters(kp, ki, lowerLimit, upperLimit);
  mPI->setInitialState(initialIntegratorState);
  mParametersSet = true;
}

void DQSymOuterController::initialize(Real timeStep) {
  if (!mParametersSet)
    throw std::logic_error("DQSym outer-controller parameters are not set.");
  mPI->initialize(timeStep);
  step();
}

void DQSymOuterController::step() {
  if (!mParametersSet)
    throw std::logic_error("DQSym outer-controller parameters are not set.");
  Real error = 0.0;
  if (**mEnable) {
    if (mType == DQSymOuterLoopType::ReactivePowerToQCurrent)
      error = (**mReference - **mMeasurement) / mNormalization;
    else
      error = (**mMeasurement - **mReference) / mNormalization;
  }
  **mError = error;
  mPI->mError->set(error);
  mPI->mEnable->set(**mEnable);
  mPI->step();
  **mIntegratorState = **mPI->mIntegratorState;
  **mUnsaturatedOutput = mOutputScale * **mPI->mUnsaturatedOutput;
  **mOutput = mOutputScale * **mPI->mOutput;
  **mSaturated = **mPI->mSaturated;
  **mUpperSaturated =
      mOutputScale >= 0.0 ? **mPI->mUpperSaturated : **mPI->mLowerSaturated;
  **mLowerSaturated =
      mOutputScale >= 0.0 ? **mPI->mLowerSaturated : **mPI->mUpperSaturated;
}

DQSymOuterController::StepTask::StepTask(DQSymOuterController &controller)
    : Task(**controller.mName + ".Step"), mController(controller) {
  mAttributeDependencies = {controller.mReference, controller.mMeasurement,
                            controller.mEnable};
  mModifiedAttributes = {controller.mError,
                         controller.mIntegratorState,
                         controller.mUnsaturatedOutput,
                         controller.mOutput,
                         controller.mSaturated,
                         controller.mUpperSaturated,
                         controller.mLowerSaturated};
}

Task::List DQSymOuterController::getTasks() {
  return {std::make_shared<StepTask>(*this)};
}

DQSymCurrentController::DQSymCurrentController(String name,
                                               Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mIdReference(mAttributes->createDynamic<Real>("id_reference")),
      mIqReference(mAttributes->createDynamic<Real>("iq_reference")),
      mId(mAttributes->createDynamic<Real>("id")),
      mIq(mAttributes->createDynamic<Real>("iq")),
      mVd(mAttributes->createDynamic<Real>("vd")),
      mVq(mAttributes->createDynamic<Real>("vq")),
      mFrequencyHz(mAttributes->createDynamic<Real>("frequency_hz")),
      mEnable(mAttributes->createDynamic<Bool>("enable")),
      mDIntegratorState(mAttributes->create<Real>("d_integrator_state", 0.0)),
      mQIntegratorState(mAttributes->create<Real>("q_integrator_state", 0.0)),
      mDUnsaturatedReference(
          mAttributes->create<Real>("d_unsaturated_reference", 0.0)),
      mQUnsaturatedReference(
          mAttributes->create<Real>("q_unsaturated_reference", 0.0)),
      mVdReference(mAttributes->create<Real>("vd_reference", 0.0)),
      mVqReference(mAttributes->create<Real>("vq_reference", 0.0)),
      mDSaturated(mAttributes->create<Bool>("d_saturated", false)),
      mQSaturated(mAttributes->create<Bool>("q_saturated", false)),
      mDUpperSaturated(mAttributes->create<Bool>("d_upper_saturated", false)),
      mDLowerSaturated(mAttributes->create<Bool>("d_lower_saturated", false)),
      mQUpperSaturated(mAttributes->create<Bool>("q_upper_saturated", false)),
      mQLowerSaturated(mAttributes->create<Bool>("q_lower_saturated", false)) {
  **mEnable = true;
}

void DQSymCurrentController::setParameters(
    const DQSymCurrentControllerParameters &parameters) {
  requireFinite(parameters.kp, "current Kp");
  requireFinite(parameters.ki, "current Ki");
  requireFinite(parameters.rFeedforwardPu, "current R feedforward");
  requireFinite(parameters.lFeedforwardPu, "current L feedforward");
  requireFinite(parameters.nominalFrequencyHz, "nominal frequency");
  if (parameters.nominalFrequencyHz <= 0.0 ||
      parameters.lowerLimitPu >= parameters.upperLimitPu)
    throw std::invalid_argument("Invalid DQSym current-controller parameters.");
  mParameters = parameters;
  mParametersSet = true;
}

void DQSymCurrentController::setInitialIntegratorStates(Real dState,
                                                        Real qState) {
  requireFinite(dState, "d integrator initial state");
  requireFinite(qState, "q integrator initial state");
  **mDIntegratorState = dState;
  **mQIntegratorState = qState;
}

void DQSymCurrentController::initialize(Real timeStep) {
  if (!mParametersSet)
    throw std::logic_error("DQSym current-controller parameters are not set.");
  if (!std::isfinite(timeStep) || timeStep <= 0.0)
    throw std::invalid_argument(
        "Current-controller time step must be finite and positive.");
  mTimeStep = timeStep;
  mInitialized = true;
}

void DQSymCurrentController::step() {
  if (!mInitialized)
    throw std::logic_error("DQSym current controller is not initialized.");
  const Real ed = **mIdReference - **mId;
  const Real eq = **mIqReference - **mIq;
  const Real omegaPu = **mFrequencyHz / mParameters.nominalFrequencyHz;
  const Real dFeedforward = **mVd + mParameters.rFeedforwardPu * **mId -
                            omegaPu * mParameters.lFeedforwardPu * **mIq;
  const Real qFeedforward = **mVq + mParameters.rFeedforwardPu * **mIq +
                            omegaPu * mParameters.lFeedforwardPu * **mId;

  auto updateAxis =
      [&](Real error, Real feedforward, Real &previousInput,
          Attribute<Real>::Ptr state, Attribute<Real>::Ptr unsaturated,
          Attribute<Real>::Ptr output, Attribute<Bool>::Ptr saturated,
          Attribute<Bool>::Ptr upper, Attribute<Bool>::Ptr lower) {
        advanceConditionalPI(mParameters.kp, mParameters.ki,
                             mParameters.lowerLimitPu, mParameters.upperLimitPu,
                             mTimeStep, true, **mEnable, error, feedforward,
                             **state, previousInput, **unsaturated, **output,
                             **saturated);
        **upper = **unsaturated > mParameters.upperLimitPu;
        **lower = **unsaturated < mParameters.lowerLimitPu;
      };

  updateAxis(ed, dFeedforward, mPreviousDIntegratorInput, mDIntegratorState,
             mDUnsaturatedReference, mVdReference, mDSaturated,
             mDUpperSaturated, mDLowerSaturated);
  updateAxis(eq, qFeedforward, mPreviousQIntegratorInput, mQIntegratorState,
             mQUnsaturatedReference, mVqReference, mQSaturated,
             mQUpperSaturated, mQLowerSaturated);
}

DQSymCurrentController::StepTask::StepTask(DQSymCurrentController &controller)
    : Task(**controller.mName + ".Step"), mController(controller) {
  mAttributeDependencies = {controller.mIdReference, controller.mIqReference,
                            controller.mId,          controller.mIq,
                            controller.mVd,          controller.mVq,
                            controller.mFrequencyHz, controller.mEnable};
  mModifiedAttributes = {controller.mDIntegratorState,
                         controller.mQIntegratorState,
                         controller.mDUnsaturatedReference,
                         controller.mQUnsaturatedReference,
                         controller.mVdReference,
                         controller.mVqReference,
                         controller.mDSaturated,
                         controller.mQSaturated,
                         controller.mDUpperSaturated,
                         controller.mDLowerSaturated,
                         controller.mQUpperSaturated,
                         controller.mQLowerSaturated};
}

Task::List DQSymCurrentController::getTasks() {
  return {std::make_shared<StepTask>(*this)};
}

DQSymModulation::DQSymModulation(String name, Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mVdCommand(mAttributes->createDynamic<Real>("vd_command")),
      mVqCommand(mAttributes->createDynamic<Real>("vq_command")),
      mDcVoltage(mAttributes->createDynamic<Real>("dc_voltage")),
      mAngle(mAttributes->createDynamic<Real>("angle")),
      mDCommand(mAttributes->create<Real>("d_command", 0.0)),
      mQCommand(mAttributes->create<Real>("q_command", 0.0)),
      mDUnsaturatedCommand(
          mAttributes->create<Real>("d_unsaturated_command", 0.0)),
      mQUnsaturatedCommand(
          mAttributes->create<Real>("q_unsaturated_command", 0.0)),
      mModulationMagnitude(
          mAttributes->create<Real>("modulation_magnitude", 0.0)),
      mAbcCommand(
          mAttributes->create<Matrix>("abc_command", Matrix::Zero(3, 1))),
      mSaturated(mAttributes->create<Bool>("saturated", false)),
      mDUpperSaturated(mAttributes->create<Bool>("d_upper_saturated", false)),
      mDLowerSaturated(mAttributes->create<Bool>("d_lower_saturated", false)),
      mQUpperSaturated(mAttributes->create<Bool>("q_upper_saturated", false)),
      mQLowerSaturated(mAttributes->create<Bool>("q_lower_saturated", false)) {}

void DQSymModulation::setParameters(Real nominalDcVoltage,
                                    Real nominalAcLineLineRms,
                                    Real lowerAxisLimit, Real upperAxisLimit) {
  if (!std::isfinite(nominalDcVoltage) ||
      !std::isfinite(nominalAcLineLineRms) || nominalDcVoltage <= 0.0 ||
      nominalAcLineLineRms <= 0.0 || lowerAxisLimit >= upperAxisLimit)
    throw std::invalid_argument("Invalid DQSym modulation parameters.");
  mNominalDcVoltage = nominalDcVoltage;
  mNominalAcLineLineRms = nominalAcLineLineRms;
  mLowerAxisLimit = lowerAxisLimit;
  mUpperAxisLimit = upperAxisLimit;
  mParametersSet = true;
}

void DQSymModulation::step() {
  if (!mParametersSet)
    throw std::logic_error("DQSym modulation parameters are not set.");
  const Real dcVoltage = clampValue(**mDcVoltage, 0.75 * mNominalDcVoltage,
                                    1.25 * mNominalDcVoltage);
  const Real scale =
      (mNominalDcVoltage / dcVoltage) /
      (mNominalDcVoltage * 0.5 * std::sqrt(3.0 / 2.0) / mNominalAcLineLineRms);
  const Real rawD = scale * **mVdCommand;
  const Real rawQ = scale * **mVqCommand;
  **mDUnsaturatedCommand = rawD;
  **mQUnsaturatedCommand = rawQ;
  **mDCommand = clampValue(rawD, mLowerAxisLimit, mUpperAxisLimit);
  **mQCommand = clampValue(rawQ, mLowerAxisLimit, mUpperAxisLimit);
  **mDUpperSaturated = rawD > mUpperAxisLimit;
  **mDLowerSaturated = rawD < mLowerAxisLimit;
  **mQUpperSaturated = rawQ > mUpperAxisLimit;
  **mQLowerSaturated = rawQ < mLowerAxisLimit;
  **mSaturated = (**mDCommand != rawD || **mQCommand != rawQ);
  **mModulationMagnitude = std::hypot(**mDCommand, **mQCommand);

  const Real theta = **mAngle;
  for (UInt phase = 0; phase < 3; ++phase) {
    const Real phaseAngle = theta - static_cast<Real>(phase) * 2.0 * PI / 3.0;
    (**mAbcCommand)(phase, 0) =
        **mDCommand * std::cos(phaseAngle) - **mQCommand * std::sin(phaseAngle);
  }
  if (!(**mAbcCommand).allFinite())
    throw std::runtime_error("DQSym modulation output contains NaN or Inf.");
}

DQSymModulation::StepTask::StepTask(DQSymModulation &modulation)
    : Task(**modulation.mName + ".Step"), mModulation(modulation) {
  mAttributeDependencies = {modulation.mVdCommand, modulation.mVqCommand,
                            modulation.mDcVoltage, modulation.mAngle};
  mModifiedAttributes = {modulation.mDCommand,
                         modulation.mQCommand,
                         modulation.mDUnsaturatedCommand,
                         modulation.mQUnsaturatedCommand,
                         modulation.mModulationMagnitude,
                         modulation.mAbcCommand,
                         modulation.mSaturated,
                         modulation.mDUpperSaturated,
                         modulation.mDLowerSaturated,
                         modulation.mQUpperSaturated,
                         modulation.mQLowerSaturated};
}

Task::List DQSymModulation::getTasks() {
  return {std::make_shared<StepTask>(*this)};
}
