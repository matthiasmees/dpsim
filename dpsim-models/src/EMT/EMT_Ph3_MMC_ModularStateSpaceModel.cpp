// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/EMT/EMT_Ph3_MMC_ModularStateSpaceModel.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include <Eigen/QR>

using namespace CPS;
using namespace CPS::EMT::Ph3;

MMC_ModularStateSpaceModel::MMC_ModularStateSpaceModel()
    : mParametersSet(false), mConfigurationRevision(0), mNominalFrequency(0.0),
      mOmegaN(0.0), mNominalAcVoltage(0.0), mNominalDcVoltage(0.0),
      mArmInductance(0.0), mArmResistance(0.0), mSubmoduleCapacitance(0.0),
      mNumberOfSubmodules(0), mReactorInductance(0.0), mReactorResistance(0.0),
      mInitialAngle(0.0), mInitialOperatingPointEnabled(false),
      mInitialActivePower(0.0), mInitialReactivePower(0.0),
      mExternalDifferentialVoltage(Matrix::Zero(2, 1)),
      mExternalCommonModeVoltage(Matrix::Zero(3, 1)), mMinimumDcVoltage(1.0),
      mJacobianRelativeStep(1e-6), mJacobianAbsoluteStep(1e-8),
      mSparseLinearizationTolerance(1e-14), mStructuredLinearization(true),
      mOperatingPointInitializationEnabled(true),
      mOperatingPointMaximumIterations(40),
      mOperatingPointNormalizedTolerance(1e-9), mFeedforwardStepCounter(0),
      mFeedforwardInitialized(false),
      mSampledFeedforwardState(Matrix::Zero(3, 1)) {}

void MMC_ModularStateSpaceModel::setParameters(
    Real nominalFrequency, Real nominalAcVoltage, Real nominalDcVoltage,
    Real armInductance, Real armResistance, Real submoduleCapacitance,
    UInt numberOfSubmodules, Real reactorInductance, Real reactorResistance) {
  if (!std::isfinite(nominalFrequency) || nominalFrequency <= 0.0 ||
      !std::isfinite(nominalAcVoltage) || nominalAcVoltage <= 0.0 ||
      !std::isfinite(nominalDcVoltage) || nominalDcVoltage <= 0.0 ||
      !std::isfinite(armInductance) || armInductance <= 0.0 ||
      !std::isfinite(armResistance) || armResistance < 0.0 ||
      !std::isfinite(submoduleCapacitance) || submoduleCapacitance <= 0.0 ||
      numberOfSubmodules == 0 || !std::isfinite(reactorInductance) ||
      reactorInductance < 0.0 || !std::isfinite(reactorResistance) ||
      reactorResistance < 0.0)
    throw std::invalid_argument("MMC electrical parameters are invalid.");

  mNominalFrequency = nominalFrequency;
  mOmegaN = 2.0 * PI * nominalFrequency;
  mNominalAcVoltage = nominalAcVoltage;
  mNominalDcVoltage = nominalDcVoltage;
  mArmInductance = armInductance;
  mArmResistance = armResistance;
  mSubmoduleCapacitance = submoduleCapacitance;
  mNumberOfSubmodules = numberOfSubmodules;
  mReactorInductance = reactorInductance;
  mReactorResistance = reactorResistance;

  const Real lEq = 0.5 * armInductance + reactorInductance;
  const Real rEq = 0.5 * armResistance + reactorResistance;
  const Real energyReference = 3.0 * submoduleCapacitance * nominalDcVoltage *
                               nominalDcVoltage /
                               static_cast<Real>(numberOfSubmodules);
  mControllers.setElectricalParameters(mOmegaN, rEq, lEq, armResistance,
                                       armInductance, nominalDcVoltage,
                                       energyReference);
  mParametersSet = true;
  mFeedforwardStepCounter = 0;
  mFeedforwardInitialized = false;
  mSampledFeedforwardState = Matrix::Zero(3, 1);

  touchConfiguration();
}

void MMC_ModularStateSpaceModel::setInitialAngle(Real angle) {
  if (!std::isfinite(angle))
    throw std::invalid_argument("Initial MMC angle must be finite.");
  mInitialAngle = angle;

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setInitialOperatingPoint(Real activePower,
                                                          Real reactivePower) {
  if (!std::isfinite(activePower) || !std::isfinite(reactivePower))
    throw std::invalid_argument("Initial MMC P/Q must be finite.");
  mInitialOperatingPointEnabled = true;
  mInitialActivePower = activePower;
  mInitialReactivePower = reactivePower;

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::clearInitialOperatingPoint() {
  mInitialOperatingPointEnabled = false;

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setPLL(Real kp, Real ki, Bool enabled) {
  mControllers.setPLL(kp, ki, enabled);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setOutputCurrentController(Real kp, Real ki) {
  mControllers.setMMCOutputCurrentController(kp, ki);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setCirculatingCurrentController(Real kp,
                                                                 Real ki) {
  mControllers.setMMCCirculatingCurrentController(kp, ki);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setZeroSequenceCurrentController(Real kp,
                                                                  Real ki) {
  mControllers.setMMCZeroSequenceCurrentController(kp, ki);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setEnergyController(Real kp, Real ki,
                                                     Bool enabled) {
  mControllers.setMMCEnergyController(kp, ki, enabled);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setActivePowerControl(Real reference, Real kp,
                                                       Real ki) {
  mControllers.setActivePowerControl(reference, kp, ki);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setActivePowerFeedforwardControl(
    Real reference, Real cutoffFrequency, Real sampleTime,
    Real minimumDaxisVoltage) {
  validateConfigured();
  if (minimumDaxisVoltage <= 0.0)
    minimumDaxisVoltage = 0.1 * std::sqrt(2.0 / 3.0) * mNominalAcVoltage;
  mSampledFeedforward.setParameters(cutoffFrequency, sampleTime,
                                    minimumDaxisVoltage);
  mControllers.setActivePowerFeedforwardControl(
      reference, cutoffFrequency, sampleTime, minimumDaxisVoltage);
  mFeedforwardInitialized = false;
  mFeedforwardStepCounter = 0;

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setActivePowerReference(Real reference) {
  mControllers.setActivePowerReference(reference);

  // In sampled-feedforward mode a new power reference enters the continuous
  // MMC equations only at the next sample instant through the held current
  // reference. Rebuilding here would apply the reference one EMT step too
  // early and would defeat affine-model reuse between controller samples.
  if (!mControllers.sampledFeedforwardEnabled())
    touchConfiguration();
}
void MMC_ModularStateSpaceModel::setDcVoltageControl(Real reference, Real kp,
                                                     Real ki) {
  mControllers.setDcVoltageControl(reference, kp, ki);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setDcDroopControl(Real activePowerReference,
                                                   Real dcVoltageReference,
                                                   Real droopGain) {
  mControllers.setDcDroopControl(activePowerReference, dcVoltageReference,
                                 droopGain);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setActiveControlOpenLoop(
    Real currentReference) {
  mControllers.setActiveOpenLoop(currentReference);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setReactivePowerControl(Real reference,
                                                         Real kp, Real ki) {
  mControllers.setReactivePowerControl(reference, kp, ki);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setAcVoltageControl(Real reference, Real kp,
                                                     Real ki) {
  mControllers.setAcVoltageControl(reference, kp, ki);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setReactiveControlOpenLoop(
    Real currentReference) {
  mControllers.setReactiveOpenLoop(currentReference);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setCirculatingCurrentReferences(
    Real dReference, Real qReference, Real zReference) {
  mControllers.setCirculatingCurrentReferences(dReference, qReference,
                                               zReference);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setMeasurementFilters(
    Real acVoltageDqTimeConstant, Real activePowerTimeConstant,
    Real reactivePowerTimeConstant, Real dcVoltageTimeConstant,
    Real acVoltageMagnitudeTimeConstant) {
  mControllers.setMeasurementFilters(
      acVoltageDqTimeConstant, activePowerTimeConstant,
      reactivePowerTimeConstant, dcVoltageTimeConstant,
      acVoltageMagnitudeTimeConstant);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setModulationDelay(Real timeDelay,
                                                    UInt padeOrder) {
  if (padeOrder != 2)
    throw std::invalid_argument(
        "MMC model supports second-order Padé delay only.");
  mControllers.setModulationDelay(timeDelay, timeDelay > 0.0);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setLimits(Real maximumAcCurrent,
                                           Real maximumCirculatingCurrent,
                                           Real maximumModulationMagnitude) {
  mControllers.setLimits(maximumAcCurrent, maximumCirculatingCurrent,
                         maximumModulationMagnitude);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setControlSource(ControlSource source) {
  mControllers.setControlSource(source);

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setExternalDifferentialVoltageCommand(
    Real dVolts, Real qVolts) {
  Matrix command(2, 1);
  command << dVolts, qVolts;
  if (command.isApprox(mExternalDifferentialVoltage, 0.0))
    return;
  mExternalDifferentialVoltage = command;
  mControllers.setExternalDifferentialVoltage(dVolts, qVolts);
  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setExternalCommonModeVoltageCommand(
    Real dVolts, Real qVolts, Real zVolts) {
  Matrix command(3, 1);
  command << dVolts, qVolts, zVolts;
  if (command.isApprox(mExternalCommonModeVoltage, 0.0))
    return;
  mExternalCommonModeVoltage = command;
  mControllers.setExternalCommonModeVoltage(dVolts, qVolts, zVolts);
  touchConfiguration();
}
MMC_ModularStateSpaceModel::ControlSource
MMC_ModularStateSpaceModel::controlSource() const {
  return mControllers.controlSource();
}
void MMC_ModularStateSpaceModel::setNumericalLinearizationParameters(
    Real relativeStep, Real absoluteStep) {
  if (!std::isfinite(relativeStep) || relativeStep <= 0.0 ||
      !std::isfinite(absoluteStep) || absoluteStep <= 0.0)
    throw std::invalid_argument("MMC Jacobian steps must be positive.");
  mJacobianRelativeStep = relativeStep;
  mJacobianAbsoluteStep = absoluteStep;

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setSparseLinearizationTolerance(
    Real tolerance) {
  if (!std::isfinite(tolerance) || tolerance < 0.0)
    throw std::invalid_argument(
        "MMC sparse-linearization tolerance must be finite and non-negative.");
  mSparseLinearizationTolerance = tolerance;

  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setStructuredLinearization(Bool enabled) {
  mStructuredLinearization = enabled;
  touchConfiguration();
}
void MMC_ModularStateSpaceModel::setOperatingPointInitialization(
    Bool enabled, UInt maximumIterations, Real normalizedTolerance) {
  if (maximumIterations == 0 || !std::isfinite(normalizedTolerance) ||
      normalizedTolerance <= 0.0)
    throw std::invalid_argument("MMC operating-point settings are invalid.");
  mOperatingPointInitializationEnabled = enabled;
  mOperatingPointMaximumIterations = maximumIterations;
  mOperatingPointNormalizedTolerance = normalizedTolerance;

  touchConfiguration();
}

void MMC_ModularStateSpaceModel::touchConfiguration() {
  ++mConfigurationRevision;
}

Bool MMC_ModularStateSpaceModel::parametersSet() const {
  return mParametersSet;
}
std::uint64_t MMC_ModularStateSpaceModel::configurationRevision() const {
  return mConfigurationRevision;
}
Real MMC_ModularStateSpaceModel::nominalFrequency() const {
  validateConfigured();
  return mNominalFrequency;
}

Real MMC_ModularStateSpaceModel::nominalDcVoltage() const {
  validateConfigured();
  return mNominalDcVoltage;
}
UInt MMC_ModularStateSpaceModel::controllerOffset() const {
  return PlantStateCount;
}
UInt MMC_ModularStateSpaceModel::stateSize() const {
  return PlantStateCount + mControllers.stateSize();
}
UInt MMC_ModularStateSpaceModel::inputSize() const { return InputCount; }
UInt MMC_ModularStateSpaceModel::outputSize() const { return OutputCount; }
std::vector<String> MMC_ModularStateSpaceModel::stateNames() const {
  std::vector<String> result = {
      "plant.iDelta_d",  "plant.iDelta_q",   "plant.iSigma_z",
      "plant.iSigma_d",  "plant.iSigma_q",   "plant.vCDelta_d",
      "plant.vCDelta_q", "plant.vCDelta_Zd", "plant.vCDelta_Zq",
      "plant.vCSigma_d", "plant.vCSigma_q",  "plant.vCSigma_z",
      "plant.grid_angle"};
  for (const auto &name : mControllers.stateNames())
    result.push_back("controller." + name);
  return result;
}

Matrix MMC_ModularStateSpaceModel::abcToDq(const Matrix &abc, Real theta) {
  if (abc.rows() != 3 || abc.cols() != 1)
    throw std::invalid_argument("abcToDq expects 3x1.");
  const Real c0 = std::cos(theta);
  const Real c1 = std::cos(theta - 2.0 * PI / 3.0);
  const Real c2 = std::cos(theta + 2.0 * PI / 3.0);
  const Real s0 = std::sin(theta);
  const Real s1 = std::sin(theta - 2.0 * PI / 3.0);
  const Real s2 = std::sin(theta + 2.0 * PI / 3.0);
  Matrix dq(2, 1);
  dq(0, 0) = (2.0 / 3.0) * (c0 * abc(0, 0) + c1 * abc(1, 0) + c2 * abc(2, 0));
  dq(1, 0) = -(2.0 / 3.0) * (s0 * abc(0, 0) + s1 * abc(1, 0) + s2 * abc(2, 0));
  return dq;
}
Matrix MMC_ModularStateSpaceModel::dqToAbc(Real d, Real q, Real theta) {
  Matrix abc(3, 1);
  abc(0, 0) = d * std::cos(theta) - q * std::sin(theta);
  abc(1, 0) = d * std::cos(theta - 2.0 * PI / 3.0) -
              q * std::sin(theta - 2.0 * PI / 3.0);
  abc(2, 0) = d * std::cos(theta + 2.0 * PI / 3.0) -
              q * std::sin(theta + 2.0 * PI / 3.0);
  return abc;
}
Real MMC_ModularStateSpaceModel::regularizeSigned(Real value,
                                                  Real minimumMagnitude) {
  if (std::abs(value) >= minimumMagnitude)
    return value;
  return value >= 0.0 ? minimumMagnitude : -minimumMagnitude;
}
void MMC_ModularStateSpaceModel::validateConfigured() const {
  if (!mParametersSet)
    throw std::logic_error(
        "MMC electrical parameters have not been configured.");
}
void MMC_ModularStateSpaceModel::validateStateInput(const Matrix &state,
                                                    const Matrix &input) const {
  validateConfigured();
  if (state.rows() != static_cast<Eigen::Index>(stateSize()) ||
      state.cols() != 1 || !state.allFinite())
    throw std::invalid_argument(
        "MMC model state has invalid dimensions or values.");
  if (input.rows() != InputCount || input.cols() != 1 || !input.allFinite())
    throw std::invalid_argument("MMC model input must be finite 5x1.");
}

Real MMC_ModularStateSpaceModel::calculateStoredEnergy(
    const Matrix &state) const {
  if (state.rows() < PlantStateCount || state.cols() != 1)
    throw std::invalid_argument(
        "MMC state is too small for energy calculation.");
  const Real scale = 3.0 * mSubmoduleCapacitance /
                     (2.0 * static_cast<Real>(mNumberOfSubmodules));
  const Real sumSquares = state(VCDeltaD, 0) * state(VCDeltaD, 0) +
                          state(VCDeltaQ, 0) * state(VCDeltaQ, 0) +
                          state(VCDeltaZd, 0) * state(VCDeltaZd, 0) +
                          state(VCDeltaZq, 0) * state(VCDeltaZq, 0) +
                          state(VCSigmaD, 0) * state(VCSigmaD, 0) +
                          state(VCSigmaQ, 0) * state(VCSigmaQ, 0) +
                          2.0 * state(VCSigmaZ, 0) * state(VCSigmaZ, 0);
  return scale * sumSquares;
}
Real MMC_ModularStateSpaceModel::activePower(const Matrix &state,
                                             const Matrix &input) const {
  validateStateInput(state, input);
  const Matrix vdq = abcToDq(input.block(0, 0, 3, 1), state(GridAngle, 0));
  return 1.5 * (vdq(0, 0) * state(IDeltaD, 0) + vdq(1, 0) * state(IDeltaQ, 0));
}
Real MMC_ModularStateSpaceModel::reactivePower(const Matrix &state,
                                               const Matrix &input) const {
  validateStateInput(state, input);
  const Matrix vdq = abcToDq(input.block(0, 0, 3, 1), state(GridAngle, 0));
  return 1.5 * (-vdq(0, 0) * state(IDeltaQ, 0) + vdq(1, 0) * state(IDeltaD, 0));
}
Real MMC_ModularStateSpaceModel::dcVoltage(const Matrix &input) const {
  if (input.rows() != InputCount || input.cols() != 1)
    throw std::invalid_argument("MMC input must be 5x1.");
  return input(Vdcp, 0) - input(Vdcn, 0);
}
Real MMC_ModularStateSpaceModel::dcCurrent(const Matrix &state) const {
  if (state.rows() != static_cast<Eigen::Index>(stateSize()) ||
      state.cols() != 1)
    throw std::invalid_argument("MMC state dimensions are invalid.");
  return 3.0 * state(ISigmaZ, 0);
}

Real MMC_ModularStateSpaceModel::powerBalanceError(const Matrix &state,
                                                   const Matrix &input) const {
  validateStateInput(state, input);
  const Real pAc = activePower(state, input);
  const Real pDc = dcVoltage(input) * dcCurrent(state);
  const Real iD = state(IDeltaD, 0);
  const Real iQ = state(IDeltaQ, 0);
  const Real differentialLoss =
      1.5 * (mArmResistance / 2.0 + mReactorResistance) * (iD * iD + iQ * iQ);
  const Real circulatingLoss =
      3.0 * mArmResistance *
          (state(ISigmaD, 0) * state(ISigmaD, 0) +
           state(ISigmaQ, 0) * state(ISigmaQ, 0)) +
      6.0 * mArmResistance * state(ISigmaZ, 0) * state(ISigmaZ, 0);
  return pAc - pDc + differentialLoss + circulatingLoss;
}

Matrix
MMC_ModularStateSpaceModel::makeControllerInput(const Matrix &state,
                                                const Matrix &input) const {
  const Matrix vdq = abcToDq(input.block(0, 0, 3, 1), state(GridAngle, 0));
  Matrix result = Matrix::Zero(Signal::MMCControllerSystem::InputCount, 1);
  result(Signal::MMCControllerSystem::VGridD, 0) = vdq(0, 0);
  result(Signal::MMCControllerSystem::VGridQ, 0) = vdq(1, 0);
  result(Signal::MMCControllerSystem::ActivePower, 0) =
      1.5 * (vdq(0, 0) * state(IDeltaD, 0) + vdq(1, 0) * state(IDeltaQ, 0));
  result(Signal::MMCControllerSystem::ReactivePower, 0) =
      1.5 * (-vdq(0, 0) * state(IDeltaQ, 0) + vdq(1, 0) * state(IDeltaD, 0));
  result(Signal::MMCControllerSystem::DcVoltage, 0) =
      regularizeSigned(input(Vdcp, 0) - input(Vdcn, 0), mMinimumDcVoltage);
  result(Signal::MMCControllerSystem::IDeltaD, 0) = state(IDeltaD, 0);
  result(Signal::MMCControllerSystem::IDeltaQ, 0) = state(IDeltaQ, 0);
  result(Signal::MMCControllerSystem::ISigmaD, 0) = state(ISigmaD, 0);
  result(Signal::MMCControllerSystem::ISigmaQ, 0) = state(ISigmaQ, 0);
  result(Signal::MMCControllerSystem::ISigmaZ, 0) = state(ISigmaZ, 0);
  result(Signal::MMCControllerSystem::StoredEnergy, 0) =
      calculateStoredEnergy(state);
  result(Signal::MMCControllerSystem::FeedforwardFilteredDaxisVoltage, 0) =
      feedforwardFilteredDaxisVoltage();
  result(Signal::MMCControllerSystem::HeldActiveCurrentReference, 0) =
      heldActiveCurrentReference();
  return result;
}

Matrix MMC_ModularStateSpaceModel::controllerOutput(const Matrix &state,
                                                    const Matrix &input) const {
  validateStateInput(state, input);
  Matrix dxController;
  Matrix output;
  mControllers.evaluate(
      state.block(controllerOffset(), 0, mControllers.stateSize(), 1),
      makeControllerInput(state, input), dxController, output);
  return output;
}

void MMC_ModularStateSpaceModel::evaluate(const Matrix &state,
                                          const Matrix &input,
                                          Matrix &stateDerivative,
                                          Matrix &output) const {
  evaluateStateDerivative(state, input, stateDerivative);
  evaluateOutput(state, input, output);
}

void MMC_ModularStateSpaceModel::evaluateStateDerivative(const Matrix &state,
                                                         const Matrix &input,
                                                         Matrix &f) const {
  validateStateInput(state, input);
  f = Matrix::Zero(stateSize(), 1);

  const Matrix controllerInput = makeControllerInput(state, input);
  Matrix dxController;
  Matrix control;
  mControllers.evaluate(
      state.block(controllerOffset(), 0, mControllers.stateSize(), 1),
      controllerInput, dxController, control);
  f.block(controllerOffset(), 0, mControllers.stateSize(), 1) = dxController;

  const Real iDeltaD = state(IDeltaD, 0);
  const Real iDeltaQ = state(IDeltaQ, 0);
  const Real iSigmaZ = state(ISigmaZ, 0);
  const Real iSigmaD = state(ISigmaD, 0);
  const Real iSigmaQ = state(ISigmaQ, 0);
  const Real vCDeltaD = state(VCDeltaD, 0);
  const Real vCDeltaQ = state(VCDeltaQ, 0);
  const Real vCDeltaZd = state(VCDeltaZd, 0);
  const Real vCDeltaZq = state(VCDeltaZq, 0);
  const Real vCSigmaD = state(VCSigmaD, 0);
  const Real vCSigmaQ = state(VCSigmaQ, 0);
  const Real vCSigmaZ = state(VCSigmaZ, 0);
  const Real vGd = controllerInput(Signal::MMCControllerSystem::VGridD, 0);
  const Real vGq = controllerInput(Signal::MMCControllerSystem::VGridQ, 0);
  const Real vDc = controllerInput(Signal::MMCControllerSystem::DcVoltage, 0);

  const Real mDeltaD = control(Signal::MMCControllerSystem::MDeltaD, 0);
  const Real mDeltaQ = control(Signal::MMCControllerSystem::MDeltaQ, 0);
  const Real mDeltaZd = 0.0;
  const Real mDeltaZq = 0.0;
  const Real mSigmaD = control(Signal::MMCControllerSystem::MSigmaD, 0);
  const Real mSigmaQ = control(Signal::MMCControllerSystem::MSigmaQ, 0);
  const Real mSigmaZ = control(Signal::MMCControllerSystem::MSigmaZ, 0);

  const Real vMDeltaD =
      (mDeltaQ * vCSigmaQ) / 4.0 - (mDeltaD * vCSigmaZ) / 2.0 -
      (mDeltaD * vCSigmaD) / 4.0 - (mDeltaZd * vCSigmaD) / 4.0 +
      (mDeltaZq * vCSigmaQ) / 4.0 - (mSigmaD * vCDeltaD) / 4.0 -
      (mSigmaZ * vCDeltaD) / 2.0 + (mSigmaQ * vCDeltaQ) / 4.0 -
      (mSigmaD * vCDeltaZd) / 4.0 + (mSigmaQ * vCDeltaZq) / 4.0;

  const Real vMDeltaQ =
      (mDeltaD * vCSigmaQ) / 4.0 + (mDeltaQ * vCSigmaD) / 4.0 -
      (mDeltaQ * vCSigmaZ) / 2.0 - (mDeltaZd * vCSigmaQ) / 4.0 -
      (mDeltaZq * vCSigmaD) / 4.0 + (mSigmaD * vCDeltaQ) / 4.0 +
      (mSigmaQ * vCDeltaD) / 4.0 - (mSigmaZ * vCDeltaQ) / 2.0 -
      (mSigmaD * vCDeltaZq) / 4.0 - (mSigmaQ * vCDeltaZd) / 4.0;

  const Real vMSigmaD =
      (mDeltaD * vCDeltaD) / 4.0 - (mDeltaQ * vCDeltaQ) / 4.0 +
      (mDeltaD * vCDeltaZd) / 4.0 + (mDeltaZd * vCDeltaD) / 4.0 +
      (mDeltaQ * vCDeltaZq) / 4.0 + (mDeltaZq * vCDeltaQ) / 4.0 +
      (mSigmaD * vCSigmaZ) / 2.0 + (mSigmaZ * vCSigmaD) / 2.0;

  const Real vMSigmaQ =
      (mDeltaQ * vCDeltaZd) / 4.0 - (mDeltaQ * vCDeltaD) / 4.0 -
      (mDeltaD * vCDeltaZq) / 4.0 - (mDeltaD * vCDeltaQ) / 4.0 +
      (mDeltaZd * vCDeltaQ) / 4.0 - (mDeltaZq * vCDeltaD) / 4.0 +
      (mSigmaQ * vCSigmaZ) / 2.0 + (mSigmaZ * vCSigmaQ) / 2.0;

  const Real vMSigmaZ =
      (mDeltaD * vCDeltaD) / 4.0 + (mDeltaQ * vCDeltaQ) / 4.0 +
      (mDeltaZd * vCDeltaZd) / 4.0 + (mDeltaZq * vCDeltaZq) / 4.0 +
      (mSigmaD * vCSigmaD) / 4.0 + (mSigmaQ * vCSigmaQ) / 4.0 +
      (mSigmaZ * vCSigmaZ) / 2.0;

  const Real lEq = 0.5 * mArmInductance + mReactorInductance;
  const Real rEq = 0.5 * mArmResistance + mReactorResistance;
  f(IDeltaD, 0) =
      (vMDeltaD - vGd - rEq * iDeltaD + lEq * iDeltaQ * mOmegaN) / lEq;
  f(IDeltaQ, 0) =
      (vMDeltaQ - vGq - rEq * iDeltaQ - lEq * iDeltaD * mOmegaN) / lEq;
  f(ISigmaD, 0) = -(vMSigmaD + mArmResistance * iSigmaD +
                    2.0 * mArmInductance * iSigmaQ * mOmegaN) /
                  mArmInductance;
  f(ISigmaQ, 0) = -(vMSigmaQ + mArmResistance * iSigmaQ -
                    2.0 * mArmInductance * iSigmaD * mOmegaN) /
                  mArmInductance;
  f(ISigmaZ, 0) =
      -(vMSigmaZ - vDc / 2.0 + mArmResistance * iSigmaZ) / mArmInductance;

  const Real cArm = mSubmoduleCapacitance;
  const Real n = static_cast<Real>(mNumberOfSubmodules);
  f(VCSigmaD, 0) = n *
                   (iSigmaD * mSigmaZ + iSigmaZ * mSigmaD +
                    iDeltaD * (mDeltaD / 4.0 + mDeltaZd / 4.0) -
                    iDeltaQ * (mDeltaQ / 4.0 - mDeltaZq / 4.0) -
                    (4.0 * cArm * vCSigmaQ * mOmegaN) / n) /
                   (2.0 * cArm);
  f(VCSigmaQ, 0) =
      -n *
      (iDeltaQ * (mDeltaD / 4.0 - mDeltaZd / 4.0) - iSigmaZ * mSigmaQ -
       iSigmaQ * mSigmaZ + iDeltaD * (mDeltaQ / 4.0 + mDeltaZq / 4.0) -
       (4.0 * cArm * vCSigmaD * mOmegaN) / n) /
      (2.0 * cArm);
  f(VCSigmaZ, 0) =
      n *
      (iDeltaD * mDeltaD + iDeltaQ * mDeltaQ + 2.0 * iSigmaD * mSigmaD +
       2.0 * iSigmaQ * mSigmaQ + 4.0 * iSigmaZ * mSigmaZ) /
      (8.0 * cArm);
  f(VCDeltaD, 0) = n *
                   (iSigmaZ * mDeltaD - (iDeltaQ * mSigmaQ) / 4.0 +
                    iSigmaD * (mDeltaD / 2.0 + mDeltaZd / 2.0) -
                    iSigmaQ * (mDeltaQ / 2.0 + mDeltaZq / 2.0) +
                    iDeltaD * (mSigmaD / 4.0 + mSigmaZ / 2.0) +
                    (2.0 * cArm * vCDeltaQ * mOmegaN) / n) /
                   (2.0 * cArm);
  f(VCDeltaQ, 0) = -n *
                   ((iDeltaD * mSigmaQ) / 4.0 - iSigmaZ * mDeltaQ +
                    iSigmaQ * (mDeltaD / 2.0 - mDeltaZd / 2.0) +
                    iSigmaD * (mDeltaQ / 2.0 - mDeltaZq / 2.0) +
                    iDeltaQ * (mSigmaD / 4.0 - mSigmaZ / 2.0) +
                    (2.0 * cArm * vCDeltaD * mOmegaN) / n) /
                   (2.0 * cArm);
  f(VCDeltaZd, 0) =
      n *
          (iDeltaD * mSigmaD + 2.0 * iSigmaD * mDeltaD + iDeltaQ * mSigmaQ +
           2.0 * iSigmaQ * mDeltaQ + 4.0 * iSigmaZ * mDeltaZd) /
          (8.0 * cArm) +
      3.0 * vCDeltaZq * mOmegaN;
  f(VCDeltaZq, 0) =
      -3.0 * vCDeltaZd * mOmegaN +
      n *
          (iDeltaQ * mSigmaD - iDeltaD * mSigmaQ + 2.0 * iSigmaD * mDeltaQ -
           2.0 * iSigmaQ * mDeltaD + 4.0 * iSigmaZ * mDeltaZq) /
          (8.0 * cArm);
  f(GridAngle, 0) = mOmegaN;

  if (!f.allFinite())
    throw std::runtime_error("MMC state derivative contains NaN or Inf.");
}

void MMC_ModularStateSpaceModel::evaluateOutput(const Matrix &state,
                                                const Matrix &input,
                                                Matrix &output) const {
  validateStateInput(state, input);
  output = Matrix::Zero(OutputCount, 1);
  const Matrix iAbc =
      dqToAbc(state(IDeltaD, 0), state(IDeltaQ, 0), state(GridAngle, 0));
  output.block(0, 0, 3, 1) = -iAbc;
  const Real idc = 3.0 * state(ISigmaZ, 0);
  output(Idcp, 0) = idc;
  output(Idcn, 0) = -idc;
}

std::vector<UInt>
MMC_ModularStateSpaceModel::activeLinearizationStateIndices() const {
  std::vector<UInt> indices;
  indices.reserve(stateSize());
  for (UInt index = 0; index < PlantStateCount; ++index)
    indices.push_back(index);
  for (const UInt local : mControllers.activeStateIndices())
    indices.push_back(controllerOffset() + local);
  return indices;
}

std::vector<UInt> MMC_ModularStateSpaceModel::equilibriumStateIndices() const {
  std::vector<UInt> indices;
  indices.reserve(stateSize());

  // All electrical plant states participate; the absolute network angle does
  // not because its derivative is the nominal angular frequency by design.
  for (UInt index = 0; index < GridAngle; ++index)
    indices.push_back(index);

  for (const UInt local : mControllers.equilibriumStateIndices())
    indices.push_back(controllerOffset() + local);
  return indices;
}

void MMC_ModularStateSpaceModel::calculateFullNumericalJacobians(
    const Matrix &state, const Matrix &input, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D) const {
  validateStateInput(state, input);
  A = Matrix::Zero(stateSize(), stateSize());
  B = Matrix::Zero(stateSize(), InputCount);
  C = Matrix::Zero(OutputCount, stateSize());
  D = Matrix::Zero(OutputCount, InputCount);

  Matrix fPlus = Matrix::Zero(stateSize(), 1);
  Matrix fMinus = Matrix::Zero(stateSize(), 1);
  Matrix gPlus = Matrix::Zero(OutputCount, 1);
  Matrix gMinus = Matrix::Zero(OutputCount, 1);

  for (UInt column = 0; column < stateSize(); ++column) {
    const Real step =
        mJacobianAbsoluteStep +
        mJacobianRelativeStep * std::max(1.0, std::abs(state(column, 0)));
    Matrix xPlus = state;
    Matrix xMinus = state;
    xPlus(column, 0) += step;
    xMinus(column, 0) -= step;
    evaluateStateDerivative(xPlus, input, fPlus);
    evaluateStateDerivative(xMinus, input, fMinus);
    evaluateOutput(xPlus, input, gPlus);
    evaluateOutput(xMinus, input, gMinus);
    A.col(column) = (fPlus - fMinus) / (2.0 * step);
    C.col(column) = (gPlus - gMinus) / (2.0 * step);
  }

  for (UInt column = 0; column < InputCount; ++column) {
    const Real step =
        mJacobianAbsoluteStep +
        mJacobianRelativeStep * std::max(1.0, std::abs(input(column, 0)));
    Matrix uPlus = input;
    Matrix uMinus = input;
    uPlus(column, 0) += step;
    uMinus(column, 0) -= step;
    evaluateStateDerivative(state, uPlus, fPlus);
    evaluateStateDerivative(state, uMinus, fMinus);
    evaluateOutput(state, uPlus, gPlus);
    evaluateOutput(state, uMinus, gMinus);
    B.col(column) = (fPlus - fMinus) / (2.0 * step);
    D.col(column) = (gPlus - gMinus) / (2.0 * step);
  }
}

void MMC_ModularStateSpaceModel::calculateNumericalJacobians(
    const Matrix &state, const Matrix &input, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D) const {
  if (mStructuredLinearization)
    calculateStructuredNumericalJacobians(state, input, A, B, C, D);
  else
    calculateFullNumericalJacobians(state, input, A, B, C, D);
}

void MMC_ModularStateSpaceModel::calculateStructuredNumericalJacobians(
    const Matrix &state, const Matrix &input, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D) const {
  validateStateInput(state, input);
  A = Matrix::Zero(stateSize(), stateSize());
  B = Matrix::Zero(stateSize(), InputCount);
  C = Matrix::Zero(OutputCount, stateSize());
  D = Matrix::Zero(OutputCount, InputCount);

  // State Jacobian: skip controller states that are structurally disconnected
  // in the selected control/filter/delay configuration.
  Matrix fPlus = Matrix::Zero(stateSize(), 1);
  Matrix fMinus = Matrix::Zero(stateSize(), 1);
  for (const UInt column : activeLinearizationStateIndices()) {
    const Real step =
        mJacobianAbsoluteStep +
        mJacobianRelativeStep * std::max(1.0, std::abs(state(column, 0)));
    Matrix xPlus = state;
    Matrix xMinus = state;
    xPlus(column, 0) += step;
    xMinus(column, 0) -= step;
    evaluateStateDerivative(xPlus, input, fPlus);
    evaluateStateDerivative(xMinus, input, fMinus);
    A.col(column) = (fPlus - fMinus) / (2.0 * step);
  }

  // The nonlinear MMC equations depend on the three AC terminal voltages only
  // through the balanced dq projection. Differentiate in the d and q
  // directions, then map those derivatives back to abc. This replaces six
  // full state-derivative evaluations by four.
  const Real theta = state(GridAngle, 0);
  const Matrix vDq = abcToDq(input.block(0, 0, 3, 1), theta);
  Matrix derivativeByD = Matrix::Zero(stateSize(), 1);
  Matrix derivativeByQ = Matrix::Zero(stateSize(), 1);
  const Matrix dDirection = dqToAbc(1.0, 0.0, theta);
  const Matrix qDirection = dqToAbc(0.0, 1.0, theta);

  const auto differentiateAcDirection =
      [&](const Matrix &direction, Real operatingValue, Matrix &derivative) {
        const Real step =
            mJacobianAbsoluteStep +
            mJacobianRelativeStep * std::max(1.0, std::abs(operatingValue));
        Matrix uPlus = input;
        Matrix uMinus = input;
        uPlus.block(0, 0, 3, 1) += step * direction;
        uMinus.block(0, 0, 3, 1) -= step * direction;
        evaluateStateDerivative(state, uPlus, fPlus);
        evaluateStateDerivative(state, uMinus, fMinus);
        derivative = (fPlus - fMinus) / (2.0 * step);
      };
  differentiateAcDirection(dDirection, vDq(0, 0), derivativeByD);
  differentiateAcDirection(qDirection, vDq(1, 0), derivativeByQ);

  const Real phaseAngles[3] = {theta, theta - 2.0 * PI / 3.0,
                               theta + 2.0 * PI / 3.0};
  for (UInt phase = 0; phase < 3; ++phase) {
    const Real dCoefficient = (2.0 / 3.0) * std::cos(phaseAngles[phase]);
    const Real qCoefficient = -(2.0 / 3.0) * std::sin(phaseAngles[phase]);
    B.col(phase) = dCoefficient * derivativeByD + qCoefficient * derivativeByQ;
  }

  // DC common-mode voltage is structurally absent. Differentiate only with
  // respect to Vdc = Vdc+ - Vdc- and map it to the two terminal columns.
  const Real vdc = input(Vdcp, 0) - input(Vdcn, 0);
  const Real dcStep = mJacobianAbsoluteStep +
                      mJacobianRelativeStep * std::max(1.0, std::abs(vdc));
  Matrix uPlus = input;
  Matrix uMinus = input;
  uPlus(Vdcp, 0) += 0.5 * dcStep;
  uPlus(Vdcn, 0) -= 0.5 * dcStep;
  uMinus(Vdcp, 0) -= 0.5 * dcStep;
  uMinus(Vdcn, 0) += 0.5 * dcStep;
  evaluateStateDerivative(state, uPlus, fPlus);
  evaluateStateDerivative(state, uMinus, fMinus);
  const Matrix derivativeByVdc = (fPlus - fMinus) / (2.0 * dcStep);
  B.col(Vdcp) = derivativeByVdc;
  B.col(Vdcn) = -derivativeByVdc;

  // Converter terminal currents have a compact analytic output Jacobian.
  // There is no direct voltage feedthrough, so D is exactly zero.
  const Real id = state(IDeltaD, 0);
  const Real iq = state(IDeltaQ, 0);
  for (UInt phase = 0; phase < 3; ++phase) {
    const Real angle = phaseAngles[phase];
    C(phase, IDeltaD) = -std::cos(angle);
    C(phase, IDeltaQ) = std::sin(angle);
    C(phase, GridAngle) = id * std::sin(angle) + iq * std::cos(angle);
  }
  C(Idcp, ISigmaZ) = 3.0;
  C(Idcn, ISigmaZ) = -3.0;
}

void MMC_ModularStateSpaceModel::buildStateSpaceModel(
    const Matrix &state, const Matrix &input, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Matrix &E, Matrix &F) const {
  calculateNumericalJacobians(state, input, A, B, C, D);

  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  Matrix output = Matrix::Zero(OutputCount, 1);
  evaluateStateDerivative(state, input, stateDerivative);
  evaluateOutput(state, input, output);

  // Local affine offsets:
  //
  //   E = f(x0,u0) - A*x0 - B*u0
  //   F = g(x0,u0) - C*x0 - D*u0
  E = stateDerivative - A * state - B * input;
  F = output - C * state - D * input;
}

Signal::MMCLinearization
MMC_ModularStateSpaceModel::getStateSpaceModel(const Matrix &state,
                                               const Matrix &input) const {
  Signal::MMCLinearization model;
  buildStateSpaceModel(state, input, model.A, model.B, model.C, model.D,
                       model.E, model.F);
  return model;
}

Signal::MMCSparseLinearization
MMC_ModularStateSpaceModel::getSparseStateSpaceModel(
    const Matrix &state, const Matrix &input) const {
  const Signal::MMCLinearization dense = getStateSpaceModel(state, input);
  Signal::MMCSparseLinearization sparse;
  sparse.A = dense.A.sparseView(0.0, mSparseLinearizationTolerance);
  sparse.B = dense.B.sparseView(0.0, mSparseLinearizationTolerance);
  sparse.C = dense.C.sparseView(0.0, mSparseLinearizationTolerance);
  sparse.D = dense.D.sparseView(0.0, mSparseLinearizationTolerance);
  sparse.E = dense.E;
  sparse.F = dense.F;
  return sparse;
}

void MMC_ModularStateSpaceModel::linearize(const Matrix &state,
                                           const Matrix &input, Matrix &A,
                                           Matrix &B, Matrix &C, Matrix &D,
                                           Matrix &E, Matrix &F) const {
  buildStateSpaceModel(state, input, A, B, C, D, E, F);
}

Matrix
MMC_ModularStateSpaceModel::initializeAnalyticalState(const Matrix &input) {
  Matrix state = Matrix::Zero(stateSize(), 1);
  state(GridAngle, 0) = mInitialAngle;
  const Matrix vdq = abcToDq(input.block(0, 0, 3, 1), mInitialAngle);
  const Real vdc =
      regularizeSigned(input(Vdcp, 0) - input(Vdcn, 0), mMinimumDcVoltage);
  const Real voltageSquared = vdq.squaredNorm();

  Real pTarget = mInitialOperatingPointEnabled
                     ? mInitialActivePower
                     : mControllers.activePowerReference();
  Real qTarget = mInitialOperatingPointEnabled
                     ? mInitialReactivePower
                     : mControllers.reactivePowerReference();

  if (voltageSquared > 0.0) {
    state(IDeltaD, 0) = (2.0 / 3.0) *
                        (vdq(0, 0) * pTarget + vdq(1, 0) * qTarget) /
                        voltageSquared;
    state(IDeltaQ, 0) = (2.0 / 3.0) *
                        (vdq(1, 0) * pTarget - vdq(0, 0) * qTarget) /
                        voltageSquared;
  }
  state(ISigmaD, 0) = mControllers.sigmaDReference();
  state(ISigmaQ, 0) = mControllers.sigmaQReference();
  // Preserve the explicitly supplied DC-current operating point whenever the
  // total-energy controller is disabled. This is essential for a loss-aware
  // P2P hot start, where P/(3*Vdc) is only an approximation to iSigma_z.
  state(ISigmaZ, 0) = mControllers.sigmaZReference();
  if (mControllers.energyControllerEnabled() &&
      (mInitialOperatingPointEnabled ||
       mControllers.activeMode() != Signal::MMCActiveMode::OpenLoop))
    state(ISigmaZ, 0) = pTarget / (3.0 * vdc);

  state(VCSigmaZ, 0) = std::abs(vdc);

  if (mControllers.sampledFeedforwardEnabled()) {
    mSampledFeedforwardState = mSampledFeedforward.initialize(
        vdq(0, 0), mControllers.activePowerReference());
    mFeedforwardInitialized = true;
  } else {
    mSampledFeedforwardState = Matrix::Zero(3, 1);
    mFeedforwardInitialized = false;
  }

  const Matrix controllerInput = makeControllerInput(state, input);
  const Real lEq = 0.5 * mArmInductance + mReactorInductance;
  const Real rEq = 0.5 * mArmResistance + mReactorResistance;
  const Real vMDeltaD =
      vdq(0, 0) + rEq * state(IDeltaD, 0) - lEq * mOmegaN * state(IDeltaQ, 0);
  const Real vMDeltaQ =
      vdq(1, 0) + rEq * state(IDeltaQ, 0) + lEq * mOmegaN * state(IDeltaD, 0);
  const Real vMSigmaD = -mArmResistance * state(ISigmaD, 0) -
                        2.0 * mArmInductance * mOmegaN * state(ISigmaQ, 0);
  const Real vMSigmaQ = -mArmResistance * state(ISigmaQ, 0) +
                        2.0 * mArmInductance * mOmegaN * state(ISigmaD, 0);
  const Real vMSigmaZ = vdc / 2.0 - mArmResistance * state(ISigmaZ, 0);
  Matrix modulation(5, 1);
  modulation << -2.0 * vMDeltaD / vdc, -2.0 * vMDeltaQ / vdc,
      2.0 * vMSigmaD / vdc, 2.0 * vMSigmaQ / vdc, 2.0 * vMSigmaZ / vdc;
  state.block(controllerOffset(), 0, mControllers.stateSize(), 1) =
      mControllers.initializeState(controllerInput, modulation);
  return state;
}

Bool MMC_ModularStateSpaceModel::solveOperatingPoint(
    Matrix &state, const Matrix &input, Real &normalizedResidual,
    Real &absoluteResidual) const {
  const std::vector<UInt> indices = equilibriumStateIndices();
  const UInt dimension = static_cast<UInt>(indices.size());
  if (dimension == 0) {
    normalizedResidual = 0.0;
    absoluteResidual = 0.0;
    return true;
  }

  const UInt activeIntegratorIndex =
      controllerOffset() + mControllers.activeIntegratorStateIndex();

  auto calculateResidual = [&](const Matrix &candidate, Matrix &rawResidual,
                               Matrix &scaledResidual, Matrix &rowScale) {
    Matrix derivative;
    evaluateStateDerivative(candidate, input, derivative);

    rawResidual.resize(dimension, 1);
    scaledResidual.resize(dimension, 1);
    rowScale.resize(dimension, 1);

    for (UInt row = 0; row < dimension; ++row) {
      const UInt stateIndex = indices[row];

      // In standalone V-type DC-voltage mode the DC terminal voltage is an
      // imposed input. The active-integrator equilibrium equation is therefore
      // redundant and would drive a minimum-norm solver toward an unloaded
      // solution. Replace it by the requested loaded active-power constraint.
      if (mInitialOperatingPointEnabled &&
          mControllers.activeMode() == Signal::MMCActiveMode::DcVoltage &&
          stateIndex == activeIntegratorIndex) {
        rawResidual(row, 0) =
            activePower(candidate, input) - mInitialActivePower;
        rowScale(row, 0) = std::max(1.0, std::abs(mInitialActivePower));
      } else {
        rawResidual(row, 0) = derivative(stateIndex, 0);
        rowScale(row, 0) = std::max(
            1.0, mOmegaN * std::max(1.0, std::abs(candidate(stateIndex, 0))));
      }
      scaledResidual(row, 0) = rawResidual(row, 0) / rowScale(row, 0);
    }
  };

  Matrix rawResidual;
  Matrix scaledResidual;
  Matrix rowScale;
  calculateResidual(state, rawResidual, scaledResidual, rowScale);
  normalizedResidual = scaledResidual.norm();
  absoluteResidual = rawResidual.norm();

  for (UInt iteration = 0; iteration < mOperatingPointMaximumIterations;
       ++iteration) {
    if (normalizedResidual <= mOperatingPointNormalizedTolerance)
      return true;

    Matrix jacobian = Matrix::Zero(dimension, dimension);
    Matrix columnScale = Matrix::Ones(dimension, 1);

    for (UInt column = 0; column < dimension; ++column) {
      const UInt stateIndex = indices[column];
      const Real step =
          mJacobianAbsoluteStep +
          mJacobianRelativeStep * std::max(1.0, std::abs(state(stateIndex, 0)));

      Matrix xPlus = state;
      Matrix xMinus = state;
      xPlus(stateIndex, 0) += step;
      xMinus(stateIndex, 0) -= step;

      Matrix rawPlus;
      Matrix scaledPlus;
      Matrix scalePlus;
      calculateResidual(xPlus, rawPlus, scaledPlus, scalePlus);

      Matrix rawMinus;
      Matrix scaledMinus;
      Matrix scaleMinus;
      calculateResidual(xMinus, rawMinus, scaledMinus, scaleMinus);

      columnScale(column, 0) = std::max(1.0, std::abs(state(stateIndex, 0)));

      for (UInt row = 0; row < dimension; ++row) {
        // Keep current-iterate row scaling fixed while differentiating.
        jacobian(row, column) = (rawPlus(row, 0) - rawMinus(row, 0)) /
                                (2.0 * step * rowScale(row, 0));
      }
    }

    Matrix scaledJacobian = jacobian;
    for (UInt column = 0; column < dimension; ++column)
      scaledJacobian.col(column) *= columnScale(column, 0);

    const Matrix scaledCorrection =
        scaledJacobian.colPivHouseholderQr().solve(-scaledResidual);
    if (!scaledCorrection.allFinite())
      return false;

    Matrix correction = Matrix::Zero(dimension, 1);
    for (UInt column = 0; column < dimension; ++column)
      correction(column, 0) =
          columnScale(column, 0) * scaledCorrection(column, 0);

    Bool accepted = false;
    Real lineSearchFactor = 1.0;
    for (UInt lineSearch = 0; lineSearch < 12; ++lineSearch) {
      Matrix candidate = state;
      for (UInt column = 0; column < dimension; ++column)
        candidate(indices[column], 0) +=
            lineSearchFactor * correction(column, 0);

      Matrix candidateRaw;
      Matrix candidateScaled;
      Matrix candidateScale;
      calculateResidual(candidate, candidateRaw, candidateScaled,
                        candidateScale);

      if (candidateScaled.norm() < normalizedResidual) {
        state = candidate;
        rawResidual = candidateRaw;
        scaledResidual = candidateScaled;
        rowScale = candidateScale;
        normalizedResidual = scaledResidual.norm();
        absoluteResidual = rawResidual.norm();
        accepted = true;
        break;
      }
      lineSearchFactor *= 0.5;
    }

    if (!accepted)
      break;
  }

  return normalizedResidual <= mOperatingPointNormalizedTolerance;
}

Matrix MMC_ModularStateSpaceModel::initializeState(const Matrix &input) {
  validateConfigured();
  if (input.rows() != InputCount || input.cols() != 1 || !input.allFinite())
    throw std::invalid_argument("MMC initialization input must be finite 5x1.");
  Matrix state = initializeAnalyticalState(input);
  if (mOperatingPointInitializationEnabled) {
    Real normalized = 0.0;
    Real absolute = 0.0;
    solveOperatingPoint(state, input, normalized, absolute);
  }
  mFeedforwardStepCounter = 0;
  return state;
}

Bool MMC_ModularStateSpaceModel::advanceSampledControllers(const Matrix &state,
                                                           const Matrix &input,
                                                           Real timeStep) {
  validateStateInput(state, input);
  if (!std::isfinite(timeStep) || timeStep <= 0.0)
    throw std::invalid_argument("MMC model time step must be positive.");
  if (!mControllers.sampledFeedforwardEnabled() ||
      mControllers.controlSource() != ControlSource::InternalControllers)
    return false;
  if (!mFeedforwardInitialized)
    throw std::logic_error(
        "MMC sampled active-power feedforward was not initialized.");

  const UInt stride = std::max<UInt>(
      1, static_cast<UInt>(
             std::llround(mControllers.sampledFeedforwardPeriod() / timeStep)));

  Bool sampled = false;
  if (mFeedforwardStepCounter == 0) {
    // Match the validated standalone implementation: the feedforward voltage
    // is measured in the nominal network dq frame, not in the PLL-deviation
    // frame.
    const Matrix vGridDq =
        abcToDq(input.block(0, 0, 3, 1), state(GridAngle, 0));
    const Real effectiveSampleTime = static_cast<Real>(stride) * timeStep;
    mSampledFeedforwardState = mSampledFeedforward.sample(
        mSampledFeedforwardState, vGridDq(0, 0),
        mControllers.activePowerReference(), effectiveSampleTime);
    sampled = true;
  }

  mFeedforwardStepCounter = (mFeedforwardStepCounter + 1) % stride;
  return sampled;
}

Real MMC_ModularStateSpaceModel::feedforwardFilteredDaxisVoltage() const {
  return mFeedforwardInitialized ? mSampledFeedforwardState(0, 0) : 0.0;
}

Real MMC_ModularStateSpaceModel::heldActiveCurrentReference() const {
  return mFeedforwardInitialized ? mSampledFeedforwardState(2, 0) : 0.0;
}
