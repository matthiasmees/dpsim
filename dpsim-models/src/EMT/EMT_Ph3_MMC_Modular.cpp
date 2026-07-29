// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/EMT/EMT_Ph3_MMC_Modular.h>

#include <array>
#include <cmath>
#include <stdexcept>

#include <dpsim-models/MathUtils.h>

using namespace CPS;

namespace {
Matrix abcToDqForLogging(const Matrix &abc, Real theta) {
  if (abc.rows() != 3 || abc.cols() != 1)
    throw std::invalid_argument("MMC logging Park transform expects 3x1.");
  Matrix dq(2, 1);
  dq(0, 0) = (2.0 / 3.0) * (std::cos(theta) * abc(0, 0) +
                            std::cos(theta - 2.0 * PI / 3.0) * abc(1, 0) +
                            std::cos(theta + 2.0 * PI / 3.0) * abc(2, 0));
  dq(1, 0) = -(2.0 / 3.0) * (std::sin(theta) * abc(0, 0) +
                             std::sin(theta - 2.0 * PI / 3.0) * abc(1, 0) +
                             std::sin(theta + 2.0 * PI / 3.0) * abc(2, 0));
  return dq;
}
} // namespace

EMT::Ph3::MMC_Modular::MMC_Modular(String uid, String name,
                                   Logger::Level logLevel)
    : VTypeVariableSSNComp(uid, name, mInputSize, mOutputSize, logLevel),
      mModelConfigured(false), mInitializationInProgress(false),
      mLinearizationUpdateInterval(1), mStepsSinceLinearization(0),
      mLinearizationTimeStep(0.0), mTimeSinceLinearization(0.0),
      mLinearizationDirty(true), mLastLinearizedConfigurationRevision(0),
      mStampZeroTolerance(1e-14), mRuntimeDiagnostics(false),
      mLinearizationRebuildCount(0), mLinearizationReuseCount(0),
      mAcTerminalVoltage(mAttributes->create<Matrix>("ac_terminal_voltage",
                                                     Matrix::Zero(3, 1))),
      mAcTerminalCurrent(mAttributes->create<Matrix>("ac_terminal_current",
                                                     Matrix::Zero(3, 1))),
      mDcPositiveVoltage(mAttributes->create<Real>("vdcp", 0.0)),
      mDcNegativeVoltage(mAttributes->create<Real>("vdcn", 0.0)),
      mDcVoltage(mAttributes->create<Real>("vdc", 0.0)),
      mDcCurrent(mAttributes->create<Real>("idc", 0.0)),
      mActivePower(mAttributes->create<Real>("p_ac", 0.0)),
      mReactivePower(mAttributes->create<Real>("q_ac", 0.0)),
      mStoredEnergy(mAttributes->create<Real>("stored_energy", 0.0)),
      mFilteredActivePower(mAttributes->create<Real>("p_filtered", 0.0)),
      mFilteredReactivePower(mAttributes->create<Real>("q_filtered", 0.0)),
      mFilteredDcVoltage(mAttributes->create<Real>("vdc_filtered", 0.0)),
      mFilteredDaxisVoltage(
          mAttributes->create<Real>("v_d_feedforward_filtered", 0.0)),
      mHeldActiveCurrentReference(
          mAttributes->create<Real>("i_delta_d_feedforward_held", 0.0)),
      mPllFrequency(mAttributes->create<Real>("pll_frequency", 0.0)),
      mPllAngleDeviation(mAttributes->create<Real>("pll_angle_deviation", 0.0)),
      mGridVoltageD(mAttributes->create<Real>("v_grid_d", 0.0)),
      mGridVoltageQ(mAttributes->create<Real>("v_grid_q", 0.0)),
      mControlVoltageD(mAttributes->create<Real>("v_control_d", 0.0)),
      mControlVoltageQ(mAttributes->create<Real>("v_control_q", 0.0)),
      mDeltaCurrentD(mAttributes->create<Real>("i_delta_d", 0.0)),
      mDeltaCurrentQ(mAttributes->create<Real>("i_delta_q", 0.0)),
      mSigmaCurrentZ(mAttributes->create<Real>("i_sigma_z", 0.0)),
      mDeltaCurrentReferenceD(mAttributes->create<Real>("i_delta_d_ref", 0.0)),
      mDeltaCurrentReferenceQ(mAttributes->create<Real>("i_delta_q_ref", 0.0)),
      mSigmaCurrentReferenceZ(mAttributes->create<Real>("i_sigma_z_ref", 0.0)),
      mDcPower(mAttributes->create<Real>("p_dc", 0.0)),
      mPowerBalanceError(mAttributes->create<Real>("power_balance_error", 0.0)),
      mAppliedModulation(mAttributes->create<Matrix>("applied_modulation",
                                                     Matrix::Zero(5, 1))),
      mDifferentialVoltageCommand(mAttributes->create<Matrix>(
          "applied_differential_voltage", Matrix::Zero(2, 1))),
      mInternalDifferentialVoltage(mAttributes->create<Matrix>(
          "internal_differential_voltage", Matrix::Zero(2, 1))),
      mCommonModeVoltageCommand(mAttributes->create<Matrix>(
          "applied_common_mode_voltage", Matrix::Zero(3, 1))),
      mExternalDifferentialVoltage(
          mAttributes->createDynamic<Matrix>("external_differential_voltage")),
      mExternalCommonModeVoltage(
          mAttributes->createDynamic<Matrix>("external_common_mode_voltage")),
      mExternalCommandActive(
          mAttributes->create<Bool>("external_command_active", false)) {
  **mExternalDifferentialVoltage = Matrix::Zero(2, 1);
  **mExternalCommonModeVoltage = Matrix::Zero(3, 1);
  mPhaseType = PhaseType::ABC;
  setTerminalNumber(3);
}

void EMT::Ph3::MMC_Modular::setParameters(
    Real nominalFrequency, Real nominalAcVoltage, Real nominalDcVoltage,
    Real armInductance, Real armResistance, Real submoduleCapacitance,
    UInt numberOfSubmodules, Real reactorInductance, Real reactorResistance) {
  mModel.setParameters(nominalFrequency, nominalAcVoltage, nominalDcVoltage,
                       armInductance, armResistance, submoduleCapacitance,
                       numberOfSubmodules, reactorInductance,
                       reactorResistance);
  markLinearizationDirty();
  VTypeVariableSSNComp::setParameters(
      Matrix::Zero(mModel.stateSize(), mModel.stateSize()),
      Matrix::Zero(mModel.stateSize(), mInputSize),
      Matrix::Zero(mOutputSize, mModel.stateSize()),
      Matrix::Zero(mOutputSize, mInputSize),
      Matrix::Zero(mModel.stateSize(), 1), Matrix::Zero(mOutputSize, 1));
  mModelConfigured = true;
}

void EMT::Ph3::MMC_Modular::setInitialAngle(Real angle) {
  mModel.setInitialAngle(angle);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setInitialOperatingPoint(Real p, Real q) {
  mModel.setInitialOperatingPoint(p, q);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setPLL(Real kp, Real ki, Bool enabled) {
  mModel.setPLL(kp, ki, enabled);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setOutputCurrentController(Real kp, Real ki) {
  mModel.setOutputCurrentController(kp, ki);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setCirculatingCurrentController(Real kp, Real ki) {
  mModel.setCirculatingCurrentController(kp, ki);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setZeroSequenceCurrentController(Real kp, Real ki) {
  mModel.setZeroSequenceCurrentController(kp, ki);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setEnergyController(Real kp, Real ki,
                                                Bool enabled) {
  mModel.setEnergyController(kp, ki, enabled);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setActivePowerControl(Real reference, Real kp,
                                                  Real ki) {
  mModel.setActivePowerControl(reference, kp, ki);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setActivePowerFeedforwardControl(
    Real reference, Real cutoffFrequency, Real sampleTime,
    Real minimumDaxisVoltage) {
  mModel.setActivePowerFeedforwardControl(reference, cutoffFrequency,
                                          sampleTime, minimumDaxisVoltage);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setActivePowerFeedforwardReference(Real reference) {
  // The sampled controller applies the new reference at its next execution
  // instant. advanceSampledControllers() then forces the required affine-model
  // rebuild; no premature dirty flag is set here.
  mModel.setActivePowerReference(reference);
}
void EMT::Ph3::MMC_Modular::setDcVoltageControl(Real reference, Real kp,
                                                Real ki) {
  mModel.setDcVoltageControl(reference, kp, ki);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setDcDroopControl(Real p, Real vdc, Real droop) {
  mModel.setDcDroopControl(p, vdc, droop);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setActiveControlOpenLoop(Real reference) {
  mModel.setActiveControlOpenLoop(reference);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setReactivePowerControl(Real reference, Real kp,
                                                    Real ki) {
  mModel.setReactivePowerControl(reference, kp, ki);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setAcVoltageControl(Real reference, Real kp,
                                                Real ki) {
  mModel.setAcVoltageControl(reference, kp, ki);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setReactiveControlOpenLoop(Real reference) {
  mModel.setReactiveControlOpenLoop(reference);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setCirculatingCurrentReferences(Real dReference,
                                                            Real qReference,
                                                            Real zReference) {
  mModel.setCirculatingCurrentReferences(dReference, qReference, zReference);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setMeasurementFilters(Real vacDq, Real p, Real q,
                                                  Real vdc, Real vacMag) {
  mModel.setMeasurementFilters(vacDq, p, q, vdc, vacMag);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setModulationDelay(Real delay, UInt padeOrder) {
  mModel.setModulationDelay(delay, padeOrder);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setLimits(Real ac, Real circulating,
                                      Real modulation) {
  mModel.setLimits(ac, circulating, modulation);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setControlSource(ControlSource source) {
  mModel.setControlSource(source);
  markLinearizationDirty();
  **mExternalCommandActive = source != ControlSource::InternalControllers;
}
void EMT::Ph3::MMC_Modular::setExternalDifferentialVoltageCommand(Real dVolts,
                                                                  Real qVolts) {
  mModel.setExternalDifferentialVoltageCommand(dVolts, qVolts);
  markLinearizationDirty();
  (**mExternalDifferentialVoltage)(0, 0) = dVolts;
  (**mExternalDifferentialVoltage)(1, 0) = qVolts;
}
void EMT::Ph3::MMC_Modular::setExternalCommonModeVoltageCommand(Real dVolts,
                                                                Real qVolts,
                                                                Real zVolts) {
  mModel.setExternalCommonModeVoltageCommand(dVolts, qVolts, zVolts);
  markLinearizationDirty();
  (**mExternalCommonModeVoltage)(0, 0) = dVolts;
  (**mExternalCommonModeVoltage)(1, 0) = qVolts;
  (**mExternalCommonModeVoltage)(2, 0) = zVolts;
}
void EMT::Ph3::MMC_Modular::setNumericalLinearizationParameters(
    Real relativeStep, Real absoluteStep) {
  mModel.setNumericalLinearizationParameters(relativeStep, absoluteStep);
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setSparseLinearizationTolerance(Real tolerance) {
  mModel.setSparseLinearizationTolerance(tolerance);
  if (!std::isfinite(tolerance) || tolerance < 0.0)
    throw std::invalid_argument(
        "MMC stamp-zero tolerance must be finite and non-negative.");
  mStampZeroTolerance = tolerance;
  markLinearizationDirty();
}

void EMT::Ph3::MMC_Modular::setStructuredLinearization(Bool enabled) {
  mModel.setStructuredLinearization(enabled);
  markLinearizationDirty();
}

void EMT::Ph3::MMC_Modular::setRuntimeDiagnostics(Bool enabled) {
  mRuntimeDiagnostics = enabled;
}

void EMT::Ph3::MMC_Modular::setLinearizationUpdateInterval(UInt intervalSteps) {
  if (intervalSteps == 0)
    throw std::invalid_argument(
        "MMC linearization update interval must be at least one step.");
  mLinearizationUpdateInterval = intervalSteps;
  mLinearizationTimeStep = 0.0;
  mStepsSinceLinearization = 0;
  mTimeSinceLinearization = 0.0;
  markLinearizationDirty();
}

void EMT::Ph3::MMC_Modular::setLinearizationTimeStep(
    Real linearizationTimeStep) {
  if (!std::isfinite(linearizationTimeStep) || linearizationTimeStep <= 0.0)
    throw std::invalid_argument(
        "MMC linearization time step must be finite and positive.");
  mLinearizationTimeStep = linearizationTimeStep;
  mStepsSinceLinearization = 0;
  mTimeSinceLinearization = 0.0;
  markLinearizationDirty();
}
void EMT::Ph3::MMC_Modular::setOperatingPointInitialization(Bool enabled,
                                                            UInt iterations,
                                                            Real tolerance) {
  mModel.setOperatingPointInitialization(enabled, iterations, tolerance);
  markLinearizationDirty();
}

EMT::Ph3::MMC_ModularStateSpaceModel &EMT::Ph3::MMC_Modular::model() {
  return mModel;
}
const EMT::Ph3::MMC_ModularStateSpaceModel &
EMT::Ph3::MMC_Modular::model() const {
  return mModel;
}

void EMT::Ph3::MMC_Modular::markLinearizationDirty() {
  mLinearizationDirty = true;
}

Bool EMT::Ph3::MMC_Modular::shouldUpdateLinearization() {
  if (mModel.configurationRevision() != mLastLinearizedConfigurationRevision)
    mLinearizationDirty = true;

  if (mInitializationInProgress || mLinearizationDirty) {
    mStepsSinceLinearization = 0;
    mTimeSinceLinearization = 0.0;
    mLinearizationDirty = false;
    return true;
  }

  if (mLinearizationTimeStep > 0.0) {
    mTimeSinceLinearization += mTimeStep;
    if (mTimeSinceLinearization + 1e-15 >= mLinearizationTimeStep) {
      mTimeSinceLinearization = 0.0;
      mStepsSinceLinearization = 0;
      return true;
    }
    return false;
  }

  ++mStepsSinceLinearization;
  if (mStepsSinceLinearization >= mLinearizationUpdateInterval) {
    mStepsSinceLinearization = 0;
    return true;
  }
  return false;
}

void EMT::Ph3::MMC_Modular::validateConfigured() const {
  if (!mModelConfigured || !mParametersSet)
    throw std::logic_error("Configure MMC_Modular before use.");
}
void EMT::Ph3::MMC_Modular::validateTerminalArrangement() const {
  const auto ac = const_cast<MMC_Modular *>(this)->node(0);
  if (ac->isGround() || ac->phaseType() != PhaseType::ABC)
    throw std::invalid_argument(
        "MMC terminal 0 requires a non-grounded ABC node.");
  for (UInt terminal = 1; terminal < 3; ++terminal) {
    const auto dc = const_cast<MMC_Modular *>(this)->node(terminal);
    if (!dc->isGround() && dc->phaseType() != PhaseType::DC)
      throw std::invalid_argument(
          "MMC DC terminals require DC nodes or ground.");
  }
}

MatrixComp EMT::Ph3::MMC_Modular::buildInitialInputFromNodes(Real) {
  validateTerminalArrangement();
  MatrixComp u = MatrixComp::Zero(mInputSize, 1);
  const Complex va = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  u(0, 0) = va;
  u(1, 0) = va * SHIFT_TO_PHASE_B;
  u(2, 0) = va * SHIFT_TO_PHASE_C;
  u(3, 0) = initialSingleVoltage(1);
  u(4, 0) = initialSingleVoltage(2);
  return u;
}

void EMT::Ph3::MMC_Modular::initializeFromNodesAndTerminals(Real frequency) {
  validateConfigured();
  validateTerminalArrangement();
  Matrix u0 = buildInitialInputFromNodes(frequency).real();
  if (std::abs(u0(3, 0) - u0(4, 0)) < 1.0) {
    const Real nominalVdc = mModel.nominalDcVoltage();
    if (!terminalNotGrounded(2)) {
      u0(3, 0) = nominalVdc;
      u0(4, 0) = 0.0;
    } else if (!terminalNotGrounded(1)) {
      u0(3, 0) = 0.0;
      u0(4, 0) = -nominalVdc;
    } else {
      u0(3, 0) = 0.5 * nominalVdc;
      u0(4, 0) = -0.5 * nominalVdc;
    }
  }
  **mIntfVoltage = u0;
  **mX = mModel.initializeState(u0);
  mInitializationInProgress = true;
  updateStateSpaceModel();
  mInitializationInProgress = false;
  mYHist = calculateHistoryVector();
  **mIntfCurrent = mW * u0 + mYHist;
  updateLogAttributes(u0);
}

Bool EMT::Ph3::MMC_Modular::updateComponentParameters() {
  const auto source = mModel.controlSource();
  if (source == ControlSource::ExternalDifferentialVoltage ||
      source == ControlSource::ExternalFullConverterVoltage) {
    if ((**mExternalDifferentialVoltage).rows() != 2 ||
        (**mExternalDifferentialVoltage).cols() != 1 ||
        !(**mExternalDifferentialVoltage).allFinite())
      throw std::runtime_error(
          "External MMC differential-voltage command must be finite 2x1.");
    mModel.setExternalDifferentialVoltageCommand(
        (**mExternalDifferentialVoltage)(0, 0),
        (**mExternalDifferentialVoltage)(1, 0));
  }
  if (source == ControlSource::ExternalFullConverterVoltage) {
    if ((**mExternalCommonModeVoltage).rows() != 3 ||
        (**mExternalCommonModeVoltage).cols() != 1 ||
        !(**mExternalCommonModeVoltage).allFinite())
      throw std::runtime_error(
          "External MMC common-mode command must be finite 3x1.");
    mModel.setExternalCommonModeVoltageCommand(
        (**mExternalCommonModeVoltage)(0, 0),
        (**mExternalCommonModeVoltage)(1, 0),
        (**mExternalCommonModeVoltage)(2, 0));
  }

  // The sampled feedforward controller is a separate hybrid state-space block.
  // Its state is intentionally not part of mX. A sample changes the held
  // current reference and therefore forces one immediate affine-model rebuild.
  if (!mInitializationInProgress && std::isfinite(mTimeStep) &&
      mTimeStep > 0.0) {
    if (mModel.advanceSampledControllers(**mX, **mIntfVoltage, mTimeStep))
      markLinearizationDirty();
  }

  if (!shouldUpdateLinearization()) {
    ++mLinearizationReuseCount;
    return false;
  }

  Matrix E;
  Matrix F;
  mModel.buildStateSpaceModel(**mX, **mIntfVoltage, mA, mB, mC, mD, E, F);
  setStateOffset(E);
  setOutputOffset(F);
  mLastLinearizedConfigurationRevision = mModel.configurationRevision();
  ++mLinearizationRebuildCount;

  if (mRuntimeDiagnostics && (mLinearizationRebuildCount <= 5 ||
                              mLinearizationRebuildCount % 1000 == 0)) {
    SPDLOG_LOGGER_DEBUG(
        mSLog,
        "MMC linearization rebuild={}, reuse={}, state_norm={}, "
        "input_norm={}",
        mLinearizationRebuildCount, mLinearizationReuseCount, (**mX).norm(),
        (**mIntfVoltage).norm());
  }

  return true;
}

void EMT::Ph3::MMC_Modular::mnaCompUpdateVoltage(const Matrix &leftVector) {
  Matrix u = Matrix::Zero(mInputSize, 1);
  if (terminalNotGrounded(0))
    for (UInt phase = 0; phase < 3; ++phase)
      u(phase, 0) =
          Math::realFromVectorElement(leftVector, matrixNodeIndex(0, phase));
  if (terminalNotGrounded(1))
    u(3, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
  if (terminalNotGrounded(2))
    u(4, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(2, 0));
  **mIntfVoltage = u;
}

void EMT::Ph3::MMC_Modular::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  std::array<Int, mInputSize> indices = {-1, -1, -1, -1, -1};
  if (terminalNotGrounded(0))
    for (UInt phase = 0; phase < 3; ++phase)
      indices[phase] = static_cast<Int>(matrixNodeIndex(0, phase));
  if (terminalNotGrounded(1))
    indices[3] = static_cast<Int>(matrixNodeIndex(1, 0));
  if (terminalNotGrounded(2))
    indices[4] = static_cast<Int>(matrixNodeIndex(2, 0));
  for (UInt row = 0; row < mOutputSize; ++row)
    for (UInt column = 0; column < mInputSize; ++column)
      if (indices[row] >= 0 && indices[column] >= 0 &&
          std::abs(mW(row, column)) > mStampZeroTolerance)
        Math::addToMatrixElement(systemMatrix, static_cast<UInt>(indices[row]),
                                 static_cast<UInt>(indices[column]),
                                 mW(row, column));
}

void EMT::Ph3::MMC_Modular::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  std::array<Int, mOutputSize> indices = {-1, -1, -1, -1, -1};
  if (terminalNotGrounded(0))
    for (UInt phase = 0; phase < 3; ++phase)
      indices[phase] = static_cast<Int>(matrixNodeIndex(0, phase));
  if (terminalNotGrounded(1))
    indices[3] = static_cast<Int>(matrixNodeIndex(1, 0));
  if (terminalNotGrounded(2))
    indices[4] = static_cast<Int>(matrixNodeIndex(2, 0));
  for (UInt row = 0; row < mOutputSize; ++row)
    if (indices[row] >= 0)
      Math::setVectorElement(rightVector, static_cast<UInt>(indices[row]),
                             -mYHist(row, 0));
}

void EMT::Ph3::MMC_Modular::updateLogAttributes(const Matrix &u) const {
  **mAcTerminalVoltage = u.block(0, 0, 3, 1);
  **mAcTerminalCurrent = (**mIntfCurrent).block(0, 0, 3, 1);
  **mDcPositiveVoltage = u(3, 0);
  **mDcNegativeVoltage = u(4, 0);
  **mDcVoltage = mModel.dcVoltage(u);
  **mDcCurrent = mModel.dcCurrent(**mX);
  **mActivePower = mModel.activePower(**mX, u);
  **mReactivePower = mModel.reactivePower(**mX, u);
  **mStoredEnergy = mModel.calculateStoredEnergy(**mX);
  **mDcPower = (**mDcVoltage) * (**mDcCurrent);
  **mPowerBalanceError = mModel.powerBalanceError(**mX, u);

  const Matrix control = mModel.controllerOutput(**mX, u);
  const Real gridAngle = (**mX)(MMC_ModularStateSpaceModel::GridAngle, 0);
  const Matrix vGridDq = abcToDqForLogging(u.block(0, 0, 3, 1), gridAngle);

  **mGridVoltageD = vGridDq(0, 0);
  **mGridVoltageQ = vGridDq(1, 0);
  **mControlVoltageD = control(Signal::MMCControllerSystem::VControlD, 0);
  **mControlVoltageQ = control(Signal::MMCControllerSystem::VControlQ, 0);
  **mDeltaCurrentD = (**mX)(MMC_ModularStateSpaceModel::IDeltaD, 0);
  **mDeltaCurrentQ = (**mX)(MMC_ModularStateSpaceModel::IDeltaQ, 0);
  **mSigmaCurrentZ = (**mX)(MMC_ModularStateSpaceModel::ISigmaZ, 0);
  **mDeltaCurrentReferenceD =
      control(Signal::MMCControllerSystem::IDeltaDReference, 0);
  **mDeltaCurrentReferenceQ =
      control(Signal::MMCControllerSystem::IDeltaQReference, 0);
  **mSigmaCurrentReferenceZ =
      control(Signal::MMCControllerSystem::ISigmaZReference, 0);
  **mPllAngleDeviation = std::remainder(
      control(Signal::MMCControllerSystem::PllAngle, 0), 2.0 * PI);
  **mPllFrequency =
      mModel.nominalFrequency() +
      control(Signal::MMCControllerSystem::DeltaOmega, 0) / (2.0 * PI);
  **mFilteredActivePower =
      control(Signal::MMCControllerSystem::FilteredActivePower, 0);
  **mFilteredReactivePower =
      control(Signal::MMCControllerSystem::FilteredReactivePower, 0);
  **mFilteredDcVoltage =
      control(Signal::MMCControllerSystem::FilteredDcVoltage, 0);
  **mFilteredDaxisVoltage = control(
      Signal::MMCControllerSystem::FeedforwardFilteredDaxisVoltageOutput, 0);
  **mHeldActiveCurrentReference = mModel.heldActiveCurrentReference();

  **mAppliedModulation = control.block(0, 0, 5, 1);
  (**mDifferentialVoltageCommand)(0, 0) =
      control(Signal::MMCControllerSystem::VMDeltaDCommand, 0);
  (**mDifferentialVoltageCommand)(1, 0) =
      control(Signal::MMCControllerSystem::VMDeltaQCommand, 0);
  (**mInternalDifferentialVoltage)(0, 0) =
      control(Signal::MMCControllerSystem::InternalVMDeltaDCommand, 0);
  (**mInternalDifferentialVoltage)(1, 0) =
      control(Signal::MMCControllerSystem::InternalVMDeltaQCommand, 0);
  (**mCommonModeVoltageCommand)(0, 0) =
      control(Signal::MMCControllerSystem::VMSigmaDCommand, 0);
  (**mCommonModeVoltageCommand)(1, 0) =
      control(Signal::MMCControllerSystem::VMSigmaQCommand, 0);
  (**mCommonModeVoltageCommand)(2, 0) =
      control(Signal::MMCControllerSystem::VMSigmaZCommand, 0);
}

std::vector<String> EMT::Ph3::MMC_Modular::getLocalStateNames() const {
  return mModel.stateNames();
}
Matrix EMT::Ph3::MMC_Modular::getState() const { return **mX; }
Matrix EMT::Ph3::MMC_Modular::getStateDerivative() const {
  Matrix dx;
  mModel.evaluateStateDerivative(**mX, **mIntfVoltage, dx);
  return dx;
}
Matrix EMT::Ph3::MMC_Modular::getInterfaceVoltage() const {
  return **mIntfVoltage;
}
Matrix EMT::Ph3::MMC_Modular::getInterfaceCurrent() const {
  return **mIntfCurrent;
}
void EMT::Ph3::MMC_Modular::getLocalLinearization(Matrix &A, Matrix &B,
                                                  Matrix &C, Matrix &D) const {
  calculateNumericalJacobians(**mX, **mIntfVoltage, A, B, C, D);
}

void EMT::Ph3::MMC_Modular::evaluateStateDerivative(const Matrix &x,
                                                    const Matrix &u,
                                                    Matrix &dx) const {
  mModel.evaluateStateDerivative(x, u, dx);
}

void EMT::Ph3::MMC_Modular::evaluateOutput(const Matrix &x, const Matrix &u,
                                           Matrix &y) const {
  mModel.evaluateOutput(x, u, y);
}

void EMT::Ph3::MMC_Modular::calculateNumericalJacobians(const Matrix &x,
                                                        const Matrix &u,
                                                        Matrix &A, Matrix &B,
                                                        Matrix &C,
                                                        Matrix &D) const {
  mModel.calculateNumericalJacobians(x, u, A, B, C, D);
}

void EMT::Ph3::MMC_Modular::buildStateSpaceModel(const Matrix &x,
                                                 const Matrix &u, Matrix &A,
                                                 Matrix &B, Matrix &C,
                                                 Matrix &D, Matrix &E,
                                                 Matrix &F) const {
  mModel.buildStateSpaceModel(x, u, A, B, C, D, E, F);
}

Signal::MMCLinearization
EMT::Ph3::MMC_Modular::getStateSpaceModel(const Matrix &x,
                                          const Matrix &u) const {
  return mModel.getStateSpaceModel(x, u);
}

Signal::MMCSparseLinearization
EMT::Ph3::MMC_Modular::getSparseStateSpaceModel(const Matrix &x,
                                                const Matrix &u) const {
  return mModel.getSparseStateSpaceModel(x, u);
}

Signal::MMCLinearization EMT::Ph3::MMC_Modular::getStateSpaceModel() const {
  return getStateSpaceModel(**mX, **mIntfVoltage);
}

Signal::MMCSparseLinearization
EMT::Ph3::MMC_Modular::getSparseStateSpaceModel() const {
  return getSparseStateSpaceModel(**mX, **mIntfVoltage);
}

std::uint64_t EMT::Ph3::MMC_Modular::linearizationRebuildCount() const {
  return mLinearizationRebuildCount;
}

std::uint64_t EMT::Ph3::MMC_Modular::linearizationReuseCount() const {
  return mLinearizationReuseCount;
}

void EMT::Ph3::MMC_Modular::resetLinearizationCounters() {
  mLinearizationRebuildCount = 0;
  mLinearizationReuseCount = 0;
}

Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::acTerminalVoltageAttribute() const {
  return mAcTerminalVoltage;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::acTerminalCurrentAttribute() const {
  return mAcTerminalCurrent;
}
Attribute<Real>::Ptr EMT::Ph3::MMC_Modular::dcPositiveVoltageAttribute() const {
  return mDcPositiveVoltage;
}
Attribute<Real>::Ptr EMT::Ph3::MMC_Modular::dcNegativeVoltageAttribute() const {
  return mDcNegativeVoltage;
}
Attribute<Real>::Ptr EMT::Ph3::MMC_Modular::dcVoltageAttribute() const {
  return mDcVoltage;
}
Attribute<Real>::Ptr EMT::Ph3::MMC_Modular::dcCurrentAttribute() const {
  return mDcCurrent;
}
Attribute<Real>::Ptr EMT::Ph3::MMC_Modular::activePowerAttribute() const {
  return mActivePower;
}
Attribute<Real>::Ptr EMT::Ph3::MMC_Modular::reactivePowerAttribute() const {
  return mReactivePower;
}
Attribute<Real>::Ptr EMT::Ph3::MMC_Modular::storedEnergyAttribute() const {
  return mStoredEnergy;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::interfaceVoltageAttribute() const {
  return mIntfVoltage;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::interfaceCurrentAttribute() const {
  return mIntfCurrent;
}
Attribute<Real>::Ptr
EMT::Ph3::MMC_Modular::filteredDaxisVoltageAttribute() const {
  return mFilteredDaxisVoltage;
}
Attribute<Real>::Ptr
EMT::Ph3::MMC_Modular::heldActiveCurrentReferenceAttribute() const {
  return mHeldActiveCurrentReference;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::appliedModulationAttribute() const {
  return mAppliedModulation;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::differentialVoltageCommandAttribute() const {
  return mDifferentialVoltageCommand;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::commonModeVoltageCommandAttribute() const {
  return mCommonModeVoltageCommand;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::appliedDifferentialVoltageAttribute() const {
  return mDifferentialVoltageCommand;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::internalDifferentialVoltageAttribute() const {
  return mInternalDifferentialVoltage;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::appliedCommonModeVoltageAttribute() const {
  return mCommonModeVoltageCommand;
}

void EMT::Ph3::MMC_Modular::addHeldControlDependencies(
    AttributeBase::List &prevStepDependencies) const {
  prevStepDependencies.push_back(mExternalDifferentialVoltage);
  prevStepDependencies.push_back(mExternalCommonModeVoltage);
}

Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::externalDifferentialVoltageAttribute() const {
  return mExternalDifferentialVoltage;
}
Attribute<Matrix>::Ptr
EMT::Ph3::MMC_Modular::externalCommonModeVoltageAttribute() const {
  return mExternalCommonModeVoltage;
}
Attribute<Bool>::Ptr
EMT::Ph3::MMC_Modular::externalCommandActiveAttribute() const {
  return mExternalCommandActive;
}
