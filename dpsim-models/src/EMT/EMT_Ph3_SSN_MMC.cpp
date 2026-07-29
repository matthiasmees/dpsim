// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

// NEW FRAME-CORRECTED VERSION: separate nominal EMT angle and Harmony PLL deviation,
// consistent dq-frame transformations, continuous/discrete eigenvalue diagnostics.
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Eigenvalues>
#include <Eigen/QR>

#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>
#include <dpsim-models/MathUtils.h>

using namespace CPS;

EMT::Ph3::SSN_MMC::SSN_MMC(String uid, String name, Logger::Level logLevel)
    : VTypeVariableSSNComp(uid, name, 5, 5, logLevel), mNominalFrequency(0.0),
      mOmegaN(0.0), mNominalAcVoltage(0.0), mNominalDcVoltage(0.0),
      mArmInductance(0.0), mArmResistance(0.0), mSubmoduleCapacitance(0.0),
      mNumberOfSubmodules(0), mReactorInductance(0.0), mReactorResistance(0.0),
      mActiveControlMode(ActiveControlMode::OpenLoop),
      mReactiveControlMode(ReactiveControlMode::OpenLoop),
      mControlSource(ControlSource::InternalControllers),
      mEnergyControllerEnabled(false), mPllEnabled(false),
      mModulationDelayEnabled(false), mActivePowerReference(0.0),
      mReactivePowerReference(0.0), mDcVoltageReference(0.0),
      mAcVoltageReference(0.0), mDroopGain(0.0), mOpenLoopIDeltaDReference(0.0),
      mOpenLoopIDeltaQReference(0.0),
      mActivePowerFeedforwardFilterCutoffFrequency(0.0),
      mActivePowerFeedforwardControlTimeStep(0.0),
      mActivePowerFeedforwardMinimumDaxisVoltage(0.0),
      mActivePowerFeedforwardFilteredVd(0.0),
      mActivePowerFeedforwardFilteredVdDerivative(0.0),
      mActivePowerFeedforwardHeldIDeltaDReference(0.0),
      mActivePowerFeedforwardStepCounter(0),
      mActivePowerFeedforwardInitialized(false),
      mActivePowerFeedforwardRuntimeEnabled(false), mISigmaDReference(0.0),
      mISigmaQReference(0.0), mISigmaZReference(0.0), mInitialAngle(0.0),
      mInitialOperatingPointEnabled(false), mInitialActivePower(0.0),
      mInitialReactivePower(0.0), mAcVoltageDqFilterTimeConstant(0.0),
      mActivePowerFilterTimeConstant(0.0),
      mReactivePowerFilterTimeConstant(0.0), mDcVoltageFilterTimeConstant(0.0),
      mAcVoltageMagnitudeFilterTimeConstant(0.0), mModulationDelay(0.0),
      mMaximumAcCurrent(std::numeric_limits<Real>::infinity()),
      mMaximumCirculatingCurrent(std::numeric_limits<Real>::infinity()),
      mMaximumModulationMagnitude(std::numeric_limits<Real>::infinity()),
      mJacobianRelativeStep(1e-6), mJacobianAbsoluteStep(1e-8),
      mOperatingPointInitializationEnabled(true),
      mOperatingPointMaximumIterations(25),
      mOperatingPointNormalizedTolerance(1e-9),
      mEigenvalueDiagnosticsEnabled(false), mEigenvalueDiagnosticsInterval(100),
      mModelUpdateCounter(0), mDiagnosticTimeStep(0.0),
      mDcVoltage(mAttributes->create<Real>("vdc")),
      mDcPositiveVoltage(mAttributes->create<Real>("vdcp")),
      mDcNegativeVoltage(mAttributes->create<Real>("vdcn")),
      mDcCurrent(mAttributes->create<Real>("idc")),
      mActivePower(mAttributes->create<Real>("p_ac")),
      mReactivePower(mAttributes->create<Real>("q_ac")),
      mAcVoltageMagnitude(mAttributes->create<Real>("v_ac_mag")),
      mStoredEnergy(mAttributes->create<Real>("stored_energy")),
      mConverterAngle(mAttributes->create<Real>("theta")),
      mPllFrequency(mAttributes->create<Real>("pll_frequency")),
      mFilteredActivePower(mAttributes->create<Real>("p_filtered")),
      mFilteredReactivePower(mAttributes->create<Real>("q_filtered")),
      mFilteredDcVoltage(mAttributes->create<Real>("vdc_filtered")),
      mFilteredDaxisVoltage(
          mAttributes->create<Real>("v_d_feedforward_filtered")),
      mHeldActiveCurrentReference(
          mAttributes->create<Real>("i_delta_d_feedforward_held")),
      mStateNorm(mAttributes->create<Real>("state_norm")),
      mStateDerivativeNorm(mAttributes->create<Real>("state_derivative_norm")),
      mEquilibriumResidualNorm(
          mAttributes->create<Real>("equilibrium_residual_norm")),
      mJacobianMaximumRealEigenvalue(
          mAttributes->create<Real>("jacobian_max_real_eigenvalue")),
      mJacobianMaximumMagnitudeEigenvalue(
          mAttributes->create<Real>("jacobian_max_abs_eigenvalue")),
      mJacobianDominantFrequency(
          mAttributes->create<Real>("jacobian_dominant_frequency")),
      mJacobianMaximumDiscreteMagnitude(
          mAttributes->create<Real>("jacobian_max_discrete_magnitude")),
      mJacobianDiscreteDominantFrequency(
          mAttributes->create<Real>("jacobian_discrete_dominant_frequency")),
      mGridAngle(mAttributes->create<Real>("grid_angle")),
      mPllAngleDeviation(mAttributes->create<Real>("pll_angle_deviation")),
      mPllError(mAttributes->create<Real>("pll_error")),
      mGridVoltageD(mAttributes->create<Real>("v_grid_d")),
      mGridVoltageQ(mAttributes->create<Real>("v_grid_q")),
      mControlVoltageD(mAttributes->create<Real>("v_control_d")),
      mControlVoltageQ(mAttributes->create<Real>("v_control_q")),
      mDeltaCurrentD(mAttributes->create<Real>("i_delta_d")),
      mDeltaCurrentQ(mAttributes->create<Real>("i_delta_q")),
      mSigmaCurrentZ(mAttributes->create<Real>("i_sigma_z")),
      mDeltaCurrentReferenceD(mAttributes->create<Real>("i_delta_d_ref")),
      mDeltaCurrentReferenceQ(mAttributes->create<Real>("i_delta_q_ref")),
      mSigmaCurrentReferenceZ(mAttributes->create<Real>("i_sigma_z_ref")),
      mDcPower(mAttributes->create<Real>("p_dc")),
      mPowerBalanceError(mAttributes->create<Real>("power_balance_error")),
      mNortonMatrixNorm(mAttributes->create<Real>("norton_matrix_norm")),
      mHistoryVectorNorm(mAttributes->create<Real>("history_vector_norm")),
      mDiagnosticsValid(mAttributes->create<Real>("diagnostics_valid")),
      mAcTerminalVoltage(mAttributes->create<Matrix>("ac_terminal_voltage",
                                                     Matrix::Zero(3, 1))),
      mAcTerminalCurrent(mAttributes->create<Matrix>("ac_terminal_current",
                                                     Matrix::Zero(3, 1))),
      mExternalDifferentialVoltage(
          mAttributes->createDynamic<Matrix>("external_differential_voltage")),
      mExternalCommonModeVoltage(
          mAttributes->createDynamic<Matrix>("external_common_mode_voltage")),
      mAppliedDifferentialVoltage(mAttributes->create<Matrix>(
          "applied_differential_voltage", Matrix::Zero(2, 1))),
      mInternalDifferentialVoltage(mAttributes->create<Matrix>(
          "internal_differential_voltage", Matrix::Zero(2, 1))),
      mAppliedCommonModeVoltage(mAttributes->create<Matrix>(
          "applied_common_mode_voltage", Matrix::Zero(3, 1))),
      mAppliedModulation(mAttributes->create<Matrix>("applied_modulation",
                                                     Matrix::Zero(5, 1))),
      mRealizedConverterVoltage(mAttributes->create<Matrix>(
          "realized_converter_voltage", Matrix::Zero(5, 1))),
      mExternalCommandActive(
          mAttributes->create<Bool>("external_command_active", false)) {

  **mDcVoltage = 0.0;
  **mDcPositiveVoltage = 0.0;
  **mDcNegativeVoltage = 0.0;
  **mDcCurrent = 0.0;
  **mActivePower = 0.0;
  **mReactivePower = 0.0;
  **mAcVoltageMagnitude = 0.0;
  **mStoredEnergy = 0.0;
  **mConverterAngle = 0.0;
  **mPllFrequency = 0.0;
  **mFilteredActivePower = 0.0;
  **mFilteredReactivePower = 0.0;
  **mFilteredDcVoltage = 0.0;
  **mFilteredDaxisVoltage = 0.0;
  **mHeldActiveCurrentReference = 0.0;
  **mStateNorm = 0.0;
  **mStateDerivativeNorm = 0.0;
  **mEquilibriumResidualNorm = 0.0;
  **mJacobianMaximumRealEigenvalue = 0.0;
  **mJacobianMaximumMagnitudeEigenvalue = 0.0;
  **mJacobianDominantFrequency = 0.0;
  **mJacobianMaximumDiscreteMagnitude = 0.0;
  **mJacobianDiscreteDominantFrequency = 0.0;
  **mGridAngle = 0.0;
  **mPllAngleDeviation = 0.0;
  **mPllError = 0.0;
  **mGridVoltageD = 0.0;
  **mGridVoltageQ = 0.0;
  **mControlVoltageD = 0.0;
  **mControlVoltageQ = 0.0;
  **mDeltaCurrentD = 0.0;
  **mDeltaCurrentQ = 0.0;
  **mSigmaCurrentZ = 0.0;
  **mDeltaCurrentReferenceD = 0.0;
  **mDeltaCurrentReferenceQ = 0.0;
  **mSigmaCurrentReferenceZ = 0.0;
  **mDcPower = 0.0;
  **mPowerBalanceError = 0.0;
  **mNortonMatrixNorm = 0.0;
  **mHistoryVectorNorm = 0.0;
  **mDiagnosticsValid = 0.0;
  **mExternalDifferentialVoltage = Matrix::Zero(2, 1);
  **mExternalCommonModeVoltage = Matrix::Zero(3, 1);

  // DPsim currently stores one phase type for the complete component. The AC
  // terminal uses all ABC phases; DC+ and DC- use only phase index 0. The
  // corresponding nodes must therefore expose a valid phase-0 MNA index.
  mPhaseType = PhaseType::ABC;
  setTerminalNumber(3);
}

MatrixComp EMT::Ph3::SSN_MMC::buildInitialInputFromNodes(Real) {
  validateTerminalArrangement();

  MatrixComp u = MatrixComp::Zero(mInputSize, 1);

  const Complex vAcA = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  u(Va, 0) = vAcA;
  u(Vb, 0) = vAcA * SHIFT_TO_PHASE_B;
  u(Vc, 0) = vAcA * SHIFT_TO_PHASE_C;

  // DC node voltages are instantaneous scalar values and must not receive the
  // RMS-to-peak AC scaling.
  u(Vdcp, 0) = initialSingleVoltage(1);
  u(Vdcn, 0) = initialSingleVoltage(2);
  return u;
}

void EMT::Ph3::SSN_MMC::validateTerminalArrangement() const {
  auto *self = const_cast<SSN_MMC *>(this);

  const auto acNode = self->node(0);
  if (acNode->isGround() || acNode->phaseType() != PhaseType::ABC)
    throw std::invalid_argument(
        "SSN_MMC terminal 0 requires a non-grounded PhaseType::ABC node.");

  for (UInt terminal = 1; terminal <= 2; ++terminal) {
    const auto dcNode = self->node(terminal);
    if (!dcNode->isGround() && dcNode->phaseType() != PhaseType::DC)
      throw std::invalid_argument(
          "SSN_MMC terminals 1 (dc+) and 2 (dc-) require PhaseType::DC "
          "nodes or ground.");
  }
}

void EMT::Ph3::SSN_MMC::validateInterfaceDimensions() const {
  if (mW.rows() != mOutputSize || mW.cols() != mInputSize)
    throw std::runtime_error("SSN_MMC Norton matrix must have size 5 x 5.");
  if (mYHist.rows() != mOutputSize || mYHist.cols() != 1)
    throw std::runtime_error("SSN_MMC history vector must have size 5 x 1.");
  if ((**mIntfVoltage).rows() != mInputSize || (**mIntfVoltage).cols() != 1)
    throw std::runtime_error(
        "SSN_MMC interface-voltage vector must have size 5 x 1.");
  if ((**mIntfCurrent).rows() != mOutputSize || (**mIntfCurrent).cols() != 1)
    throw std::runtime_error(
        "SSN_MMC interface-current vector must have size 5 x 1.");
}

void EMT::Ph3::SSN_MMC::mnaCompUpdateVoltage(const Matrix &leftVector) {
  validateInterfaceDimensions();
  Matrix u = Matrix::Zero(mInputSize, 1);

  if (terminalNotGrounded(0)) {
    for (UInt phase = 0; phase < 3; ++phase)
      u(phase, 0) =
          Math::realFromVectorElement(leftVector, matrixNodeIndex(0, phase));
  }

  if (terminalNotGrounded(1))
    u(Vdcp, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
  if (terminalNotGrounded(2))
    u(Vdcn, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(2, 0));

  **mIntfVoltage = u;
}

void EMT::Ph3::SSN_MMC::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  validateTerminalArrangement();
  validateInterfaceDimensions();
  std::array<Int, mInputSize> indices = {-1, -1, -1, -1, -1};

  if (terminalNotGrounded(0)) {
    indices[Va] = static_cast<Int>(matrixNodeIndex(0, 0));
    indices[Vb] = static_cast<Int>(matrixNodeIndex(0, 1));
    indices[Vc] = static_cast<Int>(matrixNodeIndex(0, 2));
  }
  if (terminalNotGrounded(1))
    indices[Vdcp] = static_cast<Int>(matrixNodeIndex(1, 0));
  if (terminalNotGrounded(2))
    indices[Vdcn] = static_cast<Int>(matrixNodeIndex(2, 0));

  for (UInt row = 0; row < mOutputSize; ++row) {
    if (indices[row] < 0)
      continue;
    for (UInt column = 0; column < mInputSize; ++column) {
      if (indices[column] < 0)
        continue;
      Math::addToMatrixElement(systemMatrix, static_cast<UInt>(indices[row]),
                               static_cast<UInt>(indices[column]),
                               mW(row, column));
    }
  }
}

void EMT::Ph3::SSN_MMC::mnaCompApplyRightSideVectorStamp(Matrix &rightVector) {
  validateTerminalArrangement();
  validateInterfaceDimensions();
  std::array<Int, mOutputSize> indices = {-1, -1, -1, -1, -1};

  if (terminalNotGrounded(0)) {
    indices[Ia] = static_cast<Int>(matrixNodeIndex(0, 0));
    indices[Ib] = static_cast<Int>(matrixNodeIndex(0, 1));
    indices[Ic] = static_cast<Int>(matrixNodeIndex(0, 2));
  }
  if (terminalNotGrounded(1))
    indices[Idcp] = static_cast<Int>(matrixNodeIndex(1, 0));
  if (terminalNotGrounded(2))
    indices[Idcn] = static_cast<Int>(matrixNodeIndex(2, 0));

  for (UInt row = 0; row < mOutputSize; ++row) {
    if (indices[row] < 0)
      continue;
    // rightVector is this component's private MNA contribution vector.
    // It persists between pre-steps, so overwrite this component's entry.
    // The MNA solver sums the private vectors of all components afterwards.
    Math::setVectorElement(rightVector, static_cast<UInt>(indices[row]),
                           -mYHist(row, 0));
  }
}

std::vector<String> EMT::Ph3::SSN_MMC::getLocalStateNames() const {
  return {
      "iDelta_d",         "iDelta_q",
      "iSigma_z",         "iSigma_d",
      "iSigma_q",         "vCDelta_d",
      "vCDelta_q",        "vCDelta_Zd",
      "vCDelta_Zq",       "vCSigma_d",
      "vCSigma_q",        "vCSigma_z",
      "xi_active",        "xi_reactive",
      "xi_occ_d",         "xi_occ_q",
      "xi_ccc_d",         "xi_ccc_q",
      "xi_zcc",           "xi_energy",
      "xi_pll",           "pll_angle_deviation",
      "grid_angle",       "filter_vac_d",
      "filter_vac_q",     "filter_p_1",
      "filter_p_2",       "filter_q_1",
      "filter_q_2",       "filter_vdc_1",
      "filter_vdc_2",     "filter_vac_mag_1",
      "filter_vac_mag_2", "delay_mDelta_d_1",
      "delay_mDelta_d_2", "delay_mDelta_q_1",
      "delay_mDelta_q_2", "delay_mSigma_d_1",
      "delay_mSigma_d_2", "delay_mSigma_q_1",
      "delay_mSigma_q_2", "delay_mSigma_z_1",
      "delay_mSigma_z_2",
  };
}

void EMT::Ph3::SSN_MMC::setParameters(
    Real nominalFrequency, Real nominalAcVoltage, Real nominalDcVoltage,
    Real armInductance, Real armResistance, Real submoduleCapacitance,
    UInt numberOfSubmodules, Real reactorInductance, Real reactorResistance) {

  if (nominalFrequency <= 0.0)
    throw std::invalid_argument("nominalFrequency must be positive.");
  if (nominalAcVoltage <= 0.0)
    throw std::invalid_argument("nominalAcVoltage must be positive.");
  if (nominalDcVoltage <= 0.0)
    throw std::invalid_argument("nominalDcVoltage must be positive.");
  if (armInductance <= 0.0)
    throw std::invalid_argument("armInductance must be positive.");
  if (armResistance < 0.0)
    throw std::invalid_argument("armResistance must be non-negative.");
  if (submoduleCapacitance <= 0.0)
    throw std::invalid_argument("submoduleCapacitance must be positive.");
  if (numberOfSubmodules == 0)
    throw std::invalid_argument("numberOfSubmodules must be positive.");
  if (reactorInductance < 0.0)
    throw std::invalid_argument("reactorInductance must be non-negative.");
  if (reactorResistance < 0.0)
    throw std::invalid_argument("reactorResistance must be non-negative.");

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

  mDcVoltageReference = nominalDcVoltage;
  mAcVoltageReference = nominalAcVoltage;

  VTypeVariableSSNComp::setParameters(Matrix::Zero(mStateSize, mStateSize),
                                      Matrix::Zero(mStateSize, mInputSize),
                                      Matrix::Zero(mOutputSize, mStateSize),
                                      Matrix::Zero(mOutputSize, mInputSize),
                                      Matrix::Zero(mStateSize, 1),
                                      Matrix::Zero(mOutputSize, 1));
}

void EMT::Ph3::SSN_MMC::setOutputCurrentController(Real kp, Real ki) {
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("OCC gains must be non-negative.");
  mOutputCurrentController = {kp, ki};
}

void EMT::Ph3::SSN_MMC::setCirculatingCurrentController(Real kp, Real ki) {
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("CCC gains must be non-negative.");
  mCirculatingCurrentController = {kp, ki};
}

void EMT::Ph3::SSN_MMC::setZeroSequenceCurrentController(Real kp, Real ki) {
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("ZCC gains must be non-negative.");
  mZeroSequenceCurrentController = {kp, ki};
}

void EMT::Ph3::SSN_MMC::setEnergyController(Real kp, Real ki, Bool enabled) {
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument(
        "Energy-controller gains must be non-negative.");
  mEnergyController = {kp, ki};
  mEnergyControllerEnabled = enabled;
}

void EMT::Ph3::SSN_MMC::setActivePowerControl(Real reference, Real kp,
                                              Real ki) {
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("Active-power gains must be non-negative.");
  mActiveControlMode = ActiveControlMode::ActivePower;
  mActivePowerFeedforwardRuntimeEnabled = false;
  mActivePowerReference = reference;
  mActiveController = {kp, ki};
}

void EMT::Ph3::SSN_MMC::setActivePowerFeedforwardControl(
    Real activePowerReference, Real voltageFilterCutoffFrequency,
    Real controlTimeStep, Real minimumDaxisVoltage) {
  if (!std::isfinite(activePowerReference))
    throw std::invalid_argument(
        "Active-power feedforward reference must be finite.");
  if (!std::isfinite(voltageFilterCutoffFrequency) ||
      voltageFilterCutoffFrequency <= 0.0)
    throw std::invalid_argument(
        "Active-power feedforward voltage-filter cutoff frequency must be "
        "finite and positive.");
  if (!std::isfinite(controlTimeStep) || controlTimeStep <= 0.0)
    throw std::invalid_argument(
        "Active-power feedforward control time step must be finite and "
        "positive.");
  if (!std::isfinite(minimumDaxisVoltage) || minimumDaxisVoltage < 0.0)
    throw std::invalid_argument(
        "Active-power feedforward minimum d-axis voltage must be finite and "
        "non-negative.");

  mActiveControlMode = ActiveControlMode::ActivePowerFeedforward;
  mActivePowerReference = activePowerReference;
  mActivePowerFeedforwardFilterCutoffFrequency = voltageFilterCutoffFrequency;
  mActivePowerFeedforwardControlTimeStep = controlTimeStep;
  mActivePowerFeedforwardMinimumDaxisVoltage = minimumDaxisVoltage;
  mActivePowerFeedforwardFilteredVd = 0.0;
  mActivePowerFeedforwardFilteredVdDerivative = 0.0;
  mActivePowerFeedforwardHeldIDeltaDReference = 0.0;
  mActivePowerFeedforwardStepCounter = 0;
  mActivePowerFeedforwardInitialized = false;
  mActivePowerFeedforwardRuntimeEnabled = false;
}

void EMT::Ph3::SSN_MMC::setActivePowerFeedforwardReference(
    Real activePowerReference) {
  if (!std::isfinite(activePowerReference))
    throw std::invalid_argument(
        "Active-power feedforward reference must be finite.");
  if (mActiveControlMode != ActiveControlMode::ActivePowerFeedforward)
    throw std::logic_error("setActivePowerFeedforwardReference() requires "
                           "ActivePowerFeedforward mode.");

  mActivePowerReference = activePowerReference;
}

void EMT::Ph3::SSN_MMC::setDcVoltageControl(Real reference, Real kp, Real ki) {
  if (reference <= 0.0)
    throw std::invalid_argument("DC-voltage reference must be positive.");
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("DC-voltage gains must be non-negative.");
  mActiveControlMode = ActiveControlMode::DcVoltage;
  mActivePowerFeedforwardRuntimeEnabled = false;
  mDcVoltageReference = reference;
  mActiveController = {kp, ki};
}

void EMT::Ph3::SSN_MMC::setDcDroopControl(Real activePowerReference,
                                          Real dcVoltageReference,
                                          Real droopGain) {
  if (dcVoltageReference <= 0.0)
    throw std::invalid_argument("DC-voltage reference must be positive.");
  if (droopGain == 0.0)
    throw std::invalid_argument("droopGain must be non-zero.");
  mActiveControlMode = ActiveControlMode::DcDroop;
  mActivePowerFeedforwardRuntimeEnabled = false;
  mActivePowerReference = activePowerReference;
  mDcVoltageReference = dcVoltageReference;
  mDroopGain = droopGain;
}

void EMT::Ph3::SSN_MMC::setActiveControlOpenLoop(Real reference) {
  mActiveControlMode = ActiveControlMode::OpenLoop;
  mActivePowerFeedforwardRuntimeEnabled = false;
  mOpenLoopIDeltaDReference = reference;
}

void EMT::Ph3::SSN_MMC::setReactivePowerControl(Real reference, Real kp,
                                                Real ki) {
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("Reactive-power gains must be non-negative.");
  mReactiveControlMode = ReactiveControlMode::ReactivePower;
  mReactivePowerReference = reference;
  mReactiveController = {kp, ki};
}

void EMT::Ph3::SSN_MMC::setAcVoltageControl(Real reference, Real kp, Real ki) {
  if (reference <= 0.0)
    throw std::invalid_argument("AC-voltage reference must be positive.");
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("AC-voltage gains must be non-negative.");
  mReactiveControlMode = ReactiveControlMode::AcVoltage;
  mAcVoltageReference = reference;
  mReactiveController = {kp, ki};
}

void EMT::Ph3::SSN_MMC::setReactiveControlOpenLoop(Real reference) {
  mReactiveControlMode = ReactiveControlMode::OpenLoop;
  mOpenLoopIDeltaQReference = reference;
}

void EMT::Ph3::SSN_MMC::setCirculatingCurrentReferences(Real dReference,
                                                        Real qReference,
                                                        Real zReference) {
  mISigmaDReference = dReference;
  mISigmaQReference = qReference;
  mISigmaZReference = zReference;
}

void EMT::Ph3::SSN_MMC::setPLL(Real kp, Real ki, Bool enabled) {
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("PLL gains must be non-negative.");
  mPllController = {kp, ki};
  mPllEnabled = enabled;
}

void EMT::Ph3::SSN_MMC::setMeasurementFilters(
    Real acVoltageDqTimeConstant, Real activePowerTimeConstant,
    Real reactivePowerTimeConstant, Real dcVoltageTimeConstant,
    Real acVoltageMagnitudeTimeConstant) {
  mAcVoltageDqFilterTimeConstant = std::max(0.0, acVoltageDqTimeConstant);
  mActivePowerFilterTimeConstant = std::max(0.0, activePowerTimeConstant);
  mReactivePowerFilterTimeConstant = std::max(0.0, reactivePowerTimeConstant);
  mDcVoltageFilterTimeConstant = std::max(0.0, dcVoltageTimeConstant);
  mAcVoltageMagnitudeFilterTimeConstant =
      std::max(0.0, acVoltageMagnitudeTimeConstant);
}

void EMT::Ph3::SSN_MMC::setModulationDelay(Real timeDelay, UInt padeOrder) {
  if (timeDelay < 0.0)
    throw std::invalid_argument("Modulation delay must be non-negative.");
  if (padeOrder != 2)
    throw std::invalid_argument(
        "This DPsim MMC currently supports Harmony's second-order Padé delay.");
  mModulationDelay = timeDelay;
  mModulationDelayEnabled = timeDelay > 0.0;
}

void EMT::Ph3::SSN_MMC::setInitialAngle(Real angle) { mInitialAngle = angle; }

void EMT::Ph3::SSN_MMC::setInitialOperatingPoint(Real activePower,
                                                 Real reactivePower) {
  if (!std::isfinite(activePower) || !std::isfinite(reactivePower))
    throw std::invalid_argument(
        "Initial MMC active/reactive power must be finite.");

  mInitialOperatingPointEnabled = true;
  mInitialActivePower = activePower;
  mInitialReactivePower = reactivePower;
}

void EMT::Ph3::SSN_MMC::setControlSource(ControlSource source) {
  const Bool externalDifferential =
      source == ControlSource::ExternalDifferentialVoltage ||
      source == ControlSource::ExternalFullConverterVoltage;

  if (externalDifferential && ((**mExternalDifferentialVoltage).rows() != 2 ||
                               (**mExternalDifferentialVoltage).cols() != 1 ||
                               !(**mExternalDifferentialVoltage).allFinite()))
    throw std::logic_error(
        "MMC external differential-voltage command must be finite 2x1.");

  if (source == ControlSource::ExternalFullConverterVoltage &&
      ((**mExternalCommonModeVoltage).rows() != 3 ||
       (**mExternalCommonModeVoltage).cols() != 1 ||
       !(**mExternalCommonModeVoltage).allFinite()))
    throw std::logic_error(
        "MMC external common-mode voltage command must be finite 3x1.");

  mControlSource = source;
  **mExternalCommandActive = source != ControlSource::InternalControllers;
}

void EMT::Ph3::SSN_MMC::setExternalDifferentialVoltageCommand(Real dVolts,
                                                              Real qVolts) {
  if (!std::isfinite(dVolts) || !std::isfinite(qVolts))
    throw std::invalid_argument(
        "MMC external differential-voltage command must be finite.");
  (**mExternalDifferentialVoltage)(0, 0) = dVolts;
  (**mExternalDifferentialVoltage)(1, 0) = qVolts;
}

void EMT::Ph3::SSN_MMC::setExternalCommonModeVoltageCommand(Real dVolts,
                                                            Real qVolts,
                                                            Real zVolts) {
  if (!std::isfinite(dVolts) || !std::isfinite(qVolts) ||
      !std::isfinite(zVolts))
    throw std::invalid_argument(
        "MMC external common-mode voltage command must be finite.");
  (**mExternalCommonModeVoltage)(0, 0) = dVolts;
  (**mExternalCommonModeVoltage)(1, 0) = qVolts;
  (**mExternalCommonModeVoltage)(2, 0) = zVolts;
}

void EMT::Ph3::SSN_MMC::setLimits(Real maximumAcCurrent,
                                  Real maximumCirculatingCurrent,
                                  Real maximumModulationMagnitude) {
  if (maximumAcCurrent <= 0.0 || maximumCirculatingCurrent <= 0.0 ||
      maximumModulationMagnitude <= 0.0)
    throw std::invalid_argument("MMC limits must be positive.");
  mMaximumAcCurrent = maximumAcCurrent;
  mMaximumCirculatingCurrent = maximumCirculatingCurrent;
  mMaximumModulationMagnitude = maximumModulationMagnitude;
}

void EMT::Ph3::SSN_MMC::setNumericalLinearizationParameters(Real relativeStep,
                                                            Real absoluteStep) {
  if (relativeStep <= 0.0 || absoluteStep <= 0.0)
    throw std::invalid_argument("Jacobian steps must be positive.");
  mJacobianRelativeStep = relativeStep;
  mJacobianAbsoluteStep = absoluteStep;
}

void EMT::Ph3::SSN_MMC::setOperatingPointInitialization(
    Bool enabled, UInt maximumIterations, Real normalizedTolerance) {
  if (maximumIterations == 0)
    throw std::invalid_argument(
        "Operating-point maximumIterations must be positive.");
  if (normalizedTolerance <= 0.0)
    throw std::invalid_argument(
        "Operating-point normalizedTolerance must be positive.");

  mOperatingPointInitializationEnabled = enabled;
  mOperatingPointMaximumIterations = maximumIterations;
  mOperatingPointNormalizedTolerance = normalizedTolerance;
}

void EMT::Ph3::SSN_MMC::setEigenvalueDiagnostics(Bool enabled,
                                                 UInt updateInterval) {
  if (updateInterval == 0)
    throw std::invalid_argument(
        "Eigenvalue diagnostic updateInterval must be positive.");

  mEigenvalueDiagnosticsEnabled = enabled;
  mEigenvalueDiagnosticsInterval = updateInterval;
  mModelUpdateCounter = 0;
}

void EMT::Ph3::SSN_MMC::setDiagnosticTimeStep(Real timeStep) {
  if (timeStep <= 0.0)
    throw std::invalid_argument("Diagnostic time step must be positive.");
  mDiagnosticTimeStep = timeStep;
}

Real EMT::Ph3::SSN_MMC::clamp(Real value, Real lower, Real upper) const {
  return std::max(lower, std::min(value, upper));
}

Real EMT::Ph3::SSN_MMC::regularizedDcVoltage(Real voltage) const {
  if (std::abs(voltage) >= mMinimumDcVoltage)
    return voltage;
  return voltage >= 0.0 ? mMinimumDcVoltage : -mMinimumDcVoltage;
}

Matrix EMT::Ph3::SSN_MMC::abcToDq(const Matrix &abc, Real theta) const {
  if (abc.rows() != 3 || abc.cols() != 1)
    throw std::invalid_argument("abcToDq expects a 3x1 vector.");

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

Matrix EMT::Ph3::SSN_MMC::dqToAbc(Real d, Real q, Real theta) const {
  Matrix abc(3, 1);
  abc(0, 0) = d * std::cos(theta) - q * std::sin(theta);
  abc(1, 0) = d * std::cos(theta - 2.0 * PI / 3.0) -
              q * std::sin(theta - 2.0 * PI / 3.0);
  abc(2, 0) = d * std::cos(theta + 2.0 * PI / 3.0) -
              q * std::sin(theta + 2.0 * PI / 3.0);
  return abc;
}

Matrix EMT::Ph3::SSN_MMC::rotateDq(const Matrix &dq, Real angle) const {
  if (dq.rows() != 2 || dq.cols() != 1)
    throw std::invalid_argument("rotateDq expects a 2x1 vector.");

  const Real c = std::cos(angle);
  const Real s = std::sin(angle);

  Matrix rotated(2, 1);
  rotated(0, 0) = c * dq(0, 0) - s * dq(1, 0);
  rotated(1, 0) = s * dq(0, 0) + c * dq(1, 0);
  return rotated;
}

Real EMT::Ph3::SSN_MMC::calculateStoredEnergy(const Matrix &x) const {
  // Harmony total stored-energy expression:
  //   wSigma_z = 3*Carm/(2*N) *
  //     (vCDd^2 + vCDq^2 + vCDZd^2 + vCDZq^2
  //      + vCSd^2 + vCSq^2 + 2*vCSz^2).
  const Real scale = 3.0 * mSubmoduleCapacitance /
                     (2.0 * static_cast<Real>(mNumberOfSubmodules));
  const Real sumSquares =
      x(VCDeltaD, 0) * x(VCDeltaD, 0) + x(VCDeltaQ, 0) * x(VCDeltaQ, 0) +
      x(VCDeltaZd, 0) * x(VCDeltaZd, 0) + x(VCDeltaZq, 0) * x(VCDeltaZq, 0) +
      x(VCSigmaD, 0) * x(VCSigmaD, 0) + x(VCSigmaQ, 0) * x(VCSigmaQ, 0) +
      2.0 * x(VCSigmaZ, 0) * x(VCSigmaZ, 0);
  return scale * sumSquares;
}

Real EMT::Ph3::SSN_MMC::applyFirstOrderFilter(const Matrix &x, Matrix &f,
                                              UInt stateIndex, Real input,
                                              Real timeConstant,
                                              Bool enabled) const {
  if (!enabled || timeConstant <= 0.0) {
    f(stateIndex, 0) = 0.0;
    return input;
  }
  f(stateIndex, 0) = (input - x(stateIndex, 0)) / timeConstant;
  return x(stateIndex, 0);
}

Real EMT::Ph3::SSN_MMC::applySecondOrderFilter(const Matrix &x, Matrix &f,
                                               UInt stateIndex1,
                                               UInt stateIndex2, Real input,
                                               Real timeConstant,
                                               Bool enabled) const {
  if (!enabled || timeConstant <= 0.0) {
    f(stateIndex1, 0) = 0.0;
    f(stateIndex2, 0) = 0.0;
    return input;
  }
  f(stateIndex1, 0) = (input - x(stateIndex1, 0)) / timeConstant;
  f(stateIndex2, 0) = (x(stateIndex1, 0) - x(stateIndex2, 0)) / timeConstant;
  return x(stateIndex2, 0);
}

Real EMT::Ph3::SSN_MMC::applyPadeDelayChannel(const Matrix &x, Matrix &f,
                                              UInt stateIndex1,
                                              UInt stateIndex2,
                                              Real input) const {
  if (!mModulationDelayEnabled || mModulationDelay <= 0.0) {
    f(stateIndex1, 0) = 0.0;
    f(stateIndex2, 0) = 0.0;
    return input;
  }

  // H(s) = (1 - Ts/2 + T^2 s^2/12) /
  //        (1 + Ts/2 + T^2 s^2/12).
  const Real a0 = 12.0 / (mModulationDelay * mModulationDelay);
  const Real a1 = 6.0 / mModulationDelay;
  const Real x1 = x(stateIndex1, 0);
  const Real x2 = x(stateIndex2, 0);
  f(stateIndex1, 0) = x2;
  f(stateIndex2, 0) = -a0 * x1 - a1 * x2 + input;
  return input - 2.0 * a1 * x2;
}

void EMT::Ph3::SSN_MMC::initializeSampledActivePowerFeedforward(
    const Matrix &x0, const Matrix &u0) {
  if (mActiveControlMode != ActiveControlMode::ActivePowerFeedforward)
    return;

  if (mActivePowerFeedforwardFilterCutoffFrequency <= 0.0 ||
      mActivePowerFeedforwardControlTimeStep <= 0.0)
    throw std::logic_error(
        "Active-power feedforward control must be configured before MMC "
        "initialization.");

  const Matrix vGridDq = abcToDq(u0.block(0, 0, 3, 1), x0(GridAngle, 0));
  const Real initialVd = vGridDq(0, 0);
  const Real minimumVd = mActivePowerFeedforwardMinimumDaxisVoltage > 0.0
                             ? mActivePowerFeedforwardMinimumDaxisVoltage
                             : 0.1 * std::sqrt(2.0 / 3.0) * mNominalAcVoltage;

  if (!std::isfinite(initialVd) || initialVd <= minimumVd)
    throw std::runtime_error(
        "Active-power feedforward hot start has an invalid initial d-axis "
        "voltage.");

  mActivePowerFeedforwardFilteredVd = initialVd;
  mActivePowerFeedforwardFilteredVdDerivative = 0.0;
  mActivePowerFeedforwardHeldIDeltaDReference =
      (2.0 / 3.0) * mActivePowerReference / initialVd;
  mActivePowerFeedforwardStepCounter = 0;
  mActivePowerFeedforwardInitialized = true;
  mActivePowerFeedforwardRuntimeEnabled = false;

  **mFilteredDaxisVoltage = mActivePowerFeedforwardFilteredVd;
  **mHeldActiveCurrentReference = mActivePowerFeedforwardHeldIDeltaDReference;
}

void EMT::Ph3::SSN_MMC::advanceSampledActivePowerFeedforward() {
  if (mActiveControlMode != ActiveControlMode::ActivePowerFeedforward ||
      !mActivePowerFeedforwardRuntimeEnabled ||
      mControlSource != ControlSource::InternalControllers)
    return;

  if (!mActivePowerFeedforwardInitialized)
    throw std::logic_error(
        "Active-power feedforward sampled controller was not initialized.");
  if (!std::isfinite(mTimeStep) || mTimeStep <= 0.0)
    throw std::logic_error(
        "Active-power feedforward requires a valid EMT time step.");

  const UInt stride = std::max<UInt>(
      1, static_cast<UInt>(
             std::llround(mActivePowerFeedforwardControlTimeStep / mTimeStep)));
  const Real effectiveControlTimeStep = static_cast<Real>(stride) * mTimeStep;

  if (mActivePowerFeedforwardStepCounter == 0) {
    const Real gridAngle = (**mX)(GridAngle, 0);
    const Matrix vGridDq =
        abcToDq((**mIntfVoltage).block(0, 0, 3, 1), gridAngle);
    const Real measuredVd = vGridDq(0, 0);
    if (!std::isfinite(measuredVd))
      throw std::runtime_error(
          "Active-power feedforward measured d-axis voltage is not finite.");

    // Exact zero-order-hold update of:
    //
    //   H(s) = wn^2 / (s^2 + 2*wn*s + wn^2), zeta = 1.
    //
    // For e = y - u and a = exp(-wn*T):
    //   e+    = a*((1 + wn*T)*e + T*ydot)
    //   ydot+ = a*((1 - wn*T)*ydot - wn^2*T*e).
    const Real omega = 2.0 * PI * mActivePowerFeedforwardFilterCutoffFrequency;
    const Real decay = std::exp(-omega * effectiveControlTimeStep);
    const Real filterError = mActivePowerFeedforwardFilteredVd - measuredVd;

    const Real updatedFilterError =
        decay * ((1.0 + omega * effectiveControlTimeStep) * filterError +
                 effectiveControlTimeStep *
                     mActivePowerFeedforwardFilteredVdDerivative);
    const Real updatedFilterDerivative =
        decay * ((1.0 - omega * effectiveControlTimeStep) *
                     mActivePowerFeedforwardFilteredVdDerivative -
                 omega * omega * effectiveControlTimeStep * filterError);

    mActivePowerFeedforwardFilteredVd = measuredVd + updatedFilterError;
    mActivePowerFeedforwardFilteredVdDerivative = updatedFilterDerivative;

    const Real minimumVd = mActivePowerFeedforwardMinimumDaxisVoltage > 0.0
                               ? mActivePowerFeedforwardMinimumDaxisVoltage
                               : 0.1 * std::sqrt(2.0 / 3.0) * mNominalAcVoltage;
    if (!std::isfinite(mActivePowerFeedforwardFilteredVd) ||
        mActivePowerFeedforwardFilteredVd <= minimumVd)
      throw std::runtime_error(
          "Active-power feedforward received an invalid filtered d-axis "
          "voltage.");

    mActivePowerFeedforwardHeldIDeltaDReference =
        (2.0 / 3.0) * mActivePowerReference / mActivePowerFeedforwardFilteredVd;

    **mFilteredDaxisVoltage = mActivePowerFeedforwardFilteredVd;
    **mHeldActiveCurrentReference = mActivePowerFeedforwardHeldIDeltaDReference;
  }

  mActivePowerFeedforwardStepCounter =
      (mActivePowerFeedforwardStepCounter + 1) % stride;
}

void EMT::Ph3::SSN_MMC::validateFinite(const Matrix &value,
                                       const char *name) const {
  if (value.allFinite())
    return;

  for (Eigen::Index row = 0; row < value.rows(); ++row) {
    for (Eigen::Index column = 0; column < value.cols(); ++column) {
      if (!std::isfinite(value(row, column))) {
        SPDLOG_LOGGER_CRITICAL(mSLog, "{} contains NaN/Inf at ({}, {}): {}",
                               name, row, column, value(row, column));
        throw std::runtime_error(String(name) + " contains NaN or Inf at (" +
                                 std::to_string(row) + ", " +
                                 std::to_string(column) + ").");
      }
    }
  }

  throw std::runtime_error(String(name) + " contains NaN or Inf.");
}

std::vector<UInt> EMT::Ph3::SSN_MMC::getEquilibriumStateIndices() const {
  std::vector<UInt> indices;
  indices.reserve(mStateSize);

  for (UInt index = 0; index < mPlantStateSize; ++index)
    indices.push_back(index);

  if (mActiveControlMode == ActiveControlMode::ActivePower ||
      mActiveControlMode == ActiveControlMode::DcVoltage)
    indices.push_back(XiActive);

  if (mReactiveControlMode == ReactiveControlMode::ReactivePower ||
      mReactiveControlMode == ReactiveControlMode::AcVoltage)
    indices.push_back(XiReactive);

  if (mOutputCurrentController.kp > 0.0 || mOutputCurrentController.ki > 0.0) {
    indices.push_back(XiOccD);
    indices.push_back(XiOccQ);
  }

  if (mCirculatingCurrentController.kp > 0.0 ||
      mCirculatingCurrentController.ki > 0.0) {
    indices.push_back(XiCccD);
    indices.push_back(XiCccQ);
  }

  if (mZeroSequenceCurrentController.kp > 0.0 ||
      mZeroSequenceCurrentController.ki > 0.0)
    indices.push_back(XiZcc);

  if (mEnergyControllerEnabled)
    indices.push_back(XiEnergy);

  if (mAcVoltageDqFilterTimeConstant > 0.0) {
    indices.push_back(FilterVacD);
    indices.push_back(FilterVacQ);
  }
  if (mActivePowerFilterTimeConstant > 0.0) {
    indices.push_back(FilterP1);
    indices.push_back(FilterP2);
  }
  if (mReactivePowerFilterTimeConstant > 0.0) {
    indices.push_back(FilterQ1);
    indices.push_back(FilterQ2);
  }
  if (mDcVoltageFilterTimeConstant > 0.0) {
    indices.push_back(FilterVdc1);
    indices.push_back(FilterVdc2);
  }
  if (mAcVoltageMagnitudeFilterTimeConstant > 0.0) {
    indices.push_back(FilterVacMag1);
    indices.push_back(FilterVacMag2);
  }

  if (mModulationDelayEnabled) {
    for (UInt index = DelayMDeltaD1; index <= DelayMSigmaZ2; ++index)
      indices.push_back(index);
  }

  return indices;
}

std::vector<UInt> EMT::Ph3::SSN_MMC::getDiagnosticStateIndices() const {
  std::vector<UInt> indices = getEquilibriumStateIndices();

  if (mPllEnabled) {
    indices.push_back(XiPll);
    indices.push_back(PllAngle);
  }

  return indices;
}

void EMT::Ph3::SSN_MMC::initializeAnalyticalOperatingPoint(
    Matrix &x0, const Matrix &u0) const {
  x0.setZero(mStateSize, 1);
  x0(PllAngle, 0) = 0.0;
  x0(GridAngle, 0) = mInitialAngle;

  const Matrix vDq = abcToDq(u0.block(0, 0, 3, 1), x0(GridAngle, 0));
  const Real vGd = vDq(0, 0);
  const Real vGq = vDq(1, 0);
  const Real vDc = regularizedDcVoltage(u0(Vdcp, 0) - u0(Vdcn, 0));
  const Real voltageSquared = vGd * vGd + vGq * vGq;

  Real pTarget = 0.0;
  if (mInitialOperatingPointEnabled)
    pTarget = mInitialActivePower;
  else if (mActiveControlMode == ActiveControlMode::ActivePower ||
           mActiveControlMode == ActiveControlMode::ActivePowerFeedforward)
    pTarget = mActivePowerReference;

  Real qTarget = 0.0;
  if (mInitialOperatingPointEnabled)
    qTarget = mInitialReactivePower;
  else if (mReactiveControlMode == ReactiveControlMode::ReactivePower)
    qTarget = mReactivePowerReference;

  Real iDeltaD = mOpenLoopIDeltaDReference;
  Real iDeltaQ = mOpenLoopIDeltaQReference;
  if (voltageSquared > 0.0 &&
      (mInitialOperatingPointEnabled ||
       mActiveControlMode == ActiveControlMode::ActivePower ||
       mActiveControlMode == ActiveControlMode::ActivePowerFeedforward ||
       mActiveControlMode == ActiveControlMode::DcVoltage ||
       mReactiveControlMode == ReactiveControlMode::ReactivePower)) {
    // Same P/Q-to-current operating-point mapping as Harmony::update_MMC().
    iDeltaD = (2.0 / 3.0) * (vGd * pTarget + vGq * qTarget) / voltageSquared;
    iDeltaQ = (2.0 / 3.0) * (vGq * pTarget - vGd * qTarget) / voltageSquared;
  }

  Real iSigmaZ = mISigmaZReference;
  if (mEnergyControllerEnabled &&
      (mInitialOperatingPointEnabled ||
       mActiveControlMode == ActiveControlMode::ActivePower ||
       mActiveControlMode == ActiveControlMode::ActivePowerFeedforward ||
       mActiveControlMode == ActiveControlMode::DcVoltage))
    // Positive pTarget exports AC power. Positive iSigmaZ is DC current
    // absorbed by the converter, so both have the same sign at equilibrium.
    iSigmaZ = pTarget / (3.0 * vDc);

  x0(IDeltaD, 0) = iDeltaD;
  x0(IDeltaQ, 0) = iDeltaQ;
  x0(ISigmaZ, 0) = iSigmaZ;
  x0(ISigmaD, 0) = mISigmaDReference;
  x0(ISigmaQ, 0) = mISigmaQReference;

  // This state alone gives the exact Harmony energy reference because the
  // energy expression contains 2*vCSigma_z^2.
  x0(VCSigmaZ, 0) = std::abs(vDc);

  const Real pInitial = 1.5 * (vGd * iDeltaD + vGq * iDeltaQ);
  const Real qInitial = 1.5 * (-vGd * iDeltaQ + vGq * iDeltaD);
  const Real vacControlMagnitude = std::hypot(vGd, vGq);

  x0(FilterVacD, 0) = vGd;
  x0(FilterVacQ, 0) = vGq;
  x0(FilterP1, 0) = x0(FilterP2, 0) = pInitial;
  x0(FilterQ1, 0) = x0(FilterQ2, 0) = qInitial;
  x0(FilterVdc1, 0) = x0(FilterVdc2, 0) = vDc;
  x0(FilterVacMag1, 0) = x0(FilterVacMag2, 0) = vacControlMagnitude;

  if (mActiveController.ki > 0.0) {
    if (mActiveControlMode == ActiveControlMode::ActivePower ||
        mActiveControlMode == ActiveControlMode::DcVoltage)
      x0(XiActive, 0) = iDeltaD / mActiveController.ki;
  }

  if (mReactiveController.ki > 0.0) {
    if (mReactiveControlMode == ReactiveControlMode::ReactivePower)
      x0(XiReactive, 0) = -iDeltaQ / mReactiveController.ki;
    else if (mReactiveControlMode == ReactiveControlMode::AcVoltage)
      x0(XiReactive, 0) = iDeltaQ / mReactiveController.ki;
  }

  const Real lEqAc = mArmInductance / 2.0 + mReactorInductance;
  const Real rEqAc = mArmResistance / 2.0 + mReactorResistance;

  if (mOutputCurrentController.ki > 0.0) {
    // The MATLAB current regulator explicitly feedforwards R*i_ref and
    // omega*L*i_ref, so the PI integrators are zero at a matched equilibrium.
    x0(XiOccD, 0) = 0.0;
    x0(XiOccQ, 0) = 0.0;
  }

  if (mCirculatingCurrentController.ki > 0.0) {
    x0(XiCccD, 0) =
        mArmResistance * x0(ISigmaD, 0) / mCirculatingCurrentController.ki;
    x0(XiCccQ, 0) =
        mArmResistance * x0(ISigmaQ, 0) / mCirculatingCurrentController.ki;
  }

  if (mZeroSequenceCurrentController.ki > 0.0)
    x0(XiZcc, 0) = mArmResistance * iSigmaZ / mZeroSequenceCurrentController.ki;

  if (mEnergyControllerEnabled && mEnergyController.ki > 0.0)
    x0(XiEnergy, 0) = (3.0 * vDc * iSigmaZ - pInitial) / mEnergyController.ki;

  x0(XiPll, 0) = 0.0;

  if (mModulationDelayEnabled) {
    const Real omega = mOmegaN;
    const Real vMDeltaD = vGd + rEqAc * iDeltaD - lEqAc * omega * iDeltaQ;
    const Real vMDeltaQ = vGq + rEqAc * iDeltaQ + lEqAc * omega * iDeltaD;
    const Real vMSigmaD = -mArmResistance * x0(ISigmaD, 0) -
                          2.0 * mArmInductance * omega * x0(ISigmaQ, 0);
    const Real vMSigmaQ = -mArmResistance * x0(ISigmaQ, 0) +
                          2.0 * mArmInductance * omega * x0(ISigmaD, 0);
    const Real vMSigmaZ = vDc / 2.0 - mArmResistance * iSigmaZ;

    const Real modulation[5] = {
        -2.0 * vMDeltaD / vDc, -2.0 * vMDeltaQ / vDc, 2.0 * vMSigmaD / vDc,
        2.0 * vMSigmaQ / vDc,  2.0 * vMSigmaZ / vDc,
    };
    const UInt firstState[5] = {
        DelayMDeltaD1, DelayMDeltaQ1, DelayMSigmaD1,
        DelayMSigmaQ1, DelayMSigmaZ1,
    };
    const UInt secondState[5] = {
        DelayMDeltaD2, DelayMDeltaQ2, DelayMSigmaD2,
        DelayMSigmaQ2, DelayMSigmaZ2,
    };
    const Real a0 = 12.0 / (mModulationDelay * mModulationDelay);
    for (UInt channel = 0; channel < 5; ++channel) {
      x0(firstState[channel], 0) = modulation[channel] / a0;
      x0(secondState[channel], 0) = 0.0;
    }
  }
}

Bool EMT::Ph3::SSN_MMC::solveOperatingPoint(Matrix &x0, const Matrix &u0,
                                            Real &normalizedResidual,
                                            Real &absoluteResidual) const {
  const std::vector<UInt> indices = getEquilibriumStateIndices();
  const UInt dimension = static_cast<UInt>(indices.size());
  if (dimension == 0) {
    normalizedResidual = 0.0;
    absoluteResidual = 0.0;
    return true;
  }

  auto calculateResidual = [&](const Matrix &state, Matrix &rawResidual,
                               Matrix &scaledResidual, Matrix &rowScale) {
    Matrix derivative = Matrix::Zero(mStateSize, 1);
    evaluateStateDerivative(state, u0, derivative);

    const Matrix vDq = abcToDq(u0.block(0, 0, 3, 1), state(GridAngle, 0));
    const Real activePower =
        1.5 * (vDq(0, 0) * state(IDeltaD, 0) + vDq(1, 0) * state(IDeltaQ, 0));

    rawResidual.resize(dimension, 1);
    scaledResidual.resize(dimension, 1);
    rowScale.resize(dimension, 1);

    for (UInt row = 0; row < dimension; ++row) {
      const UInt stateIndex = indices[row];

      // In V-type DC-voltage mode the terminal Vdc is an imposed input.
      // Consequently d(XiActive)/dt = Vdc_ref - Vdc is identically zero at
      // initialization and does not constrain the non-zero active-current
      // operating point. A plain least-squares solve therefore drives
      // iDeltaD and XiActive toward their minimum-norm, nearly zero solution.
      //
      // Harmony avoids this degeneracy because its DC-voltage-controlled MMC
      // changes causality: Idc is an input and Vdc is an internal state. For
      // the present DPsim V-type realization, replace the redundant
      // XiActive equilibrium equation by the requested loaded active-power
      // constraint during initialization.
      if (mInitialOperatingPointEnabled &&
          mActiveControlMode == ActiveControlMode::DcVoltage &&
          stateIndex == XiActive) {
        rawResidual(row, 0) = activePower - mInitialActivePower;
        rowScale(row, 0) = std::max(1.0, std::abs(mInitialActivePower));
      } else {
        rawResidual(row, 0) = derivative(stateIndex, 0);
        rowScale(row, 0) = std::max(
            1.0, mOmegaN * std::max(1.0, std::abs(state(stateIndex, 0))));
      }

      scaledResidual(row, 0) = rawResidual(row, 0) / rowScale(row, 0);
    }
  };

  Matrix rawResidual;
  Matrix scaledResidual;
  Matrix rowScale;
  calculateResidual(x0, rawResidual, scaledResidual, rowScale);

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
          mJacobianRelativeStep * std::max(1.0, std::abs(x0(stateIndex, 0)));

      Matrix xPlus = x0;
      Matrix xMinus = x0;
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

      columnScale(column, 0) = std::max(1.0, std::abs(x0(stateIndex, 0)));

      for (UInt row = 0; row < dimension; ++row) {
        // Keep the current-iterate row scaling fixed while differentiating.
        jacobian(row, column) = (rawPlus(row, 0) - rawMinus(row, 0)) /
                                (2.0 * step * rowScale(row, 0));
      }
    }

    Matrix scaledJacobian = jacobian;
    for (UInt column = 0; column < dimension; ++column)
      scaledJacobian.col(column) *= columnScale(column, 0);

    const Matrix scaledCorrection =
        scaledJacobian.colPivHouseholderQr().solve(-scaledResidual);
    validateFinite(scaledCorrection, "MMC operating-point correction");

    Matrix correction = Matrix::Zero(dimension, 1);
    for (UInt column = 0; column < dimension; ++column)
      correction(column, 0) =
          columnScale(column, 0) * scaledCorrection(column, 0);

    Bool accepted = false;
    Real lineSearchFactor = 1.0;
    for (UInt lineSearch = 0; lineSearch < 12; ++lineSearch) {
      Matrix candidate = x0;
      for (UInt column = 0; column < dimension; ++column)
        candidate(indices[column], 0) +=
            lineSearchFactor * correction(column, 0);

      Matrix candidateRaw;
      Matrix candidateScaled;
      Matrix candidateScale;
      calculateResidual(candidate, candidateRaw, candidateScaled,
                        candidateScale);

      if (candidateScaled.norm() < normalizedResidual) {
        x0 = candidate;
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

void EMT::Ph3::SSN_MMC::updateDiagnostics(const Matrix &A, const Matrix &x,
                                          const Matrix &u, Bool force) {
  ++mModelUpdateCounter;

  Matrix derivative = Matrix::Zero(mStateSize, 1);
  evaluateStateDerivative(x, u, derivative);
  **mStateNorm = x.norm();
  **mStateDerivativeNorm = derivative.norm();

  if (!mEigenvalueDiagnosticsEnabled)
    return;

  if (!force && (mModelUpdateCounter % mEigenvalueDiagnosticsInterval) != 0)
    return;

  const std::vector<UInt> indices = getDiagnosticStateIndices();
  const UInt dimension = static_cast<UInt>(indices.size());
  if (dimension == 0)
    return;

  Matrix activeA = Matrix::Zero(dimension, dimension);
  for (UInt row = 0; row < dimension; ++row)
    for (UInt column = 0; column < dimension; ++column)
      activeA(row, column) = A(indices[row], indices[column]);

  if (!activeA.allFinite()) {
    **mDiagnosticsValid = 0.0;
    SPDLOG_LOGGER_ERROR(mSLog, "MMC diagnostic Jacobian contains NaN or Inf.");
    return;
  }

  Eigen::EigenSolver<Matrix> continuousSolver(activeA, false);
  if (continuousSolver.info() != Eigen::Success) {
    **mDiagnosticsValid = 0.0;
    SPDLOG_LOGGER_WARN(mSLog,
                       "MMC continuous-time eigenvalue calculation failed.");
    return;
  }

  const auto continuousEigenvalues = continuousSolver.eigenvalues();

  Real maximumReal = -std::numeric_limits<Real>::infinity();
  Real maximumMagnitude = 0.0;
  Real dominantFrequency = 0.0;

  for (Int index = 0; index < continuousEigenvalues.rows(); ++index) {
    const Real realPart = continuousEigenvalues(index).real();
    const Real imaginaryPart = continuousEigenvalues(index).imag();
    const Real magnitude = std::hypot(realPart, imaginaryPart);

    maximumMagnitude = std::max(maximumMagnitude, magnitude);

    if (realPart > maximumReal) {
      maximumReal = realPart;
      dominantFrequency = std::abs(imaginaryPart) / (2.0 * PI);
    }
  }

  **mJacobianMaximumRealEigenvalue = maximumReal;
  **mJacobianMaximumMagnitudeEigenvalue = maximumMagnitude;
  **mJacobianDominantFrequency = dominantFrequency;
  **mJacobianMaximumDiscreteMagnitude = 0.0;
  **mJacobianDiscreteDominantFrequency = 0.0;
  **mDiagnosticsValid = 1.0;

  Real maximumDiscreteMagnitude = 0.0;
  Real discreteDominantFrequency = 0.0;
  Bool discreteValid = false;

  if (mDiagnosticTimeStep > 0.0) {
    const Matrix identity = Matrix::Identity(dimension, dimension);
    const Matrix left = identity - 0.5 * mDiagnosticTimeStep * activeA;
    const Matrix right = identity + 0.5 * mDiagnosticTimeStep * activeA;

    const Matrix transition = left.colPivHouseholderQr().solve(right);

    if (transition.allFinite()) {
      Eigen::EigenSolver<Matrix> discreteSolver(transition, false);

      if (discreteSolver.info() == Eigen::Success) {
        const auto discreteEigenvalues = discreteSolver.eigenvalues();

        for (Int index = 0; index < discreteEigenvalues.rows(); ++index) {
          const Real realPart = discreteEigenvalues(index).real();
          const Real imaginaryPart = discreteEigenvalues(index).imag();
          const Real magnitude = std::hypot(realPart, imaginaryPart);

          if (magnitude > maximumDiscreteMagnitude) {
            maximumDiscreteMagnitude = magnitude;
            discreteDominantFrequency =
                std::abs(std::atan2(imaginaryPart, realPart)) /
                (2.0 * PI * mDiagnosticTimeStep);
          }
        }

        **mJacobianMaximumDiscreteMagnitude = maximumDiscreteMagnitude;
        **mJacobianDiscreteDominantFrequency = discreteDominantFrequency;
        discreteValid = true;

        if (force) {
          for (Int index = 0; index < discreteEigenvalues.rows(); ++index) {
            SPDLOG_LOGGER_DEBUG(mSLog,
                                "MMC Phi eigenvalue {}: {:.9e} {:+.9e}j, "
                                "|lambda|={:.9e}",
                                index, discreteEigenvalues(index).real(),
                                discreteEigenvalues(index).imag(),
                                std::abs(discreteEigenvalues(index)));
          }
        }
      }
    }
  }

  SPDLOG_LOGGER_INFO(mSLog,
                     "MMC local diagnostics: active_states={}, "
                     "max_real_A={:.9e} 1/s, "
                     "dominant_frequency_A={:.9e} Hz, "
                     "max_abs_A={:.9e} 1/s, "
                     "max_abs_Phi={:.9e}, "
                     "dominant_frequency_Phi={:.9e} Hz, "
                     "discrete_valid={}, dx_norm={:.9e}, "
                     "W_norm={:.9e}, yhist_norm={:.9e}",
                     dimension, maximumReal, dominantFrequency,
                     maximumMagnitude, maximumDiscreteMagnitude,
                     discreteDominantFrequency, discreteValid,
                     derivative.norm(), mW.norm(), mYHist.norm());

  if (force) {
    for (Int index = 0; index < continuousEigenvalues.rows(); ++index) {
      SPDLOG_LOGGER_DEBUG(mSLog, "MMC A eigenvalue {}: {:.9e} {:+.9e}j", index,
                          continuousEigenvalues(index).real(),
                          continuousEigenvalues(index).imag());
    }
  }
}

void EMT::Ph3::SSN_MMC::evaluateStateDerivative(const Matrix &x,
                                                const Matrix &u,
                                                Matrix &f) const {
  if (x.rows() != mStateSize || x.cols() != 1)
    throw std::invalid_argument("MMC state vector has invalid dimensions.");
  if (u.rows() != mInputSize || u.cols() != 1)
    throw std::invalid_argument("MMC input vector has invalid dimensions.");

  validateFinite(x, "MMC state");
  validateFinite(u, "MMC input");
  f.setZero(mStateSize, 1);

  // -----------------------------------------------------------------------
  // EMT/network frame and Harmony PLL frame
  // -----------------------------------------------------------------------
  //
  // The physical MMC plant states are kept in a nominal dq frame rotating at
  // omega_0. GridAngle is the corresponding absolute EMT angle. Harmony's PLL
  // state is an angle deviation relative to that nominal frame. Keeping these
  // two angles separate prevents the PLL from directly changing the frequency
  // of the current injected into the EMT network.
  const Real gridAngle = x(GridAngle, 0);
  const Real pllAngle = x(PllAngle, 0);

  const Matrix vAcAbc = u.block(0, 0, 3, 1);
  const Matrix vGridDqRaw = abcToDq(vAcAbc, gridAngle);
  const Matrix vControlDqRaw = rotateDq(vGridDqRaw, -pllAngle);

  const Real pllError = vControlDqRaw(1, 0);
  Real deltaOmega = 0.0;
  if (mPllEnabled) {
    f(XiPll, 0) = pllError;
    deltaOmega = mPllController.kp * pllError + mPllController.ki * x(XiPll, 0);
    f(PllAngle, 0) = deltaOmega;
  } else {
    f(XiPll, 0) = 0.0;
    f(PllAngle, 0) = 0.0;
  }

  f(GridAngle, 0) = mOmegaN;
  const Real omegaControl = mOmegaN + deltaOmega;

  // Plant states remain in the nominal network dq frame.
  const Real iDeltaD = x(IDeltaD, 0);
  const Real iDeltaQ = x(IDeltaQ, 0);
  const Real iSigmaZ = x(ISigmaZ, 0);
  const Real iSigmaD = x(ISigmaD, 0);
  const Real iSigmaQ = x(ISigmaQ, 0);
  const Real vCDeltaD = x(VCDeltaD, 0);
  const Real vCDeltaQ = x(VCDeltaQ, 0);
  const Real vCDeltaZd = x(VCDeltaZd, 0);
  const Real vCDeltaZq = x(VCDeltaZq, 0);
  const Real vCSigmaD = x(VCSigmaD, 0);
  const Real vCSigmaQ = x(VCSigmaQ, 0);
  const Real vCSigmaZ = x(VCSigmaZ, 0);

  const Real vGd = vGridDqRaw(0, 0);
  const Real vGq = vGridDqRaw(1, 0);

  // Active and reactive power are invariant under a common dq rotation, so
  // they are evaluated in the nominal network frame.
  const Real pAcRaw = 1.5 * (vGd * iDeltaD + vGq * iDeltaQ);
  const Real qAcRaw = 1.5 * (-vGd * iDeltaQ + vGq * iDeltaD);

  Real vControlD = applyFirstOrderFilter(x, f, FilterVacD, vControlDqRaw(0, 0),
                                         mAcVoltageDqFilterTimeConstant,
                                         mAcVoltageDqFilterTimeConstant > 0.0);
  Real vControlQ = applyFirstOrderFilter(x, f, FilterVacQ, vControlDqRaw(1, 0),
                                         mAcVoltageDqFilterTimeConstant,
                                         mAcVoltageDqFilterTimeConstant > 0.0);

  const Real vDcRaw = regularizedDcVoltage(u(Vdcp, 0) - u(Vdcn, 0));
  const Real vAcMagnitudeRaw = 1.5 * std::hypot(vControlD, vControlQ);

  const Real pAc = applySecondOrderFilter(x, f, FilterP1, FilterP2, pAcRaw,
                                          mActivePowerFilterTimeConstant,
                                          mActivePowerFilterTimeConstant > 0.0);
  const Real qAc = applySecondOrderFilter(
      x, f, FilterQ1, FilterQ2, qAcRaw, mReactivePowerFilterTimeConstant,
      mReactivePowerFilterTimeConstant > 0.0);
  const Real vDc = applySecondOrderFilter(x, f, FilterVdc1, FilterVdc2, vDcRaw,
                                          mDcVoltageFilterTimeConstant,
                                          mDcVoltageFilterTimeConstant > 0.0);
  const Real vAcMagnitude = applySecondOrderFilter(
      x, f, FilterVacMag1, FilterVacMag2, vAcMagnitudeRaw,
      mAcVoltageMagnitudeFilterTimeConstant,
      mAcVoltageMagnitudeFilterTimeConstant > 0.0);

  // -----------------------------------------------------------------------
  // Harmony outer loops: references are formed in the nominal network frame
  // and then rotated into the PLL controller frame.
  //
  // The current limit is applied to the dq-vector magnitude rather than to
  // each axis independently. Conditional integration prevents the outer-loop
  // PI states from winding up while the current command is saturated.
  // -----------------------------------------------------------------------
  auto conditionalIntegratorError = [](Real error, Real unsaturatedOutput,
                                       Real saturatedOutput,
                                       Real outputPerIntegratorSign) {
    const Real saturationDirection = unsaturatedOutput - saturatedOutput;
    const Real integratorDriveDirection = outputPerIntegratorSign * error;

    if (saturationDirection * integratorDriveDirection > 0.0)
      return 0.0;
    return error;
  };

  Real activeError = 0.0;
  Real activeIntegratorError = 0.0;
  Real iDeltaDReferenceRaw = mOpenLoopIDeltaDReference;

  switch (mActiveControlMode) {
  case ActiveControlMode::ActivePower:
    activeError = mActivePowerReference - pAc;
    iDeltaDReferenceRaw = mActiveController.kp * activeError +
                          mActiveController.ki * x(XiActive, 0);
    break;
  case ActiveControlMode::ActivePowerFeedforward:
    iDeltaDReferenceRaw = mActivePowerFeedforwardHeldIDeltaDReference;
    break;
  case ActiveControlMode::DcVoltage:
    // MATLAB VDC regulator: e_vdc = Vdc_measured - Vdc_reference.
    // At the rectifier, positive error increases iDelta_d (makes it less
    // negative), thereby reducing AC absorption and DC-side charging.
    activeError = vDc - mDcVoltageReference;
    iDeltaDReferenceRaw = mActiveController.kp * activeError +
                          mActiveController.ki * x(XiActive, 0);
    break;
  case ActiveControlMode::DcDroop:
    iDeltaDReferenceRaw =
        (mActivePowerReference + (mDcVoltageReference - vDc) / mDroopGain) /
        std::max(1.0, std::abs(vControlD));
    break;
  case ActiveControlMode::OpenLoop:
    break;
  }

  Real reactiveError = 0.0;
  Real reactiveIntegratorError = 0.0;
  Real reactiveIntegratorOutputSign = 1.0;
  Real iDeltaQReferenceRaw = mOpenLoopIDeltaQReference;

  switch (mReactiveControlMode) {
  case ReactiveControlMode::ReactivePower:
    reactiveError = mReactivePowerReference - qAc;
    reactiveIntegratorOutputSign = -1.0;
    iDeltaQReferenceRaw = -(mReactiveController.kp * reactiveError +
                            mReactiveController.ki * x(XiReactive, 0));
    break;
  case ReactiveControlMode::AcVoltage:
    reactiveError = mAcVoltageReference - vAcMagnitude;
    reactiveIntegratorOutputSign = 1.0;
    iDeltaQReferenceRaw = mReactiveController.kp * reactiveError +
                          mReactiveController.ki * x(XiReactive, 0);
    break;
  case ReactiveControlMode::OpenLoop:
    break;
  }

  Real iDeltaDReferenceNetwork = iDeltaDReferenceRaw;
  Real iDeltaQReferenceNetwork = iDeltaQReferenceRaw;

  const Real currentReferenceMagnitude =
      std::hypot(iDeltaDReferenceNetwork, iDeltaQReferenceNetwork);
  if (currentReferenceMagnitude > mMaximumAcCurrent) {
    const Real scale = mMaximumAcCurrent / currentReferenceMagnitude;
    iDeltaDReferenceNetwork *= scale;
    iDeltaQReferenceNetwork *= scale;
  }

  if ((mActiveControlMode == ActiveControlMode::ActivePower ||
       mActiveControlMode == ActiveControlMode::DcVoltage) &&
      mActiveController.ki > 0.0) {
    activeIntegratorError = conditionalIntegratorError(
        activeError, iDeltaDReferenceRaw, iDeltaDReferenceNetwork, 1.0);
  }

  if ((mReactiveControlMode == ReactiveControlMode::ReactivePower ||
       mReactiveControlMode == ReactiveControlMode::AcVoltage) &&
      mReactiveController.ki > 0.0) {
    reactiveIntegratorError = conditionalIntegratorError(
        reactiveError, iDeltaQReferenceRaw, iDeltaQReferenceNetwork,
        reactiveIntegratorOutputSign);
  }

  Matrix iDeltaNetwork(2, 1);
  iDeltaNetwork << iDeltaD, iDeltaQ;
  const Matrix iDeltaControl = rotateDq(iDeltaNetwork, -pllAngle);

  Matrix iDeltaReferenceNetwork(2, 1);
  iDeltaReferenceNetwork << iDeltaDReferenceNetwork, iDeltaQReferenceNetwork;
  const Matrix iDeltaReferenceControl =
      rotateDq(iDeltaReferenceNetwork, -pllAngle);

  const Real lEqAc = mArmInductance / 2.0 + mReactorInductance;
  const Real rEqAc = mArmResistance / 2.0 + mReactorResistance;

  // -----------------------------------------------------------------------
  // MATLAB output-current regulator in PLL coordinates:
  //   vMd = vGd + PI_d + R*iD_ref - omega*L*iQ_ref
  //   vMq = vGq + PI_q + R*iQ_ref + omega*L*iD_ref
  // with PI error i_ref - i_measured.
  // -----------------------------------------------------------------------
  const Real occErrorD = iDeltaReferenceControl(0, 0) - iDeltaControl(0, 0);
  const Real occErrorQ = iDeltaReferenceControl(1, 0) - iDeltaControl(1, 0);

  Matrix vMDeltaControl(2, 1);
  vMDeltaControl(0, 0) = vControlD + mOutputCurrentController.kp * occErrorD +
                         mOutputCurrentController.ki * x(XiOccD, 0) +
                         rEqAc * iDeltaReferenceControl(0, 0) -
                         omegaControl * lEqAc * iDeltaReferenceControl(1, 0);
  vMDeltaControl(1, 0) = vControlQ + mOutputCurrentController.kp * occErrorQ +
                         mOutputCurrentController.ki * x(XiOccQ, 0) +
                         rEqAc * iDeltaReferenceControl(1, 0) +
                         omegaControl * lEqAc * iDeltaReferenceControl(0, 0);

  // Controller command back to the nominal network dq frame.
  const Matrix vMDeltaNetwork = rotateDq(vMDeltaControl, pllAngle);
  **mInternalDifferentialVoltage = vMDeltaNetwork;
  Real vMDeltaDReference = vMDeltaNetwork(0, 0);
  Real vMDeltaQReference = vMDeltaNetwork(1, 0);
  if (mControlSource == ControlSource::ExternalDifferentialVoltage ||
      mControlSource == ControlSource::ExternalFullConverterVoltage) {
    if ((**mExternalDifferentialVoltage).rows() != 2 ||
        (**mExternalDifferentialVoltage).cols() != 1 ||
        !(**mExternalDifferentialVoltage).allFinite())
      throw std::runtime_error(
          "MMC external differential-voltage command is not finite 2x1.");
    vMDeltaDReference = (**mExternalDifferentialVoltage)(0, 0);
    vMDeltaQReference = (**mExternalDifferentialVoltage)(1, 0);
  }
  (**mAppliedDifferentialVoltage)(0, 0) = vMDeltaDReference;
  (**mAppliedDifferentialVoltage)(1, 0) = vMDeltaQReference;

  // -----------------------------------------------------------------------
  // Harmony circulating-current controller. The sigma dq quantities are
  // second-harmonic variables and therefore use twice the PLL deviation.
  // -----------------------------------------------------------------------
  Matrix iSigmaNetwork(2, 1);
  iSigmaNetwork << iSigmaD, iSigmaQ;
  const Matrix iSigmaControl = rotateDq(iSigmaNetwork, -2.0 * pllAngle);

  Matrix iSigmaReferenceNetwork(2, 1);
  iSigmaReferenceNetwork << mISigmaDReference, mISigmaQReference;
  const Matrix iSigmaReferenceControl =
      rotateDq(iSigmaReferenceNetwork, -2.0 * pllAngle);

  const Real cccErrorD = iSigmaReferenceControl(0, 0) - iSigmaControl(0, 0);
  const Real cccErrorQ = iSigmaReferenceControl(1, 0) - iSigmaControl(1, 0);

  Matrix vMSigmaControl(2, 1);
  vMSigmaControl(0, 0) =
      -(mCirculatingCurrentController.kp * cccErrorD +
        mCirculatingCurrentController.ki * x(XiCccD, 0) +
        2.0 * omegaControl * mArmInductance * iSigmaControl(1, 0));
  vMSigmaControl(1, 0) =
      -(mCirculatingCurrentController.kp * cccErrorQ +
        mCirculatingCurrentController.ki * x(XiCccQ, 0) -
        2.0 * omegaControl * mArmInductance * iSigmaControl(0, 0));

  const Matrix vMSigmaNetwork = rotateDq(vMSigmaControl, 2.0 * pllAngle);
  Real vMSigmaDReference = vMSigmaNetwork(0, 0);
  Real vMSigmaQReference = vMSigmaNetwork(1, 0);

  // -----------------------------------------------------------------------
  // Harmony stored-energy controller and zero-sequence current controller
  // -----------------------------------------------------------------------
  Real energyError = 0.0;
  Real energyIntegratorError = 0.0;
  Real iSigmaZReferenceRaw = mISigmaZReference;

  if (mEnergyControllerEnabled) {
    // Harmony total stored-energy reference.
    //
    // Use the nominal DC voltage. A DC-voltage disturbance must create an
    // energy error instead of moving the energy reference together with Vdc.
    const Real energyReference = 3.0 * mSubmoduleCapacitance *
                                 mNominalDcVoltage * mNominalDcVoltage /
                                 static_cast<Real>(mNumberOfSubmodules);

    energyError = energyReference - calculateStoredEnergy(x);

    // Harmony energy controller produces a power correction. Active power is
    // used as feedforward and the resulting DC-power command is converted into
    // the zero-sequence circulating-current reference.
    const Real energyPowerCommand = mEnergyController.kp * energyError +
                                    mEnergyController.ki * x(XiEnergy, 0) + pAc;

    iSigmaZReferenceRaw =
        energyPowerCommand / (3.0 * regularizedDcVoltage(vDc));
  }

  const Real iSigmaZReference =
      clamp(iSigmaZReferenceRaw, -mMaximumCirculatingCurrent,
            mMaximumCirculatingCurrent);

  if (mEnergyControllerEnabled && mEnergyController.ki > 0.0) {
    energyIntegratorError = conditionalIntegratorError(
        energyError, iSigmaZReferenceRaw, iSigmaZReference, 1.0);
  }

  const Real zccError = iSigmaZReference - iSigmaZ;
  Real vMSigmaZReference =
      -(mZeroSequenceCurrentController.kp * zccError +
        mZeroSequenceCurrentController.ki * x(XiZcc, 0) - vDc / 2.0);

  if (mControlSource == ControlSource::ExternalFullConverterVoltage) {
    if ((**mExternalCommonModeVoltage).rows() != 3 ||
        (**mExternalCommonModeVoltage).cols() != 1 ||
        !(**mExternalCommonModeVoltage).allFinite())
      throw std::runtime_error(
          "MMC external common-mode voltage command is not finite 3x1.");
    vMSigmaDReference = (**mExternalCommonModeVoltage)(0, 0);
    vMSigmaQReference = (**mExternalCommonModeVoltage)(1, 0);
    vMSigmaZReference = (**mExternalCommonModeVoltage)(2, 0);
  }

  (**mAppliedCommonModeVoltage)(0, 0) = vMSigmaDReference;
  (**mAppliedCommonModeVoltage)(1, 0) = vMSigmaQReference;
  (**mAppliedCommonModeVoltage)(2, 0) = vMSigmaZReference;

  const Real vDcReg = regularizedDcVoltage(vDc);
  Real mDeltaD = -2.0 * vMDeltaDReference / vDcReg;
  Real mDeltaQ = -2.0 * vMDeltaQReference / vDcReg;
  Real mDeltaZd = 0.0;
  Real mDeltaZq = 0.0;
  Real mSigmaD = 2.0 * vMSigmaDReference / vDcReg;
  Real mSigmaQ = 2.0 * vMSigmaQReference / vDcReg;
  Real mSigmaZ = 2.0 * vMSigmaZReference / vDcReg;

  Real occIntegratorErrorD = occErrorD;
  Real occIntegratorErrorQ = occErrorQ;

  const Real differentialMagnitude = std::hypot(mDeltaD, mDeltaQ);
  if (differentialMagnitude > mMaximumModulationMagnitude) {
    const Real scale = mMaximumModulationMagnitude / differentialMagnitude;

    // Vector anti-windup for the internal output-current controller. The
    // modulation limiter is radial, so freeze both OCC integrators only when
    // their combined drive would push the voltage command farther outside the
    // feasible modulation circle.
    if (mControlSource == ControlSource::InternalControllers) {
      Matrix saturationExcessNetwork(2, 1);
      saturationExcessNetwork << (1.0 - scale) * vMDeltaDReference,
          (1.0 - scale) * vMDeltaQReference;

      Matrix integratorDriveControl(2, 1);
      integratorDriveControl << mOutputCurrentController.ki * occErrorD,
          mOutputCurrentController.ki * occErrorQ;
      const Matrix integratorDriveNetwork =
          rotateDq(integratorDriveControl, pllAngle);

      if (saturationExcessNetwork.cwiseProduct(integratorDriveNetwork).sum() >
          0.0) {
        occIntegratorErrorD = 0.0;
        occIntegratorErrorQ = 0.0;
      }
    }

    mDeltaD *= scale;
    mDeltaQ *= scale;
  }

  const Real sigmaMagnitude = std::hypot(mSigmaD, mSigmaQ);
  if (sigmaMagnitude > mMaximumModulationMagnitude) {
    const Real scale = mMaximumModulationMagnitude / sigmaMagnitude;
    mSigmaD *= scale;
    mSigmaQ *= scale;
  }
  mSigmaZ =
      clamp(mSigmaZ, -mMaximumModulationMagnitude, mMaximumModulationMagnitude);

  mDeltaD = applyPadeDelayChannel(x, f, DelayMDeltaD1, DelayMDeltaD2, mDeltaD);
  mDeltaQ = applyPadeDelayChannel(x, f, DelayMDeltaQ1, DelayMDeltaQ2, mDeltaQ);
  mSigmaD = applyPadeDelayChannel(x, f, DelayMSigmaD1, DelayMSigmaD2, mSigmaD);
  mSigmaQ = applyPadeDelayChannel(x, f, DelayMSigmaQ1, DelayMSigmaQ2, mSigmaQ);
  mSigmaZ = applyPadeDelayChannel(x, f, DelayMSigmaZ1, DelayMSigmaZ2, mSigmaZ);
  (**mAppliedModulation)(0, 0) = mDeltaD;
  (**mAppliedModulation)(1, 0) = mDeltaQ;
  (**mAppliedModulation)(2, 0) = mSigmaD;
  (**mAppliedModulation)(3, 0) = mSigmaQ;
  (**mAppliedModulation)(4, 0) = mSigmaZ;

  // Harmony averaged MMC modulation-voltage equations.
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
  (**mRealizedConverterVoltage)(0, 0) = vMDeltaD;
  (**mRealizedConverterVoltage)(1, 0) = vMDeltaQ;
  (**mRealizedConverterVoltage)(2, 0) = vMSigmaD;
  (**mRealizedConverterVoltage)(3, 0) = vMSigmaQ;
  (**mRealizedConverterVoltage)(4, 0) = vMSigmaZ;

  // The electrical plant states use the same q=-sin(theta) Park convention
  // as the MATLAB controller. iDelta is positive from converter to AC grid:
  //   vMd = vGd + R*iD + L*diD/dt - omega*L*iQ
  //   vMq = vGq + R*iQ + L*diQ/dt + omega*L*iD
  f(IDeltaD, 0) =
      (vMDeltaD - vGd - rEqAc * iDeltaD + lEqAc * iDeltaQ * mOmegaN) / lEqAc;
  f(IDeltaQ, 0) =
      (vMDeltaQ - vGq - rEqAc * iDeltaQ - lEqAc * iDeltaD * mOmegaN) / lEqAc;
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

  f(XiActive, 0) = (mActiveControlMode == ActiveControlMode::ActivePower ||
                    mActiveControlMode == ActiveControlMode::DcVoltage)
                       ? activeIntegratorError
                       : 0.0;
  f(XiReactive, 0) =
      (mReactiveControlMode == ReactiveControlMode::ReactivePower ||
       mReactiveControlMode == ReactiveControlMode::AcVoltage)
          ? reactiveIntegratorError
          : 0.0;
  f(XiOccD, 0) = mControlSource == ControlSource::InternalControllers
                     ? occIntegratorErrorD
                     : 0.0;
  f(XiOccQ, 0) = mControlSource == ControlSource::InternalControllers
                     ? occIntegratorErrorQ
                     : 0.0;
  const Bool internalCommonModeControl =
      mControlSource != ControlSource::ExternalFullConverterVoltage;
  f(XiCccD, 0) = internalCommonModeControl ? cccErrorD : 0.0;
  f(XiCccQ, 0) = internalCommonModeControl ? cccErrorQ : 0.0;
  f(XiZcc, 0) = internalCommonModeControl ? zccError : 0.0;
  f(XiEnergy, 0) = mEnergyControllerEnabled ? energyIntegratorError : 0.0;

  validateFinite(f, "MMC state derivative");
}

void EMT::Ph3::SSN_MMC::evaluateOutput(const Matrix &x, const Matrix &u,
                                       Matrix &output) const {
  if (x.rows() != mStateSize || x.cols() != 1)
    throw std::invalid_argument("MMC state vector has invalid dimensions.");
  if (u.rows() != mInputSize || u.cols() != 1)
    throw std::invalid_argument("MMC input vector has invalid dimensions.");

  validateFinite(x, "MMC output state");
  validateFinite(u, "MMC output input");
  output.setZero(mOutputSize, 1);

  // MATLAB iDelta is positive from converter to AC grid. DPsim's interface
  // current is positive from the MNA node into the component.
  const Matrix iAcAbc = dqToAbc(x(IDeltaD, 0), x(IDeltaQ, 0), x(GridAngle, 0));
  output.block(0, 0, 3, 1) = -iAcAbc;

  // Positive iSigmaZ is current absorbed from the DC network: current
  // enters the converter at dc+ and leaves it at dc-. DPsim interface current
  // is positive from the MNA node into the component.
  const Real iDc = 3.0 * x(ISigmaZ, 0);
  output(Idcp, 0) = iDc;
  output(Idcn, 0) = -iDc;
  validateFinite(output, "MMC output");
}

void EMT::Ph3::SSN_MMC::calculateNumericalJacobians(const Matrix &x,
                                                    const Matrix &u, Matrix &A,
                                                    Matrix &B, Matrix &C,
                                                    Matrix &D) const {
  A.setZero(mStateSize, mStateSize);
  B.setZero(mStateSize, mInputSize);
  C.setZero(mOutputSize, mStateSize);
  D.setZero(mOutputSize, mInputSize);

  Matrix fPlus = Matrix::Zero(mStateSize, 1);
  Matrix fMinus = Matrix::Zero(mStateSize, 1);
  Matrix gPlus = Matrix::Zero(mOutputSize, 1);
  Matrix gMinus = Matrix::Zero(mOutputSize, 1);

  for (UInt column = 0; column < mStateSize; ++column) {
    const Real step =
        mJacobianAbsoluteStep +
        mJacobianRelativeStep * std::max(1.0, std::abs(x(column, 0)));
    Matrix xPlus = x;
    Matrix xMinus = x;
    xPlus(column, 0) += step;
    xMinus(column, 0) -= step;
    evaluateStateDerivative(xPlus, u, fPlus);
    evaluateStateDerivative(xMinus, u, fMinus);
    evaluateOutput(xPlus, u, gPlus);
    evaluateOutput(xMinus, u, gMinus);
    A.col(column) = (fPlus - fMinus) / (2.0 * step);
    C.col(column) = (gPlus - gMinus) / (2.0 * step);
  }

  for (UInt column = 0; column < mInputSize; ++column) {
    const Real step =
        mJacobianAbsoluteStep +
        mJacobianRelativeStep * std::max(1.0, std::abs(u(column, 0)));
    Matrix uPlus = u;
    Matrix uMinus = u;
    uPlus(column, 0) += step;
    uMinus(column, 0) -= step;
    evaluateStateDerivative(x, uPlus, fPlus);
    evaluateStateDerivative(x, uMinus, fMinus);
    evaluateOutput(x, uPlus, gPlus);
    evaluateOutput(x, uMinus, gMinus);
    B.col(column) = (fPlus - fMinus) / (2.0 * step);
    D.col(column) = (gPlus - gMinus) / (2.0 * step);
  }

  validateFinite(A, "MMC Jacobian A");
  validateFinite(B, "MMC Jacobian B");
  validateFinite(C, "MMC Jacobian C");
  validateFinite(D, "MMC Jacobian D");
}

void EMT::Ph3::SSN_MMC::buildStateSpaceModel(const Matrix &x, const Matrix &u,
                                             Matrix &A, Matrix &B, Matrix &C,
                                             Matrix &D, Matrix &E,
                                             Matrix &F) const {
  calculateNumericalJacobians(x, u, A, B, C, D);
  Matrix derivative = Matrix::Zero(mStateSize, 1);
  Matrix output = Matrix::Zero(mOutputSize, 1);
  evaluateStateDerivative(x, u, derivative);
  evaluateOutput(x, u, output);
  E = derivative - A * x - B * u;
  F = output - C * x - D * u;
  validateFinite(E, "MMC affine state offset");
  validateFinite(F, "MMC affine output offset");
}

Bool EMT::Ph3::SSN_MMC::updateComponentParameters() {
  advanceSampledActivePowerFeedforward();

  Matrix stateOffset;
  Matrix outputOffset;
  buildStateSpaceModel(**mX, **mIntfVoltage, mA, mB, mC, mD, stateOffset,
                       outputOffset);
  setStateOffset(stateOffset);
  setOutputOffset(outputOffset);
  updateDiagnostics(mA, **mX, **mIntfVoltage, false);
  return true;
}

void EMT::Ph3::SSN_MMC::addHeldControlDependencies(
    AttributeBase::List &prevStepDependencies) const {
  prevStepDependencies.push_back(mExternalDifferentialVoltage);
  prevStepDependencies.push_back(mExternalCommonModeVoltage);
}

void EMT::Ph3::SSN_MMC::initializeFromNodesAndTerminals(Real frequency) {
  validateTerminalArrangement();

  if (!mParametersSet)
    throw std::logic_error(
        "setParameters() must be called before MMC initialization.");

  if (std::abs(frequency - mNominalFrequency) >
      1e-6 * std::max(1.0, mNominalFrequency)) {
    SPDLOG_LOGGER_WARN(
        mSLog,
        "Initialization frequency {:.9g} Hz differs from MMC nominal "
        "frequency {:.9g} Hz.",
        frequency, mNominalFrequency);
  }

  const MatrixComp uPhasor = buildInitialInputFromNodes(frequency);
  if (uPhasor.rows() != mInputSize || uPhasor.cols() != 1)
    throw std::logic_error("MMC initial input must be 5x1.");

  Matrix u0 = uPhasor.real();
  if (std::abs(u0(Vdcp, 0) - u0(Vdcn, 0)) < mMinimumDcVoltage) {
    // Respect grounded-pole topology instead of inventing a symmetric bipolar
    // voltage when one physical terminal is grounded.
    if (!terminalNotGrounded(2)) {
      u0(Vdcp, 0) = mNominalDcVoltage;
      u0(Vdcn, 0) = 0.0;
    } else if (!terminalNotGrounded(1)) {
      u0(Vdcp, 0) = 0.0;
      u0(Vdcn, 0) = -mNominalDcVoltage;
    } else {
      u0(Vdcp, 0) = 0.5 * mNominalDcVoltage;
      u0(Vdcn, 0) = -0.5 * mNominalDcVoltage;
    }
  }
  validateFinite(u0, "MMC initial input");

  Matrix x0 = Matrix::Zero(mStateSize, 1);
  initializeAnalyticalOperatingPoint(x0, u0);
  initializeSampledActivePowerFeedforward(x0, u0);

  Real normalizedResidual = 0.0;
  Real absoluteResidual = 0.0;
  Bool converged = true;
  if (mOperatingPointInitializationEnabled) {
    converged =
        solveOperatingPoint(x0, u0, normalizedResidual, absoluteResidual);
  } else {
    Matrix derivative = Matrix::Zero(mStateSize, 1);
    evaluateStateDerivative(x0, u0, derivative);
    derivative(PllAngle, 0) = 0.0;
    derivative(GridAngle, 0) = 0.0;
    absoluteResidual = derivative.norm();
    normalizedResidual = absoluteResidual;
  }

  validateFinite(x0, "MMC initial state");
  **mX = x0;
  **mIntfVoltage = u0;
  **mEquilibriumResidualNorm = absoluteResidual;

  updateStateSpaceModel();
  mYHist = calculateHistoryVector();
  validateFinite(mYHist, "MMC history vector");

  **mIntfCurrent = mW * (**mIntfVoltage) + mYHist;
  validateFinite(**mIntfCurrent, "MMC initial interface current");

  updateLogAttributes(**mIntfVoltage);
  updateDiagnostics(mA, **mX, **mIntfVoltage, true);

  Matrix derivative = Matrix::Zero(mStateSize, 1);
  Matrix nonlinearOutput = Matrix::Zero(mOutputSize, 1);
  evaluateStateDerivative(**mX, **mIntfVoltage, derivative);
  evaluateOutput(**mX, **mIntfVoltage, nonlinearOutput);
  derivative(PllAngle, 0) = 0.0;
  derivative(GridAngle, 0) = 0.0;
  const Matrix ssnOutput = mW * (**mIntfVoltage) + mYHist;

  if (!converged) {
    SPDLOG_LOGGER_WARN(
        mSLog,
        "MMC operating-point solve did not reach the requested normalized "
        "tolerance: normalized_residual={:.9e}, absolute_residual={:.9e}.",
        normalizedResidual, absoluteResidual);
  }

  SPDLOG_LOGGER_INFO(
      mSLog,
      "\n--- SSN MMC initialization ---"
      "\nInput u: {:s}"
      "\nState x: {:s}"
      "\nOperating-point converged: {}"
      "\nNormalized equilibrium residual: {:.9e}"
      "\nAbsolute equilibrium residual: {:.9e}"
      "\nElectrical derivative norm (angle-clock states excluded): {:.9e}"
      "\nNonlinear output: {:s}"
      "\nSSN output: {:s}"
      "\nOutput mismatch norm: {:.9e}"
      "\nLocal max-real eigenvalue: {:.9e} 1/s"
      "\nLocal dominant frequency: {:.9e} Hz"
      "\n--- SSN MMC initialization finished ---",
      Logger::matrixToString(**mIntfVoltage), Logger::matrixToString(**mX),
      converged, normalizedResidual, absoluteResidual, derivative.norm(),
      Logger::matrixToString(nonlinearOutput),
      Logger::matrixToString(ssnOutput), (nonlinearOutput - ssnOutput).norm(),
      **mJacobianMaximumRealEigenvalue, **mJacobianDominantFrequency);

  if (mActiveControlMode == ActiveControlMode::ActivePowerFeedforward) {
    mActivePowerFeedforwardStepCounter = 0;
    mActivePowerFeedforwardRuntimeEnabled = true;
  }
}

void EMT::Ph3::SSN_MMC::updateLogAttributes(const Matrix &u) const {
  validateFinite(**mX, "MMC logged state");
  validateFinite(u, "MMC logged input");

  const Real gridAngle = (**mX)(GridAngle, 0);
  const Real pllAngle = (**mX)(PllAngle, 0);
  const Matrix vGridDq = abcToDq(u.block(0, 0, 3, 1), gridAngle);
  const Matrix vControlDq = rotateDq(vGridDq, -pllAngle);

  const Real iD = (**mX)(IDeltaD, 0);
  const Real iQ = (**mX)(IDeltaQ, 0);
  const Real p = 1.5 * (vGridDq(0, 0) * iD + vGridDq(1, 0) * iQ);
  const Real q = 1.5 * (-vGridDq(0, 0) * iQ + vGridDq(1, 0) * iD);

  **mDcPositiveVoltage = u(Vdcp, 0);
  **mDcNegativeVoltage = u(Vdcn, 0);
  **mDcVoltage = u(Vdcp, 0) - u(Vdcn, 0);
  **mDcCurrent = 3.0 * (**mX)(ISigmaZ, 0);
  **mActivePower = p;
  **mReactivePower = q;
  **mAcVoltageMagnitude = std::hypot(vGridDq(0, 0), vGridDq(1, 0));
  **mStoredEnergy = calculateStoredEnergy(**mX);
  **mAcTerminalVoltage = u.block(0, 0, 3, 1);
  **mAcTerminalCurrent = (**mIntfCurrent).block(0, 0, 3, 1);

  const Real pllError = vControlDq(1, 0);
  const Real deltaOmega = mPllEnabled ? mPllController.kp * pllError +
                                            mPllController.ki * (**mX)(XiPll, 0)
                                      : 0.0;

  **mGridAngle = std::remainder(gridAngle, 2.0 * PI);
  **mPllAngleDeviation = std::remainder(pllAngle, 2.0 * PI);
  **mConverterAngle = std::remainder(gridAngle + pllAngle, 2.0 * PI);
  **mPllFrequency = (mOmegaN + deltaOmega) / (2.0 * PI);
  **mPllError = pllError;

  **mGridVoltageD = vGridDq(0, 0);
  **mGridVoltageQ = vGridDq(1, 0);
  **mControlVoltageD = vControlDq(0, 0);
  **mControlVoltageQ = vControlDq(1, 0);
  **mDeltaCurrentD = iD;
  **mDeltaCurrentQ = iQ;
  **mSigmaCurrentZ = (**mX)(ISigmaZ, 0);

  const Real pControl =
      mActivePowerFilterTimeConstant > 0.0 ? (**mX)(FilterP2, 0) : p;
  const Real qControl =
      mReactivePowerFilterTimeConstant > 0.0 ? (**mX)(FilterQ2, 0) : q;
  const Real vDcControl = mDcVoltageFilterTimeConstant > 0.0
                              ? (**mX)(FilterVdc2, 0)
                              : (**mDcVoltage);
  const Real vControlDForOuterLoop = mAcVoltageDqFilterTimeConstant > 0.0
                                         ? (**mX)(FilterVacD, 0)
                                         : vControlDq(0, 0);
  const Real vAcMagnitudeForOuterLoop =
      mAcVoltageMagnitudeFilterTimeConstant > 0.0
          ? (**mX)(FilterVacMag2, 0)
          : std::hypot(vControlDq(0, 0), vControlDq(1, 0));

  Real iDeltaDReference = mOpenLoopIDeltaDReference;
  switch (mActiveControlMode) {
  case ActiveControlMode::ActivePower:
    iDeltaDReference =
        mActiveController.kp * (mActivePowerReference - pControl) +
        mActiveController.ki * (**mX)(XiActive, 0);
    break;
  case ActiveControlMode::ActivePowerFeedforward:
    iDeltaDReference = mActivePowerFeedforwardHeldIDeltaDReference;
    break;
  case ActiveControlMode::DcVoltage:
    iDeltaDReference =
        mActiveController.kp * (vDcControl - mDcVoltageReference) +
        mActiveController.ki * (**mX)(XiActive, 0);
    break;
  case ActiveControlMode::DcDroop:
    iDeltaDReference = (mActivePowerReference +
                        (mDcVoltageReference - vDcControl) / mDroopGain) /
                       std::max(1.0, std::abs(vControlDForOuterLoop));
    break;
  case ActiveControlMode::OpenLoop:
    break;
  }
  Real iDeltaQReference = mOpenLoopIDeltaQReference;
  switch (mReactiveControlMode) {
  case ReactiveControlMode::ReactivePower:
    iDeltaQReference =
        -(mReactiveController.kp * (mReactivePowerReference - qControl) +
          mReactiveController.ki * (**mX)(XiReactive, 0));
    break;
  case ReactiveControlMode::AcVoltage:
    iDeltaQReference = mReactiveController.kp *
                           (mAcVoltageReference - vAcMagnitudeForOuterLoop) +
                       mReactiveController.ki * (**mX)(XiReactive, 0);
    break;
  case ReactiveControlMode::OpenLoop:
    break;
  }
  const Real diagnosticCurrentReferenceMagnitude =
      std::hypot(iDeltaDReference, iDeltaQReference);
  if (diagnosticCurrentReferenceMagnitude > mMaximumAcCurrent) {
    const Real scale = mMaximumAcCurrent / diagnosticCurrentReferenceMagnitude;
    iDeltaDReference *= scale;
    iDeltaQReference *= scale;
  }

  Real iSigmaZReference = mISigmaZReference;
  if (mEnergyControllerEnabled) {
    const Real energyReference = 3.0 * mSubmoduleCapacitance *
                                 mNominalDcVoltage * mNominalDcVoltage /
                                 static_cast<Real>(mNumberOfSubmodules);
    const Real energyError = energyReference - calculateStoredEnergy(**mX);
    const Real energyPowerCommand = mEnergyController.kp * energyError +
                                    mEnergyController.ki * (**mX)(XiEnergy, 0) +
                                    pControl;
    iSigmaZReference =
        energyPowerCommand / (3.0 * regularizedDcVoltage(vDcControl));
  }
  iSigmaZReference = clamp(iSigmaZReference, -mMaximumCirculatingCurrent,
                           mMaximumCirculatingCurrent);

  **mDeltaCurrentReferenceD = iDeltaDReference;
  **mDeltaCurrentReferenceQ = iDeltaQReference;
  **mSigmaCurrentReferenceZ = iSigmaZReference;

  const Real pDc = (**mDcVoltage) * (**mDcCurrent);
  **mDcPower = pDc;

  // p is AC power exported by the converter. pDc is DC power absorbed by
  // the converter because mDcCurrent = 3*iSigmaZ is positive from dc+ into
  // the converter. At a stationary operating point:
  //   p - pDc + represented copper losses = 0.
  // During transients the stored-energy derivative is intentionally omitted.
  const Real differentialLoss =
      1.5 * (mArmResistance / 2.0 + mReactorResistance) * (iD * iD + iQ * iQ);
  const Real circulatingLoss =
      3.0 * mArmResistance *
          ((**mX)(ISigmaD, 0) * (**mX)(ISigmaD, 0) +
           (**mX)(ISigmaQ, 0) * (**mX)(ISigmaQ, 0)) +
      6.0 * mArmResistance * (**mX)(ISigmaZ, 0) * (**mX)(ISigmaZ, 0);
  **mPowerBalanceError = p - pDc + differentialLoss + circulatingLoss;

  **mFilteredActivePower =
      mActivePowerFilterTimeConstant > 0.0 ? (**mX)(FilterP2, 0) : p;
  **mFilteredReactivePower =
      mReactivePowerFilterTimeConstant > 0.0 ? (**mX)(FilterQ2, 0) : q;
  **mFilteredDcVoltage = mDcVoltageFilterTimeConstant > 0.0
                             ? (**mX)(FilterVdc2, 0)
                             : u(Vdcp, 0) - u(Vdcn, 0);
  **mFilteredDaxisVoltage =
      mActiveControlMode == ActiveControlMode::ActivePowerFeedforward
          ? mActivePowerFeedforwardFilteredVd
          : 0.0;

  **mHeldActiveCurrentReference =
      mActiveControlMode == ActiveControlMode::ActivePowerFeedforward
          ? mActivePowerFeedforwardHeldIDeltaDReference
          : 0.0;

  Matrix derivative = Matrix::Zero(mStateSize, 1);
  evaluateStateDerivative(**mX, u, derivative);
  **mStateNorm = (**mX).norm();
  **mStateDerivativeNorm = derivative.norm();
  **mNortonMatrixNorm = mW.norm();
  **mHistoryVectorNorm = mYHist.norm();
}

Matrix EMT::Ph3::SSN_MMC::getState() const { return **mX; }

Matrix EMT::Ph3::SSN_MMC::getStateDerivative() const {
  Matrix derivative = Matrix::Zero(mStateSize, 1);
  evaluateStateDerivative(**mX, **mIntfVoltage, derivative);
  return derivative;
}

Matrix EMT::Ph3::SSN_MMC::getInterfaceVoltage() const { return **mIntfVoltage; }

Matrix EMT::Ph3::SSN_MMC::getInterfaceCurrent() const { return **mIntfCurrent; }

void EMT::Ph3::SSN_MMC::getLocalLinearization(Matrix &A, Matrix &B, Matrix &C,
                                              Matrix &D) const {
  calculateNumericalJacobians(**mX, **mIntfVoltage, A, B, C, D);
}

Attribute<Matrix>::Ptr EMT::Ph3::SSN_MMC::acTerminalVoltageAttribute() const {
  return mAcTerminalVoltage;
}

Attribute<Matrix>::Ptr EMT::Ph3::SSN_MMC::acTerminalCurrentAttribute() const {
  return mAcTerminalCurrent;
}

Attribute<Matrix>::Ptr EMT::Ph3::SSN_MMC::interfaceVoltageAttribute() const {
  return mIntfVoltage;
}

Attribute<Matrix>::Ptr EMT::Ph3::SSN_MMC::interfaceCurrentAttribute() const {
  return mIntfCurrent;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::dcPositiveVoltageAttribute() const {
  return mDcPositiveVoltage;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::dcNegativeVoltageAttribute() const {
  return mDcNegativeVoltage;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::dcVoltageAttribute() const {
  return mDcVoltage;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::dcCurrentAttribute() const {
  return mDcCurrent;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::activePowerAttribute() const {
  return mActivePower;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::reactivePowerAttribute() const {
  return mReactivePower;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::storedEnergyAttribute() const {
  return mStoredEnergy;
}

Attribute<Real>::Ptr EMT::Ph3::SSN_MMC::filteredDaxisVoltageAttribute() const {
  return mFilteredDaxisVoltage;
}

Attribute<Real>::Ptr
EMT::Ph3::SSN_MMC::heldActiveCurrentReferenceAttribute() const {
  return mHeldActiveCurrentReference;
}

Attribute<Matrix>::Ptr
EMT::Ph3::SSN_MMC::externalDifferentialVoltageAttribute() const {
  return mExternalDifferentialVoltage;
}

Attribute<Matrix>::Ptr
EMT::Ph3::SSN_MMC::externalCommonModeVoltageAttribute() const {
  return mExternalCommonModeVoltage;
}

Attribute<Matrix>::Ptr
EMT::Ph3::SSN_MMC::appliedDifferentialVoltageAttribute() const {
  return mAppliedDifferentialVoltage;
}

Attribute<Matrix>::Ptr
EMT::Ph3::SSN_MMC::internalDifferentialVoltageAttribute() const {
  return mInternalDifferentialVoltage;
}

Attribute<Matrix>::Ptr
EMT::Ph3::SSN_MMC::appliedCommonModeVoltageAttribute() const {
  return mAppliedCommonModeVoltage;
}

Attribute<Matrix>::Ptr EMT::Ph3::SSN_MMC::appliedModulationAttribute() const {
  return mAppliedModulation;
}

Attribute<Matrix>::Ptr
EMT::Ph3::SSN_MMC::realizedConverterVoltageAttribute() const {
  return mRealizedConverterVoltage;
}

Attribute<Bool>::Ptr EMT::Ph3::SSN_MMC::externalCommandActiveAttribute() const {
  return mExternalCommandActive;
}
