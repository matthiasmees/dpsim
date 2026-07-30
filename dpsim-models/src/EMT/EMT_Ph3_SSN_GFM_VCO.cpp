// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <cmath>
#include <stdexcept>

#include <dpsim-models/EMT/EMT_Ph3_SSN_GFM_VCO.h>

using namespace CPS;

EMT::Ph3::SSN_GFM_VCO::SSN_GFM_VCO(String uid, String name,
                                   Logger::Level logLevel)
    : TwoTerminalVTypeVariableSSNComp(uid, name, logLevel), mLf(0.0), mCf(0.0),
      mRf(0.0), mRc(0.0), mSystemOmega(0.0), mOmegaN(0.0), mVdRef(0.0),
      mVqRef(0.0), mKpVoltage(0.0), mKiVoltage(0.0), mKpCurrent(0.0),
      mKiCurrent(0.0), mGeneralParametersSet(false),
      mControllerParametersSet(false), mFilterParametersSet(false),
      mInitialStateOverrideEnabled(false), mInitialVoltageIntegratorD(0.0),
      mInitialVoltageIntegratorQ(0.0), mInitialCurrentIntegratorD(0.0),
      mInitialCurrentIntegratorQ(0.0), mInitialPowerOverrideEnabled(false),
      mInitialActivePower(0.0), mInitialReactivePower(0.0),
      mLinearizationUpdateInterval(1), mStepsSinceLinearization(0),
      mLinearizationInitialized(false),
      mPInst(mAttributes->create<Real>("p_inst")),
      mQInst(mAttributes->create<Real>("q_inst")),
      mOmegaGFM(mAttributes->create<Real>("omega_gfm")),
      mThetaGFM(mAttributes->create<Real>("theta_gfm")),
      mVoltageMagnitudeGFM(mAttributes->create<Real>("voltage_magnitude_gfm")),
      mVcD(mAttributes->create<Real>("vc_d")),
      mVcQ(mAttributes->create<Real>("vc_q")),
      mIGridD(mAttributes->create<Real>("i_grid_d")),
      mIGridQ(mAttributes->create<Real>("i_grid_q")),
      mIfD(mAttributes->create<Real>("if_d")),
      mIfQ(mAttributes->create<Real>("if_q")),
      mVoltageReferenceD(mAttributes->create<Real>("v_ref_d")),
      mVoltageReferenceQ(mAttributes->create<Real>("v_ref_q")),
      mConverterVoltageD(mAttributes->create<Real>("v_converter_d")),
      mConverterVoltageQ(mAttributes->create<Real>("v_converter_q")) {

  **mIntfVoltage = Matrix::Zero(mInputSize, 1);
  **mIntfCurrent = Matrix::Zero(mOutputSize, 1);

  **mPInst = 0.0;
  **mQInst = 0.0;
  **mOmegaGFM = 0.0;
  **mThetaGFM = 0.0;
  **mVoltageMagnitudeGFM = 0.0;
  **mVcD = 0.0;
  **mVcQ = 0.0;
  **mIGridD = 0.0;
  **mIGridQ = 0.0;
  **mIfD = 0.0;
  **mIfQ = 0.0;
  **mVoltageReferenceD = 0.0;
  **mVoltageReferenceQ = 0.0;
  **mConverterVoltageD = 0.0;
  **mConverterVoltageQ = 0.0;
}

std::vector<String> EMT::Ph3::SSN_GFM_VCO::getLocalStateNames() const {
  return {
      "theta_vco",
      "voltage_integrator_d",
      "voltage_integrator_q",
      "current_integrator_d",
      "current_integrator_q",
      "vc_a",
      "vc_b",
      "vc_c",
      "if_a",
      "if_b",
      "if_c",
  };
}

std::vector<EMT::SSNComp::LocalAbcStateBlock>
EMT::Ph3::SSN_GFM_VCO::getLocalAbcStateBlocks() const {
  return {
      {{static_cast<Int>(VcA), static_cast<Int>(VcB), static_cast<Int>(VcC)},
       "vc"},
      {{static_cast<Int>(IfA), static_cast<Int>(IfB), static_cast<Int>(IfC)},
       "if"},
  };
}

void EMT::Ph3::SSN_GFM_VCO::markLinearizationDirty() {
  mLinearizationInitialized = false;
  mStepsSinceLinearization = 0;
}

void EMT::Ph3::SSN_GFM_VCO::configureStateSpaceDimensionsIfReady() {
  if (!mGeneralParametersSet || !mControllerParametersSet ||
      !mFilterParametersSet)
    return;

  VTypeVariableSSNComp::setParameters(Matrix::Zero(mStateSize, mStateSize),
                                      Matrix::Zero(mStateSize, mInputSize),
                                      Matrix::Zero(mOutputSize, mStateSize),
                                      Matrix::Zero(mOutputSize, mInputSize),
                                      Matrix::Zero(mStateSize, 1),
                                      Matrix::Zero(mOutputSize, 1));
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::setParameters(Real sysOmega, Real vdRef,
                                          Real vqRef) {
  if (!std::isfinite(sysOmega) || sysOmega <= 0.0)
    throw std::invalid_argument(
        "Nominal angular frequency must be finite and positive.");
  if (!std::isfinite(vdRef) || !std::isfinite(vqRef))
    throw std::invalid_argument("Voltage references must be finite.");

  mSystemOmega = sysOmega;
  if (!mControllerParametersSet)
    mOmegaN = sysOmega;

  mVdRef = vdRef;
  mVqRef = vqRef;
  mGeneralParametersSet = true;

  configureStateSpaceDimensionsIfReady();
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::setControllerParameters(Real kpVoltage,
                                                    Real kiVoltage,
                                                    Real kpCurrent,
                                                    Real kiCurrent,
                                                    Real omegaNominal) {
  if (!std::isfinite(kpVoltage) || !std::isfinite(kiVoltage) ||
      !std::isfinite(kpCurrent) || !std::isfinite(kiCurrent))
    throw std::invalid_argument("Controller gains must be finite.");
  if (kiVoltage == 0.0)
    throw std::invalid_argument(
        "Voltage-controller integral gain must be non-zero.");
  if (kiCurrent == 0.0)
    throw std::invalid_argument(
        "Current-controller integral gain must be non-zero.");
  if (!std::isfinite(omegaNominal) || omegaNominal <= 0.0)
    throw std::invalid_argument(
        "VCO angular frequency must be finite and positive.");

  mKpVoltage = kpVoltage;
  mKiVoltage = kiVoltage;
  mKpCurrent = kpCurrent;
  mKiCurrent = kiCurrent;
  mOmegaN = omegaNominal;
  mControllerParametersSet = true;

  configureStateSpaceDimensionsIfReady();
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::setFilterParameters(Real lf, Real cf, Real rf,
                                                Real rc) {
  if (!std::isfinite(lf) || lf <= 0.0)
    throw std::invalid_argument("Filter inductance lf must be positive.");
  if (!std::isfinite(cf) || cf <= 0.0)
    throw std::invalid_argument("Filter capacitance cf must be positive.");
  if (!std::isfinite(rf) || rf < 0.0)
    throw std::invalid_argument(
        "Filter resistance rf must be finite and non-negative.");
  if (!std::isfinite(rc) || rc <= 0.0)
    throw std::invalid_argument("Coupling resistance rc must be positive.");

  mLf = lf;
  mCf = cf;
  mRf = rf;
  mRc = rc;
  mFilterParametersSet = true;

  configureStateSpaceDimensionsIfReady();
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::setParameters(Real lf, Real cf, Real rf, Real rc,
                                          Real omegaN, Real vdRef, Real vqRef,
                                          Real kpVoltage, Real kiVoltage,
                                          Real kpCurrent, Real kiCurrent) {
  setParameters(omegaN, vdRef, vqRef);
  setControllerParameters(kpVoltage, kiVoltage, kpCurrent, kiCurrent, omegaN);
  setFilterParameters(lf, cf, rf, rc);
}

void EMT::Ph3::SSN_GFM_VCO::setInitialStateValues(Real voltageIntegratorD,
                                                  Real voltageIntegratorQ,
                                                  Real currentIntegratorD,
                                                  Real currentIntegratorQ) {
  if (!std::isfinite(voltageIntegratorD) ||
      !std::isfinite(voltageIntegratorQ) ||
      !std::isfinite(currentIntegratorD) || !std::isfinite(currentIntegratorQ))
    throw std::invalid_argument("Initial controller states must be finite.");

  mInitialVoltageIntegratorD = voltageIntegratorD;
  mInitialVoltageIntegratorQ = voltageIntegratorQ;
  mInitialCurrentIntegratorD = currentIntegratorD;
  mInitialCurrentIntegratorQ = currentIntegratorQ;
  mInitialStateOverrideEnabled = true;
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::clearInitialStateValues() {
  mInitialStateOverrideEnabled = false;
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::setInitialPower(Real activePower,
                                            Real reactivePower) {
  if (!std::isfinite(activePower) || !std::isfinite(reactivePower))
    throw std::invalid_argument(
        "Initial active/reactive power must be finite.");

  mInitialActivePower = activePower;
  mInitialReactivePower = reactivePower;
  mInitialPowerOverrideEnabled = true;
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::clearInitialPower() {
  mInitialPowerOverrideEnabled = false;
  markLinearizationDirty();
}

void EMT::Ph3::SSN_GFM_VCO::setNumericalLinearizationParameters(
    Real relativeStep, Real absoluteStep) {
  if (relativeStep <= 0.0)
    throw std::invalid_argument(
        "Relative finite-difference step must be positive.");
  if (absoluteStep <= 0.0)
    throw std::invalid_argument(
        "Absolute finite-difference step must be positive.");

  // Compatibility no-op: exact analytical Jacobians are used.
}

void EMT::Ph3::SSN_GFM_VCO::setLinearizationUpdateInterval(
    UInt updateInterval) {
  if (updateInterval == 0)
    throw std::invalid_argument(
        "Linearization update interval must be at least one step.");

  mLinearizationUpdateInterval = updateInterval;
  markLinearizationDirty();
}

Matrix EMT::Ph3::SSN_GFM_VCO::getParkTransformMatrix(Real theta) const {
  theta = std::remainder(theta, 2.0 * PI);

  Matrix transform(2, 3);
  const Real scale = std::sqrt(2.0 / 3.0);

  transform.row(0) << scale * std::cos(theta),
      scale * std::cos(theta - 2.0 * PI / 3.0),
      scale * std::cos(theta + 2.0 * PI / 3.0);

  transform.row(1) << -scale * std::sin(theta),
      -scale * std::sin(theta - 2.0 * PI / 3.0),
      -scale * std::sin(theta + 2.0 * PI / 3.0);

  return transform;
}

Matrix EMT::Ph3::SSN_GFM_VCO::getInverseParkTransformMatrix(Real theta) const {
  theta = std::remainder(theta, 2.0 * PI);

  Matrix transform(3, 2);
  const Real scale = std::sqrt(2.0 / 3.0);

  transform << scale * std::cos(theta), -scale * std::sin(theta),
      scale * std::cos(theta - 2.0 * PI / 3.0),
      -scale * std::sin(theta - 2.0 * PI / 3.0),
      scale * std::cos(theta + 2.0 * PI / 3.0),
      -scale * std::sin(theta + 2.0 * PI / 3.0);

  return transform;
}

void EMT::Ph3::SSN_GFM_VCO::evaluateStateDerivative(
    const Matrix &x, const Matrix &u, Matrix &stateDerivative) const {
  if (x.rows() != mStateSize || x.cols() != 1)
    throw std::invalid_argument(
        "SSN_GFM_VCO state vector has an invalid dimension.");
  if (u.rows() != mInputSize || u.cols() != 1)
    throw std::invalid_argument(
        "SSN_GFM_VCO input vector has an invalid dimension.");

  stateDerivative.setZero(mStateSize, 1);

  const Real theta = x(Theta, 0);
  const Real voltageIntegratorD = x(VoltageIntegratorD, 0);
  const Real voltageIntegratorQ = x(VoltageIntegratorQ, 0);
  const Real currentIntegratorD = x(CurrentIntegratorD, 0);
  const Real currentIntegratorQ = x(CurrentIntegratorQ, 0);

  const Matrix vcAbc = x.block(VcA, 0, 3, 1);
  const Matrix ifAbc = x.block(IfA, 0, 3, 1);

  const Matrix parkTransform = getParkTransformMatrix(theta);
  const Matrix inverseParkTransform = getInverseParkTransformMatrix(theta);

  const Matrix vcDq = parkTransform * vcAbc;
  const Matrix ifDq = parkTransform * ifAbc;

  const Real vcD = vcDq(0, 0);
  const Real vcQ = vcDq(1, 0);
  const Real ifD = ifDq(0, 0);
  const Real ifQ = ifDq(1, 0);

  // Fixed-frequency VCO.
  stateDerivative(Theta, 0) = mOmegaN;

  // Outer voltage PI controller.
  const Real voltageErrorD = mVdRef - vcD;
  const Real voltageErrorQ = mVqRef - vcQ;

  stateDerivative(VoltageIntegratorD, 0) = voltageErrorD;
  stateDerivative(VoltageIntegratorQ, 0) = voltageErrorQ;

  const Real currentReferenceD =
      mKpVoltage * voltageErrorD + mKiVoltage * voltageIntegratorD;
  const Real currentReferenceQ =
      mKpVoltage * voltageErrorQ + mKiVoltage * voltageIntegratorQ;

  // Inner current PI controller.
  const Real currentErrorD = currentReferenceD - ifD;
  const Real currentErrorQ = currentReferenceQ - ifQ;

  stateDerivative(CurrentIntegratorD, 0) = currentErrorD;
  stateDerivative(CurrentIntegratorQ, 0) = currentErrorQ;

  // This is algebraically identical to VoltageControllerVSI:
  // vInv = vc + Kp_i * (iRef - if) + Ki_i * gamma.
  const Real converterVoltageD =
      vcD + mKpCurrent * currentErrorD + mKiCurrent * currentIntegratorD;
  const Real converterVoltageQ =
      vcQ + mKpCurrent * currentErrorQ + mKiCurrent * currentIntegratorQ;

  Matrix converterVoltageDq(2, 1);
  converterVoltageDq << converterVoltageD, converterVoltageQ;

  const Matrix converterVoltageAbc = inverseParkTransform * converterVoltageDq;

  // Electrical L-C-R filter.
  const Matrix vcDerivative = ifAbc / mCf + (u - vcAbc) / (mCf * mRc);
  const Matrix ifDerivative = (converterVoltageAbc - vcAbc - mRf * ifAbc) / mLf;

  stateDerivative.block(VcA, 0, 3, 1) = vcDerivative;
  stateDerivative.block(IfA, 0, 3, 1) = ifDerivative;
}

void EMT::Ph3::SSN_GFM_VCO::evaluateOutput(const Matrix &x, const Matrix &u,
                                           Matrix &output) const {
  if (x.rows() != mStateSize || x.cols() != 1)
    throw std::invalid_argument(
        "SSN_GFM_VCO state vector has an invalid dimension.");
  if (u.rows() != mInputSize || u.cols() != 1)
    throw std::invalid_argument(
        "SSN_GFM_VCO input vector has an invalid dimension.");

  const Matrix vcAbc = x.block(VcA, 0, 3, 1);
  output = (u - vcAbc) / mRc;
}

void EMT::Ph3::SSN_GFM_VCO::calculateAnalyticalJacobians(const Matrix &x,
                                                         const Matrix &u,
                                                         Matrix &A, Matrix &B,
                                                         Matrix &C,
                                                         Matrix &D) const {
  if (x.rows() != mStateSize || x.cols() != 1)
    throw std::invalid_argument(
        "SSN_GFM_VCO state vector has an invalid dimension.");
  if (u.rows() != mInputSize || u.cols() != 1)
    throw std::invalid_argument(
        "SSN_GFM_VCO input vector has an invalid dimension.");

  const Real theta = x(Theta, 0);
  const Real voltageIntegratorD = x(VoltageIntegratorD, 0);
  const Real voltageIntegratorQ = x(VoltageIntegratorQ, 0);
  const Real currentIntegratorD = x(CurrentIntegratorD, 0);
  const Real currentIntegratorQ = x(CurrentIntegratorQ, 0);

  const Matrix vcAbc = x.block(VcA, 0, 3, 1);
  const Matrix ifAbc = x.block(IfA, 0, 3, 1);

  const Matrix identity3 = Matrix::Identity(3, 3);

  const Matrix parkTransform = getParkTransformMatrix(theta);
  const Matrix tD = parkTransform.row(0);
  const Matrix tQ = parkTransform.row(1);

  const Matrix inverseParkTransform = getInverseParkTransformMatrix(theta);
  const Matrix sD = inverseParkTransform.col(0);
  const Matrix sQ = inverseParkTransform.col(1);

  // Exact transform derivatives:
  // d(tD)/dtheta = tQ, d(tQ)/dtheta = -tD,
  // d(sD)/dtheta = sQ, d(sQ)/dtheta = -sD.
  const Real vcD = (tD * vcAbc)(0, 0);
  const Real vcQ = (tQ * vcAbc)(0, 0);
  const Real ifD = (tD * ifAbc)(0, 0);
  const Real ifQ = (tQ * ifAbc)(0, 0);

  Matrix dVcDByX = Matrix::Zero(1, mStateSize);
  Matrix dVcQByX = Matrix::Zero(1, mStateSize);
  Matrix dIfDByX = Matrix::Zero(1, mStateSize);
  Matrix dIfQByX = Matrix::Zero(1, mStateSize);

  dVcDByX(0, Theta) = vcQ;
  dVcDByX.block(0, VcA, 1, 3) = tD;
  dVcQByX(0, Theta) = -vcD;
  dVcQByX.block(0, VcA, 1, 3) = tQ;

  dIfDByX(0, Theta) = ifQ;
  dIfDByX.block(0, IfA, 1, 3) = tD;
  dIfQByX(0, Theta) = -ifD;
  dIfQByX.block(0, IfA, 1, 3) = tQ;

  Matrix unitVoltageIntegratorD = Matrix::Zero(1, mStateSize);
  Matrix unitVoltageIntegratorQ = Matrix::Zero(1, mStateSize);
  Matrix unitCurrentIntegratorD = Matrix::Zero(1, mStateSize);
  Matrix unitCurrentIntegratorQ = Matrix::Zero(1, mStateSize);

  unitVoltageIntegratorD(0, VoltageIntegratorD) = 1.0;
  unitVoltageIntegratorQ(0, VoltageIntegratorQ) = 1.0;
  unitCurrentIntegratorD(0, CurrentIntegratorD) = 1.0;
  unitCurrentIntegratorQ(0, CurrentIntegratorQ) = 1.0;

  const Matrix dVoltageErrorDByX = -dVcDByX;
  const Matrix dVoltageErrorQByX = -dVcQByX;

  const Matrix dCurrentReferenceDByX =
      mKpVoltage * dVoltageErrorDByX + mKiVoltage * unitVoltageIntegratorD;
  const Matrix dCurrentReferenceQByX =
      mKpVoltage * dVoltageErrorQByX + mKiVoltage * unitVoltageIntegratorQ;

  const Matrix dCurrentErrorDByX = dCurrentReferenceDByX - dIfDByX;
  const Matrix dCurrentErrorQByX = dCurrentReferenceQByX - dIfQByX;

  const Real voltageErrorD = mVdRef - vcD;
  const Real voltageErrorQ = mVqRef - vcQ;
  const Real currentReferenceD =
      mKpVoltage * voltageErrorD + mKiVoltage * voltageIntegratorD;
  const Real currentReferenceQ =
      mKpVoltage * voltageErrorQ + mKiVoltage * voltageIntegratorQ;
  const Real currentErrorD = currentReferenceD - ifD;
  const Real currentErrorQ = currentReferenceQ - ifQ;

  const Real converterVoltageD =
      vcD + mKpCurrent * currentErrorD + mKiCurrent * currentIntegratorD;
  const Real converterVoltageQ =
      vcQ + mKpCurrent * currentErrorQ + mKiCurrent * currentIntegratorQ;

  const Matrix dConverterVoltageDByX = dVcDByX +
                                       mKpCurrent * dCurrentErrorDByX +
                                       mKiCurrent * unitCurrentIntegratorD;
  const Matrix dConverterVoltageQByX = dVcQByX +
                                       mKpCurrent * dCurrentErrorQByX +
                                       mKiCurrent * unitCurrentIntegratorQ;

  Matrix dConverterVoltageAbcByX =
      sD * dConverterVoltageDByX + sQ * dConverterVoltageQByX;
  dConverterVoltageAbcByX.col(Theta) +=
      sQ * converterVoltageD - sD * converterVoltageQ;

  A.setZero(mStateSize, mStateSize);
  B.setZero(mStateSize, mInputSize);
  C.setZero(mOutputSize, mStateSize);
  D.setZero(mOutputSize, mInputSize);

  // Fixed-frequency VCO: theta_dot is a constant, so its Jacobian row is zero.

  // Cascaded PI controller.
  A.row(VoltageIntegratorD) = dVoltageErrorDByX;
  A.row(VoltageIntegratorQ) = dVoltageErrorQByX;
  A.row(CurrentIntegratorD) = dCurrentErrorDByX;
  A.row(CurrentIntegratorQ) = dCurrentErrorQByX;

  // Electrical filter.
  A.block(VcA, VcA, 3, 3) = -1.0 / (mCf * mRc) * identity3;
  A.block(VcA, IfA, 3, 3) = 1.0 / mCf * identity3;
  B.block(VcA, 0, 3, 3) = 1.0 / (mCf * mRc) * identity3;

  A.block(IfA, 0, 3, mStateSize) = dConverterVoltageAbcByX / mLf;
  A.block(IfA, VcA, 3, 3) -= 1.0 / mLf * identity3;
  A.block(IfA, IfA, 3, 3) -= mRf / mLf * identity3;

  // SSN output y = (u - vc) / Rc.
  C.block(0, VcA, 3, 3) = -1.0 / mRc * identity3;
  D = 1.0 / mRc * identity3;
}

void EMT::Ph3::SSN_GFM_VCO::buildStateSpaceModel(const Matrix &x,
                                                 const Matrix &u, Matrix &A,
                                                 Matrix &B, Matrix &C,
                                                 Matrix &D, Matrix &E,
                                                 Matrix &F) const {
  calculateAnalyticalJacobians(x, u, A, B, C, D);

  Matrix stateDerivative = Matrix::Zero(mStateSize, 1);
  Matrix output = Matrix::Zero(mOutputSize, 1);

  evaluateStateDerivative(x, u, stateDerivative);
  evaluateOutput(x, u, output);

  E = stateDerivative - A * x - B * u;
  F = output - C * x - D * u;
}

Bool EMT::Ph3::SSN_GFM_VCO::updateComponentParameters() {
  if (mLinearizationInitialized && mLinearizationUpdateInterval > 1) {
    ++mStepsSinceLinearization;
    if (mStepsSinceLinearization < mLinearizationUpdateInterval)
      return false;
  }

  Matrix eVector;
  Matrix fVector;

  buildStateSpaceModel(**mX, **mIntfVoltage, mA, mB, mC, mD, eVector, fVector);

  setStateOffset(eVector);
  setOutputOffset(fVector);

  mLinearizationInitialized = true;
  mStepsSinceLinearization = 0;
  return true;
}

void EMT::Ph3::SSN_GFM_VCO::updateLogAttributes(const Matrix &u) const {
  const Matrix &x = **mX;
  const Real theta = x(Theta, 0);

  const Matrix parkTransform = getParkTransformMatrix(theta);
  const Matrix vcAbc = x.block(VcA, 0, 3, 1);
  const Matrix ifAbc = x.block(IfA, 0, 3, 1);
  const Matrix iGridAbc = (vcAbc - u) / mRc;

  const Matrix vcDq = parkTransform * vcAbc;
  const Matrix ifDq = parkTransform * ifAbc;
  const Matrix iGridDq = parkTransform * iGridAbc;
  const Matrix uDq = parkTransform * u;

  const Real vcD = vcDq(0, 0);
  const Real vcQ = vcDq(1, 0);
  const Real ifD = ifDq(0, 0);
  const Real ifQ = ifDq(1, 0);
  const Real iGridD = iGridDq(0, 0);
  const Real iGridQ = iGridDq(1, 0);

  const Real voltageErrorD = mVdRef - vcD;
  const Real voltageErrorQ = mVqRef - vcQ;
  const Real currentReferenceD =
      mKpVoltage * voltageErrorD + mKiVoltage * x(VoltageIntegratorD, 0);
  const Real currentReferenceQ =
      mKpVoltage * voltageErrorQ + mKiVoltage * x(VoltageIntegratorQ, 0);
  const Real converterVoltageD = vcD + mKpCurrent * (currentReferenceD - ifD) +
                                 mKiCurrent * x(CurrentIntegratorD, 0);
  const Real converterVoltageQ = vcQ + mKpCurrent * (currentReferenceQ - ifQ) +
                                 mKiCurrent * x(CurrentIntegratorQ, 0);

  **mVcD = vcD;
  **mVcQ = vcQ;
  **mIGridD = iGridD;
  **mIGridQ = iGridQ;
  **mIfD = ifD;
  **mIfQ = ifQ;

  // Power-invariant dq transform: no additional 3/2 factor is required.
  **mPInst = uDq(0, 0) * iGridD + uDq(1, 0) * iGridQ;
  **mQInst = uDq(1, 0) * iGridD - uDq(0, 0) * iGridQ;

  **mOmegaGFM = mOmegaN;
  **mThetaGFM = theta;
  **mVoltageMagnitudeGFM = std::hypot(mVdRef, mVqRef);
  **mVoltageReferenceD = mVdRef;
  **mVoltageReferenceQ = mVqRef;
  **mConverterVoltageD = converterVoltageD;
  **mConverterVoltageQ = converterVoltageQ;
}

void EMT::Ph3::SSN_GFM_VCO::initializeFromNodesAndTerminals(Real frequency) {
  if (!mParametersSet)
    throw std::logic_error(
        "setParameters(), setControllerParameters() and "
        "setFilterParameters() must be called before initialization.");

  const Real omegaInitialization = mOmegaN;
  const Complex imaginaryUnit(0.0, 1.0);

  const MatrixComp uPhasor = buildInitialInputFromNodes(frequency);

  Complex initialPower;
  if (mInitialPowerOverrideEnabled) {
    initialPower = Complex(mInitialActivePower, mInitialReactivePower);
  } else {
    // This matches VSIVoltageControlVCO. For EMT initialization from an SP
    // synchronous generator, SystemTopology::initWithPowerflow() transfers the
    // generator power to terminal 0.
    initialPower = terminal(0)->singlePower();
  }

  MatrixComp yPhasor = MatrixComp::Zero(3, 1);
  if (std::abs(uPhasor(0, 0)) > mInitializationTolerance) {
    const Complex currentA =
        -std::conj((2.0 / 3.0) * initialPower / uPhasor(0, 0));
    yPhasor << currentA, currentA * SHIFT_TO_PHASE_B,
        currentA * SHIFT_TO_PHASE_C;
  }

  // y is current entering the component. Physical grid injection is -y.
  const MatrixComp vcPhasor = uPhasor - mRc * yPhasor;
  const MatrixComp capacitorCurrentPhasor =
      imaginaryUnit * omegaInitialization * mCf * vcPhasor;
  const MatrixComp ifPhasor = capacitorCurrentPhasor - yPhasor;
  const MatrixComp converterVoltagePhasor =
      vcPhasor + (mRf + imaginaryUnit * omegaInitialization * mLf) * ifPhasor;

  const Matrix vcAbc0 = vcPhasor.real();
  const Matrix ifAbc0 = ifPhasor.real();
  const Matrix converterVoltageAbc0 = converterVoltagePhasor.real();

  const Real theta0 = std::arg(vcPhasor(0, 0));
  const Matrix parkTransform = getParkTransformMatrix(theta0);

  const Matrix vcDq0 = parkTransform * vcAbc0;
  const Matrix ifDq0 = parkTransform * ifAbc0;
  const Matrix converterVoltageDq0 = parkTransform * converterVoltageAbc0;

  const Real vcD0 = vcDq0(0, 0);
  const Real vcQ0 = vcDq0(1, 0);
  const Real ifD0 = ifDq0(0, 0);
  const Real ifQ0 = ifDq0(1, 0);

  Matrix x0 = Matrix::Zero(mStateSize, 1);
  x0(Theta, 0) = theta0;

  if (mInitialStateOverrideEnabled) {
    x0(VoltageIntegratorD, 0) = mInitialVoltageIntegratorD;
    x0(VoltageIntegratorQ, 0) = mInitialVoltageIntegratorQ;
    x0(CurrentIntegratorD, 0) = mInitialCurrentIntegratorD;
    x0(CurrentIntegratorQ, 0) = mInitialCurrentIntegratorQ;
  } else {
    const Real voltageErrorD0 = mVdRef - vcD0;
    const Real voltageErrorQ0 = mVqRef - vcQ0;

    // Set iRef = if initially.
    x0(VoltageIntegratorD, 0) =
        (ifD0 - mKpVoltage * voltageErrorD0) / mKiVoltage;
    x0(VoltageIntegratorQ, 0) =
        (ifQ0 - mKpVoltage * voltageErrorQ0) / mKiVoltage;

    // With zero initial current error, solve the current PI output for the
    // converter bridge voltage required by the phasor operating point.
    x0(CurrentIntegratorD, 0) = (converterVoltageDq0(0, 0) - vcD0) / mKiCurrent;
    x0(CurrentIntegratorQ, 0) = (converterVoltageDq0(1, 0) - vcQ0) / mKiCurrent;
  }

  x0.block(VcA, 0, 3, 1) = vcAbc0;
  x0.block(IfA, 0, 3, 1) = ifAbc0;

  **mX = x0;
  **mIntfVoltage = uPhasor.real();

  markLinearizationDirty();
  updateStateSpaceModel();

  mYHist = calculateHistoryVector();
  **mIntfCurrent = mW * (**mIntfVoltage) + mYHist;

  updateLogAttributes(**mIntfVoltage);

  Matrix stateDerivative = Matrix::Zero(mStateSize, 1);
  Matrix nonlinearOutput = Matrix::Zero(mOutputSize, 1);

  evaluateStateDerivative(**mX, **mIntfVoltage, stateDerivative);
  evaluateOutput(**mX, **mIntfVoltage, nonlinearOutput);

  const Matrix ssnOutput = mW * (**mIntfVoltage) + mYHist;

  SPDLOG_LOGGER_INFO(
      mSLog,
      "\n--- SSN GFM VCO initialization ---"
      "\nInput voltage u: {:s}"
      "\nInterface current y: {:s}"
      "\nInitial transferred/overridden power: [{:.6e}, {:.6e}]"
      "\nState x: {:s}"
      "\nState derivative norm: {:.6e}"
      "\nNonlinear output: {:s}"
      "\nSSN output: {:s}"
      "\nOutput mismatch norm: {:.6e}"
      "\nW norm: {:.6e}"
      "\nHistory-vector norm: {:.6e}"
      "\nVc dq: [{:.6e}, {:.6e}]"
      "\nIf dq: [{:.6e}, {:.6e}]"
      "\nConverter voltage dq: [{:.6e}, {:.6e}]"
      "\n--- SSN GFM VCO initialization finished ---",
      Logger::matrixToString(**mIntfVoltage),
      Logger::matrixToString(**mIntfCurrent), initialPower.real(),
      initialPower.imag(), Logger::matrixToString(**mX), stateDerivative.norm(),
      Logger::matrixToString(nonlinearOutput),
      Logger::matrixToString(ssnOutput), (nonlinearOutput - ssnOutput).norm(),
      mW.norm(), mYHist.norm(), vcD0, vcQ0, ifD0, ifQ0,
      converterVoltageDq0(0, 0), converterVoltageDq0(1, 0));
}

Matrix EMT::Ph3::SSN_GFM_VCO::getState() const { return **mX; }

Matrix EMT::Ph3::SSN_GFM_VCO::getStateDerivative() const {
  Matrix stateDerivative = Matrix::Zero(mStateSize, 1);
  evaluateStateDerivative(**mX, **mIntfVoltage, stateDerivative);
  return stateDerivative;
}

Matrix EMT::Ph3::SSN_GFM_VCO::getInterfaceVoltage() const {
  return **mIntfVoltage;
}

Matrix EMT::Ph3::SSN_GFM_VCO::getInterfaceCurrent() const {
  return **mIntfCurrent;
}
