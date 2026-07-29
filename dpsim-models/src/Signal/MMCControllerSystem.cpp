// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/Signal/MMCControllerSystem.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/QR>

using namespace CPS;
using namespace CPS::Signal;

namespace {
MMCStateSlice assignSlice(UInt &offset, UInt count) {
  MMCStateSlice slice{offset, count};
  offset += count;
  return slice;
}

void requireFinite(Real value, const char *name) {
  if (!std::isfinite(value))
    throw std::invalid_argument(String(name) + " must be finite.");
}
} // namespace

MMCControllerSystem::Slices MMCControllerSystem::makeSlices() {
  Slices s;
  UInt offset = 0;
  s.active = assignSlice(offset, 1);
  s.reactive = assignSlice(offset, 1);
  s.occ = assignSlice(offset, 2);
  s.ccc = assignSlice(offset, 2);
  s.zcc = assignSlice(offset, 1);
  s.energy = assignSlice(offset, 1);
  s.pll = assignSlice(offset, 2);
  s.vacFilter = assignSlice(offset, 2);
  s.pFilter = assignSlice(offset, 2);
  s.qFilter = assignSlice(offset, 2);
  s.vdcFilter = assignSlice(offset, 2);
  s.vacMagnitudeFilter = assignSlice(offset, 2);
  s.delay = assignSlice(offset, 10);
  s.total = offset;
  return s;
}

MMCControllerSystem::MMCControllerSystem()
    : mSlices(makeSlices()), mVacFilter(2), mPFilter(1), mQFilter(1),
      mVdcFilter(1), mVacMagnitudeFilter(1), mDelay(5), mOmegaN(0.0),
      mEquivalentAcResistance(0.0), mEquivalentAcInductance(0.0),
      mArmResistance(0.0), mArmInductance(0.0), mNominalDcVoltage(0.0),
      mEnergyReference(0.0), mEnergyKp(0.0), mEnergyKi(0.0),
      mActivePowerReference(0.0), mDcVoltageReference(0.0), mDroopGain(1.0),
      mOpenLoopIDeltaDReference(0.0), mReactivePowerReference(0.0),
      mAcVoltageReference(0.0), mOpenLoopIDeltaQReference(0.0),
      mISigmaDReference(0.0), mISigmaQReference(0.0), mISigmaZReference(0.0),
      mFeedforwardSampleTime(40e-6), mActiveMode(MMCActiveMode::OpenLoop),
      mReactiveMode(MMCReactiveMode::OpenLoop), mEnergyEnabled(false),
      mDelayEnabled(false), mFeedforwardEnabled(false), mPllEnabled(false),
      mOccEnabled(false), mCccEnabled(false), mZccEnabled(false),
      mVacFilterEnabled(false), mPFilterEnabled(false), mQFilterEnabled(false),
      mVdcFilterEnabled(false), mVacMagnitudeFilterEnabled(false),
      mControlSource(ControlSource::InternalControllers),
      mExternalDifferentialVoltage(Matrix::Zero(2, 1)),
      mExternalCommonModeVoltage(Matrix::Zero(3, 1)),
      mMaximumAcCurrent(std::numeric_limits<Real>::infinity()),
      mMaximumCirculatingCurrent(std::numeric_limits<Real>::infinity()),
      mMaximumModulationMagnitude(std::numeric_limits<Real>::infinity()) {}

void MMCControllerSystem::setElectricalParameters(
    Real nominalOmega, Real equivalentAcResistance, Real equivalentAcInductance,
    Real armResistance, Real armInductance, Real nominalDcVoltage,
    Real energyReference) {
  if (!std::isfinite(nominalOmega) || nominalOmega <= 0.0 ||
      !std::isfinite(equivalentAcResistance) || equivalentAcResistance < 0.0 ||
      !std::isfinite(equivalentAcInductance) || equivalentAcInductance <= 0.0 ||
      !std::isfinite(armResistance) || armResistance < 0.0 ||
      !std::isfinite(armInductance) || armInductance <= 0.0 ||
      !std::isfinite(nominalDcVoltage) || nominalDcVoltage <= 0.0 ||
      !std::isfinite(energyReference) || energyReference < 0.0)
    throw std::invalid_argument(
        "Controller-system electrical parameters are invalid.");
  mOmegaN = nominalOmega;
  mEquivalentAcResistance = equivalentAcResistance;
  mEquivalentAcInductance = equivalentAcInductance;
  mArmResistance = armResistance;
  mArmInductance = armInductance;
  mNominalDcVoltage = nominalDcVoltage;
  mEnergyReference = energyReference;
  mEnergy.setParameters(mEnergyKp, mEnergyKi, mEnergyEnabled, mEnergyReference,
                        1.0);
}

void MMCControllerSystem::setPLL(Real kp, Real ki, Bool enabled) {
  mPllEnabled = enabled;
  mPll.setParameters(kp, ki, enabled);
}
void MMCControllerSystem::setMMCOutputCurrentController(Real kp, Real ki) {
  mOcc.setPI(kp, ki);
  mOccEnabled = kp > 0.0 || ki > 0.0;
}
void MMCControllerSystem::setMMCCirculatingCurrentController(Real kp, Real ki) {
  mCcc.setPI(kp, ki);
  mCccEnabled = kp > 0.0 || ki > 0.0;
}
void MMCControllerSystem::setMMCZeroSequenceCurrentController(Real kp,
                                                              Real ki) {
  mZcc.setPI(kp, ki);
  mZccEnabled = kp > 0.0 || ki > 0.0;
}
void MMCControllerSystem::setMMCEnergyController(Real kp, Real ki,
                                                 Bool enabled) {
  mEnergyKp = kp;
  mEnergyKi = ki;
  mEnergyEnabled = enabled;
  mEnergy.setParameters(mEnergyKp, mEnergyKi, enabled, mEnergyReference, 1.0);
}

void MMCControllerSystem::setActivePowerControl(Real reference, Real kp,
                                                Real ki) {
  requireFinite(reference, "Active-power reference");
  mActiveMode = MMCActiveMode::ActivePower;
  mFeedforwardEnabled = false;
  mActivePowerReference = reference;
  mActive.setPI(kp, ki);
  mActive.setMode(mActiveMode);
  mActive.setReferences(mActivePowerReference, mDcVoltageReference, mDroopGain,
                        mOpenLoopIDeltaDReference);
}
void MMCControllerSystem::setActivePowerFeedforwardControl(
    Real reference, Real cutoffFrequency, Real sampleTime,
    Real minimumDaxisVoltage) {
  requireFinite(reference, "Active-power feedforward reference");
  if (!std::isfinite(cutoffFrequency) || cutoffFrequency <= 0.0 ||
      !std::isfinite(sampleTime) || sampleTime <= 0.0 ||
      !std::isfinite(minimumDaxisVoltage) || minimumDaxisVoltage < 0.0)
    throw std::invalid_argument(
        "Sampled active-power feedforward parameters are invalid.");
  mActiveMode = MMCActiveMode::SampledPowerFeedforward;
  mFeedforwardEnabled = true;
  mFeedforwardSampleTime = sampleTime;
  mActivePowerReference = reference;
  mActive.setMode(mActiveMode);
  mActive.setReferences(mActivePowerReference, mDcVoltageReference, mDroopGain,
                        mOpenLoopIDeltaDReference);
}
void MMCControllerSystem::setActivePowerReference(Real reference) {
  requireFinite(reference, "Active-power reference");
  mActivePowerReference = reference;
  mActive.setReferences(mActivePowerReference, mDcVoltageReference, mDroopGain,
                        mOpenLoopIDeltaDReference);
}
void MMCControllerSystem::setDcVoltageControl(Real reference, Real kp,
                                              Real ki) {
  if (!std::isfinite(reference) || reference <= 0.0)
    throw std::invalid_argument("DC-voltage reference must be positive.");
  mActiveMode = MMCActiveMode::DcVoltage;
  mFeedforwardEnabled = false;
  mDcVoltageReference = reference;
  mActive.setPI(kp, ki);
  mActive.setMode(mActiveMode);
  mActive.setReferences(mActivePowerReference, mDcVoltageReference, mDroopGain,
                        mOpenLoopIDeltaDReference);
}
void MMCControllerSystem::setDcDroopControl(Real activePowerReference,
                                            Real dcVoltageReference,
                                            Real droopGain) {
  if (!std::isfinite(dcVoltageReference) || dcVoltageReference <= 0.0 ||
      !std::isfinite(droopGain) || droopGain == 0.0)
    throw std::invalid_argument("DC-droop parameters are invalid.");
  mActiveMode = MMCActiveMode::DcDroop;
  mFeedforwardEnabled = false;
  mActivePowerReference = activePowerReference;
  mDcVoltageReference = dcVoltageReference;
  mDroopGain = droopGain;
  mActive.setMode(mActiveMode);
  mActive.setReferences(mActivePowerReference, mDcVoltageReference, mDroopGain,
                        mOpenLoopIDeltaDReference);
}
void MMCControllerSystem::setActiveOpenLoop(Real currentReference) {
  requireFinite(currentReference, "Open-loop d-axis current reference");
  mActiveMode = MMCActiveMode::OpenLoop;
  mFeedforwardEnabled = false;
  mOpenLoopIDeltaDReference = currentReference;
  mActive.setMode(mActiveMode);
  mActive.setReferences(mActivePowerReference, mDcVoltageReference, mDroopGain,
                        mOpenLoopIDeltaDReference);
}

void MMCControllerSystem::setReactivePowerControl(Real reference, Real kp,
                                                  Real ki) {
  requireFinite(reference, "Reactive-power reference");
  mReactiveMode = MMCReactiveMode::ReactivePower;
  mReactivePowerReference = reference;
  mReactive.setPI(kp, ki);
  mReactive.setMode(mReactiveMode);
  mReactive.setReferences(mReactivePowerReference, mAcVoltageReference,
                          mOpenLoopIDeltaQReference);
}
void MMCControllerSystem::setAcVoltageControl(Real reference, Real kp,
                                              Real ki) {
  if (!std::isfinite(reference) || reference <= 0.0)
    throw std::invalid_argument("AC-voltage reference must be positive.");
  mReactiveMode = MMCReactiveMode::AcVoltage;
  mAcVoltageReference = reference;
  mReactive.setPI(kp, ki);
  mReactive.setMode(mReactiveMode);
  mReactive.setReferences(mReactivePowerReference, mAcVoltageReference,
                          mOpenLoopIDeltaQReference);
}
void MMCControllerSystem::setReactiveOpenLoop(Real currentReference) {
  requireFinite(currentReference, "Open-loop q-axis current reference");
  mReactiveMode = MMCReactiveMode::OpenLoop;
  mOpenLoopIDeltaQReference = currentReference;
  mReactive.setMode(mReactiveMode);
  mReactive.setReferences(mReactivePowerReference, mAcVoltageReference,
                          mOpenLoopIDeltaQReference);
}

void MMCControllerSystem::setCirculatingCurrentReferences(Real dReference,
                                                          Real qReference,
                                                          Real zReference) {
  requireFinite(dReference, "Sigma-d current reference");
  requireFinite(qReference, "Sigma-q current reference");
  requireFinite(zReference, "Sigma-z current reference");
  mISigmaDReference = dReference;
  mISigmaQReference = qReference;
  mISigmaZReference = zReference;
}
void MMCControllerSystem::setMeasurementFilters(
    Real acVoltageDqTimeConstant, Real activePowerTimeConstant,
    Real reactivePowerTimeConstant, Real dcVoltageTimeConstant,
    Real acVoltageMagnitudeTimeConstant) {
  if (!std::isfinite(acVoltageDqTimeConstant) ||
      acVoltageDqTimeConstant < 0.0 ||
      !std::isfinite(activePowerTimeConstant) ||
      activePowerTimeConstant < 0.0 ||
      !std::isfinite(reactivePowerTimeConstant) ||
      reactivePowerTimeConstant < 0.0 ||
      !std::isfinite(dcVoltageTimeConstant) || dcVoltageTimeConstant < 0.0 ||
      !std::isfinite(acVoltageMagnitudeTimeConstant) ||
      acVoltageMagnitudeTimeConstant < 0.0)
    throw std::invalid_argument(
        "MMC measurement-filter time constants must be finite and "
        "non-negative.");

  const Real vacDq = acVoltageDqTimeConstant;
  const Real p = activePowerTimeConstant;
  const Real q = reactivePowerTimeConstant;
  const Real vdc = dcVoltageTimeConstant;
  const Real vacMagnitude = acVoltageMagnitudeTimeConstant;
  mVacFilterEnabled = vacDq > 0.0;
  mPFilterEnabled = p > 0.0;
  mQFilterEnabled = q > 0.0;
  mVdcFilterEnabled = vdc > 0.0;
  mVacMagnitudeFilterEnabled = vacMagnitude > 0.0;
  mVacFilter.setTimeConstant(vacDq);
  mPFilter.setTimeConstant(p);
  mQFilter.setTimeConstant(q);
  mVdcFilter.setTimeConstant(vdc);
  mVacMagnitudeFilter.setTimeConstant(vacMagnitude);
}
void MMCControllerSystem::setModulationDelay(Real timeDelay, Bool enabled) {
  mDelayEnabled = enabled && timeDelay > 0.0;
  mDelay.setDelay(timeDelay, mDelayEnabled);
}
void MMCControllerSystem::setLimits(Real maximumAcCurrent,
                                    Real maximumCirculatingCurrent,
                                    Real maximumModulationMagnitude) {
  if (!std::isfinite(maximumAcCurrent) || maximumAcCurrent <= 0.0 ||
      !std::isfinite(maximumCirculatingCurrent) ||
      maximumCirculatingCurrent <= 0.0 ||
      !std::isfinite(maximumModulationMagnitude) ||
      maximumModulationMagnitude <= 0.0)
    throw std::invalid_argument("MMC controller limits must be positive.");
  mMaximumAcCurrent = maximumAcCurrent;
  mMaximumCirculatingCurrent = maximumCirculatingCurrent;
  mMaximumModulationMagnitude = maximumModulationMagnitude;
}

void MMCControllerSystem::setControlSource(ControlSource source) {
  mControlSource = source;
}
void MMCControllerSystem::setExternalDifferentialVoltage(Real dVolts,
                                                         Real qVolts) {
  requireFinite(dVolts, "External differential d voltage");
  requireFinite(qVolts, "External differential q voltage");
  mExternalDifferentialVoltage << dVolts, qVolts;
}
void MMCControllerSystem::setExternalCommonModeVoltage(Real dVolts, Real qVolts,
                                                       Real zVolts) {
  requireFinite(dVolts, "External common-mode d voltage");
  requireFinite(qVolts, "External common-mode q voltage");
  requireFinite(zVolts, "External common-mode z voltage");
  mExternalCommonModeVoltage << dVolts, qVolts, zVolts;
}
MMCControllerSystem::ControlSource MMCControllerSystem::controlSource() const {
  return mControlSource;
}

MMCActiveMode MMCControllerSystem::activeMode() const { return mActiveMode; }
MMCReactiveMode MMCControllerSystem::reactiveMode() const {
  return mReactiveMode;
}
Real MMCControllerSystem::activePowerReference() const {
  return mActivePowerReference;
}
Real MMCControllerSystem::dcVoltageReference() const {
  return mDcVoltageReference;
}
Real MMCControllerSystem::reactivePowerReference() const {
  return mReactivePowerReference;
}
Real MMCControllerSystem::acVoltageReference() const {
  return mAcVoltageReference;
}
Real MMCControllerSystem::sigmaDReference() const { return mISigmaDReference; }
Real MMCControllerSystem::sigmaQReference() const { return mISigmaQReference; }
Real MMCControllerSystem::sigmaZReference() const { return mISigmaZReference; }
Bool MMCControllerSystem::energyControllerEnabled() const {
  return mEnergyEnabled;
}
Bool MMCControllerSystem::sampledFeedforwardEnabled() const {
  return mFeedforwardEnabled;
}
Real MMCControllerSystem::sampledFeedforwardPeriod() const {
  return mFeedforwardSampleTime;
}
UInt MMCControllerSystem::activeIntegratorStateIndex() const {
  return mSlices.active.first;
}

std::vector<UInt> MMCControllerSystem::activeStateIndices() const {
  std::vector<UInt> indices;
  indices.reserve(stateSize());
  const auto appendSlice = [&](const MMCStateSlice &slice) {
    for (UInt index = 0; index < slice.count; ++index)
      indices.push_back(slice.first + index);
  };

  if (mActiveMode == MMCActiveMode::ActivePower ||
      mActiveMode == MMCActiveMode::DcVoltage)
    appendSlice(mSlices.active);
  if (mReactiveMode == MMCReactiveMode::ReactivePower ||
      mReactiveMode == MMCReactiveMode::AcVoltage)
    appendSlice(mSlices.reactive);
  if (mControlSource == ControlSource::InternalControllers && mOccEnabled)
    appendSlice(mSlices.occ);
  if (mControlSource != ControlSource::ExternalFullConverterVoltage) {
    if (mCccEnabled)
      appendSlice(mSlices.ccc);
    if (mZccEnabled)
      appendSlice(mSlices.zcc);
  }
  if (mEnergyEnabled)
    appendSlice(mSlices.energy);
  if (mPllEnabled)
    appendSlice(mSlices.pll);
  if (mVacFilterEnabled)
    appendSlice(mSlices.vacFilter);
  if (mPFilterEnabled)
    appendSlice(mSlices.pFilter);
  if (mQFilterEnabled)
    appendSlice(mSlices.qFilter);
  if (mVdcFilterEnabled)
    appendSlice(mSlices.vdcFilter);
  if (mVacMagnitudeFilterEnabled)
    appendSlice(mSlices.vacMagnitudeFilter);
  if (mDelayEnabled)
    appendSlice(mSlices.delay);
  return indices;
}

std::vector<UInt> MMCControllerSystem::equilibriumStateIndices() const {
  std::vector<UInt> indices;
  indices.reserve(stateSize());
  const auto appendSlice = [&](const MMCStateSlice &slice) {
    for (UInt index = 0; index < slice.count; ++index)
      indices.push_back(slice.first + index);
  };

  if (mActiveMode == MMCActiveMode::ActivePower ||
      mActiveMode == MMCActiveMode::DcVoltage)
    appendSlice(mSlices.active);
  if (mReactiveMode == MMCReactiveMode::ReactivePower ||
      mReactiveMode == MMCReactiveMode::AcVoltage)
    appendSlice(mSlices.reactive);
  if (mControlSource == ControlSource::InternalControllers && mOccEnabled)
    appendSlice(mSlices.occ);
  if (mControlSource != ControlSource::ExternalFullConverterVoltage) {
    if (mCccEnabled)
      appendSlice(mSlices.ccc);
    if (mZccEnabled)
      appendSlice(mSlices.zcc);
  }
  if (mEnergyEnabled)
    appendSlice(mSlices.energy);
  if (mVacFilterEnabled)
    appendSlice(mSlices.vacFilter);
  if (mPFilterEnabled)
    appendSlice(mSlices.pFilter);
  if (mQFilterEnabled)
    appendSlice(mSlices.qFilter);
  if (mVdcFilterEnabled)
    appendSlice(mSlices.vdcFilter);
  if (mVacMagnitudeFilterEnabled)
    appendSlice(mSlices.vacMagnitudeFilter);
  if (mDelayEnabled)
    appendSlice(mSlices.delay);
  return indices;
}

UInt MMCControllerSystem::stateSize() const { return mSlices.total; }
UInt MMCControllerSystem::inputSize() const { return InputCount; }
UInt MMCControllerSystem::outputSize() const { return OutputCount; }

std::vector<String> MMCControllerSystem::stateNames() const {
  std::vector<String> result;
  const auto append = [&](const MMCStateSpaceBlock &block,
                          const String &prefix) {
    for (const auto &name : block.stateNames())
      result.push_back(prefix + name);
  };
  append(mActive, "active.");
  append(mReactive, "reactive.");
  append(mOcc, "occ.");
  append(mCcc, "ccc.");
  append(mZcc, "zcc.");
  append(mEnergy, "energy.");
  append(mPll, "pll.");
  append(mVacFilter, "vac.");
  append(mPFilter, "p.");
  append(mQFilter, "q.");
  append(mVdcFilter, "vdc.");
  append(mVacMagnitudeFilter, "vac_mag.");
  append(mDelay, "modulation.");
  return result;
}

Matrix MMCControllerSystem::rotateDq(const Matrix &dq, Real angle) {
  if (dq.rows() != 2 || dq.cols() != 1)
    throw std::invalid_argument("rotateDq expects 2x1.");
  const Real c = std::cos(angle);
  const Real s = std::sin(angle);
  Matrix result(2, 1);
  result(0, 0) = c * dq(0, 0) - s * dq(1, 0);
  result(1, 0) = s * dq(0, 0) + c * dq(1, 0);
  return result;
}
Real MMCControllerSystem::clamp(Real value, Real lower, Real upper) {
  return std::max(lower, std::min(value, upper));
}
Real MMCControllerSystem::conditionalIntegratorError(
    Real error, Real unsaturatedOutput, Real saturatedOutput,
    Real outputPerIntegratorSign) {
  const Real saturationDirection = unsaturatedOutput - saturatedOutput;
  const Real driveDirection = outputPerIntegratorSign * error;
  return saturationDirection * driveDirection > 0.0 ? 0.0 : error;
}

void MMCControllerSystem::evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                                   Matrix &y) const {
  validateDimensions(x, u);
  dx = Matrix::Zero(stateSize(), 1);
  y = Matrix::Zero(outputSize(), 1);

  Matrix vGrid(2, 1);
  vGrid << u(VGridD, 0), u(VGridQ, 0);
  const Real pllAngle = mPllEnabled ? x(mSlices.pll.first + 1, 0) : 0.0;
  const Matrix vControlRaw = rotateDq(vGrid, -pllAngle);

  Matrix dxPll, yPll;
  Matrix pllInput = Matrix::Constant(1, 1, vControlRaw(1, 0));
  mPll.evaluate(mSlices.pll.get(x), pllInput, dxPll, yPll);
  mSlices.pll.set(dx, dxPll);
  const Real deltaOmega = yPll(0, 0);
  const Real omegaControl = mOmegaN + deltaOmega;

  Matrix dxVac, yVac;
  mVacFilter.evaluate(mSlices.vacFilter.get(x), vControlRaw, dxVac, yVac);
  mSlices.vacFilter.set(dx, dxVac);
  const Real vControlD = yVac(0, 0);
  const Real vControlQ = yVac(1, 0);

  Matrix one(1, 1), dxBlock, yBlock;
  one(0, 0) = u(ActivePower, 0);
  mPFilter.evaluate(mSlices.pFilter.get(x), one, dxBlock, yBlock);
  mSlices.pFilter.set(dx, dxBlock);
  const Real p = yBlock(0, 0);

  one(0, 0) = u(ReactivePower, 0);
  mQFilter.evaluate(mSlices.qFilter.get(x), one, dxBlock, yBlock);
  mSlices.qFilter.set(dx, dxBlock);
  const Real q = yBlock(0, 0);

  one(0, 0) = u(DcVoltage, 0);
  mVdcFilter.evaluate(mSlices.vdcFilter.get(x), one, dxBlock, yBlock);
  mSlices.vdcFilter.set(dx, dxBlock);
  const Real vdc = yBlock(0, 0);

  const Real vacRaw = 1.5 * std::hypot(vControlD, vControlQ);
  one(0, 0) = vacRaw;
  mVacMagnitudeFilter.evaluate(mSlices.vacMagnitudeFilter.get(x), one, dxBlock,
                               yBlock);
  mSlices.vacMagnitudeFilter.set(dx, dxBlock);
  const Real vac = yBlock(0, 0);

  const Real heldIdReference = u(HeldActiveCurrentReference, 0);

  Matrix activeInput(4, 1);
  activeInput << p, vdc, vControlD, heldIdReference;
  Matrix dxActive, yActive;
  mActive.evaluate(mSlices.active.get(x), activeInput, dxActive, yActive);

  Matrix reactiveInput(2, 1);
  reactiveInput << q, vac;
  Matrix dxReactive, yReactive;
  mReactive.evaluate(mSlices.reactive.get(x), reactiveInput, dxReactive,
                     yReactive);

  Real idRefNetwork = yActive(0, 0);
  Real iqRefNetwork = yReactive(0, 0);
  const Real rawIdRef = idRefNetwork;
  const Real rawIqRef = iqRefNetwork;
  const Real currentMagnitude = std::hypot(idRefNetwork, iqRefNetwork);
  if (currentMagnitude > mMaximumAcCurrent) {
    const Real scale = mMaximumAcCurrent / currentMagnitude;
    idRefNetwork *= scale;
    iqRefNetwork *= scale;
  }
  if (mActiveMode == MMCActiveMode::ActivePower ||
      mActiveMode == MMCActiveMode::DcVoltage)
    dxActive(0, 0) = conditionalIntegratorError(yActive(1, 0), rawIdRef,
                                                idRefNetwork, yActive(2, 0));
  if (mReactiveMode == MMCReactiveMode::ReactivePower ||
      mReactiveMode == MMCReactiveMode::AcVoltage)
    dxReactive(0, 0) = conditionalIntegratorError(
        yReactive(1, 0), rawIqRef, iqRefNetwork, yReactive(2, 0));
  mSlices.active.set(dx, dxActive);
  mSlices.reactive.set(dx, dxReactive);

  Matrix iDeltaNetwork(2, 1);
  iDeltaNetwork << u(IDeltaD, 0), u(IDeltaQ, 0);
  const Matrix iDeltaControl = rotateDq(iDeltaNetwork, -pllAngle);
  Matrix iDeltaRefNetwork(2, 1);
  iDeltaRefNetwork << idRefNetwork, iqRefNetwork;
  const Matrix iDeltaRefControl = rotateDq(iDeltaRefNetwork, -pllAngle);

  Matrix occInput(9, 1);
  occInput << iDeltaRefControl(0, 0), iDeltaRefControl(1, 0),
      iDeltaControl(0, 0), iDeltaControl(1, 0), vControlD, vControlQ,
      omegaControl, mEquivalentAcResistance, mEquivalentAcInductance;
  Matrix dxOcc, yOcc;
  mOcc.evaluate(mSlices.occ.get(x), occInput, dxOcc, yOcc);
  Matrix vMDeltaControl(2, 1);
  vMDeltaControl << yOcc(0, 0), yOcc(1, 0);
  const Matrix vMDeltaInternal = rotateDq(vMDeltaControl, pllAngle);
  Matrix vMDeltaApplied = vMDeltaInternal;
  if (mControlSource == ControlSource::ExternalDifferentialVoltage ||
      mControlSource == ControlSource::ExternalFullConverterVoltage) {
    vMDeltaApplied = mExternalDifferentialVoltage;
    // External differential voltage bypasses the internal OCC command.
    dxOcc.setZero();
  }

  Matrix iSigmaNetwork(2, 1);
  iSigmaNetwork << u(ISigmaD, 0), u(ISigmaQ, 0);
  const Matrix iSigmaControl = rotateDq(iSigmaNetwork, -2.0 * pllAngle);
  Matrix iSigmaRefNetwork(2, 1);
  iSigmaRefNetwork << mISigmaDReference, mISigmaQReference;
  const Matrix iSigmaRefControl = rotateDq(iSigmaRefNetwork, -2.0 * pllAngle);

  Matrix cccInput(6, 1);
  cccInput << iSigmaRefControl(0, 0), iSigmaRefControl(1, 0),
      iSigmaControl(0, 0), iSigmaControl(1, 0), omegaControl, mArmInductance;
  Matrix dxCcc, yCcc;
  mCcc.evaluate(mSlices.ccc.get(x), cccInput, dxCcc, yCcc);
  Matrix vMSigmaControl(2, 1);
  vMSigmaControl << yCcc(0, 0), yCcc(1, 0);
  const Matrix vMSigmaInternal = rotateDq(vMSigmaControl, 2.0 * pllAngle);
  Matrix vMSigmaApplied = vMSigmaInternal;

  Matrix energyInput(3, 1);
  energyInput << u(StoredEnergy, 0), p, vdc;
  Matrix dxEnergy, yEnergy;
  mEnergy.evaluate(mSlices.energy.get(x), energyInput, dxEnergy, yEnergy);
  Real iSigmaZRefRaw = mEnergyEnabled ? yEnergy(0, 0) : mISigmaZReference;
  const Real iSigmaZRef = clamp(iSigmaZRefRaw, -mMaximumCirculatingCurrent,
                                mMaximumCirculatingCurrent);
  if (mEnergyEnabled)
    dxEnergy(0, 0) = conditionalIntegratorError(yEnergy(1, 0), iSigmaZRefRaw,
                                                iSigmaZRef, 1.0);
  mSlices.energy.set(dx, dxEnergy);

  Matrix zccInput(3, 1);
  zccInput << iSigmaZRef, u(ISigmaZ, 0), vdc;
  Matrix dxZcc, yZcc;
  mZcc.evaluate(mSlices.zcc.get(x), zccInput, dxZcc, yZcc);
  Real vMSigmaZApplied = yZcc(0, 0);

  if (mControlSource == ControlSource::ExternalFullConverterVoltage) {
    vMSigmaApplied(0, 0) = mExternalCommonModeVoltage(0, 0);
    vMSigmaApplied(1, 0) = mExternalCommonModeVoltage(1, 0);
    vMSigmaZApplied = mExternalCommonModeVoltage(2, 0);
    dxCcc.setZero();
    dxZcc.setZero();
  }
  mSlices.ccc.set(dx, dxCcc);
  mSlices.zcc.set(dx, dxZcc);

  const Real vdcReg = std::abs(vdc) >= 1.0 ? vdc : (vdc >= 0.0 ? 1.0 : -1.0);
  Real mDeltaD = -2.0 * vMDeltaApplied(0, 0) / vdcReg;
  Real mDeltaQ = -2.0 * vMDeltaApplied(1, 0) / vdcReg;
  Real mSigmaD = 2.0 * vMSigmaApplied(0, 0) / vdcReg;
  Real mSigmaQ = 2.0 * vMSigmaApplied(1, 0) / vdcReg;
  Real mSigmaZ = 2.0 * vMSigmaZApplied / vdcReg;

  const Real differentialMagnitude = std::hypot(mDeltaD, mDeltaQ);
  if (differentialMagnitude > mMaximumModulationMagnitude) {
    const Real scale = mMaximumModulationMagnitude / differentialMagnitude;
    if (mControlSource == ControlSource::InternalControllers) {
      Matrix saturationExcessNetwork(2, 1);
      saturationExcessNetwork << (1.0 - scale) * vMDeltaApplied(0, 0),
          (1.0 - scale) * vMDeltaApplied(1, 0);
      Matrix integratorDriveControl(2, 1);
      integratorDriveControl << yOcc(2, 0), yOcc(3, 0);
      const Matrix integratorDriveNetwork =
          rotateDq(integratorDriveControl, pllAngle);
      if (saturationExcessNetwork.cwiseProduct(integratorDriveNetwork).sum() >
          0.0)
        dxOcc.setZero();
    }
    mDeltaD *= scale;
    mDeltaQ *= scale;
  }
  mSlices.occ.set(dx, dxOcc);

  const Real sigmaMagnitude = std::hypot(mSigmaD, mSigmaQ);
  if (sigmaMagnitude > mMaximumModulationMagnitude) {
    const Real scale = mMaximumModulationMagnitude / sigmaMagnitude;
    mSigmaD *= scale;
    mSigmaQ *= scale;
  }
  mSigmaZ =
      clamp(mSigmaZ, -mMaximumModulationMagnitude, mMaximumModulationMagnitude);

  Matrix modulationInput(5, 1);
  modulationInput << mDeltaD, mDeltaQ, mSigmaD, mSigmaQ, mSigmaZ;
  Matrix dxDelay, delayedModulation;
  mDelay.evaluate(mSlices.delay.get(x), modulationInput, dxDelay,
                  delayedModulation);
  mSlices.delay.set(dx, dxDelay);

  y(MDeltaD, 0) = delayedModulation(0, 0);
  y(MDeltaQ, 0) = delayedModulation(1, 0);
  y(MSigmaD, 0) = delayedModulation(2, 0);
  y(MSigmaQ, 0) = delayedModulation(3, 0);
  y(MSigmaZ, 0) = delayedModulation(4, 0);
  y(IDeltaDReference, 0) = idRefNetwork;
  y(IDeltaQReference, 0) = iqRefNetwork;
  y(ISigmaZReference, 0) = iSigmaZRef;
  y(DeltaOmega, 0) = deltaOmega;
  y(PllAngle, 0) = pllAngle;
  y(VControlD, 0) = vControlD;
  y(VControlQ, 0) = vControlQ;
  y(FilteredActivePower, 0) = p;
  y(FilteredReactivePower, 0) = q;
  y(FilteredDcVoltage, 0) = vdc;
  y(FilteredAcMagnitude, 0) = vac;
  y(FeedforwardFilteredDaxisVoltageOutput, 0) =
      u(FeedforwardFilteredDaxisVoltage, 0);
  y(VMDeltaDCommand, 0) = vMDeltaApplied(0, 0);
  y(VMDeltaQCommand, 0) = vMDeltaApplied(1, 0);
  y(VMSigmaDCommand, 0) = vMSigmaApplied(0, 0);
  y(VMSigmaQCommand, 0) = vMSigmaApplied(1, 0);
  y(VMSigmaZCommand, 0) = vMSigmaZApplied;
  y(InternalVMDeltaDCommand, 0) = vMDeltaInternal(0, 0);
  y(InternalVMDeltaQCommand, 0) = vMDeltaInternal(1, 0);
}

Matrix MMCControllerSystem::initializeState(
    const Matrix &input, const Matrix &initialUndelayedModulation) const {
  if (input.rows() != InputCount || input.cols() != 1 || !input.allFinite())
    throw std::invalid_argument(
        "Controller initialization input must be finite.");
  if (initialUndelayedModulation.rows() != 5 ||
      initialUndelayedModulation.cols() != 1 ||
      !initialUndelayedModulation.allFinite())
    throw std::invalid_argument("Initial modulation must be finite 5x1.");

  Matrix x = Matrix::Zero(stateSize(), 1);
  // PLL aligned with the nominal network frame.
  x(mSlices.pll.first + 0, 0) = 0.0;
  x(mSlices.pll.first + 1, 0) = 0.0;

  Matrix vacState(2, 1);
  vacState << input(VGridD, 0), input(VGridQ, 0);
  mSlices.vacFilter.set(x, vacState);

  const auto initializeSecondOrder = [&](const MMCStateSlice &slice,
                                         Real value) {
    Matrix local(2, 1);
    local << value, value;
    slice.set(x, local);
  };
  initializeSecondOrder(mSlices.pFilter, input(ActivePower, 0));
  initializeSecondOrder(mSlices.qFilter, input(ReactivePower, 0));
  initializeSecondOrder(mSlices.vdcFilter, input(DcVoltage, 0));
  initializeSecondOrder(mSlices.vacMagnitudeFilter,
                        1.5 * std::hypot(input(VGridD, 0), input(VGridQ, 0)));

  // Preload each controller from its own local affine state-space model. This
  // keeps gain knowledge inside the signal class while giving the enclosing
  // MMC an exact hot-start estimate for all PI states.
  const auto stateForDesiredOutput =
      [](const MMCStateSpaceBlock &block, const Matrix &blockInput,
         const Matrix &desiredOutput, UInt outputRows) {
        const Matrix zero = Matrix::Zero(block.stateSize(), 1);
        const MMCLinearization lin = block.getStateSpaceModel(zero, blockInput);
        const Matrix c = lin.C.topRows(outputRows);
        const Matrix rhs = desiredOutput.topRows(outputRows) -
                           lin.D.topRows(outputRows) * blockInput -
                           lin.F.topRows(outputRows);
        Matrix local = c.colPivHouseholderQr().solve(rhs);
        if (!local.allFinite())
          local = zero;
        return local;
      };

  const Real held = input(HeldActiveCurrentReference, 0);
  Matrix activeInput(4, 1);
  activeInput << input(ActivePower, 0), input(DcVoltage, 0), input(VGridD, 0),
      held;
  Matrix desiredActive(1, 1);
  desiredActive(0, 0) = input(IDeltaD, 0);
  mSlices.active.set(
      x, stateForDesiredOutput(mActive, activeInput, desiredActive, 1));

  Matrix reactiveInput(2, 1);
  reactiveInput << input(ReactivePower, 0),
      1.5 * std::hypot(input(VGridD, 0), input(VGridQ, 0));
  Matrix desiredReactive(1, 1);
  desiredReactive(0, 0) = input(IDeltaQ, 0);
  mSlices.reactive.set(
      x, stateForDesiredOutput(mReactive, reactiveInput, desiredReactive, 1));

  Matrix occInput(9, 1);
  occInput << input(IDeltaD, 0), input(IDeltaQ, 0), input(IDeltaD, 0),
      input(IDeltaQ, 0), input(VGridD, 0), input(VGridQ, 0), mOmegaN,
      mEquivalentAcResistance, mEquivalentAcInductance;
  Matrix desiredOcc(2, 1);
  desiredOcc << -0.5 * input(DcVoltage, 0) * initialUndelayedModulation(0, 0),
      -0.5 * input(DcVoltage, 0) * initialUndelayedModulation(1, 0);
  mSlices.occ.set(x, stateForDesiredOutput(mOcc, occInput, desiredOcc, 2));

  Matrix cccInput(6, 1);
  cccInput << input(ISigmaD, 0), input(ISigmaQ, 0), input(ISigmaD, 0),
      input(ISigmaQ, 0), mOmegaN, mArmInductance;
  Matrix desiredCcc(2, 1);
  desiredCcc << -mArmResistance * input(ISigmaD, 0) -
                    2.0 * mOmegaN * mArmInductance * input(ISigmaQ, 0),
      -mArmResistance * input(ISigmaQ, 0) +
          2.0 * mOmegaN * mArmInductance * input(ISigmaD, 0);
  mSlices.ccc.set(x, stateForDesiredOutput(mCcc, cccInput, desiredCcc, 2));

  Matrix energyInput(3, 1);
  energyInput << input(StoredEnergy, 0), input(ActivePower, 0),
      input(DcVoltage, 0);
  Matrix desiredEnergy(1, 1);
  desiredEnergy(0, 0) = input(ISigmaZ, 0);
  mSlices.energy.set(
      x, stateForDesiredOutput(mEnergy, energyInput, desiredEnergy, 1));

  Matrix zccInput(3, 1);
  zccInput << input(ISigmaZ, 0), input(ISigmaZ, 0), input(DcVoltage, 0);
  Matrix desiredZcc(1, 1);
  desiredZcc(0, 0) =
      input(DcVoltage, 0) / 2.0 - mArmResistance * input(ISigmaZ, 0);
  mSlices.zcc.set(x, stateForDesiredOutput(mZcc, zccInput, desiredZcc, 1));

  if (mDelayEnabled) {
    // For the controllable-canonical Padé realization, solve the local block
    // equilibrium rather than duplicating its coefficients here.
    const Matrix zero = Matrix::Zero(mDelay.stateSize(), 1);
    const MMCLinearization lin =
        mDelay.getStateSpaceModel(zero, initialUndelayedModulation);
    Matrix equilibrium = lin.A.colPivHouseholderQr().solve(
        -(lin.B * initialUndelayedModulation + lin.E));
    if (equilibrium.allFinite())
      mSlices.delay.set(x, equilibrium);
  }
  return x;
}

void MMCControllerSystem::evaluateStateDerivative(const Matrix &x,
                                                  const Matrix &u,
                                                  Matrix &dx) const {
  Matrix output = Matrix::Zero(outputSize(), 1);
  evaluate(x, u, dx, output);
}

void MMCControllerSystem::evaluateOutput(const Matrix &x, const Matrix &u,
                                         Matrix &y) const {
  Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
  evaluate(x, u, stateDerivative, y);
}

void MMCControllerSystem::calculateNumericalJacobians(
    const Matrix &x, const Matrix &u, Matrix &A, Matrix &B, Matrix &C,
    Matrix &D, Real relativeStep, Real absoluteStep) const {
  calculateNumericalJacobiansGeneric(x, u, A, B, C, D, relativeStep,
                                     absoluteStep);
}

void MMCControllerSystem::buildStateSpaceModel(const Matrix &x, const Matrix &u,
                                               Matrix &A, Matrix &B, Matrix &C,
                                               Matrix &D, Matrix &E, Matrix &F,
                                               Real relativeStep,
                                               Real absoluteStep) const {
  buildStateSpaceModelGeneric(x, u, A, B, C, D, E, F, relativeStep,
                              absoluteStep);
}

MMCLinearization
MMCControllerSystem::getStateSpaceModel(const Matrix &x, const Matrix &u,
                                        Real relativeStep,
                                        Real absoluteStep) const {
  return getStateSpaceModelGeneric(x, u, relativeStep, absoluteStep);
}

MMCSparseLinearization MMCControllerSystem::getSparseStateSpaceModel(
    const Matrix &x, const Matrix &u, Real relativeStep, Real absoluteStep,
    Real sparseTolerance) const {
  return getSparseStateSpaceModelGeneric(x, u, relativeStep, absoluteStep,
                                         sparseTolerance);
}
