// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/Signal/MMCActiveOuterController.h>
#include <dpsim-models/Signal/MMCCirculatingCurrentController.h>
#include <dpsim-models/Signal/MMCEnergyController.h>
#include <dpsim-models/Signal/MMCFirstOrderFilter.h>
#include <dpsim-models/Signal/MMCOutputCurrentController.h>
#include <dpsim-models/Signal/MMCPadeDelay.h>
#include <dpsim-models/Signal/MMCReactiveOuterController.h>
#include <dpsim-models/Signal/MMCSecondOrderFilter.h>
#include <dpsim-models/Signal/MMCSrfPll.h>
#include <dpsim-models/Signal/MMCZeroSequenceCurrentController.h>

namespace CPS {
namespace Signal {

/// Complete continuous-time MMC controller stack assembled from separate
/// state-space signal blocks.
///
/// The sampled active-power feedforward block is deliberately not part of this
/// continuous state vector. Its zero-order-held values are supplied as inputs
/// by MMC_ModularStateSpaceModel. This preserves the validated 43-state MMC model and
/// avoids modifying VTypeVariableSSNComp's integrated state outside the SSN
/// integration task.
class MMCControllerSystem final : public MMCStateSpaceBlock {
public:
  enum class ControlSource {
    InternalControllers,
    ExternalDifferentialVoltage,
    ExternalFullConverterVoltage
  };

  enum Input : UInt {
    VGridD = 0,
    VGridQ,
    ActivePower,
    ReactivePower,
    DcVoltage,
    IDeltaD,
    IDeltaQ,
    ISigmaD,
    ISigmaQ,
    ISigmaZ,
    StoredEnergy,
    FeedforwardFilteredDaxisVoltage,
    HeldActiveCurrentReference,
    InputCount
  };

  enum Output : UInt {
    MDeltaD = 0,
    MDeltaQ,
    MSigmaD,
    MSigmaQ,
    MSigmaZ,
    IDeltaDReference,
    IDeltaQReference,
    ISigmaZReference,
    DeltaOmega,
    PllAngle,
    VControlD,
    VControlQ,
    FilteredActivePower,
    FilteredReactivePower,
    FilteredDcVoltage,
    FilteredAcMagnitude,
    FeedforwardFilteredDaxisVoltageOutput,
    VMDeltaDCommand,
    VMDeltaQCommand,
    VMSigmaDCommand,
    VMSigmaQCommand,
    VMSigmaZCommand,
    InternalVMDeltaDCommand,
    InternalVMDeltaQCommand,
    OutputCount
  };

  MMCControllerSystem();

  void setElectricalParameters(Real nominalOmega, Real equivalentAcResistance,
                               Real equivalentAcInductance, Real armResistance,
                               Real armInductance, Real nominalDcVoltage,
                               Real energyReference);
  void setPLL(Real kp, Real ki, Bool enabled);
  void setMMCOutputCurrentController(Real kp, Real ki);
  void setMMCCirculatingCurrentController(Real kp, Real ki);
  void setMMCZeroSequenceCurrentController(Real kp, Real ki);
  void setMMCEnergyController(Real kp, Real ki, Bool enabled);

  void setActivePowerControl(Real reference, Real kp, Real ki);
  void setActivePowerFeedforwardControl(Real reference, Real cutoffFrequency,
                                        Real sampleTime,
                                        Real minimumDaxisVoltage);
  void setActivePowerReference(Real reference);
  void setDcVoltageControl(Real reference, Real kp, Real ki);
  void setDcDroopControl(Real activePowerReference, Real dcVoltageReference,
                         Real droopGain);
  void setActiveOpenLoop(Real currentReference);

  void setReactivePowerControl(Real reference, Real kp, Real ki);
  void setAcVoltageControl(Real reference, Real kp, Real ki);
  void setReactiveOpenLoop(Real currentReference);

  void setCirculatingCurrentReferences(Real dReference, Real qReference,
                                       Real zReference);
  void setMeasurementFilters(Real acVoltageDqTimeConstant,
                             Real activePowerTimeConstant,
                             Real reactivePowerTimeConstant,
                             Real dcVoltageTimeConstant,
                             Real acVoltageMagnitudeTimeConstant);
  void setModulationDelay(Real timeDelay, Bool enabled);
  void setLimits(Real maximumAcCurrent, Real maximumCirculatingCurrent,
                 Real maximumModulationMagnitude);

  void setControlSource(ControlSource source);
  void setExternalDifferentialVoltage(Real dVolts, Real qVolts);
  void setExternalCommonModeVoltage(Real dVolts, Real qVolts, Real zVolts);
  ControlSource controlSource() const;

  MMCActiveMode activeMode() const;
  MMCReactiveMode reactiveMode() const;
  Real activePowerReference() const;
  Real dcVoltageReference() const;
  Real reactivePowerReference() const;
  Real acVoltageReference() const;
  Real sigmaDReference() const;
  Real sigmaQReference() const;
  Real sigmaZReference() const;
  Bool energyControllerEnabled() const;
  Bool sampledFeedforwardEnabled() const;
  Real sampledFeedforwardPeriod() const;

  /// Local controller-state index used by the MMC operating-point solver.
  UInt activeIntegratorStateIndex() const;

  /// Continuous states that can influence the selected nonlinear model.
  std::vector<UInt> activeStateIndices() const;

  /// Continuous controller states that participate in the steady-state solve.
  /// PLL states are excluded because the plant is represented in the nominal
  /// network dq frame and the PLL deviation is initialized to zero.
  std::vector<UInt> equilibriumStateIndices() const;

  UInt stateSize() const override;
  UInt inputSize() const override;
  UInt outputSize() const override;
  std::vector<String> stateNames() const override;
  void evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                Matrix &y) const override;
  void evaluateStateDerivative(const Matrix &x, const Matrix &u,
                               Matrix &dx) const override;
  void evaluateOutput(const Matrix &x, const Matrix &u,
                      Matrix &y) const override;
  void calculateNumericalJacobians(const Matrix &x, const Matrix &u, Matrix &A,
                                   Matrix &B, Matrix &C, Matrix &D,
                                   Real relativeStep = 1e-6,
                                   Real absoluteStep = 1e-8) const override;
  void buildStateSpaceModel(const Matrix &x, const Matrix &u, Matrix &A,
                            Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                            Matrix &F, Real relativeStep = 1e-6,
                            Real absoluteStep = 1e-8) const override;
  MMCLinearization getStateSpaceModel(const Matrix &x, const Matrix &u,
                                      Real relativeStep = 1e-6,
                                      Real absoluteStep = 1e-8) const override;
  MMCSparseLinearization
  getSparseStateSpaceModel(const Matrix &x, const Matrix &u,
                           Real relativeStep = 1e-6, Real absoluteStep = 1e-8,
                           Real sparseTolerance = 1e-14) const override;

  /// Initializes all continuous controller/filter/delay states at a supplied
  /// operating point. The sampled feedforward values are supplied in input.
  Matrix initializeState(const Matrix &input,
                         const Matrix &initialUndelayedModulation) const;

private:
  struct Slices {
    MMCStateSlice active;
    MMCStateSlice reactive;
    MMCStateSlice occ;
    MMCStateSlice ccc;
    MMCStateSlice zcc;
    MMCStateSlice energy;
    MMCStateSlice pll;
    MMCStateSlice vacFilter;
    MMCStateSlice pFilter;
    MMCStateSlice qFilter;
    MMCStateSlice vdcFilter;
    MMCStateSlice vacMagnitudeFilter;
    MMCStateSlice delay;
    UInt total = 0;
  };

  static Slices makeSlices();
  static Matrix rotateDq(const Matrix &dq, Real angle);
  static Real clamp(Real value, Real lower, Real upper);
  static Real conditionalIntegratorError(Real error, Real unsaturatedOutput,
                                         Real saturatedOutput,
                                         Real outputPerIntegratorSign);

  Slices mSlices;

  MMCFirstOrderFilter mVacFilter;
  MMCSecondOrderFilter mPFilter;
  MMCSecondOrderFilter mQFilter;
  MMCSecondOrderFilter mVdcFilter;
  MMCSecondOrderFilter mVacMagnitudeFilter;
  MMCSrfPll mPll;
  MMCActiveOuterController mActive;
  MMCReactiveOuterController mReactive;
  MMCOutputCurrentController mOcc;
  MMCCirculatingCurrentController mCcc;
  MMCEnergyController mEnergy;
  MMCZeroSequenceCurrentController mZcc;
  MMCPadeDelay mDelay;

  Real mOmegaN;
  Real mEquivalentAcResistance;
  Real mEquivalentAcInductance;
  Real mArmResistance;
  Real mArmInductance;
  Real mNominalDcVoltage;
  Real mEnergyReference;
  Real mEnergyKp;
  Real mEnergyKi;

  Real mActivePowerReference;
  Real mDcVoltageReference;
  Real mDroopGain;
  Real mOpenLoopIDeltaDReference;
  Real mReactivePowerReference;
  Real mAcVoltageReference;
  Real mOpenLoopIDeltaQReference;
  Real mISigmaDReference;
  Real mISigmaQReference;
  Real mISigmaZReference;
  Real mFeedforwardSampleTime;

  MMCActiveMode mActiveMode;
  MMCReactiveMode mReactiveMode;
  Bool mEnergyEnabled;
  Bool mDelayEnabled;
  Bool mFeedforwardEnabled;
  Bool mPllEnabled;
  Bool mOccEnabled;
  Bool mCccEnabled;
  Bool mZccEnabled;
  Bool mVacFilterEnabled;
  Bool mPFilterEnabled;
  Bool mQFilterEnabled;
  Bool mVdcFilterEnabled;
  Bool mVacMagnitudeFilterEnabled;

  ControlSource mControlSource;
  Matrix mExternalDifferentialVoltage;
  Matrix mExternalCommonModeVoltage;

  Real mMaximumAcCurrent;
  Real mMaximumCirculatingCurrent;
  Real mMaximumModulationMagnitude;
};

} // namespace Signal
} // namespace CPS
