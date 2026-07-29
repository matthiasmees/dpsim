// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <dpsim-models/Definitions.h>
#include <dpsim-models/Signal/MMCControllerSystem.h>
#include <dpsim-models/Signal/MMCSampledPowerFeedforward.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

/// Reusable nonlinear averaged MMC model without MNA ownership.
///
/// Continuous state layout:
///   13 plant/network-frame states + 30 continuous controller states = 43.
///
/// The sampled active-power feedforward controller remains a separate hybrid
/// signal class. Its three discrete states are held outside the continuous SSN
/// state vector and cause an affine-model rebuild only at sample instants.
/// This preserves the behavior of the previously validated standalone MMC.
///
/// Input:  [Va,Vb,Vc,Vdc+,Vdc-]
/// Output: [Ia,Ib,Ic,Idc+,Idc-], positive into the component.
class MMC_ModularStateSpaceModel {
public:
  using ControlSource = Signal::MMCControllerSystem::ControlSource;
  enum Input : UInt { Va = 0, Vb, Vc, Vdcp, Vdcn, InputCount };
  enum Output : UInt { Ia = 0, Ib, Ic, Idcp, Idcn, OutputCount };

  enum PlantState : UInt {
    IDeltaD = 0,
    IDeltaQ,
    ISigmaZ,
    ISigmaD,
    ISigmaQ,
    VCDeltaD,
    VCDeltaQ,
    VCDeltaZd,
    VCDeltaZq,
    VCSigmaD,
    VCSigmaQ,
    VCSigmaZ,
    GridAngle,
    PlantStateCount
  };

  MMC_ModularStateSpaceModel();

  void setParameters(Real nominalFrequency, Real nominalAcVoltage,
                     Real nominalDcVoltage, Real armInductance,
                     Real armResistance, Real submoduleCapacitance,
                     UInt numberOfSubmodules, Real reactorInductance,
                     Real reactorResistance);

  void setInitialAngle(Real angle);
  void setInitialOperatingPoint(Real activePower, Real reactivePower);
  void clearInitialOperatingPoint();

  void setPLL(Real kp, Real ki, Bool enabled);
  void setOutputCurrentController(Real kp, Real ki);
  void setCirculatingCurrentController(Real kp, Real ki);
  void setZeroSequenceCurrentController(Real kp, Real ki);
  void setEnergyController(Real kp, Real ki, Bool enabled);

  void setActivePowerControl(Real reference, Real kp, Real ki);
  void setActivePowerFeedforwardControl(Real reference, Real cutoffFrequency,
                                        Real sampleTime,
                                        Real minimumDaxisVoltage = 0.0);
  void setActivePowerReference(Real reference);
  void setDcVoltageControl(Real reference, Real kp, Real ki);
  void setDcDroopControl(Real activePowerReference, Real dcVoltageReference,
                         Real droopGain);
  void setActiveControlOpenLoop(Real currentReference);

  void setReactivePowerControl(Real reference, Real kp, Real ki);
  void setAcVoltageControl(Real reference, Real kp, Real ki);
  void setReactiveControlOpenLoop(Real currentReference);

  void setCirculatingCurrentReferences(Real dReference, Real qReference,
                                       Real zReference);
  void setMeasurementFilters(Real acVoltageDqTimeConstant,
                             Real activePowerTimeConstant,
                             Real reactivePowerTimeConstant,
                             Real dcVoltageTimeConstant,
                             Real acVoltageMagnitudeTimeConstant);
  void setModulationDelay(Real timeDelay, UInt padeOrder = 2);
  void setLimits(Real maximumAcCurrent, Real maximumCirculatingCurrent,
                 Real maximumModulationMagnitude);
  void setControlSource(ControlSource source);
  void setExternalDifferentialVoltageCommand(Real dVolts, Real qVolts);
  void setExternalCommonModeVoltageCommand(Real dVolts, Real qVolts,
                                           Real zVolts);
  ControlSource controlSource() const;
  void setNumericalLinearizationParameters(Real relativeStep,
                                           Real absoluteStep);
  void setSparseLinearizationTolerance(Real tolerance);
  void setStructuredLinearization(Bool enabled);
  void setOperatingPointInitialization(Bool enabled, UInt maximumIterations,
                                       Real normalizedTolerance);

  Bool parametersSet() const;
  std::uint64_t configurationRevision() const;
  Real nominalFrequency() const;
  Real nominalDcVoltage() const;
  UInt stateSize() const;
  UInt inputSize() const;
  UInt outputSize() const;
  std::vector<String> stateNames() const;

  void evaluate(const Matrix &state, const Matrix &input,
                Matrix &stateDerivative, Matrix &output) const;
  void evaluateStateDerivative(const Matrix &state, const Matrix &input,
                               Matrix &stateDerivative) const;
  void evaluateOutput(const Matrix &state, const Matrix &input,
                      Matrix &output) const;

  void calculateNumericalJacobians(const Matrix &state, const Matrix &input,
                                   Matrix &A, Matrix &B, Matrix &C,
                                   Matrix &D) const;
  void buildStateSpaceModel(const Matrix &state, const Matrix &input, Matrix &A,
                            Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                            Matrix &F) const;
  Signal::MMCLinearization getStateSpaceModel(const Matrix &state,
                                              const Matrix &input) const;
  Signal::MMCSparseLinearization
  getSparseStateSpaceModel(const Matrix &state, const Matrix &input) const;

  void linearize(const Matrix &state, const Matrix &input, Matrix &A, Matrix &B,
                 Matrix &C, Matrix &D, Matrix &E, Matrix &F) const;

  Matrix initializeState(const Matrix &input);

  /// Advance the separate sampled feedforward controller. Returns true when a
  /// new zero-order-held reference was produced and the aggregate affine model
  /// therefore has to be rebuilt immediately.
  Bool advanceSampledControllers(const Matrix &state, const Matrix &input,
                                 Real timeStep);

  Real calculateStoredEnergy(const Matrix &state) const;
  Real activePower(const Matrix &state, const Matrix &input) const;
  Real reactivePower(const Matrix &state, const Matrix &input) const;
  Real dcVoltage(const Matrix &input) const;
  Real dcCurrent(const Matrix &state) const;
  Real powerBalanceError(const Matrix &state, const Matrix &input) const;
  Matrix controllerOutput(const Matrix &state, const Matrix &input) const;

  Real feedforwardFilteredDaxisVoltage() const;
  Real heldActiveCurrentReference() const;

private:
  static Matrix abcToDq(const Matrix &abc, Real theta);
  static Matrix dqToAbc(Real d, Real q, Real theta);
  static Real regularizeSigned(Real value, Real minimumMagnitude);
  void validateConfigured() const;
  void validateStateInput(const Matrix &state, const Matrix &input) const;
  Matrix makeControllerInput(const Matrix &state, const Matrix &input) const;
  Matrix initializeAnalyticalState(const Matrix &input);
  std::vector<UInt> activeLinearizationStateIndices() const;
  std::vector<UInt> equilibriumStateIndices() const;
  Bool solveOperatingPoint(Matrix &state, const Matrix &input,
                           Real &normalizedResidual,
                           Real &absoluteResidual) const;

  void calculateFullNumericalJacobians(const Matrix &state, const Matrix &input,
                                       Matrix &A, Matrix &B, Matrix &C,
                                       Matrix &D) const;
  void calculateStructuredNumericalJacobians(const Matrix &state,
                                             const Matrix &input, Matrix &A,
                                             Matrix &B, Matrix &C,
                                             Matrix &D) const;

  void touchConfiguration();
  UInt controllerOffset() const;

  Bool mParametersSet;
  std::uint64_t mConfigurationRevision;
  Real mNominalFrequency;
  Real mOmegaN;
  Real mNominalAcVoltage;
  Real mNominalDcVoltage;
  Real mArmInductance;
  Real mArmResistance;
  Real mSubmoduleCapacitance;
  UInt mNumberOfSubmodules;
  Real mReactorInductance;
  Real mReactorResistance;

  Real mInitialAngle;
  Bool mInitialOperatingPointEnabled;
  Real mInitialActivePower;
  Real mInitialReactivePower;
  Matrix mExternalDifferentialVoltage;
  Matrix mExternalCommonModeVoltage;

  Real mMinimumDcVoltage;
  Real mJacobianRelativeStep;
  Real mJacobianAbsoluteStep;
  Real mSparseLinearizationTolerance;
  Bool mStructuredLinearization;
  Bool mOperatingPointInitializationEnabled;
  UInt mOperatingPointMaximumIterations;
  Real mOperatingPointNormalizedTolerance;

  UInt mFeedforwardStepCounter;
  Bool mFeedforwardInitialized;
  Signal::MMCSampledPowerFeedforward mSampledFeedforward;
  Matrix mSampledFeedforwardState;

  Signal::MMCControllerSystem mControllers;
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
