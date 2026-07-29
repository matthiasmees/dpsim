// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <cstdint>
#include <memory>

#include <dpsim-models/EMT/EMT_Ph3_MMC_ModularStateSpaceModel.h>
#include <dpsim-models/EMT/EMT_VTypeVariableSSNComp.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

/// SSN/MNA wrapper around MMC_ModularStateSpaceModel.
///
/// All controller states and equations are supplied by separate Signal MMC
/// state-space classes. This wrapper owns only terminal mapping, MNA stamping,
/// logging attributes, and the final aggregate affine linearization.
class MMC_Modular : public VTypeVariableSSNComp,
                    public SharedFactory<MMC_Modular> {
public:
  using SharedFactory<MMC_Modular>::make;
  using Ptr = std::shared_ptr<MMC_Modular>;
  using ControlSource = MMC_ModularStateSpaceModel::ControlSource;

  MMC_Modular(String uid, String name,
              Logger::Level logLevel = Logger::Level::off);
  MMC_Modular(String name, Logger::Level logLevel = Logger::Level::off)
      : MMC_Modular(name, name, logLevel) {}

  void setParameters(Real nominalFrequency, Real nominalAcVoltage,
                     Real nominalDcVoltage, Real armInductance,
                     Real armResistance, Real submoduleCapacitance,
                     UInt numberOfSubmodules, Real reactorInductance,
                     Real reactorResistance);
  void setInitialAngle(Real angle);
  void setInitialOperatingPoint(Real activePower, Real reactivePower);
  void setPLL(Real kp, Real ki, Bool enabled = true);
  void setOutputCurrentController(Real kp, Real ki);
  void setCirculatingCurrentController(Real kp, Real ki);
  void setZeroSequenceCurrentController(Real kp, Real ki);
  void setEnergyController(Real kp, Real ki, Bool enabled);
  void setActivePowerControl(Real reference, Real kp, Real ki);
  void setActivePowerFeedforwardControl(Real reference, Real cutoffFrequency,
                                        Real sampleTime,
                                        Real minimumDaxisVoltage = 0.0);
  void setActivePowerFeedforwardReference(Real reference);
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
  void setNumericalLinearizationParameters(Real relativeStep,
                                           Real absoluteStep);
  void setSparseLinearizationTolerance(Real tolerance);
  void setStructuredLinearization(Bool enabled = true);
  void setRuntimeDiagnostics(Bool enabled);
  void setLinearizationUpdateInterval(UInt intervalSteps);
  void setLinearizationTimeStep(Real linearizationTimeStep);
  void setOperatingPointInitialization(Bool enabled, UInt maximumIterations,
                                       Real normalizedTolerance);

  MMC_ModularStateSpaceModel &model();
  const MMC_ModularStateSpaceModel &model() const;

  std::vector<String> getLocalStateNames() const override;

  Matrix getState() const;
  Matrix getStateDerivative() const;
  Matrix getInterfaceVoltage() const;
  Matrix getInterfaceCurrent() const;
  void getLocalLinearization(Matrix &A, Matrix &B, Matrix &C, Matrix &D) const;
  void evaluateStateDerivative(const Matrix &x, const Matrix &u,
                               Matrix &dx) const;
  void evaluateOutput(const Matrix &x, const Matrix &u, Matrix &y) const;
  void calculateNumericalJacobians(const Matrix &x, const Matrix &u, Matrix &A,
                                   Matrix &B, Matrix &C, Matrix &D) const;
  void buildStateSpaceModel(const Matrix &x, const Matrix &u, Matrix &A,
                            Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                            Matrix &F) const;
  Signal::MMCLinearization getStateSpaceModel(const Matrix &x,
                                              const Matrix &u) const;
  Signal::MMCSparseLinearization
  getSparseStateSpaceModel(const Matrix &x, const Matrix &u) const;
  Signal::MMCLinearization getStateSpaceModel() const;
  Signal::MMCSparseLinearization getSparseStateSpaceModel() const;

  std::uint64_t linearizationRebuildCount() const;
  std::uint64_t linearizationReuseCount() const;
  void resetLinearizationCounters();

  Attribute<Matrix>::Ptr acTerminalVoltageAttribute() const;
  Attribute<Matrix>::Ptr acTerminalCurrentAttribute() const;
  Attribute<Real>::Ptr dcPositiveVoltageAttribute() const;
  Attribute<Real>::Ptr dcNegativeVoltageAttribute() const;
  Attribute<Real>::Ptr dcVoltageAttribute() const;
  Attribute<Real>::Ptr dcCurrentAttribute() const;
  Attribute<Real>::Ptr activePowerAttribute() const;
  Attribute<Real>::Ptr reactivePowerAttribute() const;
  Attribute<Real>::Ptr storedEnergyAttribute() const;
  Attribute<Matrix>::Ptr interfaceVoltageAttribute() const;
  Attribute<Matrix>::Ptr interfaceCurrentAttribute() const;
  Attribute<Real>::Ptr filteredDaxisVoltageAttribute() const;
  Attribute<Real>::Ptr heldActiveCurrentReferenceAttribute() const;
  Attribute<Matrix>::Ptr appliedModulationAttribute() const;
  Attribute<Matrix>::Ptr differentialVoltageCommandAttribute() const;
  Attribute<Matrix>::Ptr commonModeVoltageCommandAttribute() const;
  Attribute<Matrix>::Ptr appliedDifferentialVoltageAttribute() const;
  Attribute<Matrix>::Ptr internalDifferentialVoltageAttribute() const;
  Attribute<Matrix>::Ptr appliedCommonModeVoltageAttribute() const;
  Attribute<Matrix>::Ptr externalDifferentialVoltageAttribute() const;
  Attribute<Matrix>::Ptr externalCommonModeVoltageAttribute() const;
  Attribute<Bool>::Ptr externalCommandActiveAttribute() const;

  void initializeFromNodesAndTerminals(Real frequency) override;
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;

protected:
  MatrixComp buildInitialInputFromNodes(Real frequency) override;
  Bool updateComponentParameters() override;
  void updateLogAttributes(const Matrix &u) const override;

  // Optional hook used by the variable-SSN scheduler on the development
  // branch. It remains a valid virtual override even without the keyword and
  // is harmless on branches where the hook is not present.
  void
  addHeldControlDependencies(AttributeBase::List &prevStepDependencies) const;

private:
  void validateTerminalArrangement() const;
  void validateConfigured() const;
  Bool shouldUpdateLinearization();
  void markLinearizationDirty();

  static constexpr UInt mInputSize = MMC_ModularStateSpaceModel::InputCount;
  static constexpr UInt mOutputSize = MMC_ModularStateSpaceModel::OutputCount;

  MMC_ModularStateSpaceModel mModel;
  Bool mModelConfigured;
  Bool mInitializationInProgress;
  UInt mLinearizationUpdateInterval;
  UInt mStepsSinceLinearization;
  Real mLinearizationTimeStep;
  Real mTimeSinceLinearization;
  Bool mLinearizationDirty;
  std::uint64_t mLastLinearizedConfigurationRevision;
  Real mStampZeroTolerance;
  Bool mRuntimeDiagnostics;
  std::uint64_t mLinearizationRebuildCount;
  std::uint64_t mLinearizationReuseCount;

  Attribute<Matrix>::Ptr mAcTerminalVoltage;
  Attribute<Matrix>::Ptr mAcTerminalCurrent;
  Attribute<Real>::Ptr mDcPositiveVoltage;
  Attribute<Real>::Ptr mDcNegativeVoltage;
  Attribute<Real>::Ptr mDcVoltage;
  Attribute<Real>::Ptr mDcCurrent;
  Attribute<Real>::Ptr mActivePower;
  Attribute<Real>::Ptr mReactivePower;
  Attribute<Real>::Ptr mStoredEnergy;
  Attribute<Real>::Ptr mFilteredActivePower;
  Attribute<Real>::Ptr mFilteredReactivePower;
  Attribute<Real>::Ptr mFilteredDcVoltage;
  Attribute<Real>::Ptr mFilteredDaxisVoltage;
  Attribute<Real>::Ptr mHeldActiveCurrentReference;
  Attribute<Real>::Ptr mPllFrequency;
  Attribute<Real>::Ptr mPllAngleDeviation;
  Attribute<Real>::Ptr mGridVoltageD;
  Attribute<Real>::Ptr mGridVoltageQ;
  Attribute<Real>::Ptr mControlVoltageD;
  Attribute<Real>::Ptr mControlVoltageQ;
  Attribute<Real>::Ptr mDeltaCurrentD;
  Attribute<Real>::Ptr mDeltaCurrentQ;
  Attribute<Real>::Ptr mSigmaCurrentZ;
  Attribute<Real>::Ptr mDeltaCurrentReferenceD;
  Attribute<Real>::Ptr mDeltaCurrentReferenceQ;
  Attribute<Real>::Ptr mSigmaCurrentReferenceZ;
  Attribute<Real>::Ptr mDcPower;
  Attribute<Real>::Ptr mPowerBalanceError;
  Attribute<Matrix>::Ptr mAppliedModulation;
  Attribute<Matrix>::Ptr mDifferentialVoltageCommand;
  Attribute<Matrix>::Ptr mInternalDifferentialVoltage;
  Attribute<Matrix>::Ptr mCommonModeVoltageCommand;
  Attribute<Matrix>::Ptr mExternalDifferentialVoltage;
  Attribute<Matrix>::Ptr mExternalCommonModeVoltage;
  Attribute<Bool>::Ptr mExternalCommandActive;
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
