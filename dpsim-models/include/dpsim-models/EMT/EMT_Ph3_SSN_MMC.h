// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

// NEW FRAME-CORRECTED VERSION: separate nominal EMT angle and Harmony PLL deviation,
// consistent dq-frame transformations, continuous/discrete eigenvalue diagnostics.
#pragma once

#include <array>
#include <memory>
#include <vector>

#include <dpsim-models/EMT/EMT_VTypeVariableSSNComp.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

/// Averaged modular multilevel converter with one three-phase AC terminal and
/// two single-conductor DC pole terminals.
///
/// SSN interface ordering:
///   u = [v_a, v_b, v_c, v_dc+, v_dc-]^T
///   y = [i_a, i_b, i_c, i_dc+, i_dc-]^T
///
/// The MMC plant depends on v_dc = v_dc+ - v_dc-. Consequently, the ideal
/// converter is invariant to a common-mode shift of the two DC pole voltages.
///
/// Plant-state ordering follows the Harmony averaged MMC model:
///   [iDelta_d, iDelta_q, iSigma_z, iSigma_d, iSigma_q,
///    vCDelta_d, vCDelta_q, vCDelta_Zd, vCDelta_Zq,
///    vCSigma_d, vCSigma_q, vCSigma_z]
///
/// Controller integrator states are appended after the 12 plant states.
class SSN_MMC final : public VTypeVariableSSNComp {
public:
  using Ptr = std::shared_ptr<SSN_MMC>;

  enum class ActiveControlMode { OpenLoop, ActivePower, DcVoltage, DcDroop };
  enum class ReactiveControlMode { OpenLoop, ReactivePower, AcVoltage };
  enum class ControlSource { InternalControllers, ExternalDifferentialVoltage };

  static Ptr make(String name, Logger::Level logLevel = Logger::Level::off) {
    return std::make_shared<SSN_MMC>(name, name, logLevel);
  }

  static Ptr make(String uid, String name,
                  Logger::Level logLevel = Logger::Level::off) {
    return std::make_shared<SSN_MMC>(uid, name, logLevel);
  }

  SSN_MMC(String uid, String name, Logger::Level logLevel = Logger::Level::off);

  /// Physical MMC parameters.
  /// submoduleCapacitance is the capacitance of one submodule capacitor.
  void setParameters(Real nominalFrequency, Real nominalAcVoltage,
                     Real nominalDcVoltage, Real armInductance,
                     Real armResistance, Real submoduleCapacitance,
                     UInt numberOfSubmodules, Real reactorInductance,
                     Real reactorResistance);

  /// Configure the inner output-current controller.
  void setOutputCurrentController(Real kp, Real ki);

  /// Configure the circulating-current controller for iSigma_d/q.
  void setCirculatingCurrentController(Real kp, Real ki);

  /// Configure the zero-sequence circulating-current controller.
  void setZeroSequenceCurrentController(Real kp, Real ki);

  /// Optional total stored-energy controller, cascaded into iSigma_z_ref.
  void setEnergyController(Real kp, Real ki, Bool enabled = true);

  void setActivePowerControl(Real activePowerReference, Real kp, Real ki);
  void setDcVoltageControl(Real dcVoltageReference, Real kp, Real ki);
  void setDcDroopControl(Real activePowerReference, Real dcVoltageReference,
                         Real droopGain);
  void setActiveControlOpenLoop(Real iDeltaDReference = 0.0);

  void setReactivePowerControl(Real reactivePowerReference, Real kp, Real ki);
  void setAcVoltageControl(Real acVoltageReference, Real kp, Real ki);
  void setReactiveControlOpenLoop(Real iDeltaQReference = 0.0);

  void setCirculatingCurrentReferences(Real iSigmaDReference,
                                       Real iSigmaQReference,
                                       Real iSigmaZReference);

  /// Enable Harmony's synchronous-reference-frame PLL. The Harmony PLL
  /// states are the PI integrator and angle deviation. A separate nominal
  /// EMT angle state is used only for abc <-> nominal-dq conversion.
  void setPLL(Real kp, Real ki, Bool enabled = true);

  /// Measurement filters corresponding to Harmony's optional filters.
  /// acVoltageDqTimeConstant uses one first-order state per dq channel.
  /// The scalar measurements use cascaded second-order low-pass filters.
  /// A non-positive time constant disables the respective filter.
  void setMeasurementFilters(Real acVoltageDqTimeConstant,
                             Real activePowerTimeConstant,
                             Real reactivePowerTimeConstant,
                             Real dcVoltageTimeConstant,
                             Real acVoltageMagnitudeTimeConstant);

  /// Enable a second-order Padé approximation of the modulation delay for
  /// [mDelta_d, mDelta_q, mSigma_d, mSigma_q, mSigma_z].
  void setModulationDelay(Real timeDelay, UInt padeOrder = 2);

  /// Initial nominal EMT grid angle. The PLL deviation starts at zero.
  void setInitialAngle(Real angle);

  /// Preload the nonlinear MMC and controller states from a loaded AC
  /// operating point before the first SSN equivalent is constructed.
  ///
  /// This is required for non-zero-power initialization, particularly for a
  /// DC-voltage-controlled station: at steady state the voltage error is zero,
  /// while the active-controller integrator must still hold the non-zero
  /// d-axis current needed to transfer power.
  void setInitialOperatingPoint(Real activePower, Real reactivePower);

  /// Select the differential converter-voltage reference source. The external
  /// command is a 2x1 [d,q] vector in peak phase volts in the plant's nominal
  /// positive-sequence dq frame. Internally SSN_MMC converts it with
  /// mDelta_dq = -2*vMDelta_dq/(v_dc+ - v_dc-). InternalControllers is the
  /// construction default and preserves all existing examples.
  void setControlSource(ControlSource source);
  ControlSource controlSource() const { return mControlSource; }
  void setExternalDifferentialVoltageCommand(Real dVolts, Real qVolts);

  Attribute<Matrix>::Ptr acTerminalVoltageAttribute() const;
  Attribute<Matrix>::Ptr acTerminalCurrentAttribute() const;
  Attribute<Matrix>::Ptr interfaceVoltageAttribute() const;
  Attribute<Matrix>::Ptr interfaceCurrentAttribute() const;
  Attribute<Real>::Ptr dcPositiveVoltageAttribute() const;
  Attribute<Real>::Ptr dcNegativeVoltageAttribute() const;
  Attribute<Real>::Ptr dcVoltageAttribute() const;
  Attribute<Real>::Ptr dcCurrentAttribute() const;
  Attribute<Real>::Ptr activePowerAttribute() const;
  Attribute<Real>::Ptr reactivePowerAttribute() const;
  Attribute<Real>::Ptr storedEnergyAttribute() const;
  Attribute<Matrix>::Ptr externalDifferentialVoltageAttribute() const;
  Attribute<Bool>::Ptr externalCommandActiveAttribute() const;

  /// Limits applied to current references and modulation commands.
  void setLimits(Real maximumAcCurrent, Real maximumCirculatingCurrent,
                 Real maximumModulationMagnitude);

  void setNumericalLinearizationParameters(Real relativeStep,
                                           Real absoluteStep);

  /// Configure the nonlinear operating-point initialization. These are
  /// numerical solver settings, not physical MMC parameters.
  void setOperatingPointInitialization(Bool enabled, UInt maximumIterations,
                                       Real normalizedTolerance);

  /// Enable local continuous-time Jacobian diagnostics. Eigenvalues are
  /// evaluated only every updateInterval model rebuilds to limit overhead.
  void setEigenvalueDiagnostics(Bool enabled, UInt updateInterval = 100);

  /// Set the DPsim simulation time step used only for the discrete
  /// trapezoidal eigenvalue diagnostic. No default is assumed.
  void setDiagnosticTimeStep(Real timeStep);

  std::vector<String> getLocalStateNames() const override;

  Matrix getState() const;
  Matrix getStateDerivative() const;
  Matrix getInterfaceVoltage() const;
  Matrix getInterfaceCurrent() const;

protected:
  static constexpr UInt mInputSize = 5;
  static constexpr UInt mOutputSize = 5;
  static constexpr UInt mPlantStateSize = 12;
  static constexpr UInt mControllerStateSize = 8;
  // XiPll + PLL angle deviation + nominal EMT grid angle.
  static constexpr UInt mPllStateSize = 3;
  static constexpr UInt mFilterStateSize = 10;
  static constexpr UInt mDelayChannelCount = 5;
  static constexpr UInt mPadeOrder = 2;
  static constexpr UInt mDelayStateSize = mDelayChannelCount * mPadeOrder;
  static constexpr UInt mStateSize = mPlantStateSize + mControllerStateSize +
                                     mPllStateSize + mFilterStateSize +
                                     mDelayStateSize;

  enum InputIndex : UInt { Va = 0, Vb, Vc, Vdcp, Vdcn };
  enum OutputIndex : UInt { Ia = 0, Ib, Ic, Idcp, Idcn };

  enum StateIndex : UInt {
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

    XiActive,
    XiReactive,
    XiOccD,
    XiOccQ,
    XiCccD,
    XiCccQ,
    XiZcc,
    XiEnergy,

    XiPll,
    PllAngle,
    GridAngle,

    FilterVacD,
    FilterVacQ,
    FilterP1,
    FilterP2,
    FilterQ1,
    FilterQ2,
    FilterVdc1,
    FilterVdc2,
    FilterVacMag1,
    FilterVacMag2,

    DelayMDeltaD1,
    DelayMDeltaD2,
    DelayMDeltaQ1,
    DelayMDeltaQ2,
    DelayMSigmaD1,
    DelayMSigmaD2,
    DelayMSigmaQ1,
    DelayMSigmaQ2,
    DelayMSigmaZ1,
    DelayMSigmaZ2,
  };

  struct PIParameters {
    Real kp = 0.0;
    Real ki = 0.0;
  };

  MatrixComp buildInitialInputFromNodes(Real frequency);
  void validateTerminalArrangement() const;
  void validateInterfaceDimensions() const;

  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;

  void initializeFromNodesAndTerminals(Real frequency) override;
  Bool updateComponentParameters() override;
  void addHeldControlDependencies(
      AttributeBase::List &prevStepDependencies) const override;

  void evaluateStateDerivative(const Matrix &x, const Matrix &u,
                               Matrix &stateDerivative) const;
  void evaluateOutput(const Matrix &x, const Matrix &u, Matrix &output) const;

  void calculateNumericalJacobians(const Matrix &x, const Matrix &u, Matrix &A,
                                   Matrix &B, Matrix &C, Matrix &D) const;
  void buildStateSpaceModel(const Matrix &x, const Matrix &u, Matrix &A,
                            Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                            Matrix &F) const;

  Matrix abcToDq(const Matrix &abc, Real theta) const;
  Matrix dqToAbc(Real d, Real q, Real theta) const;
  Matrix rotateDq(const Matrix &dq, Real angle) const;

  Real clamp(Real value, Real lower, Real upper) const;
  Real regularizedDcVoltage(Real voltage) const;
  Real calculateStoredEnergy(const Matrix &x) const;
  Real applyFirstOrderFilter(const Matrix &x, Matrix &f, UInt stateIndex,
                             Real input, Real timeConstant, Bool enabled) const;
  Real applySecondOrderFilter(const Matrix &x, Matrix &f, UInt stateIndex1,
                              UInt stateIndex2, Real input, Real timeConstant,
                              Bool enabled) const;
  Real applyPadeDelayChannel(const Matrix &x, Matrix &f, UInt stateIndex1,
                             UInt stateIndex2, Real input) const;

  std::vector<UInt> getEquilibriumStateIndices() const;
  std::vector<UInt> getDiagnosticStateIndices() const;
  void initializeAnalyticalOperatingPoint(Matrix &x0, const Matrix &u0) const;
  Bool solveOperatingPoint(Matrix &x0, const Matrix &u0,
                           Real &normalizedResidual,
                           Real &absoluteResidual) const;
  void updateDiagnostics(const Matrix &A, const Matrix &x, const Matrix &u,
                         Bool force);
  void validateFinite(const Matrix &value, const char *name) const;

  void updateLogAttributes(const Matrix &u) const;

  // Physical parameters
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

  // Control configuration
  ActiveControlMode mActiveControlMode;
  ReactiveControlMode mReactiveControlMode;
  ControlSource mControlSource;

  PIParameters mActiveController;
  PIParameters mReactiveController;
  PIParameters mOutputCurrentController;
  PIParameters mCirculatingCurrentController;
  PIParameters mZeroSequenceCurrentController;
  PIParameters mEnergyController;

  Bool mEnergyControllerEnabled;
  Bool mPllEnabled;
  Bool mModulationDelayEnabled;

  Real mActivePowerReference;
  Real mReactivePowerReference;
  Real mDcVoltageReference;
  Real mAcVoltageReference;
  Real mDroopGain;

  Real mOpenLoopIDeltaDReference;
  Real mOpenLoopIDeltaQReference;
  Real mISigmaDReference;
  Real mISigmaQReference;
  Real mISigmaZReference;
  Real mInitialAngle;

  // Optional loaded operating point used for analytical/nonlinear
  // initialization before the first SSN model is constructed.
  Bool mInitialOperatingPointEnabled;
  Real mInitialActivePower;
  Real mInitialReactivePower;

  PIParameters mPllController;

  Real mAcVoltageDqFilterTimeConstant;
  Real mActivePowerFilterTimeConstant;
  Real mReactivePowerFilterTimeConstant;
  Real mDcVoltageFilterTimeConstant;
  Real mAcVoltageMagnitudeFilterTimeConstant;
  Real mModulationDelay;

  Real mMaximumAcCurrent;
  Real mMaximumCirculatingCurrent;
  Real mMaximumModulationMagnitude;

  Real mJacobianRelativeStep;
  Real mJacobianAbsoluteStep;

  Bool mOperatingPointInitializationEnabled;
  UInt mOperatingPointMaximumIterations;
  Real mOperatingPointNormalizedTolerance;

  Bool mEigenvalueDiagnosticsEnabled;
  UInt mEigenvalueDiagnosticsInterval;
  UInt mModelUpdateCounter;
  Real mDiagnosticTimeStep;

  static constexpr Real mMinimumDcVoltage = 1.0;

  // Logged attributes
  Attribute<Real>::Ptr mDcVoltage;
  Attribute<Real>::Ptr mDcPositiveVoltage;
  Attribute<Real>::Ptr mDcNegativeVoltage;
  Attribute<Real>::Ptr mDcCurrent;
  Attribute<Real>::Ptr mActivePower;
  Attribute<Real>::Ptr mReactivePower;
  Attribute<Real>::Ptr mAcVoltageMagnitude;
  Attribute<Real>::Ptr mStoredEnergy;
  Attribute<Real>::Ptr mConverterAngle;
  Attribute<Real>::Ptr mPllFrequency;
  Attribute<Real>::Ptr mFilteredActivePower;
  Attribute<Real>::Ptr mFilteredReactivePower;
  Attribute<Real>::Ptr mFilteredDcVoltage;

  Attribute<Real>::Ptr mStateNorm;
  Attribute<Real>::Ptr mStateDerivativeNorm;
  Attribute<Real>::Ptr mEquilibriumResidualNorm;
  Attribute<Real>::Ptr mJacobianMaximumRealEigenvalue;
  Attribute<Real>::Ptr mJacobianMaximumMagnitudeEigenvalue;
  Attribute<Real>::Ptr mJacobianDominantFrequency;
  Attribute<Real>::Ptr mJacobianMaximumDiscreteMagnitude;
  Attribute<Real>::Ptr mJacobianDiscreteDominantFrequency;
  Attribute<Real>::Ptr mGridAngle;
  Attribute<Real>::Ptr mPllAngleDeviation;
  Attribute<Real>::Ptr mPllError;
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
  Attribute<Real>::Ptr mNortonMatrixNorm;
  Attribute<Real>::Ptr mHistoryVectorNorm;
  Attribute<Real>::Ptr mDiagnosticsValid;
  Attribute<Matrix>::Ptr mAcTerminalVoltage;
  Attribute<Matrix>::Ptr mAcTerminalCurrent;
  Attribute<Matrix>::Ptr mExternalDifferentialVoltage;
  Attribute<Bool>::Ptr mExternalCommandActive;
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
