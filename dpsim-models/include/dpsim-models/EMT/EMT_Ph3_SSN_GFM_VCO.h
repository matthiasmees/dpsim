// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <vector>

#include <dpsim-models/EMT/EMT_Ph3_TwoTerminalVTypeVariableSSNComp.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

/// \brief Fixed-frequency grid-forming VSI with cascaded voltage/current PI
/// controllers, represented as a locally linearized affine SSN model.
///
/// This class is the SSN counterpart of VSIVoltageControlVCO. It retains:
///
///   - the fixed-frequency VCO angle,
///   - the outer dq voltage PI controller,
///   - the inner dq current PI controller,
///   - the three-phase L-C-R filter,
///   - the coupling resistance to the PCC.
///
/// It deliberately does not include a PLL, active/reactive-power controller,
/// virtual-synchronous-machine dynamics, excitation dynamics, active damping,
/// or a converter-delay state.
///
/// The nonlinear model
///
///   x_dot = f(x, u)
///   y     = g(x, u)
///
/// is linearized at the current operating point:
///
///   x_dot = A x + B u + E
///   y     = C x + D u + F
///
/// The inherited terminal convention is:
///
///   u = v_terminal1 - v_terminal0
///
/// and the SSN current entering the component is:
///
///   y = (u - v_c) / R_c.
///
/// The physical inverter current injected into the grid is therefore -y.
class SSN_GFM_VCO final : public TwoTerminalVTypeVariableSSNComp,
                          public SharedFactory<SSN_GFM_VCO> {

private:
  static constexpr Int mStateSize = 11;
  static constexpr Int mInputSize = 3;
  static constexpr Int mOutputSize = 3;

  enum StateIndex : Int {
    Theta = 0,

    VoltageIntegratorD = 1,
    VoltageIntegratorQ = 2,

    CurrentIntegratorD = 3,
    CurrentIntegratorQ = 4,

    VcA = 5,
    VcB = 6,
    VcC = 7,

    IfA = 8,
    IfB = 9,
    IfC = 10
  };

  // Electrical filter.
  Real mLf;
  Real mCf;
  Real mRf;
  Real mRc;

  // Fixed-frequency VCO and voltage references.
  Real mSystemOmega;
  Real mOmegaN;
  Real mVdRef;
  Real mVqRef;

  // Cascaded voltage/current PI controllers.
  Real mKpVoltage;
  Real mKiVoltage;
  Real mKpCurrent;
  Real mKiCurrent;

  // Track the legacy staged configuration API.
  Bool mGeneralParametersSet;
  Bool mControllerParametersSet;
  Bool mFilterParametersSet;

  // Optional explicit controller-state initialization.
  Bool mInitialStateOverrideEnabled;
  Real mInitialVoltageIntegratorD;
  Real mInitialVoltageIntegratorQ;
  Real mInitialCurrentIntegratorD;
  Real mInitialCurrentIntegratorQ;

  // Optional explicit initial three-phase injected power. If disabled, the
  // power transferred by initWithPowerflow() to terminal 0 is used.
  Bool mInitialPowerOverrideEnabled;
  Real mInitialActivePower;
  Real mInitialReactivePower;

  // Local-linearization scheduling. Interval one is the accurate default.
  UInt mLinearizationUpdateInterval;
  UInt mStepsSinceLinearization;
  Bool mLinearizationInitialized;

  // Logging attributes. Names intentionally match SSN_GFM where possible so
  // diagnostics can exchange the two models with minimal changes.
  const Attribute<Real>::Ptr mPInst;
  const Attribute<Real>::Ptr mQInst;
  const Attribute<Real>::Ptr mOmegaGFM;
  const Attribute<Real>::Ptr mThetaGFM;
  const Attribute<Real>::Ptr mVoltageMagnitudeGFM;

  const Attribute<Real>::Ptr mVcD;
  const Attribute<Real>::Ptr mVcQ;

  const Attribute<Real>::Ptr mIGridD;
  const Attribute<Real>::Ptr mIGridQ;

  const Attribute<Real>::Ptr mIfD;
  const Attribute<Real>::Ptr mIfQ;

  const Attribute<Real>::Ptr mVoltageReferenceD;
  const Attribute<Real>::Ptr mVoltageReferenceQ;

  const Attribute<Real>::Ptr mConverterVoltageD;
  const Attribute<Real>::Ptr mConverterVoltageQ;

  Matrix getParkTransformMatrix(Real theta) const;
  Matrix getInverseParkTransformMatrix(Real theta) const;

  void configureStateSpaceDimensionsIfReady();
  void markLinearizationDirty();

  /// \brief Evaluate the nonlinear state derivative x_dot = f(x,u).
  void evaluateStateDerivative(const Matrix &x, const Matrix &u,
                               Matrix &stateDerivative) const;

  /// \brief Evaluate the nonlinear SSN output y = g(x,u).
  void evaluateOutput(const Matrix &x, const Matrix &u, Matrix &output) const;

  /// \brief Calculate the exact local Jacobians A, B, C and D.
  void calculateAnalyticalJacobians(const Matrix &x, const Matrix &u, Matrix &A,
                                    Matrix &B, Matrix &C, Matrix &D) const;

  /// \brief Construct the complete local affine state-space model.
  void buildStateSpaceModel(const Matrix &x, const Matrix &u, Matrix &A,
                            Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                            Matrix &F) const;

protected:
  Bool updateComponentParameters() override final;
  void updateLogAttributes(const Matrix &u) const override final;

public:
  using SharedFactory<SSN_GFM_VCO>::make;

  SSN_GFM_VCO(String uid, String name,
              Logger::Level logLevel = Logger::Level::off);

  SSN_GFM_VCO(String name, Logger::Level logLevel = Logger::Level::off)
      : SSN_GFM_VCO(name, name, logLevel) {}

  std::vector<String> getLocalStateNames() const override final;

  std::vector<SSNComp::LocalAbcStateBlock>
  getLocalAbcStateBlocks() const override final;

  /// \brief Legacy-compatible general VCO/reference setter.
  ///
  /// With the power-invariant Park transform used by the reference class,
  /// VdRef is normally the line-to-line RMS voltage of a balanced system.
  void setParameters(Real sysOmega, Real vdRef, Real vqRef);

  /// \brief Legacy-compatible cascaded-controller setter.
  void setControllerParameters(Real kpVoltage, Real kiVoltage, Real kpCurrent,
                               Real kiCurrent, Real omegaNominal);

  /// \brief Legacy-compatible L-C-R filter setter.
  void setFilterParameters(Real lf, Real cf, Real rf, Real rc);

  /// \brief Combined SSN-style configuration setter.
  void setParameters(Real lf, Real cf, Real rf, Real rc, Real omegaN,
                     Real vdRef, Real vqRef, Real kpVoltage, Real kiVoltage,
                     Real kpCurrent, Real kiCurrent);

  /// \brief Override the automatically calculated PI-integrator initial states.
  void setInitialStateValues(Real voltageIntegratorD, Real voltageIntegratorQ,
                             Real currentIntegratorD, Real currentIntegratorQ);

  /// \brief Return to automatic equilibrium initialization of PI states.
  void clearInitialStateValues();

  /// \brief Explicitly set the initial injected three-phase P/Q.
  ///
  /// This is useful when the power-flow image is not a synchronous generator,
  /// because the current DPsim initWithPowerflow() path only transfers
  /// synchronous-generator terminal power automatically.
  void setInitialPower(Real activePower, Real reactivePower);

  /// \brief Return to terminal-power-based initialization.
  void clearInitialPower();

  /// \brief Compatibility no-op matching the analytical SSN_GFM API.
  void setNumericalLinearizationParameters(Real relativeStep,
                                           Real absoluteStep);

  /// \brief Rebuild the local model every updateInterval EMT steps.
  ///
  /// The default value one performs exact analytical relinearization every
  /// step. Larger values are an explicit approximation and can change results.
  void setLinearizationUpdateInterval(UInt updateInterval);

  void initializeFromNodesAndTerminals(Real frequency) override final;

  Matrix getState() const;
  Matrix getStateDerivative() const;
  Matrix getInterfaceVoltage() const;
  Matrix getInterfaceCurrent() const;
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
