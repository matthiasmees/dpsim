// SPDX-License-Identifier: MPL-2.0
#pragma once
#include <array>
#include <cstdint>
#include <dpsim-fmi/ModelExchange.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNANonlinearCompInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>
#include <functional>

namespace DPsim::FMI {
enum class MultirateMode {
  /// Restore and replay the complete interval for converged network Newton.
  RollbackNewton,
  /// Advance the FMU once after the network solve. Terminal current is lagged
  /// by one macro step and the coupling discrepancy is reported, not iterated.
  SinglePassExplicit
};
struct ElectricalPort3Ph {
  std::array<std::string, 3> voltageInputs;
  std::array<std::string, 3> currentOutputs;
  /// +1 for current entering positive terminal, -1 for generation into bus.
  double currentSignToPassive = 1;
};
struct NewtonSettings {
  double absoluteStateTolerance = 1e-11;
  double relativeStateTolerance = 1e-10;
  double absoluteCurrentTolerance = 1e-9;
  double relativeCurrentTolerance = 1e-10;
  /// Replay and finite-difference noise make the interval-map current less
  /// exact than the single-rate local residual.
  double multirateAbsoluteCurrentTolerance = 1e-6;
  double multirateRelativeCurrentTolerance = 1e-6;
  /// Permit a single-pass run to accept a discontinuous terminal current at an
  /// FMU event. This is useful for explicit co-simulation, but does not imply
  /// implicit network convergence at the event.
  bool allowInstantaneousCurrentEvents = false;
  double relativeJacobianPerturbation = 6e-6;
  int maxIterations = 20;
  /// Optional chord predictor across accepted steps. Any failed true residual
  /// check refreshes the Jacobian before the next iteration. Zero disables it.
  int maxJacobianReuseSteps = 0;
  /// Automatically replay an FMI interval when internal time events are
  /// faster than the network step. The FMU must support state rollback.
  bool automaticEventSubstepping = true;
  MultirateMode multirateMode = MultirateMode::RollbackNewton;
  /// Guard against an invalid or impractically dense event schedule.
  int maxInternalSteps = 10000;
};

/// Implicit trapezoidal FMI state residual condensed into an MNA Norton
/// tangent. The solver repeats Newton steps until true state and KCL current
/// residuals pass; attributes/history are committed only in post-step.
class NonlinearSSN final : public CPS::MNASimPowerComp<CPS::Real>,
                           public CPS::MNAVariableCompInterface,
                           public CPS::MNANonlinearCompInterface,
                           public SharedFactory<NonlinearSSN> {
  ModelExchange::Ptr mModel;
  NewtonSettings mSettings;
  std::array<std::size_t, 3> mVoltageIndices, mCurrentRows;
  double mSign, mStep = 0, mTime = 0;
  CPS::Matrix mInitialVoltage, mOldState, mOldDerivative;
  CPS::Matrix mTrialState, mTrialVoltage, mTrialDerivative, mTrialOutputs;
  CPS::Matrix mLinearizationState, mLinearizationVoltage, mCorrection,
      mSensitivity;
  CPS::Matrix mNewInputs, mW, mHistory;
  CPS::Matrix mMacroOldVoltage;
  CPS::Attribute<CPS::Matrix>::Ptr mState, mInputs, mOutputs;
  CPS::Attribute<CPS::Real>::Ptr mIterationAttribute, mResidualAttribute,
      mInternalStepAttribute, mCouplingErrorAttribute,
      mEvaluationRequestAttribute, mEvaluationTimeAttribute;
  CPS::Attribute<CPS::Matrix>::Ptr mNortonConductanceAttribute,
      mNortonHistoryAttribute, mNortonStateAttribute,
      mNortonStateCorrectionAttribute, mNortonVoltageAttribute;
  std::vector<CPS::Attribute<CPS::Real>::Ptr> mControls;
  std::vector<std::function<double(double)>> mInputFunctions;
  std::vector<std::vector<double>> mInputEvents;
  double mStateResidual = 0, mCurrentResidual = 0;
  double mCouplingError = 0;
  int mIterations = 0;
  int mInternalSteps = 0;
  std::uint64_t mEvaluationRequestBaseline = 0;
  double mEvaluationTimeBaseline = 0;
  bool mMultirateStep = false;
  bool mSinglePassStep = false;
  bool mRollbackVerified = false;
  CPS::Matrix mCachedA, mCachedB, mCachedC, mCachedD;
  int mJacobianAge = 0;
  bool mHasJacobian = false;
  CPS::Matrix mInternalCachedA;
  int mInternalJacobianAge = 0;
  bool mHasInternalJacobian = false;
  void validateTerminals();
  CPS::Matrix inputsAt(const CPS::Matrix &voltage, double time,
                       bool left) const;
  CPS::Matrix currents(const CPS::Matrix &outputs) const;
  CPS::Matrix voltageFromSolution(const CPS::Matrix &solution);
  bool eventAt(double time) const;
  void initializeValues();
  struct IntervalResult {
    CPS::Matrix state, inputs, outputs, current;
    double residual = 0;
    double nextEvent = 0;
    int steps = 0;
  };
  IntervalResult integrateInterval(const CPS::Matrix &endVoltage,
                                   bool restoreCheckpoint = true);
  void prepareMultirateIteration();

public:
  using CPS::SimPowerComp<CPS::Real>::connect;

  /// One-call import for the conventional three-phase electrical interface.
  /// Inputs and start values are read from modelDescription.xml.
  static std::shared_ptr<NonlinearSSN>
  importFMU(CPS::String name, const std::string &path,
            ModelExchangeOptions modelOptions = {},
            NewtonSettings settings = {});
  /// Import using the model name from modelDescription.xml as the component
  /// name. This is the shortest form for the conventional abc interface.
  static std::shared_ptr<NonlinearSSN>
  importFMU(const std::string &path, ModelExchangeOptions modelOptions = {},
            NewtonSettings settings = {});
  const ModelDescription &modelDescription() const {
    return mModel->description();
  }
  /// Conventional abc names (v_a/v_b/v_c or u_a/u_b/u_c with i_a/i_b/i_c)
  /// are inferred for the common case. Use the explicit overload for
  /// nonstandard port names or current polarity.
  NonlinearSSN(CPS::String name, ModelExchange::Ptr model,
               NewtonSettings settings = {});
  NonlinearSSN(CPS::String name, ModelExchange::Ptr model,
               const ElectricalPort3Ph &port, NewtonSettings settings = {});
  /// Connect the conventional three-phase terminal to a bus, with the
  /// negative terminal grounded. Use the inherited two-node connect overload
  /// for line-to-line or otherwise non-grounded connections.
  void connect(CPS::SimNode<CPS::Real>::Ptr node);
  CPS::SimPowerComp<CPS::Real>::Ptr clone(CPS::String) override;
  void initialize(CPS::Matrix frequencies) override;
  void initializeFromNodesAndTerminals(CPS::Real) override;
  /// Continuous control function with optional known discontinuities. Solver
  /// steps must land on each event; left/right endpoint handling is explicit.
  void setInputFunction(const std::string &name,
                        std::function<double(double)> function,
                        std::vector<double> discontinuityTimes = {});
  void mnaCompInitialize(CPS::Real, CPS::Real,
                         CPS::Attribute<CPS::Matrix>::Ptr) override;
  void mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &) override;
  void mnaCompApplyRightSideVectorStamp(CPS::Matrix &) override;
  void mnaCompPreStep(CPS::Real, CPS::Int) override;
  void mnaCompPostStep(CPS::Real, CPS::Int,
                       CPS::Attribute<CPS::Matrix>::Ptr &) override;
  void mnaCompAddPreStepDependencies(CPS::AttributeBase::List &,
                                     CPS::AttributeBase::List &,
                                     CPS::AttributeBase::List &) override;
  void mnaCompAddPostStepDependencies(
      CPS::AttributeBase::List &, CPS::AttributeBase::List &,
      CPS::AttributeBase::List &, CPS::Attribute<CPS::Matrix>::Ptr &) override;
  CPS::Bool hasParameterChanged() override { return true; }
  bool mnaCheckNonlinearConvergence(const CPS::Matrix &) override;
  void mnaPrepareNonlinearIteration() override;
};

/// Import a conventional three-phase Model Exchange FMU with its model name.
/// This is the concise public entry point: `FMI::importFMU(path)`.
std::shared_ptr<NonlinearSSN> importFMU(const std::string &path);
std::shared_ptr<NonlinearSSN> importFMU(const std::string &path,
                                        ModelExchangeOptions modelOptions,
                                        NewtonSettings settings = {});
std::shared_ptr<NonlinearSSN> importFMU(CPS::String name,
                                        const std::string &path,
                                        ModelExchangeOptions modelOptions = {},
                                        NewtonSettings settings = {});
} // namespace DPsim::FMI
