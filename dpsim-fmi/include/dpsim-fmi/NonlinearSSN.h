// SPDX-License-Identifier: MPL-2.0
#pragma once
#include <array>
#include <dpsim-fmi/ModelExchange.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNANonlinearCompInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>
#include <functional>

namespace DPsim::FMI {
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
  double relativeJacobianPerturbation = 6e-6;
  int maxIterations = 20;
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
  CPS::Attribute<CPS::Matrix>::Ptr mState, mInputs, mOutputs;
  CPS::Attribute<CPS::Real>::Ptr mIterationAttribute, mResidualAttribute;
  std::vector<CPS::Attribute<CPS::Real>::Ptr> mControls;
  std::vector<std::function<double(double)>> mInputFunctions;
  std::vector<std::vector<double>> mInputEvents;
  double mStateResidual = 0, mCurrentResidual = 0;
  int mIterations = 0;
  void validateTerminals();
  CPS::Matrix inputsAt(const CPS::Matrix &voltage, double time,
                       bool left) const;
  CPS::Matrix currents(const CPS::Matrix &outputs) const;
  CPS::Matrix voltageFromSolution(const CPS::Matrix &solution);
  bool eventAt(double time) const;
  void initializeValues();

public:
  NonlinearSSN(CPS::String name, ModelExchange::Ptr model,
               const ElectricalPort3Ph &port, NewtonSettings settings = {});
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
} // namespace DPsim::FMI
