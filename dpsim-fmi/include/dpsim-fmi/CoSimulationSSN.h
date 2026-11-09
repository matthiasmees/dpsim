// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-fmi/CoSimulation.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNANonlinearCompInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>
#include <functional>
#include <vector>

namespace DPsim::FMI {

struct CoSimulationElectricalPort {
  std::vector<std::string> voltageInputs;
  std::vector<std::string> currentOutputs;
  /// +1 when positive FMU current enters the passive terminal; -1 for a
  /// generator output current that is positive from the FMU into the bus.
  double currentSignToPassive = 1;
};

/// Electrical coupling policy. SinglePassExplicit remains the default and
/// advances the FMU once using the last accepted terminal voltage. RollbackNewton
/// replays trial communication steps from an FMI state checkpoint and iterates
/// the Norton coupling to the configured current tolerance.
enum class CoSimulationCouplingMode { SinglePassExplicit, RollbackNewton };

struct CoSimulationNewtonSettings {
  int maxIterations = 20;
  double absoluteCurrentTolerance = 1e-5;
  double relativeCurrentTolerance = 1e-7;
  double relativeJacobianPerturbation = 6e-6;
};

/// Explicit fixed-step electrical coupling for an FMI Co-Simulation FMU.
/// The MNA solve uses the last completed FMU current. After each accepted
/// network solve, that terminal voltage is held over the next communication
/// step. No implicit Newton convergence is claimed for this component.
class CoSimulationSSN final : public CPS::MNASimPowerComp<CPS::Real>,
                              public CPS::MNAVariableCompInterface,
                              public CPS::MNANonlinearCompInterface,
                              public SharedFactory<CoSimulationSSN> {
  CoSimulation::Ptr mModel;
  CoSimulationElectricalPort mPort;
  CoSimulationCouplingMode mCouplingMode;
  CoSimulationNewtonSettings mNewtonSettings;
  double mSign;
  double mTimeStep = 0;
  std::vector<std::size_t> mVoltageIndices, mCurrentIndices;
  std::vector<CPS::Attribute<CPS::Real>::Ptr> mControls;
  std::vector<std::function<double(double)>> mInputFunctions;
  std::vector<std::vector<double>> mInputEvents;
  CPS::Attribute<CPS::Matrix>::Ptr mInputs, mOutputs;
  CPS::Attribute<CPS::Real>::Ptr mStepAttribute, mInternalStepAttribute;
  CPS::Attribute<CPS::Real>::Ptr mIterationAttribute, mResidualAttribute;
  CPS::Matrix mStepStartVoltage, mTrialVoltage, mTrialInputs, mTrialOutputs,
      mTrialCurrent;
  CPS::Matrix mNortonConductance, mNortonHistory;
  double mCurrentTime = 0, mCurrentResidual = 0;
  double mControlTime = 0;
  CPS::Int mIterations = 0;
  bool mRollbackStep = false;

  void validateTerminals();
  void updateTerminalVoltage(const CPS::Matrix &solution);
  CPS::Matrix voltageFromSolution(const CPS::Matrix &solution);
  CPS::Matrix currents(const CPS::Matrix &outputs) const;
  CPS::Matrix inputsAt(const CPS::Matrix &voltage, double controlTime) const;
  CPS::Matrix replayAtVoltage(const CPS::Matrix &voltage);
  void prepareRollbackLinearization();

public:
  using CPS::SimPowerComp<CPS::Real>::connect;
  CoSimulationSSN(CPS::String name, CoSimulation::Ptr model,
                  CoSimulationElectricalPort port,
                  CoSimulationCouplingMode couplingMode =
                      CoSimulationCouplingMode::SinglePassExplicit,
                  CoSimulationNewtonSettings newtonSettings = {});
  static std::shared_ptr<CoSimulationSSN>
  importFMU(CPS::String name, const std::string &path,
            CoSimulationElectricalPort port, CoSimulationOptions options = {},
            CoSimulationCouplingMode couplingMode =
                CoSimulationCouplingMode::SinglePassExplicit,
            CoSimulationNewtonSettings newtonSettings = {});
  CPS::SimPowerComp<CPS::Real>::Ptr clone(CPS::String name) override;
  void connect(CPS::SimNode<CPS::Real>::Ptr node);
  void initialize(CPS::Matrix frequencies) override;
  void initializeFromNodesAndTerminals(CPS::Real frequency) override;
  void setInputFunction(const std::string &name,
                        std::function<double(double)> function,
                        std::vector<double> discontinuityTimes = {});
  void mnaCompInitialize(CPS::Real, CPS::Real,
                         CPS::Attribute<CPS::Matrix>::Ptr) override;
  void mnaCompApplyRightSideVectorStamp(CPS::Matrix &) override;
  void mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &) override;
  void mnaCompUpdateVoltage(const CPS::Matrix &) override;
  void mnaCompPreStep(CPS::Real, CPS::Int) override;
  void mnaCompPostStep(CPS::Real, CPS::Int,
                       CPS::Attribute<CPS::Matrix>::Ptr &) override;
  bool mnaCheckNonlinearConvergence(const CPS::Matrix &) override;
  void mnaPrepareNonlinearIteration() override;
  CPS::Bool hasParameterChanged() override {
    return mCouplingMode == CoSimulationCouplingMode::RollbackNewton;
  }
  void mnaCompAddPreStepDependencies(CPS::AttributeBase::List &,
                                     CPS::AttributeBase::List &,
                                     CPS::AttributeBase::List &) override;
  void mnaCompAddPostStepDependencies(
      CPS::AttributeBase::List &, CPS::AttributeBase::List &,
      CPS::AttributeBase::List &, CPS::Attribute<CPS::Matrix>::Ptr &) override;
};

} // namespace DPsim::FMI
