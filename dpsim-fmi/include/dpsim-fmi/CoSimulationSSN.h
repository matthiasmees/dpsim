// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-fmi/CoSimulation.h>
#include <dpsim-models/MNASimPowerComp.h>
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

/// Explicit fixed-step electrical coupling for an FMI Co-Simulation FMU.
/// The MNA solve uses the last completed FMU current. After each accepted
/// network solve, that terminal voltage is held over the next communication
/// step. No implicit Newton convergence is claimed for this component.
class CoSimulationSSN final : public CPS::MNASimPowerComp<CPS::Real>,
                              public SharedFactory<CoSimulationSSN> {
  CoSimulation::Ptr mModel;
  CoSimulationElectricalPort mPort;
  double mSign;
  double mTimeStep = 0;
  std::vector<std::size_t> mVoltageIndices, mCurrentIndices;
  std::vector<CPS::Attribute<CPS::Real>::Ptr> mControls;
  std::vector<std::function<double(double)>> mInputFunctions;
  std::vector<std::vector<double>> mInputEvents;
  CPS::Attribute<CPS::Matrix>::Ptr mInputs, mOutputs;
  CPS::Attribute<CPS::Real>::Ptr mStepAttribute, mInternalStepAttribute;

  void validateTerminals();
  void updateTerminalVoltage(const CPS::Matrix &solution);

public:
  using CPS::SimPowerComp<CPS::Real>::connect;
  CoSimulationSSN(CPS::String name, CoSimulation::Ptr model,
                  CoSimulationElectricalPort port);
  static std::shared_ptr<CoSimulationSSN>
  importFMU(CPS::String name, const std::string &path,
            CoSimulationElectricalPort port, CoSimulationOptions options = {});
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
  void mnaCompUpdateVoltage(const CPS::Matrix &) override;
  void mnaCompPreStep(CPS::Real, CPS::Int) override;
  void mnaCompPostStep(CPS::Real, CPS::Int,
                       CPS::Attribute<CPS::Matrix>::Ptr &) override;
  void mnaCompAddPreStepDependencies(CPS::AttributeBase::List &,
                                     CPS::AttributeBase::List &,
                                     CPS::AttributeBase::List &) override;
  void mnaCompAddPostStepDependencies(
      CPS::AttributeBase::List &, CPS::AttributeBase::List &,
      CPS::AttributeBase::List &, CPS::Attribute<CPS::Matrix>::Ptr &) override;
};

} // namespace DPsim::FMI
