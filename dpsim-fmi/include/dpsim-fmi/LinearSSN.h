// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-fmi/LinearModel.h>
#include <dpsim-models/EMT/EMT_VTypeSSNComp.h>

namespace DPsim::FMI {

/// One-channel or three-phase EMT electrical branch. Terminal order is {negative, positive};
/// v = v1-v0 and current flows from terminal 1 to terminal 0.
/// DPsim's general SSN base performs trapezoidal integration inside MNA.
class LinearSSN final : public CPS::EMT::VTypeSSNComp,
                        public SharedFactory<LinearSSN> {
  LinearModel::Ptr mModel;
  void validateTerminals() const;
  CPS::MatrixComp buildInitialInputFromNodes(CPS::Real frequency) override;

public:
  using CPS::SimPowerComp<CPS::Real>::connect;

  LinearSSN(CPS::String name, LinearModel::Ptr model,
            CPS::Logger::Level logLevel = CPS::Logger::Level::off);
  /// Ground the negative terminal and connect the positive terminal to a node.
  void connect(CPS::SimNode<CPS::Real>::Ptr node);
  CPS::SimPowerComp<CPS::Real>::Ptr clone(CPS::String name) override;
  const LinearModel &model() const { return *mModel; }
  std::vector<CPS::String> getLocalStateNames() const override;
  void initialize(CPS::Matrix frequencies) override;
  /// Preserve FMU initial states instead of replacing them with AC steady state.
  void initializeFromNodesAndTerminals(CPS::Real frequency) override;
  void mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &matrix) override;
  void mnaCompApplyRightSideVectorStamp(CPS::Matrix &vector) override;
  void mnaCompUpdateVoltage(const CPS::Matrix &solution) override;
};
} // namespace DPsim::FMI
