// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-models/EMT/EMT_LinearSSNModel.h>
#include <dpsim-models/MNASimPowerComp.h>

#include <unordered_map>

namespace CPS::EMT {

/// Paper-style SSN group generated from a collection of linear companion
/// models. Internal network nodes are eliminated once by a Schur complement;
/// only the retained terminal nodes remain in the surrounding MNA problem.
class LinearSSNGroup final : public MNASimPowerComp<Real> {
public:
  using Ptr = std::shared_ptr<LinearSSNGroup>;
  using Component = MNASimPowerComp<Real>;
  using ComponentList = std::vector<Component::Ptr>;

  static Ptr make(String name, ComponentList components, SimNode::List allNodes,
                  SimNode::List retainedNodes,
                  Logger::Level level = Logger::Level::off) {
    return std::make_shared<LinearSSNGroup>(
        std::move(name), std::move(components), std::move(allNodes),
        std::move(retainedNodes), level);
  }

  LinearSSNGroup(String name, ComponentList components, SimNode::List allNodes,
                 SimNode::List retainedNodes,
                 Logger::Level level = Logger::Level::off);

  void initializeFromNodesAndTerminals(Real frequency) override;
  Bool mnaHasCompactRightSideVectorStamp() const override { return true; }
  void mnaAddCompactRightSideVectorStamp(Matrix &rightVector) override;

protected:
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;
  void mnaCompUpdateVoltage(const Matrix &) override {}
  void mnaCompUpdateCurrent(const Matrix &) override {}
  void mnaCompPreStep(Real, Int) override;
  void mnaCompPostStep(Real, Int, Attribute<Matrix>::Ptr &leftVector) override;
  void mnaCompAddPreStepDependencies(
      AttributeBase::List &prevStepDependencies,
      AttributeBase::List &attributeDependencies,
      AttributeBase::List &modifiedAttributes) override;
  void
  mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                 AttributeBase::List &attributeDependencies,
                                 AttributeBase::List &modifiedAttributes,
                                 Attribute<Matrix>::Ptr &leftVector) override;

private:
  ComponentList mComponents;
  SimNode::List mAllNodes;
  std::vector<String> mRetainedNames;
  std::vector<LinearSSNModel> mModels;
  Attribute<Matrix>::Ptr mState;
  SparseMatrix mF, mBi, mBn, mHi, mHn, mGin, mGni;
  Matrix mDd;
  Eigen::PartialPivLU<Matrix> mInternalFactor;
  std::vector<UInt> mStateOffsets;
  std::vector<std::vector<UInt>> mModelInputIndices;
  std::unordered_map<String, UInt> mLocalNodeIndex;
  std::unordered_map<String, UInt> mRetainedNodeIndex;
  std::vector<UInt> mInternalIndices;
  std::vector<UInt> mTerminalIndices;
  // Reused per-step workspaces.  The paper realization is assembled once,
  // but its state/history vectors are updated at every EMT step.
  mutable Matrix mFullHistoryScratch;
  mutable Matrix mInternalHistoryScratch;
  mutable Matrix mInternalVoltageScratch;
  Matrix mLocalVoltageScratch;
  Matrix mTerminalInputScratch;
  mutable Matrix mTerminalHistoryScratch;
  Matrix mStateUpdateScratch;
  Bool mReady = false;

  void buildRealization(Real omega, Real timeStep,
                        Attribute<Matrix>::Ptr leftVector);
  void assembleFullHistory(Matrix &fullHistory) const;
  void currentHistory(Matrix &history) const;
  Matrix terminalVoltage(const Matrix &leftVector);
};

} // namespace CPS::EMT
