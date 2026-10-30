// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-models/MNASimPowerComp.h>

namespace CPS::EMT {

/// Composite MNA component used to benchmark a single scheduler task for the
/// fixed linear EMT network.  The contained SSN elements keep their original
/// state-space models and stamps; the aggregate only collapses their MNA
/// lifecycle and scheduler representation into one component.
class LinearSSNAggregate final : public MNASimPowerComp<Real> {
public:
  using Ptr = std::shared_ptr<LinearSSNAggregate>;
  using Component = MNASimPowerComp<Real>;
  using ComponentList = std::vector<Component::Ptr>;

  static Ptr make(String name, ComponentList components, SimNode::List nodes,
                  Logger::Level level = Logger::Level::off) {
    return std::make_shared<LinearSSNAggregate>(
        std::move(name), std::move(components), std::move(nodes), level);
  }

  LinearSSNAggregate(String name, ComponentList components, SimNode::List nodes,
                     Logger::Level level = Logger::Level::off);

  void initializeFromNodesAndTerminals(Real frequency) override;

  Bool mnaHasCompactRightSideVectorStamp() const override { return true; }
  void mnaAddCompactRightSideVectorStamp(Matrix &rightVector) override;

protected:
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  void mnaCompUpdateCurrent(const Matrix &leftVector) override;
  void mnaCompPreStep(Real time, Int timeStepCount) override;
  void mnaCompPostStep(Real time, Int timeStepCount,
                       Attribute<Matrix>::Ptr &leftVector) override;
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
  std::vector<Component *> mCompactComponents;
  std::vector<Component *> mDenseComponents;
};

} // namespace CPS::EMT
