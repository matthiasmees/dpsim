// SPDX-License-Identifier: MPL-2.0
#include <dpsim-models/EMT/EMT_LinearSSNAggregate.h>

using namespace CPS;

EMT::LinearSSNAggregate::LinearSSNAggregate(String name,
                                            ComponentList components,
                                            SimNode::List nodes,
                                            Logger::Level level)
    : MNASimPowerComp<Real>(name, true, true, level),
      mComponents(std::move(components)) {
  mPhaseType = PhaseType::ABC;
  setTerminalNumber(static_cast<UInt>(nodes.size()));
  connect(std::move(nodes));
  mCompactComponents.reserve(mComponents.size());
  mDenseComponents.reserve(mComponents.size());
  for (const auto &component : mComponents) {
    if (component->mnaHasCompactRightSideVectorStamp())
      mCompactComponents.push_back(component.get());
    else
      mDenseComponents.push_back(component.get());
  }
}

void EMT::LinearSSNAggregate::initializeFromNodesAndTerminals(Real frequency) {
  for (const auto &component : mComponents)
    component->initializeFromNodesAndTerminals(frequency);
}

void EMT::LinearSSNAggregate::mnaCompInitialize(
    Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
  for (const auto &component : mComponents)
    component->mnaInitialize(omega, timeStep, leftVector);
}

void EMT::LinearSSNAggregate::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  for (const auto &component : mComponents)
    component->mnaCompApplySystemMatrixStamp(systemMatrix);
}

void EMT::LinearSSNAggregate::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  // Child components use set-style stamps for their initial source vector.
  // They normally stamp into separate full-size RHS attributes, which the
  // solver later sums.  Passing the aggregate RHS directly would make the
  // last child touching a row overwrite earlier contributions on that row.
  // Keep the child buffers isolated, then add each complete child RHS here.
  for (const auto &component : mComponents) {
    auto childRightVector = component->getRightVector();
    childRightVector->get().setZero();
    if (component->mnaHasCompactRightSideVectorStamp())
      component->mnaAddCompactRightSideVectorStamp(childRightVector->get());
    else
      component->mnaApplyRightSideVectorStamp(childRightVector->get());
    rightVector += childRightVector->get();
  }
}

void EMT::LinearSSNAggregate::mnaAddCompactRightSideVectorStamp(
    Matrix &rightVector) {
  // Use direct compact stamps whenever a child provides them. Keep a dense
  // fallback for child types that still use set-style full-vector stamps.
  for (auto *component : mCompactComponents)
    component->mnaAddCompactRightSideVectorStamp(rightVector);
  for (auto *component : mDenseComponents) {
    auto childRightVector = component->getRightVector();
    childRightVector->get().setZero();
    component->mnaApplyRightSideVectorStamp(childRightVector->get());
    rightVector += childRightVector->get();
  }
}

void EMT::LinearSSNAggregate::mnaCompUpdateVoltage(const Matrix &leftVector) {
  for (const auto &component : mComponents)
    component->mnaCompUpdateVoltage(leftVector);
}

void EMT::LinearSSNAggregate::mnaCompUpdateCurrent(const Matrix &leftVector) {
  for (const auto &component : mComponents)
    component->mnaCompUpdateCurrent(leftVector);
}

void EMT::LinearSSNAggregate::mnaCompPreStep(Real time, Int timeStepCount) {
  for (const auto &component : mComponents)
    component->mnaCompPreStep(time, timeStepCount);
}

void EMT::LinearSSNAggregate::mnaCompPostStep(
    Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
  for (const auto &component : mComponents)
    component->mnaCompPostStep(time, timeStepCount, leftVector);
}

void EMT::LinearSSNAggregate::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes) {
  // The solver schedules its solve task from the aggregate RHS attribute.
  // Publish it as modified even though the compact stamp reads the child RHS
  // buffers directly; this preserves the pre-step -> solve ordering.
  modifiedAttributes.push_back(mRightVector);
  for (const auto &component : mComponents)
    component->mnaAddPreStepDependencies(
        prevStepDependencies, attributeDependencies, modifiedAttributes);
}

void EMT::LinearSSNAggregate::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  for (const auto &component : mComponents)
    component->mnaAddPostStepDependencies(prevStepDependencies,
                                          attributeDependencies,
                                          modifiedAttributes, leftVector);
}
