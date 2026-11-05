// SPDX-License-Identifier: MPL-2.0
#pragma once
#include <dpsim-models/MNASimPowerComp.h>
#include <functional>

/// Ideal prescribed voltage boundary, independent of the imported branch.
/// Adds a branch-current unknown and a voltage constraint to DPsim's MNA.
class VoltageBoundary final : public CPS::MNASimPowerComp<CPS::Real>,
                              public SharedFactory<VoltageBoundary> {
  std::function<CPS::Real(CPS::Real)> mFunction;

public:
  VoltageBoundary(CPS::String name,
                  std::function<CPS::Real(CPS::Real)> function)
      : MNASimPowerComp(name, true, true), mFunction(std::move(function)) {
    setTerminalNumber(2);
    setVirtualNodeNumber(1);
    **mIntfVoltage = CPS::Matrix::Constant(1, 1, mFunction(0));
    **mIntfCurrent = CPS::Matrix::Zero(1, 1);
  }
  void initializeFromNodesAndTerminals(CPS::Real) override {}
  void mnaCompInitialize(CPS::Real, CPS::Real,
                         CPS::Attribute<CPS::Matrix>::Ptr) override {
    updateMatrixNodeIndices();
  }
  void mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &matrix) override {
    const auto branch = mVirtualNodes[0]->matrixNodeIndex();
    for (CPS::UInt index = 0; index < 2; ++index) {
      if (terminalNotGrounded(index)) {
        const CPS::Real sign = index ? 1 : -1;
        CPS::Math::addToMatrixElement(matrix, matrixNodeIndex(index), branch,
                                      sign);
        CPS::Math::addToMatrixElement(matrix, branch, matrixNodeIndex(index),
                                      sign);
      }
    }
  }
  void mnaCompPreStep(CPS::Real time, CPS::Int) override {
    (**mIntfVoltage)(0, 0) = mFunction(time);
    CPS::Math::setVectorElement(**mRightVector,
                                mVirtualNodes[0]->matrixNodeIndex(),
                                (**mIntfVoltage)(0, 0));
  }
  void mnaCompPostStep(CPS::Real, CPS::Int,
                       CPS::Attribute<CPS::Matrix>::Ptr &solution) override {
    (**mIntfCurrent)(0, 0) = CPS::Math::realFromVectorElement(
        **solution, mVirtualNodes[0]->matrixNodeIndex());
  }
  void
  mnaCompAddPreStepDependencies(CPS::AttributeBase::List &,
                                CPS::AttributeBase::List &,
                                CPS::AttributeBase::List &modified) override {
    modified.push_back(mRightVector);
    modified.push_back(mIntfVoltage);
  }
  void mnaCompAddPostStepDependencies(
      CPS::AttributeBase::List &, CPS::AttributeBase::List &dependencies,
      CPS::AttributeBase::List &modified,
      CPS::Attribute<CPS::Matrix>::Ptr &solution) override {
    dependencies.push_back(solution);
    modified.push_back(mIntfCurrent);
  }
};
