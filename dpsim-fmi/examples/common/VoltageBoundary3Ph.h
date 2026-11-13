// SPDX-License-Identifier: MPL-2.0
#pragma once
#include <dpsim-models/MNASimPowerComp.h>
#include <functional>
#include <stdexcept>

/// Ideal prescribed voltage boundary, independent of the imported branch.
/// Adds three branch-current unknowns and voltage constraints to DPsim's MNA.
class VoltageBoundary3Ph final : public CPS::MNASimPowerComp<CPS::Real>,
                                 public SharedFactory<VoltageBoundary3Ph> {
  std::function<CPS::Matrix(CPS::Real)> mFunction;
  CPS::Matrix evaluate(CPS::Real time) const {
    CPS::Matrix value = mFunction(time);
    if (value.rows() != 3 || value.cols() != 1 || !value.allFinite())
      throw std::invalid_argument(
          "ABC voltage boundary requires a finite 3x1 voltage");
    return value;
  }

public:
  VoltageBoundary3Ph(CPS::String name,
                     std::function<CPS::Matrix(CPS::Real)> function)
      : MNASimPowerComp(name, true, true), mFunction(std::move(function)) {
    mPhaseType = CPS::PhaseType::ABC;
    setTerminalNumber(2);
    setVirtualNodeNumber(1);
    **mIntfVoltage = evaluate(0);
    **mIntfCurrent = CPS::Matrix::Zero(3, 1);
  }
  void initializeFromNodesAndTerminals(CPS::Real) override {}
  void mnaCompInitialize(CPS::Real, CPS::Real,
                         CPS::Attribute<CPS::Matrix>::Ptr) override {
    for (CPS::UInt index = 0; index < 2; ++index) {
      const auto terminalNode = node(index);
      if (!terminalNode || (!terminalNode->isGround() &&
                            terminalNode->phaseType() != CPS::PhaseType::ABC))
        throw std::invalid_argument(
            "ABC voltage boundary requires ABC terminal nodes or ground");
    }
    updateMatrixNodeIndices();
    **mIntfVoltage = evaluate(0);
  }
  void mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &matrix) override {
    for (const auto phase :
         {CPS::PhaseType::A, CPS::PhaseType::B, CPS::PhaseType::C}) {
      const auto branch = mVirtualNodes[0]->matrixNodeIndex(phase);
      for (CPS::UInt index = 0; index < 2; ++index) {
        if (terminalNotGrounded(index)) {
          const CPS::Real sign = index ? 1 : -1;
          CPS::Math::addToMatrixElement(
              matrix, matrixNodeIndex(index, static_cast<CPS::UInt>(phase)),
              branch, sign);
          CPS::Math::addToMatrixElement(
              matrix, branch,
              matrixNodeIndex(index, static_cast<CPS::UInt>(phase)), sign);
        }
      }
    }
  }
  void mnaCompPreStep(CPS::Real time, CPS::Int) override {
    **mIntfVoltage = evaluate(time);
    for (int channel = 0; channel < 3; ++channel) {
      const auto phase = static_cast<CPS::PhaseType>(channel);
      CPS::Math::setVectorElement(**mRightVector,
                                  mVirtualNodes[0]->matrixNodeIndex(phase),
                                  (**mIntfVoltage)(channel, 0));
    }
  }
  void mnaCompPostStep(CPS::Real, CPS::Int,
                       CPS::Attribute<CPS::Matrix>::Ptr &solution) override {
    for (int channel = 0; channel < 3; ++channel) {
      const auto phase = static_cast<CPS::PhaseType>(channel);
      (**mIntfCurrent)(channel, 0) = CPS::Math::realFromVectorElement(
          **solution, mVirtualNodes[0]->matrixNodeIndex(phase));
    }
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
