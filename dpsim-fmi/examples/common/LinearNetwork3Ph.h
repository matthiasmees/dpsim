// SPDX-License-Identifier: MPL-2.0
#pragma once
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>
#include <functional>
#include <stdexcept>

/// Six-state external grid: C*v' = i_port-i_line-G*v;
/// L*i_line' = v-v_grid-R*i_line. It is a current-input SSN branch;
/// its three branch-current unknowns participate directly in network MNA.
/// Parameters belong to the surrounding grid, not the imported FMU.
class LinearNetwork3Ph final : public CPS::MNASimPowerComp<CPS::Real>,
                               public CPS::MNAVariableCompInterface,
                               public SharedFactory<LinearNetwork3Ph> {
  double mC, mR, mL, mStep = 0, mG = -1;
  bool mChanged = true;
  CPS::Matrix mInitialState, mInitialCurrent, mAd, mBd, mForce, mW, mHistory;
  std::function<CPS::Matrix(double)> mGrid;
  std::function<double(double)> mLoad;
  CPS::Attribute<CPS::Matrix>::Ptr mState;
  CPS::Attribute<CPS::Real>::Ptr mLoadAttribute, mLoadPower, mLinePower,
      mLineLoss;
  void restore() {
    **mState = mInitialState;
    **mIntfVoltage = mInitialState.topRows(3);
    **mIntfCurrent = mInitialCurrent;
    measure(0);
  }
  void measure(double t) {
    const auto x = **mState;
    **mLoadAttribute = mLoad(t);
    **mLoadPower = mLoad(t) * x.topRows(3).squaredNorm();
    **mLinePower = x.col(0).head(3).dot(x.col(0).tail(3));
    **mLineLoss = mR * x.bottomRows(3).squaredNorm();
  }
  void prepare(double t) {
    const double g = mLoad(t - mStep / 2);
    if (!std::isfinite(g) || g < 0)
      throw std::invalid_argument("Invalid load conductance");
    mChanged = (g != mG);
    mG = g;
    CPS::Matrix a = CPS::Matrix::Zero(6, 6), b = CPS::Matrix::Zero(6, 3);
    a.topLeftCorner(3, 3) = (-g / mC) * CPS::Matrix::Identity(3, 3);
    a.topRightCorner(3, 3) = (-1 / mC) * CPS::Matrix::Identity(3, 3);
    a.bottomLeftCorner(3, 3) = (1 / mL) * CPS::Matrix::Identity(3, 3);
    a.bottomRightCorner(3, 3) = (-mR / mL) * CPS::Matrix::Identity(3, 3);
    b.topRows(3) = (1 / mC) * CPS::Matrix::Identity(3, 3);
    const CPS::Matrix identity = CPS::Matrix::Identity(6, 6);
    const CPS::Matrix left = identity - mStep * a / 2;
    auto lu = left.fullPivLu();
    mAd = lu.solve(identity + mStep * a / 2);
    mBd = lu.solve(mStep * b / 2);
    CPS::Matrix forcing = CPS::Matrix::Zero(6, 1);
    forcing.bottomRows(3) = -(mGrid(t - mStep) + mGrid(t)) / mL;
    mForce = lu.solve(mStep * forcing / 2);
    mW = mBd.topRows(3);
    const CPS::Matrix history =
        mAd * (**mState) + mBd * (**mIntfCurrent) + mForce;
    mHistory = history.topRows(3);
  }

public:
  LinearNetwork3Ph(CPS::String name, double capacitance, double resistance,
                   double inductance, CPS::Matrix initialState,
                   CPS::Matrix initialCurrent,
                   std::function<CPS::Matrix(double)> grid,
                   std::function<double(double)> load)
      : MNASimPowerComp(name, true, true), mC(capacitance), mR(resistance),
        mL(inductance), mInitialState(std::move(initialState)),
        mInitialCurrent(std::move(initialCurrent)), mGrid(std::move(grid)),
        mLoad(std::move(load)) {
    if (!std::isfinite(mC) || !std::isfinite(mR) || !std::isfinite(mL) ||
        mC <= 0 || mR < 0 || mL <= 0 || mInitialState.rows() != 6 ||
        mInitialState.cols() != 1 || mInitialCurrent.rows() != 3 ||
        mInitialCurrent.cols() != 1 || !mInitialState.allFinite() ||
        !mInitialCurrent.allFinite())
      throw std::invalid_argument(
          "Invalid external grid parameters or initialization");
    mPhaseType = CPS::PhaseType::ABC;
    setTerminalNumber(2);
    setVirtualNodeNumber(1);
    mState = mAttributes->create<CPS::Matrix>("x");
    mLoadAttribute = mAttributes->create<CPS::Real>("g_load_S");
    mLoadPower = mAttributes->create<CPS::Real>("p_load_W");
    mLinePower = mAttributes->create<CPS::Real>("p_line_bus_W");
    mLineLoss = mAttributes->create<CPS::Real>("p_line_loss_W");
    for (int p = 0; p < 3; ++p) {
      const std::string phase(1, "abc"[p]);
      auto v = mAttributes->createDynamic<CPS::Real>("v_bus_" + phase + "_V");
      auto i = mAttributes->createDynamic<CPS::Real>("i_line_" + phase + "_A");
      v->setReference(mState->deriveCoeff<CPS::Real>(p, 0));
      i->setReference(mState->deriveCoeff<CPS::Real>(p + 3, 0));
    }
    restore();
  }
  void initialize(CPS::Matrix f) override {
    CPS::SimPowerComp<CPS::Real>::initialize(std::move(f));
    restore();
  }
  void initializeFromNodesAndTerminals(CPS::Real) override { restore(); }
  void mnaCompInitialize(CPS::Real, CPS::Real step,
                         CPS::Attribute<CPS::Matrix>::Ptr) override {
    mStep = step;
    updateMatrixNodeIndices();
    for (CPS::UInt p = 0; p < 2; ++p)
      if (!node(p) ||
          (!node(p)->isGround() && node(p)->phaseType() != CPS::PhaseType::ABC))
        throw std::invalid_argument("External grid needs ABC nodes or ground");
    prepare(step);
    for (CPS::UInt i = 0; i < 3; ++i)
      for (CPS::UInt j = 0; j < 3; ++j)
        mVariableSystemMatrixEntries.emplace_back(
            mVirtualNodes[0]->matrixNodeIndex(static_cast<CPS::PhaseType>(i)),
            mVirtualNodes[0]->matrixNodeIndex(static_cast<CPS::PhaseType>(j)));
  }
  CPS::Bool hasParameterChanged() override { return mChanged; }
  void mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &matrix) override {
    for (CPS::UInt p = 0; p < 3; ++p) {
      const auto branch =
          mVirtualNodes[0]->matrixNodeIndex(static_cast<CPS::PhaseType>(p));
      for (CPS::UInt t = 0; t < 2; ++t)
        if (terminalNotGrounded(t)) {
          const double sign = t ? 1 : -1;
          CPS::Math::addToMatrixElement(matrix, matrixNodeIndex(t, p), branch,
                                        sign);
          CPS::Math::addToMatrixElement(matrix, branch, matrixNodeIndex(t, p),
                                        sign);
        }
      for (CPS::UInt q = 0; q < 3; ++q)
        CPS::Math::addToMatrixElement(
            matrix, branch,
            mVirtualNodes[0]->matrixNodeIndex(static_cast<CPS::PhaseType>(q)),
            -mW(p, q));
    }
  }
  void mnaCompApplyRightSideVectorStamp(CPS::Matrix &vector) override {
    for (CPS::UInt p = 0; p < 3; ++p)
      CPS::Math::setVectorElement(
          vector,
          mVirtualNodes[0]->matrixNodeIndex(static_cast<CPS::PhaseType>(p)),
          mHistory(p, 0));
  }
  void mnaCompPreStep(CPS::Real time, CPS::Int) override {
    prepare(time);
    mnaCompApplyRightSideVectorStamp(**mRightVector);
  }
  void mnaCompPostStep(CPS::Real time, CPS::Int,
                       CPS::Attribute<CPS::Matrix>::Ptr &solution) override {
    const CPS::Matrix oldCurrent = **mIntfCurrent;
    for (CPS::UInt p = 0; p < 3; ++p)
      (**mIntfCurrent)(p, 0) = CPS::Math::realFromVectorElement(
          **solution,
          mVirtualNodes[0]->matrixNodeIndex(static_cast<CPS::PhaseType>(p)));
    **mState = mAd * (**mState) + mBd * (oldCurrent + **mIntfCurrent) + mForce;
    **mIntfVoltage = (**mState).topRows(3);
    measure(time);
  }
  void
  mnaCompAddPreStepDependencies(CPS::AttributeBase::List &old,
                                CPS::AttributeBase::List &,
                                CPS::AttributeBase::List &modified) override {
    old.push_back(mState);
    old.push_back(mIntfCurrent);
    modified.push_back(mRightVector);
  }
  void mnaCompAddPostStepDependencies(
      CPS::AttributeBase::List &, CPS::AttributeBase::List &dependencies,
      CPS::AttributeBase::List &modified,
      CPS::Attribute<CPS::Matrix>::Ptr &solution) override {
    dependencies.push_back(solution);
    for (auto &a : {mState, mIntfVoltage, mIntfCurrent})
      modified.push_back(a);
    for (auto &a : {mLoadAttribute, mLoadPower, mLinePower, mLineLoss})
      modified.push_back(a);
  }
};
