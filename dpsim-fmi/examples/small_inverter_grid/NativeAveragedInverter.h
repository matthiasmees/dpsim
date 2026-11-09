// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <array>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/PtrFactory.h>
#include <dpsim-models/Solver/MNANonlinearCompInterface.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>

namespace SmallGrid {
/// Native DPsim implementation of the averaged inverter used by the small-grid
/// FMI. It uses the same P/Q current reference, 1.2 pu limit, passive-current
/// sign, and 5 ms first-order current response as IEEE39AvgInverter.fmu.
class NativeAveragedInverter final
    : public CPS::MNASimPowerComp<CPS::Real>,
      public CPS::MNAVariableCompInterface,
      public CPS::MNANonlinearCompInterface,
      public SharedFactory<NativeAveragedInverter> {
  double mStep = 0;
  CPS::Matrix mInitialVoltage;
  CPS::Matrix mOldState, mOldTarget;
  CPS::Matrix mTrialState, mTrialVoltage;
  CPS::Matrix mW, mHistory;
  CPS::Attribute<CPS::Matrix>::Ptr mState;
  CPS::Attribute<CPS::Real>::Ptr mIterationAttribute, mResidualAttribute;
  int mIterations = 0;
  double mResidual = 0;
  bool mUsePeriodicFilterStart = true;
  double mActivePowerW = 125e6;
  double mReactivePowerVar = 80.881e6;
  double mNominalLineVoltageV = 345e3;
  double mCurrentLimitPu = 1.20;
  double mTauSeconds = 5e-3;
  double mFrequencyHz = 60.0;
  static constexpr int kMaximumNewtonIterations = 20;
  static constexpr double kAbsoluteCurrentTolerance = 1e-9;
  static constexpr double kRelativeCurrentTolerance = 1e-10;
  static constexpr double kJacobianPerturbation = 6e-6;

  CPS::Matrix passiveCurrentReference(const CPS::Matrix &voltage) const;
  CPS::Matrix periodicSteadyCurrent(const CPS::Matrix &voltage) const;
  CPS::Matrix endpointCurrent(const CPS::Matrix &voltage) const;
  CPS::Matrix voltageFromSolution(const CPS::Matrix &solution);
  void validateTerminals();
  void prepareTangent();

public:
  using CPS::SimPowerComp<CPS::Real>::connect;
  explicit NativeAveragedInverter(CPS::String name);
  void setUsePeriodicFilterStart(bool enabled) {
    mUsePeriodicFilterStart = enabled;
  }
  void setOperatingPoint(double activePowerW, double reactivePowerVar,
                         double nominalLineVoltageV,
                         double currentLimitPu = 1.20,
                         double responseTimeS = 5e-3,
                         double frequencyHz = 60.0);
  CPS::SimPowerComp<CPS::Real>::Ptr clone(CPS::String) override;
  void connect(CPS::SimNode<CPS::Real>::Ptr node);
  void initialize(CPS::Matrix frequencies) override;
  void initializeFromNodesAndTerminals(CPS::Real frequency) override;
  void mnaCompInitialize(CPS::Real, CPS::Real,
                         CPS::Attribute<CPS::Matrix>::Ptr) override;
  void mnaCompApplySystemMatrixStamp(CPS::SparseMatrixRow &) override;
  void mnaCompApplyRightSideVectorStamp(CPS::Matrix &) override;
  void mnaCompUpdateVoltage(const CPS::Matrix &) override;
  void mnaCompUpdateCurrent(const CPS::Matrix &) override;
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
  void mnaPrepareNonlinearIteration() override { prepareTangent(); }
};
} // namespace SmallGrid
