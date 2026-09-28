// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <vector>

#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNAVariableCompInterface.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

/// Averaged MMC implemented with direct MNA companion equations.
///
/// The five inductive branch currents and seven capacitor-voltage harmonics
/// are auxiliary unknowns of the NETWORK solve. Inductor KVL, capacitor KCL,
/// modulation-dependent controlled sources, and the abc terminal mapping are
/// stamped directly; no continuous/discrete state-space matrices, SSN base
/// class, Norton reduction, or local state-transition matrices are used.
/// Algebraic controlled sources contribute incremental conductances directly
/// to the global circuit equations, including terminal-voltage feedforward.
///
/// The harmonic truncation and signs match SSN_MMC: iDelta exports AC power,
/// 3*iSigma_z absorbs DC current. Terminals are {ABC, DC+, DC-}. The internal
/// nominal dq frame is independent of the PLL controller frame.
/// PI integrators, PLL, measurement filters and optional modulation delay use
/// scalar predictor/corrector updates. The sampled P/Vd path retains its
/// exact held-input filter recurrence and sampling period.
class MMC_MNA final : public MNASimPowerComp<Real>,
                      public MNAVariableCompInterface,
                      public SharedFactory<MMC_MNA> {
public:
  using Ptr = std::shared_ptr<MMC_MNA>;
  using SharedFactory<MMC_MNA>::make;
  enum class ActiveControlMode {
    OpenLoop,
    ActivePower,
    ActivePowerFeedforward,
    DcVoltage,
    DcDroop
  };
  enum class ReactiveControlMode { OpenLoop, ReactivePower, AcVoltage };

  MMC_MNA(String uid, String name, Logger::Level level = Logger::Level::off);
  MMC_MNA(String name, Logger::Level level = Logger::Level::off)
      : MMC_MNA(name, name, level) {}

  void setParameters(Real frequency, Real acVoltage, Real dcVoltage,
                     Real armInductance, Real armResistance,
                     Real submoduleCapacitance, UInt submodules,
                     Real reactorInductance, Real reactorResistance);
  void setInitialAngle(Real angle);
  void setInitialOperatingPoint(Real activePower, Real reactivePower);
  void setOutputCurrentController(Real kp, Real ki);
  void setCirculatingCurrentController(Real kp, Real ki);
  void setZeroSequenceCurrentController(Real kp, Real ki);
  void setEnergyController(Real kp, Real ki, Bool enabled = true);
  void setPLL(Real kp, Real ki, Bool enabled = true);
  void setActivePowerControl(Real reference, Real kp, Real ki);
  void setActivePowerFeedforwardControl(Real reference, Real cutoffFrequency,
                                        Real sampleTime, Real minimumVd = 0.0);
  void setActivePowerFeedforwardReference(Real reference);
  void setDcVoltageControl(Real reference, Real kp, Real ki);
  void setDcDroopControl(Real power, Real voltage, Real droop);
  void setActiveControlOpenLoop(Real reference = 0.0);
  void setReactivePowerControl(Real reference, Real kp, Real ki);
  void setAcVoltageControl(Real reference, Real kp, Real ki);
  void setReactiveControlOpenLoop(Real reference = 0.0);
  void setCirculatingCurrentReferences(Real d, Real q, Real z);
  void setMeasurementFilters(Real acDq, Real p, Real q, Real dc,
                             Real acMagnitude);
  void setModulationDelay(Real delay, UInt padeOrder = 2);
  void setLimits(Real acCurrent, Real circulatingCurrent, Real modulation);
  void setTheta(Real theta);
  void setOperatingPointInitialization(Bool enabled, UInt iterations,
                                       Real tolerance);

  Attribute<Matrix>::Ptr acTerminalVoltageAttribute() const;
  Attribute<Matrix>::Ptr acTerminalCurrentAttribute() const;
  Attribute<Matrix>::Ptr appliedModulationAttribute() const {
    return mModulationAttribute;
  }
  Attribute<Matrix>::Ptr appliedDifferentialVoltageAttribute() const {
    return mDifferentialAttribute;
  }
  Attribute<Matrix>::Ptr appliedCommonModeVoltageAttribute() const {
    return mCommonAttribute;
  }
  Attribute<Real>::Ptr dcVoltageAttribute() const;
  Attribute<Real>::Ptr dcCurrentAttribute() const;
  Attribute<Real>::Ptr activePowerAttribute() const;
  Attribute<Real>::Ptr reactivePowerAttribute() const;
  Attribute<Real>::Ptr storedEnergyAttribute() const;

  void initializeFromNodesAndTerminals(Real frequency) override;
  void mnaCompInitialize(Real omega, Real step,
                         Attribute<Matrix>::Ptr leftVector) override;
  Bool hasParameterChanged() override { return true; }
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &matrix) override;
  void mnaCompApplyRightSideVectorStamp(Matrix &rhs) override;
  void mnaCompPreStep(Real time, Int step) override;
  void mnaCompPostStep(Real time, Int step,
                       Attribute<Matrix>::Ptr &leftVector) override;
  void mnaCompAddPreStepDependencies(AttributeBase::List &previous,
                                     AttributeBase::List &dependencies,
                                     AttributeBase::List &modified) override;
  void
  mnaCompAddPostStepDependencies(AttributeBase::List &previous,
                                 AttributeBase::List &dependencies,
                                 AttributeBase::List &modified,
                                 Attribute<Matrix>::Ptr &leftVector) override;

private:
  enum Branch : UInt {
    IDeltaD,
    IDeltaQ,
    ISigmaZ,
    ISigmaD,
    ISigmaQ,
    VCDeltaD,
    VCDeltaQ,
    VCDeltaZd,
    VCDeltaZq,
    VCSigmaD,
    VCSigmaQ,
    VCSigmaZ,
    BranchCount
  };
  enum Integrator : UInt {
    Active,
    Reactive,
    OccD,
    OccQ,
    CccD,
    CccQ,
    Zcc,
    Energy
  };
  using Electrical = std::array<Real, BranchCount>;
  using Modulation = std::array<Real, 5>;
  using Dq = std::array<Real, 2>;
  struct PIParameters {
    Real kp = 0.0, ki = 0.0;
  };
  struct Memory {
    std::array<Real, 8> integral{};
    Real pllIntegral = 0.0, pllAngle = 0.0;
    std::array<Real, 10> filter{};
    std::array<Real, 10> delay{};
  };
  struct Measurement {
    Real vd, vq, vdc;
  };
  struct Command {
    Memory rate;
    Modulation modulation{};
    Dq differential{};
    std::array<Real, 3> common{};
    Real idRef = 0.0, iqRef = 0.0, izRef = 0.0;
    Real pFiltered = 0.0, qFiltered = 0.0, vdcFiltered = 0.0;
    Real pllError = 0.0, deltaOmega = 0.0;
  };

  static Dq rotate(Dq value, Real angle);
  static Memory advance(const Memory &memory, const Memory &rate, Real step);
  static Real regularize(Real voltage);
  static PIParameters gains(Real kp, Real ki);
  Measurement measure(const Matrix &voltage, Real angle) const;
  Command control(const Electrical &electrical, const Memory &memory,
                  const Measurement &measurement) const;
  Real energy(const Electrical &electrical) const;
  /// Coefficients of physical KVL (volts) / KCL (amps), NOT derivatives.
  void
  electricalCouplings(const Modulation &modulation,
                      const std::function<void(UInt, UInt, Real)> &add) const;
  Electrical electricalBalance(const Electrical &electrical,
                               const Modulation &modulation,
                               const Measurement &measurement) const;
  void initializeOperatingPoint();
  void sampleFeedforward();
  void updateAttributes();
  void validateTerminals() const;
  void logScalar(const String &name, Real value);
  Attribute<Real>::Ptr scalar(const String &name) const;

  Real mOmega = 0.0, mAcVoltage = 0.0, mDcVoltage = 0.0;
  Real mArmL = 0.0, mArmR = 0.0, mCapacitance = 0.0;
  Real mAcL = 0.0, mAcR = 0.0;
  UInt mSubmodules = 0;
  Real mStep = 0.0, mTheta = 0.5, mAngle = 0.0, mInitialAngle = 0.0;
  Real mNextAngle = 0.0;
  std::uint64_t mStepCount = 0;
  Electrical mElectrical{}, mStorage{}, mHistory{};
  // Incremental conductances of nonlinear controlled sources with respect to
  // the 12 branch unknowns and the three measured terminal voltages (d,q,dc).
  std::array<std::array<Real, BranchCount + 3>, BranchCount>
      mSourceConductance{};
  Modulation mStepModulation{};
  Memory mMemory, mPredictedMemory, mOldRate;
  Measurement mMeasurement{0.0, 0.0, 0.0};

  ActiveControlMode mActiveMode = ActiveControlMode::OpenLoop;
  ReactiveControlMode mReactiveMode = ReactiveControlMode::OpenLoop;
  PIParameters mActive, mReactive, mOutput, mCirculating, mZero, mEnergy, mPll;
  Bool mEnergyEnabled = false, mPllEnabled = false;
  Real mPRef = 0.0, mQRef = 0.0, mVdcRef = 0.0, mVacRef = 0.0;
  Real mIdRef = 0.0, mIqRef = 0.0, mDroop = 0.0;
  std::array<Real, 3> mSigmaRef{};
  std::array<Real, 5> mFilterTime{};
  Real mDelay = 0.0;
  Real mMaxAc = std::numeric_limits<Real>::infinity();
  Real mMaxSigma = std::numeric_limits<Real>::infinity();
  Real mMaxModulation = std::numeric_limits<Real>::infinity();
  Real mSampleTime = 0.0, mCutoff = 0.0, mMinimumVd = 0.0;
  Real mFilteredVd = 0.0, mFilteredVdDerivative = 0.0, mHeldId = 0.0;
  UInt mSampleCounter = 0;
  Bool mLoadedStart = false, mSolveStart = true;
  Real mInitialP = 0.0, mInitialQ = 0.0, mInitialTolerance = 1e-9;
  UInt mInitialIterations = 100;
  std::array<Int, 5> mTerminalIndices{};
  std::array<UInt, BranchCount> mBranchIndices{};
  std::map<String, Attribute<Real>::Ptr> mScalars;
  Attribute<Matrix>::Ptr mElectricalAttribute, mAcVoltageAttribute;
  Attribute<Matrix>::Ptr mAcCurrentAttribute, mModulationAttribute;
  Attribute<Matrix>::Ptr mDifferentialAttribute, mCommonAttribute;
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
