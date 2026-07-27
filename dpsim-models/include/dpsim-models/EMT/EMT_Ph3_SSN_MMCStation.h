// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>
#include <dpsim-models/Signal/DQSymControllerBlocks.h>

namespace CPS::EMT::Ph3 {

struct SSN_MMCStationParameters {
  Real nominalPower = 0.0;         // VA
  Real nominalAcLineLineRms = 0.0; // V
  Real nominalDcVoltage = 0.0;     // pole-to-pole V
  Real nominalFrequencyHz = 0.0;   // Hz
  Real controllerTimeStep = 0.0;   // s
  Real measurementFilterFrequencyHz = 1000.0;
  Real measurementFilterDamping = 1.0;
  Real armResistancePu = 0.0;
  Real armInductancePu = 0.0;
};

/// Signal-domain composition around one existing SSN_MMC plant.
///
/// This class owns no node, terminal, MNA stamp, or MMC electrical state.
/// Measurements at step k produce a held external plant command consumed by
/// SSN_MMC's next pre-step model update (k+1).
class SSN_MMCStation final : public SimSignalComp,
                             public SharedFactory<SSN_MMCStation> {
public:
  enum class ControlMode { ActivePowerReactivePower, DCVoltageReactivePower };
  enum class State { Blocked, Ready, Enabled };
  enum class EnableDiagnostic {
    None,
    NotReady,
    NonFiniteMeasurement,
    NonFiniteReference,
    ControlErrorTooLarge,
    CommandMismatch,
    InvalidPlantMode
  };

  const Attribute<Real>::Ptr mAngle;
  const Attribute<Real>::Ptr mAngularFrequency;
  const Attribute<Real>::Ptr mActivePowerReferencePu;
  const Attribute<Real>::Ptr mReactivePowerReferencePu;
  const Attribute<Real>::Ptr mDcVoltageReference;
  const Attribute<Real>::Ptr mVdPu;
  const Attribute<Real>::Ptr mVqPu;
  const Attribute<Real>::Ptr mIdPu;
  const Attribute<Real>::Ptr mIqPu;
  const Attribute<Real>::Ptr mActivePowerPu;
  const Attribute<Real>::Ptr mReactivePowerPu;
  const Attribute<Real>::Ptr mFilteredActivePowerPu;
  const Attribute<Real>::Ptr mFilteredReactivePowerPu;
  const Attribute<Real>::Ptr mFilteredDcVoltage;
  const Attribute<Real>::Ptr mIdReferencePu;
  const Attribute<Real>::Ptr mIqReferencePu;
  const Attribute<Real>::Ptr mVdReferencePu;
  const Attribute<Real>::Ptr mVqReferencePu;
  const Attribute<Matrix>::Ptr mConverterPhaseCommand;
  const Attribute<Matrix>::Ptr mPlantDifferentialVoltageCommand;
  const Attribute<Real>::Ptr mModulationMagnitude;
  const Attribute<Real>::Ptr mAcPower;
  const Attribute<Real>::Ptr mDcPower;
  const Attribute<Real>::Ptr mPowerBalanceError;
  const Attribute<Int>::Ptr mControlMode;
  const Attribute<Int>::Ptr mState;
  const Attribute<Int>::Ptr mEnableDiagnostic;
  const Attribute<Bool>::Ptr mControllerEnabled;
  const Attribute<Bool>::Ptr mOuterLoopsEnabled;
  const Attribute<Bool>::Ptr mCurrentDSaturated;
  const Attribute<Bool>::Ptr mCurrentQSaturated;
  const Attribute<Bool>::Ptr mModulationSaturated;
  const Attribute<Bool>::Ptr mCurrentDUpperSaturated;
  const Attribute<Bool>::Ptr mCurrentDLowerSaturated;
  const Attribute<Bool>::Ptr mCurrentQUpperSaturated;
  const Attribute<Bool>::Ptr mCurrentQLowerSaturated;
  const Attribute<Bool>::Ptr mModulationDUpperSaturated;
  const Attribute<Bool>::Ptr mModulationDLowerSaturated;
  const Attribute<Bool>::Ptr mModulationQUpperSaturated;
  const Attribute<Bool>::Ptr mModulationQLowerSaturated;
  const Attribute<Bool>::Ptr mActiveOuterUpperSaturated;
  const Attribute<Bool>::Ptr mActiveOuterLowerSaturated;
  const Attribute<Bool>::Ptr mDcOuterUpperSaturated;
  const Attribute<Bool>::Ptr mDcOuterLowerSaturated;
  const Attribute<Bool>::Ptr mReactiveOuterUpperSaturated;
  const Attribute<Bool>::Ptr mReactiveOuterLowerSaturated;

  SSN_MMCStation(String name, SSN_MMC::Ptr plant,
                 Logger::Level logLevel = Logger::Level::off);

  void setParameters(const SSN_MMCStationParameters &parameters);
  void setControlMode(ControlMode mode);
  void setReferences(Real activePowerPu, Real reactivePowerPu,
                     Real dcVoltageVolts);
  void setOuterLoopsEnabled(Bool enabled);
  void setCurrentReferences(Real idReferencePu, Real iqReferencePu);
  void initializeFromOperatingPoint(Real angle, Real angularFrequency,
                                    Real activePowerPu, Real reactivePowerPu,
                                    Real dcVoltageVolts);
  Bool requestEnable(Real errorTolerancePu = 1e-6,
                     Real commandToleranceVolts = 1e-3);
  void block();
  void initialize(Real timeStep) override;
  Task::List getTasks() override;

  SSN_MMC::Ptr plant() const { return mPlant; }

private:
  SSN_MMC::Ptr mPlant;
  SSN_MMCStationParameters mParameters;
  std::shared_ptr<Signal::ExternallyAngledDQAdapter> mTransform;
  std::shared_ptr<Signal::DQSymSecondOrderFilter> mPFilter;
  std::shared_ptr<Signal::DQSymSecondOrderFilter> mQFilter;
  std::shared_ptr<Signal::DQSymSecondOrderFilter> mVdcFilter;
  std::shared_ptr<Signal::DQSymOuterController> mActiveController;
  std::shared_ptr<Signal::DQSymOuterController> mDcVoltageController;
  std::shared_ptr<Signal::DQSymOuterController> mReactiveController;
  std::shared_ptr<Signal::DQSymCurrentController> mCurrentController;
  std::shared_ptr<Signal::DQSymModulation> mModulation;
  Matrix mHeldPlantCommand = Matrix::Zero(2, 1);
  Real mEnableErrorTolerancePu = 1e-6;
  Real mEnableCommandToleranceVolts = 1e-3;
  Bool mParametersSet = false;
  Bool mInitialized = false;

  void measurementStep();
  void filterAndOuterStep();
  void currentStep();
  void commandStep();
  void setEnableDiagnostic(EnableDiagnostic diagnostic);

  class MeasurementTask : public Task {
  public:
    explicit MeasurementTask(SSN_MMCStation &station);
    void execute(Real, Int) override { mStation.measurementStep(); }

  private:
    SSN_MMCStation &mStation;
  };
  class OuterTask : public Task {
  public:
    explicit OuterTask(SSN_MMCStation &station);
    void execute(Real, Int) override { mStation.filterAndOuterStep(); }

  private:
    SSN_MMCStation &mStation;
  };
  class CurrentTask : public Task {
  public:
    explicit CurrentTask(SSN_MMCStation &station);
    void execute(Real, Int) override { mStation.currentStep(); }

  private:
    SSN_MMCStation &mStation;
  };
  class CommandTask : public Task {
  public:
    explicit CommandTask(SSN_MMCStation &station);
    void execute(Real, Int) override { mStation.commandStep(); }

  private:
    SSN_MMCStation &mStation;
  };
};

} // namespace CPS::EMT::Ph3
