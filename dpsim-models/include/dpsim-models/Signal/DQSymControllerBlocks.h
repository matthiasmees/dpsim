// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-models/SimSignalComp.h>
#include <dpsim-models/Task.h>

namespace CPS::Signal {

/// Explicit externally-angled transform boundary used by MMCStation.
///
/// Positive sequence is a-b-c, theta increases counter-clockwise, theta=0
/// aligns phase-a with the cosine-oriented d axis, and the amplitude-invariant
/// transform uses 2/3 scaling. Zero sequence is omitted. Voltage and current
/// abc inputs use peak instantaneous SI values; current is positive out of the
/// converter. P=1.5(vd*id+vq*iq), Q=1.5(vq*id-vd*iq).
/// This adapter is not claimed equivalent to MathWorks SPS abc/dq blocks.
class ExternallyAngledDQAdapter
    : public SimSignalComp,
      public SharedFactory<ExternallyAngledDQAdapter> {
public:
  const Attribute<Matrix>::Ptr mVoltageAbc;
  const Attribute<Matrix>::Ptr mCurrentAbc;
  const Attribute<Real>::Ptr mAngle;
  const Attribute<Real>::Ptr mAngularFrequency;
  const Attribute<Real>::Ptr mVd;
  const Attribute<Real>::Ptr mVq;
  const Attribute<Real>::Ptr mId;
  const Attribute<Real>::Ptr mIq;
  const Attribute<Real>::Ptr mActivePower;
  const Attribute<Real>::Ptr mReactivePower;

  ExternallyAngledDQAdapter(String name,
                            Logger::Level logLevel = Logger::Level::off);
  void step();
  Matrix dqToAbc(Real d, Real q) const;
  Task::List getTasks() override;

private:
  class StepTask : public Task {
  public:
    explicit StepTask(ExternallyAngledDQAdapter &adapter);
    void execute(Real, Int) override { mAdapter.step(); }

  private:
    ExternallyAngledDQAdapter &mAdapter;
  };
};

enum class DQSymOuterLoopType {
  ActivePowerToDcVoltage,
  DcVoltageToDCurrent,
  ReactivePowerToQCurrent
};

/// Scalar PI using DQsym's trapezoidal conditional-integration anti-windup.
/// Error, state, feedforward and outputs use one caller-selected unit system.
/// Feedforward is summed algebraically after the PI, as in the DQsym current
/// controller; it is not integrated. Disable holds the integrator state.
class DQSymPIController : public SimSignalComp,
                          public SharedFactory<DQSymPIController> {
public:
  const Attribute<Real>::Ptr mError;
  const Attribute<Real>::Ptr mFeedforward;
  const Attribute<Bool>::Ptr mEnable;
  const Attribute<Real>::Ptr mIntegratorState;
  const Attribute<Real>::Ptr mUnsaturatedOutput;
  const Attribute<Real>::Ptr mOutput;
  const Attribute<Bool>::Ptr mSaturated;
  const Attribute<Bool>::Ptr mUpperSaturated;
  const Attribute<Bool>::Ptr mLowerSaturated;

  DQSymPIController(String name, Logger::Level logLevel = Logger::Level::off);
  void setParameters(Real kp, Real ki, Real lowerLimit, Real upperLimit);
  void setInitialState(Real integratorState, Real initialError = 0.0,
                       Real initialFeedforward = 0.0);
  void initialize(Real timeStep) override;
  void step();
  Task::List getTasks() override;

private:
  Real mKp = 0.0;
  Real mKi = 0.0;
  Real mLowerLimit = 0.0;
  Real mUpperLimit = 0.0;
  Real mTimeStep = 0.0;
  Real mPreviousIntegratorInput = 0.0;
  Bool mParametersSet = false;
  Bool mInitialized = false;

  class StepTask : public Task {
  public:
    explicit StepTask(DQSymPIController &controller);
    void execute(Real, Int) override { mController.step(); }

  private:
    DQSymPIController &mController;
  };
};

class DQSymSecondOrderFilter : public SimSignalComp,
                               public SharedFactory<DQSymSecondOrderFilter> {
public:
  const Attribute<Real>::Ptr mInput;
  const Attribute<Matrix>::Ptr mState;
  const Attribute<Real>::Ptr mOutput;

  DQSymSecondOrderFilter(String name,
                         Logger::Level logLevel = Logger::Level::off);
  void setParameters(Real naturalFrequencyHz, Real dampingRatio);
  void setInitialValue(Real value);
  void initialize(Real timeStep) override;
  void step();
  Task::List getTasks() override;

private:
  Real mNaturalFrequencyHz = 0.0;
  Real mDampingRatio = 0.0;
  Real mTimeStep = 0.0;
  Matrix mPreviousState = Matrix::Zero(2, 1);
  Real mPreviousInput = 0.0;
  Bool mParametersSet = false;
  Bool mInitialized = false;

  class StepTask : public Task {
  public:
    explicit StepTask(DQSymSecondOrderFilter &filter);
    void execute(Real, Int) override { mFilter.step(); }

  private:
    DQSymSecondOrderFilter &mFilter;
  };
};

/// DQsym outer-loop cascade element.
/// Active and reactive power inputs are per unit. DC-voltage inputs are volts
/// and are divided by the configured nominal-voltage normalization. Outputs
/// are Vdc_ref in volts or id/iq references in per unit according to mType.
class DQSymOuterController : public SimSignalComp,
                             public SharedFactory<DQSymOuterController> {
public:
  const Attribute<Real>::Ptr mReference;
  const Attribute<Real>::Ptr mMeasurement;
  const Attribute<Bool>::Ptr mEnable;
  const Attribute<Real>::Ptr mError;
  const Attribute<Real>::Ptr mIntegratorState;
  const Attribute<Real>::Ptr mUnsaturatedOutput;
  const Attribute<Real>::Ptr mOutput;
  const Attribute<Bool>::Ptr mSaturated;
  const Attribute<Bool>::Ptr mUpperSaturated;
  const Attribute<Bool>::Ptr mLowerSaturated;

  DQSymOuterController(String name,
                       Logger::Level logLevel = Logger::Level::off);
  void setParameters(DQSymOuterLoopType type, Real kp, Real ki, Real lowerLimit,
                     Real upperLimit, Real normalization, Real outputScale,
                     Real initialIntegratorState);
  void setInitialIntegratorState(Real state);
  void initialize(Real timeStep) override;
  void step();
  Task::List getTasks() override;

private:
  DQSymOuterLoopType mType = DQSymOuterLoopType::DcVoltageToDCurrent;
  std::shared_ptr<DQSymPIController> mPI;
  Real mNormalization = 1.0;
  Real mOutputScale = 1.0;
  Bool mParametersSet = false;

  class StepTask : public Task {
  public:
    explicit StepTask(DQSymOuterController &controller);
    void execute(Real, Int) override { mController.step(); }

  private:
    DQSymOuterController &mController;
  };
};

struct DQSymCurrentControllerParameters {
  /// PI gains in the DQsym per-unit controller and seconds.
  Real kp = 0.0;
  Real ki = 0.0;
  /// Half-arm resistance and inductance feedforward, both per unit.
  Real rFeedforwardPu = 0.0;
  Real lFeedforwardPu = 0.0;
  /// Frequency base in hertz.
  Real nominalFrequencyHz = 0.0;
  /// Per-axis converter-voltage limits in per unit.
  Real lowerLimitPu = 0.0;
  Real upperLimitPu = 0.0;
};

/// DQ current controller. Currents and voltages are amplitude-based per-unit
/// quantities. Current is positive out of the converter. Positive id produces
/// positive generated P; negative iq produces positive generated Q.
class DQSymCurrentController : public SimSignalComp,
                               public SharedFactory<DQSymCurrentController> {
public:
  const Attribute<Real>::Ptr mIdReference;
  const Attribute<Real>::Ptr mIqReference;
  const Attribute<Real>::Ptr mId;
  const Attribute<Real>::Ptr mIq;
  const Attribute<Real>::Ptr mVd;
  const Attribute<Real>::Ptr mVq;
  const Attribute<Real>::Ptr mFrequencyHz;
  const Attribute<Bool>::Ptr mEnable;
  const Attribute<Real>::Ptr mDIntegratorState;
  const Attribute<Real>::Ptr mQIntegratorState;
  const Attribute<Real>::Ptr mDUnsaturatedReference;
  const Attribute<Real>::Ptr mQUnsaturatedReference;
  const Attribute<Real>::Ptr mVdReference;
  const Attribute<Real>::Ptr mVqReference;
  const Attribute<Bool>::Ptr mDSaturated;
  const Attribute<Bool>::Ptr mQSaturated;
  const Attribute<Bool>::Ptr mDUpperSaturated;
  const Attribute<Bool>::Ptr mDLowerSaturated;
  const Attribute<Bool>::Ptr mQUpperSaturated;
  const Attribute<Bool>::Ptr mQLowerSaturated;

  DQSymCurrentController(String name,
                         Logger::Level logLevel = Logger::Level::off);
  void setParameters(const DQSymCurrentControllerParameters &parameters);
  void setInitialIntegratorStates(Real dState, Real qState);
  void initialize(Real timeStep) override;
  void step();
  Task::List getTasks() override;

private:
  DQSymCurrentControllerParameters mParameters;
  Real mTimeStep = 0.0;
  Real mPreviousDIntegratorInput = 0.0;
  Real mPreviousQIntegratorInput = 0.0;
  Bool mParametersSet = false;
  Bool mInitialized = false;

  class StepTask : public Task {
  public:
    explicit StepTask(DQSymCurrentController &controller);
    void execute(Real, Int) override { mController.step(); }

  private:
    DQSymCurrentController &mController;
  };
};

class DQSymModulation : public SimSignalComp,
                        public SharedFactory<DQSymModulation> {
public:
  const Attribute<Real>::Ptr mVdCommand;
  const Attribute<Real>::Ptr mVqCommand;
  const Attribute<Real>::Ptr mDcVoltage;
  const Attribute<Real>::Ptr mAngle;
  const Attribute<Real>::Ptr mDCommand;
  const Attribute<Real>::Ptr mQCommand;
  const Attribute<Real>::Ptr mDUnsaturatedCommand;
  const Attribute<Real>::Ptr mQUnsaturatedCommand;
  const Attribute<Real>::Ptr mModulationMagnitude;
  const Attribute<Matrix>::Ptr mAbcCommand;
  const Attribute<Bool>::Ptr mSaturated;
  const Attribute<Bool>::Ptr mDUpperSaturated;
  const Attribute<Bool>::Ptr mDLowerSaturated;
  const Attribute<Bool>::Ptr mQUpperSaturated;
  const Attribute<Bool>::Ptr mQLowerSaturated;

  DQSymModulation(String name, Logger::Level logLevel = Logger::Level::off);
  /// DC base is pole-to-pole volts; AC base is line-line RMS volts. Axis
  /// commands and abc output are dimensionless per-unit reference quantities.
  void setParameters(Real nominalDcVoltage, Real nominalAcLineLineRms,
                     Real lowerAxisLimit = -2.0, Real upperAxisLimit = 2.0);
  void step();
  Task::List getTasks() override;

private:
  Real mNominalDcVoltage = 0.0;
  Real mNominalAcLineLineRms = 0.0;
  Real mLowerAxisLimit = -2.0;
  Real mUpperAxisLimit = 2.0;
  Bool mParametersSet = false;

  class StepTask : public Task {
  public:
    explicit StepTask(DQSymModulation &modulation);
    void execute(Real, Int) override { mModulation.step(); }

  private:
    DQSymModulation &mModulation;
  };
};

} // namespace CPS::Signal
