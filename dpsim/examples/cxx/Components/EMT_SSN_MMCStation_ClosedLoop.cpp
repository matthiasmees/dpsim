// SPDX-License-Identifier: MPL-2.0
#include "../Examples.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMCStation.h>

using namespace CPS;
using namespace DPsim;

namespace {
MatrixComp balancedVoltageReference(Real phaseVoltageAmplitude) {
  MatrixComp reference = MatrixComp::Zero(3, 1);
  const Complex phaseA(phaseVoltageAmplitude / RMS3PH_TO_PEAK1PH, 0.0);
  reference(0, 0) = phaseA;
  reference(1, 0) = phaseA * SHIFT_TO_PHASE_B;
  reference(2, 0) = phaseA * SHIFT_TO_PHASE_C;
  return reference;
}

void require(Bool condition, const String &message) {
  if (!condition)
    throw std::runtime_error(message);
}

void requireFinite(const Matrix &value, const String &name) {
  if (!value.allFinite())
    throw std::runtime_error(name + " contains NaN or Inf.");
}

struct Metrics {
  Real maxAcKcl = 0.0;
  Real maxDcPositiveKcl = 0.0;
  Real maxDcNegativeKcl = 0.0;
  Real minEnergy = std::numeric_limits<Real>::infinity();
  Real maxEnergy = -std::numeric_limits<Real>::infinity();
  Real maxCurrentPu = 0.0;
  Real maxCommandJump = 0.0;
  Real initialCurrentNorm = 0.0;
  Real postEnableCurrentNorm = 0.0;
  Real positiveIdPeak = -std::numeric_limits<Real>::infinity();
  Real negativeIdPeak = std::numeric_limits<Real>::infinity();
  Real positiveIqPeak = -std::numeric_limits<Real>::infinity();
  Real negativeIqPeak = std::numeric_limits<Real>::infinity();
  Real finalId = 0.0;
  Real finalIq = 0.0;
  Bool saturationEntered = false;
  Bool saturationRecovered = false;
  Real positiveQPeak = -std::numeric_limits<Real>::infinity();
  Real negativeQPeak = std::numeric_limits<Real>::infinity();
  Real maxReactiveActivePowerDisturbancePu = 0.0;
  Real reactiveFinalQ = 0.0;
  Real positivePPeak = -std::numeric_limits<Real>::infinity();
  Real negativePPeak = std::numeric_limits<Real>::infinity();
  Real activeFinalP = 0.0;
};
} // namespace

int main() {
  const Real timeStep = 40e-6;
  const Real finalTime = 0.24;
  const Real enableTime = 0.005;
  const Real frequency = 50.0;
  const Real omega = 2.0 * PI * frequency;
  const Real acVoltageAmplitude = 345e3;
  const Real dcVoltage = 440e3;
  const Real nominalPower = 1e9;

  auto acNode = SimNode<Real>::make("ac", PhaseType::ABC);
  auto dcPositive = SimNode<Real>::make("dcPositive", PhaseType::DC);
  auto dcNegative = SimNode<Real>::make("dcNegative", PhaseType::DC);
  acNode->setInitialVoltage(
      Complex(acVoltageAmplitude / RMS3PH_TO_PEAK1PH, 0.0));
  dcPositive->setInitialVoltage(Complex(dcVoltage / 2.0, 0.0));
  dcNegative->setInitialVoltage(Complex(-dcVoltage / 2.0, 0.0));

  auto acSource = EMT::Ph3::NetworkInjection::make("acSource");
  acSource->setParameters(balancedVoltageReference(acVoltageAmplitude),
                          frequency);
  acSource->connect({acNode});

  auto positiveSource = EMT::DC::VoltageSource::make("positiveSource");
  positiveSource->setParameters(dcVoltage / 2.0);
  positiveSource->connect({SimNode<Real>::GND, dcPositive});
  auto negativeSource = EMT::DC::VoltageSource::make("negativeSource");
  negativeSource->setParameters(dcVoltage / 2.0);
  negativeSource->connect({dcNegative, SimNode<Real>::GND});

  auto mmc = EMT::Ph3::SSN_MMC::make("mmc");
  mmc->setParameters(frequency, acVoltageAmplitude, dcVoltage, 0.05, 1.07, 0.01,
                     400, 0.0005, 0.0001);
  mmc->setInitialAngle(0.0);
  mmc->setPLL(0.0, 0.0, false);
  mmc->setActiveControlOpenLoop(0.0);
  mmc->setReactiveControlOpenLoop(0.0);
  mmc->setOutputCurrentController(117.93, 8.5e4);
  mmc->setCirculatingCurrentController(19.93, 4500.0);
  mmc->setZeroSequenceCurrentController(19.93, 4500.0);
  mmc->setCirculatingCurrentReferences(0.0, 0.0, 0.0);
  mmc->setLimits(1000.0, 1000.0, 2.0);
  mmc->setOperatingPointInitialization(true, 50, 1e-8);
  mmc->connect({acNode, dcPositive, dcNegative});

  auto station = EMT::Ph3::SSN_MMCStation::make("station", mmc);
  EMT::Ph3::SSN_MMCStationParameters parameters;
  parameters.nominalPower = nominalPower;
  parameters.nominalAcLineLineRms = acVoltageAmplitude;
  parameters.nominalDcVoltage = dcVoltage;
  parameters.nominalFrequencyHz = frequency;
  parameters.controllerTimeStep = timeStep;
  parameters.armResistancePu = 0.0015;
  parameters.armInductancePu = 0.15;
  station->setParameters(parameters);
  station->setControlMode(
      EMT::Ph3::SSN_MMCStation::ControlMode::DCVoltageReactivePower);
  station->setOuterLoopsEnabled(false);
  station->mAngle->set(0.0);
  station->mAngularFrequency->set(omega);

  SystemTopology system(frequency,
                        SystemNodeList{acNode, dcPositive, dcNegative},
                        SystemComponentList{acSource, positiveSource,
                                            negativeSource, mmc, station});
  Simulation sim("EMT_SSN_MMCStation_ClosedLoop", Logger::Level::off);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.initialize();

  require(mmc->controlSource() ==
              EMT::Ph3::SSN_MMC::ControlSource::InternalControllers,
          "MMC external command became active during initialization.");
  require(**station->mState ==
              static_cast<Int>(EMT::Ph3::SSN_MMCStation::State::Ready),
          "Station did not initialize into Ready.");

  Metrics metrics;
  metrics.initialCurrentNorm =
      mmc->getInterfaceCurrent().block(0, 0, 3, 1).norm();
  Matrix commandBeforeEnable = **station->mPlantDifferentialVoltageCommand;
  Bool enabled = false;
  Bool capturedPostEnable = false;
  Bool positiveIdApplied = false;
  Bool negativeIdApplied = false;
  Bool idRestored = false;
  Bool positiveIqApplied = false;
  Bool negativeIqApplied = false;
  Bool iqRestored = false;
  Bool saturationApplied = false;
  Bool saturationRemoved = false;
  Bool reactiveModeEnabled = false;
  Bool positiveQApplied = false;
  Bool negativeQApplied = false;
  Bool qReferenceRestored = false;
  Bool activeModeEnabled = false;
  Bool positivePApplied = false;
  Bool negativePApplied = false;
  Bool pReferenceRestored = false;

  sim.start();
  while (sim.time() < sim.finalTime()) {
    station->mAngle->set(omega * sim.time());
    if (!enabled && sim.time() >= enableTime) {
      require(station->requestEnable(2e-5, 1.0),
              "Initialized external-command enable was rejected; diagnostic=" +
                  std::to_string(**station->mEnableDiagnostic));
      enabled = true;
      commandBeforeEnable = **station->mPlantDifferentialVoltageCommand;
    }
    if (enabled && !positiveIdApplied && sim.time() >= 0.010) {
      station->setCurrentReferences(0.01, 0.0);
      positiveIdApplied = true;
    }
    if (enabled && !negativeIdApplied && sim.time() >= 0.014) {
      station->setCurrentReferences(-0.01, 0.0);
      negativeIdApplied = true;
    }
    if (enabled && !idRestored && sim.time() >= 0.018) {
      station->setCurrentReferences(0.0, 0.0);
      idRestored = true;
    }
    if (enabled && !positiveIqApplied && sim.time() >= 0.022) {
      station->setCurrentReferences(0.0, 0.01);
      positiveIqApplied = true;
    }
    if (enabled && !negativeIqApplied && sim.time() >= 0.026) {
      station->setCurrentReferences(0.0, -0.01);
      negativeIqApplied = true;
    }
    if (enabled && !iqRestored && sim.time() >= 0.030) {
      station->setCurrentReferences(0.0, 0.0);
      iqRestored = true;
    }
    if (enabled && !saturationApplied && sim.time() >= 0.034) {
      station->setCurrentReferences(3.0, 0.0);
      saturationApplied = true;
    }
    if (enabled && !saturationRemoved && sim.time() >= 0.0342) {
      station->setCurrentReferences(0.0, 0.0);
      saturationRemoved = true;
    }
    if (enabled && !reactiveModeEnabled && sim.time() >= 0.085) {
      station->block();
      station->setOuterLoopsEnabled(true);
      station->setControlMode(
          EMT::Ph3::SSN_MMCStation::ControlMode::DCVoltageReactivePower);
      station->initializeFromOperatingPoint(
          omega * sim.time(), omega, **station->mActivePowerPu,
          **station->mReactivePowerPu, **mmc->dcVoltageAttribute());
      station->initialize(timeStep);
      require(station->requestEnable(2e-5, 1.0),
              "Reactive-loop enable was rejected; diagnostic=" +
                  std::to_string(**station->mEnableDiagnostic));
      reactiveModeEnabled = true;
    }
    if (reactiveModeEnabled && !positiveQApplied && sim.time() >= 0.095) {
      station->mReactivePowerReferencePu->set(0.01);
      positiveQApplied = true;
    }
    if (reactiveModeEnabled && !negativeQApplied && sim.time() >= 0.110) {
      station->mReactivePowerReferencePu->set(-0.01);
      negativeQApplied = true;
    }
    if (reactiveModeEnabled && !qReferenceRestored && sim.time() >= 0.125) {
      station->mReactivePowerReferencePu->set(0.0);
      qReferenceRestored = true;
    }
    if (reactiveModeEnabled && !activeModeEnabled && sim.time() >= 0.155) {
      station->block();
      station->setControlMode(
          EMT::Ph3::SSN_MMCStation::ControlMode::ActivePowerReactivePower);
      station->initializeFromOperatingPoint(
          omega * sim.time(), omega, **station->mActivePowerPu,
          **station->mReactivePowerPu, **mmc->dcVoltageAttribute());
      station->initialize(timeStep);
      require(station->requestEnable(2e-5, 1.0),
              "Active-power mode enable was rejected; diagnostic=" +
                  std::to_string(**station->mEnableDiagnostic));
      station->mReactivePowerReferencePu->set(0.0);
      activeModeEnabled = true;
    }
    if (activeModeEnabled && !positivePApplied && sim.time() >= 0.165) {
      station->mActivePowerReferencePu->set(0.01);
      positivePApplied = true;
    }
    if (activeModeEnabled && !negativePApplied && sim.time() >= 0.185) {
      station->mActivePowerReferencePu->set(-0.01);
      negativePApplied = true;
    }
    if (activeModeEnabled && !pReferenceRestored && sim.time() >= 0.205) {
      station->mActivePowerReferencePu->set(0.0);
      pReferenceRestored = true;
    }

    sim.step();

    requireFinite(mmc->getState(), "MMC state");
    requireFinite(mmc->getInterfaceVoltage(), "MMC interface voltage");
    requireFinite(mmc->getInterfaceCurrent(), "MMC interface current");
    requireFinite(**station->mPlantDifferentialVoltageCommand,
                  "station plant command");

    const Matrix acKcl =
        acSource->intfCurrent() + mmc->getInterfaceCurrent().block(0, 0, 3, 1);
    const Real dcPositiveKcl =
        mmc->getInterfaceCurrent()(3, 0) + positiveSource->intfCurrent()(0, 0);
    const Real dcNegativeKcl =
        mmc->getInterfaceCurrent()(4, 0) - negativeSource->intfCurrent()(0, 0);
    metrics.maxAcKcl = std::max(metrics.maxAcKcl, acKcl.cwiseAbs().maxCoeff());
    metrics.maxDcPositiveKcl =
        std::max(metrics.maxDcPositiveKcl, std::abs(dcPositiveKcl));
    metrics.maxDcNegativeKcl =
        std::max(metrics.maxDcNegativeKcl, std::abs(dcNegativeKcl));
    const Real energy = **mmc->storedEnergyAttribute();
    metrics.minEnergy = std::min(metrics.minEnergy, energy);
    metrics.maxEnergy = std::max(metrics.maxEnergy, energy);
    metrics.maxCurrentPu = std::max(
        metrics.maxCurrentPu, std::hypot(**station->mIdPu, **station->mIqPu));
    if (sim.time() >= 0.010 && sim.time() < 0.014)
      metrics.positiveIdPeak =
          std::max(metrics.positiveIdPeak, **station->mIdPu);
    if (sim.time() >= 0.014 && sim.time() < 0.018)
      metrics.negativeIdPeak =
          std::min(metrics.negativeIdPeak, **station->mIdPu);
    if (sim.time() >= 0.022 && sim.time() < 0.026)
      metrics.positiveIqPeak =
          std::max(metrics.positiveIqPeak, **station->mIqPu);
    if (sim.time() >= 0.026 && sim.time() < 0.030)
      metrics.negativeIqPeak =
          std::min(metrics.negativeIqPeak, **station->mIqPu);
    metrics.saturationEntered =
        metrics.saturationEntered || **station->mCurrentDSaturated;
    if (saturationRemoved && !**station->mCurrentDSaturated)
      metrics.saturationRecovered = true;
    metrics.finalId = **station->mIdPu;
    metrics.finalIq = **station->mIqPu;
    if (sim.time() >= 0.095 && sim.time() < 0.110) {
      metrics.positiveQPeak =
          std::max(metrics.positiveQPeak, **station->mReactivePowerPu);
      metrics.maxReactiveActivePowerDisturbancePu =
          std::max(metrics.maxReactiveActivePowerDisturbancePu,
                   std::abs(**station->mActivePowerPu));
    }
    if (sim.time() >= 0.110 && sim.time() < 0.125) {
      metrics.negativeQPeak =
          std::min(metrics.negativeQPeak, **station->mReactivePowerPu);
      metrics.maxReactiveActivePowerDisturbancePu =
          std::max(metrics.maxReactiveActivePowerDisturbancePu,
                   std::abs(**station->mActivePowerPu));
    }
    if (qReferenceRestored)
      metrics.reactiveFinalQ = **station->mReactivePowerPu;
    if (sim.time() >= 0.165 && sim.time() < 0.185)
      metrics.positivePPeak =
          std::max(metrics.positivePPeak, **station->mActivePowerPu);
    if (sim.time() >= 0.185 && sim.time() < 0.205)
      metrics.negativePPeak =
          std::min(metrics.negativePPeak, **station->mActivePowerPu);
    if (pReferenceRestored)
      metrics.activeFinalP = **station->mActivePowerPu;

    if (enabled && !capturedPostEnable) {
      metrics.postEnableCurrentNorm =
          mmc->getInterfaceCurrent().block(0, 0, 3, 1).norm();
      metrics.maxCommandJump =
          (**station->mPlantDifferentialVoltageCommand - commandBeforeEnable)
              .norm();
      capturedPostEnable = true;
    }
  }
  sim.stop();

  SPDLOG_INFO(
      "Stage3B precheck: command_jump={} V, initial_Iabc_norm={} A, "
      "post_enable_Iabc_norm={} A, max_AC_KCL={} A, max_dc+_KCL={} A, "
      "max_dc-_KCL={} A, energy=[{},{}] J, max_current={} pu, "
      "id_peaks=[{},{}] pu, iq_peaks=[{},{}] pu, final_idq=[{},{}] pu, "
      "saturation_entered={}, saturation_recovered={}, Q_peaks=[{},{}] pu, "
      "Q_final={} pu, max_reactive_P_disturbance={} pu, P_peaks=[{},{}] pu, "
      "P_final={} pu",
      metrics.maxCommandJump, metrics.initialCurrentNorm,
      metrics.postEnableCurrentNorm, metrics.maxAcKcl, metrics.maxDcPositiveKcl,
      metrics.maxDcNegativeKcl, metrics.minEnergy, metrics.maxEnergy,
      metrics.maxCurrentPu, metrics.negativeIdPeak, metrics.positiveIdPeak,
      metrics.negativeIqPeak, metrics.positiveIqPeak, metrics.finalId,
      metrics.finalIq, metrics.saturationEntered, metrics.saturationRecovered,
      metrics.negativeQPeak, metrics.positiveQPeak, metrics.reactiveFinalQ,
      metrics.maxReactiveActivePowerDisturbancePu, metrics.negativePPeak,
      metrics.positivePPeak, metrics.activeFinalP);
  require(enabled && capturedPostEnable,
          "External command transition was not exercised.");
  require(metrics.maxAcKcl < 1e-6, "AC KCL residual exceeds 1e-6 A.");
  require(metrics.maxDcPositiveKcl < 1e-6, "dc+ KCL residual exceeds 1e-6 A.");
  require(metrics.maxDcNegativeKcl < 1e-6, "dc- KCL residual exceeds 1e-6 A.");
  require(metrics.maxCommandJump < 5.0,
          "External-command transition exceeds 5 V dq tolerance.");
  require(std::abs(metrics.postEnableCurrentNorm - metrics.initialCurrentNorm) <
              0.05,
          "External-command transition caused an AC-current discontinuity.");
  require(metrics.minEnergy > 0.0 &&
              metrics.maxEnergy / metrics.minEnergy < 1.05,
          "MMC stored energy is not bounded around the initialized point.");
  require(metrics.positiveIdPeak > 0.002,
          "Positive id reference did not produce positive plant current.");
  require(metrics.negativeIdPeak < -0.002,
          "Negative id reference did not produce negative plant current.");
  require(metrics.positiveIqPeak > 0.002,
          "Positive iq reference did not produce positive plant current.");
  require(metrics.negativeIqPeak < -0.002,
          "Negative iq reference did not produce negative plant current.");
  require(std::hypot(metrics.finalId, metrics.finalIq) < 0.002,
          "Inner current loop did not return to the original operating point.");
  require(metrics.saturationEntered && metrics.saturationRecovered,
          "Current-controller saturation entry/recovery was not observed.");
  require(metrics.positiveQPeak > 0.001,
          "Positive Q reference did not produce positive EMT reactive power.");
  require(metrics.negativeQPeak < -0.001,
          "Negative Q reference did not produce negative EMT reactive power.");
  require(std::abs(metrics.reactiveFinalQ) < 0.002,
          "Reactive-power loop did not recover after reference restoration.");
  require(metrics.maxReactiveActivePowerDisturbancePu < 0.01,
          "Reactive-power transient caused excessive active-power coupling.");
  require(metrics.positivePPeak > 0.002,
          "Positive P reference did not produce positive EMT active power.");
  require(metrics.negativePPeak < -0.002,
          "Negative P reference did not produce negative EMT active power.");
  require(
      std::abs(metrics.activeFinalP) < 0.002,
      "Active-power hierarchy did not recover after reference restoration.");

  SPDLOG_INFO(
      "Stage3B: command_jump={} V, initial_Iabc_norm={} A, "
      "post_enable_Iabc_norm={} A, max_AC_KCL={} A, max_dc+_KCL={} A, "
      "max_dc-_KCL={} A, energy=[{},{}] J, max_current={} pu, "
      "id_peaks=[{},{}] pu, iq_peaks=[{},{}] pu, final_idq=[{},{}] pu, "
      "saturation_entered={}, saturation_recovered={}",
      metrics.maxCommandJump, metrics.initialCurrentNorm,
      metrics.postEnableCurrentNorm, metrics.maxAcKcl, metrics.maxDcPositiveKcl,
      metrics.maxDcNegativeKcl, metrics.minEnergy, metrics.maxEnergy,
      metrics.maxCurrentPu, metrics.negativeIdPeak, metrics.positiveIdPeak,
      metrics.negativeIqPeak, metrics.positiveIqPeak, metrics.finalId,
      metrics.finalIq, metrics.saturationEntered, metrics.saturationRecovered);
  return 0;
}
