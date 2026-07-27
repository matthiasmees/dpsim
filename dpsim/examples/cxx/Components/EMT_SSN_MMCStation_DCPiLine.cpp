// SPDX-License-Identifier: MPL-2.0
#include "../Examples.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Capacitor.h>
#include <dpsim-models/EMT/EMT_DC_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_DC_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMCStation.h>

using namespace CPS;
using namespace DPsim;

namespace {
MatrixComp balancedVoltageReference(Real amplitude) {
  MatrixComp reference(3, 1);
  const Complex phaseA(amplitude / RMS3PH_TO_PEAK1PH, 0.0);
  reference << phaseA, phaseA * SHIFT_TO_PHASE_B, phaseA * SHIFT_TO_PHASE_C;
  return reference;
}
void require(Bool condition, const String &message) {
  if (!condition)
    throw std::runtime_error(message);
}
} // namespace

int main() {
  const Real timeStep = 40e-6;
  const Real finalTime = 0.30;
  const Real frequency = 50.0;
  const Real omega = 2.0 * PI * frequency;
  const Real acVoltage = 345e3;
  const Real nominalDcVoltage = 440e3;
  const Real nominalPower = 1e9;
  const Real localCapacitance = 150e-6;

  auto acNode = SimNode<Real>::make("ac", PhaseType::ABC);
  auto localPositive = SimNode<Real>::make("localPositive", PhaseType::DC);
  auto localNegative = SimNode<Real>::make("localNegative", PhaseType::DC);
  auto remotePositive = SimNode<Real>::make("remotePositive", PhaseType::DC);
  auto remoteNegative = SimNode<Real>::make("remoteNegative", PhaseType::DC);
  acNode->setInitialVoltage(Complex(acVoltage / RMS3PH_TO_PEAK1PH, 0.0));
  localPositive->setInitialVoltage(Complex(nominalDcVoltage / 2.0, 0.0));
  localNegative->setInitialVoltage(Complex(-nominalDcVoltage / 2.0, 0.0));
  remotePositive->setInitialVoltage(Complex(nominalDcVoltage / 2.0, 0.0));
  remoteNegative->setInitialVoltage(Complex(-nominalDcVoltage / 2.0, 0.0));

  auto acSource = EMT::Ph3::NetworkInjection::make("acSource");
  acSource->setParameters(balancedVoltageReference(acVoltage), frequency);
  acSource->connect({acNode});

  auto positiveLine = EMT::DC::SSN::PiLine::make("positiveLine");
  auto negativeLine = EMT::DC::SSN::PiLine::make("negativeLine");
  positiveLine->setParameters(0.25, 7.5e-3, 0.0, 0.0, 0.0);
  negativeLine->setParameters(0.25, 7.5e-3, 0.0, 0.0, 0.0);
  // Explicit validation-study setting. The component default remains the
  // mathematical trapezoidal reference theta=0.5.
  positiveLine->setTheta(0.55);
  negativeLine->setTheta(0.55);
  positiveLine->connect({remotePositive, localPositive});
  negativeLine->connect({localNegative, remoteNegative});

  auto localCapacitor = EMT::DC::SSN::Capacitor::make("localCapacitor");
  localCapacitor->setParameters(localCapacitance);
  localCapacitor->connect({localNegative, localPositive});
  // This single-station validation follows the validated Phase-1 boundary:
  // ideal remote pole sources represent the second DC system. Their symmetric
  // small reference steps exercise both differential-voltage directions
  // without introducing an unvalidated second converter.
  auto remotePositiveSource =
      EMT::DC::VoltageSource::make("remotePositiveSource");
  auto remoteNegativeSource =
      EMT::DC::VoltageSource::make("remoteNegativeSource");
  remotePositiveSource->setParameters(nominalDcVoltage / 2.0);
  remoteNegativeSource->setParameters(nominalDcVoltage / 2.0);
  remotePositiveSource->connect({SimNode<Real>::GND, remotePositive});
  remoteNegativeSource->connect({remoteNegative, SimNode<Real>::GND});

  auto mmc = EMT::Ph3::SSN_MMC::make("mmc");
  mmc->setParameters(frequency, acVoltage, nominalDcVoltage, 0.05, 1.07, 0.01,
                     400, 0.0005, 0.0001);
  mmc->setInitialAngle(0.0);
  mmc->setPLL(0.0, 0.0, false);
  mmc->setActiveControlOpenLoop(0.0);
  mmc->setReactiveControlOpenLoop(0.0);
  mmc->setEnergyController(120.0, 400.0, true);
  mmc->setOutputCurrentController(117.93, 8.5e4);
  mmc->setCirculatingCurrentController(19.93, 4500.0);
  mmc->setZeroSequenceCurrentController(19.93, 4500.0);
  mmc->setCirculatingCurrentReferences(0.0, 0.0, 0.0);
  mmc->setLimits(1000.0, 1000.0, 2.0);
  mmc->setOperatingPointInitialization(true, 50, 1e-8);
  mmc->connect({acNode, localPositive, localNegative});

  auto station = EMT::Ph3::SSN_MMCStation::make("station", mmc);
  EMT::Ph3::SSN_MMCStationParameters parameters;
  parameters.nominalPower = nominalPower;
  parameters.nominalAcLineLineRms = acVoltage;
  parameters.nominalDcVoltage = nominalDcVoltage;
  parameters.nominalFrequencyHz = frequency;
  parameters.controllerTimeStep = timeStep;
  parameters.armResistancePu = 0.0015;
  parameters.armInductancePu = 0.15;
  station->setParameters(parameters);
  station->setControlMode(
      EMT::Ph3::SSN_MMCStation::ControlMode::DCVoltageReactivePower);
  station->mAngle->set(0.0);
  station->mAngularFrequency->set(omega);

  SystemTopology system(
      frequency,
      SystemNodeList{acNode, localPositive, localNegative, remotePositive,
                     remoteNegative},
      SystemComponentList{acSource, positiveLine, negativeLine, localCapacitor,
                          remotePositiveSource, remoteNegativeSource, mmc,
                          station});
  Simulation sim("EMT_SSN_MMCStation_DCPiLine", Logger::Level::off);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.initialize();
  require(station->requestEnable(2e-5, 1.0),
          "PiLine station enable rejected; diagnostic=" +
              std::to_string(**station->mEnableDiagnostic));

  Real maxLocalPositiveKcl = 0.0;
  Real maxLocalNegativeKcl = 0.0;
  Real maxRemotePositiveKcl = 0.0;
  Real maxRemoteNegativeKcl = 0.0;
  Real minLocalVdc = std::numeric_limits<Real>::infinity();
  Real maxLocalVdc = -std::numeric_limits<Real>::infinity();
  Real minRemoteVdc = std::numeric_limits<Real>::infinity();
  Real maxRemoteVdc = -std::numeric_limits<Real>::infinity();
  Real minEnergy = std::numeric_limits<Real>::infinity();
  Real maxEnergy = -std::numeric_limits<Real>::infinity();
  Real finalLocalVdc = nominalDcVoltage;
  Real finalRemoteVdc = nominalDcVoltage;
  Real positiveQ = -std::numeric_limits<Real>::infinity();
  Real finalQ = 0.0;

  sim.start();
  while (sim.time() < sim.finalTime()) {
    station->mAngle->set(omega * sim.time());
    Real remotePoleVoltage = nominalDcVoltage / 2.0;
    if (sim.time() >= 0.03 && sim.time() < 0.05)
      remotePoleVoltage += 50.0;
    else if (sim.time() >= 0.08 && sim.time() < 0.10)
      remotePoleVoltage -= 50.0;
    remotePositiveSource->mVoltageRef->set(remotePoleVoltage);
    remoteNegativeSource->mVoltageRef->set(remotePoleVoltage);
    if (sim.time() >= 0.23 && sim.time() < 0.25)
      station->mReactivePowerReferencePu->set(0.01);
    else
      station->mReactivePowerReferencePu->set(0.0);
    sim.step();

    require(mmc->getState().allFinite(), "PiLine MMC state became non-finite.");
    const Real positiveCurrent = positiveLine->intfCurrent()(0, 0);
    const Real negativeCurrent = negativeLine->intfCurrent()(0, 0);
    const Real localCapCurrent = localCapacitor->intfCurrent()(0, 0);
    const Real remotePositiveSourceCurrent =
        remotePositiveSource->intfCurrent()(0, 0);
    const Real remoteNegativeSourceCurrent =
        remoteNegativeSource->intfCurrent()(0, 0);
    const Real localPositiveKcl =
        positiveCurrent + localCapCurrent + mmc->getInterfaceCurrent()(3, 0);
    const Real localNegativeKcl =
        -negativeCurrent - localCapCurrent + mmc->getInterfaceCurrent()(4, 0);
    const Real remotePositiveKcl =
        -positiveCurrent + remotePositiveSourceCurrent;
    const Real remoteNegativeKcl =
        negativeCurrent - remoteNegativeSourceCurrent;
    maxLocalPositiveKcl =
        std::max(maxLocalPositiveKcl, std::abs(localPositiveKcl));
    maxLocalNegativeKcl =
        std::max(maxLocalNegativeKcl, std::abs(localNegativeKcl));
    maxRemotePositiveKcl =
        std::max(maxRemotePositiveKcl, std::abs(remotePositiveKcl));
    maxRemoteNegativeKcl =
        std::max(maxRemoteNegativeKcl, std::abs(remoteNegativeKcl));

    const Real localVdc = **mmc->dcVoltageAttribute();
    const Real remoteVdc =
        remotePositive->voltage()(0, 0) - remoteNegative->voltage()(0, 0);
    minLocalVdc = std::min(minLocalVdc, localVdc);
    maxLocalVdc = std::max(maxLocalVdc, localVdc);
    minRemoteVdc = std::min(minRemoteVdc, remoteVdc);
    maxRemoteVdc = std::max(maxRemoteVdc, remoteVdc);
    const Real energy = **mmc->storedEnergyAttribute();
    minEnergy = std::min(minEnergy, energy);
    maxEnergy = std::max(maxEnergy, energy);
    if (sim.time() >= 0.23 && sim.time() < 0.25)
      positiveQ = std::max(positiveQ, **station->mReactivePowerPu);
    finalLocalVdc = localVdc;
    finalRemoteVdc = remoteVdc;
    finalQ = **station->mReactivePowerPu;
  }
  sim.stop();

  SPDLOG_INFO("Stage3F PiLine: local_Vdc=[{},{}] final={}, remote_Vdc=[{},{}] "
              "final={}, KCL=[{},{},{},{}] A, energy=[{},{}] J, Q_peak={} pu, "
              "Q_final={} pu, line_currents=[{},{}] A",
              minLocalVdc, maxLocalVdc, finalLocalVdc, minRemoteVdc,
              maxRemoteVdc, finalRemoteVdc, maxLocalPositiveKcl,
              maxLocalNegativeKcl, maxRemotePositiveKcl, maxRemoteNegativeKcl,
              minEnergy, maxEnergy, positiveQ, finalQ,
              positiveLine->intfCurrent()(0, 0),
              negativeLine->intfCurrent()(0, 0));

  require(maxLocalPositiveKcl < 1e-6 && maxLocalNegativeKcl < 1e-6 &&
              maxRemotePositiveKcl < 1e-6 && maxRemoteNegativeKcl < 1e-6,
          "PiLine KCL residual exceeds 1e-6 A.");
  require(maxRemoteVdc > nominalDcVoltage + 1.0 &&
              minRemoteVdc < nominalDcVoltage - 1.0,
          "PiLine disturbances did not exercise both Vdc directions.");
  require(std::abs(finalLocalVdc - nominalDcVoltage) < 200.0 &&
              std::abs(finalRemoteVdc - nominalDcVoltage) < 200.0,
          "PiLine DC voltage did not restore.");
  require(positiveQ > 0.001 && std::abs(finalQ) < 0.002,
          "PiLine reactive-power step or recovery failed.");
  require(minEnergy > 0.0 && maxEnergy / minEnergy < 1.1,
          "PiLine MMC energy is unbounded.");
  require(std::abs(positiveLine->intfCurrent()(0, 0) -
                   negativeLine->intfCurrent()(0, 0)) < 0.1,
          "PiLine pole-current continuity failed.");
  return 0;
}
