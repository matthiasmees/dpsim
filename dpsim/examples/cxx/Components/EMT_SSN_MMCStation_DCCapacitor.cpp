// SPDX-License-Identifier: MPL-2.0
#include "../Examples.h"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_CurrentSource.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Capacitor.h>
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
  const Real dcCapacitance = 1.5e-3;

  auto acNode = SimNode<Real>::make("ac", PhaseType::ABC);
  auto dcPositive = SimNode<Real>::make("dcPositive", PhaseType::DC);
  acNode->setInitialVoltage(Complex(acVoltage / RMS3PH_TO_PEAK1PH, 0.0));
  dcPositive->setInitialVoltage(Complex(nominalDcVoltage, 0.0));

  auto acSource = EMT::Ph3::NetworkInjection::make("acSource");
  acSource->setParameters(balancedVoltageReference(acVoltage), frequency);
  acSource->connect({acNode});

  auto dcCapacitor = EMT::DC::SSN::Capacitor::make("dcCapacitor");
  dcCapacitor->setParameters(dcCapacitance);
  dcCapacitor->connect({SimNode<Real>::GND, dcPositive});

  auto disturbance = EMT::DC::CurrentSource::make("disturbance");
  disturbance->setParameters(0.0);
  disturbance->connect({dcPositive, SimNode<Real>::GND});

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
  mmc->connect({acNode, dcPositive, SimNode<Real>::GND});

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
      frequency, SystemNodeList{acNode, dcPositive},
      SystemComponentList{acSource, dcCapacitor, disturbance, mmc, station});
  Simulation sim("EMT_SSN_MMCStation_DCCapacitor", Logger::Level::off);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.initialize();
  require(station->requestEnable(2e-5, 1.0),
          "DC-capacitor station enable rejected; diagnostic=" +
              std::to_string(**station->mEnableDiagnostic));

  Real maximumKcl = 0.0;
  Real minimumVdc = std::numeric_limits<Real>::infinity();
  Real maximumVdc = -std::numeric_limits<Real>::infinity();
  Real minimumEnergy = std::numeric_limits<Real>::infinity();
  Real maximumEnergy = -std::numeric_limits<Real>::infinity();
  Real positiveDeviation = 0.0;
  Real negativeDeviation = 0.0;
  Real positiveQ = -std::numeric_limits<Real>::infinity();
  Real finalVdc = nominalDcVoltage;
  Real finalQ = 0.0;

  sim.start();
  while (sim.time() < sim.finalTime()) {
    station->mAngle->set(omega * sim.time());
    if (sim.time() >= 0.03 && sim.time() < 0.05)
      disturbance->mCurrentRef->set(10.0);
    else if (sim.time() >= 0.07 && sim.time() < 0.09)
      disturbance->mCurrentRef->set(-10.0);
    else
      disturbance->mCurrentRef->set(0.0);
    if (sim.time() >= 0.23 && sim.time() < 0.25)
      station->mReactivePowerReferencePu->set(0.01);
    else
      station->mReactivePowerReferencePu->set(0.0);

    sim.step();
    require(mmc->getState().allFinite(), "MMC state became non-finite.");
    const Real kcl = mmc->getInterfaceCurrent()(3, 0) -
                     disturbance->intfCurrent()(0, 0) +
                     dcCapacitor->intfCurrent()(0, 0);
    maximumKcl = std::max(maximumKcl, std::abs(kcl));
    const Real vdc = **mmc->dcVoltageAttribute();
    minimumVdc = std::min(minimumVdc, vdc);
    maximumVdc = std::max(maximumVdc, vdc);
    if (sim.time() >= 0.03 && sim.time() < 0.07)
      positiveDeviation = std::max(positiveDeviation, vdc - nominalDcVoltage);
    if (sim.time() >= 0.07 && sim.time() < 0.11)
      negativeDeviation = std::min(negativeDeviation, vdc - nominalDcVoltage);
    if (sim.time() >= 0.23 && sim.time() < 0.25)
      positiveQ = std::max(positiveQ, **station->mReactivePowerPu);
    const Real energy = **mmc->storedEnergyAttribute();
    minimumEnergy = std::min(minimumEnergy, energy);
    maximumEnergy = std::max(maximumEnergy, energy);
    finalVdc = vdc;
    finalQ = **station->mReactivePowerPu;
  }
  sim.stop();

  SPDLOG_INFO(
      "Stage3F capacitor: Vdc=[{},{}] V, deviations=[{},{}] V, final={} V, "
      "max_KCL={} A, energy=[{},{}] J, positive_Q={} pu, final_Q={} pu, "
      "final_id_ref={} pu, final_id={} pu, final_P={} pu, final_idc={} A",
      minimumVdc, maximumVdc, negativeDeviation, positiveDeviation, finalVdc,
      maximumKcl, minimumEnergy, maximumEnergy, positiveQ, finalQ,
      **station->mIdReferencePu, **station->mIdPu, **station->mActivePowerPu,
      **mmc->dcCurrentAttribute());

  require(maximumKcl < 1e-6, "DC-capacitor KCL residual exceeds 1e-6 A.");
  require(positiveDeviation > 1.0,
          "Positive injected current did not raise DC voltage.");
  require(negativeDeviation < -1.0,
          "Negative injected current did not lower DC voltage.");
  require(std::abs(finalVdc - nominalDcVoltage) < 200.0,
          "DC-voltage loop did not restore the capacitor voltage.");
  require(positiveQ > 0.001,
          "Positive Q reference did not produce positive reactive power.");
  require(std::abs(finalQ) < 0.002,
          "Reactive power did not recover after the capacitor test.");
  require(minimumEnergy > 0.0 && maximumEnergy / minimumEnergy < 1.1,
          "MMC energy is unbounded in the DC-capacitor test.");
  return 0;
}
