// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <cmath>
#include <complex>
#include <stdexcept>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Resistor.h>
#include <dpsim-models/EMT/EMT_DC_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>

using namespace DPsim;
using namespace CPS;

namespace {

MatrixComp balancedVoltageReference(Real phaseVoltageAmplitude) {
  MatrixComp reference = MatrixComp::Zero(3, 1);

  const Real lineToLineRms = phaseVoltageAmplitude / RMS3PH_TO_PEAK1PH;
  const Complex phaseA(lineToLineRms, 0.0);

  reference(0, 0) = phaseA;
  reference(1, 0) = phaseA * SHIFT_TO_PHASE_B;
  reference(2, 0) = phaseA * SHIFT_TO_PHASE_C;
  return reference;
}

void requireFinite(const Matrix &value, const String &name) {
  if (!value.allFinite())
    throw std::runtime_error(name + " contains NaN or Inf.");
}

} // namespace

int main(int argc, char *argv[]) {
  const String simName = "EMT_SSN_MMC_ActivePowerController_Isolation_Test";

  // -----------------------------------------------------------------------
  // Purpose
  // -----------------------------------------------------------------------
  //
  // Isolate the only controller that has not yet been validated independently:
  // the active-power controller used by the second MMC in the P2P case.
  //
  // Deliberately excluded:
  // - DC PiLine,
  // - external DC capacitor,
  // - second MMC,
  // - AC series RL,
  // - PLL,
  // - reactive-power controller.
  //
  // Retained from the Harmony controller stack:
  // - active-power PI,
  // - energy PI,
  // - output-current controller,
  // - circulating-current controller,
  // - zero-sequence-current controller.
  //
  // Both AC and DC terminal voltages are stiff. Therefore, divergence in this
  // test is inside the active-power/energy/ZCC/OCC MMC model and cannot be
  // attributed to the DC line or the voltage-controlled converter.

  const Real timeStep = 50.0e-6;
  const Real finalTime = 0.80;

  const Real positiveStepTime = 0.10;
  const Real negativeStepTime = 0.35;
  const Real zeroRestoreTime = 0.60;

  const Real positivePowerReference = 10.0e6;
  const Real negativePowerReference = -10.0e6;

  const Real nominalFrequency = 50.0;
  const Real acVoltageAmplitude = 345.0e3;
  const Real nominalDcVoltage = 440.0e3;
  const Real dcLoadResistance = 44.0e3;

  const Real armInductance = 0.05;
  const Real armResistance = 1.07;
  const Real submoduleCapacitance = 0.01;
  const UInt numberOfSubmodules = 400;
  const Real reactorInductance = 0.0005;
  const Real reactorResistance = 0.0001;

  // Harmony controller gains used by the active-power-controlled converter.
  const Real kpActivePower = 6.6667e-7;
  const Real kiActivePower = 3.3333e-4;

  const Real kpEnergy = 120.0;
  const Real kiEnergy = 400.0;

  const Real kpOutputCurrent = 117.93;
  const Real kiOutputCurrent = 8.5e4;

  const Real kpCirculatingCurrent = 19.93;
  const Real kiCirculatingCurrent = 4500.0;

  const Real kpZeroSequenceCurrent = 19.93;
  const Real kiZeroSequenceCurrent = 4500.0;

  const Real maximumAcCurrent = 100.0;
  const Real maximumCirculatingCurrent = 100.0;
  const Real maximumModulationMagnitude = 2.0;

  Logger::setLogDir("logs/" + simName);

  // -----------------------------------------------------------------------
  // Nodes and stiff terminal voltages
  // -----------------------------------------------------------------------
  auto nAc = SimNode<Real>::make("nAc", PhaseType::ABC);

  auto nDc = SimNode<Real>::make("nDc", PhaseType::DC);

  const Complex initialAcPhasor(acVoltageAmplitude / RMS3PH_TO_PEAK1PH, 0.0);

  nAc->setInitialVoltage(initialAcPhasor);
  nDc->setInitialVoltage(Complex(nominalDcVoltage, 0.0));

  auto acSource =
      EMT::Ph3::NetworkInjection::make("AcSource", Logger::Level::info);

  acSource->setParameters(balancedVoltageReference(acVoltageAmplitude),
                          nominalFrequency);

  acSource->connect({nAc});

  auto dcSource = EMT::DC::VoltageSource::make("DcSource", Logger::Level::info);

  dcSource->setParameters(nominalDcVoltage);

  // Positive DC voltage: V(nDc)-V(GND) = 440 kV.
  dcSource->connect({
      SimNode<Real>::GND,
      nDc,
  });

  // Scalar-DC resistor regression. With the stiff 440 kV bus, the passive
  // terminal convention v=v1-v0 and i:1->0 gives exactly +10 A and +4.4 MW.
  auto dcLoad = EMT::DC::SSN::Resistor::make("DcLoad", Logger::Level::info);
  dcLoad->setParameters(dcLoadResistance);
  dcLoad->connect({SimNode<Real>::GND, nDc});

  // -----------------------------------------------------------------------
  // Active-power-controlled MMC
  // -----------------------------------------------------------------------
  auto mmc = EMT::Ph3::SSN_MMC::make("MMC", "MMC", Logger::Level::debug);

  mmc->setParameters(nominalFrequency, acVoltageAmplitude, nominalDcVoltage,
                     armInductance, armResistance, submoduleCapacitance,
                     numberOfSubmodules, reactorInductance, reactorResistance);

  // Keep frames and AC terminal dynamics out of this first isolation test.
  mmc->setInitialAngle(0.0);
  mmc->setPLL(0.0, 0.0, false);
  mmc->setReactiveControlOpenLoop(0.0);

  mmc->setActivePowerControl(0.0, kpActivePower, kiActivePower);

  mmc->setEnergyController(kpEnergy, kiEnergy, true);

  mmc->setOutputCurrentController(kpOutputCurrent, kiOutputCurrent);

  mmc->setCirculatingCurrentController(kpCirculatingCurrent,
                                       kiCirculatingCurrent);

  mmc->setZeroSequenceCurrentController(kpZeroSequenceCurrent,
                                        kiZeroSequenceCurrent);

  mmc->setCirculatingCurrentReferences(0.0, 0.0, 0.0);

  mmc->setLimits(maximumAcCurrent, maximumCirculatingCurrent,
                 maximumModulationMagnitude);

  mmc->setOperatingPointInitialization(true, 50, 1.0e-8);

  mmc->setEigenvalueDiagnostics(true, 100);

  mmc->setDiagnosticTimeStep(timeStep);

  mmc->connect({
      nAc,
      nDc,
      SimNode<Real>::GND,
  });

  auto powerReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(0.0));

  // -----------------------------------------------------------------------
  // System and logging
  // -----------------------------------------------------------------------
  auto system = SystemTopology(nominalFrequency,
                               SystemNodeList{
                                   nAc,
                                   nDc,
                               },
                               SystemComponentList{
                                   acSource,
                                   dcSource,
                                   dcLoad,
                                   mmc,
                               });

  auto logger = DataLogger::make(simName);

  logger->logAttribute("Voltage_AC", nAc->attribute("v"));

  logger->logAttribute("Voltage_DC", nDc->attribute("v"));
  logger->logAttribute("Current_DC_Resistor", dcLoad->attribute("i_intf"));

  logger->logAttribute("ActivePowerReference", powerReferenceSignal);

  logger->logAttribute("MMC_InterfaceCurrent", mmc->attribute("i_intf"));

  logger->logAttribute("MMC_Vdc", mmc->attribute("vdc"));

  logger->logAttribute("MMC_Idc", mmc->attribute("idc"));

  logger->logAttribute("MMC_Pac", mmc->attribute("p_ac"));

  logger->logAttribute("MMC_Qac", mmc->attribute("q_ac"));

  logger->logAttribute("MMC_Pdc", mmc->attribute("p_dc"));

  logger->logAttribute("MMC_PowerBalanceError",
                       mmc->attribute("power_balance_error"));

  logger->logAttribute("MMC_IDeltaD", mmc->attribute("i_delta_d"));

  logger->logAttribute("MMC_IDeltaDReference", mmc->attribute("i_delta_d_ref"));

  logger->logAttribute("MMC_IDeltaQ", mmc->attribute("i_delta_q"));

  logger->logAttribute("MMC_ISigmaZ", mmc->attribute("i_sigma_z"));

  logger->logAttribute("MMC_ISigmaZReference", mmc->attribute("i_sigma_z_ref"));

  logger->logAttribute("MMC_StoredEnergy", mmc->attribute("stored_energy"));

  logger->logAttribute("MMC_StateNorm", mmc->attribute("state_norm"));

  logger->logAttribute("MMC_StateDerivativeNorm",
                       mmc->attribute("state_derivative_norm"));

  logger->logAttribute("MMC_JacobianMaxRealEigenvalue",
                       mmc->attribute("jacobian_max_real_eigenvalue"));

  logger->logAttribute("MMC_JacobianMaxDiscreteMagnitude",
                       mmc->attribute("jacobian_max_discrete_magnitude"));

  logger->logAttribute("MMC_DiagnosticsValid",
                       mmc->attribute("diagnostics_valid"));

  // -----------------------------------------------------------------------
  // Simulation
  // -----------------------------------------------------------------------
  Simulation sim(simName, Logger::Level::info);

  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);

  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.addLogger(logger);

  Bool applyPositiveStep = true;
  Bool applyNegativeStep = true;
  Bool restoreZero = true;

  auto setPowerReference = [&](Real reference, const char *description) {
    SPDLOG_INFO("{}: setting active-power reference to {:.3f} MW "
                "at t={:.6f} s.",
                description, reference / 1.0e6, sim.time());

    mmc->setActivePowerControl(reference, kpActivePower, kiActivePower);

    powerReferenceSignal->set(reference);
  };

  sim.initialize();
  sim.start();

  while (sim.time() < sim.finalTime()) {
    if (applyPositiveStep && sim.time() >= positiveStepTime) {
      setPowerReference(positivePowerReference,
                        "Applying positive active-power step");
      applyPositiveStep = false;
    }

    if (applyNegativeStep && sim.time() >= negativeStepTime) {
      setPowerReference(negativePowerReference,
                        "Applying negative active-power step");
      applyNegativeStep = false;
    }

    if (restoreZero && sim.time() >= zeroRestoreTime) {
      setPowerReference(0.0, "Restoring zero active power");
      restoreZero = false;
    }

    sim.step();
  }

  sim.stop();

  requireFinite(mmc->getState(), "MMC state");
  requireFinite(mmc->getInterfaceVoltage(), "MMC interface voltage");
  requireFinite(mmc->getInterfaceCurrent(), "MMC interface current");

  const Real loadCurrent = dcLoad->intfCurrent()(0, 0);
  const Real dcNodeKclResidual = loadCurrent + dcSource->intfCurrent()(0, 0) +
                                 mmc->getInterfaceCurrent()(3, 0);
  const Real finalPac = mmc->attributeTyped<Real>("p_ac")->get();
  const Real finalPdc = mmc->attributeTyped<Real>("p_dc")->get();

  if (!std::isfinite(dcNodeKclResidual) || std::abs(dcNodeKclResidual) > 1.0e-6)
    throw std::runtime_error(
        fmt::format("scalar DC-node KCL residual is {} A", dcNodeKclResidual));
  if (std::abs(loadCurrent - nominalDcVoltage / dcLoadResistance) > 1.0e-9)
    throw std::runtime_error(
        fmt::format("scalar DC resistor current is {} A, expected {} A",
                    loadCurrent, nominalDcVoltage / dcLoadResistance));

  SPDLOG_INFO("Scalar-DC MMC resistor validation: Vdc={} V, Iload={} A, "
              "KCL residual={} A, Pac={} W, Pdc={} W, Pac-Pdc={} W",
              nDc->voltage()(0, 0), loadCurrent, dcNodeKclResidual, finalPac,
              finalPdc, finalPac - finalPdc);

  return 0;
}
