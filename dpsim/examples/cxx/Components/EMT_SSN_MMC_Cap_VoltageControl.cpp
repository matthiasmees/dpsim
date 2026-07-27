// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <cmath>
#include <complex>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_CurrentSource.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>

using namespace DPsim;
using namespace CPS;

namespace {

MatrixComp balancedVoltageReference(Real phaseVoltageAmplitude) {
  MatrixComp reference = MatrixComp::Zero(3, 1);

  // DPsim's three-phase EMT voltage source expects a line-to-line RMS phasor.
  // Harmony's V_m is the phase-to-neutral peak amplitude.
  const Real lineToLineRms = phaseVoltageAmplitude / RMS3PH_TO_PEAK1PH;
  const Complex phaseA(lineToLineRms, 0.0);

  reference(0, 0) = phaseA;
  reference(1, 0) = phaseA * SHIFT_TO_PHASE_B;
  reference(2, 0) = phaseA * SHIFT_TO_PHASE_C;
  return reference;
}

} // namespace

int main(int argc, char *argv[]) {
  const String simName = "EMT_SSN_MMC_DcVoltageController_CurrentPulse_Test";

  // -----------------------------------------------------------------------
  // Test objective
  // -----------------------------------------------------------------------
  //
  // Validate the complete DC-voltage-control path:
  //
  //   DC-current disturbance
  //       -> dynamic DC voltage
  //       -> DC-voltage PI
  //       -> iDeltaD reference
  //       -> output-current controller
  //       -> AC active power
  //       -> MMC energy controller
  //       -> iSigmaZ / DC current
  //       -> DC-voltage regulation.
  //
  // The AC terminal is connected directly to an ideal three-phase source.
  // There is deliberately no external AC RL branch in this test.
  //
  // Harmony implements Vdc as an internal state in DC-voltage-control mode.
  // The current DPsim MMC remains a V-type component, so this test represents
  // Harmony's effective capacitance as an external, correctly initialized
  // SSN capacitor at the MMC DC terminal.

  const Real timeStep = 50.0e-6;
  const Real finalTime = 1.0;

  // Positive current is injected into the DC bus from 0.10 s to 0.25 s.
  const Real currentPulseStart = 0.10;
  const Real currentPulseEnd = 0.25;
  const Real dcCurrentPulse = 20.0;

  // -----------------------------------------------------------------------
  // Converter and operating-point parameters
  // -----------------------------------------------------------------------
  const Real nominalFrequency = 50.0;
  const Real acVoltageAmplitude = 345.0e3;
  const Real nominalDcVoltage = 440.0e3;

  const Real armInductance = 0.05;
  const Real armResistance = 1.07;
  const Real submoduleCapacitance = 0.01;
  const UInt numberOfSubmodules = 400;
  const Real reactorInductance = 0.0005;
  const Real reactorResistance = 0.0001;

  // Harmony-equivalent DC capacitance:
  //
  //   Ce = 6 * Carm / N
  //
  // For Carm=0.01 F and N=400:
  //
  //   Ce = 150 uF.
  const Real effectiveDcCapacitance =
      6.0 * submoduleCapacitance / static_cast<Real>(numberOfSubmodules);

  // -----------------------------------------------------------------------
  // Controller parameters
  // -----------------------------------------------------------------------
  const Real kpOutputCurrent = 117.93;
  const Real kiOutputCurrent = 8.5e4;

  const Real kpCirculatingCurrent = 19.93;
  const Real kiCirculatingCurrent = 4500.0;

  const Real kpZeroSequenceCurrent = 19.93;
  const Real kiZeroSequenceCurrent = 4500.0;

  const Real kpEnergy = 120.0;
  const Real kiEnergy = 400.0;

  // Around the zero-power operating point:
  //
  //   Ceff * dVdc/dt =
  //       Idc_source - kId * iDeltaD
  //
  // with:
  //
  //   kId = 1.5 * Vac_peak / Vdc.
  //
  // The implemented outer controller uses:
  //
  //   eV = Vdc - Vdc_ref
  //   iDeltaD_ref = Kp * eV + Ki * integral(eV).
  //
  // Gains below place the idealized outer loop at 3 Hz with zeta=1.
  const Real dcVoltageLoopBandwidth = 2.0 * PI * 3.0;
  const Real dcVoltageLoopDamping = 1.0;
  const Real acCurrentToDcCurrentGain =
      1.5 * acVoltageAmplitude / nominalDcVoltage;

  const Real kpDcVoltage = 2.0 * dcVoltageLoopDamping * dcVoltageLoopBandwidth *
                           effectiveDcCapacitance / acCurrentToDcCurrentGain;

  const Real kiDcVoltage = dcVoltageLoopBandwidth * dcVoltageLoopBandwidth *
                           effectiveDcCapacitance / acCurrentToDcCurrentGain;

  const Real maximumAcCurrent = 100.0;
  const Real maximumCirculatingCurrent = 100.0;
  const Real maximumModulationMagnitude = 2.0;

  Logger::setLogDir("logs/" + simName);

  // -----------------------------------------------------------------------
  // Nodes
  // -----------------------------------------------------------------------
  auto nAc = SimNode<Real>::make("nAc", PhaseType::ABC);
  auto nDc = SimNode<Real>::make("nDc", PhaseType::DC);

  const Complex initialAcPhasor(acVoltageAmplitude / RMS3PH_TO_PEAK1PH, 0.0);

  nAc->setInitialVoltage(initialAcPhasor);
  nDc->setInitialVoltage(Complex(nominalDcVoltage, 0.0));

  // -----------------------------------------------------------------------
  // Ideal AC source
  // -----------------------------------------------------------------------
  auto acSource =
      EMT::Ph3::NetworkInjection::make("AcSource", Logger::Level::info);

  acSource->setParameters(balancedVoltageReference(acVoltageAmplitude),
                          nominalFrequency);

  acSource->connect({nAc});

  // -----------------------------------------------------------------------
  // Controlled DC current source
  // -----------------------------------------------------------------------
  //
  // EMT::DC::CurrentSource defines positive current from terminal 1 to
  // terminal 0. Connecting {nDc, GND} therefore makes a positive reference
  // inject current into the positive DC bus.
  auto dcCurrentSource =
      EMT::DC::CurrentSource::make("DcCurrentSource", Logger::Level::info);

  dcCurrentSource->setParameters(0.0);

  dcCurrentSource->connect({
      nDc,
      SimNode<Real>::GND,
  });

  auto dcCurrentReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(0.0));

  // -----------------------------------------------------------------------
  // Harmony-equivalent dynamic DC capacitance
  // -----------------------------------------------------------------------
  //
  auto effectiveDcCapacitor = EMT::DC::SSN::Capacitor::make(
      "EffectiveDcCapacitor", Logger::Level::info);

  effectiveDcCapacitor->setParameters(effectiveDcCapacitance);

  // The scalar DC capacitor uses voltage terminal1-terminal0. This gives
  // Vdc = V(nDc)-V(GND).
  effectiveDcCapacitor->connect({
      SimNode<Real>::GND,
      nDc,
  });

  // -----------------------------------------------------------------------
  // MMC
  // -----------------------------------------------------------------------
  auto mmc = EMT::Ph3::SSN_MMC::make("MMC1", "MMC1", Logger::Level::debug);

  mmc->setParameters(nominalFrequency, acVoltageAmplitude, nominalDcVoltage,
                     armInductance, armResistance, submoduleCapacitance,
                     numberOfSubmodules, reactorInductance, reactorResistance);

  // Keep the AC/control-frame test as simple as possible.
  mmc->setInitialAngle(0.0);
  mmc->setPLL(0.0, 0.0, false);
  mmc->setReactiveControlOpenLoop(0.0);

  mmc->setDcVoltageControl(nominalDcVoltage, kpDcVoltage, kiDcVoltage);

  // The energy controller is required in this external-capacitor realization.
  // It links the AC power command to the zero-sequence circulating current and
  // therefore to the current drawn from the external DC terminal.
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

  mmc->setEigenvalueDiagnostics(true, 200);

  mmc->setDiagnosticTimeStep(timeStep);

  mmc->connect({
      nAc,
      nDc,
      SimNode<Real>::GND,
  });

  auto dcVoltageReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(nominalDcVoltage));

  // -----------------------------------------------------------------------
  // System
  // -----------------------------------------------------------------------
  auto system = SystemTopology(nominalFrequency,
                               SystemNodeList{
                                   nAc,
                                   nDc,
                               },
                               SystemComponentList{
                                   acSource,
                                   dcCurrentSource,
                                   effectiveDcCapacitor,
                                   mmc,
                               });

  // -----------------------------------------------------------------------
  // Logging
  // -----------------------------------------------------------------------
  auto logger = DataLogger::make(simName);

  logger->logAttribute("Voltage_AC", nAc->attribute("v"));

  logger->logAttribute("Voltage_DC", nDc->attribute("v"));

  logger->logAttribute("DcCurrentReference", dcCurrentReferenceSignal);

  logger->logAttribute("Current_DC_Source",
                       dcCurrentSource->attribute("i_intf"));

  logger->logAttribute("Current_DC_Capacitor",
                       effectiveDcCapacitor->attribute("i_intf"));

  logger->logAttribute("VdcReference", dcVoltageReferenceSignal);

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

  Bool applyCurrentPulse = true;
  Bool removeCurrentPulse = true;

  auto setDcCurrent = [&](Real current, const char *description) {
    SPDLOG_INFO("{}: setting injected DC current to {:.3f} A "
                "at t={:.6f} s.",
                description, current, sim.time());

    dcCurrentSource->setParameters(current);

    dcCurrentReferenceSignal->set(current);
  };

  sim.initialize();

  sim.start();

  while (sim.time() < sim.finalTime()) {
    if (applyCurrentPulse && sim.time() >= currentPulseStart) {
      setDcCurrent(dcCurrentPulse, "Applying positive DC-current pulse");

      applyCurrentPulse = false;
    }

    if (removeCurrentPulse && sim.time() >= currentPulseEnd) {
      setDcCurrent(0.0, "Removing DC-current pulse");

      removeCurrentPulse = false;
    }

    sim.step();
  }

  sim.stop();
  return 0;
}
