// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <vector>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_Ph1_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph1_CurrentSource.h>
#include <dpsim-models/EMT/EMT_Ph1_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph1_SSNTypeI2T.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_Full_Serial_RLC.h>
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
  const String simName = "EMT_SSN_MMC_PiLine";

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
  // The AC terminal remains connected through the already validated
  // SSN series RL branch. The new element in this test is a single-pole
  // EMT::Ph1::PiLine on the DC side. The current disturbance is applied at
  // the remote end of that DC line.
  //
  // Harmony implements Vdc as an internal state in DC-voltage-control mode.
  // The current DPsim MMC remains a V-type component, so this test represents
  // Harmony's effective capacitance as an external, correctly initialized
  // SSN capacitor at the MMC DC terminal.

  const Real timeStep = 50.0e-6;
  const Real finalTime = 0.70;

  // Positive current is injected into the DC bus from 0.10 s to 0.25 s.
  const Real currentPulseStart = 0.10;
  const Real currentPulseEnd = 0.25;
  const Real dcCurrentPulse = 20.0;

  // -----------------------------------------------------------------------
  // Converter and operating-point parameters
  // -----------------------------------------------------------------------
  const Real nominalFrequency = 50.0;
  const Real nominalOmega = 2.0 * PI * nominalFrequency;
  const Real acVoltageAmplitude = 345.0e3;
  const Real nominalDcVoltage = 440.0e3;

  // Harmony source plus AC branch:
  //
  //   Zsource + Zbranch = 5 ohm + (5 + j140) ohm
  //                     = 10 + j140 ohm.
  const Real acSeriesResistance = 10.0;
  const Real acSeriesReactance = 140.0;
  const Real acSeriesInductance = acSeriesReactance / nominalOmega;

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
  // DC PiLine test-bench parameters
  // -----------------------------------------------------------------------
  //
  // These are deliberately moderate first-test values rather than claimed
  // Harmony data. They introduce a remote DC node, cable current dynamics and
  // distributed shunt capacitance without changing the proven converter
  // tuning.
  const Real dcLineSeriesResistance = 5.0;        // ohm
  const Real dcLineSeriesInductance = 50.0e-3;    // H
  const Real dcLineParallelCapacitance = 20.0e-6; // F total
  const Real dcLineParallelConductance = 1.0e-8;  // S total

  // -----------------------------------------------------------------------
  // Controller parameters
  // -----------------------------------------------------------------------
  const Real kpPll = 0.001103374;
  const Real kiPll = 0.00073;

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
  //   eV = Vdc_ref - Vdc
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
  auto nAcSource = SimNode<Real>::make("nAcSource", PhaseType::ABC);

  auto nMmcAc = SimNode<Real>::make("nMmcAc", PhaseType::ABC);

  auto nDcMmc = SimNode<Real>::make("nDcMmc", PhaseType::Single);

  auto nDcRemote = SimNode<Real>::make("nDcRemote", PhaseType::Single);

  const Complex initialAcPhasor(acVoltageAmplitude / RMS3PH_TO_PEAK1PH, 0.0);

  // The initial operating point is zero power, so the SSN RL current is zero
  // and both AC nodes start at the same phasor.
  nAcSource->setInitialVoltage(initialAcPhasor);
  nMmcAc->setInitialVoltage(initialAcPhasor);
  nDcMmc->setInitialVoltage(Complex(nominalDcVoltage, 0.0));
  nDcRemote->setInitialVoltage(Complex(nominalDcVoltage, 0.0));

  // -----------------------------------------------------------------------
  // Ideal AC source
  // -----------------------------------------------------------------------
  auto acSource =
      EMT::Ph3::NetworkInjection::make("AcSource", Logger::Level::info);

  acSource->setParameters(balancedVoltageReference(acVoltageAmplitude),
                          nominalFrequency);

  acSource->connect({nAcSource});

  // -----------------------------------------------------------------------
  // SSN source-and-branch series RL
  // -----------------------------------------------------------------------
  auto acSeriesRl =
      EMT::Ph3::SSN::Full_Serial_RLC::make("AcSeriesRL", Logger::Level::info);

  acSeriesRl->setIntegrationTheta(0.55);

  acSeriesRl->setParameters(
      Math::singlePhaseParameterToThreePhase(acSeriesResistance),
      Math::singlePhaseParameterToThreePhase(acSeriesInductance),
      Matrix::Zero(3, 3));

  // TwoTerminalVTypeSSNComp defines:
  //
  //   u = v_terminal1 - v_terminal0.
  //
  // This orientation gives u = Vsource - Vmmc and positive branch current
  // from the source toward the MMC.
  acSeriesRl->connect({
      nMmcAc,
      nAcSource,
  });

  // -----------------------------------------------------------------------
  // Single-pole DC PiLine
  // -----------------------------------------------------------------------
  auto dcLine = EMT::Ph1::PiLine::make("DcPiLine", Logger::Level::info);

  dcLine->setParameters(dcLineSeriesResistance, dcLineSeriesInductance,
                        dcLineParallelCapacitance, dcLineParallelConductance);

  dcLine->connect({
      nDcMmc,
      nDcRemote,
  });

  // Construct the composite subcomponents now so that the two internal shunt
  // capacitors can be logged and corrected after MNA initialization.
  dcLine->createSubComponents();

  std::vector<EMT::Ph1::Capacitor::Ptr> dcLineShuntCapacitors;
  for (const auto &subComponent : dcLine->mnaSubComponents()) {
    auto capacitor =
        std::dynamic_pointer_cast<EMT::Ph1::Capacitor>(subComponent);
    if (capacitor)
      dcLineShuntCapacitors.push_back(capacitor);
  }

  if (dcLineShuntCapacitors.size() != 2)
    throw std::runtime_error(
        "Expected exactly two DC PiLine shunt capacitors.");

  auto dcLineCapMmc = dcLineShuntCapacitors[0];
  auto dcLineCapRemote = dcLineShuntCapacitors[1];

  // -----------------------------------------------------------------------
  // Controlled remote-end DC current source
  // -----------------------------------------------------------------------
  //
  // EMT::Ph1::CurrentSource defines positive current from terminal 0 to
  // terminal 1. Connecting {GND, nDcRemote} therefore injects a positive
  // current into the remote DC-line node.
  auto dcCurrentSource = EMT::Ph1::CurrentSource::make("DcRemoteCurrentSource",
                                                       Logger::Level::info);

  dcCurrentSource->setParameters(Complex(0.0, 0.0), -1.0);

  dcCurrentSource->connect({
      SimNode<Real>::GND,
      nDcRemote,
  });

  auto dcCurrentReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(0.0));

  // -----------------------------------------------------------------------
  // Harmony-equivalent dynamic DC capacitance
  // -----------------------------------------------------------------------
  //
  // Use SSNTypeI2T instead of EMT::Ph1::Capacitor to avoid interpreting
  // 440 kV as an AC RMS phasor during initialization.
  //
  // State-space model:
  //
  //   x = Vdc
  //   dx/dt = iC / Ce
  //   y = Vdc.
  Matrix capacitorA = Matrix::Zero(1, 1);
  Matrix capacitorB = Matrix::Zero(1, 1);
  Matrix capacitorC = Matrix::Zero(1, 1);
  Matrix capacitorD = Matrix::Zero(1, 1);

  capacitorB(0, 0) = 1.0 / effectiveDcCapacitance;
  capacitorC(0, 0) = 1.0;

  auto effectiveDcCapacitor =
      EMT::Ph1::SSNTypeI2T::make("EffectiveDcCapacitor", Logger::Level::info);

  effectiveDcCapacitor->setParameters(capacitorA, capacitorB, capacitorC,
                                      capacitorD);

  // SSNTypeI2T uses the voltage terminal1-terminal0. This connection gives
  // Vdc = V(nDcMmc)-V(GND).
  effectiveDcCapacitor->connect({
      SimNode<Real>::GND,
      nDcMmc,
  });

  // -----------------------------------------------------------------------
  // MMC
  // -----------------------------------------------------------------------
  auto mmc = EMT::Ph3::SSN_MMC::make("MMC1", "MMC1", Logger::Level::debug);

  mmc->setParameters(nominalFrequency, acVoltageAmplitude, nominalDcVoltage,
                     armInductance, armResistance, submoduleCapacitance,
                     numberOfSubmodules, reactorInductance, reactorResistance);

  // Keep reactive control open-loop, but retain the already validated PLL.
  mmc->setInitialAngle(0.0);
  mmc->setPLL(kpPll, kiPll, true);
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
      nMmcAc,
      nDcMmc,
      SimNode<Real>::GND,
  });

  auto dcVoltageReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(nominalDcVoltage));

  // -----------------------------------------------------------------------
  // System
  // -----------------------------------------------------------------------
  auto system = SystemTopology(nominalFrequency,
                               SystemNodeList{
                                   nAcSource,
                                   nMmcAc,
                                   nDcMmc,
                                   nDcRemote,
                               },
                               SystemComponentList{
                                   acSource,
                                   acSeriesRl,
                                   dcLine,
                                   dcCurrentSource,
                                   effectiveDcCapacitor,
                                   mmc,
                               });

  // -----------------------------------------------------------------------
  // Logging
  // -----------------------------------------------------------------------
  auto logger = DataLogger::make(simName);

  logger->logAttribute("Voltage_AC_Source", nAcSource->attribute("v"));

  logger->logAttribute("Voltage_AC_MMC", nMmcAc->attribute("v"));

  logger->logAttribute("Current_AC_SeriesRL", acSeriesRl->attribute("i_intf"));

  logger->logAttribute("MMC_InterfaceCurrent", mmc->attribute("i_intf"));

  logger->logAttribute("Voltage_DC_MMC", nDcMmc->attribute("v"));

  logger->logAttribute("Voltage_DC_Remote", nDcRemote->attribute("v"));

  logger->logAttribute("Voltage_DC_Line", dcLine->attribute("v_intf"));

  logger->logAttribute("Current_DC_Line", dcLine->attribute("i_intf"));

  logger->logAttribute("DcCurrentReference", dcCurrentReferenceSignal);

  logger->logAttribute("Current_DC_RemoteSource",
                       dcCurrentSource->attribute("i_intf"));

  logger->logAttribute("Current_DC_EffectiveCapacitor",
                       effectiveDcCapacitor->attribute("i_intf"));

  logger->logAttribute("Current_DC_LineCap_MMC",
                       dcLineCapMmc->attribute("i_intf"));

  logger->logAttribute("Current_DC_LineCap_Remote",
                       dcLineCapRemote->attribute("i_intf"));

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

  logger->logAttribute("MMC_PLLFrequency", mmc->attribute("pll_frequency"));

  logger->logAttribute("MMC_PLLError", mmc->attribute("pll_error"));

  logger->logAttribute("MMC_PLLAngleDeviation",
                       mmc->attribute("pll_angle_deviation"));

  logger->logAttribute("MMC_ControlVoltageD", mmc->attribute("v_control_d"));

  logger->logAttribute("MMC_ControlVoltageQ", mmc->attribute("v_control_q"));

  logger->logAttribute("MMC_IDeltaQReference", mmc->attribute("i_delta_q_ref"));

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

    dcCurrentSource->setParameters(Complex(current, 0.0), -1.0);

    dcCurrentReferenceSignal->set(current);
  };

  sim.initialize();

  // SSNTypeI2T initializes its state to zero. Restore the charged converter
  // DC-bus operating point after initialization and before the first step.
  effectiveDcCapacitor->manualInit(Matrix::Constant(1, 1, nominalDcVoltage),
                                   Matrix::Zero(1, 1), Matrix::Zero(1, 1), 0.0,
                                   nominalDcVoltage);

  // EMT::Ph1::PiLine internally uses EMT::Ph1::Capacitor, whose generic
  // initialization interprets node values as AC RMS phasors and multiplies
  // them by RMS3PH_TO_PEAK1PH. For a scalar DC line this would initialize
  // 440 kV as 359.258 kV. Correct the two private shunt-capacitor states here.
  for (const auto &capacitor : dcLineShuntCapacitors) {
    capacitor->setIntfVoltage(Matrix::Constant(1, 1, nominalDcVoltage));
    capacitor->setIntfCurrent(Matrix::Zero(1, 1));
  }

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
