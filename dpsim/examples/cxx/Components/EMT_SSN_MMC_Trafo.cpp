// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <cmath>
#include <complex>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_Ph1_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>

using namespace DPsim;
using namespace CPS;

namespace {

MatrixComp balancedVoltageReferenceFromLineToLineRms(Real lineToLineRms) {
  MatrixComp reference = MatrixComp::Zero(3, 1);

  const Complex phaseA(lineToLineRms, 0.0);
  reference(0, 0) = phaseA;
  reference(1, 0) = phaseA * SHIFT_TO_PHASE_B;
  reference(2, 0) = phaseA * SHIFT_TO_PHASE_C;
  return reference;
}

Complex initialNodePhasorFromLineToLineRms(Real lineToLineRms,
                                           Real angle = 0.0) {
  return std::polar(lineToLineRms, angle);
}

} // namespace

int main(int argc, char *argv[]) {
  const String simName = "EMT_SSN_MMC_DQsym_Trafo";

  // -----------------------------------------------------------------------
  // Scope of this test
  // -----------------------------------------------------------------------
  //
  // This test validates only the new AC topology taken from the DQsym/SPS
  // P2P case:
  //
  //   ideal 400-kV grid source
  //       -> grid series RL
  //       -> 400/333-kV transformer
  //       -> SSN MMC
  //       -> stiff bipolar +/-320-kV DC source
  //
  // It intentionally does NOT reproduce the DQsym start-up sequence yet.
  // The present SSN_MMC has no blocked-diode precharge mode, so the converter
  // is initialized in its already-energized operating state.
  //
  // Once this passes, the same grid/transformer assembly can be duplicated
  // for the second station and connected through the bipolar DC cable.

  const Real timeStep = 20.0e-6;
  const Real finalTime = 0.40;

  const Real powerStepTime = 0.10;
  const Real powerRestoreTime = 0.30;
  const Real activePowerStep = 50.0e6;

  // -----------------------------------------------------------------------
  // DQsym P2P base data
  // -----------------------------------------------------------------------
  const Real nominalFrequency = 50.0;
  const Real nominalOmega = 2.0 * PI * nominalFrequency;

  const Real converterRatedPower = 1000.0e6;
  const Real primaryVoltageRmsLl = 400.0e3;
  const Real secondaryVoltageRmsLl = 333.0e3;
  const Real nominalDcVoltage = 640.0e3;

  // Power-invariant d-axis / phase-peak voltage used by SSN_MMC.
  const Real converterAcVoltageAmplitude =
      RMS3PH_TO_PEAK1PH * secondaryVoltageRmsLl;

  // DQsym converter hardware.
  const UInt numberOfSubmodules = 36;
  const Real submoduleCapacitance = 1.758e-3;

  const Real armInductancePu = 0.15;
  const Real armResistancePu = armInductancePu / 100.0;

  const Real converterAcImpedanceBase =
      secondaryVoltageRmsLl * secondaryVoltageRmsLl / converterRatedPower;

  const Real armInductance =
      armInductancePu * converterAcImpedanceBase / nominalOmega;

  const Real armResistance = armResistancePu * converterAcImpedanceBase;

  // The transformer is modelled explicitly, so no additional external
  // converter reactor is included in the MMC parameters.
  const Real reactorInductance = 0.0;
  const Real reactorResistance = 0.0;

  // -----------------------------------------------------------------------
  // DQsym transformer data
  // -----------------------------------------------------------------------
  //
  // SPS uses total winding values:
  //   Rxfo = 0.003 pu
  //   Lxfo = 0.12 pu
  //
  // EMT::Ph3::Transformer places its series R/L on the higher-voltage side,
  // therefore the values are referred to the 400-kV side.
  const Real transformerResistancePu = 0.003;
  const Real transformerReactancePu = 0.12;

  const Real transformerPrimaryBaseImpedance =
      primaryVoltageRmsLl * primaryVoltageRmsLl / converterRatedPower;

  const Real transformerResistance =
      transformerResistancePu * transformerPrimaryBaseImpedance;

  const Real transformerInductance =
      transformerReactancePu * transformerPrimaryBaseImpedance / nominalOmega;

  const Real transformerRatio = primaryVoltageRmsLl / secondaryVoltageRmsLl;

  // Explicit RL grid branch from the uploaded SPS model's right-hand grid.
  const Real gridResistance = 0.8929;
  const Real gridInductance = 16.58e-3;

  // -----------------------------------------------------------------------
  // Existing, previously validated MMC controller values
  // -----------------------------------------------------------------------
  //
  // These are deliberately retained for the first topology test. They are
  // NOT the per-unit DQsym controller gains. The DQsym control architecture
  // will be implemented only after the physical grid/transformer topology
  // has been validated.
  const Real kpPll = 0.001103374;
  const Real kiPll = 0.00073;

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

  Logger::setLogDir("logs/" + simName);

  // -----------------------------------------------------------------------
  // Nodes
  // -----------------------------------------------------------------------
  auto nGridSource = SimNode<Real>::make("nGridSource", PhaseType::ABC);

  auto nGridAfterResistance =
      SimNode<Real>::make("nGridAfterResistance", PhaseType::ABC);

  auto nTransformerPrimary =
      SimNode<Real>::make("nTransformerPrimary", PhaseType::ABC);

  auto nMmcAc = SimNode<Real>::make("nMmcAc", PhaseType::ABC);

  auto nDcPositive = SimNode<Real>::make("nDcPositive", PhaseType::Single);

  auto nDcNegative = SimNode<Real>::make("nDcNegative", PhaseType::Single);

  // Zero-power initial state:
  // no grid/transformer voltage drop and nominal transformer ratio.
  const Complex primaryInitialPhasor =
      initialNodePhasorFromLineToLineRms(primaryVoltageRmsLl);

  const Complex secondaryInitialPhasor =
      initialNodePhasorFromLineToLineRms(secondaryVoltageRmsLl);

  nGridSource->setInitialVoltage(primaryInitialPhasor);
  nGridAfterResistance->setInitialVoltage(primaryInitialPhasor);
  nTransformerPrimary->setInitialVoltage(primaryInitialPhasor);
  nMmcAc->setInitialVoltage(secondaryInitialPhasor);

  nDcPositive->setInitialVoltage(Complex(+0.5 * nominalDcVoltage, 0.0));
  nDcNegative->setInitialVoltage(Complex(-0.5 * nominalDcVoltage, 0.0));

  // -----------------------------------------------------------------------
  // 400-kV AC source and series RL grid equivalent
  // -----------------------------------------------------------------------
  auto gridSource =
      EMT::Ph3::NetworkInjection::make("GridSource", Logger::Level::info);

  gridSource->setParameters(
      balancedVoltageReferenceFromLineToLineRms(primaryVoltageRmsLl),
      nominalFrequency);

  gridSource->connect({nGridSource});

  auto gridResistor =
      EMT::Ph3::Resistor::make("GridResistance", Logger::Level::info);

  gridResistor->setParameters(
      Math::singlePhaseParameterToThreePhase(gridResistance));

  gridResistor->connect({
      nGridSource,
      nGridAfterResistance,
  });

  auto gridInductor =
      EMT::Ph3::Inductor::make("GridInductance", Logger::Level::info);

  gridInductor->setParameters(
      Math::singlePhaseParameterToThreePhase(gridInductance));

  gridInductor->connect({
      nGridAfterResistance,
      nTransformerPrimary,
  });

  // -----------------------------------------------------------------------
  // 400/333-kV, 1000-MVA transformer
  // -----------------------------------------------------------------------
  auto transformer = EMT::Ph3::Transformer::make(
      "GridTransformer", "GridTransformer", Logger::Level::info, true);

  transformer->setParameters(
      primaryVoltageRmsLl, secondaryVoltageRmsLl, converterRatedPower,
      transformerRatio, 0.0,
      Math::singlePhaseParameterToThreePhase(transformerResistance),
      Math::singlePhaseParameterToThreePhase(transformerInductance));

  transformer->connect({
      nTransformerPrimary,
      nMmcAc,
  });

  // -----------------------------------------------------------------------
  // Stiff bipolar +/-320-kV DC terminal
  // -----------------------------------------------------------------------
  auto positiveDcSource =
      EMT::Ph1::VoltageSource::make("PositiveDcSource", Logger::Level::info);

  positiveDcSource->setParameters(Complex(+0.5 * nominalDcVoltage, 0.0), 0.0);

  positiveDcSource->connect({
      SimNode<Real>::GND,
      nDcPositive,
  });

  auto negativeDcSource =
      EMT::Ph1::VoltageSource::make("NegativeDcSource", Logger::Level::info);

  negativeDcSource->setParameters(Complex(+0.5 * nominalDcVoltage, 0.0), 0.0);

  // V(GND)-V(nDcNegative) = +320 kV.
  negativeDcSource->connect({
      nDcNegative,
      SimNode<Real>::GND,
  });

  // -----------------------------------------------------------------------
  // Existing SSN MMC
  // -----------------------------------------------------------------------
  auto mmc = EMT::Ph3::SSN_MMC::make("MMC", "MMC", Logger::Level::debug);

  mmc->setParameters(nominalFrequency, converterAcVoltageAmplitude,
                     nominalDcVoltage, armInductance, armResistance,
                     submoduleCapacitance, numberOfSubmodules,
                     reactorInductance, reactorResistance);

  mmc->setInitialAngle(0.0);
  mmc->setPLL(kpPll, kiPll, true);
  mmc->setReactiveControlOpenLoop(0.0);

  mmc->setActivePowerControl(0.0, kpActivePower, kiActivePower);

  mmc->setEnergyController(kpEnergy, kiEnergy, true);

  mmc->setOutputCurrentController(kpOutputCurrent, kiOutputCurrent);

  mmc->setCirculatingCurrentController(kpCirculatingCurrent,
                                       kiCirculatingCurrent);

  mmc->setZeroSequenceCurrentController(kpZeroSequenceCurrent,
                                        kiZeroSequenceCurrent);

  mmc->setCirculatingCurrentReferences(0.0, 0.0, 0.0);

  mmc->setLimits(5000.0, 2500.0, 2.0);

  mmc->setOperatingPointInitialization(true, 75, 1.0e-8);

  mmc->setEigenvalueDiagnostics(true, 200);

  mmc->setDiagnosticTimeStep(timeStep);

  mmc->connect({
      nMmcAc,
      nDcPositive,
      nDcNegative,
  });

  auto powerReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(0.0));

  // -----------------------------------------------------------------------
  // System
  // -----------------------------------------------------------------------
  auto system = SystemTopology(nominalFrequency,
                               SystemNodeList{
                                   nGridSource,
                                   nGridAfterResistance,
                                   nTransformerPrimary,
                                   nMmcAc,
                                   nDcPositive,
                                   nDcNegative,
                               },
                               SystemComponentList{
                                   gridSource,
                                   gridResistor,
                                   gridInductor,
                                   transformer,
                                   positiveDcSource,
                                   negativeDcSource,
                                   mmc,
                               });

  auto logger = DataLogger::make(simName);

  logger->logAttribute("Voltage_GridSource", nGridSource->attribute("v"));

  logger->logAttribute("Voltage_TransformerPrimary",
                       nTransformerPrimary->attribute("v"));

  logger->logAttribute("Voltage_MMC_AC", nMmcAc->attribute("v"));

  logger->logAttribute("Current_GridRL", gridInductor->attribute("i_intf"));

  logger->logAttribute("Current_Transformer", transformer->attribute("i_intf"));

  logger->logAttribute("Voltage_DC_Positive", nDcPositive->attribute("v"));

  logger->logAttribute("Voltage_DC_Negative", nDcNegative->attribute("v"));

  logger->logAttribute("ActivePowerReference", powerReferenceSignal);

  logger->logAttribute("MMC_InterfaceCurrent", mmc->attribute("i_intf"));

  logger->logAttribute("MMC_Pac", mmc->attribute("p_ac"));

  logger->logAttribute("MMC_Qac", mmc->attribute("q_ac"));

  logger->logAttribute("MMC_Pdc", mmc->attribute("p_dc"));

  logger->logAttribute("MMC_Vdc", mmc->attribute("vdc"));

  logger->logAttribute("MMC_Idc", mmc->attribute("idc"));

  logger->logAttribute("MMC_IDeltaD", mmc->attribute("i_delta_d"));

  logger->logAttribute("MMC_IDeltaDReference", mmc->attribute("i_delta_d_ref"));

  logger->logAttribute("MMC_PLLFrequency", mmc->attribute("pll_frequency"));

  logger->logAttribute("MMC_PLLAngleDeviation",
                       mmc->attribute("pll_angle_deviation"));

  logger->logAttribute("MMC_StoredEnergy", mmc->attribute("stored_energy"));

  logger->logAttribute("MMC_StateNorm", mmc->attribute("state_norm"));

  logger->logAttribute("MMC_StateDerivativeNorm",
                       mmc->attribute("state_derivative_norm"));

  logger->logAttribute("MMC_JacobianMaxRealEigenvalue",
                       mmc->attribute("jacobian_max_real_eigenvalue"));

  logger->logAttribute("MMC_JacobianMaxDiscreteMagnitude",
                       mmc->attribute("jacobian_max_discrete_magnitude"));

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

  Bool applyPowerStep = true;
  Bool restoreZeroPower = true;

  sim.initialize();
  sim.start();

  while (sim.time() < sim.finalTime()) {
    if (applyPowerStep && sim.time() >= powerStepTime) {
      SPDLOG_INFO("Applying +{:.3f} MW active-power reference at t={:.6f} s.",
                  activePowerStep / 1.0e6, sim.time());

      mmc->setActivePowerControl(activePowerStep, kpActivePower, kiActivePower);

      powerReferenceSignal->set(activePowerStep);

      applyPowerStep = false;
    }

    if (restoreZeroPower && sim.time() >= powerRestoreTime) {
      SPDLOG_INFO("Restoring zero active-power reference at t={:.6f} s.",
                  sim.time());

      mmc->setActivePowerControl(0.0, kpActivePower, kiActivePower);

      powerReferenceSignal->set(0.0);

      restoreZeroPower = false;
    }

    sim.step();
  }

  sim.stop();
  return 0;
}
