// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_Ph3_GFL_Siemens.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/SP/SP_Ph1_ControlledCurrentSource.h>
#include <dpsim-models/SP/SP_Ph1_Inductor.h>
#include <dpsim-models/SP/SP_Ph1_Resistor.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

namespace GflSiemensPqSetpointStep {

// =============================================================================
// Parameters
// =============================================================================

struct SimulationParameters {
  String name = "EMT_GFL_Siemens_PQ_SetpointStep";
  Real frequency = 50.0;
  Real timeStep = 20e-6;
  Real finalTime = 3.0;

  Real activePowerStepTime = 1.0;
  Real reactivePowerStepTime = 2.0;

  UInt logDownsampling = 10;
  Bool recomputeSystemMatrix = false;
};

/// Parameters transcribed from the MATLAB/Simulink CSI model.
struct SiemensControllerParameters {
  Real ratedApparentPower = 50e3; // total three-phase S_base [VA]
  Real ratedVoltage = 257.0;      // V_base line-to-line RMS [V]
  Real nominalFrequency = 50.0;   // f_base [Hz]

  // CSI inverse droop:
  //   P_cmd = P_ref - K_fP (f_PLL - f_ref)
  //   Q_cmd = Q_ref - K_VQ (V_PCC - V_ref)
  Real frequencyToActivePowerGainPu = 1.0 / 0.05;
  Real voltageToReactivePowerGainPu = 1.0 / 0.05;

  // SRF PLL parameters from the CSI subsystem mask.
  Real pllKp = 0.449 / (9.1e-3 + 250e-6);
  Real pllKi = 250e-6 / (5.9 * (9.1e-3 + 250e-6));

  // CSI dq current controller.
  Real currentControllerKp = 0.74 / 10.0;
  Real currentControllerKi = 11.9 / 10.0; // 1.19 / 10.0 from Simulink
  Real pccVoltageFeedforwardGain = 1.0;

  // Measurement filters. PT3/PT4 are reassigned to measured P/Q. The CSI
  // source contains no current-sensor time constant, so tau_I is an explicit
  // DPsim test parameter.
  Real activePowerMeasurementTimeConstant = 0.05;
  Real reactivePowerMeasurementTimeConstant = 0.05;
  Real currentMeasurementTimeConstant = 0.5e-3;

  // Physical converter-side Rf-Lf-Cf filter.
  Real filterResistance = 0.0279;
  Real filterInductance = 2.72e-4;
  Real filterCapacitance = 3.5e-6;

  // Protective test limits. These are not claimed as Siemens parameters.
  Real minimumFrequencyPu = 0.80;
  Real maximumFrequencyPu = 1.20;
  Real maximumCurrentReferencePu = 1.50;
  Real maximumVoltageCommandPu = 1.20;
  Real minimumVoltageForCurrentReferencePu = 0.05;

  Real baseImpedance() const {
    return ratedVoltage * ratedVoltage / ratedApparentPower;
  }

  Real baseOmega() const { return 2.0 * PI * nominalFrequency; }
};

/// Test commands:
///
///   t < t_P:              P_ref=0.50 pu, Q_ref=0.00 pu
///   t_P <= t < t_Q:       P_ref=0.75 pu, Q_ref=0.00 pu
///   t >= t_Q:             P_ref=0.75 pu, Q_ref=0.10 pu
///
/// The GFL is initialized from the first operating point.
struct CommandParameters {
  Real initialActivePowerPu = 0.50;
  Real initialReactivePowerPu = 0.00;

  Real activePowerAfterStepPu = 0.75;
  Real reactivePowerAfterStepPu = 0.10;
};

/// Stiff-grid connection:
///
///   ideal grid -- R_grid + L_grid -- PCC -- GFL_Siemens
///
/// The impedance is large enough to make the PLL and decoupling visible, but
/// sufficiently small for a first validation of the GFL controller.
struct GridParameters {
  Real voltagePu = 1.0;
  Real resistancePu = 0.01;
  Real reactancePu = 0.05;

  Real resistance(const SiemensControllerParameters &gfl) const {
    return resistancePu * gfl.baseImpedance();
  }

  Real inductance(const SiemensControllerParameters &gfl) const {
    return reactancePu * gfl.baseImpedance() / gfl.baseOmega();
  }
};

struct Parameters {
  SimulationParameters simulation;
  SiemensControllerParameters controller;
  CommandParameters commands;
  GridParameters grid;
};

// =============================================================================
// Helpers and validation
// =============================================================================

struct PowerFlowResult {
  SystemTopology system;
  Complex gflPower;
  Complex gridBusVoltage;
  Complex pccVoltage;
  Complex gridPower;
  Complex gflCurrent;
};

Complex
generatorPositivePower(const std::shared_ptr<SP::Ph1::NetworkInjection> &source,
                       const SimNode<Complex>::Ptr &sourceBus) {
  const Complex voltage = sourceBus->singleVoltage();
  const Complex consumerPositiveCurrent = (**source->mIntfCurrent)(0, 0);
  return -voltage * std::conj(consumerPositiveCurrent);
}

void validateParameters(const Parameters &p) {
  if (!(p.simulation.frequency > 0.0) || !(p.simulation.timeStep > 0.0) ||
      !(p.simulation.finalTime > 0.0) ||
      !(p.simulation.activePowerStepTime > 0.0) ||
      !(p.simulation.reactivePowerStepTime >
        p.simulation.activePowerStepTime) ||
      !(p.simulation.reactivePowerStepTime < p.simulation.finalTime) ||
      !(p.simulation.logDownsampling > 0)) {
    throw std::invalid_argument(
        "Require f>0, dt>0, finalTime>0, "
        "0<P_STEP_TIME<Q_STEP_TIME<finalTime and logDownsampling>0");
  }

  const auto &c = p.controller;
  if (!(c.ratedApparentPower > 0.0) || !(c.ratedVoltage > 0.0) ||
      !(c.nominalFrequency > 0.0) || !(c.filterInductance > 0.0) ||
      !(c.filterCapacitance > 0.0) || !(c.filterResistance >= 0.0)) {
    throw std::invalid_argument("Invalid GFL_Siemens parameters");
  }

  if (std::abs(p.simulation.frequency - c.nominalFrequency) > 1e-9) {
    throw std::invalid_argument(
        "Simulation frequency and GFL nominal frequency must match");
  }

  const auto finite = [](Real value) { return std::isfinite(value); };
  if (!finite(p.commands.initialActivePowerPu) ||
      !finite(p.commands.initialReactivePowerPu) ||
      !finite(p.commands.activePowerAfterStepPu) ||
      !finite(p.commands.reactivePowerAfterStepPu)) {
    throw std::invalid_argument("Non-finite GFL command value");
  }

  const Real initialMagnitude = std::hypot(p.commands.initialActivePowerPu,
                                           p.commands.initialReactivePowerPu);
  const Real pStepMagnitude = std::hypot(p.commands.activePowerAfterStepPu,
                                         p.commands.initialReactivePowerPu);
  const Real pqStepMagnitude = std::hypot(p.commands.activePowerAfterStepPu,
                                          p.commands.reactivePowerAfterStepPu);

  const Real commandMargin = 0.98 * c.maximumCurrentReferencePu;
  if (!(initialMagnitude < commandMargin) ||
      !(pStepMagnitude < commandMargin) || !(pqStepMagnitude < commandMargin)) {
    throw std::invalid_argument(
        "Requested P/Q operating point is too close to the configured "
        "current-reference limit");
  }

  if (!(p.grid.voltagePu > 0.0) || !(p.grid.resistancePu >= 0.0) ||
      !(p.grid.reactancePu > 0.0)) {
    throw std::invalid_argument("Invalid stiff-grid parameters");
  }
}

void logParameterSummary(const Parameters &p) {
  const auto &c = p.controller;
  const auto &g = p.grid;
  const auto &u = p.commands;

  SPDLOG_INFO(
      "GFL_Siemens P/Q setpoint-step test:"
      "\n  topology: ideal grid -- GRID_R -- GRID_L -- PCC -- GFL_Siemens"
      "\n  simulation: dt={} s, finalTime={} s, P-step={} s, Q-step={} s"
      "\n  base: S={} VA, V_LL={} V RMS, f={} Hz, Zbase={} Ohm"
      "\n  initial command: P={} pu ({} W), Q={} pu ({} var)"
      "\n  after P step: P={} pu ({} W), Q={} pu"
      "\n  after Q step: P={} pu, Q={} pu ({} var)"
      "\n  grid: V={} pu, R={} pu ({} Ohm), X={} pu, L={} H"
      "\n  inverse droop: K_fP={}, K_VQ={}"
      "\n  PLL: Kp={}, Ki={}"
      "\n  current PI: Kp={}, Ki={}, Kff_v={}"
      "\n  measurement filters: tauP={} s, tauQ={} s, tauI={} s"
      "\n  converter filter: Rf={} Ohm, Lf={} H, Cf={} F",
      p.simulation.timeStep, p.simulation.finalTime,
      p.simulation.activePowerStepTime, p.simulation.reactivePowerStepTime,
      c.ratedApparentPower, c.ratedVoltage, c.nominalFrequency,
      c.baseImpedance(), u.initialActivePowerPu,
      u.initialActivePowerPu * c.ratedApparentPower, u.initialReactivePowerPu,
      u.initialReactivePowerPu * c.ratedApparentPower, u.activePowerAfterStepPu,
      u.activePowerAfterStepPu * c.ratedApparentPower, u.initialReactivePowerPu,
      u.activePowerAfterStepPu, u.reactivePowerAfterStepPu,
      u.reactivePowerAfterStepPu * c.ratedApparentPower, g.voltagePu,
      g.resistancePu, g.resistance(c), g.reactancePu, g.inductance(c),
      c.frequencyToActivePowerGainPu, c.voltageToReactivePowerGainPu, c.pllKp,
      c.pllKi, c.currentControllerKp, c.currentControllerKi,
      c.pccVoltageFeedforwardGain, c.activePowerMeasurementTimeConstant,
      c.reactivePowerMeasurementTimeConstant, c.currentMeasurementTimeConstant,
      c.filterResistance, c.filterInductance, c.filterCapacitance);
}

Real alignedEventTime(Real requestedTime, Real timeStep) {
  return std::round(requestedTime / timeStep) * timeStep;
}

// =============================================================================
// SP power-flow initialization
// =============================================================================

PowerFlowResult buildAndRunPowerFlow(const Parameters &p) {
  const String simulationName = p.simulation.name + "_PF";
  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  // The standalone SP inductor is an MNA component rather than an NRP branch.
  // Therefore, this example first solves the exact two-bus constant-PQ
  // operating point and then validates/transfers it through an SP-MNA network
  // built from the same conventional R and L components as the EMT network.
  const Complex gridVoltage =
      Math::polar(p.grid.voltagePu * p.controller.ratedVoltage, 0.0);
  const Complex gridImpedance(p.grid.resistance(p.controller),
                              p.controller.baseOmega() *
                                  p.grid.inductance(p.controller));
  const Complex targetGflPower(
      p.commands.initialActivePowerPu * p.controller.ratedApparentPower,
      p.commands.initialReactivePowerPu * p.controller.ratedApparentPower);

  Complex pccVoltage = gridVoltage;
  constexpr UInt maximumIterations = 200;
  constexpr Real voltageTolerance = 1e-12;
  Bool converged = false;

  for (UInt iteration = 0; iteration < maximumIterations; ++iteration) {
    if (!(std::abs(pccVoltage) > 1e-9)) {
      throw std::runtime_error(
          "Two-bus power-flow iteration reached zero PCC voltage");
    }

    // Generator-positive convention:
    //   S_GFL = V_PCC * conj(I_GFL)
    //   V_PCC = V_GRID + Z_GRID * I_GFL
    const Complex gflCurrent = std::conj(targetGflPower / pccVoltage);
    const Complex updatedVoltage = gridVoltage + gridImpedance * gflCurrent;

    if (std::abs(updatedVoltage - pccVoltage) <=
        voltageTolerance * std::max<Real>(1.0, std::abs(updatedVoltage))) {
      pccVoltage = updatedVoltage;
      converged = true;
      break;
    }

    pccVoltage = updatedVoltage;
  }

  if (!converged) {
    throw std::runtime_error(
        "Two-bus constant-PQ power-flow initialization did not converge");
  }

  const Complex gflCurrent = std::conj(targetGflPower / pccVoltage);
  const Complex reconstructedPower = pccVoltage * std::conj(gflCurrent);
  const Complex reconstructedVoltage = gridVoltage + gridImpedance * gflCurrent;

  if (std::abs(reconstructedPower - targetGflPower) >
          1e-9 * p.controller.ratedApparentPower ||
      std::abs(reconstructedVoltage - pccVoltage) >
          1e-9 * p.controller.ratedVoltage) {
    throw std::runtime_error("Two-bus power-flow consistency check failed");
  }

  // Node names intentionally match the EMT topology.
  auto nGrid = SimNode<Complex>::make("BUS_GRID", PhaseType::Single);
  auto nAfterGridR =
      SimNode<Complex>::make("BUS_AFTER_GRID_R", PhaseType::Single);
  auto nPcc = SimNode<Complex>::make("BUS_PCC", PhaseType::Single);

  auto grid = SP::Ph1::NetworkInjection::make("GRID", Logger::Level::debug);
  grid->setParameters(gridVoltage, p.simulation.frequency);

  auto gridResistor = SP::Ph1::Resistor::make("GRID_R", Logger::Level::debug);
  gridResistor->setParameters(p.grid.resistance(p.controller));

  auto gridInductor = SP::Ph1::Inductor::make("GRID_L", Logger::Level::debug);
  gridInductor->setParameters(p.grid.inductance(p.controller));

  // Positive current flows from terminal 0 to terminal 1. Connecting the
  // source from ground to PCC therefore represents generator-positive current
  // injection by the GFL at the solved operating point.
  auto gflEquivalent = SP::Ph1::ControlledCurrentSource::make(
      "GFL_PF_EQUIVALENT", Logger::Level::debug);
  gflEquivalent->setParameters(gflCurrent);

  grid->connect({nGrid});
  gridResistor->connect({nGrid, nAfterGridR});
  gridInductor->connect({nAfterGridR, nPcc});
  gflEquivalent->connect({SimNode<Complex>::GND, nPcc});

  SystemTopology system(
      p.simulation.frequency, SystemNodeList{nGrid, nAfterGridR, nPcc},
      SystemComponentList{grid, gridResistor, gridInductor, gflEquivalent});

  auto logger = DataLogger::make(simulationName);
  logger->logAttribute("BUS_GRID.v", nGrid->attribute("v"));
  logger->logAttribute("BUS_AFTER_GRID_R.v", nAfterGridR->attribute("v"));
  logger->logAttribute("BUS_PCC.v", nPcc->attribute("v"));
  logger->logAttribute("GRID.i", grid->attribute("i_intf"));
  logger->logAttribute("GRID_R.i", gridResistor->attribute("i_intf"));
  logger->logAttribute("GRID_L.i", gridInductor->attribute("i_intf"));
  logger->logAttribute("GFL_PF_EQUIVALENT.i",
                       gflEquivalent->attribute("i_intf"));

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(1.0);
  simulation.setFinalTime(1.0);
  simulation.setDomain(Domain::SP);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.addLogger(logger);
  simulation.run();

  const Complex solvedGridVoltage = nGrid->singleVoltage();
  const Complex solvedPccVoltage = nPcc->singleVoltage();
  const Complex solvedGflPower = solvedPccVoltage * std::conj(gflCurrent);
  const Complex gridPower = generatorPositivePower(grid, nGrid);

  const auto finiteComplex = [](const Complex &value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  };

  if (!finiteComplex(solvedGridVoltage) || !finiteComplex(solvedPccVoltage) ||
      !finiteComplex(solvedGflPower) || !finiteComplex(gridPower) ||
      !(std::abs(solvedPccVoltage) > 0.0)) {
    throw std::runtime_error("Invalid GFL SP-MNA initialization result");
  }

  const Real powerErrorPu = std::abs(solvedGflPower - targetGflPower) /
                            p.controller.ratedApparentPower;
  const Real voltageErrorPu =
      std::abs(solvedPccVoltage - pccVoltage) / p.controller.ratedVoltage;

  if (powerErrorPu > 1e-8 || voltageErrorPu > 1e-8) {
    throw std::runtime_error(
        "SP-MNA operating point differs from the solved two-bus power flow");
  }

  SPDLOG_INFO("SP initialization solved with conventional series R and L:"
              "\n  GFL injection: P={} W = {} pu, Q={} var = {} pu"
              "\n  injection current: magnitude={} A, angle={} deg"
              "\n  grid injection: P={} W, Q={} var"
              "\n  |V_GRID|={} V_LL RMS = {} pu, angle={} deg"
              "\n  |V_PCC|={} V_LL RMS = {} pu, angle={} deg"
              "\n  PCC-grid angle difference={} deg"
              "\n  PF/SP-MNA errors: |dS|={} pu, |dV|={} pu",
              solvedGflPower.real(),
              solvedGflPower.real() / p.controller.ratedApparentPower,
              solvedGflPower.imag(),
              solvedGflPower.imag() / p.controller.ratedApparentPower,
              std::abs(gflCurrent), std::arg(gflCurrent) * 180.0 / PI,
              gridPower.real(), gridPower.imag(), std::abs(solvedGridVoltage),
              std::abs(solvedGridVoltage) / p.controller.ratedVoltage,
              std::arg(solvedGridVoltage) * 180.0 / PI,
              std::abs(solvedPccVoltage),
              std::abs(solvedPccVoltage) / p.controller.ratedVoltage,
              std::arg(solvedPccVoltage) * 180.0 / PI,
              std::arg(solvedPccVoltage / solvedGridVoltage) * 180.0 / PI,
              powerErrorPu, voltageErrorPu);

  return {system,           solvedGflPower, solvedGridVoltage,
          solvedPccVoltage, gridPower,      gflCurrent};
}

// =============================================================================
// GFL configuration
// =============================================================================

void configureGfl(const std::shared_ptr<EMT::Ph3::GFL_Siemens> &gfl,
                  const Parameters &p, const PowerFlowResult &powerFlow) {
  const auto &c = p.controller;

  const Real initialVoltagePu = std::abs(powerFlow.pccVoltage) / c.ratedVoltage;
  const Real initialActivePowerPu =
      powerFlow.gflPower.real() / c.ratedApparentPower;
  const Real initialReactivePowerPu =
      powerFlow.gflPower.imag() / c.ratedApparentPower;

  gfl->setBaseParameters(c.ratedApparentPower, c.ratedVoltage,
                         c.nominalFrequency);

  // V_ref is initialized to the solved PCC voltage. Consequently the inverse
  // Q-V droop contributes zero bias at t=0, while the GFL injects the exact
  // initial P/Q operating point from the power flow.
  gfl->setReferencesPerUnit(1.0, initialVoltagePu, initialActivePowerPu,
                            initialReactivePowerPu);

  gfl->setDroopParametersPerUnit(c.frequencyToActivePowerGainPu,
                                 c.voltageToReactivePowerGainPu);
  gfl->setPllParameters(c.pllKp, c.pllKi);
  gfl->setCurrentControllerParameters(c.currentControllerKp,
                                      c.currentControllerKi,
                                      c.pccVoltageFeedforwardGain);
  gfl->setMeasurementFilterTimeConstants(c.activePowerMeasurementTimeConstant,
                                         c.reactivePowerMeasurementTimeConstant,
                                         c.currentMeasurementTimeConstant);
  gfl->setControllerLimitsPerUnit(
      c.minimumFrequencyPu, c.maximumFrequencyPu, c.maximumCurrentReferencePu,
      c.maximumVoltageCommandPu, c.minimumVoltageForCurrentReferencePu);
  gfl->setFilterParameters(c.filterInductance, c.filterCapacitance,
                           c.filterResistance);
  gfl->withControl(true);

  SPDLOG_INFO("Configured GFL_Siemens from SP operating point:"
              "\n  f_ref=1 pu, V_ref={} pu, P_ref={} pu, Q_ref={} pu",
              initialVoltagePu, initialActivePowerPu, initialReactivePowerPu);
}

// =============================================================================
// EMT simulation
// =============================================================================

void runEmt(const Parameters &p, const PowerFlowResult &powerFlow) {
  const String simulationName = p.simulation.name + "_EMT";
  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto nGrid = SimNode<Real>::make("BUS_GRID", PhaseType::ABC);
  auto nAfterGridR = SimNode<Real>::make("BUS_AFTER_GRID_R", PhaseType::ABC);
  auto nPcc = SimNode<Real>::make("BUS_PCC", PhaseType::ABC);

  auto grid = EMT::Ph3::NetworkInjection::make("GRID", Logger::Level::debug);

  // EMT::Ph3::NetworkInjection / EMT::Ph3::VoltageSource expects a
  // three-phase vector whose entries are line-to-line RMS phasors. The
  // VoltageSource converts each entry internally with RMS3PH_TO_PEAK1PH.
  // Reuse the exact SP grid phasor so the EMT source starts at precisely the
  // same magnitude and angle as the initialization network.
  MatrixComp gridVoltageRef = MatrixComp::Zero(3, 1);
  gridVoltageRef(0, 0) = powerFlow.gridBusVoltage;
  gridVoltageRef(1, 0) = gridVoltageRef(0, 0) * SHIFT_TO_PHASE_B;
  gridVoltageRef(2, 0) = gridVoltageRef(0, 0) * SHIFT_TO_PHASE_C;
  grid->setParameters(gridVoltageRef, p.simulation.frequency);

  SPDLOG_INFO("EMT grid source initialized from SP phasor: |V|={} V_LL RMS, "
              "angle={} deg",
              std::abs(powerFlow.gridBusVoltage),
              Math::phase(powerFlow.gridBusVoltage) * 180.0 / PI);

  auto gridResistor = EMT::Ph3::Resistor::make("GRID_R", Logger::Level::debug);
  gridResistor->setParameters(
      Math::singlePhaseParameterToThreePhase(p.grid.resistance(p.controller)));

  auto gridInductor = EMT::Ph3::Inductor::make("GRID_L", Logger::Level::debug);
  gridInductor->setParameters(
      Math::singlePhaseParameterToThreePhase(p.grid.inductance(p.controller)));

  auto gfl = EMT::Ph3::GFL_Siemens::make("GFL", "GFL", Logger::Level::debug);
  configureGfl(gfl, p, powerFlow);

  grid->connect({nGrid});
  gridResistor->connect({nGrid, nAfterGridR});
  gridInductor->connect({nAfterGridR, nPcc});
  gfl->connect({nPcc});

  SystemTopology system(
      p.simulation.frequency, SystemNodeList{nGrid, nAfterGridR, nPcc},
      SystemComponentList{grid, gridResistor, gridInductor, gfl});

  system.initWithPowerflow(powerFlow.system, Domain::EMT);
  gfl->terminal(0)->setPower(powerFlow.gflPower);

  auto logger =
      DataLogger::make(simulationName, true, p.simulation.logDownsampling);

  // Physical network.
  logger->logAttribute("BUS_GRID.v", nGrid->attribute("v"));
  logger->logAttribute("BUS_AFTER_GRID_R.v", nAfterGridR->attribute("v"));
  logger->logAttribute("BUS_PCC.v", nPcc->attribute("v"));
  logger->logAttribute("GRID.i", grid->attribute("i_intf"));
  logger->logAttribute("GRID_R.i", gridResistor->attribute("i_intf"));
  logger->logAttribute("GRID_L.i", gridInductor->attribute("i_intf"));

  // References, inverse-droop commands and measured power.
  logger->logAttribute("GFL.P_ref_pu", gfl->attribute("P_ref_pu"));
  logger->logAttribute("GFL.Q_ref_pu", gfl->attribute("Q_ref_pu"));
  logger->logAttribute("GFL.V_ref_pu", gfl->attribute("V_ref_pu"));
  logger->logAttribute("GFL.P_command_pu", gfl->attribute("P_command_pu"));
  logger->logAttribute("GFL.Q_command_pu", gfl->attribute("Q_command_pu"));
  logger->logAttribute("GFL.P_elec_pu", gfl->attribute("P_elec_pu"));
  logger->logAttribute("GFL.Q_elec_pu", gfl->attribute("Q_elec_pu"));
  logger->logAttribute("GFL.P_filtered_pu", gfl->attribute("P_filtered_pu"));
  logger->logAttribute("GFL.Q_filtered_pu", gfl->attribute("Q_filtered_pu"));

  // PLL.
  logger->logAttribute("GFL.frequency_pu", gfl->attribute("frequency_pu"));
  logger->logAttribute("GFL.V_magnitude_pu", gfl->attribute("V_magnitude_pu"));
  logger->logAttribute("GFL.pll_vq_error_pu",
                       gfl->attribute("pll_vq_error_pu"));
  logger->logAttribute("GFL.pll_integrator", gfl->attribute("pll_integrator"));
  logger->logAttribute("GFL.theta", gfl->attribute("theta"));

  // dq current control.
  logger->logAttribute("GFL.v_pcc_dq_pu", gfl->attribute("v_pcc_dq_pu"));
  logger->logAttribute("GFL.i_pcc_dq_pu", gfl->attribute("i_pcc_dq_pu"));
  logger->logAttribute("GFL.i_pcc_filtered_dq_pu",
                       gfl->attribute("i_pcc_filtered_dq_pu"));
  logger->logAttribute("GFL.i_ref_dq_pu", gfl->attribute("i_ref_dq_pu"));
  logger->logAttribute("GFL.e_i_dq_pu", gfl->attribute("e_i_dq_pu"));
  logger->logAttribute("GFL.xi_i_dq_pu", gfl->attribute("xi_i_dq_pu"));
  logger->logAttribute("GFL.v_cmd_dq_pu", gfl->attribute("v_cmd_dq_pu"));

  // Internal converter/filter quantities.
  logger->logAttribute("GFL.i_pcc", gfl->attribute("i_pcc"));
  logger->logAttribute("GFL.i_lf", gfl->attribute("i_lf"));
  logger->logAttribute("GFL.i_cf", gfl->attribute("i_cf"));
  logger->logAttribute("GFL.Vsref_pu", gfl->attribute("Vsref_pu"));
  logger->logAttribute("GFL.Vsref", gfl->attribute("Vsref"));

  const Real pStepTime =
      alignedEventTime(p.simulation.activePowerStepTime, p.simulation.timeStep);
  const Real qStepTime = alignedEventTime(p.simulation.reactivePowerStepTime,
                                          p.simulation.timeStep);

  auto activePowerStep = AttributeEvent<Real>::make(
      pStepTime, gfl->mActivePowerRefPu, p.commands.activePowerAfterStepPu);
  auto reactivePowerStep = AttributeEvent<Real>::make(
      qStepTime, gfl->mReactivePowerRefPu, p.commands.reactivePowerAfterStepPu);

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(p.simulation.timeStep);
  simulation.setFinalTime(p.simulation.finalTime);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doSystemMatrixRecomputation(p.simulation.recomputeSystemMatrix);
  simulation.addEvent(activePowerStep);
  simulation.addEvent(reactivePowerStep);
  simulation.addLogger(logger);

  SPDLOG_INFO("Starting GFL EMT run: P_ref steps {} -> {} pu at t={} s; "
              "Q_ref steps {} -> {} pu at t={} s.",
              p.commands.initialActivePowerPu,
              p.commands.activePowerAfterStepPu, pStepTime,
              p.commands.initialReactivePowerPu,
              p.commands.reactivePowerAfterStepPu, qStepTime);

  simulation.run();

  SPDLOG_INFO("Simulation completed. Results: logs/{}/{}.csv", simulationName,
              simulationName);
}

} // namespace GflSiemensPqSetpointStep

int main(int argc, char *argv[]) {
  try {
    GflSiemensPqSetpointStep::Parameters parameters;
    CommandLineArgs args(argc, argv);

    if (argc > 1) {
      parameters.simulation.timeStep = args.timeStep;
      parameters.simulation.finalTime = args.duration;

      if (args.name != "dpsim")
        parameters.simulation.name = args.name;

      if (args.options.find("P_STEP_TIME") != args.options.end()) {
        parameters.simulation.activePowerStepTime =
            args.getOptionReal("P_STEP_TIME");
      }

      if (args.options.find("Q_STEP_TIME") != args.options.end()) {
        parameters.simulation.reactivePowerStepTime =
            args.getOptionReal("Q_STEP_TIME");
      }

      if (args.options.find("INITIAL_P_PU") != args.options.end()) {
        parameters.commands.initialActivePowerPu =
            args.getOptionReal("INITIAL_P_PU");
      }

      if (args.options.find("INITIAL_Q_PU") != args.options.end()) {
        parameters.commands.initialReactivePowerPu =
            args.getOptionReal("INITIAL_Q_PU");
      }

      if (args.options.find("P_AFTER_STEP_PU") != args.options.end()) {
        parameters.commands.activePowerAfterStepPu =
            args.getOptionReal("P_AFTER_STEP_PU");
      }

      if (args.options.find("Q_AFTER_STEP_PU") != args.options.end()) {
        parameters.commands.reactivePowerAfterStepPu =
            args.getOptionReal("Q_AFTER_STEP_PU");
      }

      if (args.options.find("CURRENT_FILTER_TAU") != args.options.end()) {
        parameters.controller.currentMeasurementTimeConstant =
            args.getOptionReal("CURRENT_FILTER_TAU");
      }

      if (args.options.find("GRID_R_PU") != args.options.end()) {
        parameters.grid.resistancePu = args.getOptionReal("GRID_R_PU");
      }

      if (args.options.find("GRID_X_PU") != args.options.end()) {
        parameters.grid.reactancePu = args.getOptionReal("GRID_X_PU");
      }

      if (args.options.find("LOG_DOWNSAMPLING") != args.options.end()) {
        const Real value = args.getOptionReal("LOG_DOWNSAMPLING");
        if (!(value >= 1.0) || std::floor(value) != value ||
            value > static_cast<Real>(std::numeric_limits<UInt>::max())) {
          throw std::invalid_argument(
              "LOG_DOWNSAMPLING must be an integer >= 1");
        }
        parameters.simulation.logDownsampling = static_cast<UInt>(value);
      }
    }

    GflSiemensPqSetpointStep::validateParameters(parameters);
    GflSiemensPqSetpointStep::logParameterSummary(parameters);

    const auto powerFlow =
        GflSiemensPqSetpointStep::buildAndRunPowerFlow(parameters);

    GflSiemensPqSetpointStep::runEmt(parameters, powerFlow);
    return EXIT_SUCCESS;
  } catch (const std::exception &exception) {
    SPDLOG_ERROR("GFL_Siemens P/Q setpoint-step example failed:\n{}",
                 exception.what());
    return EXIT_FAILURE;
  }
}
