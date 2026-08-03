// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_Ph3_GFM_Siemens.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

namespace GfmSiemensInitialLoadStep {

// =============================================================================
// Parameters
// =============================================================================

struct SimulationParameters {
  String name = "EMT_GFM_Siemens_LoadStep_Initial_0_5pu";
  Real frequency = 50.0;
  Real timeStep = 20e-6;
  Real finalTime = 3.0;
  Real loadStepTime = 1.0;
  UInt logDownsampling = 10;
  Bool recomputeSystemMatrix = true;
};

/// Parameters transcribed from the MATLAB/Simulink New VSI Model.
struct SiemensControllerParameters {
  Real ratedApparentPower = 50e3; // total three-phase S_base [VA]
  Real ratedVoltage = 257.0;      // V_base line-to-line RMS [V]
  Real nominalFrequency = 50.0;   // f_base [Hz]

  // P-f and Q-V droop.
  Real activePowerDroopPu = 0.02;
  Real reactivePowerDroopPu = 0.0311;

  // P/Q measurement filters.
  Real activePowerMeasurementTimeConstant = 0.1;
  Real reactivePowerMeasurementTimeConstant = 0.1;

  // dq capacitor-voltage controller.
  Real voltageControllerKp = 0.52;
  Real voltageControllerKi = 1.16;
  Real outputCurrentFeedforwardGain = 1.0;

  // dq filter-current controller.
  Real currentControllerKp = 0.74;
  Real currentControllerKi = 1.19;
  Real pccVoltageFeedforwardGain = 1.0;

  // Original Simulink PWM delay. Bypassed by default for the EMT time step.
  Real pwmDelayTimeConstantSlx = 1e-6;
  Bool usePwmDelay = false;

  // Physical converter-side Rf-Lf-Cf filter.
  Real filterResistance = 0.0279;
  Real filterInductance = 2.72e-4;
  Real filterCapacitance = 3.5e-6;

  // Protective test limits, not claimed as Siemens/Simulink parameters.
  Real minimumFrequencyPu = 0.80;
  Real maximumFrequencyPu = 1.20;
  Real maximumCurrentReferencePu = 2.0;
  Real maximumVoltageCommandPu = 1.20;

  Real baseImpedance() const {
    return ratedVoltage * ratedVoltage / ratedApparentPower;
  }

  Real baseOmega() const { return 2.0 * PI * nominalFrequency; }

  Real appliedPwmDelayTimeConstant() const {
    return usePwmDelay ? pwmDelayTimeConstantSlx : 0.0;
  }
};

/// External test network:
///
///                           +-- R_base -- L_base -- GND
///                           |
///   GFM -- R_grid -- L_grid +-- breaker -- R_step -- L_step -- GND
///
/// The SP initialization contains only the permanently connected base load.
/// Its impedance is selected so that the GFM terminal supplies exactly the
/// requested initial P+jQ at 1 pu GFM voltage, including the grid R-L losses.
struct NetworkParameters {
  Real gridResistancePu = 0.01;
  Real gridReactancePu = 0.05;

  // Desired initial GFM terminal power before the load step.
  Real initialGfmActivePowerPu = 0.50;
  Real initialGfmReactivePowerPu = 0.15;

  // Additional switched constant-impedance load, specified at 1 pu load-bus
  // voltage. The actual step power depends on the post-step bus voltage.
  Real stepLoadActivePowerPu = 0.25;
  Real stepLoadReactivePowerPu = 0.075;

  Real breakerOpenResistance = 1e12;
  Real breakerClosedResistancePu = 1e-4;

  static Real impedanceResistancePuForPower(Real activePowerPu,
                                            Real reactivePowerPu) {
    const Real denominator =
        activePowerPu * activePowerPu + reactivePowerPu * reactivePowerPu;
    if (!(denominator > 0.0))
      throw std::invalid_argument("RL load P and Q cannot both be zero");
    return activePowerPu / denominator;
  }

  static Real impedanceReactancePuForPower(Real activePowerPu,
                                           Real reactivePowerPu) {
    const Real denominator =
        activePowerPu * activePowerPu + reactivePowerPu * reactivePowerPu;
    if (!(denominator > 0.0))
      throw std::invalid_argument("RL load P and Q cannot both be zero");
    return reactivePowerPu / denominator;
  }

  /// Total series impedance seen from the ideal 1 pu source required for
  /// S_source = P_initial + j Q_initial:
  ///
  ///   Z_total,pu = 1 / conj(S_initial,pu)
  ///              = (P_initial + j Q_initial) / |S_initial|^2.
  Real initialTotalResistancePu() const {
    return impedanceResistancePuForPower(initialGfmActivePowerPu,
                                         initialGfmReactivePowerPu);
  }

  Real initialTotalReactancePu() const {
    return impedanceReactancePuForPower(initialGfmActivePowerPu,
                                        initialGfmReactivePowerPu);
  }

  /// The permanent load is the total required impedance minus the shared grid
  /// impedance. This makes the SP source power exactly the requested initial
  /// power at a 1 pu source voltage.
  Real baseLoadResistancePu() const {
    return initialTotalResistancePu() - gridResistancePu;
  }

  Real baseLoadReactancePu() const {
    return initialTotalReactancePu() - gridReactancePu;
  }

  Real stepLoadResistancePu() const {
    return impedanceResistancePuForPower(stepLoadActivePowerPu,
                                         stepLoadReactivePowerPu);
  }

  Real stepLoadReactancePu() const {
    return impedanceReactancePuForPower(stepLoadActivePowerPu,
                                        stepLoadReactivePowerPu);
  }

  Real gridResistance(const SiemensControllerParameters &gfm) const {
    return gridResistancePu * gfm.baseImpedance();
  }

  Real gridInductance(const SiemensControllerParameters &gfm) const {
    return gridReactancePu * gfm.baseImpedance() / gfm.baseOmega();
  }

  Real baseLoadResistance(const SiemensControllerParameters &gfm) const {
    return baseLoadResistancePu() * gfm.baseImpedance();
  }

  Real baseLoadInductance(const SiemensControllerParameters &gfm) const {
    return baseLoadReactancePu() * gfm.baseImpedance() / gfm.baseOmega();
  }

  Real stepLoadResistance(const SiemensControllerParameters &gfm) const {
    return stepLoadResistancePu() * gfm.baseImpedance();
  }

  Real stepLoadInductance(const SiemensControllerParameters &gfm) const {
    return stepLoadReactancePu() * gfm.baseImpedance() / gfm.baseOmega();
  }

  Real breakerClosedResistance(const SiemensControllerParameters &gfm) const {
    return breakerClosedResistancePu * gfm.baseImpedance();
  }
};

struct Parameters {
  SimulationParameters simulation;
  SiemensControllerParameters controller;
  NetworkParameters network;
};

// =============================================================================
// Helpers and validation
// =============================================================================

struct PowerFlowResult {
  SystemTopology system;
  Complex gfmPower;
  Complex gfmBusVoltage;
  Complex loadBusVoltage;
  Complex stepBreakerLoadSideVoltage;
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
      !(p.simulation.finalTime > 0.0) || !(p.simulation.loadStepTime > 0.0) ||
      !(p.simulation.loadStepTime < p.simulation.finalTime) ||
      !(p.simulation.logDownsampling > 0)) {
    throw std::invalid_argument(
        "Require f>0, dt>0, finalTime>0, "
        "0<loadStepTime<finalTime and logDownsampling>0");
  }

  if (!(p.controller.ratedApparentPower > 0.0) ||
      !(p.controller.ratedVoltage > 0.0) ||
      !(p.controller.nominalFrequency > 0.0)) {
    throw std::invalid_argument("Invalid GFM_Siemens base parameters");
  }

  if (std::abs(p.simulation.frequency - p.controller.nominalFrequency) > 1e-9) {
    throw std::invalid_argument(
        "Simulation frequency and GFM nominal frequency must match");
  }

  if (!(p.network.initialGfmActivePowerPu > 0.0) ||
      !(p.network.initialGfmReactivePowerPu >= 0.0) ||
      !(p.network.stepLoadActivePowerPu >= 0.0) ||
      !(p.network.stepLoadReactivePowerPu >= 0.0)) {
    throw std::invalid_argument(
        "This RL-load example requires positive initial P and non-negative Q");
  }

  if (!(p.network.baseLoadResistancePu() > 0.0) ||
      !(p.network.baseLoadReactancePu() > 0.0)) {
    throw std::invalid_argument(
        "Initial target power and grid impedance yield a non-positive "
        "base-load R or X");
  }

  if (!(p.network.stepLoadResistancePu() > 0.0) ||
      !(p.network.stepLoadReactancePu() >= 0.0)) {
    throw std::invalid_argument("Invalid switched RL load impedance");
  }
}

void logParameterSummary(const Parameters &p) {
  const auto &c = p.controller;
  const auto &n = p.network;

  SPDLOG_INFO(
      "GFM_Siemens nonzero-initial-power load step:"
      "\n  topology: GFM--Rg--Lg--BUS_LOAD with permanent and switched RL "
      "branches"
      "\n  simulation: dt={} s, finalTime={} s, stepTime={} s"
      "\n  base: S={} VA, V_LL={} V RMS, f={} Hz, Z={} Ohm"
      "\n  initial GFM target: P={} pu ({} W), Q={} pu ({} var)"
      "\n  grid: R={} pu ({} Ohm), X={} pu, L={} H"
      "\n  permanent load: R={} pu ({} Ohm), X={} pu, L={} H"
      "\n  switched load nominal at 1 pu load-bus voltage: "
      "P={} pu, Q={} pu, R={} pu ({} Ohm), X={} pu, L={} H"
      "\n  controller: kp={}, kq={}, tauP={} s, tauQ={} s"
      "\n  voltage PI: Kp={}, Ki={}, Kff_i={}"
      "\n  current PI: Kp={}, Ki={}, Kff_v={}"
      "\n  converter filter: Rf={} Ohm, Lf={} H, Cf={} F",
      p.simulation.timeStep, p.simulation.finalTime, p.simulation.loadStepTime,
      c.ratedApparentPower, c.ratedVoltage, c.nominalFrequency,
      c.baseImpedance(), n.initialGfmActivePowerPu,
      n.initialGfmActivePowerPu * c.ratedApparentPower,
      n.initialGfmReactivePowerPu,
      n.initialGfmReactivePowerPu * c.ratedApparentPower, n.gridResistancePu,
      n.gridResistance(c), n.gridReactancePu, n.gridInductance(c),
      n.baseLoadResistancePu(), n.baseLoadResistance(c),
      n.baseLoadReactancePu(), n.baseLoadInductance(c), n.stepLoadActivePowerPu,
      n.stepLoadReactivePowerPu, n.stepLoadResistancePu(),
      n.stepLoadResistance(c), n.stepLoadReactancePu(), n.stepLoadInductance(c),
      c.activePowerDroopPu, c.reactivePowerDroopPu,
      c.activePowerMeasurementTimeConstant,
      c.reactivePowerMeasurementTimeConstant, c.voltageControllerKp,
      c.voltageControllerKi, c.outputCurrentFeedforwardGain,
      c.currentControllerKp, c.currentControllerKi, c.pccVoltageFeedforwardGain,
      c.filterResistance, c.filterInductance, c.filterCapacitance);
}

// =============================================================================
// SP initialization
// =============================================================================

PowerFlowResult buildAndRunPowerFlow(const Parameters &p) {
  const String simulationName = p.simulation.name + "_PF";
  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  // Node names intentionally match the EMT topology.
  auto nGfm = SimNode<Complex>::make("BUS_GFM", PhaseType::Single);
  auto nAfterGridR =
      SimNode<Complex>::make("BUS_AFTER_GRID_R", PhaseType::Single);
  auto nLoadBus = SimNode<Complex>::make("BUS_LOAD", PhaseType::Single);
  auto nAfterBaseR =
      SimNode<Complex>::make("BUS_AFTER_BASE_R", PhaseType::Single);
  auto nStepBreakerLoad =
      SimNode<Complex>::make("BUS_STEP_BREAKER_LOAD", PhaseType::Single);
  auto nAfterStepR =
      SimNode<Complex>::make("BUS_AFTER_STEP_R", PhaseType::Single);

  auto source = SP::Ph1::NetworkInjection::make("GFM", Logger::Level::debug);
  source->setParameters(Math::polar(p.controller.ratedVoltage, 0.0),
                        p.simulation.frequency);

  auto gridResistor = SP::Ph1::Resistor::make("GRID_R", Logger::Level::debug);
  gridResistor->setParameters(p.network.gridResistance(p.controller));

  auto gridInductor = SP::Ph1::Inductor::make("GRID_L", Logger::Level::debug);
  gridInductor->setParameters(p.network.gridInductance(p.controller));

  auto baseLoadResistor =
      SP::Ph1::Resistor::make("BASE_LOAD_R", Logger::Level::debug);
  baseLoadResistor->setParameters(p.network.baseLoadResistance(p.controller));

  auto baseLoadInductor =
      SP::Ph1::Inductor::make("BASE_LOAD_L", Logger::Level::debug);
  baseLoadInductor->setParameters(p.network.baseLoadInductance(p.controller));

  auto stepBreaker =
      SP::Ph1::Switch::make("STEP_BREAKER", Logger::Level::debug);
  stepBreaker->setParameters(p.network.breakerOpenResistance,
                             p.network.breakerClosedResistance(p.controller),
                             false);

  auto stepLoadResistor =
      SP::Ph1::Resistor::make("STEP_LOAD_R", Logger::Level::debug);
  stepLoadResistor->setParameters(p.network.stepLoadResistance(p.controller));

  auto stepLoadInductor =
      SP::Ph1::Inductor::make("STEP_LOAD_L", Logger::Level::debug);
  stepLoadInductor->setParameters(p.network.stepLoadInductance(p.controller));

  source->connect({nGfm});
  gridResistor->connect({nGfm, nAfterGridR});
  gridInductor->connect({nAfterGridR, nLoadBus});

  // Permanently connected base load.
  baseLoadResistor->connect({nLoadBus, nAfterBaseR});
  baseLoadInductor->connect({nAfterBaseR, SimNode<Complex>::GND});

  // Initially disconnected step load.
  stepBreaker->connect({nLoadBus, nStepBreakerLoad});
  stepLoadResistor->connect({nStepBreakerLoad, nAfterStepR});
  stepLoadInductor->connect({nAfterStepR, SimNode<Complex>::GND});

  SystemTopology system(p.simulation.frequency,
                        SystemNodeList{nGfm, nAfterGridR, nLoadBus, nAfterBaseR,
                                       nStepBreakerLoad, nAfterStepR},
                        SystemComponentList{source, gridResistor, gridInductor,
                                            baseLoadResistor, baseLoadInductor,
                                            stepBreaker, stepLoadResistor,
                                            stepLoadInductor});

  auto logger = DataLogger::make(simulationName);
  logger->logAttribute("BUS_GFM.v", nGfm->attribute("v"));
  logger->logAttribute("BUS_LOAD.v", nLoadBus->attribute("v"));
  logger->logAttribute("GFM.i", source->attribute("i_intf"));
  logger->logAttribute("BASE_LOAD_R.i", baseLoadResistor->attribute("i_intf"));
  logger->logAttribute("STEP_BREAKER.i", stepBreaker->attribute("i_intf"));

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(1.0);
  simulation.setFinalTime(1.0);
  simulation.setDomain(Domain::SP);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.addLogger(logger);
  simulation.run();

  const Complex gfmPower = generatorPositivePower(source, nGfm);
  const Complex gfmVoltage = nGfm->singleVoltage();
  const Complex loadBusVoltage = nLoadBus->singleVoltage();
  const Complex stepLoadSideVoltage = nStepBreakerLoad->singleVoltage();

  const Real initialPPu = gfmPower.real() / p.controller.ratedApparentPower;
  const Real initialQPu = gfmPower.imag() / p.controller.ratedApparentPower;

  if (!std::isfinite(initialPPu) || !std::isfinite(initialQPu) ||
      !(std::abs(gfmVoltage) > 0.0)) {
    throw std::runtime_error("Invalid SP initialization result");
  }

  SPDLOG_INFO("SP initialization solved with permanent base load:"
              "\n  GFM P={} W = {} pu (target {} pu)"
              "\n  GFM Q={} var = {} pu (target {} pu)"
              "\n  |V_GFM|={} V_LL RMS = {} pu"
              "\n  |V_LOAD|={} V_LL RMS = {} pu"
              "\n  open step-branch load-side voltage={} V_LL RMS",
              gfmPower.real(), initialPPu, p.network.initialGfmActivePowerPu,
              gfmPower.imag(), initialQPu, p.network.initialGfmReactivePowerPu,
              std::abs(gfmVoltage),
              std::abs(gfmVoltage) / p.controller.ratedVoltage,
              std::abs(loadBusVoltage),
              std::abs(loadBusVoltage) / p.controller.ratedVoltage,
              std::abs(stepLoadSideVoltage));

  const Real powerTolerancePu = 1e-6;
  if (std::abs(initialPPu - p.network.initialGfmActivePowerPu) >
          powerTolerancePu ||
      std::abs(initialQPu - p.network.initialGfmReactivePowerPu) >
          powerTolerancePu) {
    SPDLOG_WARN(
        "SP source power differs from its analytical target by more than {} "
        "pu. "
        "This can indicate a parameter/sign convention difference in the "
        "local DPsim branch.",
        powerTolerancePu);
  }

  return {system, gfmPower, gfmVoltage, loadBusVoltage, stepLoadSideVoltage};
}

// =============================================================================
// GFM configuration
// =============================================================================

void configureGfm(const std::shared_ptr<EMT::Ph3::GFM_Siemens> &gfm,
                  const SiemensControllerParameters &parameters,
                  const PowerFlowResult &powerFlow) {
  const Real initialVoltagePu =
      std::abs(powerFlow.gfmBusVoltage) / parameters.ratedVoltage;
  const Real initialActivePowerPu =
      powerFlow.gfmPower.real() / parameters.ratedApparentPower;
  const Real initialReactivePowerPu =
      powerFlow.gfmPower.imag() / parameters.ratedApparentPower;

  gfm->setBaseParameters(parameters.ratedApparentPower, parameters.ratedVoltage,
                         parameters.nominalFrequency);

  // The droop references equal the nonzero SP operating point, so the EMT run
  // starts at 1 pu frequency while delivering approximately 0.5 pu active
  // power. The load step is a disturbance around this operating point.
  gfm->setReferencesPerUnit(1.0, initialVoltagePu, initialActivePowerPu,
                            initialReactivePowerPu);

  gfm->setDroopParametersPerUnit(parameters.activePowerDroopPu,
                                 parameters.reactivePowerDroopPu);

  gfm->setPowerMeasurementFilterTimeConstants(
      parameters.activePowerMeasurementTimeConstant,
      parameters.reactivePowerMeasurementTimeConstant);

  gfm->setVoltageControllerParameters(parameters.voltageControllerKp,
                                      parameters.voltageControllerKi,
                                      parameters.outputCurrentFeedforwardGain);

  gfm->setCurrentControllerParameters(parameters.currentControllerKp,
                                      parameters.currentControllerKi,
                                      parameters.pccVoltageFeedforwardGain);

  gfm->setPwmDelayTimeConstant(parameters.appliedPwmDelayTimeConstant());

  gfm->setControllerLimitsPerUnit(
      parameters.minimumFrequencyPu, parameters.maximumFrequencyPu,
      parameters.maximumCurrentReferencePu, parameters.maximumVoltageCommandPu);

  gfm->setFilterParameters(parameters.filterInductance,
                           parameters.filterCapacitance,
                           parameters.filterResistance);

  gfm->withControl(true);

  SPDLOG_INFO("Configured GFM_Siemens nonzero operating point:"
              "\n  P_ref={} pu, Q_ref={} pu, V_ref={} pu, f_ref=1 pu",
              initialActivePowerPu, initialReactivePowerPu, initialVoltagePu);
}

// =============================================================================
// EMT load-step simulation
// =============================================================================

void runEmt(const Parameters &p, const PowerFlowResult &powerFlow) {
  const String simulationName = p.simulation.name + "_EMT";
  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto nGfm = SimNode<Real>::make("BUS_GFM", PhaseType::ABC);
  auto nAfterGridR = SimNode<Real>::make("BUS_AFTER_GRID_R", PhaseType::ABC);
  auto nLoadBus = SimNode<Real>::make("BUS_LOAD", PhaseType::ABC);
  auto nAfterBaseR = SimNode<Real>::make("BUS_AFTER_BASE_R", PhaseType::ABC);
  auto nStepBreakerLoad =
      SimNode<Real>::make("BUS_STEP_BREAKER_LOAD", PhaseType::ABC);
  auto nAfterStepR = SimNode<Real>::make("BUS_AFTER_STEP_R", PhaseType::ABC);

  auto gfm = EMT::Ph3::GFM_Siemens::make("GFM", "GFM", Logger::Level::debug);
  configureGfm(gfm, p.controller, powerFlow);

  auto gridResistor = EMT::Ph3::Resistor::make("GRID_R", Logger::Level::debug);
  gridResistor->setParameters(Math::singlePhaseParameterToThreePhase(
      p.network.gridResistance(p.controller)));

  auto gridInductor = EMT::Ph3::Inductor::make("GRID_L", Logger::Level::debug);
  gridInductor->setParameters(Math::singlePhaseParameterToThreePhase(
      p.network.gridInductance(p.controller)));

  auto baseLoadResistor =
      EMT::Ph3::Resistor::make("BASE_LOAD_R", Logger::Level::debug);
  baseLoadResistor->setParameters(Math::singlePhaseParameterToThreePhase(
      p.network.baseLoadResistance(p.controller)));

  auto baseLoadInductor =
      EMT::Ph3::Inductor::make("BASE_LOAD_L", Logger::Level::debug);
  baseLoadInductor->setParameters(Math::singlePhaseParameterToThreePhase(
      p.network.baseLoadInductance(p.controller)));

  auto stepBreaker =
      EMT::Ph3::Switch::make("STEP_BREAKER", Logger::Level::debug);
  stepBreaker->setParameters(
      Math::singlePhaseParameterToThreePhase(p.network.breakerOpenResistance),
      Math::singlePhaseParameterToThreePhase(
          p.network.breakerClosedResistance(p.controller)));
  stepBreaker->openSwitch();

  auto stepLoadResistor =
      EMT::Ph3::Resistor::make("STEP_LOAD_R", Logger::Level::debug);
  stepLoadResistor->setParameters(Math::singlePhaseParameterToThreePhase(
      p.network.stepLoadResistance(p.controller)));

  auto stepLoadInductor =
      EMT::Ph3::Inductor::make("STEP_LOAD_L", Logger::Level::debug);
  stepLoadInductor->setParameters(Math::singlePhaseParameterToThreePhase(
      p.network.stepLoadInductance(p.controller)));

  gfm->connect({nGfm});
  gridResistor->connect({nGfm, nAfterGridR});
  gridInductor->connect({nAfterGridR, nLoadBus});

  baseLoadResistor->connect({nLoadBus, nAfterBaseR});
  baseLoadInductor->connect({nAfterBaseR, SimNode<Real>::GND});

  stepBreaker->connect({nLoadBus, nStepBreakerLoad});
  stepLoadResistor->connect({nStepBreakerLoad, nAfterStepR});
  stepLoadInductor->connect({nAfterStepR, SimNode<Real>::GND});

  SystemTopology system(p.simulation.frequency,
                        SystemNodeList{nGfm, nAfterGridR, nLoadBus, nAfterBaseR,
                                       nStepBreakerLoad, nAfterStepR},
                        SystemComponentList{gfm, gridResistor, gridInductor,
                                            baseLoadResistor, baseLoadInductor,
                                            stepBreaker, stepLoadResistor,
                                            stepLoadInductor});

  system.initWithPowerflow(powerFlow.system, Domain::EMT);
  gfm->terminal(0)->setPower(powerFlow.gfmPower);

  auto logger =
      DataLogger::make(simulationName, true, p.simulation.logDownsampling);

  // External network.
  logger->logAttribute("BUS_GFM.v", nGfm->attribute("v"));
  logger->logAttribute("BUS_LOAD.v", nLoadBus->attribute("v"));
  logger->logAttribute("BUS_STEP_BREAKER_LOAD.v",
                       nStepBreakerLoad->attribute("v"));
  logger->logAttribute("GRID_R.i", gridResistor->attribute("i_intf"));
  logger->logAttribute("GRID_L.i", gridInductor->attribute("i_intf"));
  logger->logAttribute("BASE_LOAD_R.i", baseLoadResistor->attribute("i_intf"));
  logger->logAttribute("BASE_LOAD_L.i", baseLoadInductor->attribute("i_intf"));
  logger->logAttribute("STEP_BREAKER.v", stepBreaker->attribute("v_intf"));
  logger->logAttribute("STEP_BREAKER.i", stepBreaker->attribute("i_intf"));
  logger->logAttribute("STEP_LOAD_R.i", stepLoadResistor->attribute("i_intf"));
  logger->logAttribute("STEP_LOAD_L.i", stepLoadInductor->attribute("i_intf"));

  // Outer-loop quantities.
  logger->logAttribute("GFM.P_elec_pu", gfm->attribute("P_elec_pu"));
  logger->logAttribute("GFM.Q_elec_pu", gfm->attribute("Q_elec_pu"));
  logger->logAttribute("GFM.P_filtered_pu", gfm->attribute("P_filtered_pu"));
  logger->logAttribute("GFM.Q_filtered_pu", gfm->attribute("Q_filtered_pu"));
  logger->logAttribute("GFM.frequency_pu", gfm->attribute("frequency_pu"));
  logger->logAttribute("GFM.V_magnitude_pu", gfm->attribute("V_magnitude_pu"));
  logger->logAttribute("GFM.f_droop_pu", gfm->attribute("f_droop_pu"));
  logger->logAttribute("GFM.V_droop_pu", gfm->attribute("V_droop_pu"));
  logger->logAttribute("GFM.theta", gfm->attribute("theta"));

  // Cascaded dq controllers.
  logger->logAttribute("GFM.v_pcc_dq_pu", gfm->attribute("v_pcc_dq_pu"));
  logger->logAttribute("GFM.i_pcc_dq_pu", gfm->attribute("i_pcc_dq_pu"));
  logger->logAttribute("GFM.i_lf_dq_pu", gfm->attribute("i_lf_dq_pu"));
  logger->logAttribute("GFM.v_ref_dq_pu", gfm->attribute("v_ref_dq_pu"));
  logger->logAttribute("GFM.i_ref_dq_pu", gfm->attribute("i_ref_dq_pu"));
  logger->logAttribute("GFM.xi_v_dq_pu", gfm->attribute("xi_v_dq_pu"));
  logger->logAttribute("GFM.xi_i_dq_pu", gfm->attribute("xi_i_dq_pu"));
  logger->logAttribute("GFM.v_cmd_pre_pwm_dq_pu",
                       gfm->attribute("v_cmd_pre_pwm_dq_pu"));
  logger->logAttribute("GFM.v_cmd_dq_pu", gfm->attribute("v_cmd_dq_pu"));

  // Internal converter/filter quantities.
  logger->logAttribute("GFM.i_pcc", gfm->attribute("i_pcc"));
  logger->logAttribute("GFM.i_lf", gfm->attribute("i_lf"));
  logger->logAttribute("GFM.i_cf", gfm->attribute("i_cf"));
  logger->logAttribute("GFM.Vsref_pu", gfm->attribute("Vsref_pu"));
  logger->logAttribute("GFM.Vsref", gfm->attribute("Vsref"));

  const Real alignedLoadStepTime =
      std::round(p.simulation.loadStepTime / p.simulation.timeStep) *
      p.simulation.timeStep;

  auto loadStepEvent =
      SwitchEvent3Ph::make(alignedLoadStepTime, stepBreaker, true);

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(p.simulation.timeStep);
  simulation.setFinalTime(p.simulation.finalTime);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doSystemMatrixRecomputation(p.simulation.recomputeSystemMatrix);
  simulation.addEvent(loadStepEvent);
  simulation.addLogger(logger);

  SPDLOG_INFO("Starting EMT run at P_ref={} pu. STEP_BREAKER closes at t={} s.",
              powerFlow.gfmPower.real() / p.controller.ratedApparentPower,
              alignedLoadStepTime);

  simulation.run();

  SPDLOG_INFO("Simulation completed. Results: logs/{}/{}.csv", simulationName,
              simulationName);
}

} // namespace GfmSiemensInitialLoadStep

int main(int argc, char *argv[]) {
  try {
    GfmSiemensInitialLoadStep::Parameters parameters;
    CommandLineArgs args(argc, argv);

    if (argc > 1) {
      parameters.simulation.timeStep = args.timeStep;
      parameters.simulation.finalTime = args.duration;

      if (args.name != "dpsim")
        parameters.simulation.name = args.name;

      if (args.options.find("LOAD_STEP_TIME") != args.options.end()) {
        parameters.simulation.loadStepTime =
            args.getOptionReal("LOAD_STEP_TIME");
      }

      if (args.options.find("INITIAL_P_PU") != args.options.end()) {
        parameters.network.initialGfmActivePowerPu =
            args.getOptionReal("INITIAL_P_PU");
      }

      if (args.options.find("INITIAL_Q_PU") != args.options.end()) {
        parameters.network.initialGfmReactivePowerPu =
            args.getOptionReal("INITIAL_Q_PU");
      }

      if (args.options.find("STEP_P_PU") != args.options.end()) {
        parameters.network.stepLoadActivePowerPu =
            args.getOptionReal("STEP_P_PU");
      }

      if (args.options.find("STEP_Q_PU") != args.options.end()) {
        parameters.network.stepLoadReactivePowerPu =
            args.getOptionReal("STEP_Q_PU");
      }

      if (args.options.find("USE_PWM_DELAY") != args.options.end()) {
        parameters.controller.usePwmDelay = args.getOptionBool("USE_PWM_DELAY");
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

    GfmSiemensInitialLoadStep::validateParameters(parameters);
    GfmSiemensInitialLoadStep::logParameterSummary(parameters);

    const auto powerFlow =
        GfmSiemensInitialLoadStep::buildAndRunPowerFlow(parameters);

    GfmSiemensInitialLoadStep::runEmt(parameters, powerFlow);
    return EXIT_SUCCESS;
  } catch (const std::exception &exception) {
    SPDLOG_ERROR("GFM_Siemens nonzero-power load-step example failed:\n{}",
                 exception.what());
    return EXIT_FAILURE;
  }
}
