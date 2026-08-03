// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <cmath>
#include <filesystem>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

namespace TwoGfmBreakerScenario {

// =============================================================================
// Parameters
// =============================================================================

struct SimulationParameters {
  String name = "EMT_Two_GFM_Droop_Breaker_PF_Init";
  Real frequency = 50.0;
  Real timeStep = 100e-6;
  Real finalTime = 10.0;

  // The steady-state initialization uses the breaker open, matching the EMT
  // initial topology. The breaker then closes during the EMT simulation.
  Real breakerCloseTime = 2.00;

  Bool recomputeSystemMatrix = true;
};

struct GfmParameters {
  Real ratedPower = 50e6;     // VA
  Real ratedVoltage = 10.5e3; // V RMS line-to-line

  // The working 20 kVA / 381 V GFM_Droop filter is transferred to each
  // machine base by preserving its per-unit R, X_L, B_C, and Rd.
  //
  // Original working values:
  //   S_base = 20 kVA
  //   V_LL   = sqrt(3) * 220 V
  //   Lf     = 3 mH
  //   Cf     = 20 uF
  //   Rf     = 0.05 Ohm
  //   Rd     = 4.0 Ohm
  Real sourceRatedPower = 20e3;
  Real sourceRatedVoltage = std::sqrt(3.0) * 220.0;
  Real sourceFilterInductance = 3.0e-3;
  Real sourceFilterCapacitance = 20.0e-6;
  Real sourceFilterResistance = 0.05;
  Real sourceCapacitorDampingResistance = 4.0;

  Real activePowerDroopFraction = 0.05;
  Real reactivePowerDroopFraction = 0.05;
  Real voltageIntegralGain = 40.0; // 1/s

  Real powerMeasurementFilterCutoffFrequency = 10.0; // Hz

  Real minimumFrequencyFraction = 0.90;
  Real maximumFrequencyFraction = 1.10;
  Real minimumVoltageFraction = 0.70;
  Real maximumVoltageFraction = 1.30;

  Real baseImpedance() const {
    return ratedVoltage * ratedVoltage / ratedPower;
  }

  Real sourceBaseImpedance() const {
    return sourceRatedVoltage * sourceRatedVoltage / sourceRatedPower;
  }

  Real filterInductance(Real frequency) const {
    const Real omega = 2.0 * PI * frequency;
    const Real sourceReactancePu =
        omega * sourceFilterInductance / sourceBaseImpedance();
    return sourceReactancePu * baseImpedance() / omega;
  }

  Real filterCapacitance(Real frequency) const {
    const Real omega = 2.0 * PI * frequency;
    const Real sourceSusceptancePu =
        omega * sourceFilterCapacitance * sourceBaseImpedance();
    return sourceSusceptancePu / (omega * baseImpedance());
  }

  Real filterResistance() const {
    const Real resistancePu = sourceFilterResistance / sourceBaseImpedance();
    return resistancePu * baseImpedance();
  }

  Real capacitorDampingResistance() const {
    const Real resistancePu =
        sourceCapacitorDampingResistance / sourceBaseImpedance();
    return resistancePu * baseImpedance();
  }

  Real activePowerDroop(Real frequency) const {
    return activePowerDroopFraction * frequency / ratedPower;
  }

  Real reactivePowerDroop(Real voltagePeakPhase) const {
    return reactivePowerDroopFraction * voltagePeakPhase / ratedPower;
  }

  Real powerFilterTimeConstant() const {
    return 1.0 / (2.0 * PI * powerMeasurementFilterCutoffFrequency);
  }
};

struct GasTransformerParameters {
  Real nominalVoltageLow = 10.5e3;
  Real nominalVoltageHigh = 220e3;
  Real ratedPower = 50e6;
  Real ratioMagnitude = nominalVoltageLow / nominalVoltageHigh;
  Real ratioPhase = 0.0;

  Real resistance() const {
    const Real baseImpedance =
        nominalVoltageHigh * nominalVoltageHigh / ratedPower;
    return 0.00376196 * baseImpedance;
  }

  Real inductance(Real frequency) const {
    const Real baseImpedance =
        nominalVoltageHigh * nominalVoltageHigh / ratedPower;
    return 0.1007298 * baseImpedance / (2.0 * PI * frequency);
  }
};

struct PshTransformerParameters {
  Real nominalVoltageLow = 18.0e3;
  Real nominalVoltageHigh = 220e3;
  Real ratedPower = 200e6;
  Real ratioMagnitude = nominalVoltageLow / nominalVoltageHigh;
  Real ratioPhase = 0.0;

  // Values from the original scenario:
  //   R = 0.41745 * 2
  //   L = 0.0481260775594524 * 2
  Real resistance() const { return 0.41745 * 2.0; }

  Real inductance(Real /*frequency*/) const { return 0.0481260775594524 * 2.0; }
};

struct Line1Parameters {
  // Parameters supplied for the cable / LINE_1 branch.
  Real baseVoltage = 220e3;
  Real resistance = 5.04669;
  Real inductance = 0.13523123641;
  Real capacitance = 1.93522865e-6;
  Real conductance = 1e-15;
};

struct Line3Parameters {
  Real baseVoltage = 220e3;
  Real length = 5.0;

  Real resistance() const { return 0.0749 * length; }
  Real inductance() const { return 1.270693e-3 * length; }
  Real capacitance() const { return 0.00466961e-6 * length; }

  Real conductance = 1e-15;
};

struct BreakerParameters {
  // Original scenario values.
  Real openResistance = 1e12;   // Ohm
  Real closedResistance = 1e-3; // Ohm
};

struct Parameters {
  SimulationParameters simulation;

  GfmParameters gasGfm;
  GfmParameters pshaGfm;

  GasTransformerParameters gasTransformer;
  PshTransformerParameters pshaTransformer;

  Line1Parameters line1;
  Line3Parameters line3;
  BreakerParameters breaker;

  // Both SP voltage sources start at the same phase angle. Change this value
  // deliberately to study out-of-phase breaker closing.
  Real pshaInitialAngleDegrees = 0.0;

  Parameters() {
    // GEN_gas data from the original scenario.
    gasGfm.ratedPower = 50e6;
    gasGfm.ratedVoltage = 10.5e3;

    // GEN_psh data from the original scenario.
    pshaGfm.ratedPower = 200e6;
    pshaGfm.ratedVoltage = 18.0e3;
  }
};

// =============================================================================
// Power-flow / steady-state initialization
// =============================================================================

struct PowerFlowResult {
  SystemTopology system;

  Complex gasPower;
  Complex pshaPower;

  Complex gasBusVoltage;
  Complex pshaBusVoltage;

  Complex breakerGridSideVoltage;
  Complex breakerPshaSideVoltage;
};

Complex
generatorPositivePower(const std::shared_ptr<SP::Ph1::NetworkInjection> &source,
                       const SimNode<Complex>::Ptr &sourceBus) {
  const Complex voltage = sourceBus->singleVoltage();

  // NetworkInjection follows the component-consumer convention:
  // positive interface current enters the voltage source from the network.
  // Generator-positive injected power therefore has the opposite sign.
  const Complex consumerCurrent = (**source->mIntfCurrent)(0, 0);
  return -voltage * std::conj(consumerCurrent);
}

PowerFlowResult buildAndRunPowerFlow(const Parameters &p) {
  const String simulationName = p.simulation.name + "_PF";

  SPDLOG_INFO("Scenario bases:"
              "\n  GEN_gas: S_rated={} VA, V_rated={} V_LL RMS"
              "\n  GEN_psh: S_rated={} VA, V_rated={} V_LL RMS"
              "\n  TR_gas : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
              "\n  TR_psh : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
              "\n  breaker: R_open={} Ohm, R_closed={} Ohm",
              p.gasGfm.ratedPower, p.gasGfm.ratedVoltage, p.pshaGfm.ratedPower,
              p.pshaGfm.ratedVoltage, p.gasTransformer.nominalVoltageLow,
              p.gasTransformer.nominalVoltageHigh, p.gasTransformer.ratedPower,
              p.gasTransformer.resistance(),
              p.gasTransformer.inductance(p.simulation.frequency),
              p.pshaTransformer.nominalVoltageLow,
              p.pshaTransformer.nominalVoltageHigh,
              p.pshaTransformer.ratedPower, p.pshaTransformer.resistance(),
              p.pshaTransformer.inductance(p.simulation.frequency),
              p.breaker.openResistance, p.breaker.closedResistance);

  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  // Node names intentionally match the EMT node names so that
  // SystemTopology::initWithPowerflow() can transfer every solved phasor.
  auto busGas = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busB = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busA = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busPshaGrid = SimNode<Complex>::make("BUS_psha_grid", PhaseType::Single);
  auto busPshaHv = SimNode<Complex>::make("BUS_psha_hv", PhaseType::Single);
  auto busPsha = SimNode<Complex>::make("BUS_psha", PhaseType::Single);

  auto gasInjection =
      SP::Ph1::NetworkInjection::make("GEN_gas_PF", Logger::Level::debug);
  gasInjection->setParameters(Math::polar(p.gasGfm.ratedVoltage, 0.0),
                              p.simulation.frequency);

  auto pshaInjection =
      SP::Ph1::NetworkInjection::make("GEN_psha_PF", Logger::Level::debug);
  pshaInjection->setParameters(
      Math::polar(p.pshaGfm.ratedVoltage,
                  p.pshaInitialAngleDegrees * PI / 180.0),
      p.simulation.frequency);

  auto gasTransformer = SP::Ph1::Transformer::make("TR_gas_PF", "TR_gas_PF",
                                                   Logger::Level::debug, true);
  gasTransformer->setParameters(
      p.gasTransformer.nominalVoltageLow, p.gasTransformer.nominalVoltageHigh,
      p.gasTransformer.ratedPower, p.gasTransformer.ratioMagnitude,
      p.gasTransformer.ratioPhase, p.gasTransformer.resistance(),
      p.gasTransformer.inductance(p.simulation.frequency));

  auto line1 = SP::Ph1::PiLine::make("LINE_1_PF", Logger::Level::debug);
  line1->setParameters(p.line1.resistance, p.line1.inductance,
                       p.line1.capacitance, p.line1.conductance);

  auto line3 = SP::Ph1::PiLine::make("LINE_3_PF", Logger::Level::debug);
  line3->setParameters(p.line3.resistance(), p.line3.inductance(),
                       p.line3.capacitance(), p.line3.conductance);

  auto breaker = SP::Ph1::Switch::make("BREAKER_psha_PF", Logger::Level::debug);
  breaker->setParameters(p.breaker.openResistance, p.breaker.closedResistance,
                         false);

  auto pshaTransformer = SP::Ph1::Transformer::make("TR_psha_PF", "TR_psha_PF",
                                                    Logger::Level::debug, true);
  pshaTransformer->setParameters(
      p.pshaTransformer.nominalVoltageLow, p.pshaTransformer.nominalVoltageHigh,
      p.pshaTransformer.ratedPower, p.pshaTransformer.ratioMagnitude,
      p.pshaTransformer.ratioPhase, p.pshaTransformer.resistance(),
      p.pshaTransformer.inductance(p.simulation.frequency));

  gasInjection->connect({busGas});
  gasTransformer->connect({busGas, busB});
  line1->connect({busB, busA});
  line3->connect({busA, busPshaGrid});
  breaker->connect({busPshaGrid, busPshaHv});
  pshaTransformer->connect({busPsha, busPshaHv});
  pshaInjection->connect({busPsha});

  SystemTopology system(p.simulation.frequency,
                        SystemNodeList{
                            busGas,
                            busB,
                            busA,
                            busPshaGrid,
                            busPshaHv,
                            busPsha,
                        },
                        SystemComponentList{
                            gasInjection,
                            gasTransformer,
                            line1,
                            line3,
                            breaker,
                            pshaTransformer,
                            pshaInjection,
                        });

  auto logger = DataLogger::make(simulationName);

  logger->logAttribute("V_BUS_gas_PF", busGas->attribute("v"));
  logger->logAttribute("V_BUS_b_PF", busB->attribute("v"));
  logger->logAttribute("V_BUS_a_PF", busA->attribute("v"));
  logger->logAttribute("V_BUS_psha_grid_PF", busPshaGrid->attribute("v"));
  logger->logAttribute("V_BUS_psha_hv_PF", busPshaHv->attribute("v"));
  logger->logAttribute("V_BUS_psha_PF", busPsha->attribute("v"));

  logger->logAttribute("I_GEN_gas_PF", gasInjection->attribute("i_intf"));
  logger->logAttribute("I_GEN_psha_PF", pshaInjection->attribute("i_intf"));
  logger->logAttribute("I_BREAKER_psha_PF", breaker->attribute("i_intf"));

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(1.0);
  simulation.setFinalTime(1.0);
  simulation.setDomain(Domain::SP);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.addLogger(logger);
  simulation.run();

  const Complex gasPower = generatorPositivePower(gasInjection, busGas);
  const Complex pshaPower = generatorPositivePower(pshaInjection, busPsha);

  const Complex gasVoltage = busGas->singleVoltage();
  const Complex pshaVoltage = busPsha->singleVoltage();
  const Complex breakerGridVoltage = busPshaGrid->singleVoltage();
  const Complex breakerPshaVoltage = busPshaHv->singleVoltage();

  const auto finiteComplex = [](const Complex &value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  };

  if (!(std::abs(gasVoltage) > 0.0) || !(std::abs(pshaVoltage) > 0.0) ||
      !finiteComplex(gasVoltage) || !finiteComplex(pshaVoltage) ||
      !finiteComplex(gasPower) || !finiteComplex(pshaPower) ||
      !finiteComplex(breakerGridVoltage) ||
      !finiteComplex(breakerPshaVoltage)) {
    throw std::runtime_error(
        "Invalid SP steady-state result for the two-GFM breaker case.");
  }

  SPDLOG_INFO(
      "Solved open-breaker SP initialization:"
      "\n  GEN_gas : P={} W, Q={} var, |V|={} V_LL RMS, angle={} deg"
      "\n  GEN_psha: P={} W, Q={} var, |V|={} V_LL RMS, angle={} deg"
      "\n  breaker grid side: |V|={} V_LL RMS, angle={} deg"
      "\n  breaker psha side: |V|={} V_LL RMS, angle={} deg"
      "\n  breaker mismatch: |dV|={} V_LL RMS, dAngle={} deg",
      gasPower.real(), gasPower.imag(), std::abs(gasVoltage),
      std::arg(gasVoltage) * 180.0 / PI, pshaPower.real(), pshaPower.imag(),
      std::abs(pshaVoltage), std::arg(pshaVoltage) * 180.0 / PI,
      std::abs(breakerGridVoltage), std::arg(breakerGridVoltage) * 180.0 / PI,
      std::abs(breakerPshaVoltage), std::arg(breakerPshaVoltage) * 180.0 / PI,
      std::abs(breakerGridVoltage - breakerPshaVoltage),
      std::arg(breakerGridVoltage / breakerPshaVoltage) * 180.0 / PI);

  return {
      system,      gasPower,           pshaPower,          gasVoltage,
      pshaVoltage, breakerGridVoltage, breakerPshaVoltage,
  };
}

// =============================================================================
// EMT model
// =============================================================================

void configureGfm(const std::shared_ptr<EMT::Ph3::GFM_Droop> &gfm,
                  const GfmParameters &parameters, Real frequency,
                  const Complex &initialBusVoltage,
                  const Complex &initialGeneratorPower) {
  const Real voltagePeakPhase = RMS3PH_TO_PEAK1PH * std::abs(initialBusVoltage);

  const Real lf = parameters.filterInductance(frequency);
  const Real cf = parameters.filterCapacitance(frequency);
  const Real rf = parameters.filterResistance();
  const Real rd = parameters.capacitorDampingResistance();

  gfm->setParameters(frequency, voltagePeakPhase, initialGeneratorPower.real(),
                     initialGeneratorPower.imag());

  gfm->setDroopParameters(parameters.activePowerDroop(frequency),
                          parameters.reactivePowerDroop(voltagePeakPhase),
                          parameters.voltageIntegralGain);

  gfm->setPowerFilterTimeConstant(parameters.powerFilterTimeConstant());

  gfm->setControllerLimits(parameters.minimumFrequencyFraction * frequency,
                           parameters.maximumFrequencyFraction * frequency,
                           parameters.minimumVoltageFraction * voltagePeakPhase,
                           parameters.maximumVoltageFraction *
                               voltagePeakPhase);

  gfm->setFilterParameters(lf, cf, rf, rd);
  gfm->withControl(true);

  SPDLOG_INFO("{} parameters:"
              "\n  P_ref={} W, Q_ref={} var"
              "\n  V_ref={} V phase peak"
              "\n  k_p={} Hz/W, k_q={} V/var, k_iv={} 1/s"
              "\n  Lf={} H, Cf={} F, Rf={} Ohm, Rd={} Ohm",
              gfm->name(), initialGeneratorPower.real(),
              initialGeneratorPower.imag(), voltagePeakPhase,
              parameters.activePowerDroop(frequency),
              parameters.reactivePowerDroop(voltagePeakPhase),
              parameters.voltageIntegralGain, lf, cf, rf, rd);
}

void runEmt(const Parameters &p, const PowerFlowResult &powerFlow) {
  const String simulationName = p.simulation.name + "_EMT";

  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busA = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto busPshaGrid = SimNode<Real>::make("BUS_psha_grid", PhaseType::ABC);
  auto busPshaHv = SimNode<Real>::make("BUS_psha_hv", PhaseType::ABC);
  auto busPsha = SimNode<Real>::make("BUS_psha", PhaseType::ABC);

  auto gasGfm = EMT::Ph3::GFM_Droop::make("GEN_gas", "GEN_gas",
                                          Logger::Level::debug, false);
  auto pshaGfm = EMT::Ph3::GFM_Droop::make("GEN_psha", "GEN_psha",
                                           Logger::Level::debug, false);

  configureGfm(gasGfm, p.gasGfm, p.simulation.frequency,
               powerFlow.gasBusVoltage, powerFlow.gasPower);

  configureGfm(pshaGfm, p.pshaGfm, p.simulation.frequency,
               powerFlow.pshaBusVoltage, powerFlow.pshaPower);

  gasGfm->withControl(true);
  pshaGfm->withControl(true);

  auto gasTransformer = EMT::Ph3::Transformer::make("TR_gas", "TR_gas",
                                                    Logger::Level::debug, true);
  gasTransformer->setParameters(
      p.gasTransformer.nominalVoltageLow, p.gasTransformer.nominalVoltageHigh,
      p.gasTransformer.ratedPower, p.gasTransformer.ratioMagnitude,
      p.gasTransformer.ratioPhase,
      Math::singlePhaseParameterToThreePhase(p.gasTransformer.resistance()),
      Math::singlePhaseParameterToThreePhase(
          p.gasTransformer.inductance(p.simulation.frequency)));

  auto line1 = EMT::Ph3::PiLine::make("LINE_1", Logger::Level::debug);
  line1->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line1.resistance),
      Math::singlePhaseParameterToThreePhase(p.line1.inductance),
      Math::singlePhaseParameterToThreePhase(p.line1.capacitance),
      Math::singlePhaseParameterToThreePhase(p.line1.conductance));

  auto line3 = EMT::Ph3::PiLine::make("LINE_3", Logger::Level::debug);
  line3->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line3.resistance()),
      Math::singlePhaseParameterToThreePhase(p.line3.inductance()),
      Math::singlePhaseParameterToThreePhase(p.line3.capacitance()),
      Math::singlePhaseParameterToThreePhase(p.line3.conductance));

  auto breaker = EMT::Ph3::Switch::make("BREAKER_psha", Logger::Level::debug);
  breaker->setParameters(
      Math::singlePhaseParameterToThreePhase(p.breaker.openResistance),
      Math::singlePhaseParameterToThreePhase(p.breaker.closedResistance));
  breaker->openSwitch();

  auto pshaTransformer = EMT::Ph3::Transformer::make(
      "TR_psha", "TR_psha", Logger::Level::debug, true);
  pshaTransformer->setParameters(
      p.pshaTransformer.nominalVoltageLow, p.pshaTransformer.nominalVoltageHigh,
      p.pshaTransformer.ratedPower, p.pshaTransformer.ratioMagnitude,
      p.pshaTransformer.ratioPhase,
      Math::singlePhaseParameterToThreePhase(p.pshaTransformer.resistance()),
      Math::singlePhaseParameterToThreePhase(
          p.pshaTransformer.inductance(p.simulation.frequency)));

  // GFM_Droop is a one-terminal composite whose external terminal is its PCC.
  gasGfm->connect({busGas});
  gasTransformer->connect({busGas, busB});
  line1->connect({busB, busA});
  line3->connect({busA, busPshaGrid});
  breaker->connect({busPshaGrid, busPshaHv});
  pshaTransformer->connect({busPsha, busPshaHv});
  pshaGfm->connect({busPsha});

  SystemTopology system(p.simulation.frequency,
                        SystemNodeList{
                            busGas,
                            busB,
                            busA,
                            busPshaGrid,
                            busPshaHv,
                            busPsha,
                        },
                        SystemComponentList{
                            gasGfm,
                            gasTransformer,
                            line1,
                            line3,
                            breaker,
                            pshaTransformer,
                            pshaGfm,
                        });

  // Every EMT node has an identically named SP node. The open-breaker
  // steady-state phasors therefore initialize both electrical islands,
  // including both sides of the breaker.
  system.initWithPowerflow(powerFlow.system, Domain::EMT);

  // GFM_Droop is not a synchronous-generator PF component, so its solved
  // generator-positive terminal power must be assigned explicitly.
  gasGfm->terminal(0)->setPower(powerFlow.gasPower);
  pshaGfm->terminal(0)->setPower(powerFlow.pshaPower);

  auto logger = DataLogger::make(simulationName);

  logger->logAttribute("BUS_gas_v", busGas->attribute("v"));
  logger->logAttribute("BUS_b_v", busB->attribute("v"));
  logger->logAttribute("BUS_a_v", busA->attribute("v"));
  logger->logAttribute("BUS_psha_grid_v", busPshaGrid->attribute("v"));
  logger->logAttribute("BUS_psha_hv_v", busPshaHv->attribute("v"));
  logger->logAttribute("BUS_psha_v", busPsha->attribute("v"));

  logger->logAttribute("TR_gas_i", gasTransformer->attribute("i_intf"));
  logger->logAttribute("LINE_1_i", line1->attribute("i_intf"));
  logger->logAttribute("LINE_3_i", line3->attribute("i_intf"));
  logger->logAttribute("BREAKER_psha_i", breaker->attribute("i_intf"));
  logger->logAttribute("BREAKER_psha_v", breaker->attribute("v_intf"));
  logger->logAttribute("TR_psha_i", pshaTransformer->attribute("i_intf"));

  logger->logAttribute("GEN_gas_P_elec", gasGfm->attribute("P_elec"));
  logger->logAttribute("GEN_gas_Q_elec", gasGfm->attribute("Q_elec"));
  logger->logAttribute("GEN_gas_P_filtered", gasGfm->attribute("P_filtered"));
  logger->logAttribute("GEN_gas_Q_filtered", gasGfm->attribute("Q_filtered"));
  logger->logAttribute("GEN_gas_frequency", gasGfm->attribute("frequency"));
  logger->logAttribute("GEN_gas_theta", gasGfm->attribute("theta"));
  logger->logAttribute("GEN_gas_voltage_magnitude",
                       gasGfm->attribute("V_magnitude"));
  logger->logAttribute("GEN_gas_i_pcc", gasGfm->attribute("i_pcc"));
  logger->logAttribute("GEN_gas_v_source", gasGfm->attribute("Vs"));

  logger->logAttribute("GEN_psha_P_elec", pshaGfm->attribute("P_elec"));
  logger->logAttribute("GEN_psha_Q_elec", pshaGfm->attribute("Q_elec"));
  logger->logAttribute("GEN_psha_P_filtered", pshaGfm->attribute("P_filtered"));
  logger->logAttribute("GEN_psha_Q_filtered", pshaGfm->attribute("Q_filtered"));
  logger->logAttribute("GEN_psha_frequency", pshaGfm->attribute("frequency"));
  logger->logAttribute("GEN_psha_theta", pshaGfm->attribute("theta"));
  logger->logAttribute("GEN_psha_voltage_magnitude",
                       pshaGfm->attribute("V_magnitude"));
  logger->logAttribute("GEN_psha_i_pcc", pshaGfm->attribute("i_pcc"));
  logger->logAttribute("GEN_psha_v_source", pshaGfm->attribute("Vs"));

  Simulation simulation(simulationName, Logger::Level::info);

  simulation.addEvent(
      SwitchEvent3Ph::make(p.simulation.breakerCloseTime, breaker, true));

  simulation.setSystem(system);
  simulation.setTimeStep(p.simulation.timeStep);
  simulation.setFinalTime(p.simulation.finalTime);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doSystemMatrixRecomputation(p.simulation.recomputeSystemMatrix);
  simulation.addLogger(logger);
  simulation.run();
}

} // namespace TwoGfmBreakerScenario

int main(int argc, char *argv[]) {
  TwoGfmBreakerScenario::Parameters parameters;

  CommandLineArgs args(argc, argv);

  if (argc > 1) {
    parameters.simulation.timeStep = args.timeStep;
    parameters.simulation.finalTime = args.duration;

    if (args.name != "dpsim")
      parameters.simulation.name = args.name;

    if (args.options.find("BREAKER_CLOSE_TIME") != args.options.end()) {
      parameters.simulation.breakerCloseTime =
          args.getOptionReal("BREAKER_CLOSE_TIME");
    }

    if (args.options.find("PSHA_INITIAL_ANGLE_DEG") != args.options.end()) {
      parameters.pshaInitialAngleDegrees =
          args.getOptionReal("PSHA_INITIAL_ANGLE_DEG");
    }
  }

  if (!(parameters.simulation.timeStep > 0.0) ||
      !(parameters.simulation.finalTime > 0.0) ||
      !(parameters.simulation.breakerCloseTime >= 0.0) ||
      !(parameters.simulation.breakerCloseTime <
        parameters.simulation.finalTime)) {
    throw std::invalid_argument("Require dt>0, finalTime>0, and "
                                "0<=breakerCloseTime<finalTime.");
  }

  const auto powerFlow =
      TwoGfmBreakerScenario::buildAndRunPowerFlow(parameters);

  TwoGfmBreakerScenario::runEmt(parameters, powerFlow);

  return 0;
}
