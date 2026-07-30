// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"
#include "../GeneratorFactory.h"

#include <DPsim.h>

#include <cmath>
#include <filesystem>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

namespace ScenarioA {

// -----------------------------------------------------------------------------
// Simulation and component parameters
// -----------------------------------------------------------------------------

struct SimulationParameters {
  String name = "EMT_Scenario_A_SSN_GFM_SwingEquation";
  Real frequency = 50.0;
  Real timeStep = 100e-6;
  Real finalTime = 60.0;

  Bool recomputeSystemMatrix = true;
  UInt linearizationUpdateInterval = 1;
};

struct GridFormingInverterParameters {
  // Electrical base inherited from the synchronous generator being replaced.
  Real ratedPower = 50e6;
  Real ratedVoltage = 10.5e3;
  Real inertiaConstantH = 5.0;

  // Active-power/frequency droop used to derive the damping coefficient.
  //
  // R = 0.05 means a 5% frequency deviation for a 1 pu active-power change
  // in the corresponding steady-state droop relation.
  Real activePowerFrequencyDroop = 0.05;

  // LCL-filter parameters from the supplied controller design.
  Real Lf = 100e-6;
  Real Cf = 2.5e-3;
  Real Rf = 2.07e-3;
  Real Rc = 1e-5;
  Real converterDelayTimeConstant = 0.5e-3;

  // Inner-controller design specifications.
  Real currentControllerBandwidth = 2000.0;
  Real voltageControllerNaturalFrequency = 387.1;
  Real voltageControllerDamping = 0.864;

  // The supplied controller design does not define a separate outer Q/V loop.
  // Keeping both gains at zero holds the initialized internal voltage
  // magnitude while the cascaded voltage/current PI loops remain active.
  Real voltageDroopGain = 0.0;
  Real reactiveIntegralGain = 0.0;

  // Optional detailed GFM features.
  Real activeDampingGain = 0.0;
  Real gridCurrentFeedforward = 0.0;
  Real virtualResistance = 0.0;
  Real virtualReactance = 0.0;

  Real omegaNominal(Real systemFrequency) const {
    return 2.0 * PI * systemFrequency;
  }

  Real kpCurrent() const { return Lf * currentControllerBandwidth; }

  Real kiCurrent() const { return Rf * currentControllerBandwidth; }

  Real kpVoltage() const {
    return 2.0 * voltageControllerDamping * voltageControllerNaturalFrequency *
           Cf;
  }

  Real kiVoltage() const {
    return voltageControllerNaturalFrequency *
           voltageControllerNaturalFrequency * Cf;
  }

  Real delayBandwidth() const {
    if (converterDelayTimeConstant <= 0.0)
      throw std::invalid_argument(
          "The converter delay time constant must be positive.");

    return 1.0 / converterDelayTimeConstant;
  }

  Real powerFilterCutoff(Real systemFrequency) const {
    // Equivalent to OmegaCutoff = OmegaNull in the supplied design.
    return omegaNominal(systemFrequency);
  }

  Real virtualInertia(Real systemFrequency) const {
    const Real omegaN = omegaNominal(systemFrequency);

    // SSN_GFM uses:
    //
    //   J * domega/dt = (Pref - P)/omega - D*(omega - omegaN)
    //
    // Conversion from the synchronous-machine inertia constant H:
    //
    //   J = 2*H*Sbase/omegaN^2
    return 2.0 * inertiaConstantH * ratedPower / (omegaN * omegaN);
  }

  Real dampingCoefficient(Real systemFrequency) const {
    if (activePowerFrequencyDroop <= 0.0)
      throw std::invalid_argument(
          "The active-power/frequency droop must be positive.");

    const Real omegaN = omegaNominal(systemFrequency);

    // Matching the steady-state swing-equation relation to a pu droop R:
    //
    //   D = Sbase / (R*omegaN^2)
    return ratedPower / (activePowerFrequencyDroop * omegaN * omegaN);
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

    // Preserve the factor used in the original scenario.
    return 0.00376196 * baseImpedance;
  }

  Real inductance(Real systemFrequency) const {
    const Real baseImpedance =
        nominalVoltageHigh * nominalVoltageHigh / ratedPower;

    // Preserve the factor used in the original scenario.
    return 0.1007298 * baseImpedance / (2.0 * PI * systemFrequency);
  }
};

struct CableParameters {
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

struct FaultParameters {
  // Retained from the original command-line interface. The original scenario
  // did not connect or schedule the fault switch, so this file does not claim
  // to apply a fault.
  Bool startFaultEvent = true;
  Bool endFaultEvent = true;
  Bool useVariableResistanceSwitch = false;

  Real startTime = 100.0;
  Real endTime = 11100.0;

  Real switchOpenResistance = 1e12;
  Real switchClosedResistance = 1e-8;
};

struct Parameters {
  SimulationParameters simulation;
  GridFormingInverterParameters gfm;
  GasTransformerParameters gasTransformer;
  CableParameters cable;
  Line3Parameters line3;
  FaultParameters fault;
};

// -----------------------------------------------------------------------------
// Power-flow result needed by the EMT GFM initialization
// -----------------------------------------------------------------------------

struct PowerFlowResult {
  SystemTopology system;

  Complex generatorPower;
  Complex generatorBusVoltage;
};

// -----------------------------------------------------------------------------
// Power-flow initialization
// -----------------------------------------------------------------------------

PowerFlowResult buildAndRunPowerFlow(const Parameters &parameters) {
  const auto &simulation = parameters.simulation;
  const auto &gfm = parameters.gfm;
  const auto &transformer = parameters.gasTransformer;
  const auto &cable = parameters.cable;
  const auto &line3 = parameters.line3;

  const String simulationName = simulation.name + "_PF";
  const std::filesystem::path logDirectory =
      std::filesystem::path("logs") / simulationName;

  std::filesystem::create_directories(logDirectory);
  Logger::setLogDir(logDirectory.string());

  auto busGas = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busB = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busA = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busPsha = SimNode<Complex>::make("BUS_psha", PhaseType::Single);

  auto generator =
      SP::Ph1::SynchronGenerator::make("GEN_gas", Logger::Level::off);

  generator->setParameters(gfm.ratedPower, gfm.ratedVoltage, 0.0,
                           gfm.ratedVoltage, PowerflowBusType::VD);

  generator->setBaseVoltage(gfm.ratedVoltage);

  auto gasTransformer =
      SP::Ph1::Transformer::make("TR_gas", Logger::Level::off);

  gasTransformer->setParameters(
      transformer.nominalVoltageLow, transformer.nominalVoltageHigh,
      transformer.ratedPower, transformer.ratioMagnitude,
      transformer.ratioPhase, transformer.resistance(),
      transformer.inductance(simulation.frequency));

  gasTransformer->setBaseVoltage(transformer.nominalVoltageHigh);

  auto cableComponent = SP::Ph1::PiLine::make("cable", Logger::Level::off);

  cableComponent->setParameters(cable.resistance, cable.inductance,
                                cable.capacitance, cable.conductance);

  cableComponent->setBaseVoltage(cable.baseVoltage);

  auto line3Component = SP::Ph1::PiLine::make("line_3", Logger::Level::off);

  line3Component->setParameters(line3.resistance(), line3.inductance(),
                                line3.capacitance(), line3.conductance);

  line3Component->setBaseVoltage(line3.baseVoltage);

  generator->connect({busGas});
  gasTransformer->connect({busGas, busB});
  cableComponent->connect({busB, busA});
  line3Component->connect({busA, busPsha});

  SystemTopology system(simulation.frequency,
                        SystemNodeList{busGas, busB, busA, busPsha},
                        SystemComponentList{
                            generator,
                            gasTransformer,
                            cableComponent,
                            line3Component,
                        });

  auto logger = DataLogger::make(simulationName, Logger::Level::off);

  logger->logAttribute("V_BUS_gas_PF", busGas->attribute("v"));
  logger->logAttribute("V_BUS_b_PF", busB->attribute("v"));
  logger->logAttribute("V_BUS_a_PF", busA->attribute("v"));
  logger->logAttribute("V_BUS_psha_PF", busPsha->attribute("v"));

  Simulation powerFlowSimulation(simulationName, Logger::Level::info);

  powerFlowSimulation.setSystem(system);
  powerFlowSimulation.setTimeStep(simulation.timeStep);
  powerFlowSimulation.setFinalTime(simulation.timeStep);
  powerFlowSimulation.setDomain(Domain::SP);
  powerFlowSimulation.setSolverType(Solver::Type::NRP);
  powerFlowSimulation.setSolverAndComponentBehaviour(
      Solver::Behaviour::Initialization);
  powerFlowSimulation.doInitFromNodesAndTerminals(false);
  powerFlowSimulation.addLogger(logger);
  powerFlowSimulation.run();

  const Complex solvedPower = generator->getApparentPower();
  const Complex solvedVoltage = busGas->singleVoltage();

  auto log = Logger::get(simulationName + "_result", Logger::Level::info,
                         Logger::Level::info);

  log->info("Solved GEN_gas operating point: "
            "P={:.9e} W, Q={:.9e} var, "
            "|V_BUS_gas|={:.9e} V RMS line-to-line, angle={:.6f} deg",
            solvedPower.real(), solvedPower.imag(), std::abs(solvedVoltage),
            std::arg(solvedVoltage) * 180.0 / PI);

  return {
      system,
      solvedPower,
      solvedVoltage,
  };
}

// -----------------------------------------------------------------------------
// EMT scenario with swing-equation SSN GFM
// -----------------------------------------------------------------------------

void runScenario(const Parameters &parameters, Real inertiaScale,
                 Real dampingScale) {
  const auto &simulation = parameters.simulation;
  const auto &gfm = parameters.gfm;
  const auto &transformer = parameters.gasTransformer;
  const auto &cable = parameters.cable;
  const auto &line3 = parameters.line3;

  if (inertiaScale <= 0.0)
    throw std::invalid_argument("SCALEINERTIA_G1 must be positive.");

  if (dampingScale < 0.0)
    throw std::invalid_argument("SCALEDAMPING_G1 must be non-negative.");

  const PowerFlowResult powerFlow = buildAndRunPowerFlow(parameters);

  const String simulationName = simulation.name + "_EMT";
  const std::filesystem::path logDirectory =
      std::filesystem::path("logs") / simulationName;

  std::filesystem::create_directories(logDirectory);
  Logger::setLogDir(logDirectory.string());

  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busA = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto busPsha = SimNode<Real>::make("BUS_psha", PhaseType::ABC);

  const Real omegaNominal = gfm.omegaNominal(simulation.frequency);

  // SSN_GFM uses an amplitude-invariant Park transform. Its nominal voltage
  // is therefore the phase-to-neutral peak value, derived from the solved
  // line-to-line RMS power-flow voltage.
  const Real nominalVoltagePeakPhase =
      RMS3PH_TO_PEAK1PH * std::abs(powerFlow.generatorBusVoltage);

  const Real virtualInertia =
      gfm.virtualInertia(simulation.frequency) * inertiaScale;

  const Real dampingCoefficient =
      gfm.dampingCoefficient(simulation.frequency) * dampingScale;

  auto generator =
      EMT::Ph3::SSN_GFM::make("GEN_gas", "GEN_gas", Logger::Level::off);

  generator->setParameters(
      gfm.Lf, gfm.Cf, gfm.Rf, gfm.Rc, nominalVoltagePeakPhase, omegaNominal,
      powerFlow.generatorPower.real(), powerFlow.generatorPower.imag(),
      virtualInertia, dampingCoefficient, gfm.voltageDroopGain,
      gfm.reactiveIntegralGain, gfm.kpVoltage(), gfm.kiVoltage(),
      gfm.kpCurrent(), gfm.kiCurrent(), gfm.activeDampingGain,
      gfm.powerFilterCutoff(simulation.frequency), gfm.delayBandwidth());

  generator->setGridCurrentFeedforward(gfm.gridCurrentFeedforward);

  generator->setVirtualImpedance(gfm.virtualResistance, gfm.virtualReactance);

  generator->setLinearizationUpdateInterval(
      simulation.linearizationUpdateInterval);

  auto gasTransformer =
      EMT::Ph3::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);

  gasTransformer->setParameters(
      transformer.nominalVoltageLow, transformer.nominalVoltageHigh,
      transformer.ratedPower, transformer.ratioMagnitude,
      transformer.ratioPhase,
      Math::singlePhaseParameterToThreePhase(transformer.resistance()),
      Math::singlePhaseParameterToThreePhase(
          transformer.inductance(simulation.frequency)));

  auto cableComponent = EMT::Ph3::PiLine::make("cable", Logger::Level::off);

  cableComponent->setParameters(
      Math::singlePhaseParameterToThreePhase(cable.resistance),
      Math::singlePhaseParameterToThreePhase(cable.inductance),
      Math::singlePhaseParameterToThreePhase(cable.capacitance),
      Math::singlePhaseParameterToThreePhase(cable.conductance));

  auto line3Component = EMT::Ph3::PiLine::make("line_3", Logger::Level::off);

  line3Component->setParameters(
      Math::singlePhaseParameterToThreePhase(line3.resistance()),
      Math::singlePhaseParameterToThreePhase(line3.inductance()),
      Math::singlePhaseParameterToThreePhase(line3.capacitance()),
      Math::singlePhaseParameterToThreePhase(line3.conductance));

  // SSN components use terminal 0 = ground and terminal 1 = PCC.
  generator->connect({SimNode<Real>::GND, busGas});
  gasTransformer->connect({busGas, busB});
  cableComponent->connect({busB, busA});
  line3Component->connect({busA, busPsha});

  SystemTopology system(simulation.frequency,
                        SystemNodeList{busGas, busB, busA, busPsha},
                        SystemComponentList{
                            generator,
                            gasTransformer,
                            cableComponent,
                            line3Component,
                        });

  system.initWithPowerflow(powerFlow.system, Domain::EMT);

  auto logger = DataLogger::make(simulationName, Logger::Level::off);

  logger->logAttribute("BUS_gas_EMT_v", busGas->attribute("v"));
  logger->logAttribute("BUS_b_EMT_v", busB->attribute("v"));
  logger->logAttribute("BUS_a_EMT_v", busA->attribute("v"));
  logger->logAttribute("BUS_psha_EMT_v", busPsha->attribute("v"));

  logger->logAttribute("TR_gas_EMT_i", gasTransformer->attribute("i_intf"));
  logger->logAttribute("cable_EMT_i", cableComponent->attribute("i_intf"));
  logger->logAttribute("line3_EMT_i", line3Component->attribute("i_intf"));

  logger->logAttribute("GEN_gas_EMT_omega", generator->attribute("omega_gfm"));
  logger->logAttribute("GEN_gas_EMT_theta", generator->attribute("theta_gfm"));
  logger->logAttribute("GEN_gas_EMT_P_elec", generator->attribute("p_inst"));
  logger->logAttribute("GEN_gas_EMT_Q_elec", generator->attribute("q_inst"));
  logger->logAttribute("GEN_gas_EMT_voltage_magnitude",
                       generator->attribute("voltage_magnitude_gfm"));
  logger->logAttribute("GEN_gas_EMT_vc_d", generator->attribute("vc_d"));
  logger->logAttribute("GEN_gas_EMT_vc_q", generator->attribute("vc_q"));
  logger->logAttribute("GEN_gas_EMT_i_grid_d",
                       generator->attribute("i_grid_d"));
  logger->logAttribute("GEN_gas_EMT_i_grid_q",
                       generator->attribute("i_grid_q"));
  logger->logAttribute("GEN_gas_EMT_if_d", generator->attribute("if_d"));
  logger->logAttribute("GEN_gas_EMT_if_q", generator->attribute("if_q"));
  logger->logAttribute("GEN_gas_EMT_i_intf", generator->attribute("i_intf"));
  logger->logAttribute("GEN_gas_EMT_v_intf", generator->attribute("v_intf"));
  logger->logAttribute("GEN_gas_EMT_x", generator->attribute("x"));

  auto log = Logger::get(simulationName + "_parameters", Logger::Level::info,
                         Logger::Level::info);

  log->info("SSN_GFM parameters:"
            "\n  f_n={:.6f} Hz, omega_n={:.9e} rad/s"
            "\n  solved P_ref={:.9e} W, Q_ref={:.9e} var"
            "\n  solved nominal voltage={:.9e} V phase peak"
            "\n  Lf={:.9e} H, Cf={:.9e} F, Rf={:.9e} ohm, Rc={:.9e} ohm"
            "\n  KpV={:.9e}, KiV={:.9e}, KpI={:.9e}, KiI={:.9e}"
            "\n  J={:.9e}, D={:.9e}"
            "\n  inertia scale={:.6f}, damping scale={:.6f}"
            "\n  power-filter cutoff={:.9e} rad/s"
            "\n  converter-delay bandwidth={:.9e} rad/s"
            "\n  linearization interval={} step(s)",
            simulation.frequency, omegaNominal, powerFlow.generatorPower.real(),
            powerFlow.generatorPower.imag(), nominalVoltagePeakPhase, gfm.Lf,
            gfm.Cf, gfm.Rf, gfm.Rc, gfm.kpVoltage(), gfm.kiVoltage(),
            gfm.kpCurrent(), gfm.kiCurrent(), virtualInertia,
            dampingCoefficient, inertiaScale, dampingScale,
            gfm.powerFilterCutoff(simulation.frequency), gfm.delayBandwidth(),
            simulation.linearizationUpdateInterval);

  Simulation emtSimulation(simulationName, Logger::Level::info);

  emtSimulation.setSystem(system);
  emtSimulation.setTimeStep(simulation.timeStep);
  emtSimulation.setFinalTime(simulation.finalTime);
  emtSimulation.setDomain(Domain::EMT);
  emtSimulation.setSolverType(Solver::Type::MNA);
  emtSimulation.doInitFromNodesAndTerminals(true);
  emtSimulation.doSystemMatrixRecomputation(simulation.recomputeSystemMatrix);
  emtSimulation.addLogger(logger);
  emtSimulation.run();
}

} // namespace ScenarioA

int main(int argc, char *argv[]) {
  ScenarioA::Parameters parameters;

  Real inertiaScale = 1.0;
  Real dampingScale = 1.0;

  CommandLineArgs args(argc, argv);

  if (argc > 1) {
    parameters.simulation.timeStep = args.timeStep;
    parameters.simulation.finalTime = args.duration;

    if (args.name != "dpsim")
      parameters.simulation.name = args.name;

    if (args.options.find("SCALEINERTIA_G1") != args.options.end())
      inertiaScale = args.getOptionReal("SCALEINERTIA_G1");

    if (args.options.find("SCALEDAMPING_G1") != args.options.end())
      dampingScale = args.getOptionReal("SCALEDAMPING_G1");

    if (args.options.find("LINEARIZATION_INTERVAL") != args.options.end()) {
      const Int interval = args.getOptionInt("LINEARIZATION_INTERVAL");

      if (interval <= 0)
        throw std::invalid_argument(
            "LINEARIZATION_INTERVAL must be at least one.");

      parameters.simulation.linearizationUpdateInterval =
          static_cast<UInt>(interval);
    }
  }

  ScenarioA::runScenario(parameters, inertiaScale, dampingScale);

  return 0;
}
