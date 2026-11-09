// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
// Closed-loop IEEE 9-bus comparison using native, FMI ME, or FMI CS inverter.
#include "../../../dpsim/examples/cxx/Examples.h"
#include "../small_inverter_grid/NativeAveragedInverter.h"

#include <DPsim.h>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim-models/Signal/ExciterDC1Simp.h>
#include <dpsim-models/Signal/TurbineGovernorType1.h>
#include <dpsim/SequentialScheduler.h>

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace CPS;
using namespace DPsim;

namespace {
using Scenario = CPS::CIM::Examples::Grids::IEEE9::ScenarioConfig;
using IEEE9Line = CPS::CIM::Examples::Grids::IEEE9::Line;
using IEEE9Transformer = CPS::CIM::Examples::Grids::IEEE9::Transformer;
constexpr Real frequency = 60.0;
constexpr Real stepPower = 25e6;
constexpr Real currentLimitPu = 1.20;
constexpr Real responseTime = 5e-3;
constexpr Real phasePeakPerLineVoltage = 0.8164965809277260;

Matrix parameter3(Real value) {
  return Math::singlePhaseParameterToThreePhase(value);
}

SystemTopology makePowerflowSystem(const Scenario &grid) {
  std::vector<SP::SimNode::Ptr> nodes;
  SystemNodeList systemNodes;
  nodes.reserve(9);
  for (int bus = 1; bus <= 9; ++bus) {
    auto node =
        SP::SimNode::make("BUS" + std::to_string(bus), PhaseType::Single);
    nodes.push_back(node);
    systemNodes.push_back(node);
  }

  SystemComponentList components;
  for (const auto *generator : {&grid.gen1, &grid.gen2}) {
    auto machine = SP::Ph1::SynchronGenerator::make(generator->Name);
    machine->setParameters(generator->RatedPower, generator->RatedVoltage,
                           generator->InitialPower, generator->InitialVoltage,
                           generator->BusType, generator->InitialPowerReactive,
                           2e8, -2e8);
    machine->setBaseVoltage(generator->RatedVoltage);
    const int bus = generator == &grid.gen1 ? 1 : 2;
    machine->connect({nodes.at(bus - 1)});
    components.push_back(machine);
  }

  // The averaged inverter replaces the standard bus-3 synchronous generator.
  auto inverterPowerflow = SP::Ph1::Load::make("GEN3_PF");
  inverterPowerflow->setParameters(-grid.gen3.InitialPower,
                                   -grid.gen3.InitialPowerReactive,
                                   grid.gen3.RatedVoltage);
  inverterPowerflow->modifyPowerFlowBusType(PowerflowBusType::PQ);
  inverterPowerflow->connect({nodes.at(2)});
  components.push_back(inverterPowerflow);

  for (const auto *load : {&grid.load5, &grid.load6, &grid.load8}) {
    auto component = SP::Ph1::Load::make(load->Name);
    component->setParameters(load->RealPower, load->ReactivePower,
                             load->BaseVoltage);
    component->modifyPowerFlowBusType(PowerflowBusType::PQ);
    const int bus = load == &grid.load5 ? 5 : (load == &grid.load6 ? 6 : 8);
    component->connect({nodes.at(bus - 1)});
    components.push_back(component);
  }

  struct LineData {
    const IEEE9Line *line;
    int from;
    int to;
  };
  const std::array<LineData, 6> lines{{{&grid.line54, 5, 4},
                                       {&grid.line64, 6, 4},
                                       {&grid.line75, 7, 5},
                                       {&grid.line96, 9, 6},
                                       {&grid.line78, 7, 8},
                                       {&grid.line89, 8, 9}}};
  for (const auto &entry : lines) {
    auto line = SP::Ph1::PiLine::make(entry.line->Name);
    line->setParameters(entry.line->Resistance, entry.line->Inductance,
                        entry.line->Capacitance, entry.line->Conductance);
    line->setBaseVoltage(entry.line->BaseVoltage);
    line->connect({nodes.at(entry.from - 1), nodes.at(entry.to - 1)});
    components.push_back(line);
  }

  struct TransformerData {
    const IEEE9Transformer *transformer;
    int lowBus;
    int highBus;
  };
  const std::array<TransformerData, 3> transformers{
      {{&grid.transf14, 1, 4}, {&grid.transf27, 2, 7}, {&grid.transf39, 3, 9}}};
  for (const auto &entry : transformers) {
    const auto &data = *entry.transformer;
    auto transformer = SP::Ph1::Transformer::make(data.Name);
    transformer->setParameters(data.VoltageLVSide, data.VoltageHVSide,
                               data.RatedPower, data.Ratio, 0.0,
                               data.Resistance, data.Inductance);
    transformer->setBaseVoltage(data.VoltageHVSide);
    transformer->connect(
        {nodes.at(entry.lowBus - 1), nodes.at(entry.highBus - 1)});
    components.push_back(transformer);
  }
  return SystemTopology(frequency, systemNodes, components);
}

void runCase(const Scenario &grid, const SystemTopology &powerflow,
             const std::string &mode, const std::string &fmuPath, Real timeStep,
             Real finalTime, Real loadStepTime) {
  const String simName = "IEEE9_Inverter_" + mode;
  const String logDirectory = "logs/" + simName;
  Logger::setLogDir(logDirectory);

  std::vector<EMT::SimNode::Ptr> nodes;
  SystemNodeList systemNodes;
  nodes.reserve(9);
  for (int bus = 1; bus <= 9; ++bus) {
    auto node = EMT::SimNode::make("BUS" + std::to_string(bus), PhaseType::ABC);
    nodes.push_back(node);
    systemNodes.push_back(node);
  }
  SystemComponentList components;

  // Keep the standard bus-1 and bus-2 synchronous machines.
  for (const auto *generator : {&grid.gen1, &grid.gen2}) {
    auto machine = EMT::Ph3::SynchronGenerator4OrderVBR::make(generator->Name);
    machine->setOperationalParametersPerUnit(
        generator->RatedPower, generator->RatedVoltage, grid.nomFreq,
        generator->H, generator->Xd, generator->Xq, generator->Xa,
        generator->XdPrime, generator->XqPrime, generator->TdoPrime,
        generator->TqoPrime);

    // Keep the IEEE9 excitation and turbine-governor dynamics on the two
    // synchronous machines.  Without them, the 25 MW load step leaves a
    // sustained generation deficit that causes frequency/voltage collapse in
    // long EMT runs, even though the powerflow initial state is valid.
    const auto &exciterData = generator == &grid.gen1 ? grid.exc1 : grid.exc2;
    const auto &governorData = generator == &grid.gen1 ? grid.gov1 : grid.gov2;
    auto exciterParams = std::make_shared<Signal::ExciterDC1SimpParameters>();
    exciterParams->Ta = exciterData.TA;
    exciterParams->Ka = exciterData.KA;
    exciterParams->Tef = exciterData.TE;
    exciterParams->Kef = exciterData.KE;
    exciterParams->Tf = exciterData.TF;
    exciterParams->Kf = exciterData.KF;
    exciterParams->Tr = 0.01;
    exciterParams->MaxVa = exciterData.VRmax;
    exciterParams->MinVa = exciterData.VRmin;
    exciterParams->Bef = std::log(exciterData.S_EX2 / exciterData.S_EX1) /
                         (exciterData.EX2 - exciterData.EX1);
    exciterParams->Aef =
        exciterData.S_EX1 / std::exp(exciterParams->Bef * exciterData.EX1);
    auto exciter = Signal::ExciterDC1Simp::make(generator->Name + "_Exciter",
                                                Logger::Level::off);
    exciter->setParameters(exciterParams);
    machine->addExciter(exciter);

    auto governor = Signal::TurbineGovernorType1::make(
        generator->Name + "_TurbineGovernor", Logger::Level::off);
    governor->setParameters(governorData.T2, 1.0, 1.0, governorData.T3,
                            governorData.T1, governorData.R, governorData.Vmin,
                            governorData.Vmax, 1.0);
    machine->addGovernor(governor);
    const int bus = generator == &grid.gen1 ? 1 : 2;
    machine->connect({nodes.at(bus - 1)});
    components.push_back(machine);
  }

  for (const auto *load : {&grid.load5, &grid.load6, &grid.load8}) {
    auto component = EMT::Ph3::RXLoad::make(load->Name);
    component->setParameters(
        Math::singlePhasePowerToThreePhase(load->RealPower),
        Math::singlePhasePowerToThreePhase(load->ReactivePower),
        load->BaseVoltage);
    const int bus = load == &grid.load5 ? 5 : (load == &grid.load6 ? 6 : 8);
    component->connect({nodes.at(bus - 1)});
    components.push_back(component);
  }

  struct LineData {
    const IEEE9Line *line;
    int from;
    int to;
  };
  const std::array<LineData, 6> lines{{{&grid.line54, 5, 4},
                                       {&grid.line64, 6, 4},
                                       {&grid.line75, 7, 5},
                                       {&grid.line96, 9, 6},
                                       {&grid.line78, 7, 8},
                                       {&grid.line89, 8, 9}}};
  for (const auto &entry : lines) {
    auto line = EMT::Ph3::PiLine::make(entry.line->Name);
    line->setParameters(parameter3(entry.line->Resistance),
                        parameter3(entry.line->Inductance),
                        parameter3(entry.line->Capacitance),
                        parameter3(entry.line->Conductance));
    line->connect({nodes.at(entry.from - 1), nodes.at(entry.to - 1)});
    components.push_back(line);
  }

  struct TransformerData {
    const IEEE9Transformer *transformer;
    int lowBus;
    int highBus;
  };
  const std::array<TransformerData, 3> transformers{
      {{&grid.transf14, 1, 4}, {&grid.transf27, 2, 7}, {&grid.transf39, 3, 9}}};
  for (const auto &entry : transformers) {
    const auto &data = *entry.transformer;
    auto transformer = EMT::Ph3::Transformer::make(data.Name);
    transformer->setParameters(
        data.VoltageLVSide, data.VoltageHVSide, data.RatedPower, data.Ratio,
        0.0, parameter3(data.Resistance), parameter3(data.Inductance));
    transformer->connect(
        {nodes.at(entry.lowBus - 1), nodes.at(entry.highBus - 1)});
    components.push_back(transformer);
  }

  FMI::NewtonSettings newtonSettings;
  newtonSettings.multirateMode = FMI::MultirateMode::SinglePassExplicit;
  FMI::CoSimulationElectricalPort csPort{
      {"v_a", "v_b", "v_c"}, {"i_a", "i_b", "i_c"}, 1.0};
  std::shared_ptr<SimPowerComp<Real>> inverter;
  if (mode == "native") {
    auto native = SmallGrid::NativeAveragedInverter::make("GEN3_Native");
    native->setOperatingPoint(
        grid.gen3.InitialPower, grid.gen3.InitialPowerReactive,
        grid.gen3.RatedVoltage, currentLimitPu, responseTime, grid.nomFreq);
    native->setUsePeriodicFilterStart(false);
    inverter = native;
  } else if (mode == "me") {
    FMI::ModelExchangeOptions options;
    auto model = std::make_shared<FMI::ModelExchange>(fmuPath, options);
    inverter = FMI::NonlinearSSN::make(
        "GEN3_AvgInverter", model,
        FMI::ElectricalPort3Ph{
            {"v_a", "v_b", "v_c"}, {"i_a", "i_b", "i_c"}, 1.0},
        newtonSettings);
  } else if (mode == "cs" || mode == "cs-newton" || mode == "cs-explicit") {
    const bool rollbackNewton = mode != "cs-explicit";
    FMI::CoSimulationOptions options;
    options.enableRollback = rollbackNewton;
    auto model = std::make_shared<FMI::CoSimulation>(fmuPath, options);
    const auto coupling =
        rollbackNewton ? FMI::CoSimulationCouplingMode::RollbackNewton
                       : FMI::CoSimulationCouplingMode::SinglePassExplicit;
    inverter =
        FMI::CoSimulationSSN::make("GEN3_AvgInverter", model, csPort, coupling);
  } else {
    throw std::invalid_argument("mode must be native, me, cs, or cs-newton");
  }
  inverter->connect({EMT::SimNode::GND, nodes.at(2)});
  components.push_back(inverter);

  // Add a 25 MW resistive load step at bus 8 at the run midpoint.
  const Real stepBusPeak =
      1.02 * grid.load8.BaseVoltage * phasePeakPerLineVoltage;
  const Real stepResistance = 1.5 * stepBusPeak * stepBusPeak / stepPower;
  auto stepNode = EMT::SimNode::make("BUS8_StepNode", PhaseType::ABC);
  auto stepLoad = EMT::Ph3::Resistor::make("BUS8_StepLoad");
  stepLoad->setParameters(parameter3(stepResistance));
  stepLoad->connect({EMT::SimNode::GND, stepNode});
  auto stepSwitch = EMT::Ph3::Switch::make("BUS8_StepSwitch");
  stepSwitch->setParameters(parameter3(1e12), parameter3(1e-3), false);
  stepSwitch->connect({stepNode, nodes.at(7)});
  systemNodes.push_back(stepNode);
  components.push_back(stepLoad);
  components.push_back(stepSwitch);

  auto system = SystemTopology(grid.nomFreq, systemNodes, components);
  system.initWithPowerflow(powerflow, Domain::EMT);
  auto logger = DataLogger::make(simName, true, 1);
  for (int bus : {1, 3, 4, 5, 8, 9})
    logger->logAttribute("BUS" + std::to_string(bus) + ".V",
                         nodes.at(bus - 1)->attribute("v"));
  logger->logAttribute("inverter.v", inverter->attribute("v_intf"));
  logger->logAttribute("inverter.i", inverter->attribute("i_intf"));
  if (mode != "native") {
    for (const auto *phase : {"a", "b", "c"}) {
      logger->logAttribute(String("fmi.input.v_") + phase,
                           inverter->attribute(String("fmi.input.v_") + phase));
      logger->logAttribute(
          String("fmi.output.i_") + phase,
          inverter->attribute(String("fmi.output.i_") + phase));
    }
  }
  logger->logAttribute("step_load.i", stepLoad->attribute("i_intf"));
  if (mode == "native" || mode == "me" || mode == "cs" || mode == "cs-newton") {
    logger->logAttribute("newton_iterations",
                         inverter->attribute("newton_iterations"));
    logger->logAttribute("newton_residual",
                         inverter->attribute("newton_residual"));
  } else {
    logger->logAttribute("fmi_do_steps", inverter->attribute("cosim_do_steps"));
  }

  Simulation sim(simName, Logger::Level::off);
  sim.addEvent(SwitchEvent3Ph::make(loadStepTime, stepSwitch, true));
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.doInitFromNodesAndTerminals(true);
  sim.doSystemMatrixRecomputation(true);
  sim.addLogger(logger);
  sim.setScheduler(std::make_shared<SequentialScheduler>(
      logDirectory + "/" + simName + "_TaskTimes.csv", Logger::Level::off));
  sim.setLogStepTimes(true);

  const auto initializeStart = std::chrono::steady_clock::now();
  sim.initialize();
  const std::chrono::duration<double> initializeElapsed =
      std::chrono::steady_clock::now() - initializeStart;
  const auto simulationStart = std::chrono::steady_clock::now();
  sim.run();
  const std::chrono::duration<double> simulationElapsed =
      std::chrono::steady_clock::now() - simulationStart;
  sim.logLUTimes();

  std::ofstream stepTiming(logDirectory + "/" + simName +
                           "_SimulationStepTimes.csv");
  stepTiming << "step,simulation_time_s,step_wall_s\n";
  const auto &stepTimes = sim.stepTimes();
  for (std::size_t index = 0; index < stepTimes.size(); ++index)
    stepTiming << index + 1 << "," << (index + 1) * timeStep << ","
               << stepTimes[index] << "\n";

  std::cout << simName << " runtime_s="
            << initializeElapsed.count() + simulationElapsed.count()
            << " initialization_s=" << initializeElapsed.count()
            << " simulation_s=" << simulationElapsed.count()
            << " steps=" << std::llround(finalTime / timeStep)
            << " dt_s=" << timeStep << " stop_s=" << finalTime
            << " load_step_s=" << loadStepTime << " load_step_MW=25\n";
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cerr << "Usage: IEEE9_AvgInverter native|me|cs|cs-explicit "
                 "inverter.fmu stop_s "
                 "step_s load_step_s\n";
    return 2;
  }
  try {
    const std::string mode = argv[1];
    if (mode != "native" && mode != "me" && mode != "cs" &&
        mode != "cs-newton" && mode != "cs-explicit")
      throw std::invalid_argument(
          "mode must be native, me, cs, cs-newton, or cs-explicit");
    const Real finalTime = std::stod(argv[3]);
    const Real timeStep = std::stod(argv[4]);
    const Real loadStepTime = std::stod(argv[5]);
    if (finalTime <= 0.0 || timeStep <= 0.0 || loadStepTime <= 0.0 ||
        loadStepTime >= finalTime ||
        std::abs(finalTime / timeStep - std::round(finalTime / timeStep)) >
            1e-8 ||
        std::abs(loadStepTime / timeStep -
                 std::round(loadStepTime / timeStep)) > 1e-8)
      throw std::invalid_argument("invalid time grid");

    Scenario grid(frequency);
    auto powerflow = makePowerflowSystem(grid);
    Simulation pf("IEEE9_AvgInverter_SP_PF", Logger::Level::off);
    pf.setSystem(powerflow);
    pf.setDomain(Domain::SP);
    pf.setSolverType(Solver::Type::NRP);
    pf.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
    pf.setTimeStep(0.01);
    pf.setFinalTime(0.01);
    pf.doInitFromNodesAndTerminals(false);
    pf.run();
    runCase(grid, powerflow, mode, argv[2], timeStep, finalTime, loadStepTime);
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
