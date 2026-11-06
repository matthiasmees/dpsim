// SPDX-License-Identifier: MPL-2.0
#include "NativeAveragedInverter.h"
#include "SmallGridNetwork.h"
#include <chrono>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Event.h>
#include <dpsim/SequentialScheduler.h>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc != 6 && argc != 7) {
    std::cerr << "Usage: FMI_SmallGrid_Inverter me|cs inverter.fmu stop_s "
                 "step_s load_step_s [manual|powerflow-unadjusted|powerflow]\n";
    return 2;
  }
  try {
    const std::string mode = argv[1];
    const std::string fmuPath = argv[2];
    const Real stopTime = std::stod(argv[3]);
    const Real timeStep = std::stod(argv[4]);
    const Real loadStep = std::stod(argv[5]);
    const std::string initialization = argc == 7 ? argv[6] : "manual";
    if ((mode != "native" && mode != "me" && mode != "cs") ||
        (initialization != "manual" &&
         initialization != "powerflow-unadjusted" &&
         initialization != "powerflow") ||
        stopTime <= 0 || timeStep <= 0 || loadStep <= 0 ||
        loadStep >= stopTime ||
        std::abs(stopTime / timeStep - std::round(stopTime / timeStep)) >
            1e-8 ||
        std::abs(loadStep / timeStep - std::round(loadStep / timeStep)) > 1e-8)
      throw std::invalid_argument("invalid mode or simulation time grid");

    const bool usePowerflow = initialization != "manual";
    const bool reinitializeFMU = initialization == "powerflow";
    const String simulationName = "SmallGrid_Inverter_" + mode +
                                  (usePowerflow ? "_" + initialization : "");
    Logger::setLogDir("logs/" + simulationName);
    std::unique_ptr<Simulation> powerflowSimulation;
    std::unique_ptr<SmallGrid::PowerflowParts> powerflow;
    std::vector<Complex> sourceStart = SmallGrid::sourcePhasors();
    std::vector<Complex> pccStart = SmallGrid::pccPhasors();
    if (usePowerflow) {
      powerflow = std::make_unique<SmallGrid::PowerflowParts>(
          SmallGrid::makePowerflowParts());
      powerflowSimulation = std::make_unique<Simulation>(
          "SmallGrid_Inverter_PF", Logger::Level::off);
      powerflowSimulation->setSystem(*powerflow->system);
      powerflowSimulation->setDomain(Domain::SP);
      powerflowSimulation->setSolverType(Solver::Type::NRP);
      powerflowSimulation->setSolverAndComponentBehaviour(
          Solver::Behaviour::Initialization);
      powerflowSimulation->setTimeStep(0.01);
      powerflowSimulation->setFinalTime(0.01);
      powerflowSimulation->doInitFromNodesAndTerminals(false);
      powerflowSimulation->run();
      sourceStart =
          SmallGrid::threePhasePhasors(powerflow->sourceNode->singleVoltage());
      pccStart =
          SmallGrid::threePhasePhasors(powerflow->pccNode->singleVoltage());
    }
    SmallGrid::NetworkParts network =
        SmallGrid::makeNetworkParts(sourceStart, pccStart);
    FMI::NewtonSettings newton;
    // Use residual-based Newton stopping (default ceiling: 20 iterations).
    // The solver exits as soon as the state and terminal-current residuals
    // meet their tolerances; this is an upper bound, not a fixed count.
    newton.multirateMode = FMI::MultirateMode::SinglePassExplicit;
    FMI::ElectricalPort3Ph mePort{
        {"v_a", "v_b", "v_c"}, {"i_a", "i_b", "i_c"}, 1.0};
    FMI::CoSimulationElectricalPort csPort{
        {"v_a", "v_b", "v_c"}, {"i_a", "i_b", "i_c"}, 1.0};
    SimPowerComp<Real>::Ptr inverter;
    if (mode == "native") {
      inverter = SmallGrid::NativeAveragedInverter::make("native_inverter");
    } else if (mode == "me") {
      FMI::ModelExchangeOptions options;
      options.initializeFromPowerflow = reinitializeFMU;
      auto model = std::make_shared<FMI::ModelExchange>(fmuPath, options);
      inverter = FMI::NonlinearSSN::make("fmi_inverter", model, mePort, newton);
    } else {
      FMI::CoSimulationOptions options;
      options.initializeFromPowerflow = reinitializeFMU;
      auto model = std::make_shared<FMI::CoSimulation>(fmuPath, options);
      inverter = FMI::CoSimulationSSN::make("fmi_inverter", model, csPort);
    }
    inverter->connect({EMT::SimNode::GND, network.pccNode});
    auto system = SmallGrid::makeSystem(network, inverter);
    if (usePowerflow)
      system.initWithPowerflow(*powerflow->system, Domain::EMT);

    auto logger = DataLogger::make(simulationName, true, 1);
    logger->logAttribute("pcc.v", network.pccNode->attribute("v"));
    logger->logAttribute("inverter.v", inverter->attribute("v_intf"));
    logger->logAttribute("inverter.i", inverter->attribute("i_intf"));
    if (mode == "native" || mode == "me") {
      logger->logAttribute("newton_iterations",
                           inverter->attribute("newton_iterations"));
      logger->logAttribute("newton_residual",
                           inverter->attribute("newton_residual"));
    } else {
      logger->logAttribute("fmi_do_steps",
                           inverter->attribute("cosim_do_steps"));
    }

    Simulation simulation(simulationName, Logger::Level::off);
    simulation.addEvent(
        SwitchEvent3Ph::make(loadStep, network.loadSwitch, true));
    simulation.setSystem(system);
    simulation.setDomain(Domain::EMT);
    simulation.setSolverType(Solver::Type::MNA);
    simulation.setTimeStep(timeStep);
    // The fixed-step clock accumulates round-off. Extend the internal stop by
    // a nanosecond so the data logger captures the sample at requested stop_s
    // for step sizes where repeated addition lands a few ulps below/above it.
    simulation.setFinalTime(stopTime + 1e-9);
    simulation.doInitFromNodesAndTerminals(true);
    simulation.doSystemMatrixRecomputation(true);
    simulation.addLogger(logger);
    const std::string logDir = "logs/" + simulationName;
    simulation.setScheduler(std::make_shared<SequentialScheduler>(
        logDir + "/" + simulationName + "_TaskTimes.csv", Logger::Level::off));
    simulation.setLogStepTimes(true);

    const auto start = std::chrono::steady_clock::now();
    const auto initStart = start;
    simulation.initialize();
    const auto initialized = std::chrono::steady_clock::now();
    simulation.run();
    const auto finish = std::chrono::steady_clock::now();
    const auto initSeconds =
        std::chrono::duration<double>(initialized - initStart).count();
    const auto simulationSeconds =
        std::chrono::duration<double>(finish - initialized).count();
    const auto wallSeconds =
        std::chrono::duration<double>(finish - start).count();
    simulation.logLUTimes();
    std::ofstream stepTiming(logDir + "/" + simulationName +
                             "_SimulationStepTimes.csv");
    stepTiming << "step,simulation_time_s,step_wall_s\n";
    for (std::size_t index = 0; index < simulation.stepTimes().size(); ++index)
      stepTiming << index + 1 << "," << (index + 1) * timeStep << ","
                 << simulation.stepTimes()[index] << '\n';
    std::cout << "small_grid mode=" << mode
              << " initialization=" << initialization
              << " steps=" << std::llround(stopTime / timeStep)
              << " dt_s=" << timeStep << " stop_s=" << stopTime
              << " load_step_s=" << loadStep << " init_s=" << initSeconds
              << " simulation_s=" << simulationSeconds
              << " wall_s=" << wallSeconds << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
