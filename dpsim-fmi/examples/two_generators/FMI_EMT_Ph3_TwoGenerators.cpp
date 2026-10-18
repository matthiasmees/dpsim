// SPDX-License-Identifier: MPL-2.0
#include "LinearNetwork3Ph.h"
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc < 4 || argc > 7) {
    std::cerr << "Usage: FMI_EMT_Ph3_TwoGenerators first.fmu second.fmu "
                 "benchmark.json [step_s [tolerance_scale [stop_time_s]]]\n";
    return 2;
  }
  try {
    // ----- PARAMETERS AND OPERATING POINT -----
    nlohmann::json config;
    std::ifstream file(argv[3]);
    file >> config;
    const Real frequency = config.at("frequency_Hz");
    const Real phaseVoltageRms = config.at("phase_rms_V");
    const Real timeStep = argc >= 5 ? std::stod(argv[4]) : 50e-6;
    const Real toleranceScale = argc >= 6 ? std::stod(argv[5]) : 1.0;
    const Real finalTime = argc >= 7 ? std::stod(argv[6]) : 0.1;
    const auto voltages =
        config.at("voltage_inputs").get<std::array<std::string, 3>>();
    const auto currents =
        config.at("current_outputs").get<std::array<std::string, 3>>();
    const auto networkValues =
        config.at("initial_network_state").get<std::vector<double>>();
    const auto generatorCurrents =
        config.at("initial_generator_current").get<std::vector<double>>();
    Matrix initialNetworkState(6, 1), initialTotalCurrent(3, 1);
    std::vector<Complex> initialNodeVoltages(3);
    std::vector<FMI::ScalarInput> initialInputs;
    for (int phase = 0; phase < 3; ++phase) {
      initialNetworkState(phase, 0) = networkValues[phase];
      initialNetworkState(phase + 3, 0) = 2 * networkValues[phase + 3];
      initialTotalCurrent(phase, 0) = 2 * generatorCurrents[phase];
      initialNodeVoltages[phase] = networkValues[phase] / RMS3PH_TO_PEAK1PH;
      initialInputs.push_back({voltages[phase], networkValues[phase]});
    }
    initialInputs.push_back({"vref_pu", 1.0});
    initialInputs.push_back({"pmref_pu", config.at("pmref_pu")});

    // ----- TWO INDEPENDENT FMI INSTANCES -----
    // Separate model identities respect MATLAB's singleton export capability.
    auto firstModel =
        std::make_shared<FMI::ModelExchange>(argv[1], initialInputs);
    auto secondModel =
        std::make_shared<FMI::ModelExchange>(argv[2], initialInputs);
    FMI::NewtonSettings settings;
    settings.absoluteStateTolerance *= toleranceScale;
    settings.relativeStateTolerance *= toleranceScale;
    settings.absoluteCurrentTolerance *= toleranceScale;
    settings.relativeCurrentTolerance *= toleranceScale;
    const FMI::ElectricalPort3Ph port{voltages, currents, -1.0};
    auto firstGenerator =
        FMI::NonlinearSSN::make("first_generator", firstModel, port, settings);
    auto secondGenerator = FMI::NonlinearSSN::make("second_generator",
                                                   secondModel, port, settings);

    // ----- DISTINCT CONTROL DISTURBANCES -----
    firstGenerator->setInputFunction(
        "vref_pu", [](Real time) { return time >= 0.05 ? 1.02 : 1.0; }, {0.05});
    secondGenerator->setInputFunction(
        "vref_pu", [](Real time) { return time >= 0.075 ? 0.99 : 1.0; },
        {0.075});

    // ----- SHARED EXTERNAL NETWORK -----
    const auto gridVoltage = [=](Real time) -> Matrix {
      Matrix voltage(3, 1);
      for (int phase = 0; phase < 3; ++phase) {
        voltage(phase, 0) =
            std::sqrt(2.0) * phaseVoltageRms *
            std::cos(2 * PI * frequency * time - phase * 2 * PI / 3);
      }
      return voltage;
    };
    const Real loadConductance = 2 * config.at("G0_S").get<Real>();
    auto network = LinearNetwork3Ph::make(
        "shared_network", 2 * config.at("network")[0].get<Real>(),
        config.at("network")[1].get<Real>() / 2,
        config.at("network")[2].get<Real>() / 2, initialNetworkState,
        initialTotalCurrent, gridVoltage,
        [=](Real) { return loadConductance; });

    // ----- NODES AND TOPOLOGY -----
    auto pcc = SimNode<Real>::make("pcc", PhaseType::ABC, initialNodeVoltages);
    firstGenerator->connect({SimNode<Real>::GND, pcc});
    secondGenerator->connect({SimNode<Real>::GND, pcc});
    network->connect({SimNode<Real>::GND, pcc});
    SystemTopology system(frequency, {pcc},
                          {firstGenerator, secondGenerator, network});

    // ----- LOGGING -----
    auto logger = DataLogger::make(String("fmi_two_generators"));
    for (const auto &name : firstModel->description().outputNames) {
      logger->logAttribute("first_" + name,
                           firstGenerator->attribute("fmi.output." + name));
      logger->logAttribute("second_" + name,
                           secondGenerator->attribute("fmi.output." + name));
    }
    logger->logAttribute("first_vref",
                         firstGenerator->attribute("fmi.input.vref_pu"));
    logger->logAttribute("second_vref",
                         secondGenerator->attribute("fmi.input.vref_pu"));
    logger->logAttribute("first_current", firstGenerator->attribute("i_intf"));
    logger->logAttribute("second_current",
                         secondGenerator->attribute("i_intf"));
    logger->logAttribute("network_current", network->attribute("i_intf"));
    logger->logAttribute("pcc_voltage", pcc->attribute("v"));
    logger->logAttribute("network_state", network->attribute("x"));
    logger->logAttribute("first_residual",
                         firstGenerator->attribute("newton_residual"));
    logger->logAttribute("second_residual",
                         secondGenerator->attribute("newton_residual"));

    // ----- SIMULATION -----
    Simulation simulation("fmi_two_generators", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Two generators completed: independent controls, shared MNA "
                 "network\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
