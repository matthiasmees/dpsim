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
  if (argc < 3 || argc > 6) {
    std::cerr << "Usage: FMI_EMT_Ph3_Nonlinear model.fmu benchmark.json "
                 "[step_s [tolerance_multiplier [stop_time_s]]]\n";
    return 2;
  }
  try {

    // ----- PARAMETERS AND FMU INTERFACE -----

    nlohmann::json config;
    std::ifstream input(argv[2]);
    input >> config;
    const auto voltageInputs =
        config.at("voltage_inputs").get<std::array<std::string, 3>>();
    const auto currentOutputs =
        config.at("current_outputs").get<std::array<std::string, 3>>();
    const auto controlNames =
        config.at("control_names").get<std::vector<std::string>>();
    const auto controlValues =
        config.at("control_initial_values").get<std::vector<double>>();
    const auto initialNetworkValues =
        config.at("initial_network_state").get<std::vector<double>>();
    const auto initialGeneratorCurrents =
        config.at("initial_generator_current").get<std::vector<double>>();
    if (initialNetworkValues.size() != 6 ||
        initialGeneratorCurrents.size() != 3 ||
        controlNames.size() != controlValues.size()) {
      throw std::invalid_argument("Invalid benchmark configuration dimensions");
    }

    // ----- INITIAL CONDITIONS -----

    Matrix initialNetworkState(6, 1), initialGeneratorCurrent(3, 1);
    std::vector<Complex> initialNodeVoltages(3);
    std::vector<FMI::ScalarInput> inputs;
    for (int k = 0; k < 6; ++k) {
      initialNetworkState(k, 0) = initialNetworkValues[k];
    }
    for (int phase = 0; phase < 3; ++phase) {
      initialGeneratorCurrent(phase, 0) = initialGeneratorCurrents[phase];
      initialNodeVoltages[phase] =
          initialNetworkValues[phase] / RMS3PH_TO_PEAK1PH;
      inputs.push_back({voltageInputs[phase], initialNetworkValues[phase]});
    }
    for (std::size_t k = 0; k < controlNames.size(); ++k) {
      inputs.push_back({controlNames[k], controlValues[k]});
    }

    // ----- IMPORTED GENERATOR -----

    auto model = std::make_shared<FMI::ModelExchange>(argv[1], inputs);
    FMI::NewtonSettings newtonSettings;
    if (argc == 5) {
      const double toleranceScale = std::stod(argv[4]);
      newtonSettings.absoluteStateTolerance *= toleranceScale;
      newtonSettings.relativeStateTolerance *= toleranceScale;
      newtonSettings.absoluteCurrentTolerance *= toleranceScale;
      newtonSettings.relativeCurrentTolerance *= toleranceScale;
    }
    auto generator = FMI::NonlinearSSN::make(
        "imported", model,
        FMI::ElectricalPort3Ph{
            voltageInputs, currentOutputs,
            config.at("current_sign_to_passive").get<double>()},
        newtonSettings);

    // ----- REFERENCE AND LOAD DISTURBANCES -----

    const double voltageReferenceStepTime = config.at("vref_step_s"),
                 voltageReferenceStep = config.at("vref_step_pu");
    generator->setInputFunction(
        "vref_pu",
        [=](double time) {
          return 1 + voltageReferenceStep * (time >= voltageReferenceStepTime);
        },
        {voltageReferenceStepTime});
    const double frequency = config.at("frequency_Hz"),
                 phaseVoltageRms = config.at("phase_rms_V");
    const double initialLoadConductance = config.at("G0_S"),
                 loadConductanceStep = config.at("Gstep_S"),
                 loadStepTime = config.at("load_step_s");
    const auto gridVoltage = [=](double time) -> Matrix {
      Matrix phaseVoltages(3, 1);
      const double angle = 2 * PI * frequency * time;
      phaseVoltages << std::sqrt(2.0) * phaseVoltageRms * std::cos(angle),
          std::sqrt(2.0) * phaseVoltageRms * std::cos(angle - 2 * PI / 3),
          std::sqrt(2.0) * phaseVoltageRms * std::cos(angle + 2 * PI / 3);
      return phaseVoltages;
    };
    const auto loadConductance = [=](double time) {
      return initialLoadConductance +
             loadConductanceStep * (time >= loadStepTime - 1e-12);
    };

    // ----- EXTERNAL NETWORK COMPONENT -----

    // The reference network retains its coupled R/L/C state equations.
    auto network = LinearNetwork3Ph::make(
        "external_grid", config.at("network")[0].get<double>(),
        config.at("network")[1].get<double>(),
        config.at("network")[2].get<double>(), initialNetworkState,
        initialGeneratorCurrent, gridVoltage, loadConductance);

    // ----- NODES AND TOPOLOGY -----

    auto bus = SimNode<Real>::make("bus", PhaseType::ABC, initialNodeVoltages);
    generator->connect(bus);
    network->connect({SimNode<Real>::GND, bus});
    SystemTopology system(frequency, {bus}, {generator, network});

    // ----- LOGGING -----

    auto logger = DataLogger::make(String("fmi_sg_nonlinear"));
    for (auto &name : model->description().outputNames) {
      logger->logAttribute(name, generator->attribute("fmi.output." + name));
    }
    for (auto &name : model->description().inputNames) {
      logger->logAttribute("input_" + name,
                           generator->attribute("fmi.input." + name));
    }
    logger->logAttribute("newton_iterations",
                         generator->attribute("newton_iterations"));
    logger->logAttribute("newton_residual",
                         generator->attribute("newton_residual"));
    logger->logAttribute("network_current", network->attribute("i_intf"));
    for (const std::string name :
         {"v_bus_a_V", "v_bus_b_V", "v_bus_c_V", "i_line_a_A", "i_line_b_A",
          "i_line_c_A", "g_load_S", "p_load_W", "p_line_bus_W",
          "p_line_loss_W"}) {
      logger->logAttribute(name, network->attribute(name));
    }

    // ----- SIMULATION -----

    const double timeStep =
        argc >= 4 ? std::stod(argv[3]) : config.at("step_s").get<double>();
    if (std::abs(loadStepTime / timeStep -
                 std::round(loadStepTime / timeStep)) > 1e-7) {
      throw std::invalid_argument("Load event must align with h");
    }
    Simulation simulation("fmi_sg_nonlinear", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    const double finalTime =
        argc == 6 ? std::stod(argv[5]) : config.at("stop_time_s").get<double>();
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed nonlinear FMI/MNA trapezoidal simulation with "
              << model->description().initialState.rows()
              << " imported states, h=" << timeStep << ".\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
