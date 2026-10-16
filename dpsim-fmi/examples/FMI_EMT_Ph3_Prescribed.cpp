// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

// Prescribed voltage benchmark: all imported branch dynamics come from the FMU.
int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    std::cerr
        << "Usage: FMI_EMT_Ph3_Prescribed model.fmu config.json [step_s]\n";
    return 2;
  }
  try {
    using namespace CPS;
    using namespace DPsim;
    nlohmann::json config;
    std::ifstream file(argv[2]);
    file >> config;
    const auto voltages =
        config.at("voltage_inputs").get<std::array<std::string, 3>>();
    const auto currents =
        config.at("current_outputs").get<std::array<std::string, 3>>();
    const double omega = config.at("angular_frequency_rad_s");
    const double amplitude = config.at("voltage_peak_V");
    const auto voltage = [=](double t) -> Matrix {
      Matrix v(3, 1);
      for (int p = 0; p < 3; ++p) {
        const double phase = p == 0 ? 0 : (p == 1 ? -2 * PI / 3 : 2 * PI / 3);
        v(p, 0) = amplitude * std::cos(omega * t + phase);
      }
      return v;
    };
    std::vector<FMI::ScalarInput> inputs;
    const Matrix initial = voltage(0);
    for (int p = 0; p < 3; ++p)
      inputs.push_back({voltages[p], initial(p, 0)});
    FMI::ModelExchangeOptions runtimeOptions;
    runtimeOptions.endpointStateEvents =
        config.value("endpoint_state_events", false);
    runtimeOptions.allowContinuousStateResets =
        config.value("allow_continuous_state_resets", false);
    auto model =
        std::make_shared<FMI::ModelExchange>(argv[1], inputs, runtimeOptions);
    auto branch = FMI::NonlinearSSN::make(
        "imported", model,
        FMI::ElectricalPort3Ph{
            voltages, currents,
            config.at("current_sign_to_passive").get<double>()});
    auto source = VoltageBoundary3Ph::make("source", voltage);
    auto node = SimNode<Real>::make("bus", PhaseType::ABC);
    branch->connect({SimNode<Real>::GND, node});
    source->connect({SimNode<Real>::GND, node});
    auto logger = DataLogger::make(String("fmi_prescribed"));
    for (const auto &name : model->description().outputNames)
      logger->logAttribute(name, branch->attribute("fmi.output." + name));
    logger->logAttribute("source_current", source->attribute("i_intf"));
    logger->logAttribute("newton_iterations",
                         branch->attribute("newton_iterations"));
    logger->logAttribute("newton_residual",
                         branch->attribute("newton_residual"));
    SystemTopology topology(omega / (2 * PI), {node}, {source, branch});
    Simulation simulation("fmi_prescribed", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(topology);
    const double step =
        argc == 4 ? std::stod(argv[3]) : config.at("step_s").get<double>();
    const double stop = config.at("stop_time_s").get<double>();
    if (!std::isfinite(step) || step <= 0 || !std::isfinite(stop) ||
        stop <= 0 || std::abs(stop / step - std::round(stop / step)) > 1e-7)
      throw std::invalid_argument(
          "Stop time must be a positive integer number of steps");
    simulation.setTimeStep(step);
    simulation.setFinalTime(stop);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    // Integer step count avoids accumulated-time rounding omitting the last sample.
    simulation.start();
    for (long long k = 0; k < std::llround(stop / step); ++k)
      simulation.step();
    simulation.stop();
    std::cout << "Imported " << model->description().initialState.rows()
              << " continuous states.\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
