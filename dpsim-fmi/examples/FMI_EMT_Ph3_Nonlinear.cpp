// SPDX-License-Identifier: MPL-2.0
#include "LinearNetwork3Ph.h"
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "Usage: FMI_EMT_Ph3_Nonlinear model.fmu benchmark.json "
                 "[step_s [tolerance_multiplier]]\n";
    return 2;
  }
  try {
    using namespace CPS;
    using namespace DPsim;
    nlohmann::json config;
    std::ifstream input(argv[2]);
    input >> config;
    const auto voltages =
        config.at("voltage_inputs").get<std::array<std::string, 3>>();
    const auto currents =
        config.at("current_outputs").get<std::array<std::string, 3>>();
    const auto controlNames =
        config.at("control_names").get<std::vector<std::string>>();
    const auto controlValues =
        config.at("control_initial_values").get<std::vector<double>>();
    const auto initial =
        config.at("initial_network_state").get<std::vector<double>>();
    const auto initialI =
        config.at("initial_generator_current").get<std::vector<double>>();
    if (initial.size() != 6 || initialI.size() != 3 ||
        controlNames.size() != controlValues.size())
      throw std::invalid_argument("Invalid benchmark configuration dimensions");
    Matrix x0(6, 1), i0(3, 1);
    std::vector<Complex> node0(3);
    std::vector<FMI::ScalarInput> inputs;
    for (int k = 0; k < 6; ++k)
      x0(k, 0) = initial[k];
    for (int p = 0; p < 3; ++p) {
      i0(p, 0) = initialI[p];
      node0[p] = initial[p] / RMS3PH_TO_PEAK1PH;
      inputs.push_back({voltages[p], initial[p]});
    }
    for (std::size_t k = 0; k < controlNames.size(); ++k)
      inputs.push_back({controlNames[k], controlValues[k]});
    auto model = std::make_shared<FMI::ModelExchange>(argv[1], inputs);
    FMI::NewtonSettings settings;
    if (argc == 5) {
      const double factor = std::stod(argv[4]);
      settings.absoluteStateTolerance *= factor;
      settings.relativeStateTolerance *= factor;
      settings.absoluteCurrentTolerance *= factor;
      settings.relativeCurrentTolerance *= factor;
    }
    auto generator = FMI::NonlinearSSN::make(
        "imported", model,
        FMI::ElectricalPort3Ph{
            voltages, currents,
            config.at("current_sign_to_passive").get<double>()},
        settings);
    const double event = config.at("vref_step_s"),
                 increment = config.at("vref_step_pu");
    generator->setInputFunction(
        "vref_pu", [=](double t) { return 1 + increment * (t >= event); },
        {event});
    const double frequency = config.at("frequency_Hz"),
                 vrms = config.at("phase_rms_V");
    const double G0 = config.at("G0_S"), Gstep = config.at("Gstep_S"),
                 loadTime = config.at("load_step_s");
    const auto grid = [=](double time) -> Matrix {
      Matrix v(3, 1);
      const double a = 2 * PI * frequency * time;
      v << std::sqrt(2.0) * vrms * std::cos(a),
          std::sqrt(2.0) * vrms * std::cos(a - 2 * PI / 3),
          std::sqrt(2.0) * vrms * std::cos(a + 2 * PI / 3);
      return v;
    };
    const auto load = [=](double t) {
      return G0 + Gstep * (t >= loadTime - 1e-12);
    };
    auto network = LinearNetwork3Ph::make(
        "external_grid", config.at("network")[0].get<double>(),
        config.at("network")[1].get<double>(),
        config.at("network")[2].get<double>(), x0, i0, grid, load);
    auto bus = SimNode<Real>::make("bus", PhaseType::ABC, node0);
    generator->connect({SimNode<Real>::GND, bus});
    network->connect({SimNode<Real>::GND, bus});
    auto logger = DataLogger::make(String("fmi_sg_nonlinear"));
    for (auto &name : model->description().outputNames)
      logger->logAttribute(name, generator->attribute("fmi.output." + name));
    for (auto &name : model->description().inputNames)
      logger->logAttribute("input_" + name,
                           generator->attribute("fmi.input." + name));
    logger->logAttribute("newton_iterations",
                         generator->attribute("newton_iterations"));
    logger->logAttribute("newton_residual",
                         generator->attribute("newton_residual"));
    logger->logAttribute("network_current", network->attribute("i_intf"));
    for (const std::string name :
         {"v_bus_a_V", "v_bus_b_V", "v_bus_c_V", "i_line_a_A", "i_line_b_A",
          "i_line_c_A", "g_load_S", "p_load_W", "p_line_bus_W",
          "p_line_loss_W"})
      logger->logAttribute(name, network->attribute(name));
    const double step =
        argc >= 4 ? std::stod(argv[3]) : config.at("step_s").get<double>();
    if (std::abs(loadTime / step - std::round(loadTime / step)) > 1e-7)
      throw std::invalid_argument("Load event must align with h");
    SystemTopology system(frequency, {bus}, {generator, network});
    Simulation simulation("fmi_sg_nonlinear", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(step);
    simulation.setFinalTime(config.at("stop_time_s").get<double>());
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed nonlinear FMI/MNA trapezoidal simulation with "
              << model->description().initialState.rows()
              << " imported states, h=" << step << ".\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
