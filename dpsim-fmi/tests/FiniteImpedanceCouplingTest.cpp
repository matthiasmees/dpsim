// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_SeriesResistor.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace CPS;
using namespace DPsim;

int main(int argc, char **argv) {
  if (argc != 4)
    return 2;
  try {
    const std::string mode = argv[2];
    const Real step = std::stod(argv[3]);
    const bool explicitPass = mode == "single_pass";
    if (!explicitPass && mode != "rollback_newton")
      throw std::invalid_argument(
          "mode must be rollback_newton or single_pass");

    constexpr Real frequency = 50.0;
    constexpr Real finalTime = 0.003;
    constexpr Real gridResistance = 5.0;
    constexpr Real gridInductance = 0.005;
    const auto sourceVoltage = [=](Real time) {
      Matrix voltage(3, 1);
      for (int phase = 0; phase < 3; ++phase) {
        const Real angle =
            phase == 0 ? 0 : (phase == 1 ? -2 * PI / 3 : 2 * PI / 3);
        voltage(phase, 0) =
            std::sqrt(2.0) * 23.0 * std::sin(2 * PI * frequency * time + angle);
      }
      return voltage;
    };

    std::vector<Complex> zero(3, Complex(0, 0));
    auto grid = SimNode<Real>::make("grid", PhaseType::ABC, zero);
    auto middle = SimNode<Real>::make("rl_mid", PhaseType::ABC, zero);
    auto pcc = SimNode<Real>::make("pcc", PhaseType::ABC, zero);
    auto source = VoltageBoundary3Ph::make("grid_source", sourceVoltage);
    auto resistor = EMT::Ph3::SeriesResistor::make("grid_R");
    resistor->setParameters(gridResistance);
    auto inductor = EMT::Ph3::Inductor::make("grid_L");
    inductor->setParameters(Matrix::Identity(3, 3) * gridInductance);

    FMI::NewtonSettings settings;
    if (explicitPass)
      settings.multirateMode = FMI::MultirateMode::SinglePassExplicit;
    auto fmu = std::make_shared<FMI::ModelExchange>(argv[1]);
    auto converter = FMI::NonlinearSSN::make(
        "nonlinear_fmu", fmu,
        FMI::ElectricalPort3Ph{{"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, 1},
        settings);

    source->connect({SimNode<Real>::GND, grid});
    resistor->connect({grid, middle});
    inductor->connect({middle, pcc});
    converter->connect({SimNode<Real>::GND, pcc});
    SystemTopology system(frequency, {grid, middle, pcc},
                          {source, resistor, inductor, converter});
    Simulation simulation("finite_impedance_coupling", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(step);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);

    auto logger = DataLogger::make(String("finite_impedance_" + mode));
    logger->logAttribute("pcc_voltage", pcc->attribute("v"));
    logger->logAttribute("converter_current", converter->attribute("i_intf"));
    logger->logAttribute("newton_iterations",
                         converter->attribute("newton_iterations"));
    logger->logAttribute("newton_residual",
                         converter->attribute("newton_residual"));
    logger->logAttribute("fmi_coupling_error_A",
                         converter->attribute("fmi_coupling_error_A"));
    logger->logAttribute("fmi_internal_steps",
                         converter->attribute("fmi_internal_steps"));
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "mode=" << mode << " step_s=" << step
              << " R_grid_ohm=" << gridResistance
              << " L_grid_H=" << gridInductance << " final_time_s=" << finalTime
              << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
