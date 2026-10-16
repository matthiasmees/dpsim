// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <cmath>
#include <dpsim-fmi/LinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>

int main(int argc, char **argv) {
  if (argc != 2 && argc != 8) {
    std::cerr
        << "Usage: FMI_EMT_Ph3_LinearSSN model.fmu [u_a u_b u_c i_a i_b i_c]\n";
    return 2;
  }
  try {
    using namespace CPS;
    using namespace DPsim;
    // Excitation describes the surrounding grid. All branch equations and
    // the state count come from the FMU, including any cross-phase coupling.
    const auto voltage = [](Real t) -> Matrix {
      Matrix u(3, 1);
      for (int p = 0; p < 3; ++p) {
        const Real angle = p == 1 ? -2 * PI / 3 : p == 2 ? 2 * PI / 3 : 0;
        u(p, 0) = std::sqrt(2.0) * 230 * std::sin(2 * PI * 50 * t + angle);
      }
      return u;
    };
    FMI::VoltagePorts ports{{"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, {}};
    if (argc == 8)
      for (int p = 0; p < 3; ++p) {
        ports.voltageInputs[p] = argv[2 + p];
        ports.currentOutputs[p] = argv[5 + p];
      }
    const Matrix u0 = voltage(0);
    ports.initialVoltages.assign(u0.data(), u0.data() + u0.size());
    const auto model = FMI::LinearModel::loadMulti(argv[1], ports);
    std::cout << "Imported " << model->modelName << " with " << model->A.rows()
              << " states and three electrical channels.\nA:\n"
              << model->A << "\nB:\n"
              << model->B << "\nC:\n"
              << model->C << "\nD:\n"
              << model->D << '\n';
    std::vector<Complex> initial(3);
    for (int p = 0; p < 3; ++p)
      initial[p] = u0(p, 0) / RMS3PH_TO_PEAK1PH;
    auto node = SimNode<Real>::make("terminal", PhaseType::ABC, initial);
    auto branch = FMI::LinearSSN::make("imported", model);
    branch->connect({SimNode<Real>::GND, node});
    auto source = VoltageBoundary3Ph::make("voltage", voltage);
    source->connect({SimNode<Real>::GND, node});
    SystemTopology system(50, {node}, {source, branch});
    auto logger = DataLogger::make(String("fmi_rlc_3ph"));
    for (const auto &input : model->voltageInputs)
      logger->logAttribute(input, branch->attribute("fmi.input." + input));
    for (const auto &output : model->outputNames)
      logger->logAttribute(output, branch->attribute("fmi.output." + output));
    logger->logAttribute("x", branch->attribute("x"));
    logger->logAttribute("source_current", source->attribute("i_intf"));
    Simulation simulation("fmi_rlc_3ph");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(1e-5);
    simulation.setFinalTime(0.1);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed inside DPsim MNA: logs/fmi_rlc_3ph.csv\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
