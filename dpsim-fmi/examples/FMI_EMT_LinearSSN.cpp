// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary.h"
#include <cmath>
#include <dpsim-fmi/LinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
  if (argc != 2 && argc != 4) {
    std::cerr << "Usage: FMI_EMT_LinearSSN model.fmu [voltage_input "
                 "current_output]\n";
    return 2;
  }
  try {
    using namespace CPS;
    using namespace DPsim;
    const FMI::VoltagePort port{argc == 4 ? argv[2] : "u",
                                argc == 4 ? argv[3] : "i"};
    const auto model = FMI::LinearModel::load(argv[1], port);
    std::cout << "Imported " << model->modelName << " with " << model->A.rows()
              << " states. Matrices recovered from FMI evaluations:\n"
              << "A:\n"
              << model->A << "\nB:\n"
              << model->B << "\nC (all outputs):\n"
              << model->C << "\nD:\n"
              << model->D << '\n';
    // The FMU has a voltage input, so excitation belongs to the surrounding grid.
    // No R/L/C values or dedicated resistor/inductor/capacitor classes occur here.
    auto node = SimNode<Real>::make("terminal", PhaseType::Single);
    auto ground = SimNode<Real>::GND;
    auto branch = FMI::LinearSSN::make("imported", model);
    branch->connect({ground, node});
    auto boundary = VoltageBoundary::make("voltage", [](Real t) {
      return std::sqrt(2.0) * 230.0 * std::sin(2 * PI * 50.0 * t);
    });
    boundary->connect({ground, node});
    SystemTopology system(50, {node}, {boundary, branch});
    auto logger = DataLogger::make(String("fmi_rlc"));
    logger->logAttribute("u_V", branch->attribute("v_intf"));
    for (const auto &output : model->outputNames)
      logger->logAttribute(output, branch->attribute("fmi.output." + output));
    logger->logAttribute("x", branch->attribute("x"));
    Simulation simulation("fmi_rlc");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(1e-5);
    simulation.setFinalTime(0.1);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed inside DPsim MNA; results: logs/fmi_rlc.csv\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
