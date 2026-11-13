// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary.h"
#include <cmath>
#include <dpsim-fmi/LinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <fstream>
#include <iostream>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc != 2 && argc != 4 && argc != 5) {
    std::cerr << "Usage: FMI_EMT_LinearSSN model.fmu [voltage_input "
                 "current_output [stop_time_s]]\n";
    return 2;
  }
  try {

    // ----- PARAMETERS AND FMU INTERFACE -----

    const Real systemFrequency = 50.0;
    const Real phaseVoltageRms = 230.0;
    const Real timeStep = 10e-6;
    const Real finalTime = argc == 5 ? std::stod(argv[4]) : 0.1;

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

    // ----- NODES AND COMPONENTS -----

    // The FMU has a voltage input, so excitation belongs to the surrounding
    // grid. No R/L/C values or dedicated resistor/inductor/capacitor classes
    // occur here.
    auto node = SimNode<Real>::make("terminal", PhaseType::Single);
    auto ground = SimNode<Real>::GND;
    auto branch = FMI::LinearSSN::make("imported", model);
    auto boundary = VoltageBoundary::make("voltage", [=](Real time) {
      return std::sqrt(2.0) * phaseVoltageRms *
             std::sin(2 * PI * systemFrequency * time);
    });

    // ----- TOPOLOGY -----

    branch->connect(node);
    boundary->connect({ground, node});
    SystemTopology system(systemFrequency, {node}, {boundary, branch});

    // ----- LOGGING -----

    auto logger = DataLogger::make(String("fmi_rlc"));
    logger->logAttribute("u_V", branch->attribute("v_intf"));
    for (const auto &output : model->outputNames) {
      logger->logAttribute(output, branch->attribute("fmi.output." + output));
    }
    logger->logAttribute("x", branch->attribute("x"));

    // ----- SIMULATION -----

    Simulation simulation("fmi_rlc");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
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
