// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <cmath>
#include <dpsim-fmi/LinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc != 2 && argc != 8 && argc != 9) {
    std::cerr
        << "Usage: FMI_EMT_Ph3_LinearSSN model.fmu [u_a u_b u_c i_a i_b i_c "
           "stop_time_s]\n";
    return 2;
  }
  try {

    // ----- PARAMETERS AND FMU INTERFACE -----

    const Real systemFrequency = 50.0;
    const Real phaseVoltageRms = 230.0;
    const Real timeStep = 10e-6;
    const Real finalTime = argc == 9 ? std::stod(argv[8]) : 0.1;

    // Excitation describes the surrounding grid. All branch equations and
    // the state count come from the FMU, including any cross-phase coupling.
    const auto voltage = [=](Real time) -> Matrix {
      Matrix phaseVoltages(3, 1);
      for (int phase = 0; phase < 3; ++phase) {
        const Real angle = phase == 1   ? -2 * PI / 3
                           : phase == 2 ? 2 * PI / 3
                                        : 0;
        phaseVoltages(phase, 0) =
            std::sqrt(2.0) * phaseVoltageRms *
            std::sin(2 * PI * systemFrequency * time + angle);
      }
      return phaseVoltages;
    };
    FMI::VoltagePorts ports{{"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, {}};
    if (argc == 8) {
      for (int phase = 0; phase < 3; ++phase) {
        ports.voltageInputs[phase] = argv[2 + phase];
        ports.currentOutputs[phase] = argv[5 + phase];
      }
    }
    const Matrix initialVoltage = voltage(0);
    ports.initialVoltages.assign(initialVoltage.data(),
                                 initialVoltage.data() + initialVoltage.size());
    const auto model = FMI::LinearModel::loadMulti(argv[1], ports);
    std::cout << "Imported " << model->modelName << " with " << model->A.rows()
              << " states and three electrical channels.\nA:\n"
              << model->A << "\nB:\n"
              << model->B << "\nC:\n"
              << model->C << "\nD:\n"
              << model->D << '\n';

    // ----- INITIAL VALUES, NODES AND COMPONENTS -----

    std::vector<Complex> initialNodeVoltages(3);
    for (int phase = 0; phase < 3; ++phase) {
      initialNodeVoltages[phase] = initialVoltage(phase, 0) / RMS3PH_TO_PEAK1PH;
    }
    auto node =
        SimNode<Real>::make("terminal", PhaseType::ABC, initialNodeVoltages);
    auto branch = FMI::LinearSSN::make("imported", model);
    auto source = VoltageBoundary3Ph::make("voltage", voltage);

    // ----- TOPOLOGY -----

    branch->connect(node);
    source->connect({SimNode<Real>::GND, node});
    SystemTopology system(systemFrequency, {node}, {source, branch});

    // ----- LOGGING -----

    auto logger = DataLogger::make(String("fmi_rlc_3ph"));
    for (const auto &input : model->voltageInputs) {
      logger->logAttribute(input, branch->attribute("fmi.input." + input));
    }
    for (const auto &output : model->outputNames) {
      logger->logAttribute(output, branch->attribute("fmi.output." + output));
    }
    logger->logAttribute("x", branch->attribute("x"));
    logger->logAttribute("source_current", source->attribute("i_intf"));

    // ----- SIMULATION -----

    Simulation simulation("fmi_rlc_3ph");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed inside DPsim MNA: logs/fmi_rlc_3ph.csv\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
