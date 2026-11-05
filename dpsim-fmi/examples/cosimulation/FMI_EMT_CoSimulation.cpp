// SPDX-License-Identifier: MPL-2.0
#include "../common/VoltageBoundary.h"
#include <cmath>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc != 2 && argc != 4 && argc != 5) {
    std::cerr << "Usage: FMI_EMT_CoSimulation model.fmu [voltage_input "
                 "current_output [stop_time_s]]\n";
    return 2;
  }
  try {
    const Real frequency = 50.0;
    const Real phaseVoltageRms = 230.0;
    const Real timeStep = 10e-6;
    const Real finalTime = argc == 5 ? std::stod(argv[4]) : 0.1;
    FMI::CoSimulationElectricalPort port{
        {argc >= 4 ? argv[2] : "u"}, {argc >= 4 ? argv[3] : "i"}, -1};
    auto model = std::make_shared<FMI::CoSimulation>(argv[1]);
    auto node = SimNode<Real>::make("terminal", PhaseType::Single);
    auto fmu = FMI::CoSimulationSSN::make("imported", model, port);
    auto source = VoltageBoundary::make("voltage", [=](Real time) {
      return std::sqrt(2.0) * phaseVoltageRms *
             std::sin(2 * PI * frequency * time);
    });
    fmu->connect(node);
    source->connect({SimNode<Real>::GND, node});
    SystemTopology system(frequency, {node}, {source, fmu});

    auto logger = DataLogger::make(String("fmi_cosimulation"));
    logger->logAttribute("u_V", fmu->attribute("v_intf"));
    logger->logAttribute("i_A", fmu->attribute("i_intf"));
    logger->logAttribute("fmi_do_steps", fmu->attribute("cosim_do_steps"));
    Simulation simulation("fmi_cosimulation");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout
        << "Completed explicit FMI Co-Simulation; logs/fmi_cosimulation.csv\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
