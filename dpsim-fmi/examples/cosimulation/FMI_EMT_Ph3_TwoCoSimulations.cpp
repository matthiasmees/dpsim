// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "Usage: FMI_EMT_Ph3_TwoCoSimulations first.fmu second.fmu "
                 "[stop_time_s [step_s]]\n";
    return 2;
  }
  try {
    const Real frequency = 50.0;
    const Real phaseVoltagePeak = std::sqrt(2.0) * 230.0;
    const Real finalTime = argc >= 4 ? std::stod(argv[3]) : 0.1;
    const Real timeStep = argc == 5 ? std::stod(argv[4]) : 10e-6;
    const FMI::CoSimulationElectricalPort port{
        {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, -1};
    auto firstModel = std::make_shared<FMI::CoSimulation>(argv[1]);
    auto secondModel = std::make_shared<FMI::CoSimulation>(argv[2]);
    auto first = FMI::CoSimulationSSN::make("first_fmu", firstModel, port);
    auto second = FMI::CoSimulationSSN::make("second_fmu", secondModel, port);

    auto gridNode = SimNode<Real>::make("grid", PhaseType::ABC);
    auto pccNode = SimNode<Real>::make("pcc", PhaseType::ABC);
    auto lineNode = SimNode<Real>::make("line", PhaseType::ABC);
    auto source = EMT::Ph3::VoltageSource::make("grid_source");
    MatrixComp phasors(3, 1);
    for (int phase = 0; phase < 3; ++phase)
      phasors(phase, 0) =
          std::polar(phaseVoltagePeak / RMS3PH_TO_PEAK1PH, -phase * 2 * PI / 3);
    source->setParameters(phasors, frequency);
    auto resistor = EMT::Ph3::Resistor::make("grid_resistor");
    resistor->setParameters(Matrix::Identity(3, 3));
    auto inductor = EMT::Ph3::Inductor::make("grid_inductor");
    inductor->setParameters(2e-3 * Matrix::Identity(3, 3));

    source->connect({SimNode<Real>::GND, gridNode});
    resistor->connect({gridNode, lineNode});
    inductor->connect({lineNode, pccNode});
    first->connect(pccNode);
    second->connect(pccNode);
    SystemTopology system(frequency, {gridNode, lineNode, pccNode},
                          {source, resistor, inductor, first, second});

    auto logger = DataLogger::make(String("fmi_cosimulation_two"));
    logger->logAttribute("pcc_voltage", pccNode->attribute("v"));
    logger->logAttribute("first_current", first->attribute("i_intf"));
    logger->logAttribute("second_current", second->attribute("i_intf"));
    logger->logAttribute("first_do_steps", first->attribute("cosim_do_steps"));
    logger->logAttribute("second_do_steps",
                         second->attribute("cosim_do_steps"));

    Simulation simulation("fmi_cosimulation_two");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed two-FMU explicit FMI Co-Simulation; "
                 "logs/fmi_cosimulation_two.csv\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
