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
  if (argc != 2 && argc != 8 && argc != 11) {
    std::cerr << "Usage: FMI_EMT_Ph3_CoSimulation model.fmu [u_a u_b u_c "
                 "i_a i_b i_c [current_sign stop_time_s step_s]]\n";
    return 2;
  }
  try {
    const Real frequency = 50.0;
    const Real phaseVoltagePeak = std::sqrt(2.0) * 230.0;
    const Real finalTime = argc == 11 ? std::stod(argv[9]) : 0.1;
    const Real timeStep = argc == 11 ? std::stod(argv[10]) : 10e-6;
    const Real lineResistance = 1.0;
    const Real lineInductance = 2e-3;

    FMI::CoSimulationElectricalPort port{
        {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, -1};
    if (argc >= 8) {
      for (int phase = 0; phase < 3; ++phase) {
        port.voltageInputs[phase] = argv[2 + phase];
        port.currentOutputs[phase] = argv[5 + phase];
      }
    }
    if (argc == 11)
      port.currentSignToPassive = std::stod(argv[8]);

    auto model = std::make_shared<FMI::CoSimulation>(argv[1]);
    auto node = SimNode<Real>::make("pcc", PhaseType::ABC);
    auto gridNode = SimNode<Real>::make("grid", PhaseType::ABC);
    auto lineNode = SimNode<Real>::make("line", PhaseType::ABC);
    auto fmu = FMI::CoSimulationSSN::make("fmu", model, port);
    auto source = EMT::Ph3::VoltageSource::make("grid_source");
    MatrixComp phasors(3, 1);
    for (int phase = 0; phase < 3; ++phase)
      phasors(phase, 0) =
          std::polar(phaseVoltagePeak / RMS3PH_TO_PEAK1PH, -phase * 2 * PI / 3);
    source->setParameters(phasors, frequency);
    auto resistor = EMT::Ph3::Resistor::make("grid_resistor");
    resistor->setParameters(lineResistance * Matrix::Identity(3, 3));
    auto inductor = EMT::Ph3::Inductor::make("grid_inductor");
    inductor->setParameters(lineInductance * Matrix::Identity(3, 3));

    source->connect({SimNode<Real>::GND, gridNode});
    inductor->connect({lineNode, node});
    fmu->connect(node);
    resistor->connect({gridNode, lineNode});
    SystemTopology system(frequency, {gridNode, lineNode, node},
                          {source, resistor, inductor, fmu});

    auto logger = DataLogger::make(String("fmi_cosimulation_3ph"));
    logger->logAttribute("u_pcc", fmu->attribute("v_intf"));
    logger->logAttribute("i_fmu", fmu->attribute("i_intf"));
    logger->logAttribute("fmi_do_steps", fmu->attribute("cosim_do_steps"));

    Simulation simulation("fmi_cosimulation_3ph");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed explicit FMI Co-Simulation; "
                 "logs/fmi_cosimulation_3ph.csv\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
