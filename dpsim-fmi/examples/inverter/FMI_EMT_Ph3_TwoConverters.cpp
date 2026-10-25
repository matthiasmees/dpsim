// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "Usage: FMI_EMT_Ph3_TwoConverters first.fmu second.fmu "
                 "[step_s [stop_time_s]]\n";
    return 2;
  }

  try {
    const double frequency = 50.0;
    const double phaseVoltagePeak = 325.269119;
    const double timeStep = argc >= 4 ? std::stod(argv[3]) : 50e-6;
    const double finalTime = argc >= 5 ? std::stod(argv[4]) : 5e-3;
    if (timeStep <= 0 || finalTime <= 0 ||
        std::abs(finalTime / timeStep - std::round(finalTime / timeStep)) >
            1e-8)
      throw std::invalid_argument(
          "step and stop time must form an integer grid");

    // The two paths are intentionally separate FMI instances. This is a
    // useful smoke test for shared-network coupling and does not clone an FMU.
    FMI::ModelExchangeOptions modelOptions;
    modelOptions.endpointStateEvents = true;
    FMI::NewtonSettings settings;
    settings.multirateMode = FMI::MultirateMode::SinglePassExplicit;
    auto firstConverter =
        FMI::importFMU("first_converter", argv[1], modelOptions, settings);
    auto secondConverter =
        FMI::importFMU("second_converter", argv[2], modelOptions, settings);

    auto gridSource = EMT::Ph3::VoltageSource::make("grid_source");
    MatrixComp phasors(3, 1);
    for (int phase = 0; phase < 3; ++phase)
      phasors(phase, 0) =
          std::polar(phaseVoltagePeak / RMS3PH_TO_PEAK1PH, -phase * 2 * PI / 3);
    std::vector<Complex> initialNodeVoltage(3);
    for (int phase = 0; phase < 3; ++phase)
      initialNodeVoltage[phase] = phasors(phase, 0).real();
    auto gridNode =
        SimNode<Real>::make("grid", PhaseType::ABC, initialNodeVoltage);
    auto pccNode =
        SimNode<Real>::make("pcc", PhaseType::ABC, initialNodeVoltage);
    gridSource->setParameters(phasors, frequency);
    auto lineResistor = EMT::Ph3::Resistor::make("line_resistor");
    lineResistor->setParameters(0.5 * Matrix::Identity(3, 3));
    auto nativeLoad = EMT::Ph3::Resistor::make("native_load");
    nativeLoad->setParameters(20.0 * Matrix::Identity(3, 3));

    gridSource->connect({SimNode<Real>::GND, gridNode});
    lineResistor->connect({gridNode, pccNode});
    nativeLoad->connect({SimNode<Real>::GND, pccNode});
    firstConverter->connect(pccNode);
    secondConverter->connect(pccNode);
    SystemTopology system(frequency, {gridNode, pccNode},
                          {gridSource, lineResistor, nativeLoad, firstConverter,
                           secondConverter});

    auto logger = DataLogger::make(String("fmi_two_converters"));
    logger->logAttribute("pcc_voltage", pccNode->attribute("v"));
    logger->logAttribute("first_current", firstConverter->attribute("i_intf"));
    logger->logAttribute("second_current",
                         secondConverter->attribute("i_intf"));
    logger->logAttribute("first_residual",
                         firstConverter->attribute("newton_residual"));
    logger->logAttribute("second_residual",
                         secondConverter->attribute("newton_residual"));
    logger->logAttribute("first_internal_steps",
                         firstConverter->attribute("fmi_internal_steps"));
    logger->logAttribute("second_internal_steps",
                         secondConverter->attribute("fmi_internal_steps"));

    Simulation simulation("fmi_two_converters", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();

    std::cout << "Two FMI converters completed on one native EMT network: "
              << std::llround(finalTime / timeStep) << " steps\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
