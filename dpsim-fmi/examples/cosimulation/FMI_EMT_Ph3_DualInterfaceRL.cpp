// SPDX-License-Identifier: MPL-2.0
#include "../common/VoltageBoundary3Ph.h"
#include <algorithm>
#include <cmath>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>
#include <memory>

using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cerr << "Usage: FMI_EMT_Ph3_DualInterfaceRL me|cs model.fmu "
                 "stop_time_s step_s log_downsampling\n";
    return 2;
  }
  try {
    const std::string mode = argv[1];
    if (mode != "me" && mode != "cs")
      throw std::invalid_argument("mode must be 'me' or 'cs'");
    const Real stopTime = std::stod(argv[3]);
    const Real timeStep = std::stod(argv[4]);
    const int downsampling = std::stoi(argv[5]);
    if (stopTime <= 0 || timeStep <= 0 || downsampling <= 0 ||
        std::abs(stopTime / timeStep - std::round(stopTime / timeStep)) > 1e-8)
      throw std::invalid_argument("invalid simulation time grid");

    constexpr Real frequency = 50;
    const Real phasePeak = std::sqrt(2.0) * 230.0;
    constexpr Real sagStart = 0.04;
    constexpr Real sagEnd = 0.06;
    const auto gridVoltage = [=](Real time) {
      Matrix v(3, 1);
      const Real envelope = std::min(1.0, std::max(0.0, time / 0.01));
      const Real sagScale = time >= sagStart && time < sagEnd ? 0.75 : 1.0;
      const Real scale = envelope * sagScale;
      for (int phase = 0; phase < 3; ++phase)
        v(phase, 0) = scale * phasePeak *
                      std::sin(2 * PI * frequency * time -
                               phase * 2 * PI / 3);
      return v;
    };

    const Matrix initialVoltage = gridVoltage(0);
    std::vector<Complex> initialNodeVoltage(3);
    for (int phase = 0; phase < 3; ++phase)
      initialNodeVoltage[phase] =
          initialVoltage(phase, 0) / RMS3PH_TO_PEAK1PH;
    auto sourceNode = SimNode<Real>::make("source", PhaseType::ABC,
                                          initialNodeVoltage);
    auto lineNode = SimNode<Real>::make("line", PhaseType::ABC,
                                        initialNodeVoltage);
    auto pccNode = SimNode<Real>::make("pcc", PhaseType::ABC,
                                       initialNodeVoltage);

    auto source = VoltageBoundary3Ph::make("grid_voltage", gridVoltage);
    auto sourceResistor = EMT::Ph3::Resistor::make("grid_resistor");
    sourceResistor->setParameters(2.0 * Matrix::Identity(3, 3));
    auto sourceInductor = EMT::Ph3::Inductor::make("grid_inductor");
    sourceInductor->setParameters(0.01 * Matrix::Identity(3, 3));
    auto nativeLoad = EMT::Ph3::Resistor::make("native_load");
    nativeLoad->setParameters(100.0 * Matrix::Identity(3, 3));

    auto fmu = [&]() -> std::shared_ptr<SimPowerComp<Real>> {
      if (mode == "me") {
        auto component = FMI::NonlinearSSN::importFMU(std::string(argv[2]));
        component->connect(pccNode);
        return component;
      }
      FMI::CoSimulationElectricalPort cosimPort{
          {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, 1.0};
      auto component = FMI::CoSimulationSSN::importFMU(
          "fmu", argv[2], cosimPort);
      component->connect(pccNode);
      return component;
    }();

    source->connect({SimNode<Real>::GND, sourceNode});
    sourceResistor->connect({sourceNode, lineNode});
    sourceInductor->connect({lineNode, pccNode});
    nativeLoad->connect({SimNode<Real>::GND, pccNode});
    SystemTopology system(frequency, {sourceNode, lineNode, pccNode},
                          {source, sourceResistor, sourceInductor, nativeLoad,
                           fmu});

    auto logger = DataLogger::make(String("dual_interface_rl_" + mode), true,
                                   static_cast<UInt>(downsampling));
    logger->logAttribute("pcc_v", fmu->attribute("v_intf"));
    logger->logAttribute("fmu_i", fmu->attribute("i_intf"));
    logger->logAttribute("source_i", source->attribute("i_intf"));
    if (mode == "me") {
      logger->logAttribute("newton_iterations",
                           fmu->attribute("newton_iterations"));
      logger->logAttribute("newton_residual", fmu->attribute("newton_residual"));
    } else {
      logger->logAttribute("fmi_do_steps", fmu->attribute("cosim_do_steps"));
    }

    Simulation simulation("dual_interface_rl_" + mode);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(stopTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.run();
    std::cout << "Completed " << mode << " run; logs/dual_interface_rl_"
              << mode << ".csv\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
