// SPDX-License-Identifier: MPL-2.0
#include "../common/VoltageBoundary.h"
#include <DPsim.h>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim-fmi/LinearSSN.h>
#include <dpsim-models/EMT/EMT_Ph1_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph1_Resistor.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>

using namespace DPsim;
using namespace CPS;

namespace {
constexpr Real kSourceVoltage = 10.0;
constexpr Real kGridResistance = 1.0;
constexpr Real kFmuResistance = 0.1;
constexpr Real kInductance = 20e-6;
constexpr Real kFrequency = 50.0;

void run(const String &mode, const String &fmuPath, Real timeStep,
         Real finalTime) {
  auto sourceNode = SimNode<Real>::make("source_node", PhaseType::Single);
  auto pcc = SimNode<Real>::make("pcc", PhaseType::Single);
  auto source =
      VoltageBoundary::make("source", [](Real) { return kSourceVoltage; });
  auto gridResistance = EMT::Ph1::Resistor::make("Rg");
  gridResistance->setParameters(kGridResistance);
  source->connect({SimNode<Real>::GND, sourceNode});
  gridResistance->connect({sourceNode, pcc});

  SystemComponentList components{source, gridResistance};
  SystemNodeList nodes{sourceNode, pcc};
  AttributeBase::Ptr current;
  AttributeBase::Ptr voltage = pcc->attribute("v");
  std::shared_ptr<FMI::CoSimulationSSN> coSimulation;

  if (mode == "native") {
    auto fmuResistance = EMT::Ph1::Resistor::make("Rf");
    fmuResistance->setParameters(kFmuResistance);
    auto inductance = EMT::Ph1::Inductor::make("L");
    inductance->setParameters(kInductance);
    auto inductorNode = SimNode<Real>::make("inductor_node", PhaseType::Single);
    fmuResistance->connect({pcc, inductorNode});
    // Orient the reported inductor current from source to ground, matching
    // the FMI passive current convention at its positive terminal.
    inductance->connect({SimNode<Real>::GND, inductorNode});
    current = inductance->attribute("i_intf");
    components.push_back(fmuResistance);
    components.push_back(inductance);
    nodes.push_back(inductorNode);
  } else if (mode == "me") {
    auto model = FMI::LinearModel::load(fmuPath, FMI::VoltagePort{"u", "i"});
    auto branch = FMI::LinearSSN::make("linear_rl_fmu", model);
    branch->connect(pcc);
    current = branch->attribute("i_intf");
    components.push_back(branch);
  } else if (mode == "cs") {
    auto model = std::make_shared<FMI::CoSimulation>(fmuPath);
    FMI::CoSimulationElectricalPort port{{"u"}, {"i"}, 1};
    coSimulation = FMI::CoSimulationSSN::make("linear_rl_fmu", model, port);
    coSimulation->connect(pcc);
    current = coSimulation->attribute("i_intf");
    components.push_back(coSimulation);
  } else {
    throw std::invalid_argument("mode must be native, me, or cs");
  }

  SystemTopology system(kFrequency, nodes, components);
  auto logger = DataLogger::make(String("linear_rl_" + mode));
  logger->logAttribute("u_V", voltage);
  logger->logAttribute("i_A", current);
  if (coSimulation)
    logger->logAttribute("fmi_do_steps",
                         coSimulation->attribute("cosim_do_steps"));

  Simulation simulation("linear_rl_" + mode, Logger::Level::warn);
  simulation.setDomain(Domain::EMT);
  simulation.setSystem(system);
  simulation.setTimeStep(timeStep);
  simulation.setFinalTime(finalTime);
  simulation.doSteadyStateInit(false);
  simulation.addLogger(logger);
  simulation.run();
  std::cout << "Completed " << mode << " at dt=" << timeStep
            << " s; results in logs/linear_rl_" << mode << ".csv\n";
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 5 ||
      (std::string(argv[1]) != "native" && std::string(argv[1]) != "me" &&
       std::string(argv[1]) != "cs")) {
    std::cerr << "Usage: FMI_LinearRLStability native|me|cs model.fmu|- "
                 "step_s stop_time_s\n";
    return 2;
  }
  try {
    run(argv[1], argv[2], std::stod(argv[3]), std::stod(argv[4]));
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
