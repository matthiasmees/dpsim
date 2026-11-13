// SPDX-License-Identifier: MPL-2.0
#include "../common/VoltageBoundary3Ph.h"
#include <DPsim.h>
#include <algorithm>
#include <cmath>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim-fmi/LinearSSN.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>

using namespace DPsim;
using namespace CPS;

namespace {
constexpr Real kPhaseVoltageRms = 10.0;
constexpr Real kFrequency = 50.0;
constexpr Real kGridResistance = 1.0;
constexpr Real kFmuResistance = 0.1;
constexpr Real kInductance = 20e-6;

Matrix sourceVoltage(Real time) {
  Matrix voltage(3, 1);
  constexpr Real phaseAngles[] = {0.0, -2.0 * PI / 3.0, 2.0 * PI / 3.0};
  const Real ramp = std::min(1.0, std::max(0.0, time / 1e-3));
  for (int phase = 0; phase < 3; ++phase)
    voltage(phase, 0) =
        ramp * std::sqrt(2.0) * kPhaseVoltageRms *
        std::sin(2.0 * PI * kFrequency * time + phaseAngles[phase]);
  return voltage;
}

void run(const String &mode, const String &fmuPath, Real step, Real finalTime) {
  const Matrix initialVoltage = sourceVoltage(0.0);
  auto sourceNode = SimNode<Real>::make("source_node", PhaseType::ABC);
  auto pcc = SimNode<Real>::make("pcc", PhaseType::ABC);
  auto source = VoltageBoundary3Ph::make("source", sourceVoltage);
  source->connect({SimNode<Real>::GND, sourceNode});

  auto gridResistance = EMT::Ph3::Resistor::make("Rg");
  gridResistance->setParameters(kGridResistance * Matrix::Identity(3, 3));
  gridResistance->connect({sourceNode, pcc});

  SystemComponentList components{source, gridResistance};
  SystemNodeList nodes{sourceNode, pcc};
  AttributeBase::Ptr current;
  std::shared_ptr<FMI::CoSimulationSSN> coSimulation;

  if (mode == "native") {
    auto fmuResistance = EMT::Ph3::Resistor::make("Rf");
    fmuResistance->setParameters(kFmuResistance * Matrix::Identity(3, 3));
    auto inductor = EMT::Ph3::Inductor::make("L");
    inductor->setParameters(kInductance * Matrix::Identity(3, 3));
    auto branchNode = SimNode<Real>::make("branch_node", PhaseType::ABC);
    fmuResistance->connect({pcc, branchNode});
    inductor->connect({SimNode<Real>::GND, branchNode});
    current = inductor->attribute("i_intf");
    components.push_back(fmuResistance);
    components.push_back(inductor);
    nodes.push_back(branchNode);
  } else if (mode == "me") {
    FMI::VoltagePorts ports{{"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, {}};
    ports.initialVoltages.assign(initialVoltage.data(),
                                 initialVoltage.data() + initialVoltage.size());
    auto model = FMI::LinearModel::loadMulti(fmuPath, ports);
    auto branch = FMI::LinearSSN::make("linear_rl_3ph_fmu", model);
    branch->connect(pcc);
    current = branch->attribute("i_intf");
    components.push_back(branch);
  } else if (mode == "cs" || mode == "cs-rollback") {
    FMI::CoSimulationOptions options;
    const bool useRollback = mode == "cs-rollback";
    options.enableRollback = useRollback;
    auto model = std::make_shared<FMI::CoSimulation>(fmuPath, options);
    FMI::CoSimulationElectricalPort port{
        {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, 1};
    const auto coupling =
        useRollback ? FMI::CoSimulationCouplingMode::RollbackNewton
                    : FMI::CoSimulationCouplingMode::SinglePassExplicit;
    coSimulation =
        FMI::CoSimulationSSN::make("linear_rl_3ph_fmu", model, port, coupling);
    coSimulation->connect(pcc);
    current = coSimulation->attribute("i_intf");
    components.push_back(coSimulation);
  } else {
    throw std::invalid_argument("mode must be native, me, cs, or cs-rollback");
  }

  SystemTopology system(kFrequency, nodes, components);
  auto logger = DataLogger::make(String("linear_rl_3ph_" + mode));
  logger->logAttribute("u_abc_V", pcc->attribute("v"));
  logger->logAttribute("i_abc_A", current);
  if (coSimulation)
    logger->logAttribute("fmi_do_steps",
                         coSimulation->attribute("cosim_do_steps"));
  if (coSimulation && mode == "cs-rollback") {
    logger->logAttribute("newton_iterations",
                         coSimulation->attribute("newton_iterations"));
    logger->logAttribute("newton_residual_A",
                         coSimulation->attribute("newton_residual"));
  }

  Simulation simulation("linear_rl_3ph_" + mode, Logger::Level::warn);
  simulation.setDomain(Domain::EMT);
  simulation.setSystem(system);
  simulation.setTimeStep(step);
  simulation.setFinalTime(finalTime);
  simulation.doSteadyStateInit(false);
  simulation.addLogger(logger);
  simulation.run();
  std::cout << "Completed 3-phase EMT " << mode << " at dt=" << step
            << " s; results in logs/linear_rl_3ph_" << mode << ".csv\n";
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 5 ||
      (std::string(argv[1]) != "native" && std::string(argv[1]) != "me" &&
       std::string(argv[1]) != "cs" && std::string(argv[1]) != "cs-rollback")) {
    std::cerr << "Usage: FMI_EMT_Ph3_LinearRLStability "
                 "native|me|cs|cs-rollback model.fmu|- "
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
