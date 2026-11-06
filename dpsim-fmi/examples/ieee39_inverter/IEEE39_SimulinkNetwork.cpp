// SPDX-License-Identifier: MPL-2.0
// C API used by the Simulink S-function. It deliberately builds the same
// dynamic IEEE39 topology and initialization as IEEE39_AvgInverter.cpp.
#define main ieee39_avg_inverter_example_main
#include "IEEE39_AvgInverter.cpp"
#undef main

#include <chrono>
#include <dpsim-models/EMT/EMT_Ph3_ControlledCurrentSource.h>
#include <memory>
#include <stdexcept>

namespace {
struct Network {
  std::vector<EMT::SimNode::Ptr> nodes;
  std::shared_ptr<EMT::Ph3::ControlledCurrentSource> inverterSource;
  std::shared_ptr<EMT::Ph3::Switch> loadSwitch;
  std::unique_ptr<Simulation> sim;
  Real dt;
};

void readVoltages(Network &network, double *out) {
  for (std::size_t bus = 0; bus < network.nodes.size(); ++bus)
    for (int phase = 0; phase < 3; ++phase)
      out[3 * bus + phase] = network.nodes[bus]->voltage()(phase, 0);
}

SystemTopology makeDynamicSystem(Network &network, Real loadStepTime) {
  SystemNodeList systemNodes;
  SystemComponentList components;
  for (const auto &bus : buses) {
    auto node =
        EMT::SimNode::make("N" + std::to_string(bus.id), PhaseType::ABC);
    network.nodes.push_back(node);
    systemNodes.push_back(node);
  }

  for (const auto &gen : generators) {
    auto machine = EMT::Ph3::SynchronGenerator4OrderVBR::make(
        "Gen_N" + std::to_string(gen.bus));
    machine->setOperationalParametersPerUnit(
        gen.baseS * 1e6, baseVoltage, frequency, gen.H, gen.Xd, gen.Xq, gen.Xl,
        gen.Xd_t, gen.Xq_t, gen.Td0_t, gen.Tq0_t);
    machine->connect({network.nodes.at(gen.bus - 1)});
    components.push_back(machine);
  }

  int loadIndex = 0;
  for (const auto &bus : buses) {
    if (bus.pd == 0.0 && bus.qd == 0.0)
      continue;
    const Real g = bus.gs * 1e6 / (baseVoltage * baseVoltage);
    const Real b = bus.bs * 1e6 / (baseVoltage * baseVoltage);
    auto load = EMT::Ph3::RXLoad::make("load" + std::to_string(++loadIndex));
    load->setParameters(parameter3(bus.pd * 1e6 / 3.0),
                        parameter3(bus.qd * 1e6 / 3.0), bus.vm * baseVoltage);
    load->connect({network.nodes.at(bus.id - 1)});
    components.push_back(load);
    if (g != 0.0 || b != 0.0) {
      auto shunt = EMT::Ph3::Shunt::make("Shunt_N" + std::to_string(bus.id));
      shunt->setParameters(g, b);
      shunt->connect({network.nodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }
  for (const auto &bus : buses) {
    const Real g = bus.gs * 1e6 / (baseVoltage * baseVoltage);
    const Real b = bus.bs * 1e6 / (baseVoltage * baseVoltage);
    if ((g != 0.0 || b != 0.0) && bus.pd == 0.0 && bus.qd == 0.0) {
      auto shunt = EMT::Ph3::Shunt::make("Shunt_N" + std::to_string(bus.id));
      shunt->setParameters(g, b);
      shunt->connect({network.nodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }
  int lineIndex = 0, transformerIndex = 0;
  for (const auto &branch : branches) {
    const Real tap = branch.tap == 0.0 ? 1.0 : branch.tap;
    const Real r = branch.r * zBase * tap * tap;
    const Real l = branch.x * zBase * tap * tap / omega;
    if (branch.tap == 0.0) {
      auto line = EMT::Ph3::PiLine::make("line" + std::to_string(++lineIndex));
      line->setParameters(parameter3(r), parameter3(l),
                          parameter3(branch.b / zBase / omega),
                          parameter3(0.0));
      line->connect(
          {network.nodes.at(branch.from - 1), network.nodes.at(branch.to - 1)});
      components.push_back(line);
    } else {
      auto trafo = EMT::Ph3::Transformer::make(
          "transformer" + std::to_string(++transformerIndex));
      trafo->setParameters(baseVoltage, baseVoltage, baseMVA * 1e6, tap, 0.0,
                           parameter3(r), parameter3(l));
      trafo->connect(
          {network.nodes.at(branch.from - 1), network.nodes.at(branch.to - 1)});
      components.push_back(trafo);
    }
  }

  network.inverterSource =
      std::dynamic_pointer_cast<EMT::Ph3::ControlledCurrentSource>(
          EMT::Ph3::ControlledCurrentSource::make("Gen30_AveragedInverter"));
  network.inverterSource->setParameters(Matrix::Zero(3, 1));
  network.inverterSource->connect({EMT::SimNode::GND, network.nodes.at(29)});
  components.push_back(network.inverterSource);

  const Real stepPower = 100e6;
  const Real bus39Peak = buses[38].vm * baseVoltage * dpEnvelopeScale;
  const Real stepResistance = 1.5 * bus39Peak * bus39Peak / stepPower;
  auto stepNode = EMT::SimNode::make("N39_LoadStep", PhaseType::ABC);
  auto stepLoad = EMT::Ph3::Resistor::make("N39_LoadStep_R");
  stepLoad->setParameters(parameter3(stepResistance));
  stepLoad->connect({EMT::SimNode::GND, stepNode});
  network.loadSwitch = EMT::Ph3::Switch::make("N39_LoadStep_Switch");
  network.loadSwitch->setParameters(parameter3(1e12), parameter3(1e-3), false);
  network.loadSwitch->connect({stepNode, network.nodes.at(38)});
  systemNodes.push_back(stepNode);
  components.push_back(stepLoad);
  components.push_back(network.loadSwitch);

  auto system = SystemTopology(frequency, systemNodes, components);
  auto systemPF = makePowerflowSystem();
  Simulation pf("IEEE39_Simulink_PF", Logger::Level::off);
  pf.setSystem(systemPF);
  pf.setDomain(Domain::SP);
  pf.setSolverType(Solver::Type::NRP);
  pf.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  pf.setTimeStep(.01);
  pf.setFinalTime(.01);
  pf.doInitFromNodesAndTerminals(false);
  pf.run();
  system.initWithPowerflow(systemPF, Domain::EMT);
  return system;
}
} // namespace

extern "C" {
void *ieee39_network_create(double dt, double finalTime, double loadStepTime,
                            double *initialVoltages,
                            double *initializationSeconds, char *error,
                            int errorCapacity) {
  try {
    if (dt <= 0.0 || finalTime <= dt || loadStepTime <= 0.0 ||
        loadStepTime >= finalTime ||
        std::abs(finalTime / dt - std::round(finalTime / dt)) > 1e-8 ||
        std::abs(loadStepTime / dt - std::round(loadStepTime / dt)) > 1e-8)
      throw std::invalid_argument("invalid simulation time grid");
    const auto begin = std::chrono::steady_clock::now();
    auto network = std::make_unique<Network>();
    network->dt = dt;
    auto system = makeDynamicSystem(*network, loadStepTime);
    network->sim = std::make_unique<Simulation>("IEEE39_Simulink_ClosedLoop",
                                                Logger::Level::off);
    network->sim->addEvent(
        SwitchEvent3Ph::make(loadStepTime, network->loadSwitch, true));
    network->sim->setSystem(system);
    network->sim->setDomain(Domain::EMT);
    network->sim->setSolverType(Solver::Type::MNA);
    network->sim->setTimeStep(dt);
    network->sim->setFinalTime(finalTime);
    network->sim->doInitFromNodesAndTerminals(true);
    network->sim->doSystemMatrixRecomputation(true);
    network->sim->initialize();
    network->sim->start();
    readVoltages(*network, initialVoltages);
    *initializationSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - begin)
            .count();
    return network.release();
  } catch (const std::exception &e) {
    if (error && errorCapacity > 0)
      std::snprintf(error, static_cast<std::size_t>(errorCapacity), "%s",
                    e.what());
    return nullptr;
  }
}

int ieee39_network_step(void *opaque, const double *current, double *voltages,
                        double *stepSeconds, char *error, int errorCapacity) {
  try {
    auto &network = *static_cast<Network *>(opaque);
    Matrix currentMatrix(3, 1);
    for (int phase = 0; phase < 3; ++phase)
      currentMatrix(phase, 0) = current[phase];
    network.inverterSource->setParameters(currentMatrix);
    const auto begin = std::chrono::steady_clock::now();
    network.sim->next();
    *stepSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - begin)
            .count();
    readVoltages(network, voltages);
    return 0;
  } catch (const std::exception &e) {
    if (error && errorCapacity > 0)
      std::snprintf(error, static_cast<std::size_t>(errorCapacity), "%s",
                    e.what());
    return 1;
  }
}

void ieee39_network_destroy(void *opaque) {
  delete static_cast<Network *>(opaque);
}
}
