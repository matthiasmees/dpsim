// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
// IEEE 118-bus three-phase DP / EMT SG4 example, based on MATPOWER case118.
#include "IEEE118_Data.h"
#include <DPsim.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace CPS;
using namespace DPsim;
using namespace IEEE118Data;

namespace {
constexpr Real frequency = 60.;
constexpr Real omega = 2. * PI * frequency;
constexpr Real baseMVA = 100.;
constexpr Real H = 4.;
constexpr Real Ra = .002;
constexpr Real Xl = .15;
constexpr Real Xd = 1.75;
constexpr Real XdTransient = .425;
constexpr Real Xq = 1.75;
constexpr Real XqTransient = .425;
constexpr Real Td0Transient = 5.;
constexpr Real Tq0Transient = .5;
constexpr Real dpEnvelopeScale = 0.8164965809277260;

Real baseV(const Bus &bus) { return bus.baseKV * 1e3; }
Real zBase(const Bus &bus) { return baseV(bus) * baseV(bus) / (baseMVA * 1e6); }
Real dynamicMachineBase(const Generator &gen) {
  return std::max(gen.mbase, gen.pmax) * 1e6;
}
Matrix parameter3(Real x) { return Math::singlePhaseParameterToThreePhase(x); }

const Bus &busFor(int id) { return buses.at(static_cast<size_t>(id - 1)); }
PowerflowBusType pfBusType(int type) {
  if (type == 3)
    return PowerflowBusType::VD;
  if (type == 2)
    return PowerflowBusType::PV;
  return PowerflowBusType::PQ;
}

SystemTopology makePowerflowSystem() {
  SystemNodeList nodes;
  SystemComponentList components;
  std::vector<SP::SimNode::Ptr> busNodes;
  busNodes.reserve(buses.size());
  for (const auto &bus : buses) {
    auto node =
        SP::SimNode::make("N" + std::to_string(bus.id), PhaseType::Single);
    busNodes.push_back(node);
    nodes.push_back(node);
  }

  for (const auto &gen : generators) {
    const auto &bus = busFor(gen.bus);
    auto machine =
        SP::Ph1::SynchronGenerator::make("Gen" + std::to_string(gen.bus));
    machine->setParameters(gen.mbase * 1e6, baseV(bus), gen.pg * 1e6,
                           gen.vg * baseV(bus), pfBusType(bus.type),
                           gen.qg * 1e6, gen.qmax * 1e6, gen.qmin * 1e6);
    machine->setBaseVoltage(baseV(bus));
    machine->connect({busNodes.at(gen.bus - 1)});
    components.push_back(machine);
  }

  int loadIndex = 0;
  for (const auto &bus : buses) {
    if (bus.pd != 0. || bus.qd != 0.) {
      auto load = SP::Ph1::Load::make("PFLoad" + std::to_string(++loadIndex));
      load->setParameters(bus.pd * 1e6, bus.qd * 1e6, baseV(bus));
      load->connect({busNodes.at(bus.id - 1)});
      components.push_back(load);
    }
    if (bus.gs != 0. || bus.bs != 0.) {
      auto shunt = SP::Ph1::Shunt::make("PFShunt" + std::to_string(bus.id));
      shunt->setParameters(bus.gs * 1e6 / (baseV(bus) * baseV(bus)),
                           bus.bs * 1e6 / (baseV(bus) * baseV(bus)));
      shunt->setBaseVoltage(baseV(bus));
      shunt->connect({busNodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }

  int lineIndex = 0;
  int transformerIndex = 0;
  for (const auto &branch : branches) {
    const auto &from = busFor(branch.from);
    const auto &to = busFor(branch.to);
    const Real vf = baseV(from), vt = baseV(to);
    const Real tap = branch.ratio == 0. ? 1. : branch.ratio;
    const Real ratio = tap * vf / vt;
    const Real zTo = zBase(to);
    if (branch.ratio == 0. && vf == vt) {
      auto line = SP::Ph1::PiLine::make("PFLine" + std::to_string(++lineIndex));
      line->setParameters(branch.r * zBase(to), branch.x * zBase(to) / omega,
                          branch.b / zBase(to) / omega, 0.);
      line->setBaseVoltage(vt);
      line->connect({busNodes.at(branch.from - 1), busNodes.at(branch.to - 1)});
      components.push_back(line);
    } else {
      auto transformer = SP::Ph1::Transformer::make(
          "PFTrafo" + std::to_string(++transformerIndex));
      const Real highV = std::max(vf, vt);
      const Real transformerBaseZ = zTo;
      const Real transformerScale = (vf > vt) ? ratio : tap;
      const Real transformerR =
          branch.r * transformerBaseZ * transformerScale * transformerScale;
      const Real transformerL = branch.x * transformerBaseZ * transformerScale *
                                transformerScale / omega;
      const Real rate = branch.rateA > 0. ? branch.rateA * 1e6 : baseMVA * 1e6;
      transformer->setParameters(vf, vt, rate, std::abs(ratio),
                                 branch.angle * PI / 180., transformerR,
                                 transformerL);
      transformer->setBaseVoltage(highV);
      transformer->connect(
          {busNodes.at(branch.from - 1), busNodes.at(branch.to - 1)});
      components.push_back(transformer);
    }
  }
  return SystemTopology(frequency, nodes, components);
}

template <typename NodePtr>
void addLoadsAndShuntsDP(const std::vector<NodePtr> &nodes,
                         SystemComponentList &components) {
  int count = 0;
  for (const auto &bus : buses) {
    const auto &node = nodes.at(bus.id - 1);
    const Real peak = bus.vm * baseV(bus) * dpEnvelopeScale;
    const Real p = bus.pd * 1e6, q = bus.qd * 1e6;
    if (p > 0.) {
      auto r = DP::Ph3::Resistor::make("Load" + std::to_string(++count) + "_R");
      r->setParameters(parameter3(1.5 * peak * peak / p));
      r->connect({node, DP::SimNode::GND});
      components.push_back(r);
    }
    if (q > 0.) {
      auto l = DP::Ph3::Inductor::make("Load" + std::to_string(count) + "_L");
      l->setParameters(parameter3(1.5 * peak * peak / (omega * q)));
      l->connect({node, DP::SimNode::GND});
      components.push_back(l);
    } else if (q < 0.) {
      auto c = DP::Ph3::Capacitor::make("Load" + std::to_string(count) + "_C");
      c->setParameters(parameter3(-q / (1.5 * peak * peak * omega)));
      c->connect({node, DP::SimNode::GND});
      components.push_back(c);
    }
    const Real g = bus.gs * 1e6 / std::pow(baseV(bus), 2);
    const Real b = bus.bs * 1e6 / std::pow(baseV(bus), 2);
    if (g > 0.) {
      auto r = DP::Ph3::Resistor::make("Shunt" + std::to_string(bus.id) + "_R");
      r->setParameters(parameter3(1. / g));
      r->connect({node, DP::SimNode::GND});
      components.push_back(r);
    }
    if (b > 0.) {
      auto c =
          DP::Ph3::Capacitor::make("Shunt" + std::to_string(bus.id) + "_C");
      c->setParameters(parameter3(b / omega));
      c->connect({node, DP::SimNode::GND});
      components.push_back(c);
    } else if (b < 0.) {
      auto l = DP::Ph3::Inductor::make("Shunt" + std::to_string(bus.id) + "_L");
      l->setParameters(parameter3(-1. / (omega * b)));
      l->connect({node, DP::SimNode::GND});
      components.push_back(l);
    }
  }
}

template <typename Nodes>
void addBranches(const Nodes &nodes, SystemComponentList &components) {
  constexpr bool emt =
      std::is_same_v<typename Nodes::value_type, EMT::SimNode::Ptr>;
  int lines = 0, transformers = 0;
  for (const auto &branch : branches) {
    const auto &from = busFor(branch.from);
    const auto &to = busFor(branch.to);
    const Real vf = baseV(from), vt = baseV(to),
               tap = branch.ratio == 0. ? 1. : branch.ratio;
    const Real ratio = tap * vf / vt;
    const Real zTo = zBase(to);
    if (branch.ratio == 0. && vf == vt) {
      const Real r = branch.r * zTo, l = branch.x * zTo / omega;
      const Real c = branch.b / zTo / omega;
      if constexpr (emt) {
        auto line = EMT::Ph3::PiLine::make("Line" + std::to_string(++lines));
        line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                            parameter3(0.));
        line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(line);
      } else {
        auto line = DP::Ph3::PiLine::make("Line" + std::to_string(++lines));
        line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                            parameter3(0.));
        line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(line);
      }
    } else {
      const Real transformerScale = (vf > vt) ? ratio : tap;
      const Real r = branch.r * zTo * transformerScale * transformerScale;
      const Real l =
          branch.x * zTo * transformerScale * transformerScale / omega;
      const Real ratedPower =
          branch.rateA > 0. ? branch.rateA * 1e6 : baseMVA * 1e6;
      if constexpr (emt) {
        auto trafo = EMT::Ph3::Transformer::make(
            "Trafo" + std::to_string(++transformers));
        trafo->setParameters(vf, vt, ratedPower, std::abs(ratio),
                             branch.angle * PI / 180., parameter3(r),
                             parameter3(l));
        trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(trafo);
      } else {
        auto trafo = DP::Ph3::Transformer::make("Trafo" +
                                                std::to_string(++transformers));
        trafo->setParameters(vf, vt, ratedPower, std::abs(ratio),
                             branch.angle * PI / 180., parameter3(r),
                             parameter3(l));
        trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
        components.push_back(trafo);
      }
    }
  }
}

void runDP(const SystemTopology &systemPF) {
  const String name = "IEEE118_DP3Ph_CPP";
  Logger::setLogDir("logs/" + name);
  std::vector<DP::SimNode::Ptr> nodes;
  SystemNodeList systemNodes;
  for (const auto &bus : buses) {
    auto node = DP::SimNode::make("N" + std::to_string(bus.id), PhaseType::ABC);
    nodes.push_back(node);
    systemNodes.push_back(node);
  }
  SystemComponentList components;
  for (const auto &gen : generators) {
    auto machine = DP::Ph3::SynchronGenerator4OrderVBR::make(
        "Gen" + std::to_string(gen.bus));
    const auto &bus = busFor(gen.bus);
    machine->setParametersOperationalPerUnit(
        dynamicMachineBase(gen), baseV(bus), frequency, 2, 1., Ra, Xd, Xq,
        XdTransient, XqTransient, Xl, Td0Transient, Tq0Transient, H);
    machine->connect({nodes.at(gen.bus - 1)});
    components.push_back(machine);
  }
  addLoadsAndShuntsDP(nodes, components);
  addBranches(nodes, components);
  SystemTopology system(frequency, systemNodes, components);
  system.initWithPowerflow(systemPF, Domain::DP);
  auto logger = DataLogger::make(name);
  for (int id : {1, 10, 30, 69, 100, 118})
    logger->logAttribute("N" + std::to_string(id) + ".V",
                         nodes.at(id - 1)->attribute("v"));
  Simulation sim(name);
  sim.setSystem(system);
  sim.setDomain(Domain::DP);
  sim.setSolverType(Solver::Type::MNA);
  sim.setTimeStep(1e-3);
  sim.setFinalTime(30.);
  sim.doInitFromNodesAndTerminals(true);
  sim.addLogger(logger);
  sim.run();
}

void runEMT(const SystemTopology &systemPF) {
  const String name = "IEEE118_EMT3Ph_CPP";
  Logger::setLogDir("logs/" + name);
  auto logger = DataLogger::make(name, true, 2);
  std::vector<EMT::SimNode::Ptr> nodes;
  SystemNodeList systemNodes;
  for (const auto &bus : buses) {
    auto node =
        EMT::SimNode::make("N" + std::to_string(bus.id), PhaseType::ABC);
    nodes.push_back(node);
    systemNodes.push_back(node);
  }
  SystemComponentList components;
  for (const auto &gen : generators) {
    auto machine = EMT::Ph3::SynchronGenerator4OrderVBR::make(
        "Gen" + std::to_string(gen.bus));
    const auto &bus = busFor(gen.bus);
    machine->setOperationalParametersPerUnit(
        dynamicMachineBase(gen), baseV(bus), frequency, H, Xd, Xq, Xl,
        XdTransient, XqTransient, Td0Transient, Tq0Transient);
    if (gen.bus == 10 || gen.bus == 69 || gen.bus == 89) {
      const auto prefix = "Gen" + std::to_string(gen.bus) + ".";
      logger->logAttribute(prefix + "w_r", machine->attribute("w_r"));
      logger->logAttribute(prefix + "Tm", machine->attribute("Tm"));
      logger->logAttribute(prefix + "Te", machine->attribute("Te"));
    }
    machine->connect({nodes.at(gen.bus - 1)});
    components.push_back(machine);
  }
  int loadIndex = 0;
  for (const auto &bus : buses) {
    if (bus.pd != 0. || bus.qd != 0.) {
      auto load = EMT::Ph3::RXLoad::make("Load" + std::to_string(++loadIndex));
      load->setParameters(parameter3(bus.pd * 1e6 / 3.),
                          parameter3(bus.qd * 1e6 / 3.), bus.vm * baseV(bus));
      load->connect({nodes.at(bus.id - 1)});
      components.push_back(load);
    }
    const Real g = bus.gs * 1e6 / std::pow(baseV(bus), 2);
    const Real b = bus.bs * 1e6 / std::pow(baseV(bus), 2);
    if (g != 0. || b != 0.) {
      auto shunt = EMT::Ph3::Shunt::make("Shunt" + std::to_string(bus.id));
      shunt->setParameters(g, b);
      shunt->connect({nodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }
  addBranches(nodes, components);
  SystemTopology system(frequency, systemNodes, components);
  system.initWithPowerflow(systemPF, Domain::EMT);
  for (int id : {1, 10, 30, 69, 100, 118})
    logger->logAttribute("N" + std::to_string(id) + ".V",
                         nodes.at(id - 1)->attribute("v"));
  Simulation sim(name);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);
  sim.setTimeStep(5e-4);
  sim.setFinalTime(30.);
  sim.doInitFromNodesAndTerminals(true);
  sim.addLogger(logger);
  sim.run();
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 2 ||
      (std::string(argv[1]) != "dp" && std::string(argv[1]) != "emt")) {
    std::cerr << "Usage: IEEE118_Ph3 dp|emt\n";
    return 2;
  }
  auto systemPF = makePowerflowSystem();
  Simulation pf("IEEE118_SP_PF", Logger::Level::off);
  pf.setSystem(systemPF);
  pf.setTimeStep(.01);
  pf.setFinalTime(.01);
  pf.setDomain(Domain::SP);
  pf.setSolverType(Solver::Type::NRP);
  pf.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  pf.doInitFromNodesAndTerminals(false);
  pf.run();
  if (std::string(argv[1]) == "dp")
    runDP(systemPF);
  else
    runEMT(systemPF);
  return 0;
}
