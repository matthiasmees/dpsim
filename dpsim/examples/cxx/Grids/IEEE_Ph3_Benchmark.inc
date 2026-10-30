// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
// Shared three-phase DP / EMT SG4 benchmark implementation, based on MATPOWER.
#ifdef IEEE300_CASE
#include "IEEE300_Data.h"
namespace IEEEGridData = IEEE300Data;
#else
#include "IEEE118_Data.h"
namespace IEEEGridData = IEEE118Data;
#endif
#include <DPsim.h>
#include <dpsim-models/DP/DP_Ph3_SSN_ImpedanceLoad.h>
#include <dpsim-models/DP/DP_Ph3_SSN_PiLine.h>
#include <dpsim-models/DP/DP_Ph3_SSN_SynchronGenerator4Order.h>
#include <dpsim-models/DP/DP_Ph3_SSN_Transformer.h>
#include <dpsim-models/EMT/EMT_LinearSSNFused.h>
#include <dpsim-models/EMT/EMT_LinearSSNGroup.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_ImpedanceLoad.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_SynchronGenerator4Order.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_Transformer.h>
#include <dpsim/SequentialScheduler.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

using namespace CPS;
using namespace DPsim;
using namespace IEEEGridData;

namespace {
#ifdef IEEE300_CASE
constexpr char gridName[] = "IEEE300";
#else
constexpr char gridName[] = "IEEE118";
#endif
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
Complex powerflowGeneratorPower(const SystemTopology &system,
                                const Generator &gen) {
  const String name = "Gen" + std::to_string(gen.bus);
  for (const auto &object : system.mComponents) {
    if (object->name() == name) {
      auto machine =
          std::dynamic_pointer_cast<SP::Ph1::SynchronGenerator>(object);
      if (machine)
        return machine->getApparentPower();
    }
  }
  return Complex(gen.pg, gen.qg) * 1e6;
}
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
                         SystemComponentList &components,
                         bool useSSNImpedanceLoad) {
  int count = 0;
  for (const auto &bus : buses) {
    const auto &node = nodes.at(bus.id - 1);
    const Real peak = bus.vm * baseV(bus) * dpEnvelopeScale;
    const Real p = bus.pd * 1e6, q = bus.qd * 1e6;
    const Real g = bus.gs * 1e6 / std::pow(baseV(bus), 2);
    const Real b = bus.bs * 1e6 / std::pow(baseV(bus), 2);
    if (useSSNImpedanceLoad) {
      if (p != 0. || q != 0. || g != 0. || b != 0.) {
        auto load = DP::Ph3::SSN::ImpedanceLoad::make("ImpedanceLoad" +
                                                      std::to_string(bus.id));
        load->setParameters(p, q, bus.vm * baseV(bus), g, b);
        load->connect({node});
        components.push_back(load);
      }
      continue;
    }
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
void addBranches(const Nodes &nodes, SystemComponentList &components,
                 bool useSSNPiLine, bool usePiLineNoSub,
                 bool useTransformerNoSub, bool useSSNTransformer) {
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
        const auto lineName = "Line" + std::to_string(++lines);
        if (usePiLineNoSub) {
          auto line = EMT::Ph3::PiLine::make(lineName, false);
          line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                              parameter3(0.));
          line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(line);
        } else if (useSSNPiLine) {
          auto line = EMT::Ph3::SSN::PiLine::make(lineName);
          line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                              parameter3(0.));
          line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(line);
        } else {
          auto line = EMT::Ph3::PiLine::make(lineName);
          line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                              parameter3(0.));
          line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(line);
        }
      } else {
        const auto lineName = "Line" + std::to_string(++lines);
        if (useSSNPiLine) {
          auto line = DP::Ph3::SSN::PiLine::make(lineName);
          line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                              parameter3(0.));
          line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(line);
        } else {
          auto line = DP::Ph3::PiLine::make(lineName);
          line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                              parameter3(0.));
          line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(line);
        }
      }
    } else {
      const Real transformerScale = (vf > vt) ? ratio : tap;
      const Real r = branch.r * zTo * transformerScale * transformerScale;
      const Real l =
          branch.x * zTo * transformerScale * transformerScale / omega;
      const Real ratedPower =
          branch.rateA > 0. ? branch.rateA * 1e6 : baseMVA * 1e6;
      if constexpr (emt) {
        const auto name = "Trafo" + std::to_string(++transformers);
        if (useSSNTransformer) {
          auto trafo = EMT::Ph3::SSN::Transformer::make(name);
          trafo->setParameters(vf, vt, ratedPower, std::abs(ratio),
                               branch.angle * PI / 180., parameter3(r),
                               parameter3(l));
          trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(trafo);
        } else {
          auto trafo = EMT::Ph3::Transformer::make(name, !useTransformerNoSub);
          trafo->setParameters(vf, vt, ratedPower, std::abs(ratio),
                               branch.angle * PI / 180., parameter3(r),
                               parameter3(l));
          trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(trafo);
        }
      } else {
        const auto name = "Trafo" + std::to_string(++transformers);
        if (useSSNTransformer) {
          auto trafo = DP::Ph3::SSN::Transformer::make(name);
          trafo->setParameters(vf, vt, ratedPower, std::abs(ratio),
                               branch.angle * PI / 180., parameter3(r),
                               parameter3(l));
          trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(trafo);
        } else {
          auto trafo = DP::Ph3::Transformer::make(name);
          trafo->setParameters(vf, vt, ratedPower, std::abs(ratio),
                               branch.angle * PI / 180., parameter3(r),
                               parameter3(l));
          trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
          components.push_back(trafo);
        }
      }
    }
  }
}

String simulationName(const String &domain, bool useSSNPiLine,
                      bool useSSNImpedanceLoad, bool useSSNTransformer,
                      bool useSSNGenerator) {
  String name = String(gridName) + "_" + domain + "3Ph_CPP";
  if (useSSNPiLine)
    name += "_SSNLine";
  if (useSSNImpedanceLoad)
    name += "_SSNLoad";
  if (useSSNTransformer)
    name += "_SSNTrafo";
  if (useSSNGenerator)
    name += "_SSNGen4Order";
  return name;
}

void runDP(const SystemTopology &systemPF, bool useSSNPiLine,
           bool useSSNImpedanceLoad, bool useSSNTransformer,
           bool useSSNGenerator, bool profile, bool shortRun, bool noLogger) {
  String name = simulationName("DP", useSSNPiLine, useSSNImpedanceLoad,
                               useSSNTransformer, useSSNGenerator);
  if (profile)
    name += "_DetailedProfile";
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
    const auto &bus = busFor(gen.bus);
    if (useSSNGenerator) {
      auto machine = DP::Ph3::SSN::SynchronGenerator4Order::make(
          "Gen" + std::to_string(gen.bus));
      const Real machineBase = dynamicMachineBase(gen);
      machine->setParametersOperationalPerUnit(
          machineBase, baseV(bus), frequency, H, Ra, XdTransient, XqTransient,
          Xd, Xq, Td0Transient, Tq0Transient);
      const Complex power = powerflowGeneratorPower(systemPF, gen);
      machine->setPowerflowInitialValues(power.real() / machineBase,
                                         power.imag() / machineBase, bus.vm,
                                         bus.va * PI / 180.);
      machine->connect({DP::SimNode::GND, nodes.at(gen.bus - 1)});
      components.push_back(machine);
    } else {
      auto machine = DP::Ph3::SynchronGenerator4OrderVBR::make(
          "Gen" + std::to_string(gen.bus));
      machine->setParametersOperationalPerUnit(
          dynamicMachineBase(gen), baseV(bus), frequency, 2, 1., Ra, Xd, Xq,
          XdTransient, XqTransient, Xl, Td0Transient, Tq0Transient, H);
      machine->connect({nodes.at(gen.bus - 1)});
      components.push_back(machine);
    }
  }
  addLoadsAndShuntsDP(nodes, components, useSSNImpedanceLoad);
  addBranches(nodes, components, useSSNPiLine, false, false,
              useSSNTransformer);
  SystemTopology system(frequency, systemNodes, components);
  system.initWithPowerflow(systemPF, Domain::DP);
  auto logger = DataLogger::make(name);
  for (int id : {1, 10, 30, 69, 100, 118})
    logger->logAttribute("N" + std::to_string(id) + ".V",
                         nodes.at(id - 1)->attribute("v"));
  Simulation sim(name, noLogger ? Logger::Level::off : Logger::Level::info);
  sim.setSystem(system);
  sim.setDomain(Domain::DP);
  sim.setSolverType(Solver::Type::MNA);
  sim.setTimeStep(1e-3);
  sim.setFinalTime(shortRun ? 1. : 30.);
  sim.doInitFromNodesAndTerminals(true);
  sim.setLogStepTimes(profile);
  if (profile) {
    sim.setScheduler(std::make_shared<SequentialScheduler>(
        "logs/" + name + "/" + name + "_TaskTimes.csv", Logger::Level::off));
  } else if (!noLogger) {
    sim.addLogger(logger);
  }
  sim.run();
  if (profile) {
    sim.logStepTimes(name + "_StepTimes");
    sim.logLUTimes();
  }
}

void runEMT(const SystemTopology &systemPF, bool useSSNPiLine,
            bool usePiLineNoSub, bool useRXLoadNoSub,
            bool useTransformerNoSub, bool useSSNImpedanceLoad,
            bool useSSNTransformer, bool useSSNGenerator, bool profile,
            bool shortRun, bool noLogger, bool useKluPartial, bool use50usStep,
            bool useLinearAggregate, bool usePaperSSNGroup, bool use10sRun,
            bool useSparseLU, bool useExactDenseLU,
            bool usePrecomputedGroupMaps, bool useFullyReducedGroupMaps) {
  String name = simulationName("EMT", useSSNPiLine, useSSNImpedanceLoad,
                               useSSNTransformer, useSSNGenerator);
  if (usePiLineNoSub || useRXLoadNoSub || useTransformerNoSub) {
    name = String(gridName) + "_EMT3Ph_CPP_NoSub";
    if (useSSNImpedanceLoad)
      name += "_SSNLoad";
    if (useSSNTransformer)
      name += "_SSNTrafo";
    if (useSSNGenerator)
      name += "_SSNGen4Order";
  }
  if (profile)
    name += "_DetailedProfile";
  if (useLinearAggregate)
    name += "_LinearAggregate";
  if (usePaperSSNGroup)
    name += "_SSNGroup";
  if (usePaperSSNGroup && !usePrecomputedGroupMaps)
    name += "_Baseline";
  else if (usePaperSSNGroup && !useFullyReducedGroupMaps)
    name += "_SharedMaps";
  else if (usePaperSSNGroup)
    name += "_FullyReduced";
  Logger::setLogDir("logs/" + name);
  auto logger = DataLogger::make(name, true, 2);
  std::vector<EMT::SimNode::Ptr> nodes;
  SystemNodeList systemNodes;
  for (const auto &bus : buses) {
    auto node =
        EMT::SimNode::make("N" + std::to_string(bus.id), PhaseType::ABC);
    nodes.push_back(node);
  }
  std::unordered_map<int, EMT::SimNode::Ptr> retainedByBus;
  if (usePaperSSNGroup) {
    for (const auto &gen : generators)
      retainedByBus.emplace(
          gen.bus,
          EMT::SimNode::make("N" + std::to_string(gen.bus), PhaseType::ABC));
    for (const auto &entry : retainedByBus)
      systemNodes.push_back(entry.second);
    for (const auto &bus : buses) {
      const auto pfNode =
          const_cast<SystemTopology &>(systemPF).node<SP::SimNode>(
              "N" + std::to_string(bus.id));
      if (!pfNode)
        throw std::runtime_error("Missing power-flow node for paper SSN group");
      // Use the exact complex PF voltage copied by SystemTopology for EMT
      // initialization, including its RMS-to-peak scaling and angle.
      nodes.at(bus.id - 1)->setInitialVoltage(pfNode->singleVoltage());
    }
  } else {
    for (const auto &node : nodes)
      systemNodes.push_back(node);
  }
  auto simulationNode = [&](int bus) -> EMT::SimNode::Ptr {
    if (usePaperSSNGroup)
      return retainedByBus.at(bus);
    return nodes.at(bus - 1);
  };
  SystemComponentList components;
  for (const auto &gen : generators) {
    const auto &bus = busFor(gen.bus);
    if (useSSNGenerator) {
      auto machine = EMT::Ph3::SSN::SynchronGenerator4Order::make(
          "Gen" + std::to_string(gen.bus));
      const Real machineBase = dynamicMachineBase(gen);
      machine->setParametersOperationalPerUnit(
          machineBase, baseV(bus), frequency, H, Ra, XdTransient, XqTransient,
          Xd, Xq, Td0Transient, Tq0Transient);
      const Complex power = powerflowGeneratorPower(systemPF, gen);
      machine->setPowerflowInitialValues(power.real() / machineBase,
                                         power.imag() / machineBase, bus.vm,
                                         bus.va * PI / 180.);
      if (gen.bus == 10) {
        logger->logAttribute("Gen10.SSNState", machine->attribute("x"));
        logger->logAttribute("Gen10.SSNCurrent", machine->attribute("i_intf"));
      }
      machine->connect({EMT::SimNode::GND, simulationNode(gen.bus)});
      components.push_back(machine);
    } else {
      auto machine = EMT::Ph3::SynchronGenerator4OrderVBR::make(
          "Gen" + std::to_string(gen.bus));
      machine->setOperationalParametersPerUnit(
          dynamicMachineBase(gen), baseV(bus), frequency, H, Xd, Xq, Xl,
          XdTransient, XqTransient, Td0Transient, Tq0Transient);
      if (gen.bus == 10 || gen.bus == 69 || gen.bus == 89) {
        const auto prefix = "Gen" + std::to_string(gen.bus) + ".";
        logger->logAttribute(prefix + "w_r", machine->attribute("w_r"));
        logger->logAttribute(prefix + "Tm", machine->attribute("Tm"));
        logger->logAttribute(prefix + "Te", machine->attribute("Te"));
      }
      machine->connect({simulationNode(gen.bus)});
      components.push_back(machine);
    }
  }
  int loadIndex = 0;
  SystemComponentList linearComponents;
  auto &linearTarget =
      (useLinearAggregate || usePaperSSNGroup) ? linearComponents : components;
  for (const auto &bus : buses) {
    const Real g = bus.gs * 1e6 / std::pow(baseV(bus), 2);
    const Real b = bus.bs * 1e6 / std::pow(baseV(bus), 2);
    if (useSSNImpedanceLoad) {
      if (bus.pd != 0. || bus.qd != 0. || g != 0. || b != 0.) {
        auto load = EMT::Ph3::SSN::ImpedanceLoad::make("ImpedanceLoad" +
                                                       std::to_string(bus.id));
        load->setParameters(bus.pd * 1e6, bus.qd * 1e6, bus.vm * baseV(bus), g,
                            b);
        load->connect({nodes.at(bus.id - 1)});
        linearTarget.push_back(load);
      }
      continue;
    }
    if (bus.pd != 0. || bus.qd != 0.) {
      auto load = EMT::Ph3::RXLoad::make("Load" + std::to_string(++loadIndex),
                                         !useRXLoadNoSub);
      load->setParameters(parameter3(bus.pd * 1e6 / 3.),
                          parameter3(bus.qd * 1e6 / 3.), bus.vm * baseV(bus));
      load->connect({nodes.at(bus.id - 1)});
      linearTarget.push_back(load);
    }
    if (g != 0. || b != 0.) {
      auto shunt = EMT::Ph3::Shunt::make("Shunt" + std::to_string(bus.id));
      shunt->setParameters(g, b);
      shunt->connect({nodes.at(bus.id - 1)});
      linearTarget.push_back(shunt);
    }
  }
  addBranches(nodes, linearTarget, useSSNPiLine, usePiLineNoSub,
              useTransformerNoSub, useSSNTransformer);
  if (useLinearAggregate) {
    EMT::LinearSSNFused::ComponentList aggregateComponents;
    aggregateComponents.reserve(linearComponents.size());
    for (const auto &component : linearComponents)
      aggregateComponents.push_back(
          std::dynamic_pointer_cast<EMT::LinearSSNFused::Component>(
              component));
    components.push_back(EMT::LinearSSNFused::make(
        "LinearNetwork", std::move(aggregateComponents), nodes));
  } else if (usePaperSSNGroup) {
    EMT::LinearSSNGroup::ComponentList groupComponents;
    groupComponents.reserve(linearComponents.size());
    for (const auto &component : linearComponents)
      groupComponents.push_back(
          std::dynamic_pointer_cast<EMT::LinearSSNGroup::Component>(component));
    EMT::SimNode::List retainedNodes;
    retainedNodes.reserve(retainedByBus.size());
    for (const auto &gen : generators)
      retainedNodes.push_back(retainedByBus.at(gen.bus));
    auto group = EMT::LinearSSNGroup::make(
        "LinearNetworkGroup", std::move(groupComponents), nodes,
        std::move(retainedNodes));
    group->setPrecomputedInternalMaps(usePrecomputedGroupMaps);
    group->setFullyReducedInternalMaps(useFullyReducedGroupMaps);
    components.push_back(std::move(group));
  }
  SystemTopology system(frequency, systemNodes, components);
  system.initWithPowerflow(systemPF, Domain::EMT);
  for (int id : {1, 10, 30, 69, 100, 118}) {
    if (!usePaperSSNGroup || retainedByBus.count(id) != 0)
      logger->logAttribute("N" + std::to_string(id) + ".V",
                           usePaperSSNGroup
                               ? retainedByBus.at(id)->attribute("v")
                               : nodes.at(id - 1)->attribute("v"));
  }
  Simulation sim(name, noLogger ? Logger::Level::off : Logger::Level::info);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);
  if (usePaperSSNGroup && useSparseLU) {
    sim.setDirectLinearSolverImplementation(DirectLinearSolverImpl::SparseLU);
  } else if (usePaperSSNGroup && !useKluPartial) {
    // The paper reduction produces a dense Schur-complement block on the
    // retained generator buses. DenseLU is the default for this model; pass
    // klupartial to benchmark the sparse KLU alternative on the same system.
    sim.setDirectLinearSolverImplementation(DirectLinearSolverImpl::DenseLU);
    if (!useExactDenseLU) {
      DirectLinearSolverConfiguration denseConfig;
      denseConfig.setStaticPreconditioner(true);
      sim.setDirectLinearSolverConfiguration(denseConfig);
    }
  } else if (useKluPartial) {
    DirectLinearSolverConfiguration kluConfig;
    kluConfig.setPartialRefactorizationMethod(
        PARTIAL_REFACTORIZATION_METHOD::FACTORIZATION_PATH);
    kluConfig.setFillInReductionMethod(FILL_IN_REDUCTION_METHOD::AMD_NV);
    sim.setDirectLinearSolverConfiguration(kluConfig);
  } else if (useLinearAggregate) {
    // Keep the measured fastest default explicit: standard KLU numeric
    // refactorization, with the optional partial-refactorization path off.
    DirectLinearSolverConfiguration standardKluConfig;
    standardKluConfig.setPartialRefactorizationMethod(
        PARTIAL_REFACTORIZATION_METHOD::NO_PARTIAL_REFACTORIZATION);
    sim.setDirectLinearSolverConfiguration(standardKluConfig);
  }
  sim.setTimeStep(use50usStep ? 5e-5 : 5e-4);
  sim.setFinalTime(use10sRun ? 10. : (shortRun ? 1. : 30.));
  sim.setLogStepTimes(profile);
  sim.doInitFromNodesAndTerminals(true);
  if (profile) {
    sim.setScheduler(std::make_shared<SequentialScheduler>(
        "logs/" + name + "/" + name + "_TaskTimes.csv", Logger::Level::off));
  } else if (!noLogger) {
    sim.addLogger(logger);
  }
  sim.run();
  if (profile) {
    sim.logStepTimes(name + "_StepTimes");
    sim.logLUTimes();
  }
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 2 || argc > 12 ||
      (std::string(argv[1]) != "dp" && std::string(argv[1]) != "emt")) {
    std::cerr << "Usage: IEEE118_Ph3/IEEE300_Ph3 dp|emt "
                 "[allssn|ssn|nosub] [ssnload] "
                 "[ssntrafo] [ssngen] [conventional] [short] [nolog] "
                 "[profile] [dt50us] [aggregate] [paper] [sparselu] "
                 "[denseexact] [groupfull|groupshared|groupbaseline]\n";
    return 2;
  }
#ifdef IEEE300_CASE
  // The measured fastest IEEE300 EMT formulation is the fused linear SSN
  // aggregate. Generators remain conventional VBR components.
  bool useSSNPiLine = true;
  bool useSSNImpedanceLoad = true;
  bool useSSNTransformer = true;
#elif defined(IEEE118_DEFAULT_SSN_MODE)
  bool useSSNPiLine = true;
  bool useSSNImpedanceLoad = true;
  bool useSSNTransformer = true;
#else
  bool useSSNPiLine = false;
  bool useSSNImpedanceLoad = false;
  bool useSSNTransformer = false;
#endif
  bool useSSNGenerator = false;
  bool usePiLineNoSub = false;
  bool useRXLoadNoSub = false;
  bool useTransformerNoSub = false;
  bool profile = false;
  bool shortRun = false;
  bool noLogger = false;
  bool useKluPartial = false;
  bool use50usStep = false;
  bool forceConventional = false;
  bool useLinearAggregate =
#ifdef IEEE300_CASE
      true;
#else
      false;
#endif
  bool usePaperSSNGroup = false;
  bool use10sRun = false;
  bool useSparseLU = false;
  bool useExactDenseLU = false;
  bool usePrecomputedGroupMaps = true;
  bool useFullyReducedGroupMaps = false;
  for (int arg = 2; arg < argc; ++arg) {
    const std::string mode(argv[arg]);
    if (mode == "allssn") {
      useSSNPiLine = true;
      useSSNImpedanceLoad = true;
      useSSNTransformer = true;
    } else if (mode == "ssn")
      useSSNPiLine = true;
    else if (mode == "nosub")
      usePiLineNoSub = useRXLoadNoSub = useTransformerNoSub = true;
    else if (mode == "nosubline")
      usePiLineNoSub = true;
    else if (mode == "nosubload")
      useRXLoadNoSub = true;
    else if (mode == "nosubtrafo")
      useTransformerNoSub = true;
    else if (mode == "ssnload")
      useSSNImpedanceLoad = true;
    else if (mode == "ssntrafo")
      useSSNTransformer = true;
    else if (mode == "ssngen")
      useSSNGenerator = true;
    else if (mode == "conventional")
      forceConventional = true;
    else if (mode == "profile") {
      profile = true;
      noLogger = true;
    } else if (mode == "short")
      shortRun = true;
    else if (mode == "10s")
      use10sRun = true;
    else if (mode == "nolog")
      noLogger = true;
    else if (mode == "klupartial")
      useKluPartial = true;
    else if (mode == "dt50us")
      use50usStep = true;
    else if (mode == "aggregate")
      useLinearAggregate = true;
    else if (mode == "paper")
      usePaperSSNGroup = true;
    else if (mode == "sparselu")
      useSparseLU = true;
    else if (mode == "denseexact")
      useExactDenseLU = true;
    else if (mode == "groupbaseline")
      usePrecomputedGroupMaps = useFullyReducedGroupMaps = false;
    else if (mode == "groupshared") {
      usePrecomputedGroupMaps = true;
      useFullyReducedGroupMaps = false;
    } else if (mode == "groupfull") {
      usePrecomputedGroupMaps = true;
      useFullyReducedGroupMaps = true;
    }
    else {
      std::cerr << "Optional modes are 'allssn', 'ssn', 'nosub', 'ssnload', "
                   "'ssntrafo', 'ssngen', 'short', 'nolog', 'klupartial', "
                   "'profile', 'dt50us', '10s', 'aggregate', 'paper', "
                   "'nosubline', 'nosubload', 'nosubtrafo', and "
                   "'conventional', 'sparselu', 'denseexact', and "
                   "'groupfull', 'groupshared', and 'groupbaseline'.\n";
      return 2;
    }
  }
  if (forceConventional) {
    useSSNPiLine = false;
    useSSNImpedanceLoad = false;
    useSSNTransformer = false;
    useSSNGenerator = false;
    useLinearAggregate = false;
    usePaperSSNGroup = false;
  }
  if (useLinearAggregate && (!useSSNPiLine || !useSSNImpedanceLoad ||
                             !useSSNTransformer || usePiLineNoSub ||
                             useRXLoadNoSub || useTransformerNoSub)) {
    std::cerr << "aggregate requires allssn and cannot be combined with "
                 "nosub or conventional.\n";
    return 2;
  }
  if (usePaperSSNGroup &&
      (!useSSNPiLine || !useSSNImpedanceLoad || !useSSNTransformer ||
       usePiLineNoSub || useRXLoadNoSub || useTransformerNoSub ||
       useLinearAggregate)) {
    std::cerr << "paper requires allssn and cannot be combined with nosub or "
                 "aggregate.\n";
    return 2;
  }
  auto systemPF = makePowerflowSystem();
  Simulation pf(String(gridName) + "_SP_PF", Logger::Level::off);
  pf.setSystem(systemPF);
  pf.setTimeStep(.01);
  pf.setFinalTime(.01);
  pf.setDomain(Domain::SP);
  pf.setSolverType(Solver::Type::NRP);
  pf.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  pf.doInitFromNodesAndTerminals(false);
  pf.run();
  if (std::string(argv[1]) == "dp")
    runDP(systemPF, useSSNPiLine, useSSNImpedanceLoad, useSSNTransformer,
          useSSNGenerator, profile, shortRun, noLogger);
  else
    runEMT(systemPF, useSSNPiLine, usePiLineNoSub, useRXLoadNoSub,
           useTransformerNoSub, useSSNImpedanceLoad, useSSNTransformer,
           useSSNGenerator, profile, shortRun, noLogger,
           useKluPartial, use50usStep, useLinearAggregate, usePaperSSNGroup,
           use10sRun, useSparseLU, useExactDenseLU,
           usePrecomputedGroupMaps, useFullyReducedGroupMaps);
  return 0;
}
