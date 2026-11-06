// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
// Hardcoded IEEE39 counterpart to examples/Notebooks/Grids/IEEE39_DP3Ph_EMT3Ph_MATPOWER.ipynb.
// Run with `IEEE39_Ph3 dp` or `IEEE39_Ph3 emt`.
#include <DPsim.h>
#include <dpsim-fmi/CoSimulationSSN.h>
#include <dpsim-fmi/NonlinearSSN.h>

#include <cmath>
#include <dpsim/SequentialScheduler.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace CPS;
using namespace DPsim;

namespace {
constexpr Real frequency = 60.0;
constexpr Real baseMVA = 100.0;
constexpr Real baseKV = 345.0;
constexpr Real baseVoltage = baseKV * 1e3;
constexpr Real omega = 2.0 * PI * frequency;
constexpr Real zBase = baseVoltage * baseVoltage / (baseMVA * 1e6);
constexpr Real dpEnvelopeScale = 0.8164965809277260;

struct Bus {
  int id;
  int type;
  Real pd;
  Real qd;
  Real gs;
  Real bs;
  Real vm;
  Real va;
};
struct Branch {
  int from;
  int to;
  Real r;
  Real x;
  Real b;
  Real tap;
};
struct Generator {
  int bus;
  Real pg, qg;
  Real vg;
  Real baseS, H, Ra, Xl, Td0_t, Td0_s, Tq0_t, Tq0_s;
  Real Xd, Xd_t, Xd_s, Xq, Xq_t, Xq_s;
};
const Bus buses[] = {
    {1, 1, 97.6, 44.2, 0, 0, 1.0393836, -13.536602},
    {2, 1, 0, 0, 0, 0, 1.0484941, -9.7852666},
    {3, 1, 322, 2.4, 0, 0, 1.0307077, -12.276384},
    {4, 1, 500, 184, 0, 0, 1.00446, -12.626734},
    {5, 1, 0, 0, 0, 0, 1.0060063, -11.192339},
    {6, 1, 0, 0, 0, 0, 1.0082256, -10.40833},
    {7, 1, 233.8, 84, 0, 0, 0.99839728, -12.755626},
    {8, 1, 522, 176.6, 0, 0, 0.99787232, -13.335844},
    {9, 1, 6.5, -66.6, 0, 0, 1.038332, -14.178442},
    {10, 1, 0, 0, 0, 0, 1.0178431, -8.170875},
    {11, 1, 0, 0, 0, 0, 1.0133858, -8.9369663},
    {12, 1, 8.53, 88, 0, 0, 1.000815, -8.9988236},
    {13, 1, 0, 0, 0, 0, 1.014923, -8.9299272},
    {14, 1, 0, 0, 0, 0, 1.012319, -10.715295},
    {15, 1, 320, 153, 0, 0, 1.0161854, -11.345399},
    {16, 1, 329, 32.3, 0, 0, 1.0325203, -10.033348},
    {17, 1, 0, 0, 0, 0, 1.0342365, -11.116436},
    {18, 1, 158, 30, 0, 0, 1.0315726, -11.986168},
    {19, 1, 0, 0, 0, 0, 1.0501068, -5.4100729},
    {20, 1, 680, 103, 0, 0, 0.99101054, -6.8211783},
    {21, 1, 274, 115, 0, 0, 1.0323192, -7.6287461},
    {22, 1, 0, 0, 0, 0, 1.0501427, -3.1831199},
    {23, 1, 247.5, 84.6, 0, 0, 1.0451451, -3.3812763},
    {24, 1, 308.6, -92.2, 0, 0, 1.038001, -9.9137585},
    {25, 1, 224, 47.2, 0, 0, 1.0576827, -8.3692354},
    {26, 1, 139, 17, 0, 0, 1.0525613, -9.4387696},
    {27, 1, 281, 75.5, 0, 0, 1.0383449, -11.362152},
    {28, 1, 206, 27.6, 0, 0, 1.0503737, -5.9283592},
    {29, 1, 283.5, 26.9, 0, 0, 1.0501149, -3.1698741},
    {30, 2, 0, 0, 0, 0, 1.0499, -7.3704746},
    {31, 3, 9.2, 4.6, 0, 0, 0.982, 0},
    {32, 2, 0, 0, 0, 0, 0.9841, -0.1884374},
    {33, 2, 0, 0, 0, 0, 0.9972, -0.19317445},
    {34, 2, 0, 0, 0, 0, 1.0123, -1.631119},
    {35, 2, 0, 0, 0, 0, 1.0494, 1.7765069},
    {36, 2, 0, 0, 0, 0, 1.0636, 4.4684374},
    {37, 2, 0, 0, 0, 0, 1.0275, -1.5828988},
    {38, 2, 0, 0, 0, 0, 1.0265, 3.8928177},
    {39, 2, 1104, 250, 0, 0, 1.03, -14.535256},
};
const Branch branches[] = {
    {1, 2, 0.0035, 0.0411, 0.6987, 0},   {1, 39, 0.001, 0.025, 0.75, 0},
    {2, 3, 0.0013, 0.0151, 0.2572, 0},   {2, 25, 0.007, 0.0086, 0.146, 0},
    {2, 30, 0, 0.0181, 0, 1.025},        {3, 4, 0.0013, 0.0213, 0.2214, 0},
    {3, 18, 0.0011, 0.0133, 0.2138, 0},  {4, 5, 0.0008, 0.0128, 0.1342, 0},
    {4, 14, 0.0008, 0.0129, 0.1382, 0},  {5, 6, 0.0002, 0.0026, 0.0434, 0},
    {5, 8, 0.0008, 0.0112, 0.1476, 0},   {6, 7, 0.0006, 0.0092, 0.113, 0},
    {6, 11, 0.0007, 0.0082, 0.1389, 0},  {6, 31, 0, 0.025, 0, 1.07},
    {7, 8, 0.0004, 0.0046, 0.078, 0},    {8, 9, 0.0023, 0.0363, 0.3804, 0},
    {9, 39, 0.001, 0.025, 1.2, 0},       {10, 11, 0.0004, 0.0043, 0.0729, 0},
    {10, 13, 0.0004, 0.0043, 0.0729, 0}, {10, 32, 0, 0.02, 0, 1.07},
    {12, 11, 0.0016, 0.0435, 0, 1.006},  {12, 13, 0.0016, 0.0435, 0, 1.006},
    {13, 14, 0.0009, 0.0101, 0.1723, 0}, {14, 15, 0.0018, 0.0217, 0.366, 0},
    {15, 16, 0.0009, 0.0094, 0.171, 0},  {16, 17, 0.0007, 0.0089, 0.1342, 0},
    {16, 19, 0.0016, 0.0195, 0.304, 0},  {16, 21, 0.0008, 0.0135, 0.2548, 0},
    {16, 24, 0.0003, 0.0059, 0.068, 0},  {17, 18, 0.0007, 0.0082, 0.1319, 0},
    {17, 27, 0.0013, 0.0173, 0.3216, 0}, {19, 20, 0.0007, 0.0138, 0, 1.06},
    {19, 33, 0.0007, 0.0142, 0, 1.07},   {20, 34, 0.0009, 0.018, 0, 1.009},
    {21, 22, 0.0008, 0.014, 0.2565, 0},  {22, 23, 0.0006, 0.0096, 0.1846, 0},
    {22, 35, 0, 0.0143, 0, 1.025},       {23, 24, 0.0022, 0.035, 0.361, 0},
    {23, 36, 0.0005, 0.0272, 0, 1},      {25, 26, 0.0032, 0.0323, 0.531, 0},
    {25, 37, 0.0006, 0.0232, 0, 1.025},  {26, 27, 0.0014, 0.0147, 0.2396, 0},
    {26, 28, 0.0043, 0.0474, 0.7802, 0}, {26, 29, 0.0057, 0.0625, 1.029, 0},
    {28, 29, 0.0014, 0.0151, 0.249, 0},  {29, 38, 0.0008, 0.0156, 0, 1.025},
};
const Generator generators[] = {
    {30, 250, 161.762, 1.0499, 10000, 500, 0, .003, 7.0, .05, .70, .035, .020,
     .006, .004, .019, .008, .004},
    {31, 677.871, 221.574, .982, 1000, 30.3, 0, .035, 6.56, .05, 1.5, .035,
     .295, .0697, .05, .282, .17, .05},
    {32, 650, 206.965, .9841, 1000, 35.8, 0, .0304, 5.7, .05, 1.5, .035, .2495,
     .0531, .045, .237, .0876, .045},
    {33, 632, 108.293, .9972, 1000, 28.6, 0, .0295, 5.69, .05, 1.5, .035, .262,
     .0436, .035, .258, .166, .035},
    {34, 508, 166.688, 1.0123, 600, 26.0, 0, .054, 5.4, .05, .44, .035, .67,
     .132, .089, .62, .166, .089},
    {35, 650, 210.661, 1.0494, 1000, 34.8, 0, .0224, 7.3, .05, .4, .035, .254,
     .05, .04, .241, .0814, .04},
    {36, 560, 100.165, 1.0636, 1000, 26.4, 0, .0322, 5.66, .05, 1.5, .035, .295,
     .049, .044, .292, .186, .044},
    {37, 540, -1.36945, 1.0275, 1000, 24.3, 0, .028, 6.7, .05, .41, .035, .29,
     .057, .045, .28, .0911, .045},
    {38, 830, 21.7327, 1.0265, 1000, 34.5, 0, .0298, 4.79, .05, 1.96, .035,
     .2106, .057, .045, .205, .0587, .045},
    {39, 1000, 78.4674, 1.03, 1000, 42.0, 0, .0125, 10.2, .05, .035, .035, .1,
     .031, .025, .069, .025, .015},
};

Matrix parameter3(Real value) {
  return Math::singlePhaseParameterToThreePhase(value);
}

PowerflowBusType pfBusType(int type) {
  if (type == 3)
    return PowerflowBusType::VD;
  if (type == 2)
    return PowerflowBusType::PV;
  return PowerflowBusType::PQ;
}

SystemTopology makePowerflowSystem() {
  SystemNodeList systemNodes;
  SystemComponentList components;
  std::vector<SP::SimNode::Ptr> nodes;
  for (const auto &bus : buses) {
    auto node =
        SP::SimNode::make("N" + std::to_string(bus.id), PhaseType::Single);
    nodes.push_back(node);
    systemNodes.push_back(node);
  }
  for (const auto &gen : generators) {
    const auto &bus = buses[gen.bus - 1];
    auto machine =
        SP::Ph1::SynchronGenerator::make("Gen_N" + std::to_string(gen.bus));
    const Real pg = gen.bus == 30 ? 0.5 * gen.pg : gen.pg;
    const Real qg = gen.bus == 30 ? 0.5 * gen.qg : gen.qg;
    machine->setParameters(gen.baseS * 1e6, baseVoltage, pg * 1e6,
                           gen.vg * baseVoltage, pfBusType(bus.type), qg * 1e6,
                           2000e6, -2000e6);
    machine->setBaseVoltage(baseVoltage);
    machine->connect({nodes.at(gen.bus - 1)});
    components.push_back(machine);
  }
  int loadIndex = 0;
  for (const auto &bus : buses) {
    if (bus.pd != 0. || bus.qd != 0.) {
      auto load = SP::Ph1::Load::make("PFLoad" + std::to_string(++loadIndex));
      load->setParameters(bus.pd * 1e6, bus.qd * 1e6, baseVoltage);
      load->connect({nodes.at(bus.id - 1)});
      components.push_back(load);
    }
    if (bus.gs != 0. || bus.bs != 0.) {
      auto shunt = SP::Ph1::Shunt::make("PFShunt_N" + std::to_string(bus.id));
      shunt->setParameters(bus.gs * 1e6 / (baseVoltage * baseVoltage),
                           bus.bs * 1e6 / (baseVoltage * baseVoltage));
      shunt->setBaseVoltage(baseVoltage);
      shunt->connect({nodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }
  int lineIndex = 0, transformerIndex = 0;
  for (const auto &branch : branches) {
    const Real tap = branch.tap == 0. ? 1. : branch.tap;
    const Real r = branch.r * zBase * tap * tap;
    const Real l = branch.x * zBase * tap * tap / omega;
    if (branch.tap == 0.) {
      auto line = SP::Ph1::PiLine::make("PFLine" + std::to_string(++lineIndex));
      line->setParameters(r, l, branch.b / zBase / omega, 0.);
      line->setBaseVoltage(baseVoltage);
      line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
      components.push_back(line);
    } else {
      auto trafo = SP::Ph1::Transformer::make(
          "PFTrafo" + std::to_string(++transformerIndex));
      trafo->setParameters(baseVoltage, baseVoltage, baseMVA * 1e6, tap, 0., r,
                           l);
      trafo->setBaseVoltage(baseVoltage);
      trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
      components.push_back(trafo);
    }
  }
  return SystemTopology(frequency, systemNodes, components);
}

void runInverter(const SystemTopology &systemPF, const std::string &mode,
                 const std::string &fmuPath, Real timeStep, Real finalTime,
                 Real loadStepTime) {
  const String simName = "IEEE39_Inverter_" + mode;
  Logger::setLogDir("logs/" + simName);
  std::vector<EMT::SimNode::Ptr> nodes;
  SystemNodeList systemNodes;
  for (const auto &bus : buses) {
    auto node =
        EMT::SimNode::make("N" + std::to_string(bus.id), PhaseType::ABC);
    nodes.push_back(node);
    systemNodes.push_back(node);
  }
  auto inverterNode = nodes.at(29);
  SystemComponentList components;

  // Retain half of the generator-30 operating point as a native machine; the
  // dual-interface converter supplies the other half at the same bus.
  for (const auto &gen : generators) {
    auto machine = EMT::Ph3::SynchronGenerator4OrderVBR::make(
        "Gen_N" + std::to_string(gen.bus));
    machine->setOperationalParametersPerUnit(
        gen.baseS * 1e6, baseVoltage, frequency, gen.H, gen.Xd, gen.Xq, gen.Xl,
        gen.Xd_t, gen.Xq_t, gen.Td0_t, gen.Tq0_t);
    machine->connect({nodes.at(gen.bus - 1)});
    components.push_back(machine);
  }

  int loadIndex = 0;
  for (const auto &bus : buses) {
    if (bus.pd == 0.0 && bus.qd == 0.0)
      continue;
    auto load = EMT::Ph3::RXLoad::make("load" + std::to_string(++loadIndex));
    load->setParameters(parameter3(bus.pd * 1e6 / 3.0),
                        parameter3(bus.qd * 1e6 / 3.0), bus.vm * baseVoltage);
    load->connect({nodes.at(bus.id - 1)});
    components.push_back(load);
    const Real g = bus.gs * 1e6 / (baseVoltage * baseVoltage);
    const Real b = bus.bs * 1e6 / (baseVoltage * baseVoltage);
    if (g != 0.0 || b != 0.0) {
      auto shunt = EMT::Ph3::Shunt::make("Shunt_N" + std::to_string(bus.id));
      shunt->setParameters(g, b);
      shunt->connect({nodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }
  for (const auto &bus : buses) {
    const Real g = bus.gs * 1e6 / (baseVoltage * baseVoltage);
    const Real b = bus.bs * 1e6 / (baseVoltage * baseVoltage);
    if ((g != 0.0 || b != 0.0) && bus.pd == 0.0 && bus.qd == 0.0) {
      auto shunt = EMT::Ph3::Shunt::make("Shunt_N" + std::to_string(bus.id));
      shunt->setParameters(g, b);
      shunt->connect({nodes.at(bus.id - 1)});
      components.push_back(shunt);
    }
  }

  int lineIndex = 0;
  int transformerIndex = 0;
  for (const auto &branch : branches) {
    const Real tap = branch.tap == 0.0 ? 1.0 : branch.tap;
    const Real r = branch.r * zBase * tap * tap;
    const Real l = branch.x * zBase * tap * tap / omega;
    if (branch.tap == 0.0) {
      const Real c = branch.b / zBase / omega;
      auto line = EMT::Ph3::PiLine::make("line" + std::to_string(++lineIndex));
      line->setParameters(parameter3(r), parameter3(l), parameter3(c),
                          parameter3(0.0));
      line->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
      components.push_back(line);
    } else {
      auto trafo = EMT::Ph3::Transformer::make(
          "transformer" + std::to_string(++transformerIndex));
      trafo->setParameters(baseVoltage, baseVoltage, baseMVA * 1e6, tap, 0.0,
                           parameter3(r), parameter3(l));
      trafo->connect({nodes.at(branch.from - 1), nodes.at(branch.to - 1)});
      components.push_back(trafo);
    }
  }

  FMI::NewtonSettings newtonSettings;
  // Keep the ME network coupling to a maximum of two Newton iterations.
  // NonlinearSSN uses the implicit trapezoidal state residual.
  newtonSettings.maxIterations = 2;
  newtonSettings.multirateMode = FMI::MultirateMode::SinglePassExplicit;
  FMI::CoSimulationElectricalPort port{
      {"v_a", "v_b", "v_c"}, {"i_a", "i_b", "i_c"}, 1.0};
  std::shared_ptr<SimPowerComp<Real>> inverter;
  if (mode == "me") {
    FMI::ModelExchangeOptions options;
    auto model = std::make_shared<FMI::ModelExchange>(fmuPath, options);
    inverter = FMI::NonlinearSSN::make(
        "Gen30_AveragedInverter", model,
        FMI::ElectricalPort3Ph{
            {"v_a", "v_b", "v_c"}, {"i_a", "i_b", "i_c"}, 1.0},
        newtonSettings);
  } else if (mode == "cs") {
    auto model = std::make_shared<FMI::CoSimulation>(fmuPath);
    inverter =
        FMI::CoSimulationSSN::make("Gen30_AveragedInverter", model, port);
  } else {
    throw std::invalid_argument("mode must be me or cs");
  }
  inverter->connect({EMT::SimNode::GND, inverterNode});
  components.push_back(inverter);

  // A 100 MW balanced load step at bus 39 is applied to all FMI modes at the
  // same absolute time. Its per-phase resistance is calculated from the
  // initialized bus voltage and three-phase P = 3*Vphase,rms^2/R.
  const Real stepPower = 100e6;
  const Real bus39Peak = buses[38].vm * baseVoltage * dpEnvelopeScale;
  const Real stepResistance = 1.5 * bus39Peak * bus39Peak / stepPower;
  auto stepNode = EMT::SimNode::make("N39_LoadStep", PhaseType::ABC);
  auto stepLoad = EMT::Ph3::Resistor::make("N39_LoadStep_R");
  stepLoad->setParameters(parameter3(stepResistance));
  stepLoad->connect({EMT::SimNode::GND, stepNode});
  auto loadSwitch = EMT::Ph3::Switch::make("N39_LoadStep_Switch");
  loadSwitch->setParameters(parameter3(1e12), parameter3(1e-3), false);
  loadSwitch->connect({stepNode, nodes.at(38)});
  systemNodes.push_back(stepNode);
  components.push_back(stepLoad);
  components.push_back(loadSwitch);

  auto system = SystemTopology(frequency, systemNodes, components);
  system.initWithPowerflow(systemPF, Domain::EMT);
  auto logger = DataLogger::make(simName, true, 1);
  for (int id : {1, 10, 20, 30, 31, 39})
    logger->logAttribute("N" + std::to_string(id) + ".V",
                         nodes.at(id - 1)->attribute("v"));
  logger->logAttribute("inverter.v", inverter->attribute("v_intf"));
  logger->logAttribute("inverter.i", inverter->attribute("i_intf"));
  logger->logAttribute("step_load.i", stepLoad->attribute("i_intf"));
  if (mode == "me") {
    logger->logAttribute("newton_iterations",
                         inverter->attribute("newton_iterations"));
    logger->logAttribute("newton_residual",
                         inverter->attribute("newton_residual"));
    logger->logAttribute("fmi_internal_steps",
                         inverter->attribute("fmi_internal_steps"));
    logger->logAttribute("fmi_evaluation_requests",
                         inverter->attribute("fmi_evaluation_requests"));
    logger->logAttribute("fmi_evaluation_wall_time_s",
                         inverter->attribute("fmi_evaluation_wall_time_s"));
  } else {
    logger->logAttribute("fmi_do_steps", inverter->attribute("cosim_do_steps"));
  }

  Simulation sim(simName, Logger::Level::off);
  sim.addEvent(SwitchEvent3Ph::make(loadStepTime, loadSwitch, true));
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.doInitFromNodesAndTerminals(true);
  sim.doSystemMatrixRecomputation(true);
  sim.addLogger(logger);
  const std::string logDirectory = "logs/" + simName;
  sim.setScheduler(std::make_shared<SequentialScheduler>(
      logDirectory + "/" + simName + "_TaskTimes.csv", Logger::Level::off));
  sim.setLogStepTimes(true);

  const auto initializeStart = std::chrono::steady_clock::now();
  sim.initialize();
  const std::chrono::duration<double> initializeElapsed =
      std::chrono::steady_clock::now() - initializeStart;
  const auto simulationStart = std::chrono::steady_clock::now();
  sim.run();
  const std::chrono::duration<double> simulationElapsed =
      std::chrono::steady_clock::now() - simulationStart;
  sim.logLUTimes();

  std::ofstream stepTiming(logDirectory + "/" + simName +
                           "_SimulationStepTimes.csv");
  stepTiming << "step,simulation_time_s,step_wall_s\n";
  const auto &stepTimes = sim.stepTimes();
  for (std::size_t index = 0; index < stepTimes.size(); ++index)
    stepTiming << index + 1 << "," << (index + 1) * timeStep << ","
               << stepTimes[index] << "\n";

  const Real totalRuntime =
      initializeElapsed.count() + simulationElapsed.count();
  std::cout << simName << " runtime_s=" << totalRuntime
            << " initialization_s=" << initializeElapsed.count()
            << " simulation_s=" << simulationElapsed.count()
            << " steps=" << std::llround(finalTime / timeStep)
            << " dt_s=" << timeStep << " stop_s=" << finalTime
            << " load_step_s=" << loadStepTime << " load_step_MW=100"
            << "\n";
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cerr << "Usage: IEEE39_AvgInverter me|cs inverter.fmu stop_s step_s "
                 "load_step_s\n";
    return 2;
  }
  try {
    const std::string mode = argv[1];
    const Real finalTime = std::stod(argv[3]);
    const Real timeStep = std::stod(argv[4]);
    const Real loadStepTime = std::stod(argv[5]);
    if (finalTime <= 0 || timeStep <= 0 || loadStepTime <= 0 ||
        loadStepTime >= finalTime ||
        std::abs(finalTime / timeStep - std::round(finalTime / timeStep)) >
            1e-8 ||
        std::abs(loadStepTime / timeStep -
                 std::round(loadStepTime / timeStep)) > 1e-8)
      throw std::invalid_argument("invalid time grid");

    auto systemPF = makePowerflowSystem();
    Simulation pf("IEEE39_Inverter_SP_PF", Logger::Level::off);
    pf.setSystem(systemPF);
    pf.setDomain(Domain::SP);
    pf.setSolverType(Solver::Type::NRP);
    pf.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
    pf.setTimeStep(.01);
    pf.setFinalTime(.01);
    pf.doInitFromNodesAndTerminals(false);
    pf.run();
    runInverter(systemPF, mode, argv[2], timeStep, finalTime, loadStepTime);
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
