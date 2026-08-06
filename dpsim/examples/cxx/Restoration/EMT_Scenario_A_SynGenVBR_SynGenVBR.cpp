// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorVBR.h>
#include <dpsim-models/Signal/Exciter.h>

#include <cmath>
#include <filesystem>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

namespace ScenarioASynGenVBRSynGenVBR {

// =============================================================================
// Parameters
// =============================================================================

struct SimulationParameters {
  String name = "EMT_Scenario_A_SynGenVBR_SynGenVBR";
  Real frequency = 50.0;
  Real timeStep = 100e-6;
  Real finalTime = 10.0;

  // Both generator islands are initialized with the breaker open. The breaker
  // is closed during the EMT simulation.
  Real breakerCloseTime = 2.0;

  Bool recomputeSystemMatrix = true;
};

struct SynGenVbrParameters {
  Real ratedPower = 0.0;
  Real ratedVoltage = 0.0;
  Real nominalFrequency = 50.0;
  Int poleNumber = 2;
  Real nominalFieldCurrent = 1300.0;

  // Operational synchronous-machine parameters in per unit.
  Real statorResistance = 0.0;
  Real ld = 0.0;
  Real lq = 0.0;
  Real ldTransient = 0.0;
  Real lqTransient = 0.0;
  Real ldSubtransient = 0.0;
  Real lqSubtransient = 0.0;
  Real leakageInductance = 0.0;
  Real td0Transient = 0.0;
  Real tq0Transient = 0.0;
  Real td0Subtransient = 0.0;
  Real tq0Subtransient = 0.0;
  Real inertia = 0.0;

  // Turbine-governor parameters.
  Real governorTa = 0.0;
  Real governorTb = 0.0;
  Real governorTc = 0.0;
  Real governorFa = 0.0;
  Real governorFb = 0.0;
  Real governorFc = 0.0;
  Real governorGain = 0.0;
  Real governorSpeedRelayTimeConstant = 0.0;
  Real governorServoMotorTimeConstant = 0.0;

  // Exciter parameters.
  Real exciterTa = 0.0;
  Real exciterKa = 0.0;
  Real exciterTe = 0.0;
  Real exciterKe = 0.0;
  Real exciterTf = 0.0;
  Real exciterKf = 0.0;
  Real exciterTr = 0.0;
  Real exciterMaximumRegulatorVoltage = 0.0;
  Real exciterMinimumRegulatorVoltage = 0.0;

  static SynGenVbrParameters gas() {
    SynGenVbrParameters p;

    p.ratedPower = 50e6;
    p.ratedVoltage = 10.5e3;
    p.nominalFrequency = 50.0;
    p.poleNumber = 2;
    p.nominalFieldCurrent = 1300.0;

    p.statorResistance = 0.002;
    p.ld = 2.4;
    p.lq = 1.33;
    p.ldTransient = 0.31;
    p.lqTransient = 1.2;
    p.ldSubtransient = 0.24;
    p.lqSubtransient = 0.35;
    p.leakageInductance = 0.135;
    p.td0Transient = 1.45;
    p.tq0Transient = 1.0e-6;
    p.td0Subtransient = 0.022;
    p.tq0Subtransient = 0.0095;
    p.inertia = 5.0;

    p.governorTa = 1.0;
    p.governorTb = 0.5;
    p.governorTc = 0.2;
    p.governorFa = 0.3;
    p.governorFb = 0.25;
    p.governorFc = 0.3;
    p.governorGain = 20.0;
    p.governorSpeedRelayTimeConstant = 0.1;
    p.governorServoMotorTimeConstant = 0.1;

    p.exciterTa = 0.005;
    p.exciterKa = 200.0;
    p.exciterTe = 0.05;
    p.exciterKe = 0.5;
    p.exciterTf = 0.3;
    p.exciterKf = 0.01;
    p.exciterTr = 0.02;
    p.exciterMaximumRegulatorVoltage = 9.0;
    p.exciterMinimumRegulatorVoltage = -5.0;

    return p;
  }

  static SynGenVbrParameters psh() {
    SynGenVbrParameters p;

    p.ratedPower = 200e6;
    p.ratedVoltage = 18.0e3;
    p.nominalFrequency = 50.0;
    p.poleNumber = 2;
    p.nominalFieldCurrent = 1300.0;

    p.statorResistance = 0.002;
    p.ld = 2.0;
    p.lq = 2.0;
    p.ldTransient = 0.3;
    p.lqTransient = 0.3;
    p.ldSubtransient = 0.2;
    p.lqSubtransient = 0.2;
    p.leakageInductance = 0.1;
    p.td0Transient = 1.0;
    p.tq0Transient = 1.0;
    p.td0Subtransient = 0.05;
    p.tq0Subtransient = 0.05;
    p.inertia = 4.0;

    p.governorTa = 0.3;
    p.governorTb = 7.0;
    p.governorTc = 0.5;
    p.governorFa = 0.25;
    p.governorFb = 0.3;
    p.governorFc = 0.15;
    p.governorGain = 25.0;
    p.governorSpeedRelayTimeConstant = 0.1;
    p.governorServoMotorTimeConstant = 0.1;

    p.exciterTa = 0.005;
    p.exciterKa = 200.0;
    p.exciterTe = 0.05;
    p.exciterKe = 0.5;
    p.exciterTf = 0.3;
    p.exciterKf = 0.01;
    p.exciterTr = 0.02;
    p.exciterMaximumRegulatorVoltage = 9.0;
    p.exciterMinimumRegulatorVoltage = -5.0;

    return p;
  }
};

struct GasTransformerParameters {
  Real nominalVoltageLow = 10.5e3;
  Real nominalVoltageHigh = 220e3;
  Real ratedPower = 50e6;
  Real ratioMagnitude = nominalVoltageLow / nominalVoltageHigh;
  Real ratioPhase = 0.0;

  Real resistance() const {
    const Real baseImpedance =
        nominalVoltageHigh * nominalVoltageHigh / ratedPower;
    return 0.00376196 * baseImpedance;
  }

  Real inductance(Real frequency) const {
    const Real baseImpedance =
        nominalVoltageHigh * nominalVoltageHigh / ratedPower;
    return 0.1007298 * baseImpedance / (2.0 * PI * frequency);
  }
};

struct PshTransformerParameters {
  Real nominalVoltageLow = 18.0e3;
  Real nominalVoltageHigh = 220e3;
  Real ratedPower = 200e6;
  Real ratioMagnitude = nominalVoltageLow / nominalVoltageHigh;
  Real ratioPhase = 0.0;

  Real resistance() const { return 0.41745 * 2.0; }

  Real inductance() const { return 0.0481260775594524 * 2.0; }
};

struct CableParameters {
  Real baseVoltage = 220e3;
  Real resistance = 5.04669;
  Real inductance = 0.13523123641;
  Real capacitance = 1.93522865e-6;
  Real conductance = 1e-15;
};

struct Line3Parameters {
  Real baseVoltage = 220e3;
  Real length = 5.0;

  Real resistance() const { return 0.0749 * length; }
  Real inductance() const { return 1.270693e-3 * length; }
  Real capacitance() const { return 0.00466961e-6 * length; }

  Real conductance = 1e-15;
};

struct BreakerParameters {
  Real openResistance = 1e12;
  // Retain the numerically robust value from the already working breaker
  // example. The old Scenario A draft declared 1e-8 Ohm, but that breaker
  // branch was commented out and therefore not validated there.
  Real closedResistance = 1e-3;
};

struct Parameters {
  SimulationParameters simulation;

  SynGenVbrParameters gasGenerator = SynGenVbrParameters::gas();
  SynGenVbrParameters pshGenerator = SynGenVbrParameters::psh();

  GasTransformerParameters gasTransformer;
  PshTransformerParameters pshTransformer;
  CableParameters cable;
  Line3Parameters line3;
  BreakerParameters breaker;

  // The initial angle is applied to the isolated PSH island in the SP
  // initialization and can be changed to study out-of-phase closing.
  Real pshInitialAngleDegrees = 0.0;
};

// =============================================================================
// Power-flow / steady-state initialization
// =============================================================================

struct PowerFlowResult {
  SystemTopology system;

  Complex gasPower;
  Complex pshPower;

  Complex gasBusVoltage;
  Complex pshBusVoltage;

  Complex breakerGridSideVoltage;
  Complex breakerPshSideVoltage;
};

Complex
generatorPositivePower(const std::shared_ptr<SP::Ph1::NetworkInjection> &source,
                       const SimNode<Complex>::Ptr &sourceBus) {
  const Complex voltage = sourceBus->singleVoltage();

  // NetworkInjection uses the component-consumer convention. Generator-positive
  // injected power therefore has the opposite sign.
  const Complex consumerCurrent = (**source->mIntfCurrent)(0, 0);
  return -voltage * std::conj(consumerCurrent);
}

PowerFlowResult buildAndRunPowerFlow(const Parameters &p) {
  const String simulationName = p.simulation.name + "_PF";

  SPDLOG_INFO("Scenario A bases:"
              "\n  GEN_gas: S_rated={} VA, V_rated={} V_LL RMS, H={} s"
              "\n  GEN_psh: S_rated={} VA, V_rated={} V_LL RMS, H={} s"
              "\n  TR_gas : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
              "\n  TR_psh : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
              "\n  breaker: R_open={} Ohm, R_closed={} Ohm",
              p.gasGenerator.ratedPower, p.gasGenerator.ratedVoltage,
              p.gasGenerator.inertia, p.pshGenerator.ratedPower,
              p.pshGenerator.ratedVoltage, p.pshGenerator.inertia,
              p.gasTransformer.nominalVoltageLow,
              p.gasTransformer.nominalVoltageHigh, p.gasTransformer.ratedPower,
              p.gasTransformer.resistance(),
              p.gasTransformer.inductance(p.simulation.frequency),
              p.pshTransformer.nominalVoltageLow,
              p.pshTransformer.nominalVoltageHigh, p.pshTransformer.ratedPower,
              p.pshTransformer.resistance(), p.pshTransformer.inductance(),
              p.breaker.openResistance, p.breaker.closedResistance);

  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  // Node names intentionally match the EMT topology so that
  // SystemTopology::initWithPowerflow() can transfer all phasors.
  auto busGas = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busB = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busA = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busPsha = SimNode<Complex>::make("BUS_psha", PhaseType::Single);
  auto busTrPsh = SimNode<Complex>::make("BUS_tr_psh", PhaseType::Single);
  auto busPsh = SimNode<Complex>::make("BUS_psh", PhaseType::Single);

  auto gasInjection =
      SP::Ph1::NetworkInjection::make("GEN_gas_PF", Logger::Level::debug);
  gasInjection->setParameters(Math::polar(p.gasGenerator.ratedVoltage, 0.0),
                              p.simulation.frequency);

  auto pshInjection =
      SP::Ph1::NetworkInjection::make("GEN_psh_PF", Logger::Level::debug);
  pshInjection->setParameters(
      Math::polar(p.pshGenerator.ratedVoltage,
                  p.pshInitialAngleDegrees * PI / 180.0),
      p.simulation.frequency);

  auto gasTransformer = SP::Ph1::Transformer::make("TR_gas_PF", "TR_gas_PF",
                                                   Logger::Level::debug, true);
  gasTransformer->setParameters(
      p.gasTransformer.nominalVoltageLow, p.gasTransformer.nominalVoltageHigh,
      p.gasTransformer.ratedPower, p.gasTransformer.ratioMagnitude,
      p.gasTransformer.ratioPhase, p.gasTransformer.resistance(),
      p.gasTransformer.inductance(p.simulation.frequency));

  auto cable = SP::Ph1::PiLine::make("cable_PF", Logger::Level::debug);
  cable->setParameters(p.cable.resistance, p.cable.inductance,
                       p.cable.capacitance, p.cable.conductance);

  auto line3 = SP::Ph1::PiLine::make("line_3_PF", Logger::Level::debug);
  line3->setParameters(p.line3.resistance(), p.line3.inductance(),
                       p.line3.capacitance(), p.line3.conductance);

  auto breaker =
      SP::Ph1::Switch::make("breaker_GEN_psh_PF", Logger::Level::debug);
  breaker->setParameters(p.breaker.openResistance, p.breaker.closedResistance,
                         false);

  auto pshTransformer = SP::Ph1::Transformer::make("TR_psh_PF", "TR_psh_PF",
                                                   Logger::Level::debug, true);
  pshTransformer->setParameters(
      p.pshTransformer.nominalVoltageLow, p.pshTransformer.nominalVoltageHigh,
      p.pshTransformer.ratedPower, p.pshTransformer.ratioMagnitude,
      p.pshTransformer.ratioPhase, p.pshTransformer.resistance(),
      p.pshTransformer.inductance());

  gasInjection->connect({busGas});
  gasTransformer->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPsha});
  breaker->connect({busPsha, busTrPsh});
  pshTransformer->connect({busPsh, busTrPsh});
  pshInjection->connect({busPsh});

  SystemTopology system(
      p.simulation.frequency,
      SystemNodeList{busGas, busB, busA, busPsha, busTrPsh, busPsh},
      SystemComponentList{gasInjection, gasTransformer, cable, line3, breaker,
                          pshTransformer, pshInjection});

  auto logger = DataLogger::make(simulationName);

  logger->logAttribute("V_BUS_gas_PF", busGas->attribute("v"));
  logger->logAttribute("V_BUS_b_PF", busB->attribute("v"));
  logger->logAttribute("V_BUS_a_PF", busA->attribute("v"));
  logger->logAttribute("V_BUS_psha_PF", busPsha->attribute("v"));
  logger->logAttribute("V_BUS_tr_psh_PF", busTrPsh->attribute("v"));
  logger->logAttribute("V_BUS_psh_PF", busPsh->attribute("v"));

  logger->logAttribute("I_GEN_gas_PF", gasInjection->attribute("i_intf"));
  logger->logAttribute("I_GEN_psh_PF", pshInjection->attribute("i_intf"));
  logger->logAttribute("I_breaker_GEN_psh_PF", breaker->attribute("i_intf"));

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(1.0);
  simulation.setFinalTime(1.0);
  simulation.setDomain(Domain::SP);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.addLogger(logger);
  simulation.run();

  const Complex gasPower = generatorPositivePower(gasInjection, busGas);
  const Complex pshPower = generatorPositivePower(pshInjection, busPsh);

  const Complex gasVoltage = busGas->singleVoltage();
  const Complex pshVoltage = busPsh->singleVoltage();
  const Complex breakerGridVoltage = busPsha->singleVoltage();
  const Complex breakerPshVoltage = busTrPsh->singleVoltage();

  const auto finiteComplex = [](const Complex &value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  };

  if (!(std::abs(gasVoltage) > 0.0) || !(std::abs(pshVoltage) > 0.0) ||
      !finiteComplex(gasVoltage) || !finiteComplex(pshVoltage) ||
      !finiteComplex(gasPower) || !finiteComplex(pshPower) ||
      !finiteComplex(breakerGridVoltage) || !finiteComplex(breakerPshVoltage)) {
    throw std::runtime_error(
        "Invalid SP steady-state result for Scenario A with two "
        "SynchronGeneratorVBR components.");
  }

  SPDLOG_INFO(
      "Solved open-breaker SP initialization:"
      "\n  GEN_gas: P={} W, Q={} var, |V|={} V_LL RMS, angle={} deg"
      "\n  GEN_psh: P={} W, Q={} var, |V|={} V_LL RMS, angle={} deg"
      "\n  breaker grid side: |V|={} V_LL RMS, angle={} deg"
      "\n  breaker PSH side : |V|={} V_LL RMS, angle={} deg"
      "\n  breaker mismatch : |dV|={} V_LL RMS, dAngle={} deg",
      gasPower.real(), gasPower.imag(), std::abs(gasVoltage),
      std::arg(gasVoltage) * 180.0 / PI, pshPower.real(), pshPower.imag(),
      std::abs(pshVoltage), std::arg(pshVoltage) * 180.0 / PI,
      std::abs(breakerGridVoltage), std::arg(breakerGridVoltage) * 180.0 / PI,
      std::abs(breakerPshVoltage), std::arg(breakerPshVoltage) * 180.0 / PI,
      std::abs(breakerGridVoltage - breakerPshVoltage),
      std::arg(breakerGridVoltage / breakerPshVoltage) * 180.0 / PI);

  return {system,     gasPower,           pshPower,         gasVoltage,
          pshVoltage, breakerGridVoltage, breakerPshVoltage};
}

// =============================================================================
// EMT model
// =============================================================================

void configureGenerator(
    const std::shared_ptr<EMT::Ph3::SynchronGeneratorVBR> &generator,
    const SynGenVbrParameters &parameters, const Complex &initialGeneratorPower,
    const String &exciterName) {
  generator->setBaseAndOperationalPerUnitParameters(
      parameters.ratedPower, parameters.ratedVoltage,
      parameters.nominalFrequency, parameters.poleNumber,
      parameters.nominalFieldCurrent, parameters.statorResistance,
      parameters.ld, parameters.lq, parameters.ldTransient,
      parameters.lqTransient, parameters.ldSubtransient,
      parameters.lqSubtransient, parameters.leakageInductance,
      parameters.td0Transient, parameters.tq0Transient,
      parameters.td0Subtransient, parameters.tq0Subtransient,
      parameters.inertia);

  // The dynamic coefficients are taken from the original Scenario A model.
  // The operating point is taken from the solved open-breaker SP system so
  // that neither governor receives an artificial power step at t = 0.
  const Real mechanicalPowerPerUnit =
      initialGeneratorPower.real() / parameters.ratedPower;

  generator->addGovernor(
      parameters.governorTa, parameters.governorTb, parameters.governorTc,
      parameters.governorFa, parameters.governorFb, parameters.governorFc,
      parameters.governorGain, parameters.governorSpeedRelayTimeConstant,
      parameters.governorServoMotorTimeConstant, mechanicalPowerPerUnit,
      mechanicalPowerPerUnit);

  auto exciterParameters = std::make_shared<Signal::ExciterParameters>();
  exciterParameters->Ta = parameters.exciterTa;
  exciterParameters->Ka = parameters.exciterKa;
  exciterParameters->Te = parameters.exciterTe;
  exciterParameters->Ke = parameters.exciterKe;
  exciterParameters->Tf = parameters.exciterTf;
  exciterParameters->Kf = parameters.exciterKf;
  exciterParameters->Tr = parameters.exciterTr;
  exciterParameters->maxVr = parameters.exciterMaximumRegulatorVoltage;
  exciterParameters->minVr = parameters.exciterMinimumRegulatorVoltage;

  auto exciter = Signal::Exciter::make(exciterName, Logger::Level::debug);
  generator->addExciter(exciter, exciterParameters);

  SPDLOG_INFO(
      "{} parameters:"
      "\n  S_rated={} VA, V_rated={} V_LL RMS, f_nom={} Hz"
      "\n  P_init={} W, Q_init={} var, Pm_init={} pu"
      "\n  Rs={} pu, Ld={} pu, Lq={} pu, Ld'={} pu, Lq'={} pu"
      "\n  Ld''={} pu, Lq''={} pu, Ll={} pu, H={} s"
      "\n  Td0'={} s, Tq0'={} s, Td0''={} s, Tq0''={} s"
      "\n  governor: Ta={}, Tb={}, Tc={}, Fa={}, Fb={}, Fc={}, K={}, "
      "Tsr={}, Tsm={}"
      "\n  exciter: Ta={}, Ka={}, Te={}, Ke={}, Tf={}, Kf={}, Tr={}, "
      "Vr=[{},{}]",
      generator->name(), parameters.ratedPower, parameters.ratedVoltage,
      parameters.nominalFrequency, initialGeneratorPower.real(),
      initialGeneratorPower.imag(), mechanicalPowerPerUnit,
      parameters.statorResistance, parameters.ld, parameters.lq,
      parameters.ldTransient, parameters.lqTransient, parameters.ldSubtransient,
      parameters.lqSubtransient, parameters.leakageInductance,
      parameters.inertia, parameters.td0Transient, parameters.tq0Transient,
      parameters.td0Subtransient, parameters.tq0Subtransient,
      parameters.governorTa, parameters.governorTb, parameters.governorTc,
      parameters.governorFa, parameters.governorFb, parameters.governorFc,
      parameters.governorGain, parameters.governorSpeedRelayTimeConstant,
      parameters.governorServoMotorTimeConstant, parameters.exciterTa,
      parameters.exciterKa, parameters.exciterTe, parameters.exciterKe,
      parameters.exciterTf, parameters.exciterKf, parameters.exciterTr,
      parameters.exciterMinimumRegulatorVoltage,
      parameters.exciterMaximumRegulatorVoltage);
}

void runEmt(const Parameters &p, const PowerFlowResult &powerFlow) {
  const String simulationName = p.simulation.name + "_EMT";

  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busA = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto busPsha = SimNode<Real>::make("BUS_psha", PhaseType::ABC);
  auto busTrPsh = SimNode<Real>::make("BUS_tr_psh", PhaseType::ABC);
  auto busPsh = SimNode<Real>::make("BUS_psh", PhaseType::ABC);

  auto gasGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_gas", Logger::Level::debug);
  auto pshGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_psh", Logger::Level::debug);

  configureGenerator(gasGenerator, p.gasGenerator, powerFlow.gasPower,
                     "GEN_gas_exciter");
  configureGenerator(pshGenerator, p.pshGenerator, powerFlow.pshPower,
                     "GEN_psh_exciter");

  auto gasTransformer = EMT::Ph3::Transformer::make("TR_gas", "TR_gas",
                                                    Logger::Level::debug, true);
  gasTransformer->setParameters(
      p.gasTransformer.nominalVoltageLow, p.gasTransformer.nominalVoltageHigh,
      p.gasTransformer.ratedPower, p.gasTransformer.ratioMagnitude,
      p.gasTransformer.ratioPhase,
      Math::singlePhaseParameterToThreePhase(p.gasTransformer.resistance()),
      Math::singlePhaseParameterToThreePhase(
          p.gasTransformer.inductance(p.simulation.frequency)));

  auto cable = EMT::Ph3::PiLine::make("cable", Logger::Level::debug);
  cable->setParameters(
      Math::singlePhaseParameterToThreePhase(p.cable.resistance),
      Math::singlePhaseParameterToThreePhase(p.cable.inductance),
      Math::singlePhaseParameterToThreePhase(p.cable.capacitance),
      Math::singlePhaseParameterToThreePhase(p.cable.conductance));

  auto line3 = EMT::Ph3::PiLine::make("line_3", Logger::Level::debug);
  line3->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line3.resistance()),
      Math::singlePhaseParameterToThreePhase(p.line3.inductance()),
      Math::singlePhaseParameterToThreePhase(p.line3.capacitance()),
      Math::singlePhaseParameterToThreePhase(p.line3.conductance));

  auto breaker =
      EMT::Ph3::Switch::make("breaker_GEN_psh", Logger::Level::debug);
  breaker->setParameters(
      Math::singlePhaseParameterToThreePhase(p.breaker.openResistance),
      Math::singlePhaseParameterToThreePhase(p.breaker.closedResistance));
  breaker->openSwitch();

  auto pshTransformer = EMT::Ph3::Transformer::make("TR_psh", "TR_psh",
                                                    Logger::Level::debug, true);
  pshTransformer->setParameters(
      p.pshTransformer.nominalVoltageLow, p.pshTransformer.nominalVoltageHigh,
      p.pshTransformer.ratedPower, p.pshTransformer.ratioMagnitude,
      p.pshTransformer.ratioPhase,
      Math::singlePhaseParameterToThreePhase(p.pshTransformer.resistance()),
      Math::singlePhaseParameterToThreePhase(p.pshTransformer.inductance()));

  gasGenerator->connect({busGas});
  gasTransformer->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPsha});
  breaker->connect({busPsha, busTrPsh});
  pshTransformer->connect({busPsh, busTrPsh});
  pshGenerator->connect({busPsh});

  SystemTopology system(
      p.simulation.frequency,
      SystemNodeList{busGas, busB, busA, busPsha, busTrPsh, busPsh},
      SystemComponentList{gasGenerator, gasTransformer, cable, line3, breaker,
                          pshTransformer, pshGenerator});

  // The open-breaker SP topology contains both islands and identically named
  // nodes. This initializes both synchronous machines and both breaker sides.
  system.initWithPowerflow(powerFlow.system, Domain::EMT);

  // SynchronGeneratorVBR expects terminal power in DPsim's consumer-positive
  // convention and negates it during initializeFromNodesAndTerminals().
  gasGenerator->terminal(0)->setPower(-powerFlow.gasPower);
  pshGenerator->terminal(0)->setPower(-powerFlow.pshPower);

  auto logger = DataLogger::make(simulationName);

  logger->logAttribute("BUS_gas_v", busGas->attribute("v"));
  logger->logAttribute("BUS_b_v", busB->attribute("v"));
  logger->logAttribute("BUS_a_v", busA->attribute("v"));
  logger->logAttribute("BUS_psha_v", busPsha->attribute("v"));
  logger->logAttribute("BUS_tr_psh_v", busTrPsh->attribute("v"));
  logger->logAttribute("BUS_psh_v", busPsh->attribute("v"));

  logger->logAttribute("TR_gas_i", gasTransformer->attribute("i_intf"));
  logger->logAttribute("cable_i", cable->attribute("i_intf"));
  logger->logAttribute("line_3_i", line3->attribute("i_intf"));
  logger->logAttribute("breaker_GEN_psh_i", breaker->attribute("i_intf"));
  logger->logAttribute("breaker_GEN_psh_v", breaker->attribute("v_intf"));
  logger->logAttribute("TR_psh_i", pshTransformer->attribute("i_intf"));

  logger->logAttribute("GEN_gas_w_r", gasGenerator->attribute("w_r"));
  logger->logAttribute("GEN_gas_P_elec", gasGenerator->attribute("P_elec"));
  logger->logAttribute("GEN_gas_Q_elec", gasGenerator->attribute("Q_elec"));
  logger->logAttribute("GEN_gas_P_mech", gasGenerator->attribute("P_mech"));
  logger->logAttribute("GEN_gas_i_intf", gasGenerator->attribute("i_intf"));
  logger->logAttribute("GEN_gas_v_intf", gasGenerator->attribute("v_intf"));
  logger->logAttribute("GEN_gas_delta_r", gasGenerator->attribute("delta_r"));
  logger->logAttribute("GEN_gas_T_e", gasGenerator->attribute("T_e"));

  logger->logAttribute("GEN_psh_w_r", pshGenerator->attribute("w_r"));
  logger->logAttribute("GEN_psh_P_elec", pshGenerator->attribute("P_elec"));
  logger->logAttribute("GEN_psh_Q_elec", pshGenerator->attribute("Q_elec"));
  logger->logAttribute("GEN_psh_P_mech", pshGenerator->attribute("P_mech"));
  logger->logAttribute("GEN_psh_i_intf", pshGenerator->attribute("i_intf"));
  logger->logAttribute("GEN_psh_v_intf", pshGenerator->attribute("v_intf"));
  logger->logAttribute("GEN_psh_delta_r", pshGenerator->attribute("delta_r"));
  logger->logAttribute("GEN_psh_T_e", pshGenerator->attribute("T_e"));

  Simulation simulation(simulationName, Logger::Level::info);

  simulation.addEvent(
      SwitchEvent3Ph::make(p.simulation.breakerCloseTime, breaker, true));

  simulation.setSystem(system);
  simulation.setTimeStep(p.simulation.timeStep);
  simulation.setFinalTime(p.simulation.finalTime);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.doSystemMatrixRecomputation(p.simulation.recomputeSystemMatrix);
  simulation.addLogger(logger);
  simulation.run();
}

} // namespace ScenarioASynGenVBRSynGenVBR

int main(int argc, char *argv[]) {
  ScenarioASynGenVBRSynGenVBR::Parameters parameters;

  CommandLineArgs args(argc, argv);

  if (argc > 1) {
    parameters.simulation.timeStep = args.timeStep;
    parameters.simulation.finalTime = args.duration;

    if (args.name != "dpsim")
      parameters.simulation.name = args.name;

    if (args.options.find("BREAKER_CLOSE_TIME") != args.options.end()) {
      parameters.simulation.breakerCloseTime =
          args.getOptionReal("BREAKER_CLOSE_TIME");
    }

    if (args.options.find("PSH_INITIAL_ANGLE_DEG") != args.options.end()) {
      parameters.pshInitialAngleDegrees =
          args.getOptionReal("PSH_INITIAL_ANGLE_DEG");
    } else if (args.options.find("PSHA_INITIAL_ANGLE_DEG") !=
               args.options.end()) {
      // Backward-compatible option name used by the GFM/SynGen example.
      parameters.pshInitialAngleDegrees =
          args.getOptionReal("PSHA_INITIAL_ANGLE_DEG");
    }
  }

  if (!(parameters.simulation.timeStep > 0.0) ||
      !(parameters.simulation.finalTime > 0.0) ||
      !(parameters.simulation.breakerCloseTime >= 0.0) ||
      !(parameters.simulation.breakerCloseTime <
        parameters.simulation.finalTime)) {
    throw std::invalid_argument("Require dt>0, finalTime>0, and "
                                "0<=breakerCloseTime<finalTime.");
  }

  const auto powerFlow =
      ScenarioASynGenVBRSynGenVBR::buildAndRunPowerFlow(parameters);

  ScenarioASynGenVBRSynGenVBR::runEmt(parameters, powerFlow);

  return 0;
}
