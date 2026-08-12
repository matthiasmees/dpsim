// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <dpsim-models/EMT/EMT_Ph3_RXLoad.h>
#include <dpsim-models/EMT/EMT_Ph3_Switch.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorVBR.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/SP/SP_Ph1_Load.h>
#include <dpsim-models/SP/SP_Ph1_Switch.h>
#include <dpsim-models/SP/SP_Ph1_SynchronGenerator.h>
#include <dpsim-models/SP/SP_Ph1_Transformer.h>
#include <dpsim-models/Signal/Exciter.h>

#include <memory>

using namespace DPsim;
using namespace CPS;

struct Parameters {
  // Simulation
  String name = "EMT_Load_loss_SynGenVBR_SynGenVBR";
  Real frequency = 50.0;
  Real pfTimeStep = 1.0;
  Real pfFinalTime = 1.0;
  Real timeStep = 100e-6;
  Real finalTime = 15.0;
  Real breakerOpenTime = 5.0;
  Bool recomputeSystemMatrix = true;
  Int logDownSampling = 10;

  // Synchronous generator GEN_gas
  Real genRatedPower = 50e6;
  Real genRatedVoltage = 10.5e3;
  Real genNominalFrequency = 50.0;
  Int genPoleNumber = 2;
  Real genNominalFieldCurrent = 1300.0;

  Real genStatorResistance = 0.002;
  Real genLd = 2.4;
  Real genLq = 1.33;
  Real genLdTransient = 0.31;
  Real genLqTransient = 1.2;
  Real genLdSubtransient = 0.24;
  Real genLqSubtransient = 0.35;
  Real genLeakageInductance = 0.135;
  Real genTd0Transient = 1.45;
  Real genTq0Transient = 1.0;
  Real genTd0Subtransient = 0.022;
  Real genTq0Subtransient = 0.0095;
  Real genInertia = 5.0;

  // Turbine governor
  Real govTa = 1.0;
  Real govTb = 0.5;
  Real govTc = 0.2;
  Real govF1a = 0.3;
  Real govFa = 0.3;
  Real govFb = 0.25;
  Real govFc = 0.3;
  Real govGain = 20.0;
  Real govSpeedRelayTimeConstant = 0.1;
  Real govServoMotorTimeConstant = 0.1;

  // Exciter
  Real excTa = 0.005;
  Real excKa = 200.0;
  Real excTe = 0.05;
  Real excKe = 0.5;
  Real excTf = 0.3;
  Real excKf = 0.01;
  Real excTr = 0.02;
  Real excMaximumRegulatorVoltage = 19.0;
  Real excMinimumRegulatorVoltage = -5.0;

  // Transformer TR_gas
  Real trGasVoltageLow = 10.5e3;
  Real trGasVoltageHigh = 220e3;
  Real trGasRatedPower = 50e6;
  Real trGasRatioMagnitude = 10.5e3 / 220e3;
  Real trGasRatioPhase = 0.0;
  Real trGasResistance = 3.64157728;
  Real trGasInductance = 0.31037265855769886;

  // Transformer TR_ld2
  Real trLoadVoltageLow = 10.0e3;
  Real trLoadVoltageHigh = 220e3;
  Real trLoadRatedPower = 50e6;
  Real trLoadRatioMagnitude = 10.0e3 / 220e3;
  Real trLoadRatioPhase = 0.0;
  Real trLoadResistance = 3.64157728;
  Real trLoadInductance = 0.31037292;

  // LOAD_2a
  Real load2aActivePower = 5e6;
  Real load2aReactivePower = 3.098721e6;
  Real load2aNominalVoltage = 10e3;

  // LOAD_2b
  Real load2bActivePower = 15e6;
  Real load2bReactivePower = 0.0;
  Real load2bNominalVoltage = 10e3;

  // BR_load_2b
  Real breakerOpenResistance = 1e12;
  Real breakerClosedResistance = 1e-3;
};

int main() {
  const Parameters p;

  // ==========================================================================
  // 1. POWER FLOW
  //
  // PF and EMT have exactly the same graph:
  //
  // GEN_gas -- BUS_gas -- TR_gas -- BUS_b -- TR_ld2 -- BUS_load_2
  //                                                       |          |
  //                                                    LOAD_2a   BR_load_2b
  //                                                                  |
  //                                                           BUS_load_2b
  //                                                                  |
  //                                                               LOAD_2b
  //
  // PF and EMT use identical component and node names. The two systems are
  // separate SystemTopology objects, so this is safe and lets the stock
  // initWithPowerflow() exact-name lookup initialize GEN_gas correctly.
  // ==========================================================================

  const String pfSimName = p.name + "_PF";
  Logger::setLogDir("logs/" + pfSimName);

  // PF nodes
  auto busGasPF = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busBPF = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busLoadPF = SimNode<Complex>::make("BUS_load_2", PhaseType::Single);
  auto busLoad2bPF = SimNode<Complex>::make("BUS_load_2b", PhaseType::Single);

  // PF synchronous generator: VD/slack bus. P is only an initial estimate;
  // the converged slack P/Q are solved by the NRP solver.
  auto gasGeneratorPF =
      SP::Ph1::SynchronGenerator::make("GEN_gas", Logger::Level::off);
  gasGeneratorPF->setParameters(p.genRatedPower, p.genRatedVoltage,
                                p.load2aActivePower + p.load2bActivePower,
                                p.genRatedVoltage, PowerflowBusType::VD);
  gasGeneratorPF->setBaseVoltage(p.genRatedVoltage);

  // PF TR_gas
  auto trGasPF =
      SP::Ph1::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  trGasPF->setParameters(p.trGasVoltageLow, p.trGasVoltageHigh,
                         p.trGasRatedPower, p.trGasRatioMagnitude,
                         p.trGasRatioPhase, p.trGasResistance,
                         p.trGasInductance);
  trGasPF->setBaseVoltage(p.trGasVoltageHigh);

  // PF TR_ld2
  auto trLoadPF =
      SP::Ph1::Transformer::make("TR_ld2", "TR_ld2", Logger::Level::off, true);
  trLoadPF->setParameters(p.trLoadVoltageLow, p.trLoadVoltageHigh,
                          p.trLoadRatedPower, p.trLoadRatioMagnitude,
                          p.trLoadRatioPhase, p.trLoadResistance,
                          p.trLoadInductance);
  trLoadPF->setBaseVoltage(p.trLoadVoltageHigh);

  // PF LOAD_2a
  auto load2aPF = SP::Ph1::Load::make("LOAD_2a", Logger::Level::off);
  load2aPF->setParameters(p.load2aActivePower, p.load2aReactivePower,
                          p.load2aNominalVoltage);
  load2aPF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // PF breaker: actual SP::Ph1::Switch, initially closed
  auto breakerPF = SP::Ph1::Switch::make("BR_load_2b", Logger::Level::off);
  breakerPF->setParameters(p.breakerOpenResistance, p.breakerClosedResistance,
                           true);
  breakerPF->setBaseVoltage(p.load2bNominalVoltage);

  // PF LOAD_2b behind the breaker
  auto load2bPF = SP::Ph1::Load::make("LOAD_2b", Logger::Level::off);
  load2bPF->setParameters(p.load2bActivePower, p.load2bReactivePower,
                          p.load2bNominalVoltage);
  load2bPF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // PF connectivity: identical graph to EMT
  gasGeneratorPF->connect({busGasPF});
  trGasPF->connect({busGasPF, busBPF});
  trLoadPF->connect({busLoadPF, busBPF});
  load2aPF->connect({busLoadPF});
  breakerPF->connect({busLoadPF, busLoad2bPF});
  load2bPF->connect({busLoad2bPF});

  SystemTopology systemPF(
      p.frequency, SystemNodeList{busGasPF, busBPF, busLoadPF, busLoad2bPF},
      SystemComponentList{gasGeneratorPF, trGasPF, trLoadPF, load2aPF,
                          breakerPF, load2bPF});

  auto loggerPF = DataLogger::make(pfSimName);
  loggerPF->logAttribute("BUS_gas_v", busGasPF->attribute("v"));
  loggerPF->logAttribute("BUS_b_v", busBPF->attribute("v"));
  loggerPF->logAttribute("BUS_load_2_v", busLoadPF->attribute("v"));
  loggerPF->logAttribute("BUS_load_2b_v", busLoad2bPF->attribute("v"));

  Simulation simulationPF(pfSimName, Logger::Level::debug);
  simulationPF.setSystem(systemPF);
  simulationPF.setTimeStep(p.pfTimeStep);
  simulationPF.setFinalTime(p.pfFinalTime);
  simulationPF.setDomain(Domain::SP);
  simulationPF.setSolverType(Solver::Type::NRP);
  simulationPF.setSolverAndComponentBehaviour(
      Solver::Behaviour::Initialization);
  simulationPF.doInitFromNodesAndTerminals(false);
  simulationPF.addLogger(loggerPF);
  simulationPF.run();

  SPDLOG_INFO("PF operating point:"
              "\n  V_GEN        = {} V RMS, angle = {} deg"
              "\n  V_BUS_b      = {} V RMS, angle = {} deg"
              "\n  V_BUS_load2  = {} V RMS, angle = {} deg"
              "\n  V_BUS_load2b = {} V RMS, angle = {} deg"
              "\n  GEN_gas      = P={} MW, Q={} Mvar",
              Math::abs(busGasPF->singleVoltage()),
              Math::phase(busGasPF->singleVoltage()) * 180.0 / PI,
              Math::abs(busBPF->singleVoltage()),
              Math::phase(busBPF->singleVoltage()) * 180.0 / PI,
              Math::abs(busLoadPF->singleVoltage()),
              Math::phase(busLoadPF->singleVoltage()) * 180.0 / PI,
              Math::abs(busLoad2bPF->singleVoltage()),
              Math::phase(busLoad2bPF->singleVoltage()) * 180.0 / PI,
              gasGeneratorPF->getApparentPower().real() / 1e6,
              gasGeneratorPF->getApparentPower().imag() / 1e6);

  // ==========================================================================
  // 2. EMT
  // ==========================================================================

  const String emtSimName = p.name + "_EMT";
  Logger::setLogDir("logs/" + emtSimName);

  // EMT nodes. No manual PF-result phasors are passed here.
  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busLoad = SimNode<Real>::make("BUS_load_2", PhaseType::ABC);
  auto busLoad2b = SimNode<Real>::make("BUS_load_2b", PhaseType::ABC);

  // EMT synchronous generator
  auto gasGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_gas", Logger::Level::off);

  gasGenerator->setBaseAndOperationalPerUnitParameters(
      p.genRatedPower, p.genRatedVoltage, p.genNominalFrequency,
      p.genPoleNumber, p.genNominalFieldCurrent, p.genStatorResistance, p.genLd,
      p.genLq, p.genLdTransient, p.genLqTransient, p.genLdSubtransient,
      p.genLqSubtransient, p.genLeakageInductance, p.genTd0Transient,
      p.genTq0Transient, p.genTd0Subtransient, p.genTq0Subtransient,
      p.genInertia);

  // Keep the existing legacy Signal::TurbineGovernor used by this project,
  // but derive its equilibrium Pm/Tm automatically from the PF operating point.
  gasGenerator->addGovernor(
      p.govTa, p.govTb, p.govTc, p.govF1a, p.govFa, p.govFb, p.govFc, p.govGain,
      p.govSpeedRelayTimeConstant, p.govServoMotorTimeConstant,
      Base::SynchronGenerator::LegacyGovernorInitialization::FromPowerflow);

  auto exciterParameters = std::make_shared<Signal::ExciterParameters>();
  exciterParameters->Ta = p.excTa;
  exciterParameters->Ka = p.excKa;
  exciterParameters->Te = p.excTe;
  exciterParameters->Ke = p.excKe;
  exciterParameters->Tf = p.excTf;
  exciterParameters->Kf = p.excKf;
  exciterParameters->Tr = p.excTr;
  exciterParameters->maxVr = p.excMaximumRegulatorVoltage;
  exciterParameters->minVr = p.excMinimumRegulatorVoltage;

  auto exciter = Signal::Exciter::make("GEN_gas_exciter", Logger::Level::off);
  gasGenerator->addExciter(exciter, exciterParameters);

  // EMT TR_gas
  auto trGas =
      EMT::Ph3::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  trGas->setParameters(
      p.trGasVoltageLow, p.trGasVoltageHigh, p.trGasRatedPower,
      p.trGasRatioMagnitude, p.trGasRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trGasResistance),
      Math::singlePhaseParameterToThreePhase(p.trGasInductance));

  // EMT TR_ld2
  auto trLoad =
      EMT::Ph3::Transformer::make("TR_ld2", "TR_ld2", Logger::Level::off, true);
  trLoad->setParameters(
      p.trLoadVoltageLow, p.trLoadVoltageHigh, p.trLoadRatedPower,
      p.trLoadRatioMagnitude, p.trLoadRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trLoadResistance),
      Math::singlePhaseParameterToThreePhase(p.trLoadInductance));

  // EMT loads intentionally have NO setParameters() call.
  // initWithPowerflow() transfers their PF P/Q to the matching EMT terminals.
  // RXLoad then derives R/X from that P/Q and the actual solved PF bus voltage
  // in doInitFromNodesAndTerminals(true).
  auto load2a = EMT::Ph3::RXLoad::make("LOAD_2a", Logger::Level::off);
  auto load2b = EMT::Ph3::RXLoad::make("LOAD_2b", Logger::Level::off);

  // EMT breaker, same topological position and same open/closed resistance
  auto breaker = EMT::Ph3::Switch::make("BR_load_2b", Logger::Level::off);
  breaker->setParameters(
      Math::singlePhaseParameterToThreePhase(p.breakerOpenResistance),
      Math::singlePhaseParameterToThreePhase(p.breakerClosedResistance));
  breaker->closeSwitch();

  // EMT connectivity: exactly the same graph as PF
  gasGenerator->connect({busGas});
  trGas->connect({busGas, busB});
  trLoad->connect({busLoad, busB});
  load2a->connect({busLoad});
  breaker->connect({busLoad, busLoad2b});
  load2b->connect({busLoad2b});

  SystemTopology systemEMT(p.frequency,
                           SystemNodeList{busGas, busB, busLoad, busLoad2b},
                           SystemComponentList{gasGenerator, trGas, trLoad,
                                               load2a, breaker, load2b});

  // Central PF -> EMT operating-point transfer:
  //   - solved node voltages by exact node name
  //   - synchronous-generator terminal P/Q with generator sign conversion
  //   - load terminal P/Q with consumer-positive sign convention
  //
  // No example-side power extraction or terminal assignment is required.
  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  auto loggerEMT = DataLogger::make(emtSimName, true, p.logDownSampling);

  loggerEMT->logAttribute("BUS_gas_v", busGas->attribute("v"));
  loggerEMT->logAttribute("BUS_b_v", busB->attribute("v"));
  loggerEMT->logAttribute("BUS_load_2_v", busLoad->attribute("v"));
  loggerEMT->logAttribute("BUS_load_2b_v", busLoad2b->attribute("v"));

  loggerEMT->logAttribute("GEN_gas_v", gasGenerator->attribute("v_intf"));
  loggerEMT->logAttribute("GEN_gas_i", gasGenerator->attribute("i_intf"));
  loggerEMT->logAttribute("GEN_gas_w_r", gasGenerator->attribute("w_r"));
  loggerEMT->logAttribute("GEN_gas_delta_r",
                          gasGenerator->attribute("delta_r"));
  loggerEMT->logAttribute("GEN_gas_P_elec", gasGenerator->attribute("P_elec"));
  loggerEMT->logAttribute("GEN_gas_Q_elec", gasGenerator->attribute("Q_elec"));
  loggerEMT->logAttribute("GEN_gas_P_mech", gasGenerator->attribute("P_mech"));

  loggerEMT->logAttribute("BR_load_2b_v", breaker->attribute("v_intf"));
  loggerEMT->logAttribute("BR_load_2b_i", breaker->attribute("i_intf"));
  loggerEMT->logAttribute("LOAD_2a_i", load2a->attribute("i_intf"));
  loggerEMT->logAttribute("LOAD_2b_i", load2b->attribute("i_intf"));

  Simulation simulationEMT(emtSimName, Logger::Level::info);
  simulationEMT.setSystem(systemEMT);
  simulationEMT.setTimeStep(p.timeStep);
  simulationEMT.setFinalTime(p.finalTime);
  simulationEMT.setDomain(Domain::EMT);

  // Dynamic models derive all internal initial states from the operating point
  // written by initWithPowerflow(). SynchronGeneratorVBR also initializes the
  // attached legacy TurbineGovernor from its PF-derived mechanical equilibrium.
  simulationEMT.doInitFromNodesAndTerminals(true);
  simulationEMT.doSystemMatrixRecomputation(p.recomputeSystemMatrix);
  simulationEMT.addLogger(loggerEMT);

  // false = open three-phase breaker -> disconnect LOAD_2b
  simulationEMT.addEvent(
      SwitchEvent3Ph::make(p.breakerOpenTime, breaker, false));

  SPDLOG_INFO("\n  Starting EMT load-loss simulation:"
              "\n  dt={} s"
              "\n  finalTime={} s"
              "\n  BR_load_2b opens at t={} s"
              "\n  disconnected load: P={} MW, Q={} Mvar",
              p.timeStep, p.finalTime, p.breakerOpenTime,
              p.load2bActivePower / 1e6, p.load2bReactivePower / 1e6);

  simulationEMT.run();

  return 0;
}
