// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_Switch.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorVBR.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/SP/SP_Ph1_PiLine.h>
#include <dpsim-models/SP/SP_Ph1_Switch.h>
#include <dpsim-models/SP/SP_Ph1_SynchronGenerator.h>
#include <dpsim-models/SP/SP_Ph1_Transformer.h>
#include <dpsim-models/Signal/Exciter.h>

#include <memory>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

// =============================================================================
// Parameters
// =============================================================================
enum class BenchmarkMode {
  Ideal,
  CurrentZero,
  ExponentialZCSEmulation,
};

EMT::Ph3::Switch::SwitchingMode toEMTMode(BenchmarkMode mode) {
  switch (mode) {
  case BenchmarkMode::Ideal:
    return EMT::Ph3::Switch::SwitchingMode::Ideal;

  case BenchmarkMode::CurrentZero:
    return EMT::Ph3::Switch::SwitchingMode::CurrentZero;

  case BenchmarkMode::ExponentialZCSEmulation:
    return EMT::Ph3::Switch::SwitchingMode::ExponentialZCSEmulation;
  }

  return EMT::Ph3::Switch::SwitchingMode::Ideal;
}
struct Parameters {
  // Simulation
  String name = "EMT_Inrush_undervoltage_SynGenVBR";
  Real frequency = 50.0;
  Real pfTimeStep = 1.0;
  Real pfFinalTime = 1.0;
  Real timeStep = 10e-6;
  Real finalTime = 10.0;
  Real breakerCloseTime = 2.0;
  Bool recomputeSystemMatrix = true;
  Int logDownSampling = 1;

  // ---------------------------------------------------------------------------
  // Synchronous generator GEN_gas
  // ---------------------------------------------------------------------------
  Real gasRatedPower = 50e6;
  Real gasRatedVoltage = 10.5e3;
  Real gasNominalFrequency = 50.0;
  Int gasPoleNumber = 2;
  Real gasNominalFieldCurrent = 1300.0;

  Real gasStatorResistance = 0.002;
  Real gasLd = 2.4;
  Real gasLq = 1.33;
  Real gasLdTransient = 0.31;
  Real gasLqTransient = 1.2;
  Real gasLdSubtransient = 0.24;
  Real gasLqSubtransient = 0.35;
  Real gasLeakageInductance = 0.135;
  Real gasTd0Transient = 1.45;
  Real gasTq0Transient = 1.0;
  Real gasTd0Subtransient = 0.022;
  Real gasTq0Subtransient = 0.0095;
  Real gasInertia = 5.0;

  // GEN_gas turbine governor
  Real gasGovTa = 1.0;
  Real gasGovTb = 0.5;
  Real gasGovTc = 0.2;
  Real gasGovF1a = 0.3;
  Real gasGovFa = 0.3;
  Real gasGovFb = 0.25;
  Real gasGovFc = 0.3;
  Real gasGovGain = 20.0;
  Real gasGovSpeedRelayTimeConstant = 0.1;
  Real gasGovServoMotorTimeConstant = 0.1;

  // GEN_gas exciter
  Real gasExcTa = 0.005;
  Real gasExcKa = 200.0;
  Real gasExcTe = 0.05;
  Real gasExcKe = 0.5;
  Real gasExcTf = 0.3;
  Real gasExcKf = 0.01;
  Real gasExcTr = 0.02;
  Real gasExcMaximumRegulatorVoltage = 9.0;
  Real gasExcMinimumRegulatorVoltage = -5.0;

  // ---------------------------------------------------------------------------
  // Synchronous generator GEN_psh
  // ---------------------------------------------------------------------------
  Real pshRatedPower = 200e6;
  Real pshRatedVoltage = 18.0e3;
  Real pshNominalFrequency = 50.0;
  Int pshPoleNumber = 2;
  Real pshNominalFieldCurrent = 1300.0;

  Real pshStatorResistance = 0.002;
  Real pshLd = 2.0;
  Real pshLq = 2.0;
  Real pshLdTransient = 0.3;
  Real pshLqTransient = 0.3;
  Real pshLdSubtransient = 0.2;
  Real pshLqSubtransient = 0.2;
  Real pshLeakageInductance = 0.1;
  Real pshTd0Transient = 1.0;
  Real pshTq0Transient = 1.0;
  Real pshTd0Subtransient = 0.05;
  Real pshTq0Subtransient = 0.05;
  Real pshInertia = 4.0;

  // GEN_psh turbine governor
  Real pshGovTa = 0.3;
  Real pshGovTb = 7.0;
  Real pshGovTc = 0.5;
  Real pshGovF1a = 0.3;
  Real pshGovFa = 0.25;
  Real pshGovFb = 0.3;
  Real pshGovFc = 0.15;
  Real pshGovGain = 25.0;
  Real pshGovSpeedRelayTimeConstant = 0.1;
  Real pshGovServoMotorTimeConstant = 0.1;

  // GEN_psh exciter
  Real pshExcTa = 0.005;
  Real pshExcKa = 200.0;
  Real pshExcTe = 0.05;
  Real pshExcKe = 0.5;
  Real pshExcTf = 0.3;
  Real pshExcKf = 0.01;
  Real pshExcTr = 0.02;
  Real pshExcMaximumRegulatorVoltage = 9.0;
  Real pshExcMinimumRegulatorVoltage = -5.0;

  // ---------------------------------------------------------------------------
  // Transformer TR_gas
  // ---------------------------------------------------------------------------
  Real trGasVoltageLow = 10.5e3;
  Real trGasVoltageHigh = 220e3;
  Real trGasRatedPower = 50e6;
  Real trGasRatioMagnitude = 10.5e3 / 220e3;
  Real trGasRatioPhase = 0.0;
  Real trGasResistance = 3.64157728;
  Real trGasInductance = 0.31037265855769886;

  // ---------------------------------------------------------------------------
  // Line
  // ---------------------------------------------------------------------------
  Real lineBaseVoltage = 220e3;
  Real lineResistance = 5.04669;
  Real lineInductance = 0.13523123641;
  Real lineCapacitance = 1.93522865e-6;
  Real lineConductance = 1e-15;

  // ---------------------------------------------------------------------------
  // Line 3
  // ---------------------------------------------------------------------------
  Real line3BaseVoltage = 220e3;
  Real line3Length = 5.0;
  Real line3ResistancePerLength = 0.0749;
  Real line3InductancePerLength = 1.270693e-3;
  Real line3CapacitancePerLength = 0.00466961e-6;
  Real line3Conductance = 1e-15;

  // ---------------------------------------------------------------------------
  // Breaker GEN_psh
  // ---------------------------------------------------------------------------
  Real breakerOpenResistance = 1e12;
  Real breakerClosedResistance = 1e-3;

  // ---------------------------------------------------------------------------
  // Transformer TR_psh
  // ---------------------------------------------------------------------------
  Real trPshVoltageLow = 18.0e3;
  Real trPshVoltageHigh = 220e3;
  Real trPshRatedPower = 200e6;
  Real trPshRatioMagnitude = 18.0e3 / 220e3;
  Real trPshRatioPhase = 0.0;
  Real trPshResistance = 0.8349;
  Real trPshInductance = 0.0962521551189048;

  // Transformer TR_load
  Real trLoadVoltageLow = 10.0e3;
  Real trLoadVoltageHigh = 220e3;
  Real trLoadRatedPower = 50e6;
  Real trLoadRatioMagnitude = 10.0e3 / 220e3;
  Real trLoadRatioPhase = 0.0;
  Real trLoadResistance = 3.64157728;
  Real trLoadInductance = 0.31037292;

  // Transformer TR_1
  Real tr1VoltageLow = 220e3;
  Real tr1VoltageHigh = 400e3;
  Real tr1RatedPower = 300e6;
  Real tr1RatioMagnitude = 220e3 / 400e3;
  Real tr1RatioPhase = 0.0;
  Real tr1Resistance = 2.6666666;
  Real tr1Inductance = 0.2009929154769824;

  

  // LOAD_2a
  Real load2aActivePower = 5e6;
  Real load2aReactivePower = 3.098721e6;
  Real load2aNominalVoltage = 10e3;

  // LOAD_2b
  Real load2bActivePower = 15e6;
  Real load2bReactivePower = 0.0;
  Real load2bNominalVoltage = 10e3;
};

int main(int argc, char *argv[]) {
  Parameters p;

  // Keep the same command-line control style as the previous Scenario A
  // example for dt, duration, simulation name, and breaker closing time.
  CommandLineArgs args(argc, argv);

  if (argc > 1) {
    p.timeStep = args.timeStep;
    p.finalTime = args.duration;

    if (args.name != "dpsim")
      p.name = args.name;

    if (args.options.find("BREAKER_CLOSE_TIME") != args.options.end())
      p.breakerCloseTime = args.getOptionReal("BREAKER_CLOSE_TIME");
  }

  if (!(p.timeStep > 0.0) || !(p.finalTime > 0.0) ||
      !(p.breakerCloseTime >= 0.0) || !(p.breakerCloseTime < p.finalTime)) {
    throw std::invalid_argument(
        "Require dt>0, finalTime>0, and 0<=breakerCloseTime<finalTime.");
  }

  const Real line3Resistance = p.line3ResistancePerLength * p.line3Length;
  const Real line3Inductance = p.line3InductancePerLength * p.line3Length;
  const Real line3Capacitance = p.line3CapacitancePerLength * p.line3Length;

  SPDLOG_INFO("Scenario A bases:"
              "\n  GEN_gas: S_rated={} VA, V_rated={} V_LL RMS, H={} s"
              "\n  GEN_psh: S_rated={} VA, V_rated={} V_LL RMS, H={} s"
              "\n  TR_gas : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
              "\n  TR_psh : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
              "\n  breaker: R_open={} Ohm, R_closed={} Ohm",
              p.gasRatedPower, p.gasRatedVoltage, p.gasInertia, p.pshRatedPower,
              p.pshRatedVoltage, p.pshInertia, p.trGasVoltageLow,
              p.trGasVoltageHigh, p.trGasRatedPower, p.trGasResistance,
              p.trGasInductance, p.trPshVoltageLow, p.trPshVoltageHigh,
              p.trPshRatedPower, p.trPshResistance, p.trPshInductance,
              p.breakerOpenResistance, p.breakerClosedResistance);

  // ==========================================================================
  // 1. POWER FLOW
  //
  // PF and EMT have exactly the same graph:
  //
  // GEN_gas -- BUS_gas -- TR_gas -- BUS_b -- cable -- BUS_a
  //                                                      |
  //                                                   line_3
  //                                                      |
  //                                                 BUS_psha
  //                                                      |
  //                                             breaker_GEN_psh
  //                                                      |
  //                                              BUS_tr_psh
  //                                                      |
  //                                                   TR_psh
  //                                                      |
  //                                                 BUS_psh
  //                                                      |
  //                                                  GEN_psh
  //
  // The breaker is OPEN in both PF and EMT initialization. Each side therefore
  // has its own VD generator. PF and EMT use identical component and node names,
  // allowing SystemTopology::initWithPowerflow() to transfer the complete
  // operating point without example-side terminal-power assignments.
  // ==========================================================================

  const String pfSimName = p.name + "_PF";
  Logger::setLogDir("logs/" + pfSimName);

  // PF nodes
  auto BUS_gas_PF = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto BUS_b_PF = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto BUS_a_PF = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto BUS_psha_PF = SimNode<Complex>::make("BUS_psha", PhaseType::Single);
  auto BUS_psh_PF = SimNode<Complex>::make("BUS_psh", PhaseType::Single); 
  auto BUS_load_PF = SimNode<Complex>::make("BUS_load", PhaseType::Single);

  auto BUS_tr_PF = SimNode<Complex>::make("BUS_tr", PhaseType::Single);
  auto BUS_c_PF = SimNode<Complex>::make("BUS_c", PhaseType::Single);
  

  // PF GEN_gas: VD/slack generator for the gas-side island.
  auto GEN_gas_PF =
      SP::Ph1::SynchronGenerator::make("GEN_gas", Logger::Level::off);
  GEN_gas_PF->setParameters(
                            p.gasRatedPower, p.gasRatedVoltage,
                            0.0, // initial P estimate; VD P/Q are solved by NRP
                            p.gasRatedVoltage, PowerflowBusType::VD);
  GEN_gas_PF->setBaseVoltage(p.gasRatedVoltage);

   
  // PF GEN_psh: VD/slack generator for the isolated PSH-side island.
  auto GEN_psh_PF = SP::Ph1::SynchronGenerator::make("GEN_psh", Logger::Level::off);
  GEN_psh_PF->setParameters(
      p.pshRatedPower, p.pshRatedVoltage,
      0.0, // initial P estimate; VD P/Q are solved by NRP
      p.pshRatedVoltage, PowerflowBusType::PV);
  GEN_psh_PF->setBaseVoltage(p.pshRatedVoltage);


  // PF TR_gas
  auto TR_gas_PF =
      SP::Ph1::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  TR_gas_PF->setParameters(p.trGasVoltageLow, p.trGasVoltageHigh,
                         p.trGasRatedPower, p.trGasRatioMagnitude,
                         p.trGasRatioPhase, p.trGasResistance,
                         p.trGasInductance);
  TR_gas_PF->setBaseVoltage(p.trGasVoltageHigh);

  // PF TR_psh
  auto TR_psh_PF =
      SP::Ph1::Transformer::make("TR_psh", "TR_psh", Logger::Level::off, true);
  TR_psh_PF->setParameters(p.trPshVoltageLow, p.trPshVoltageHigh,
                         p.trPshRatedPower, p.trPshRatioMagnitude,
                         p.trPshRatioPhase, p.trPshResistance,
                         p.trPshInductance);
  TR_psh_PF->setBaseVoltage(p.trPshVoltageHigh);

  // PF TR_ld2
  auto TR_load_PF =
      SP::Ph1::Transformer::make("TR_load", "TR_load", Logger::Level::off, true);
  TR_load_PF->setParameters(p.trLoadVoltageLow, p.trLoadVoltageHigh,
                          p.trLoadRatedPower, p.trLoadRatioMagnitude,
                          p.trLoadRatioPhase, p.trLoadResistance,
                          p.trLoadInductance);
  TR_load_PF->setBaseVoltage(p.trLoadVoltageHigh);

  // PF line 1
  auto LINE_1_PF = SP::Ph1::PiLine::make("LINE_1", Logger::Level::off);
  LINE_1_PF->setParameters(p.lineResistance, p.lineInductance, p.lineCapacitance, p.lineConductance);
  LINE_1_PF->setBaseVoltage(p.lineBaseVoltage);

  // PF line 2
  auto LINE_2_PF = SP::Ph1::PiLine::make("LINE_2", Logger::Level::off);
  LINE_2_PF->setParameters(p.lineResistance, p.lineInductance, p.lineCapacitance, p.lineConductance);
  LINE_2_PF->setBaseVoltage(p.lineBaseVoltage);

  // PF line_3
  auto LINE_3_PF = SP::Ph1::PiLine::make("LINE_3", Logger::Level::off);
  LINE_3_PF->setParameters(line3Resistance, line3Inductance, line3Capacitance, p.line3Conductance);
  LINE_3_PF->setBaseVoltage(p.line3BaseVoltage);

  // PF LOAD_2a
  auto LOAD_a_PF = SP::Ph1::Load::make("LOAD_a", Logger::Level::off);
  LOAD_a_PF->setParameters(p.load2aActivePower, p.load2aReactivePower,
                          p.load2aNominalVoltage);
  LOAD_a_PF->modifyPowerFlowBusType(PowerflowBusType::PQ);


  // PF LOAD_2b behind the breaker
  auto LOAD_b_PF = SP::Ph1::Load::make("LOAD_b", Logger::Level::off);
  LOAD_b_PF->setParameters(p.load2bActivePower, p.load2bReactivePower,
                          p.load2bNominalVoltage);
  LOAD_b_PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // PF breaker: actual SP switch, initially open.
  auto BR_tr_PF = SP::Ph1::Switch::make("BR_tr", Logger::Level::off);
  BR_tr_PF->setParameters(p.breakerOpenResistance, p.breakerClosedResistance, false);
  BR_tr_PF->setBaseVoltage(p.trPshVoltageHigh);

  // PF TR_ld2
  auto TR_1_PF =
      SP::Ph1::Transformer::make("TR_1", "TR_1", Logger::Level::off, true);
  TR_1_PF->setParameters(p.tr1VoltageLow, p.tr1VoltageHigh,
                          p.tr1RatedPower, p.tr1RatioMagnitude,
                          p.tr1RatioPhase, p.trLoadResistance,
                          p.tr1Inductance);
  TR_1_PF->setBaseVoltage(p.tr1VoltageHigh);

  
  // PF connectivity: identical graph to EMT.
  GEN_gas_PF->connect({BUS_gas_PF});
  TR_gas_PF->connect({BUS_gas_PF, BUS_b_PF});
  LINE_1_PF->connect({BUS_b_PF, BUS_a_PF});
  LINE_2_PF->connect({BUS_b_PF, BUS_a_PF});

  LINE_3_PF->connect({BUS_a_PF, BUS_psha_PF});
  TR_psh_PF->connect({BUS_psh_PF, BUS_psha_PF});
  GEN_psh_PF->connect({BUS_psh_PF});
  
  TR_load_PF->connect({BUS_load_PF, BUS_b_PF});
  LOAD_a_PF->connect({BUS_load_PF});
  LOAD_b_PF->connect({BUS_load_PF});

  BR_tr_PF->connect({BUS_b_PF, BUS_tr_PF});
  TR_1_PF->connect({BUS_tr_PF, BUS_c_PF});

  SystemTopology systemPF(
      p.frequency,
      SystemNodeList{BUS_gas_PF, BUS_b_PF, BUS_a_PF, BUS_psha_PF, BUS_psh_PF, BUS_load_PF, BUS_tr_PF, BUS_c_PF},
      SystemComponentList{GEN_gas_PF, TR_gas_PF, LINE_1_PF, LINE_2_PF, LINE_3_PF,
                          TR_psh_PF, GEN_psh_PF, TR_load_PF, LOAD_a_PF, LOAD_b_PF, BR_tr_PF, TR_1_PF});

  auto loggerPF = DataLogger::make(pfSimName);
  loggerPF->logAttribute("BUS_gas_v", BUS_gas_PF->attribute("v"));
  loggerPF->logAttribute("BUS_b_v", BUS_b_PF->attribute("v"));
  loggerPF->logAttribute("BUS_a_v", BUS_a_PF->attribute("v"));
  loggerPF->logAttribute("BUS_psha_v", BUS_psha_PF->attribute("v"));
  loggerPF->logAttribute("BUS_tr_psh_v", BUS_psh_PF->attribute("v"));
  loggerPF->logAttribute("BUS_psh_v", BUS_load_PF->attribute("v"));
  loggerPF->logAttribute("BUS_psh_v", BUS_tr_PF->attribute("v"));
  loggerPF->logAttribute("BUS_psh_v", BUS_c_PF->attribute("v"));
  loggerPF->logAttribute("breaker_GEN_psh_i", BR_tr_PF->attribute("i_intf"));

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

  // ==========================================================================
  // 2. EMT
  // ==========================================================================

  const String emtSimName = p.name + "_EMT";
  Logger::setLogDir("logs/" + emtSimName);

  // EMT nodes. No manual PF-result phasors are passed here.
  auto BUS_gas_EMT = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto BUS_b_EMT = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto BUS_a_EMT = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto BUS_psha_EMT = SimNode<Real>::make("BUS_psha", PhaseType::ABC);
  auto BUS_psh_EMT = SimNode<Real>::make("BUS_psh", PhaseType::ABC);
  auto BUS_load_EMT = SimNode<Real>::make("BUS_load", PhaseType::ABC);

  auto BUS_tr_EMT = SimNode<Real>::make("BUS_tr", PhaseType::ABC);
  auto BUS_c_EMT = SimNode<Real>::make("BUS_c", PhaseType::ABC);

  
  // ---------------------------------------------------------------------------
  // EMT GEN_gas
  // ---------------------------------------------------------------------------
  auto GEN_gas_EMT =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_gas", Logger::Level::off);

  GEN_gas_EMT->setBaseAndOperationalPerUnitParameters(
      p.gasRatedPower, p.gasRatedVoltage, p.gasNominalFrequency,
      p.gasPoleNumber, p.gasNominalFieldCurrent, p.gasStatorResistance, p.gasLd,
      p.gasLq, p.gasLdTransient, p.gasLqTransient, p.gasLdSubtransient,
      p.gasLqSubtransient, p.gasLeakageInductance, p.gasTd0Transient,
      p.gasTq0Transient, p.gasTd0Subtransient, p.gasTq0Subtransient,
      p.gasInertia);

  // Initialize the legacy governor directly from the PF operating point.
  GEN_gas_EMT->addGovernor(
      p.gasGovTa, p.gasGovTb, p.gasGovTc, p.gasGovF1a, p.gasGovFa, p.gasGovFb,
      p.gasGovFc, p.gasGovGain, p.gasGovSpeedRelayTimeConstant,
      p.gasGovServoMotorTimeConstant,
      Base::SynchronGenerator::LegacyGovernorInitialization::FromPowerflow);

  auto gasExciterParameters = std::make_shared<Signal::ExciterParameters>();
  gasExciterParameters->Ta = p.gasExcTa;
  gasExciterParameters->Ka = p.gasExcKa;
  gasExciterParameters->Te = p.gasExcTe;
  gasExciterParameters->Ke = p.gasExcKe;
  gasExciterParameters->Tf = p.gasExcTf;
  gasExciterParameters->Kf = p.gasExcKf;
  gasExciterParameters->Tr = p.gasExcTr;
  gasExciterParameters->maxVr = p.gasExcMaximumRegulatorVoltage;
  gasExciterParameters->minVr = p.gasExcMinimumRegulatorVoltage;

  auto gasExciter =
      Signal::Exciter::make("GEN_gas_exciter", Logger::Level::off);
  GEN_gas_EMT->addExciter(gasExciter, gasExciterParameters);

  // ---------------------------------------------------------------------------
  // EMT GEN_psh
  // ---------------------------------------------------------------------------
  auto GEN_psh_EMT =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_psh", Logger::Level::off);

  GEN_psh_EMT->setBaseAndOperationalPerUnitParameters(
      p.pshRatedPower, p.pshRatedVoltage, p.pshNominalFrequency,
      p.pshPoleNumber, p.pshNominalFieldCurrent, p.pshStatorResistance, p.pshLd,
      p.pshLq, p.pshLdTransient, p.pshLqTransient, p.pshLdSubtransient,
      p.pshLqSubtransient, p.pshLeakageInductance, p.pshTd0Transient,
      p.pshTq0Transient, p.pshTd0Subtransient, p.pshTq0Subtransient,
      p.pshInertia);

  GEN_psh_EMT->addGovernor(
      p.pshGovTa, p.pshGovTb, p.pshGovTc, p.pshGovF1a, p.pshGovFa, p.pshGovFb,
      p.pshGovFc, p.pshGovGain, p.pshGovSpeedRelayTimeConstant,
      p.pshGovServoMotorTimeConstant,
      Base::SynchronGenerator::LegacyGovernorInitialization::FromPowerflow);

  auto pshExciterParameters = std::make_shared<Signal::ExciterParameters>();
  pshExciterParameters->Ta = p.pshExcTa;
  pshExciterParameters->Ka = p.pshExcKa;
  pshExciterParameters->Te = p.pshExcTe;
  pshExciterParameters->Ke = p.pshExcKe;
  pshExciterParameters->Tf = p.pshExcTf;
  pshExciterParameters->Kf = p.pshExcKf;
  pshExciterParameters->Tr = p.pshExcTr;
  pshExciterParameters->maxVr = p.pshExcMaximumRegulatorVoltage;
  pshExciterParameters->minVr = p.pshExcMinimumRegulatorVoltage;

  auto pshExciter =
      Signal::Exciter::make("GEN_psh_exciter", Logger::Level::off);
  GEN_psh_EMT->addExciter(pshExciter, pshExciterParameters);


  // EMT TR_gas
  auto TR_gas_EMT =
      EMT::Ph3::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  TR_gas_EMT->setParameters(
      p.trGasVoltageLow, p.trGasVoltageHigh, p.trGasRatedPower,
      p.trGasRatioMagnitude, p.trGasRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trGasResistance),
      Math::singlePhaseParameterToThreePhase(p.trGasInductance));

  // EMT TR_psh
  auto TR_psh_EMT =
      EMT::Ph3::Transformer::make("TR_psh", "TR_psh", Logger::Level::off, true);
  TR_psh_EMT->setParameters(
      p.trPshVoltageLow, p.trPshVoltageHigh, p.trPshRatedPower,
      p.trPshRatioMagnitude, p.trPshRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trPshResistance),
      Math::singlePhaseParameterToThreePhase(p.trPshInductance));

  // EMT TR_1
  auto TR_1_EMT =
      EMT::Ph3::Transformer::make("TR_1", "TR_1", Logger::Level::off, true);
  TR_1_EMT->setParameters(
      p.tr1VoltageLow, p.tr1VoltageHigh, p.tr1RatedPower,
      p.tr1RatioMagnitude, p.tr1RatioPhase,
      Math::singlePhaseParameterToThreePhase(p.tr1Resistance),
      Math::singlePhaseParameterToThreePhase(p.tr1Inductance));
  
  // EMT TR_ld2
  auto TR_load_EMT =
      EMT::Ph3::Transformer::make("TR_load", "TR_load", Logger::Level::off, true);
  TR_load_EMT->setParameters(
      p.trLoadVoltageLow, p.trLoadVoltageHigh, p.trLoadRatedPower,
      p.trLoadRatioMagnitude, p.trLoadRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trLoadResistance),
      Math::singlePhaseParameterToThreePhase(p.trLoadInductance));

  // EMT line 1
  auto LINE_1_EMT = EMT::Ph3::PiLine::make("LINE_1", Logger::Level::off);
  LINE_1_EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(p.lineResistance),
      Math::singlePhaseParameterToThreePhase(p.lineInductance),
      Math::singlePhaseParameterToThreePhase(p.lineCapacitance),
      Math::singlePhaseParameterToThreePhase(p.lineConductance));

  // EMT line 2
  auto LINE_2_EMT = EMT::Ph3::PiLine::make("LINE_2", Logger::Level::off);
  LINE_2_EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(p.lineResistance),
      Math::singlePhaseParameterToThreePhase(p.lineInductance),
      Math::singlePhaseParameterToThreePhase(p.lineCapacitance),
      Math::singlePhaseParameterToThreePhase(p.lineConductance));

  // EMT line_3
  auto LINE_3_EMT = EMT::Ph3::PiLine::make("LINE_3", Logger::Level::off);
  LINE_3_EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(line3Resistance),
      Math::singlePhaseParameterToThreePhase(line3Inductance),
      Math::singlePhaseParameterToThreePhase(line3Capacitance),
      Math::singlePhaseParameterToThreePhase(p.line3Conductance));


  // EMT breaker: same topological position and same initial OPEN state as PF.
  auto BR_tr_EMT = EMT::Ph3::Switch::make("BR_tr", Logger::Level::off);
  BR_tr_EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(p.breakerOpenResistance),
      Math::singlePhaseParameterToThreePhase(p.breakerClosedResistance));
  BR_tr_EMT->openSwitch();


  auto LOAD_a_EMT = EMT::Ph3::RXLoad::make("LOAD_a", Logger::Level::off);
  auto LOAD_b_EMT = EMT::Ph3::RXLoad::make("LOAD_b", Logger::Level::off);

  // Piecewise-linear inductor
  ///////////////////////////////////////////////////////////////////////////////////
  std::vector<Real> mFluxBreakpoints;
  std::vector<Real> mCurrentBreakpoints;

  // Real SourceVoltagePeak = (220e3 * sqrt(2)) / sqrt(3);
  // const Real baseFlux = std::sqrt(2.0) * SourceVoltagePeak / (2.0 * PI * 50);
  // const std::vector<Real> fluxInPerUnit{0.0, 0.5, 0.9, 1.0};
  // const std::vector<Real> currentInPerUnit{0.0, 0.01, 0.1, 10.0};

  const std::vector<Real> fluxInPerUnit{0.0, 1.2, 1.2 + 0.2*(1 - 1e-3)};
  const std::vector<Real> currentInPerUnit{0.0, 1.1e-3, 1};
  
  Real voltage = p.tr1VoltageHigh;
  Real tr1BaseCurrent_onePhase_RMS =  p.tr1RatedPower / (voltage * sqrt(3.0));
  Real tr1BaseCurrent_onePhase_Peak = sqrt(2.0) * (tr1BaseCurrent_onePhase_RMS);
  Real tr1_flux_LG_Peak = (sqrt(2.0/3.0) * voltage) / (2 * PI *50);

    for (size_t k = 0; k < fluxInPerUnit.size(); ++k) {
      mFluxBreakpoints.push_back(fluxInPerUnit[k] * tr1_flux_LG_Peak);
      mCurrentBreakpoints.push_back(currentInPerUnit[k] * tr1BaseCurrent_onePhase_Peak);
    }
  

  auto lPiecewiseLinear = EMT::Ph3::PiecewiseLinearInductor::make(
        "LPiecewiseLinear", Logger::Level::debug);
    lPiecewiseLinear->setParameters(mFluxBreakpoints, mCurrentBreakpoints);
  //////////////////////////////////////////////////////////////////////////////////////////////

  // nonlinear inductor old
  //Base values taken from CoPilot for this circuit first inductor
	// Real V_LL_RMS = 400e3;
	// Real V_LG_Peak = sqrt(2.0/3.0) * V_LL_RMS; 
	// Real S_three_phase = 300e6;
	// Real I_one_phase_RMS = S_three_phase / (V_LL_RMS * sqrt(3.0));
	// Real I_one_phase_Peak = sqrt(2.0) * I_one_phase_RMS;    // base value for characteristics - current
	// Real flux_LG_Peak = V_LG_Peak / (2 * M_PI * 50);		// base value for characteristics - flux

	// // Real VB   = sqrt(2.0/3.0) * 220e3;
	// // Real SB   = 200e6;
	// // Real wB   = 2*M_PI*50;
	// //Real IB   = SB/(sqrt(Real(3))*VB);
	// // Real IB   = SB/VB;
	// // Real psiB = VB/wB;

	// std::vector<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::Point> Pts = {
	// { 0.0 * I_one_phase_Peak,    0.0 * flux_LG_Peak  },
	// { 1.1e-3 * I_one_phase_Peak, 1.2 * flux_LG_Peak },
	// };

	// // std::vector<CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic::Point> Pts = {
	// // { 0.0 * IB,    0.0 * psiB  },
	// // { 8.8e-4 * IB, 0.88 * psiB },
	// // { 0.04 * IB,   1.04 * psiB },
	// // { 0.2 *IB,    1.10 * psiB },
	// // };

	// auto characteristic = std::make_shared<
	// CPS::EMT::Ph1::PieceWiseNonLinearCharacteristic
	// >(Pts, 0.2 * flux_LG_Peak/I_one_phase_Peak);

	// auto nl_ind = CPS::EMT::Ph3::NonLinearInductor::make("nl_ind", Logger::Level::info);
	// nl_ind->setPieceWiseCharacteristic(characteristic);

	// // Set time step for the nonlinear inductor
	// nl_ind->setTimeStep(p.timeStep);
	

	// auto resistor_EMT = EMT::Ph3::Resistor::make("parallel_resistor", Logger::Level::info);
	// resistor_EMT->setParameters(Math::singlePhaseParameterToThreePhase(1e3));
  ///////////////////////////////////////////////////////////////////////////////////////


  // EMT connectivity: exactly the same graph and names as PF.
  GEN_gas_EMT->connect({BUS_gas_EMT});
  TR_gas_EMT->connect({BUS_gas_EMT, BUS_b_EMT});
  LINE_1_EMT->connect({BUS_b_EMT, BUS_a_EMT});
  LINE_2_EMT->connect({BUS_b_EMT, BUS_a_EMT});

  LINE_3_EMT->connect({BUS_a_EMT, BUS_psha_EMT});
  TR_psh_EMT->connect({BUS_psh_EMT, BUS_psha_EMT});
  GEN_psh_EMT->connect({BUS_psh_EMT});
  
  TR_load_EMT->connect({BUS_load_EMT, BUS_b_EMT});
  LOAD_a_EMT->connect({BUS_load_EMT});
  LOAD_b_EMT->connect({BUS_load_EMT});

  BR_tr_EMT->connect({BUS_b_EMT, BUS_tr_EMT});
  TR_1_EMT->connect({BUS_tr_EMT, BUS_c_EMT});

   
  lPiecewiseLinear->connect({BUS_c_EMT, EMT::SimNode::GND});
  // nl_ind->connect({BUS_c_EMT, EMT::SimNode::GND});
  // resistor_EMT->connect({ BUS_c_EMT, EMT::SimNode::GND });


  SystemTopology systemEMT(
      p.frequency,
      SystemNodeList{BUS_gas_EMT, BUS_b_EMT, BUS_a_EMT, BUS_psha_EMT, BUS_psh_EMT, BUS_load_EMT, BUS_tr_EMT, BUS_c_EMT},
      SystemComponentList{GEN_gas_EMT, TR_gas_EMT, LINE_1_EMT, LINE_2_EMT, LINE_3_EMT, TR_psh_EMT,
                          GEN_psh_EMT, TR_load_EMT, LOAD_a_EMT, LOAD_b_EMT, BR_tr_EMT, TR_1_EMT, lPiecewiseLinear});

  // Central PF -> EMT operating-point transfer:
  //   - solved node voltages by exact node name
  //   - synchronous-generator terminal P/Q with generator sign conversion
  //
  // No manual generatorPositivePower(), terminal()->setPower(), or manually
  // supplied initial phasors are required in the example.
  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  auto loggerEMT = DataLogger::make(emtSimName, true, p.logDownSampling);

  loggerEMT->logAttribute("BUS_gas_v", BUS_gas_EMT->attribute("v"));
  loggerEMT->logAttribute("BUS_b_v", BUS_b_EMT->attribute("v"));
  loggerEMT->logAttribute("BUS_a_v", BUS_a_EMT->attribute("v"));
  loggerEMT->logAttribute("BUS_psha_v", BUS_psha_EMT->attribute("v"));
  loggerEMT->logAttribute("BUS_c_v", BUS_c_EMT->attribute("v"));
  loggerEMT->logAttribute("BUS_psh_v", BUS_psh_EMT->attribute("v"));

  loggerEMT->logAttribute("TR_gas_i", TR_gas_EMT->attribute("i_intf"));
  loggerEMT->logAttribute("line_1_i", LINE_1_EMT->attribute("i_intf"));
  loggerEMT->logAttribute("line_3_i", LINE_3_EMT->attribute("i_intf"));
  loggerEMT->logAttribute("breaker_tr_i", BR_tr_EMT->attribute("i_intf"));
  loggerEMT->logAttribute("breaker_tr_v", BR_tr_EMT->attribute("v_intf"));
  loggerEMT->logAttribute("TR_psh_i", TR_psh_EMT->attribute("i_intf"));

  loggerEMT->logAttribute("GEN_gas_w_r", GEN_gas_EMT->attribute("w_r"));
  loggerEMT->logAttribute("GEN_gas_P_elec", GEN_gas_EMT->attribute("P_elec"));
  loggerEMT->logAttribute("GEN_gas_Q_elec", GEN_gas_EMT->attribute("Q_elec"));
  loggerEMT->logAttribute("GEN_gas_P_mech", GEN_gas_EMT->attribute("P_mech"));
  loggerEMT->logAttribute("GEN_gas_i_intf", GEN_gas_EMT->attribute("i_intf"));
  loggerEMT->logAttribute("GEN_gas_v_intf", GEN_gas_EMT->attribute("v_intf"));
  loggerEMT->logAttribute("GEN_gas_delta_r", GEN_gas_EMT->attribute("delta_r"));
  loggerEMT->logAttribute("GEN_gas_T_e", GEN_gas_EMT->attribute("T_e"));

  loggerEMT->logAttribute("GEN_psh_w_r", GEN_psh_EMT->attribute("w_r"));
  loggerEMT->logAttribute("GEN_psh_P_elec", GEN_psh_EMT->attribute("P_elec"));
  loggerEMT->logAttribute("GEN_psh_Q_elec", GEN_psh_EMT->attribute("Q_elec"));
  loggerEMT->logAttribute("GEN_psh_P_mech", GEN_psh_EMT->attribute("P_mech"));
  loggerEMT->logAttribute("GEN_psh_i_intf", GEN_psh_EMT->attribute("i_intf"));
  loggerEMT->logAttribute("GEN_psh_v_intf", GEN_psh_EMT->attribute("v_intf"));
  loggerEMT->logAttribute("GEN_psh_delta_r", GEN_psh_EMT->attribute("delta_r"));
  loggerEMT->logAttribute("GEN_psh_T_e", GEN_psh_EMT->attribute("T_e"));

  loggerEMT->logAttribute("voltage", BUS_tr_EMT->attribute("v"));
  loggerEMT->logAttribute("current", lPiecewiseLinear->attribute("i_intf"));
  loggerEMT->logAttribute("flux", lPiecewiseLinear->attribute("x"));

  // loggerEMT->logAttribute("current", nl_ind->attribute("i_intf"));
  // loggerEMT->logAttribute("flux", nl_ind->attribute("flux"));

  Simulation simulationEMT(emtSimName, Logger::Level::info);
  simulationEMT.setSystem(systemEMT);
  simulationEMT.setTimeStep(p.timeStep);
  simulationEMT.setFinalTime(p.finalTime);
  simulationEMT.setDomain(Domain::EMT);
  simulationEMT.setSolverType(Solver::Type::MNA);

  // Dynamic models derive all internal states from the operating point written
  // by initWithPowerflow(). Both legacy governors use their PF-derived
  // mechanical equilibrium, avoiding an artificial power step at t = 0.
  simulationEMT.doInitFromNodesAndTerminals(true);
  simulationEMT.doSystemMatrixRecomputation(p.recomputeSystemMatrix);
  simulationEMT.addLogger(loggerEMT);

  // true = close the three-phase breaker and connect the two generator islands.
  simulationEMT.addEvent(
      SwitchEvent3Ph::make(p.breakerCloseTime, BR_tr_EMT, true));

  // SPDLOG_INFO("\nStarting Scenario A EMT simulation:"
  //             "\n  dt={} s"
  //             "\n  finalTime={} s"
  //             "\n  breaker closes at t={} s"
  //             "\n  initial breaker state=open"
  //             "\n  matrix recomputation={}",
  //             p.timeStep, p.finalTime, p.breakerCloseTime,
  //             p.recomputeSystemMatrix ? "ON" : "OFF");

  simulationEMT.run();

  return 0;
}
