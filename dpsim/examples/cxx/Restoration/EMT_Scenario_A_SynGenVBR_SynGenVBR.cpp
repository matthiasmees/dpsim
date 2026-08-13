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

struct Parameters {
  // Simulation
  String name = "EMT_Scenario_A_SynGenVBR_SynGenVBR";
  Real frequency = 50.0;
  Real pfTimeStep = 1.0;
  Real pfFinalTime = 1.0;
  Real timeStep = 100e-6;
  Real finalTime = 10.0;
  Real breakerCloseTime = 5.0;
  Bool recomputeSystemMatrix = true;
  Int logDownSampling = 10;

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
  // Cable
  // ---------------------------------------------------------------------------
  Real cableBaseVoltage = 220e3;
  Real cableResistance = 5.04669;
  Real cableInductance = 0.13523123641;
  Real cableCapacitance = 1.93522865e-6;
  Real cableConductance = 1e-15;

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
  auto busGasPF = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busBPF = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busAPF = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busPshaPF = SimNode<Complex>::make("BUS_psha", PhaseType::Single);
  auto busTrPshPF = SimNode<Complex>::make("BUS_tr_psh", PhaseType::Single);
  auto busPshPF = SimNode<Complex>::make("BUS_psh", PhaseType::Single);

  // PF GEN_gas: VD/slack generator for the gas-side island.
  auto gasGeneratorPF =
      SP::Ph1::SynchronGenerator::make("GEN_gas", Logger::Level::off);
  gasGeneratorPF->setParameters(
      p.gasRatedPower, p.gasRatedVoltage,
      0.0, // initial P estimate; VD P/Q are solved by NRP
      p.gasRatedVoltage, PowerflowBusType::VD);
  gasGeneratorPF->setBaseVoltage(p.gasRatedVoltage);

  // PF TR_gas
  auto trGasPF =
      SP::Ph1::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  trGasPF->setParameters(p.trGasVoltageLow, p.trGasVoltageHigh,
                         p.trGasRatedPower, p.trGasRatioMagnitude,
                         p.trGasRatioPhase, p.trGasResistance,
                         p.trGasInductance);
  trGasPF->setBaseVoltage(p.trGasVoltageHigh);

  // PF cable
  auto cablePF = SP::Ph1::PiLine::make("cable", Logger::Level::off);
  cablePF->setParameters(p.cableResistance, p.cableInductance,
                         p.cableCapacitance, p.cableConductance);
  cablePF->setBaseVoltage(p.cableBaseVoltage);

  // PF line_3
  auto line3PF = SP::Ph1::PiLine::make("line_3", Logger::Level::off);
  line3PF->setParameters(line3Resistance, line3Inductance, line3Capacitance,
                         p.line3Conductance);
  line3PF->setBaseVoltage(p.line3BaseVoltage);

  // PF breaker: actual SP switch, initially open.
  auto breakerPF = SP::Ph1::Switch::make("breaker_GEN_psh", Logger::Level::off);
  breakerPF->setParameters(p.breakerOpenResistance, p.breakerClosedResistance,
                           false);
  breakerPF->setBaseVoltage(p.trPshVoltageHigh);

  // PF TR_psh
  auto trPshPF =
      SP::Ph1::Transformer::make("TR_psh", "TR_psh", Logger::Level::off, true);
  trPshPF->setParameters(p.trPshVoltageLow, p.trPshVoltageHigh,
                         p.trPshRatedPower, p.trPshRatioMagnitude,
                         p.trPshRatioPhase, p.trPshResistance,
                         p.trPshInductance);
  trPshPF->setBaseVoltage(p.trPshVoltageHigh);

  // PF GEN_psh: VD/slack generator for the isolated PSH-side island.
  auto pshGeneratorPF =
      SP::Ph1::SynchronGenerator::make("GEN_psh", Logger::Level::off);
  pshGeneratorPF->setParameters(
      p.pshRatedPower, p.pshRatedVoltage,
      0.0, // initial P estimate; VD P/Q are solved by NRP
      p.pshRatedVoltage, PowerflowBusType::VD);
  pshGeneratorPF->setBaseVoltage(p.pshRatedVoltage);

  // PF connectivity: identical graph to EMT.
  gasGeneratorPF->connect({busGasPF});
  trGasPF->connect({busGasPF, busBPF});
  cablePF->connect({busBPF, busAPF});
  line3PF->connect({busAPF, busPshaPF});
  breakerPF->connect({busPshaPF, busTrPshPF});
  trPshPF->connect({busPshPF, busTrPshPF});
  pshGeneratorPF->connect({busPshPF});

  SystemTopology systemPF(
      p.frequency,
      SystemNodeList{busGasPF, busBPF, busAPF, busPshaPF, busTrPshPF, busPshPF},
      SystemComponentList{gasGeneratorPF, trGasPF, cablePF, line3PF, breakerPF,
                          trPshPF, pshGeneratorPF});

  auto loggerPF = DataLogger::make(pfSimName);
  loggerPF->logAttribute("BUS_gas_v", busGasPF->attribute("v"));
  loggerPF->logAttribute("BUS_b_v", busBPF->attribute("v"));
  loggerPF->logAttribute("BUS_a_v", busAPF->attribute("v"));
  loggerPF->logAttribute("BUS_psha_v", busPshaPF->attribute("v"));
  loggerPF->logAttribute("BUS_tr_psh_v", busTrPshPF->attribute("v"));
  loggerPF->logAttribute("BUS_psh_v", busPshPF->attribute("v"));
  loggerPF->logAttribute("breaker_GEN_psh_i", breakerPF->attribute("i_intf"));

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

  const Complex gasPowerPF = gasGeneratorPF->getApparentPower();
  const Complex pshPowerPF = pshGeneratorPF->getApparentPower();

  SPDLOG_INFO(
      "PF operating point:"
      "\n  GEN_gas: P={} MW, Q={} Mvar, |V|={} V RMS, angle={} deg"
      "\n  GEN_psh: P={} MW, Q={} Mvar, |V|={} V RMS, angle={} deg"
      "\n  breaker grid side: |V|={} V RMS, angle={} deg"
      "\n  breaker PSH side : |V|={} V RMS, angle={} deg"
      "\n  breaker mismatch : |dV|={} V RMS, dAngle={} deg",
      gasPowerPF.real() / 1e6, gasPowerPF.imag() / 1e6,
      Math::abs(busGasPF->singleVoltage()),
      Math::phase(busGasPF->singleVoltage()) * 180.0 / PI,
      pshPowerPF.real() / 1e6, pshPowerPF.imag() / 1e6,
      Math::abs(busPshPF->singleVoltage()),
      Math::phase(busPshPF->singleVoltage()) * 180.0 / PI,
      Math::abs(busPshaPF->singleVoltage()),
      Math::phase(busPshaPF->singleVoltage()) * 180.0 / PI,
      Math::abs(busTrPshPF->singleVoltage()),
      Math::phase(busTrPshPF->singleVoltage()) * 180.0 / PI,
      Math::abs(busPshaPF->singleVoltage() - busTrPshPF->singleVoltage()),
      Math::phase(busPshaPF->singleVoltage() / busTrPshPF->singleVoltage()) *
          180.0 / PI);

  // ==========================================================================
  // 2. EMT
  // ==========================================================================

  const String emtSimName = p.name + "_EMT";
  Logger::setLogDir("logs/" + emtSimName);

  // EMT nodes. No manual PF-result phasors are passed here.
  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busA = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto busPsha = SimNode<Real>::make("BUS_psha", PhaseType::ABC);
  auto busTrPsh = SimNode<Real>::make("BUS_tr_psh", PhaseType::ABC);
  auto busPsh = SimNode<Real>::make("BUS_psh", PhaseType::ABC);

  // ---------------------------------------------------------------------------
  // EMT GEN_gas
  // ---------------------------------------------------------------------------
  auto gasGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_gas", Logger::Level::off);

  gasGenerator->setBaseAndOperationalPerUnitParameters(
      p.gasRatedPower, p.gasRatedVoltage, p.gasNominalFrequency,
      p.gasPoleNumber, p.gasNominalFieldCurrent, p.gasStatorResistance, p.gasLd,
      p.gasLq, p.gasLdTransient, p.gasLqTransient, p.gasLdSubtransient,
      p.gasLqSubtransient, p.gasLeakageInductance, p.gasTd0Transient,
      p.gasTq0Transient, p.gasTd0Subtransient, p.gasTq0Subtransient,
      p.gasInertia);

  // Initialize the legacy governor directly from the PF operating point.
  gasGenerator->addGovernor(
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
  gasGenerator->addExciter(gasExciter, gasExciterParameters);

  // EMT TR_gas
  auto trGas =
      EMT::Ph3::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  trGas->setParameters(
      p.trGasVoltageLow, p.trGasVoltageHigh, p.trGasRatedPower,
      p.trGasRatioMagnitude, p.trGasRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trGasResistance),
      Math::singlePhaseParameterToThreePhase(p.trGasInductance));

  // EMT cable
  auto cable = EMT::Ph3::PiLine::make("cable", Logger::Level::off);
  cable->setParameters(
      Math::singlePhaseParameterToThreePhase(p.cableResistance),
      Math::singlePhaseParameterToThreePhase(p.cableInductance),
      Math::singlePhaseParameterToThreePhase(p.cableCapacitance),
      Math::singlePhaseParameterToThreePhase(p.cableConductance));

  // EMT line_3
  auto line3 = EMT::Ph3::PiLine::make("line_3", Logger::Level::off);
  line3->setParameters(
      Math::singlePhaseParameterToThreePhase(line3Resistance),
      Math::singlePhaseParameterToThreePhase(line3Inductance),
      Math::singlePhaseParameterToThreePhase(line3Capacitance),
      Math::singlePhaseParameterToThreePhase(p.line3Conductance));

  // EMT breaker: same topological position and same initial OPEN state as PF.
  auto breaker = EMT::Ph3::Switch::make("breaker_GEN_psh", Logger::Level::off);
  breaker->setParameters(
      Math::singlePhaseParameterToThreePhase(p.breakerOpenResistance),
      Math::singlePhaseParameterToThreePhase(p.breakerClosedResistance));
  breaker->openSwitch();

  // EMT TR_psh
  auto trPsh =
      EMT::Ph3::Transformer::make("TR_psh", "TR_psh", Logger::Level::off, true);
  trPsh->setParameters(
      p.trPshVoltageLow, p.trPshVoltageHigh, p.trPshRatedPower,
      p.trPshRatioMagnitude, p.trPshRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trPshResistance),
      Math::singlePhaseParameterToThreePhase(p.trPshInductance));

  // ---------------------------------------------------------------------------
  // EMT GEN_psh
  // ---------------------------------------------------------------------------
  auto pshGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_psh", Logger::Level::off);

  pshGenerator->setBaseAndOperationalPerUnitParameters(
      p.pshRatedPower, p.pshRatedVoltage, p.pshNominalFrequency,
      p.pshPoleNumber, p.pshNominalFieldCurrent, p.pshStatorResistance, p.pshLd,
      p.pshLq, p.pshLdTransient, p.pshLqTransient, p.pshLdSubtransient,
      p.pshLqSubtransient, p.pshLeakageInductance, p.pshTd0Transient,
      p.pshTq0Transient, p.pshTd0Subtransient, p.pshTq0Subtransient,
      p.pshInertia);

  pshGenerator->addGovernor(
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
  pshGenerator->addExciter(pshExciter, pshExciterParameters);

  // EMT connectivity: exactly the same graph and names as PF.
  gasGenerator->connect({busGas});
  trGas->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPsha});
  breaker->connect({busPsha, busTrPsh});
  trPsh->connect({busPsh, busTrPsh});
  pshGenerator->connect({busPsh});

  SystemTopology systemEMT(
      p.frequency,
      SystemNodeList{busGas, busB, busA, busPsha, busTrPsh, busPsh},
      SystemComponentList{gasGenerator, trGas, cable, line3, breaker, trPsh,
                          pshGenerator});

  // Central PF -> EMT operating-point transfer:
  //   - solved node voltages by exact node name
  //   - synchronous-generator terminal P/Q with generator sign conversion
  //
  // No manual generatorPositivePower(), terminal()->setPower(), or manually
  // supplied initial phasors are required in the example.
  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  auto loggerEMT = DataLogger::make(emtSimName, true, p.logDownSampling);

  loggerEMT->logAttribute("BUS_gas_v", busGas->attribute("v"));
  loggerEMT->logAttribute("BUS_b_v", busB->attribute("v"));
  loggerEMT->logAttribute("BUS_a_v", busA->attribute("v"));
  loggerEMT->logAttribute("BUS_psha_v", busPsha->attribute("v"));
  loggerEMT->logAttribute("BUS_tr_psh_v", busTrPsh->attribute("v"));
  loggerEMT->logAttribute("BUS_psh_v", busPsh->attribute("v"));

  loggerEMT->logAttribute("TR_gas_i", trGas->attribute("i_intf"));
  loggerEMT->logAttribute("cable_i", cable->attribute("i_intf"));
  loggerEMT->logAttribute("line_3_i", line3->attribute("i_intf"));
  loggerEMT->logAttribute("breaker_GEN_psh_i", breaker->attribute("i_intf"));
  loggerEMT->logAttribute("breaker_GEN_psh_v", breaker->attribute("v_intf"));
  loggerEMT->logAttribute("TR_psh_i", trPsh->attribute("i_intf"));

  loggerEMT->logAttribute("GEN_gas_w_r", gasGenerator->attribute("w_r"));
  loggerEMT->logAttribute("GEN_gas_P_elec", gasGenerator->attribute("P_elec"));
  loggerEMT->logAttribute("GEN_gas_Q_elec", gasGenerator->attribute("Q_elec"));
  loggerEMT->logAttribute("GEN_gas_P_mech", gasGenerator->attribute("P_mech"));
  loggerEMT->logAttribute("GEN_gas_i_intf", gasGenerator->attribute("i_intf"));
  loggerEMT->logAttribute("GEN_gas_v_intf", gasGenerator->attribute("v_intf"));
  loggerEMT->logAttribute("GEN_gas_delta_r",
                          gasGenerator->attribute("delta_r"));
  loggerEMT->logAttribute("GEN_gas_T_e", gasGenerator->attribute("T_e"));

  loggerEMT->logAttribute("GEN_psh_w_r", pshGenerator->attribute("w_r"));
  loggerEMT->logAttribute("GEN_psh_P_elec", pshGenerator->attribute("P_elec"));
  loggerEMT->logAttribute("GEN_psh_Q_elec", pshGenerator->attribute("Q_elec"));
  loggerEMT->logAttribute("GEN_psh_P_mech", pshGenerator->attribute("P_mech"));
  loggerEMT->logAttribute("GEN_psh_i_intf", pshGenerator->attribute("i_intf"));
  loggerEMT->logAttribute("GEN_psh_v_intf", pshGenerator->attribute("v_intf"));
  loggerEMT->logAttribute("GEN_psh_delta_r",
                          pshGenerator->attribute("delta_r"));
  loggerEMT->logAttribute("GEN_psh_T_e", pshGenerator->attribute("T_e"));

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
      SwitchEvent3Ph::make(p.breakerCloseTime, breaker, true));

  SPDLOG_INFO("\nStarting Scenario A EMT simulation:"
              "\n  dt={} s"
              "\n  finalTime={} s"
              "\n  breaker closes at t={} s"
              "\n  initial breaker state=open"
              "\n  matrix recomputation={}",
              p.timeStep, p.finalTime, p.breakerCloseTime,
              p.recomputeSystemMatrix ? "ON" : "OFF");

  simulationEMT.run();

  return 0;
}
