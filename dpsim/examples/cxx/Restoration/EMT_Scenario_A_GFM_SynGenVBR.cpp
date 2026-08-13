// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <dpsim-models/EMT/EMT_Ph3_GFM_Droop.h>
#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_Switch.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorVBR.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/SP/SP_Ph1_PiLine.h>
#include <dpsim-models/SP/SP_Ph1_Switch.h>
#include <dpsim-models/SP/SP_Ph1_SynchronGenerator.h>
#include <dpsim-models/SP/SP_Ph1_Transformer.h>
#include <dpsim-models/Signal/Exciter.h>

#include <cmath>
#include <memory>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

// =============================================================================
// Parameters
// =============================================================================

struct Parameters {
  // Simulation
  String name = "EMT_GFM_SynGenVBR";
  Real frequency = 50.0;
  Real pfTimeStep = 1.0;
  Real pfFinalTime = 1.0;
  Real timeStep = 100e-6;
  Real finalTime = 10.0;
  Real breakerCloseTime = 2.0;
  Int logDownSampling = 10;

  // ---------------------------------------------------------------------------
  // Grid-forming converter GEN_gas
  // ---------------------------------------------------------------------------
  Real gfmRatedPower = 50e6;
  Real gfmRatedVoltage = 10.5e3;

  // GFM P-f / Q-V droop controller
  Real gfmActivePowerDroopPu = 0.05;
  Real gfmReactivePowerDroopPu = 0.05;
  Real gfmVoltageIntegralGain = 40.0;

  // P/Q measurement filter
  Real gfmPowerMeasurementFilterCutoffFrequency = 10.0;

  // Controller limits
  Real gfmMinimumFrequencyPu = 0.90;
  Real gfmMaximumFrequencyPu = 1.10;
  Real gfmMinimumVoltagePu = 0.70;
  Real gfmMaximumVoltagePu = 1.30;

  // ---------------------------------------------------------------------------
  // GFM validation filter
  //
  // The original working filter is defined on a 20 kVA / 381 V base.
  // Its per-unit R, X_L, B_C and Rd are preserved when used on the
  // 50 MVA / 10.5 kV GFM base.
  // ---------------------------------------------------------------------------
  Real gfmFilterReferencePower = 20e3;
  Real gfmFilterReferenceVoltage = std::sqrt(3.0) * 220.0;
  Real gfmFilterReferenceFrequency = 50.0;

  Real gfmFilterInductance = 3.0e-3;
  Real gfmFilterCapacitance = 20.0e-6;
  Real gfmFilterResistance = 0.05;
  Real gfmFilterDampingResistance = 4.0;

  Real gfmFilterBaseImpedance() const {
    return gfmFilterReferenceVoltage * gfmFilterReferenceVoltage /
           gfmFilterReferencePower;
  }

  Real gfmFilterBaseOmega() const {
    return 2.0 * PI * gfmFilterReferenceFrequency;
  }

  Real gfmFilterInductiveReactancePu() const {
    return gfmFilterBaseOmega() * gfmFilterInductance /
           gfmFilterBaseImpedance();
  }

  Real gfmFilterCapacitiveSusceptancePu() const {
    return gfmFilterBaseOmega() * gfmFilterCapacitance *
           gfmFilterBaseImpedance();
  }

  Real gfmFilterResistancePu() const {
    return gfmFilterResistance / gfmFilterBaseImpedance();
  }

  Real gfmFilterDampingResistancePu() const {
    return gfmFilterDampingResistance / gfmFilterBaseImpedance();
  }

  Real gfmPowerFilterTimeConstant() const {
    return 1.0 / (2.0 * PI * gfmPowerMeasurementFilterCutoffFrequency);
  }

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
  // Cable / LINE_1
  // ---------------------------------------------------------------------------
  Real line1BaseVoltage = 220e3;
  Real line1Resistance = 5.04669;
  Real line1Inductance = 0.13523123641;
  Real line1Capacitance = 1.93522865e-6;
  Real line1Conductance = 1e-15;

  // ---------------------------------------------------------------------------
  // LINE_3
  // ---------------------------------------------------------------------------
  Real line3BaseVoltage = 220e3;
  Real line3Length = 5.0;
  Real line3ResistancePerLength = 0.0749;
  Real line3InductancePerLength = 1.270693e-3;
  Real line3CapacitancePerLength = 0.00466961e-6;
  Real line3Conductance = 1e-15;

  Real line3Resistance() const {
    return line3ResistancePerLength * line3Length;
  }

  Real line3Inductance() const {
    return line3InductancePerLength * line3Length;
  }

  Real line3Capacitance() const {
    return line3CapacitancePerLength * line3Length;
  }

  // ---------------------------------------------------------------------------
  // BREAKER_psh
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

  CommandLineArgs args(argc, argv);

  if (argc > 1) {
    p.timeStep = args.timeStep;
    p.finalTime = args.duration;

    if (args.name != "dpsim")
      p.name = args.name;

    if (args.options.find("BREAKER_CLOSE_TIME") != args.options.end())
      p.breakerCloseTime = args.getOptionReal("BREAKER_CLOSE_TIME");
  }

  if (!(p.frequency > 0.0) || !(p.timeStep > 0.0) || !(p.finalTime > 0.0) ||
      !(p.logDownSampling > 0) || !(p.breakerCloseTime >= 0.0) ||
      !(p.breakerCloseTime < p.finalTime)) {
    throw std::invalid_argument(
        "Require frequency>0, dt>0, finalTime>0, logDownSampling>0 and "
        "0<=breakerCloseTime<finalTime.");
  }

  SPDLOG_INFO(
      "GFM-SynGen breaker scenario:"
      "\n  GFM_gas: S_rated={} VA, V_rated={} V_LL RMS"
      "\n  GEN_psh: S_rated={} VA, V_rated={} V_LL RMS, H={} s"
      "\n  TR_gas : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
      "\n  TR_psh : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
      "\n  breaker: R_open={} Ohm, R_closed={} Ohm"
      "\n  GFM droop: k_p={} pu, k_q={} pu, k_iv={} 1/s"
      "\n  GFM filter: Rf={} pu, X_Lf={} pu, B_Cf={} pu, Rd={} pu",
      p.gfmRatedPower, p.gfmRatedVoltage, p.pshRatedPower, p.pshRatedVoltage,
      p.pshInertia, p.trGasVoltageLow, p.trGasVoltageHigh, p.trGasRatedPower,
      p.trGasResistance, p.trGasInductance, p.trPshVoltageLow,
      p.trPshVoltageHigh, p.trPshRatedPower, p.trPshResistance,
      p.trPshInductance, p.breakerOpenResistance, p.breakerClosedResistance,
      p.gfmActivePowerDroopPu, p.gfmReactivePowerDroopPu,
      p.gfmVoltageIntegralGain, p.gfmFilterResistancePu(),
      p.gfmFilterInductiveReactancePu(), p.gfmFilterCapacitiveSusceptancePu(),
      p.gfmFilterDampingResistancePu());

  // ==========================================================================
  // 1. POWER FLOW
  //
  // PF and EMT have exactly the same electrical graph:
  //
  // GEN_gas -- BUS_gas -- TR_gas -- BUS_b -- LINE_1 -- BUS_a
  //                                                        |
  //                                                     LINE_3
  //                                                        |
  //                                                BUS_psha_grid
  //                                                        |
  //                                                 BREAKER_psh
  //                                                        |
  //                                                 BUS_psha_hv
  //                                                        |
  //                                                    TR_psh
  //                                                        |
  //                                                  BUS_psha
  //                                                        |
  //                                                   GEN_psh
  //
  // The breaker is OPEN during PF and EMT initialization, so the gas/GFM and
  // PSH sides are two independent islands. Each island has one VD source.
  //
  // Node names are identical in PF and EMT, allowing initWithPowerflow() to
  // transfer all solved bus phasors directly.
  //
  // GEN_psh also uses the same component name in PF and EMT, so its P/Q is
  // transferred automatically. The gas-side PF source is deliberately named
  // GEN_gas_PF because GFM_Droop currently uses generator-positive terminal
  // power for initialization, whereas the generic SynchronGenerator PF
  // transfer applies the standard consumer-positive dynamic sign convention.
  // ==========================================================================

  const String pfSimName = p.name + "_PF";
  Logger::setLogDir("logs/" + pfSimName);

  // PF nodes
  auto busGasPF = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busBPF = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busAPF = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busPshaGridPF =
      SimNode<Complex>::make("BUS_psha_grid", PhaseType::Single);
  auto busPshaHvPF = SimNode<Complex>::make("BUS_psha_hv", PhaseType::Single);
  auto busPshPF = SimNode<Complex>::make("BUS_psha", PhaseType::Single);

  // PF gas-side VD source.
  auto gasGeneratorPF =
      SP::Ph1::SynchronGenerator::make("GEN_gas_PF", Logger::Level::off);
  gasGeneratorPF->setParameters(
      p.gfmRatedPower, p.gfmRatedVoltage,
      0.0, // initial P estimate; VD P/Q are solved by NRP
      p.gfmRatedVoltage, PowerflowBusType::VD);
  gasGeneratorPF->setBaseVoltage(p.gfmRatedVoltage);

  // PF TR_gas
  auto trGasPF =
      SP::Ph1::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  trGasPF->setParameters(p.trGasVoltageLow, p.trGasVoltageHigh,
                         p.trGasRatedPower, p.trGasRatioMagnitude,
                         p.trGasRatioPhase, p.trGasResistance,
                         p.trGasInductance);
  trGasPF->setBaseVoltage(p.trGasVoltageHigh);

  // PF LINE_1 / cable
  auto line1PF = SP::Ph1::PiLine::make("LINE_1", Logger::Level::off);
  line1PF->setParameters(p.line1Resistance, p.line1Inductance,
                         p.line1Capacitance, p.line1Conductance);
  line1PF->setBaseVoltage(p.line1BaseVoltage);

  // PF LINE_3
  auto line3PF = SP::Ph1::PiLine::make("LINE_3", Logger::Level::off);
  line3PF->setParameters(p.line3Resistance(), p.line3Inductance(),
                         p.line3Capacitance(), p.line3Conductance);
  line3PF->setBaseVoltage(p.line3BaseVoltage);

  // PF breaker: actual SP switch, initially OPEN.
  auto breakerPF = SP::Ph1::Switch::make("BREAKER_psh", Logger::Level::off);
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

  // PF GEN_psh: VD/slack generator for the isolated PSH island.
  auto pshGeneratorPF =
      SP::Ph1::SynchronGenerator::make("GEN_psh", Logger::Level::off);
  pshGeneratorPF->setParameters(
      p.pshRatedPower, p.pshRatedVoltage,
      0.0, // initial P estimate; VD P/Q are solved by NRP
      p.pshRatedVoltage, PowerflowBusType::VD);
  pshGeneratorPF->setBaseVoltage(p.pshRatedVoltage);

  // PF connectivity
  gasGeneratorPF->connect({busGasPF});
  trGasPF->connect({busGasPF, busBPF});
  line1PF->connect({busBPF, busAPF});
  line3PF->connect({busAPF, busPshaGridPF});
  breakerPF->connect({busPshaGridPF, busPshaHvPF});
  trPshPF->connect({busPshPF, busPshaHvPF});
  pshGeneratorPF->connect({busPshPF});

  SystemTopology systemPF(p.frequency,
                          SystemNodeList{busGasPF, busBPF, busAPF,
                                         busPshaGridPF, busPshaHvPF, busPshPF},
                          SystemComponentList{gasGeneratorPF, trGasPF, line1PF,
                                              line3PF, breakerPF, trPshPF,
                                              pshGeneratorPF});

  auto loggerPF = DataLogger::make(pfSimName);
  loggerPF->logAttribute("BUS_gas_v", busGasPF->attribute("v"));
  loggerPF->logAttribute("BUS_b_v", busBPF->attribute("v"));
  loggerPF->logAttribute("BUS_a_v", busAPF->attribute("v"));
  loggerPF->logAttribute("BUS_psha_grid_v", busPshaGridPF->attribute("v"));
  loggerPF->logAttribute("BUS_psha_hv_v", busPshaHvPF->attribute("v"));
  loggerPF->logAttribute("BUS_psha_v", busPshPF->attribute("v"));

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
      "\n  GFM_gas: P={} MW, Q={} Mvar, |V|={} V RMS, angle={} deg"
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
      Math::abs(busPshaGridPF->singleVoltage()),
      Math::phase(busPshaGridPF->singleVoltage()) * 180.0 / PI,
      Math::abs(busPshaHvPF->singleVoltage()),
      Math::phase(busPshaHvPF->singleVoltage()) * 180.0 / PI,
      Math::abs(busPshaGridPF->singleVoltage() - busPshaHvPF->singleVoltage()),
      Math::phase(busPshaGridPF->singleVoltage() /
                  busPshaHvPF->singleVoltage()) *
          180.0 / PI);

  // ==========================================================================
  // 2. EMT
  // ==========================================================================

  const String emtSimName = p.name + "_EMT";
  Logger::setLogDir("logs/" + emtSimName);

  // EMT nodes: exact same names as PF.
  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busA = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto busPshaGrid = SimNode<Real>::make("BUS_psha_grid", PhaseType::ABC);
  auto busPshaHv = SimNode<Real>::make("BUS_psha_hv", PhaseType::ABC);
  auto busPsh = SimNode<Real>::make("BUS_psha", PhaseType::ABC);

  // ---------------------------------------------------------------------------
  // EMT GFM GEN_gas
  // ---------------------------------------------------------------------------
  auto gasGfm = EMT::Ph3::GFM_Droop::make("GEN_gas", "GEN_gas",
                                          Logger::Level::off, false);

  const Real gfmInitialVoltagePu =
      Math::abs(busGasPF->singleVoltage()) / p.gfmRatedVoltage;

  gasGfm->setBaseParameters(p.gfmRatedPower, p.gfmRatedVoltage, p.frequency);

  gasGfm->setParametersPerUnit(1.0, gfmInitialVoltagePu,
                               gasPowerPF.real() / p.gfmRatedPower,
                               gasPowerPF.imag() / p.gfmRatedPower);

  gasGfm->setDroopParametersPerUnit(p.gfmActivePowerDroopPu,
                                    p.gfmReactivePowerDroopPu,
                                    p.gfmVoltageIntegralGain);

  gasGfm->setPowerFilterTimeConstant(p.gfmPowerFilterTimeConstant());

  gasGfm->setControllerLimitsPerUnit(
      p.gfmMinimumFrequencyPu, p.gfmMaximumFrequencyPu, p.gfmMinimumVoltagePu,
      p.gfmMaximumVoltagePu);

  gasGfm->setFilterParametersPerUnit(
      p.gfmFilterInductiveReactancePu(), p.gfmFilterCapacitiveSusceptancePu(),
      p.gfmFilterResistancePu(), p.gfmFilterDampingResistancePu());

  gasGfm->withControl(true);

  // EMT TR_gas
  auto trGas =
      EMT::Ph3::Transformer::make("TR_gas", "TR_gas", Logger::Level::off, true);
  trGas->setParameters(
      p.trGasVoltageLow, p.trGasVoltageHigh, p.trGasRatedPower,
      p.trGasRatioMagnitude, p.trGasRatioPhase,
      Math::singlePhaseParameterToThreePhase(p.trGasResistance),
      Math::singlePhaseParameterToThreePhase(p.trGasInductance));

  // EMT LINE_1 / cable
  auto line1 = EMT::Ph3::PiLine::make("LINE_1", Logger::Level::off);
  line1->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line1Resistance),
      Math::singlePhaseParameterToThreePhase(p.line1Inductance),
      Math::singlePhaseParameterToThreePhase(p.line1Capacitance),
      Math::singlePhaseParameterToThreePhase(p.line1Conductance));

  // EMT LINE_3
  auto line3 = EMT::Ph3::PiLine::make("LINE_3", Logger::Level::off);
  line3->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line3Resistance()),
      Math::singlePhaseParameterToThreePhase(p.line3Inductance()),
      Math::singlePhaseParameterToThreePhase(p.line3Capacitance()),
      Math::singlePhaseParameterToThreePhase(p.line3Conductance));

  // EMT breaker: same position and same initial OPEN state as PF.
  auto breaker = EMT::Ph3::Switch::make("BREAKER_psh", Logger::Level::off);
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

  // EMT connectivity: exactly the same electrical graph as PF.
  gasGfm->connect({busGas});
  trGas->connect({busGas, busB});
  line1->connect({busB, busA});
  line3->connect({busA, busPshaGrid});
  breaker->connect({busPshaGrid, busPshaHv});
  trPsh->connect({busPsh, busPshaHv});
  pshGenerator->connect({busPsh});

  SystemTopology systemEMT(
      p.frequency,
      SystemNodeList{busGas, busB, busA, busPshaGrid, busPshaHv, busPsh},
      SystemComponentList{gasGfm, trGas, line1, line3, breaker, trPsh,
                          pshGenerator});

  // Central PF -> EMT initialization:
  //   - all solved node voltages are transferred by exact node name
  //   - GEN_psh P/Q are transferred automatically because the PF and EMT
  //     synchronous-generator components have the same name
  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  // GFM_Droop currently initializes from generator-positive terminal power,
  // unlike the generic consumer-positive dynamic terminal convention used by
  // initWithPowerflow() for SP::Ph1::SynchronGenerator. Therefore the gas-side
  // PF source is intentionally named GEN_gas_PF and its solved generator power
  // is assigned explicitly here.
  gasGfm->terminal(0)->setPower(gasPowerPF);

  auto loggerEMT = DataLogger::make(emtSimName, true, p.logDownSampling);

  loggerEMT->logAttribute("BUS_gas_v", busGas->attribute("v"));
  loggerEMT->logAttribute("BUS_b_v", busB->attribute("v"));
  loggerEMT->logAttribute("BUS_a_v", busA->attribute("v"));
  loggerEMT->logAttribute("BUS_psha_grid_v", busPshaGrid->attribute("v"));
  loggerEMT->logAttribute("BUS_psha_hv_v", busPshaHv->attribute("v"));
  loggerEMT->logAttribute("BUS_psha_v", busPsh->attribute("v"));

  loggerEMT->logAttribute("TR_gas_i", trGas->attribute("i_intf"));
  loggerEMT->logAttribute("LINE_1_i", line1->attribute("i_intf"));
  loggerEMT->logAttribute("LINE_3_i", line3->attribute("i_intf"));
  loggerEMT->logAttribute("BREAKER_psh_i", breaker->attribute("i_intf"));
  loggerEMT->logAttribute("BREAKER_psh_v", breaker->attribute("v_intf"));
  loggerEMT->logAttribute("TR_psh_i", trPsh->attribute("i_intf"));

  loggerEMT->logAttribute("GEN_gas_P_elec", gasGfm->attribute("P_elec"));
  loggerEMT->logAttribute("GEN_gas_Q_elec", gasGfm->attribute("Q_elec"));
  loggerEMT->logAttribute("GEN_gas_P_filtered",
                          gasGfm->attribute("P_filtered"));
  loggerEMT->logAttribute("GEN_gas_Q_filtered",
                          gasGfm->attribute("Q_filtered"));
  loggerEMT->logAttribute("GEN_gas_frequency", gasGfm->attribute("frequency"));
  loggerEMT->logAttribute("GEN_gas_theta", gasGfm->attribute("theta"));
  loggerEMT->logAttribute("GEN_gas_voltage_magnitude",
                          gasGfm->attribute("V_magnitude"));
  loggerEMT->logAttribute("GEN_gas_i_pcc", gasGfm->attribute("i_pcc"));
  loggerEMT->logAttribute("GEN_gas_v_source", gasGfm->attribute("Vs"));

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

  // Dynamic states are initialized from the PF operating point.
  // Matrix-recomputation mode is intentionally left at Auto; the breaker
  // supports precomputed open/closed system matrices.
  simulationEMT.doInitFromNodesAndTerminals(true);
  simulationEMT.addLogger(loggerEMT);

  // Close the breaker and connect the GFM and synchronous-generator islands.
  simulationEMT.addEvent(
      SwitchEvent3Ph::make(p.breakerCloseTime, breaker, true));

  SPDLOG_INFO("\nStarting GFM-SynGen breaker EMT simulation:"
              "\n  dt={} s"
              "\n  finalTime={} s"
              "\n  breaker closes at t={} s"
              "\n  initial breaker state=open"
              "\n  GFM initial P={} MW, Q={} Mvar"
              "\n  GEN_psh initial P={} MW, Q={} Mvar",
              p.timeStep, p.finalTime, p.breakerCloseTime,
              gasPowerPF.real() / 1e6, gasPowerPF.imag() / 1e6,
              pshPowerPF.real() / 1e6, pshPowerPF.imag() / 1e6);

  simulationEMT.run();

  return 0;
}
