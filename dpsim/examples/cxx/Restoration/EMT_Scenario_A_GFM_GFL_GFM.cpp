// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <dpsim-models/EMT/EMT_Ph3_GFL_Siemens.h>
#include <dpsim-models/EMT/EMT_Ph3_GFM_Siemens.h>
#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_Switch.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/SP/SP_Ph1_Load.h>
#include <dpsim-models/SP/SP_Ph1_PiLine.h>
#include <dpsim-models/SP/SP_Ph1_Switch.h>
#include <dpsim-models/SP/SP_Ph1_SynchronGenerator.h>
#include <dpsim-models/SP/SP_Ph1_Transformer.h>

#include <cmath>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

// =============================================================================
// Parameters
// =============================================================================

struct Parameters {
  // Simulation
  String name = "EMT_Scenario_A_GFM_GFL_GFM";
  Real frequency = 50.0;
  Real pfTimeStep = 1.0;
  Real pfFinalTime = 1.0;
  Real timeStep = 100e-6;
  Real finalTime = 15.0;
  Real breakerCloseTime = 10.0;
  Real gflReactivePowerStepTime = 5.0;
  Real gflReactivePowerAfterStep = -25e6;
  Int logDownSampling = 10;

  // ---------------------------------------------------------------------------
  // Converter ratings
  // ---------------------------------------------------------------------------
  Real gasRatedPower = 50e6;
  Real gasRatedVoltage = 10.5e3;

  Real pshRatedPower = 200e6;
  Real pshRatedVoltage = 18.0e3;

  Real pvRatedPower = 50e6;
  Real pvRatedVoltage = 11.0e3;

  Real networkVoltage = 220e3;

  // ---------------------------------------------------------------------------
  // Siemens GFM controller
  // ---------------------------------------------------------------------------
  Real gfmActivePowerDroopPu = 0.02;
  Real gfmReactivePowerDroopPu = 0.0311;

  Real gfmActivePowerMeasurementTimeConstant = 0.1;
  Real gfmReactivePowerMeasurementTimeConstant = 0.1;

  Real gfmVoltageControllerKp = 0.52;
  Real gfmVoltageControllerKi = 1.16;
  Real gfmOutputCurrentFeedforwardGain = 1.0;

  Real gfmCurrentControllerKp = 0.74;
  Real gfmCurrentControllerKi = 1.19;
  Real gfmPccVoltageFeedforwardGain = 1.0;

  Real gfmPwmDelayTimeConstant = 0.0;

  Real gfmMinimumFrequencyPu = 0.80;
  Real gfmMaximumFrequencyPu = 1.20;
  Real gfmMaximumCurrentReferencePu = 2.0;
  Real gfmMaximumVoltageCommandPu = 1.20;

  // ---------------------------------------------------------------------------
  // Siemens GFL controller
  // ---------------------------------------------------------------------------
  Real gflFrequencyToActivePowerGainPu = 4.0;
  Real gflVoltageToReactivePowerGainPu = 2.0;

  Real gflPllKp = 0.449 / (9.1e-3 + 250e-6);
  Real gflPllKi = 250e-6 / (5.9 * (9.1e-3 + 250e-6));

  Real gflCurrentControllerKp = 0.74 / 10.0;
  // Retain the value that gave the well-damped response in the load-loss test.
  Real gflCurrentControllerKi = 11.9 / 10.0;
  Real gflPccVoltageFeedforwardGain = 1.0;

  Real gflActivePowerMeasurementTimeConstant = 0.05;
  Real gflReactivePowerMeasurementTimeConstant = 0.05;
  Real gflCurrentMeasurementTimeConstant = 0.5e-3;
  Real gflVoltageFeedforwardMeasurementTimeConstant = 1.0e-3;

  Real gflMinimumFrequencyPu = 0.80;
  Real gflMaximumFrequencyPu = 1.20;
  Real gflMaximumCurrentReferencePu = 1.50;
  Real gflMaximumVoltageCommandPu = 1.20;
  Real gflMinimumVoltageForCurrentReferencePu = 0.05;

  // ---------------------------------------------------------------------------
  // Common Siemens validation filter
  // ---------------------------------------------------------------------------
  Real filterReferencePower = 50e3;
  Real filterReferenceVoltage = 257.0;
  Real filterReferenceFrequency = 50.0;

  Real filterResistance = 0.0279;
  Real filterInductance = 2.72e-4;
  Real filterCapacitance = 3.5e-6;

  Real filterBaseImpedance() const {
    return filterReferenceVoltage * filterReferenceVoltage /
           filterReferencePower;
  }

  Real filterBaseOmega() const { return 2.0 * PI * filterReferenceFrequency; }

  Real filterResistancePu() const {
    return filterResistance / filterBaseImpedance();
  }

  Real filterInductiveReactancePu() const {
    return filterBaseOmega() * filterInductance / filterBaseImpedance();
  }

  Real filterCapacitiveSusceptancePu() const {
    return filterBaseOmega() * filterCapacitance * filterBaseImpedance();
  }

  // ---------------------------------------------------------------------------
  // TR_gas
  // ---------------------------------------------------------------------------
  Real trGasVoltageLow = 10.5e3;
  Real trGasVoltageHigh = 220e3;
  Real trGasRatedPower = 50e6;
  Real trGasRatioMagnitude = 10.5e3 / 220e3;
  Real trGasRatioPhase = 0.0;
  Real trGasResistance = 3.64157728;
  Real trGasInductance = 0.31037265855769886;

  // ---------------------------------------------------------------------------
  // Cable BUS_b -- BUS_a
  // ---------------------------------------------------------------------------
  Real cableBaseVoltage = 220e3;
  Real cableResistance = 5.04669;
  Real cableInductance = 0.13523123641;
  Real cableCapacitance = 1.93522865e-6;
  Real cableConductance = 1e-15;

  // ---------------------------------------------------------------------------
  // LINE_3 BUS_a -- BUS_psha_grid
  // ---------------------------------------------------------------------------
  Real line3BaseVoltage = 220e3;
  Real line3Resistance = 0.0749 * 5.0;
  Real line3Inductance = 1.270693e-3 * 5.0;
  Real line3Capacitance = 0.00466961e-6 * 5.0;
  Real line3Conductance = 1e-15;

  // ---------------------------------------------------------------------------
  // PSH breaker
  // ---------------------------------------------------------------------------
  Real breakerOpenResistance = 1e12;
  Real breakerClosedResistance = 1e-3;

  // ---------------------------------------------------------------------------
  // TR_psh
  // ---------------------------------------------------------------------------
  Real trPshVoltageLow = 18.0e3;
  Real trPshVoltageHigh = 220e3;
  Real trPshRatedPower = 200e6;
  Real trPshRatioMagnitude = 18.0e3 / 220e3;
  Real trPshRatioPhase = 0.0;
  Real trPshResistance = 0.8349;
  Real trPshInductance = 0.0962521551189048;

  // ---------------------------------------------------------------------------
  // TR_pv
  // ---------------------------------------------------------------------------
  Real trPvVoltageHigh = 220e3;
  Real trPvVoltageLow = 11.0e3;
  Real trPvRatedPower = 50e6;
  Real trPvRatioMagnitude = 220e3 / 11.0e3;
  Real trPvRatioPhase = 0.0;
  Real trPvResistance = 0.0;
  Real trPvInductance = 0.928e-3;
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

    if (args.options.find("GFL_Q_STEP_TIME") != args.options.end())
      p.gflReactivePowerStepTime = args.getOptionReal("GFL_Q_STEP_TIME");

    if (args.options.find("GFL_Q_AFTER_STEP") != args.options.end())
      p.gflReactivePowerAfterStep = args.getOptionReal("GFL_Q_AFTER_STEP");
  }

  if (!(p.frequency > 0.0) || !(p.timeStep > 0.0) || !(p.finalTime > 0.0) ||
      !(p.logDownSampling > 0) || !(p.breakerCloseTime >= 0.0) ||
      !(p.breakerCloseTime < p.finalTime) ||
      !(p.gflReactivePowerStepTime >= 0.0) ||
      !(p.gflReactivePowerStepTime < p.finalTime)) {
    throw std::invalid_argument(
        "Require frequency>0, dt>0, finalTime>0, logDownSampling>0, "
        "0<=breakerCloseTime<finalTime and "
        "0<=gflReactivePowerStepTime<finalTime.");
  }

  SPDLOG_INFO("Scenario A Siemens GFM/GFL:"
              "\n  gas GFM: S={} VA, V={} V_LL RMS"
              "\n  PSH GFM: S={} VA, V={} V_LL RMS"
              "\n  PV GFL : S={} VA, V={} V_LL RMS"
              "\n  GFL droop: K_fP={}, K_VQ={}"
              "\n  GFL voltage-feedforward filter: tau={} s"
              "\n  filter: Rf={} pu, X_Lf={} pu, B_Cf={} pu"
              "\n  breaker closes at t={} s"
              "\n  GFL Q_ref steps to {} Mvar at t={} s",
              p.gasRatedPower, p.gasRatedVoltage, p.pshRatedPower,
              p.pshRatedVoltage, p.pvRatedPower, p.pvRatedVoltage,
              p.gflFrequencyToActivePowerGainPu,
              p.gflVoltageToReactivePowerGainPu,
              p.gflVoltageFeedforwardMeasurementTimeConstant,
              p.filterResistancePu(), p.filterInductiveReactancePu(),
              p.filterCapacitiveSusceptancePu(), p.breakerCloseTime,
              p.gflReactivePowerAfterStep / 1e6, p.gflReactivePowerStepTime);

  // ==========================================================================
  // 1. POWER FLOW
  //
  // PF and EMT use the same electrical graph:
  //
  // GEN_gas -- BUS_gas -- TR_gas -- BUS_b -- cable -- BUS_a
  //                                                    |    |
  //                                                    |   TR_pv
  //                                                    |    |
  //                                                  line_3 BUS_pv -- pv
  //                                                    |
  //                                             BUS_psha_grid
  //                                                    |
  //                                             breaker_GEN_psh
  //                                                    |
  //                                              BUS_tr_psh
  //                                                    |
  //                                                TR_psh
  //                                                    |
  //                                                BUS_psh
  //                                                    |
  //                                                GEN_psh
  //
  // The breaker is OPEN during PF and EMT initialization. Therefore the gas
  // and PSH sides are independent islands and each island uses one VD source.
  //
  // GFM_Siemens initializes from generator-positive terminal power. The generic
  // SynchronGenerator PF transfer uses the dynamic consumer-positive
  // convention, so the PF GFM surrogates deliberately have "_PF" names and
  // their solved P/Q are assigned explicitly after initWithPowerflow().
  //
  // The PV GFL starts at P=Q=0. Its PF surrogate is an SP load with the same
  // name "pv"; zero power is sign-convention neutral.
  // ==========================================================================

  const String pfSimName = p.name + "_PF";
  Logger::setLogDir("logs/" + pfSimName);

  // PF nodes
  auto busGasPF = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busBPF = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busAPF = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busPshaGridPF =
      SimNode<Complex>::make("BUS_psha_grid", PhaseType::Single);
  auto busTrPshPF = SimNode<Complex>::make("BUS_tr_psh", PhaseType::Single);
  auto busPshPF = SimNode<Complex>::make("BUS_psh", PhaseType::Single);
  auto busPvPF = SimNode<Complex>::make("BUS_pv", PhaseType::Single);

  // PF gas-side VD source
  auto gasGeneratorPF =
      SP::Ph1::SynchronGenerator::make("GEN_gas_PF", Logger::Level::off);
  gasGeneratorPF->setParameters(p.gasRatedPower, p.gasRatedVoltage, 0.0,
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

  // PF LINE_3
  auto line3PF = SP::Ph1::PiLine::make("line_3", Logger::Level::off);
  line3PF->setParameters(p.line3Resistance, p.line3Inductance,
                         p.line3Capacitance, p.line3Conductance);
  line3PF->setBaseVoltage(p.line3BaseVoltage);

  // PF breaker: actual SP switch, initially OPEN
  auto breakerPF = SP::Ph1::Switch::make("breaker_GEN_psh", Logger::Level::off);
  breakerPF->setParameters(p.breakerOpenResistance, p.breakerClosedResistance,
                           false);
  breakerPF->setBaseVoltage(p.networkVoltage);

  // PF TR_psh
  auto trPshPF =
      SP::Ph1::Transformer::make("TR_psh", "TR_psh", Logger::Level::off, true);
  trPshPF->setParameters(p.trPshVoltageLow, p.trPshVoltageHigh,
                         p.trPshRatedPower, p.trPshRatioMagnitude,
                         p.trPshRatioPhase, p.trPshResistance,
                         p.trPshInductance);
  trPshPF->setBaseVoltage(p.trPshVoltageHigh);

  // PF PSH-side VD source because the breaker is open.
  auto pshGeneratorPF =
      SP::Ph1::SynchronGenerator::make("GEN_psh_PF", Logger::Level::off);
  pshGeneratorPF->setParameters(p.pshRatedPower, p.pshRatedVoltage, 0.0,
                                p.pshRatedVoltage, PowerflowBusType::VD);
  pshGeneratorPF->setBaseVoltage(p.pshRatedVoltage);

  // PF TR_pv: system side -> converter side.
  auto trPvPF =
      SP::Ph1::Transformer::make("TR_pv", "TR_pv", Logger::Level::off, false);
  trPvPF->setParameters(p.trPvVoltageHigh, p.trPvVoltageLow, p.trPvRatedPower,
                        p.trPvRatioMagnitude, p.trPvRatioPhase,
                        p.trPvResistance, p.trPvInductance);
  trPvPF->setBaseVoltage(p.trPvVoltageHigh);

  // Zero-power PF surrogate for the GFL.
  auto pvPF = SP::Ph1::Load::make("pv", Logger::Level::off);
  pvPF->setParameters(0.0, 0.0, p.pvRatedVoltage);
  pvPF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // PF connectivity
  gasGeneratorPF->connect({busGasPF});
  trGasPF->connect({busGasPF, busBPF});
  cablePF->connect({busBPF, busAPF});
  line3PF->connect({busAPF, busPshaGridPF});
  breakerPF->connect({busPshaGridPF, busTrPshPF});
  trPshPF->connect({busPshPF, busTrPshPF});
  pshGeneratorPF->connect({busPshPF});
  trPvPF->connect({busAPF, busPvPF});
  pvPF->connect({busPvPF});

  SystemTopology systemPF(
      p.frequency,
      SystemNodeList{busGasPF, busBPF, busAPF, busPshaGridPF, busTrPshPF,
                     busPshPF, busPvPF},
      SystemComponentList{gasGeneratorPF, trGasPF, cablePF, line3PF, breakerPF,
                          trPshPF, pshGeneratorPF, trPvPF, pvPF});

  auto loggerPF = DataLogger::make(pfSimName);
  loggerPF->logAttribute("BUS_gas_v", busGasPF->attribute("v"));
  loggerPF->logAttribute("BUS_a_v", busAPF->attribute("v"));
  loggerPF->logAttribute("BUS_psha_grid_v", busPshaGridPF->attribute("v"));
  loggerPF->logAttribute("BUS_tr_psh_v", busTrPshPF->attribute("v"));
  loggerPF->logAttribute("BUS_psh_v", busPshPF->attribute("v"));
  loggerPF->logAttribute("BUS_pv_v", busPvPF->attribute("v"));

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
  const Complex pvPowerPF(0.0, 0.0);

  SPDLOG_INFO(
      "PF operating point:"
      "\n  gas GFM: P={} MW, Q={} Mvar, |V|={} V RMS"
      "\n  PSH GFM: P={} MW, Q={} Mvar, |V|={} V RMS"
      "\n  PV GFL : P={} MW, Q={} Mvar, |V|={} V RMS"
      "\n  breaker grid side: |V|={} V RMS, angle={} deg"
      "\n  breaker PSH side : |V|={} V RMS, angle={} deg"
      "\n  breaker mismatch : |dV|={} V RMS, dAngle={} deg",
      gasPowerPF.real() / 1e6, gasPowerPF.imag() / 1e6,
      Math::abs(busGasPF->singleVoltage()), pshPowerPF.real() / 1e6,
      pshPowerPF.imag() / 1e6, Math::abs(busPshPF->singleVoltage()),
      pvPowerPF.real() / 1e6, pvPowerPF.imag() / 1e6,
      Math::abs(busPvPF->singleVoltage()),
      Math::abs(busPshaGridPF->singleVoltage()),
      Math::phase(busPshaGridPF->singleVoltage()) * 180.0 / PI,
      Math::abs(busTrPshPF->singleVoltage()),
      Math::phase(busTrPshPF->singleVoltage()) * 180.0 / PI,
      Math::abs(busPshaGridPF->singleVoltage() - busTrPshPF->singleVoltage()),
      Math::phase(busPshaGridPF->singleVoltage() /
                  busTrPshPF->singleVoltage()) *
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
  auto busTrPsh = SimNode<Real>::make("BUS_tr_psh", PhaseType::ABC);
  auto busPsh = SimNode<Real>::make("BUS_psh", PhaseType::ABC);
  auto busPv = SimNode<Real>::make("BUS_pv", PhaseType::ABC);

  // ---------------------------------------------------------------------------
  // EMT gas GFM
  // ---------------------------------------------------------------------------
  auto gasGfm =
      EMT::Ph3::GFM_Siemens::make("GEN_gas", "GEN_gas", Logger::Level::off);

  gasGfm->setBaseParameters(p.gasRatedPower, p.gasRatedVoltage, p.frequency);

  gasGfm->setReferencesPerUnit(
      1.0, Math::abs(busGasPF->singleVoltage()) / p.gasRatedVoltage,
      gasPowerPF.real() / p.gasRatedPower, gasPowerPF.imag() / p.gasRatedPower);

  gasGfm->setDroopParametersPerUnit(p.gfmActivePowerDroopPu,
                                    p.gfmReactivePowerDroopPu);

  gasGfm->setPowerMeasurementFilterTimeConstants(
      p.gfmActivePowerMeasurementTimeConstant,
      p.gfmReactivePowerMeasurementTimeConstant);

  gasGfm->setVoltageControllerParameters(p.gfmVoltageControllerKp,
                                         p.gfmVoltageControllerKi,
                                         p.gfmOutputCurrentFeedforwardGain);

  gasGfm->setCurrentControllerParameters(p.gfmCurrentControllerKp,
                                         p.gfmCurrentControllerKi,
                                         p.gfmPccVoltageFeedforwardGain);

  gasGfm->setPwmDelayTimeConstant(p.gfmPwmDelayTimeConstant);

  gasGfm->setControllerLimitsPerUnit(
      p.gfmMinimumFrequencyPu, p.gfmMaximumFrequencyPu,
      p.gfmMaximumCurrentReferencePu, p.gfmMaximumVoltageCommandPu);

  gasGfm->setFilterParametersPerUnit(p.filterInductiveReactancePu(),
                                     p.filterCapacitiveSusceptancePu(),
                                     p.filterResistancePu());

  gasGfm->withControl(true);

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

  // EMT LINE_3
  auto line3 = EMT::Ph3::PiLine::make("line_3", Logger::Level::off);
  line3->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line3Resistance),
      Math::singlePhaseParameterToThreePhase(p.line3Inductance),
      Math::singlePhaseParameterToThreePhase(p.line3Capacitance),
      Math::singlePhaseParameterToThreePhase(p.line3Conductance));

  // EMT breaker: same initial OPEN state as PF.
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
  // EMT PSH GFM
  // ---------------------------------------------------------------------------
  auto pshGfm =
      EMT::Ph3::GFM_Siemens::make("GEN_psh", "GEN_psh", Logger::Level::off);

  pshGfm->setBaseParameters(p.pshRatedPower, p.pshRatedVoltage, p.frequency);

  pshGfm->setReferencesPerUnit(
      1.0, Math::abs(busPshPF->singleVoltage()) / p.pshRatedVoltage,
      pshPowerPF.real() / p.pshRatedPower, pshPowerPF.imag() / p.pshRatedPower);

  pshGfm->setDroopParametersPerUnit(p.gfmActivePowerDroopPu,
                                    p.gfmReactivePowerDroopPu);

  pshGfm->setPowerMeasurementFilterTimeConstants(
      p.gfmActivePowerMeasurementTimeConstant,
      p.gfmReactivePowerMeasurementTimeConstant);

  pshGfm->setVoltageControllerParameters(p.gfmVoltageControllerKp,
                                         p.gfmVoltageControllerKi,
                                         p.gfmOutputCurrentFeedforwardGain);

  pshGfm->setCurrentControllerParameters(p.gfmCurrentControllerKp,
                                         p.gfmCurrentControllerKi,
                                         p.gfmPccVoltageFeedforwardGain);

  pshGfm->setPwmDelayTimeConstant(p.gfmPwmDelayTimeConstant);

  pshGfm->setControllerLimitsPerUnit(
      p.gfmMinimumFrequencyPu, p.gfmMaximumFrequencyPu,
      p.gfmMaximumCurrentReferencePu, p.gfmMaximumVoltageCommandPu);

  pshGfm->setFilterParametersPerUnit(p.filterInductiveReactancePu(),
                                     p.filterCapacitiveSusceptancePu(),
                                     p.filterResistancePu());

  pshGfm->withControl(true);

  // EMT TR_pv
  auto trPv =
      EMT::Ph3::Transformer::make("TR_pv", "TR_pv", Logger::Level::off, false);
  trPv->setParameters(p.trPvVoltageHigh, p.trPvVoltageLow, p.trPvRatedPower,
                      p.trPvRatioMagnitude, p.trPvRatioPhase,
                      Math::singlePhaseParameterToThreePhase(p.trPvResistance),
                      Math::singlePhaseParameterToThreePhase(p.trPvInductance));

  // ---------------------------------------------------------------------------
  // EMT PV GFL
  // ---------------------------------------------------------------------------
  auto pvGfl = EMT::Ph3::GFL_Siemens::make("pv", "pv", Logger::Level::off);

  pvGfl->setBaseParameters(p.pvRatedPower, p.pvRatedVoltage, p.frequency);

  pvGfl->setReferencesPerUnit(
      1.0, Math::abs(busPvPF->singleVoltage()) / p.pvRatedVoltage, 0.0, 0.0);

  pvGfl->setDroopParametersPerUnit(p.gflFrequencyToActivePowerGainPu,
                                   p.gflVoltageToReactivePowerGainPu);

  pvGfl->setPllParameters(p.gflPllKp, p.gflPllKi);

  pvGfl->setCurrentControllerParameters(p.gflCurrentControllerKp,
                                        p.gflCurrentControllerKi,
                                        p.gflPccVoltageFeedforwardGain);

  pvGfl->setMeasurementFilterTimeConstants(
      p.gflActivePowerMeasurementTimeConstant,
      p.gflReactivePowerMeasurementTimeConstant,
      p.gflCurrentMeasurementTimeConstant,
      p.gflVoltageFeedforwardMeasurementTimeConstant);

  pvGfl->setControllerLimitsPerUnit(
      p.gflMinimumFrequencyPu, p.gflMaximumFrequencyPu,
      p.gflMaximumCurrentReferencePu, p.gflMaximumVoltageCommandPu,
      p.gflMinimumVoltageForCurrentReferencePu);

  pvGfl->setFilterParametersPerUnit(p.filterInductiveReactancePu(),
                                    p.filterCapacitiveSusceptancePu(),
                                    p.filterResistancePu());

  pvGfl->withControl(true);

  // EMT connectivity: exact same graph as PF.
  gasGfm->connect({busGas});
  trGas->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPshaGrid});
  breaker->connect({busPshaGrid, busTrPsh});
  trPsh->connect({busPsh, busTrPsh});
  pshGfm->connect({busPsh});
  trPv->connect({busA, busPv});
  pvGfl->connect({busPv});

  SystemTopology systemEMT(
      p.frequency,
      SystemNodeList{busGas, busB, busA, busPshaGrid, busTrPsh, busPsh, busPv},
      SystemComponentList{gasGfm, trGas, cable, line3, breaker, trPsh, pshGfm,
                          trPv, pvGfl});

  // Central PF -> EMT voltage initialization.
  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  // GFM_Siemens uses generator-positive terminal power during controller
  // initialization. Assign the solved PF injections explicitly.
  gasGfm->terminal(0)->setPower(gasPowerPF);
  pshGfm->terminal(0)->setPower(pshPowerPF);

  // GFL starts at P=Q=0; keep the terminal operating point explicit.
  pvGfl->terminal(0)->setPower(pvPowerPF);

  auto loggerEMT = DataLogger::make(emtSimName, true, p.logDownSampling);

  // Network
  loggerEMT->logAttribute("BUS_gas_v", busGas->attribute("v"));
  loggerEMT->logAttribute("BUS_a_v", busA->attribute("v"));
  loggerEMT->logAttribute("BUS_b_v", busB->attribute("v"));
  loggerEMT->logAttribute("BUS_psha_grid_v", busPshaGrid->attribute("v"));
  loggerEMT->logAttribute("BUS_tr_psh_v", busTrPsh->attribute("v"));
  loggerEMT->logAttribute("BUS_psh_v", busPsh->attribute("v"));
  loggerEMT->logAttribute("BUS_pv_v", busPv->attribute("v"));
  loggerEMT->logAttribute("breaker_GEN_psh_i", breaker->attribute("i_intf"));

  // Gas GFM
  loggerEMT->logAttribute("GFM_gas_P_elec_pu", gasGfm->attribute("P_elec_pu"));
  loggerEMT->logAttribute("GFM_gas_Q_elec_pu", gasGfm->attribute("Q_elec_pu"));
  loggerEMT->logAttribute("GFM_gas_frequency_pu",
                          gasGfm->attribute("frequency_pu"));
  loggerEMT->logAttribute("GFM_gas_V_magnitude_pu",
                          gasGfm->attribute("V_magnitude_pu"));

  // PSH GFM
  loggerEMT->logAttribute("GFM_psh_P_elec_pu", pshGfm->attribute("P_elec_pu"));
  loggerEMT->logAttribute("GFM_psh_Q_elec_pu", pshGfm->attribute("Q_elec_pu"));
  loggerEMT->logAttribute("GFM_psh_frequency_pu",
                          pshGfm->attribute("frequency_pu"));
  loggerEMT->logAttribute("GFM_psh_V_magnitude_pu",
                          pshGfm->attribute("V_magnitude_pu"));

  // PV GFL
  loggerEMT->logAttribute("GFL_pv_P_ref_pu", pvGfl->attribute("P_ref_pu"));
  loggerEMT->logAttribute("GFL_pv_Q_ref_pu", pvGfl->attribute("Q_ref_pu"));
  loggerEMT->logAttribute("GFL_pv_P_command_pu",
                          pvGfl->attribute("P_command_pu"));
  loggerEMT->logAttribute("GFL_pv_Q_command_pu",
                          pvGfl->attribute("Q_command_pu"));
  loggerEMT->logAttribute("GFL_pv_P_elec_pu", pvGfl->attribute("P_elec_pu"));
  loggerEMT->logAttribute("GFL_pv_Q_elec_pu", pvGfl->attribute("Q_elec_pu"));
  loggerEMT->logAttribute("GFL_pv_frequency_pu",
                          pvGfl->attribute("frequency_pu"));
  loggerEMT->logAttribute("GFL_pv_V_magnitude_pu",
                          pvGfl->attribute("V_magnitude_pu"));
  loggerEMT->logAttribute("GFL_pv_i_pcc_dq_pu",
                          pvGfl->attribute("i_pcc_dq_pu"));
  loggerEMT->logAttribute("GFL_pv_i_ref_dq_pu",
                          pvGfl->attribute("i_ref_dq_pu"));
  loggerEMT->logAttribute("GFL_pv_v_cmd_dq_pu",
                          pvGfl->attribute("v_cmd_dq_pu"));
  loggerEMT->logAttribute("GFL_pv_v_pcc_filtered_dq_pu",
                          pvGfl->attribute("v_pcc_filtered_dq_pu"));

  Simulation simulationEMT(emtSimName, Logger::Level::info);
  simulationEMT.setSystem(systemEMT);
  simulationEMT.setTimeStep(p.timeStep);
  simulationEMT.setFinalTime(p.finalTime);
  simulationEMT.setDomain(Domain::EMT);

  // Dynamic states are initialized from the PF operating point.
  // Matrix recomputation remains in the solver's default Auto mode; the EMT
  // switch supports precomputed open/closed system matrices.
  simulationEMT.doInitFromNodesAndTerminals(true);
  simulationEMT.addLogger(loggerEMT);

  // GFL reactive-power reference step.
  const Real qStepTime =
      std::round(p.gflReactivePowerStepTime / p.timeStep) * p.timeStep;
  const Real qAfterStepPu = p.gflReactivePowerAfterStep / p.pvRatedPower;

  simulationEMT.addEvent(AttributeEvent<Real>::make(
      qStepTime, pvGfl->mReactivePowerRefPu, qAfterStepPu));

  // Close the PSH breaker and connect both GFM islands.
  simulationEMT.addEvent(
      SwitchEvent3Ph::make(p.breakerCloseTime, breaker, true));

  SPDLOG_INFO("\nStarting Scenario A Siemens EMT simulation:"
              "\n  dt={} s"
              "\n  finalTime={} s"
              "\n  GFL Q_ref step={} pu at t={} s"
              "\n  PSH breaker closes at t={} s"
              "\n  initial breaker state=open",
              p.timeStep, p.finalTime, qAfterStepPu, qStepTime,
              p.breakerCloseTime);

  simulationEMT.run();

  return 0;
}
