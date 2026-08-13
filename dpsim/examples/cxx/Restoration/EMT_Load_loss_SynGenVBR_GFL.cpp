// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <dpsim-models/EMT/EMT_Ph3_GFL_Siemens.h>
#include <dpsim-models/EMT/EMT_Ph3_RXLoad.h>
#include <dpsim-models/EMT/EMT_Ph3_Switch.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorVBR.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/SP/SP_Ph1_Load.h>
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
  String name = "EMT_Load_loss_SynGenVBR_GFL";
  Real frequency = 50.0;
  Real pfTimeStep = 1.0;
  Real pfFinalTime = 1.0;
  Real timeStep = 100e-6;
  Real finalTime = 5.0;
  Real breakerOpenTime = 2.0;
  Real gflActivePowerStepTime = 2.2;
  Int logDownSampling = 10;

  // ---------------------------------------------------------------------------
  // Synchronous generator GEN_gas
  // ---------------------------------------------------------------------------
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

  // GEN_gas turbine governor
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

  // GEN_gas exciter
  Real excTa = 0.005;
  Real excKa = 200.0;
  Real excTe = 0.05;
  Real excKe = 0.5;
  Real excTf = 0.3;
  Real excKf = 0.01;
  Real excTr = 0.02;
  Real excMaximumRegulatorVoltage = 19.0;
  Real excMinimumRegulatorVoltage = -5.0;

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
  // Transformer TR_ld2
  // ---------------------------------------------------------------------------
  Real trLoadVoltageLow = 10.0e3;
  Real trLoadVoltageHigh = 220e3;
  Real trLoadRatedPower = 50e6;
  Real trLoadRatioMagnitude = 10.0e3 / 220e3;
  Real trLoadRatioPhase = 0.0;
  Real trLoadResistance = 3.64157728;
  Real trLoadInductance = 0.31037292;

  // ---------------------------------------------------------------------------
  // LOAD_2a
  // ---------------------------------------------------------------------------
  Real load2aActivePower = 5e6;
  Real load2aReactivePower = 3.098721e6;
  Real load2aNominalVoltage = 10e3;

  // ---------------------------------------------------------------------------
  // LOAD_2b
  // ---------------------------------------------------------------------------
  Real load2bActivePower = 15e6;
  Real load2bReactivePower = 0.0;
  Real load2bNominalVoltage = 10e3;

  // ---------------------------------------------------------------------------
  // BR_load_2b
  // ---------------------------------------------------------------------------
  Real breakerOpenResistance = 1e12;
  Real breakerClosedResistance = 1e-3;

  // ---------------------------------------------------------------------------
  // PV GFL and transformer
  // ---------------------------------------------------------------------------
  Real gflRatedPower = 50e6;
  Real gflRatedVoltage = 11.0e3;
  Real gflNetworkVoltage = 220e3;

  Real trPvResistance = 0.0;
  Real trPvInductance = 0.928e-3;

  // The GFL starts at zero P/Q and reacts to the load-loss disturbance
  // through its inverse frequency/voltage droop.
  Real gflInitialActivePower = 0.0;
  Real gflInitialReactivePower = 0.0;

  // GFL active-power reference step
  Real gflActivePowerAfterStep = -5e6;

  // GFL inverse-droop gains
  Real gflFrequencyToActivePowerGainPu = 4.0;
  Real gflVoltageToReactivePowerGainPu = 2.0;

  // PLL
  Real gflPllKp = 0.449 / (9.1e-3 + 250e-6);
  Real gflPllKi = 250e-6 / (5.9 * (9.1e-3 + 250e-6));

  // Current controller
  Real gflCurrentControllerKp = 0.74 / 10.0;
  Real gflCurrentControllerKi = 11.9 / 10.0;
  Real gflPccVoltageFeedforwardGain = 1.0;

  // Measurement filters
  Real gflActivePowerMeasurementTimeConstant = 0.05;
  Real gflReactivePowerMeasurementTimeConstant = 0.05;
  Real gflCurrentMeasurementTimeConstant = 0.5e-3;
  Real gflVoltageFeedforwardMeasurementTimeConstant = 1.0e-3;

  // Controller limits
  Real gflMinimumFrequencyPu = 0.80;
  Real gflMaximumFrequencyPu = 1.20;
  Real gflMaximumCurrentReferencePu = 1.50;
  Real gflMaximumVoltageCommandPu = 1.20;
  Real gflMinimumVoltageForCurrentReferencePu = 0.05;

  // ---------------------------------------------------------------------------
  // Siemens validation filter
  // ---------------------------------------------------------------------------
  Real gflFilterReferencePower = 50e3;
  Real gflFilterReferenceVoltage = 257.0;
  Real gflFilterReferenceFrequency = 50.0;

  Real gflFilterResistance = 0.0279;
  Real gflFilterInductance = 2.72e-4;
  Real gflFilterCapacitance = 3.5e-6;

  Real gflFilterBaseImpedance() const {
    return gflFilterReferenceVoltage * gflFilterReferenceVoltage /
           gflFilterReferencePower;
  }

  Real gflFilterBaseOmega() const {
    return 2.0 * PI * gflFilterReferenceFrequency;
  }

  Real gflFilterResistancePu() const {
    return gflFilterResistance / gflFilterBaseImpedance();
  }

  Real gflFilterInductiveReactancePu() const {
    return gflFilterBaseOmega() * gflFilterInductance /
           gflFilterBaseImpedance();
  }

  Real gflFilterCapacitiveSusceptancePu() const {
    return gflFilterBaseOmega() * gflFilterCapacitance *
           gflFilterBaseImpedance();
  }
};

int main(int argc, char *argv[]) {
  Parameters p;

  CommandLineArgs args(argc, argv);

  if (argc > 1) {
    p.timeStep = args.timeStep;
    p.finalTime = args.duration;

    if (args.name != "dpsim")
      p.name = args.name;

    if (args.options.find("BREAKER_OPEN_TIME") != args.options.end())
      p.breakerOpenTime = args.getOptionReal("BREAKER_OPEN_TIME");

    if (args.options.find("GFL_P_STEP_TIME") != args.options.end())
      p.gflActivePowerStepTime = args.getOptionReal("GFL_P_STEP_TIME");

    if (args.options.find("GFL_P_AFTER_STEP") != args.options.end())
      p.gflActivePowerAfterStep = args.getOptionReal("GFL_P_AFTER_STEP");
  }

  if (!(p.frequency > 0.0) || !(p.timeStep > 0.0) || !(p.finalTime > 0.0) ||
      !(p.logDownSampling > 0) || !(p.breakerOpenTime >= 0.0) ||
      !(p.breakerOpenTime < p.finalTime) ||
      !(p.gflActivePowerStepTime >= 0.0) ||
      !(p.gflActivePowerStepTime < p.finalTime)) {
    throw std::invalid_argument(
        "Require frequency>0, dt>0, finalTime>0, logDownSampling>0, "
        "0<=breakerOpenTime<finalTime and "
        "0<=gflActivePowerStepTime<finalTime.");
  }

  SPDLOG_INFO("Load-loss scenario:"
              "\n  GEN_gas: S_rated={} VA, V_rated={} V_LL RMS, H={} s"
              "\n  GFL_pv : S_rated={} VA, V_rated={} V_LL RMS"
              "\n  TR_gas : {}/{} V, S_rated={} VA, R={} Ohm, L={} H"
              "\n  load retained   : P={} MW, Q={} Mvar"
              "\n  load disconnected: P={} MW, Q={} Mvar"
              "\n  GFL P step: {} MW at t={} s"
              "\n  GFL droop: K_fP={} pu, K_VQ={} pu"
              "\n  GFL filter: Rf={} pu, X_Lf={} pu, B_Cf={} pu",
              p.genRatedPower, p.genRatedVoltage, p.genInertia, p.gflRatedPower,
              p.gflRatedVoltage, p.trGasVoltageLow, p.trGasVoltageHigh,
              p.trGasRatedPower, p.trGasResistance, p.trGasInductance,
              p.load2aActivePower / 1e6, p.load2aReactivePower / 1e6,
              p.load2bActivePower / 1e6, p.load2bReactivePower / 1e6,
              p.gflActivePowerAfterStep / 1e6, p.gflActivePowerStepTime,
              p.gflFrequencyToActivePowerGainPu,
              p.gflVoltageToReactivePowerGainPu, p.gflFilterResistancePu(),
              p.gflFilterInductiveReactancePu(),
              p.gflFilterCapacitiveSusceptancePu());

  // ==========================================================================
  // 1. POWER FLOW
  //
  // PF and EMT use exactly the same node/component names and topology:
  //
  // GEN_gas -- BUS_gas -- TR_gas -- BUS_b
  //                                  |    |
  //                                  |   TR_pv
  //                                  |    |
  //                                  |  BUS_pv
  //                                  |    |
  //                                  |    pv
  //                                  |
  //                                TR_ld2
  //                                  |
  //                            BUS_load_2
  //                              |       |
  //                           LOAD_2a  BR_load_2b
  //                                      |
  //                                BUS_load_2b
  //                                      |
  //                                   LOAD_2b
  //
  // BR_load_2b is CLOSED during PF and EMT initialization. At t = 5 s it
  // opens and disconnects LOAD_2b.
  //
  // The GFL is represented in PF by a zero-power SP::Ph1::Load named "pv".
  // Since the initial GFL operating point is P=Q=0, the generic load-power
  // transfer is sign-convention neutral and can initialize the matching EMT
  // GFL component directly by name.
  // ==========================================================================

  const String pfSimName = p.name + "_PF";
  Logger::setLogDir("logs/" + pfSimName);

  // PF nodes
  auto busGasPF = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busBPF = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busLoadPF = SimNode<Complex>::make("BUS_load_2", PhaseType::Single);
  auto busLoad2bPF = SimNode<Complex>::make("BUS_load_2b", PhaseType::Single);
  auto busPvPF = SimNode<Complex>::make("BUS_pv", PhaseType::Single);

  // PF GEN_gas: VD/slack generator.
  // P is only an initial estimate; NRP solves the converged slack P/Q.
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

  // PF breaker: actual SP switch, initially CLOSED.
  auto breakerPF = SP::Ph1::Switch::make("BR_load_2b", Logger::Level::off);
  breakerPF->setParameters(p.breakerOpenResistance, p.breakerClosedResistance,
                           true);
  breakerPF->setBaseVoltage(p.load2bNominalVoltage);

  // PF LOAD_2b behind the breaker
  auto load2bPF = SP::Ph1::Load::make("LOAD_2b", Logger::Level::off);
  load2bPF->setParameters(p.load2bActivePower, p.load2bReactivePower,
                          p.load2bNominalVoltage);
  load2bPF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // PF TR_pv: network side -> converter side.
  auto trPvPF =
      SP::Ph1::Transformer::make("TR_pv", "TR_pv", Logger::Level::off, false);
  trPvPF->setParameters(p.gflNetworkVoltage, p.gflRatedVoltage, p.gflRatedPower,
                        p.gflNetworkVoltage / p.gflRatedVoltage, 0.0,
                        p.trPvResistance, p.trPvInductance);
  trPvPF->setBaseVoltage(p.gflNetworkVoltage);

  // PF surrogate for the GFL.
  auto pvPF = SP::Ph1::Load::make("pv", Logger::Level::off);
  pvPF->setParameters(p.gflInitialActivePower, p.gflInitialReactivePower,
                      p.gflRatedVoltage);
  pvPF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // PF connectivity
  gasGeneratorPF->connect({busGasPF});
  trGasPF->connect({busGasPF, busBPF});
  trLoadPF->connect({busLoadPF, busBPF});
  load2aPF->connect({busLoadPF});
  breakerPF->connect({busLoadPF, busLoad2bPF});
  load2bPF->connect({busLoad2bPF});
  trPvPF->connect({busBPF, busPvPF});
  pvPF->connect({busPvPF});

  SystemTopology systemPF(
      p.frequency,
      SystemNodeList{busGasPF, busBPF, busLoadPF, busLoad2bPF, busPvPF},
      SystemComponentList{gasGeneratorPF, trGasPF, trLoadPF, load2aPF,
                          breakerPF, load2bPF, trPvPF, pvPF});

  auto loggerPF = DataLogger::make(pfSimName);
  loggerPF->logAttribute("BUS_gas_v", busGasPF->attribute("v"));
  loggerPF->logAttribute("BUS_b_v", busBPF->attribute("v"));
  loggerPF->logAttribute("BUS_load_2_v", busLoadPF->attribute("v"));
  loggerPF->logAttribute("BUS_load_2b_v", busLoad2bPF->attribute("v"));
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

  SPDLOG_INFO("PF operating point:"
              "\n  GEN_gas      : P={} MW, Q={} Mvar"
              "\n  V_GEN        : |V|={} V RMS, angle={} deg"
              "\n  V_BUS_b      : |V|={} V RMS, angle={} deg"
              "\n  V_BUS_load2  : |V|={} V RMS, angle={} deg"
              "\n  V_BUS_load2b : |V|={} V RMS, angle={} deg"
              "\n  V_BUS_pv     : |V|={} V RMS, angle={} deg"
              "\n  GFL_pv       : P={} MW, Q={} Mvar",
              gasGeneratorPF->getApparentPower().real() / 1e6,
              gasGeneratorPF->getApparentPower().imag() / 1e6,
              Math::abs(busGasPF->singleVoltage()),
              Math::phase(busGasPF->singleVoltage()) * 180.0 / PI,
              Math::abs(busBPF->singleVoltage()),
              Math::phase(busBPF->singleVoltage()) * 180.0 / PI,
              Math::abs(busLoadPF->singleVoltage()),
              Math::phase(busLoadPF->singleVoltage()) * 180.0 / PI,
              Math::abs(busLoad2bPF->singleVoltage()),
              Math::phase(busLoad2bPF->singleVoltage()) * 180.0 / PI,
              Math::abs(busPvPF->singleVoltage()),
              Math::phase(busPvPF->singleVoltage()) * 180.0 / PI,
              p.gflInitialActivePower / 1e6, p.gflInitialReactivePower / 1e6);

  // ==========================================================================
  // 2. EMT
  // ==========================================================================

  const String emtSimName = p.name + "_EMT";
  Logger::setLogDir("logs/" + emtSimName);

  // EMT nodes: exact same names as PF.
  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busLoad = SimNode<Real>::make("BUS_load_2", PhaseType::ABC);
  auto busLoad2b = SimNode<Real>::make("BUS_load_2b", PhaseType::ABC);
  auto busPv = SimNode<Real>::make("BUS_pv", PhaseType::ABC);

  // ---------------------------------------------------------------------------
  // EMT GEN_gas
  // ---------------------------------------------------------------------------
  auto gasGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_gas", Logger::Level::off);

  gasGenerator->setBaseAndOperationalPerUnitParameters(
      p.genRatedPower, p.genRatedVoltage, p.genNominalFrequency,
      p.genPoleNumber, p.genNominalFieldCurrent, p.genStatorResistance, p.genLd,
      p.genLq, p.genLdTransient, p.genLqTransient, p.genLdSubtransient,
      p.genLqSubtransient, p.genLeakageInductance, p.genTd0Transient,
      p.genTq0Transient, p.genTd0Subtransient, p.genTq0Subtransient,
      p.genInertia);

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

  // EMT loads intentionally have no setParameters() call.
  // initWithPowerflow() transfers their PF terminal P/Q and RXLoad derives
  // R/X from that operating point and the solved PF bus voltage.
  auto load2a = EMT::Ph3::RXLoad::make("LOAD_2a", Logger::Level::off);
  auto load2b = EMT::Ph3::RXLoad::make("LOAD_2b", Logger::Level::off);

  // EMT breaker: same position and same initial CLOSED state as PF.
  auto breaker = EMT::Ph3::Switch::make("BR_load_2b", Logger::Level::off);
  breaker->setParameters(
      Math::singlePhaseParameterToThreePhase(p.breakerOpenResistance),
      Math::singlePhaseParameterToThreePhase(p.breakerClosedResistance));
  breaker->closeSwitch();

  // ---------------------------------------------------------------------------
  // EMT TR_pv + GFL
  // ---------------------------------------------------------------------------
  auto trPv =
      EMT::Ph3::Transformer::make("TR_pv", "TR_pv", Logger::Level::off, false);
  trPv->setParameters(p.gflNetworkVoltage, p.gflRatedVoltage, p.gflRatedPower,
                      p.gflNetworkVoltage / p.gflRatedVoltage, 0.0,
                      Math::singlePhaseParameterToThreePhase(p.trPvResistance),
                      Math::singlePhaseParameterToThreePhase(p.trPvInductance));

  auto pvGfl = EMT::Ph3::GFL_Siemens::make("pv", "pv", Logger::Level::off);

  const Real gflInitialVoltagePu =
      Math::abs(busPvPF->singleVoltage()) / p.gflRatedVoltage;

  pvGfl->setBaseParameters(p.gflRatedPower, p.gflRatedVoltage, p.frequency);

  pvGfl->setReferencesPerUnit(1.0, gflInitialVoltagePu,
                              p.gflInitialActivePower / p.gflRatedPower,
                              p.gflInitialReactivePower / p.gflRatedPower);

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

  pvGfl->setFilterParametersPerUnit(p.gflFilterInductiveReactancePu(),
                                    p.gflFilterCapacitiveSusceptancePu(),
                                    p.gflFilterResistancePu());

  pvGfl->withControl(true);

  // EMT connectivity: exactly the same graph and component names as PF.
  gasGenerator->connect({busGas});
  trGas->connect({busGas, busB});
  trLoad->connect({busLoad, busB});
  load2a->connect({busLoad});
  breaker->connect({busLoad, busLoad2b});
  load2b->connect({busLoad2b});
  trPv->connect({busB, busPv});
  pvGfl->connect({busPv});

  SystemTopology systemEMT(
      p.frequency, SystemNodeList{busGas, busB, busLoad, busLoad2b, busPv},
      SystemComponentList{gasGenerator, trGas, trLoad, load2a, breaker, load2b,
                          trPv, pvGfl});

  // Central PF -> EMT operating-point transfer:
  //   - solved node voltages by exact node name
  //   - synchronous-generator terminal P/Q with generator sign conversion
  //   - load terminal P/Q with consumer-positive sign convention
  //
  // For the GFL, the PF surrogate and EMT model are both named "pv".
  // Its initial P/Q are zero, so the generic SP-load transfer is
  // sign-convention neutral.
  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  auto loggerEMT = DataLogger::make(emtSimName, true, p.logDownSampling);

  loggerEMT->logAttribute("BUS_gas_v", busGas->attribute("v"));
  loggerEMT->logAttribute("BUS_b_v", busB->attribute("v"));
  loggerEMT->logAttribute("BUS_load_2_v", busLoad->attribute("v"));
  loggerEMT->logAttribute("BUS_load_2b_v", busLoad2b->attribute("v"));
  loggerEMT->logAttribute("BUS_pv_v", busPv->attribute("v"));

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

  loggerEMT->logAttribute("TR_pv_i", trPv->attribute("i_intf"));
  loggerEMT->logAttribute("GFL_pv_P_ref_pu", pvGfl->attribute("P_ref_pu"));
  loggerEMT->logAttribute("GFL_pv_Q_ref_pu", pvGfl->attribute("Q_ref_pu"));
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
  loggerEMT->logAttribute("GFL_pv_v_pcc_dq_pu",
                          pvGfl->attribute("v_pcc_dq_pu"));

  Simulation simulationEMT(emtSimName, Logger::Level::info);
  simulationEMT.setSystem(systemEMT);
  simulationEMT.setTimeStep(p.timeStep);
  simulationEMT.setFinalTime(p.finalTime);
  simulationEMT.setDomain(Domain::EMT);

  // All dynamic states are derived from the PF operating point.
  // Matrix-recomputation mode is intentionally left at Auto. The breaker
  // supports precomputed open/closed matrices, and the GFL RLC network stamp
  // is fixed while its controller acts through the source vector.
  simulationEMT.doInitFromNodesAndTerminals(true);
  simulationEMT.addLogger(loggerEMT);

  // Open the load breaker -> disconnect LOAD_2b.
  simulationEMT.addEvent(
      SwitchEvent3Ph::make(p.breakerOpenTime, breaker, false));

  // Change the GFL active-power setpoint after the load-loss event.
  // The GFL uses per-unit references, therefore -1 MW on a 50 MVA base is
  // -0.02 pu.
  const Real gflActivePowerStepTime =
      std::round(p.gflActivePowerStepTime / p.timeStep) * p.timeStep;
  const Real gflActivePowerAfterStepPu =
      p.gflActivePowerAfterStep / p.gflRatedPower;

  simulationEMT.addEvent(AttributeEvent<Real>::make(gflActivePowerStepTime,
                                                    pvGfl->mActivePowerRefPu,
                                                    gflActivePowerAfterStepPu));

  SPDLOG_INFO("\nStarting EMT SynGen + GFL load-loss simulation:"
              "\n  dt={} s"
              "\n  finalTime={} s"
              "\n  BR_load_2b opens at t={} s"
              "\n  disconnected load: P={} MW, Q={} Mvar"
              "\n  GFL initial P={} MW, Q={} Mvar"
              "\n  GFL P reference changes to {} MW ({} pu) at t={} s"
              "\n  GFL droop: K_fP={} pu, K_VQ={} pu",
              p.timeStep, p.finalTime, p.breakerOpenTime,
              p.load2bActivePower / 1e6, p.load2bReactivePower / 1e6,
              p.gflInitialActivePower / 1e6, p.gflInitialReactivePower / 1e6,
              p.gflActivePowerAfterStep / 1e6, gflActivePowerAfterStepPu,
              gflActivePowerStepTime, p.gflFrequencyToActivePowerGainPu,
              p.gflVoltageToReactivePowerGainPu);

  simulationEMT.run();

  return 0;
}
