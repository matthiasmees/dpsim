// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_Ph3_GFL_Siemens.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorVBR.h>
#include <dpsim-models/Signal/Exciter.h>

#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;

namespace ScenarioASynGenVBRSynGenVBRGflSiemens {

// =============================================================================
// Parameters
// =============================================================================

struct SimulationParameters {
  String name = "EMT_Scenario_A_SynGenVBR_GFL_SynGenVBR";
  Real frequency = 50.0;
  Real timeStep = 100e-6;
  Real finalTime = 15.0;

  // Both generator islands are initialized with the breaker open. The breaker
  // is closed during the EMT simulation.
  Real breakerCloseTime = 10.0;

  // The PV GFL starts from the zero-power PF operating point and receives the
  // same reactive-power reference step as in the Siemens GFL scenario.
  Real gflReactivePowerStepTime = 5.0;
  Real gflReactivePowerAfterStep = -25e6;

  // Diagnostic event switches. Defaults reproduce the requested sequence:
  // GFL Q step first, PSH breaker closing second.
  Bool enableGflReactivePowerStep = true;
  Bool enableBreakerClosing = true;

  UInt logDownsampling = 10;
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
    p.tq0Transient = 1.0;
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

struct GflValidationParameters {
  // Siemens GFL validation/controller parameters copied from the working
  // Scenario A GFM/GFL example.
  Real frequencyToActivePowerGainPu = 4.0;
  Real voltageToReactivePowerGainPu = 2.0;

  Real pllKp = 0.449 / (9.1e-3 + 250e-6);
  Real pllKi = 250e-6 / (5.9 * (9.1e-3 + 250e-6));
  Real pllGainScale = 1.0;

  Real currentControllerKp = 0.74 / 10.0;
  Real currentControllerKi = 11.9 / 10.0;
  Real currentControllerGainScale = 1.0;
  Real pccVoltageFeedforwardGain = 1.0;

  Real activePowerMeasurementTimeConstant = 0.05;
  Real reactivePowerMeasurementTimeConstant = 0.05;
  Real currentMeasurementTimeConstant = 0.5e-3;

  Real minimumFrequencyPu = 0.80;
  Real maximumFrequencyPu = 1.20;
  Real maximumCurrentReferencePu = 1.50;
  Real maximumVoltageCommandPu = 1.20;
  Real minimumVoltageForCurrentReferencePu = 0.05;
};

struct SiemensValidationFilter {
  // Physical Siemens validation filter: 50 kVA / 257 V. Its R, X_L and B_C
  // are transferred to the 50 MVA / 11 kV GFL base in per unit.
  Real referencePower = 50e3;
  Real referenceVoltage = 257.0;
  Real referenceFrequency = 50.0;

  Real resistance = 0.0279;
  Real inductance = 2.72e-4;
  Real capacitance = 3.5e-6;

  Real baseImpedance() const {
    return referenceVoltage * referenceVoltage / referencePower;
  }

  Real baseOmega() const { return 2.0 * PI * referenceFrequency; }

  Real resistancePu() const { return resistance / baseImpedance(); }

  Real inductiveReactancePu() const {
    return baseOmega() * inductance / baseImpedance();
  }

  Real capacitiveSusceptancePu() const {
    return baseOmega() * capacitance * baseImpedance();
  }
};

struct PvGflParameters {
  Real ratedPower = 50e6;
  Real ratedVoltage = 11.0e3;
  Real networkVoltage = 220e3;

  // Transformer orientation is network side -> converter side, matching the
  // former transformer embedded in AvVoltageSourceInverterDQ.
  Real transformerResistance = 0.0;
  Real transformerInductance = 0.928e-3;

  Real initialActivePower = 0.0;
  Real initialReactivePower = 0.0;
};

struct DiagnosticParameters {
  // Controller switches isolate whether the synchronous-machine outer
  // controls participate in the unstable mode.
  Bool enableGasGovernor = true;
  Bool enableGasExciter = true;
  Bool enablePshGovernor = true;
  Bool enablePshExciter = true;
  Bool enableGflControl = true;

  // Selected with --option DIAGNOSTIC_CASE=<0..12>.
  Int caseId = 0;
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

  PvGflParameters pvGfl;
  GflValidationParameters gfl;
  SiemensValidationFilter gflFilter;
  DiagnosticParameters diagnostics;

  // The initial angle is applied to the isolated PSH island in the SP
  // initialization and can be changed to study out-of-phase closing.
  Real pshInitialAngleDegrees = 0.0;
};

// =============================================================================
// Diagnostics
// =============================================================================

void applyDiagnosticCase(Parameters &p, Int caseId) {
  p.diagnostics.caseId = caseId;

  switch (caseId) {
  case 0:
    // Full requested case: Q step at 5 s, breaker closing at 10 s.
    break;
  case 1:
    // Initialization only: no external events.
    p.simulation.enableGflReactivePowerStep = false;
    p.simulation.enableBreakerClosing = false;
    break;
  case 2:
    // Isolate the GFL Q-step response on the gas-machine island.
    p.simulation.enableGflReactivePowerStep = true;
    p.simulation.enableBreakerClosing = false;
    break;
  case 3:
    // Isolate the breaker-closing response.
    p.simulation.enableGflReactivePowerStep = false;
    p.simulation.enableBreakerClosing = true;
    break;
  case 4:
    // Small Q step on the gas-machine island.
    p.simulation.enableGflReactivePowerStep = true;
    p.simulation.enableBreakerClosing = false;
    p.simulation.gflReactivePowerAfterStep = -5e6;
    break;
  case 5:
    // Direct P/Q control: remove inverse-droop feedback.
    p.simulation.enableGflReactivePowerStep = true;
    p.simulation.enableBreakerClosing = false;
    p.gfl.frequencyToActivePowerGainPu = 0.0;
    p.gfl.voltageToReactivePowerGainPu = 0.0;
    break;
  case 6:
    // Reduce PLL bandwidth while keeping the original Q outer loop.
    p.simulation.enableGflReactivePowerStep = true;
    p.simulation.enableBreakerClosing = false;
    p.gfl.pllGainScale = 0.25;
    break;
  case 7:
    // Hold gas-machine mechanical torque and field voltage fixed.
    p.simulation.enableGflReactivePowerStep = true;
    p.simulation.enableBreakerClosing = false;
    p.diagnostics.enableGasGovernor = false;
    p.diagnostics.enableGasExciter = false;
    break;
  case 8: {
    // Artificial diagnostic only: approximate the GFM filter reactance with
    // the synchronous-machine subtransient reactances. Xd'' cannot be reduced
    // to this value with the original Xl=0.135 pu, so Xl is reduced jointly.
    const Real targetSubtransientReactance = p.gflFilter.inductiveReactancePu();
    p.simulation.enableGflReactivePowerStep = true;
    p.simulation.enableBreakerClosing = false;
    p.gasGenerator.leakageInductance = 0.5 * targetSubtransientReactance;
    p.gasGenerator.ldSubtransient = targetSubtransientReactance;
    p.gasGenerator.lqSubtransient = targetSubtransientReactance;
    break;
  }
  case 9:
    // Initialization only with direct P/Q control.
    p.simulation.enableGflReactivePowerStep = false;
    p.simulation.enableBreakerClosing = false;
    p.gfl.frequencyToActivePowerGainPu = 0.0;
    p.gfl.voltageToReactivePowerGainPu = 0.0;
    break;
  case 10:
    // Initialization only with direct P/Q and no current-sensor lag.
    p.simulation.enableGflReactivePowerStep = false;
    p.simulation.enableBreakerClosing = false;
    p.gfl.frequencyToActivePowerGainPu = 0.0;
    p.gfl.voltageToReactivePowerGainPu = 0.0;
    p.gfl.currentMeasurementTimeConstant = 0.0;
    break;
  case 11:
    // Initialization only with reduced PLL gains.
    p.simulation.enableGflReactivePowerStep = false;
    p.simulation.enableBreakerClosing = false;
    p.gfl.pllGainScale = 0.25;
    break;
  case 12:
    // Initialization only with the GFL current controller bypassed.
    // The initialized source-voltage dq command is retained and rotated by
    // the PLL/open-loop path.
    p.simulation.enableGflReactivePowerStep = false;
    p.simulation.enableBreakerClosing = false;
    p.diagnostics.enableGflControl = false;
    break;
  default:
    throw std::invalid_argument(
        "DIAGNOSTIC_CASE must be an integer from 0 to 12");
  }
}

Bool boolFromOption(Real value, const String &name) {
  if (value == 0.0)
    return false;
  if (value == 1.0)
    return true;
  throw std::invalid_argument(name + " must be either 0 or 1");
}

void validateMachineParameters(const SynGenVbrParameters &m,
                               const String &name) {
  const Bool dAxisValid = m.ld > m.ldTransient &&
                          m.ldTransient > m.ldSubtransient &&
                          m.ldSubtransient > m.leakageInductance;
  const Bool qAxisValid = m.lq > m.lqTransient &&
                          m.lqTransient > m.lqSubtransient &&
                          m.lqSubtransient > m.leakageInductance;

  if (!dAxisValid || !qAxisValid || !(m.leakageInductance > 0.0)) {
    throw std::invalid_argument(name + " requires Xd > Xd' > Xd'' > Xl and "
                                       "Xq > Xq' > Xq'' > Xl");
  }
}

void validateParameters(const Parameters &p) {
  if (!(p.simulation.frequency > 0.0) || !(p.simulation.timeStep > 0.0) ||
      !(p.simulation.finalTime > 0.0) || !(p.simulation.logDownsampling > 0)) {
    throw std::invalid_argument(
        "Require frequency>0, dt>0, finalTime>0 and LOG_DOWNSAMPLING>0");
  }

  if (p.simulation.enableGflReactivePowerStep &&
      (!(p.simulation.gflReactivePowerStepTime >= 0.0) ||
       !(p.simulation.gflReactivePowerStepTime < p.simulation.finalTime))) {
    throw std::invalid_argument(
        "Enabled GFL Q step must satisfy 0<=PV_Q_STEP_TIME<finalTime");
  }

  if (p.simulation.enableBreakerClosing &&
      (!(p.simulation.breakerCloseTime >= 0.0) ||
       !(p.simulation.breakerCloseTime < p.simulation.finalTime))) {
    throw std::invalid_argument(
        "Enabled breaker event must satisfy 0<=BREAKER_CLOSE_TIME<finalTime");
  }

  if (p.simulation.enableGflReactivePowerStep &&
      p.simulation.enableBreakerClosing &&
      !(p.simulation.gflReactivePowerStepTime <
        p.simulation.breakerCloseTime)) {
    throw std::invalid_argument("This diagnostic example requires the GFL Q "
                                "step before breaker closing");
  }

  if (!std::isfinite(p.simulation.gflReactivePowerAfterStep) ||
      !(p.pvGfl.ratedPower > 0.0) || !(p.pvGfl.ratedVoltage > 0.0) ||
      !(p.gflFilter.inductiveReactancePu() > 0.0) ||
      !(p.gflFilter.capacitiveSusceptancePu() > 0.0) ||
      !(p.gflFilter.resistancePu() >= 0.0) || !(p.gfl.pllGainScale > 0.0) ||
      !(p.gfl.currentControllerGainScale > 0.0) ||
      !(p.gfl.currentMeasurementTimeConstant >= 0.0)) {
    throw std::invalid_argument(
        "Invalid GFL ratings, filter, references, gain scales, or filter time");
  }

  validateMachineParameters(p.gasGenerator, "GEN_gas");
  validateMachineParameters(p.pshGenerator, "GEN_psh");
}

void logDiagnosticSummary(const Parameters &p) {
  const Real omega = 2.0 * PI * p.simulation.frequency;
  const Real zBaseHv = p.gasTransformer.nominalVoltageHigh *
                       p.gasTransformer.nominalVoltageHigh /
                       p.gasGenerator.ratedPower;

  const Real transformerReactancePu =
      omega * p.gasTransformer.inductance(p.simulation.frequency) / zBaseHv;
  const Real cableReactancePu = omega * p.cable.inductance / zBaseHv;

  const Real synGenPathReactancePu =
      p.gasGenerator.ldSubtransient + transformerReactancePu + cableReactancePu;
  const Real gfmReferencePathReactancePu = p.gflFilter.inductiveReactancePu() +
                                           transformerReactancePu +
                                           cableReactancePu;

  const Real synGenScr =
      synGenPathReactancePu > 0.0 ? 1.0 / synGenPathReactancePu : 0.0;
  const Real gfmReferenceScr = gfmReferencePathReactancePu > 0.0
                                   ? 1.0 / gfmReferencePathReactancePu
                                   : 0.0;

  SPDLOG_INFO(
      "Diagnostic configuration:"
      "\n  case={} dt={} s, Q-step enabled={}, breaker enabled={}"
      "\n  Q step: {} var at {} s; breaker close: {} s"
      "\n  gas controls: governor={}, exciter={}"
      "\n  GFL: control={}, KfP={}, KVQ={}, PLL scale={}, "
      "current-PI scale={}, tau_i={} s"
      "\n  gas machine: Xl={}, Xd''={}, Xq''={}"
      "\n  path estimate on 50 MVA base:"
      " Xtr={} pu, Xcable={} pu, Xsyn={} pu (SCR~{}),"
      " Xgfm-reference={} pu (SCR~{})"
      "\n  exact Xd'' target for equal scalar path reactance={} pu",
      p.diagnostics.caseId, p.simulation.timeStep,
      p.simulation.enableGflReactivePowerStep,
      p.simulation.enableBreakerClosing, p.simulation.gflReactivePowerAfterStep,
      p.simulation.gflReactivePowerStepTime, p.simulation.breakerCloseTime,
      p.diagnostics.enableGasGovernor, p.diagnostics.enableGasExciter,
      p.diagnostics.enableGflControl, p.gfl.frequencyToActivePowerGainPu,
      p.gfl.voltageToReactivePowerGainPu, p.gfl.pllGainScale,
      p.gfl.currentControllerGainScale, p.gfl.currentMeasurementTimeConstant,
      p.gasGenerator.leakageInductance, p.gasGenerator.ldSubtransient,
      p.gasGenerator.lqSubtransient, transformerReactancePu, cableReactancePu,
      synGenPathReactancePu, synGenScr, gfmReferencePathReactancePu,
      gfmReferenceScr, p.gflFilter.inductiveReactancePu());

  if (p.gflFilter.inductiveReactancePu() <= p.gasGenerator.leakageInductance) {
    SPDLOG_WARN(
        "Matching only GEN_gas Xd'' to the GFM filter value is invalid: "
        "target Xd''={} pu is not greater than Xl={} pu. Use diagnostic "
        "case 8 only as an artificial joint Xl/Xd''/Xq'' sensitivity test.",
        p.gflFilter.inductiveReactancePu(), p.gasGenerator.leakageInductance);
  }
}

// =============================================================================
// Power-flow / steady-state initialization
// =============================================================================

struct PowerFlowResult {
  SystemTopology system;

  Complex gasPower;
  Complex pshPower;
  Complex pvPower;

  Complex gasBusVoltage;
  Complex pshBusVoltage;
  Complex pvBusVoltage;

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
  auto busPv = SimNode<Complex>::make("BUS_pv", PhaseType::Single);

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

  auto pvTransformer = SP::Ph1::Transformer::make("TR_pv_PF", "TR_pv_PF",
                                                  Logger::Level::debug, false);
  pvTransformer->setParameters(
      p.pvGfl.networkVoltage, p.pvGfl.ratedVoltage, p.pvGfl.ratedPower,
      p.pvGfl.networkVoltage / p.pvGfl.ratedVoltage, 0.0,
      p.pvGfl.transformerResistance, p.pvGfl.transformerInductance);

  // Zero-power PF surrogate for the initial GFL terminal operating point.
  auto pvEquivalent = SP::Ph1::Load::make("pv_PF", Logger::Level::debug);
  pvEquivalent->setParameters(p.pvGfl.initialActivePower,
                              p.pvGfl.initialReactivePower,
                              p.pvGfl.ratedVoltage);

  gasInjection->connect({busGas});
  gasTransformer->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPsha});
  breaker->connect({busPsha, busTrPsh});
  pshTransformer->connect({busPsh, busTrPsh});
  pshInjection->connect({busPsh});

  pvTransformer->connect({busA, busPv});
  pvEquivalent->connect({busPv});

  SystemTopology system(
      p.simulation.frequency,
      SystemNodeList{busGas, busB, busA, busPsha, busTrPsh, busPsh, busPv},
      SystemComponentList{gasInjection, gasTransformer, cable, line3, breaker,
                          pshTransformer, pshInjection, pvTransformer,
                          pvEquivalent});

  auto logger = DataLogger::make(simulationName);

  logger->logAttribute("V_BUS_gas_PF", busGas->attribute("v"));
  logger->logAttribute("V_BUS_b_PF", busB->attribute("v"));
  logger->logAttribute("V_BUS_a_PF", busA->attribute("v"));
  logger->logAttribute("V_BUS_psha_PF", busPsha->attribute("v"));
  logger->logAttribute("V_BUS_tr_psh_PF", busTrPsh->attribute("v"));
  logger->logAttribute("V_BUS_psh_PF", busPsh->attribute("v"));
  logger->logAttribute("V_BUS_pv_PF", busPv->attribute("v"));

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
  const Complex pvPower(p.pvGfl.initialActivePower,
                        p.pvGfl.initialReactivePower);

  const Complex gasVoltage = busGas->singleVoltage();
  const Complex pshVoltage = busPsh->singleVoltage();
  const Complex pvVoltage = busPv->singleVoltage();
  const Complex breakerGridVoltage = busPsha->singleVoltage();
  const Complex breakerPshVoltage = busTrPsh->singleVoltage();

  const auto finiteComplex = [](const Complex &value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  };

  if (!(std::abs(gasVoltage) > 0.0) || !(std::abs(pshVoltage) > 0.0) ||
      !(std::abs(pvVoltage) > 0.0) || !finiteComplex(gasVoltage) ||
      !finiteComplex(pshVoltage) || !finiteComplex(pvVoltage) ||
      !finiteComplex(gasPower) || !finiteComplex(pshPower) ||
      !finiteComplex(pvPower) || !finiteComplex(breakerGridVoltage) ||
      !finiteComplex(breakerPshVoltage)) {
    throw std::runtime_error(
        "Invalid SP steady-state result for Scenario A with two "
        "SynchronGeneratorVBR components.");
  }

  SPDLOG_INFO(
      "Solved open-breaker SP initialization:"
      "\n  GEN_gas: P={} W, Q={} var, |V|={} V_LL RMS, angle={} deg"
      "\n  GEN_psh: P={} W, Q={} var, |V|={} V_LL RMS, angle={} deg"
      "\n  PV GFL : P={} W, Q={} var, |V|={} V_LL RMS, angle={} deg"
      "\n  breaker grid side: |V|={} V_LL RMS, angle={} deg"
      "\n  breaker PSH side : |V|={} V_LL RMS, angle={} deg"
      "\n  breaker mismatch : |dV|={} V_LL RMS, dAngle={} deg",
      gasPower.real(), gasPower.imag(), std::abs(gasVoltage),
      std::arg(gasVoltage) * 180.0 / PI, pshPower.real(), pshPower.imag(),
      std::abs(pshVoltage), std::arg(pshVoltage) * 180.0 / PI, pvPower.real(),
      pvPower.imag(), std::abs(pvVoltage), std::arg(pvVoltage) * 180.0 / PI,
      std::abs(breakerGridVoltage), std::arg(breakerGridVoltage) * 180.0 / PI,
      std::abs(breakerPshVoltage), std::arg(breakerPshVoltage) * 180.0 / PI,
      std::abs(breakerGridVoltage - breakerPshVoltage),
      std::arg(breakerGridVoltage / breakerPshVoltage) * 180.0 / PI);

  return {system,    gasPower,           pshPower,
          pvPower,   gasVoltage,         pshVoltage,
          pvVoltage, breakerGridVoltage, breakerPshVoltage};
}

// =============================================================================
// EMT model
// =============================================================================

void configureGenerator(
    const std::shared_ptr<EMT::Ph3::SynchronGeneratorVBR> &generator,
    const SynGenVbrParameters &parameters, const Complex &initialGeneratorPower,
    const String &exciterName, Bool enableGovernor, Bool enableExciter) {
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

  if (enableGovernor) {
    generator->addGovernor(
        parameters.governorTa, parameters.governorTb, parameters.governorTc,
        parameters.governorFa, parameters.governorFb, parameters.governorFc,
        parameters.governorGain, parameters.governorSpeedRelayTimeConstant,
        parameters.governorServoMotorTimeConstant, mechanicalPowerPerUnit,
        mechanicalPowerPerUnit);
  }

  if (enableExciter) {
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
  }

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

void configureGfl(const std::shared_ptr<EMT::Ph3::GFL_Siemens> &gfl,
                  const Parameters &p, const Complex &initialPower,
                  const Complex &initialVoltage) {
  const Real initialVoltagePu = std::abs(initialVoltage) / p.pvGfl.ratedVoltage;
  const Real initialActivePowerPu = initialPower.real() / p.pvGfl.ratedPower;
  const Real initialReactivePowerPu = initialPower.imag() / p.pvGfl.ratedPower;

  gfl->setBaseParameters(p.pvGfl.ratedPower, p.pvGfl.ratedVoltage,
                         p.simulation.frequency);
  gfl->setReferencesPerUnit(1.0, initialVoltagePu, initialActivePowerPu,
                            initialReactivePowerPu);

  gfl->setDroopParametersPerUnit(p.gfl.frequencyToActivePowerGainPu,
                                 p.gfl.voltageToReactivePowerGainPu);
  gfl->setPllParameters(p.gfl.pllGainScale * p.gfl.pllKp,
                        p.gfl.pllGainScale * p.gfl.pllKi);
  gfl->setCurrentControllerParameters(
      p.gfl.currentControllerGainScale * p.gfl.currentControllerKp,
      p.gfl.currentControllerGainScale * p.gfl.currentControllerKi,
      p.gfl.pccVoltageFeedforwardGain);
  gfl->setMeasurementFilterTimeConstants(
      p.gfl.activePowerMeasurementTimeConstant,
      p.gfl.reactivePowerMeasurementTimeConstant,
      p.gfl.currentMeasurementTimeConstant);
  gfl->setControllerLimitsPerUnit(
      p.gfl.minimumFrequencyPu, p.gfl.maximumFrequencyPu,
      p.gfl.maximumCurrentReferencePu, p.gfl.maximumVoltageCommandPu,
      p.gfl.minimumVoltageForCurrentReferencePu);

  gfl->setFilterParametersPerUnit(p.gflFilter.inductiveReactancePu(),
                                  p.gflFilter.capacitiveSusceptancePu(),
                                  p.gflFilter.resistancePu());
  gfl->withControl(p.diagnostics.enableGflControl);

  SPDLOG_INFO(
      "{} parameters: S_base={} VA, V_base_LL={} V, P_init={} pu, "
      "Q_init={} pu, V_init={} pu, Rf={} pu, X_Lf={} pu, B_Cf={} pu, "
      "control={}, KfP={}, KVQ={}, PLL scale={}, current PI scale={}",
      gfl->name(), p.pvGfl.ratedPower, p.pvGfl.ratedVoltage,
      initialActivePowerPu, initialReactivePowerPu, initialVoltagePu,
      p.gflFilter.resistancePu(), p.gflFilter.inductiveReactancePu(),
      p.gflFilter.capacitiveSusceptancePu(), p.diagnostics.enableGflControl,
      p.gfl.frequencyToActivePowerGainPu, p.gfl.voltageToReactivePowerGainPu,
      p.gfl.pllGainScale, p.gfl.currentControllerGainScale);
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
  auto busPv = SimNode<Real>::make("BUS_pv", PhaseType::ABC);

  auto gasGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_gas", Logger::Level::debug);
  auto pshGenerator =
      EMT::Ph3::SynchronGeneratorVBR::make("GEN_psh", Logger::Level::debug);

  configureGenerator(gasGenerator, p.gasGenerator, powerFlow.gasPower,
                     "GEN_gas_exciter", p.diagnostics.enableGasGovernor,
                     p.diagnostics.enableGasExciter);
  configureGenerator(pshGenerator, p.pshGenerator, powerFlow.pshPower,
                     "GEN_psh_exciter", p.diagnostics.enablePshGovernor,
                     p.diagnostics.enablePshExciter);

  auto pvGfl = EMT::Ph3::GFL_Siemens::make("pv", "pv", Logger::Level::debug);
  configureGfl(pvGfl, p, powerFlow.pvPower, powerFlow.pvBusVoltage);

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

  // Keep the exact transformer orientation and zero-resistance construction
  // from the working Siemens GFL scenario.
  auto pvTransformer = EMT::Ph3::Transformer::make("TR_pv", "TR_pv",
                                                   Logger::Level::debug, false);
  pvTransformer->setParameters(
      p.pvGfl.networkVoltage, p.pvGfl.ratedVoltage, p.pvGfl.ratedPower,
      p.pvGfl.networkVoltage / p.pvGfl.ratedVoltage, 0.0,
      Math::singlePhaseParameterToThreePhase(p.pvGfl.transformerResistance),
      Math::singlePhaseParameterToThreePhase(p.pvGfl.transformerInductance));

  gasGenerator->connect({busGas});
  gasTransformer->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPsha});
  breaker->connect({busPsha, busTrPsh});
  pshTransformer->connect({busPsh, busTrPsh});
  pshGenerator->connect({busPsh});

  pvTransformer->connect({busA, busPv});
  pvGfl->connect({busPv});

  SystemTopology system(
      p.simulation.frequency,
      SystemNodeList{busGas, busB, busA, busPsha, busTrPsh, busPsh, busPv},
      SystemComponentList{gasGenerator, gasTransformer, cable, line3, breaker,
                          pshTransformer, pshGenerator, pvTransformer, pvGfl});

  // The open-breaker SP topology contains both islands and identically named
  // nodes. This initializes both synchronous machines and both breaker sides.
  system.initWithPowerflow(powerFlow.system, Domain::EMT);

  // SynchronGeneratorVBR expects terminal power in DPsim's consumer-positive
  // convention and negates it during initializeFromNodesAndTerminals().
  gasGenerator->terminal(0)->setPower(-powerFlow.gasPower);
  pshGenerator->terminal(0)->setPower(-powerFlow.pshPower);

  // GFL_Siemens uses generator-positive P/Q references internally.
  pvGfl->terminal(0)->setPower(powerFlow.pvPower);

  SPDLOG_INFO(
      "EMT initialization transfer:"
      "\n  GEN_gas terminal power (consumer convention)=[{}, {}] VA"
      "\n  GEN_psh terminal power (consumer convention)=[{}, {}] VA"
      "\n  GFL terminal power (generator convention used by GFL)=[{}, {}] VA"
      "\n  BUS_gas initial phasor: |V|={} V_LL RMS, angle={} deg"
      "\n  BUS_psh initial phasor: |V|={} V_LL RMS, angle={} deg"
      "\n  BUS_pv initial phasor : |V|={} V_LL RMS, angle={} deg",
      -powerFlow.gasPower.real(), -powerFlow.gasPower.imag(),
      -powerFlow.pshPower.real(), -powerFlow.pshPower.imag(),
      powerFlow.pvPower.real(), powerFlow.pvPower.imag(),
      std::abs(busGas->initialSingleVoltage()),
      std::arg(busGas->initialSingleVoltage()) * 180.0 / PI,
      std::abs(busPsh->initialSingleVoltage()),
      std::arg(busPsh->initialSingleVoltage()) * 180.0 / PI,
      std::abs(busPv->initialSingleVoltage()),
      std::arg(busPv->initialSingleVoltage()) * 180.0 / PI);

  auto logger =
      DataLogger::make(simulationName, true, p.simulation.logDownsampling);

  logger->logAttribute("BUS_gas_v", busGas->attribute("v"));
  logger->logAttribute("BUS_b_v", busB->attribute("v"));
  logger->logAttribute("BUS_a_v", busA->attribute("v"));
  logger->logAttribute("BUS_psha_v", busPsha->attribute("v"));
  logger->logAttribute("BUS_tr_psh_v", busTrPsh->attribute("v"));
  logger->logAttribute("BUS_psh_v", busPsh->attribute("v"));
  logger->logAttribute("BUS_pv_v", busPv->attribute("v"));

  logger->logAttribute("TR_gas_i", gasTransformer->attribute("i_intf"));
  logger->logAttribute("cable_i", cable->attribute("i_intf"));
  logger->logAttribute("line_3_i", line3->attribute("i_intf"));
  logger->logAttribute("breaker_GEN_psh_i", breaker->attribute("i_intf"));
  logger->logAttribute("breaker_GEN_psh_v", breaker->attribute("v_intf"));
  logger->logAttribute("TR_psh_i", pshTransformer->attribute("i_intf"));
  logger->logAttribute("TR_pv_i", pvTransformer->attribute("i_intf"));

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

  logger->logAttribute("GFL_pv_P_ref_pu", pvGfl->attribute("P_ref_pu"));
  logger->logAttribute("GFL_pv_Q_ref_pu", pvGfl->attribute("Q_ref_pu"));
  logger->logAttribute("GFL_pv_P_command_pu", pvGfl->attribute("P_command_pu"));
  logger->logAttribute("GFL_pv_Q_command_pu", pvGfl->attribute("Q_command_pu"));
  logger->logAttribute("GFL_pv_P_elec_pu", pvGfl->attribute("P_elec_pu"));
  logger->logAttribute("GFL_pv_Q_elec_pu", pvGfl->attribute("Q_elec_pu"));
  logger->logAttribute("GFL_pv_frequency_pu", pvGfl->attribute("frequency_pu"));
  logger->logAttribute("GFL_pv_V_magnitude_pu",
                       pvGfl->attribute("V_magnitude_pu"));
  logger->logAttribute("GFL_pv_pll_vq_error_pu",
                       pvGfl->attribute("pll_vq_error_pu"));
  logger->logAttribute("GFL_pv_i_pcc_dq_pu", pvGfl->attribute("i_pcc_dq_pu"));
  logger->logAttribute("GFL_pv_i_ref_dq_pu", pvGfl->attribute("i_ref_dq_pu"));
  logger->logAttribute("GFL_pv_v_cmd_dq_pu", pvGfl->attribute("v_cmd_dq_pu"));
  logger->logAttribute("GFL_pv_v_pcc_dq_pu", pvGfl->attribute("v_pcc_dq_pu"));
  logger->logAttribute("GFL_pv_i_lf_dq_pu", pvGfl->attribute("i_lf_dq_pu"));
  logger->logAttribute("GFL_pv_i_pcc_filtered_dq_pu",
                       pvGfl->attribute("i_pcc_filtered_dq_pu"));
  logger->logAttribute("GFL_pv_i_error_dq_pu", pvGfl->attribute("e_i_dq_pu"));
  logger->logAttribute("GFL_pv_current_integrator_dq_pu",
                       pvGfl->attribute("xi_i_dq_pu"));
  logger->logAttribute("GFL_pv_P_filtered_pu",
                       pvGfl->attribute("P_filtered_pu"));
  logger->logAttribute("GFL_pv_Q_filtered_pu",
                       pvGfl->attribute("Q_filtered_pu"));
  logger->logAttribute("GFL_pv_pll_integrator",
                       pvGfl->attribute("pll_integrator"));
  logger->logAttribute("GFL_pv_theta", pvGfl->attribute("theta"));
  logger->logAttribute("GFL_pv_i_lf_abc", pvGfl->attribute("i_lf"));
  logger->logAttribute("GFL_pv_i_cf_abc", pvGfl->attribute("i_cf"));

  Simulation simulation(simulationName, Logger::Level::info);

  if (p.simulation.enableBreakerClosing) {
    simulation.addEvent(
        SwitchEvent3Ph::make(p.simulation.breakerCloseTime, breaker, true));
  }

  if (p.simulation.enableGflReactivePowerStep &&
      p.simulation.gflReactivePowerStepTime <= p.simulation.finalTime) {
    const Real eventTime = std::round(p.simulation.gflReactivePowerStepTime /
                                      p.simulation.timeStep) *
                           p.simulation.timeStep;
    const Real reactivePowerReferencePu =
        p.simulation.gflReactivePowerAfterStep / p.pvGfl.ratedPower;

    simulation.addEvent(AttributeEvent<Real>::make(
        eventTime, pvGfl->mReactivePowerRefPu, reactivePowerReferencePu));
  }

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

} // namespace ScenarioASynGenVBRSynGenVBRGflSiemens

int main(int argc, char *argv[]) {
  try {
    ScenarioASynGenVBRSynGenVBRGflSiemens::Parameters parameters;
    CommandLineArgs args(argc, argv);

    if (argc > 1) {
      parameters.simulation.timeStep = args.timeStep;
      parameters.simulation.finalTime = args.duration;

      if (args.name != "dpsim")
        parameters.simulation.name = args.name;

      // Apply a diagnostic preset first. Explicit options below override it.
      if (args.options.find("DIAGNOSTIC_CASE") != args.options.end()) {
        const Real value = args.getOptionReal("DIAGNOSTIC_CASE");
        if (std::floor(value) != value)
          throw std::invalid_argument("DIAGNOSTIC_CASE must be an integer");
        ScenarioASynGenVBRSynGenVBRGflSiemens::applyDiagnosticCase(
            parameters, static_cast<Int>(value));
      }

      if (args.options.find("BREAKER_CLOSE_TIME") != args.options.end()) {
        parameters.simulation.breakerCloseTime =
            args.getOptionReal("BREAKER_CLOSE_TIME");
      }
      if (args.options.find("PV_Q_STEP_TIME") != args.options.end()) {
        parameters.simulation.gflReactivePowerStepTime =
            args.getOptionReal("PV_Q_STEP_TIME");
      }
      if (args.options.find("PV_Q_AFTER_STEP") != args.options.end()) {
        parameters.simulation.gflReactivePowerAfterStep =
            args.getOptionReal("PV_Q_AFTER_STEP");
      }

      if (args.options.find("ENABLE_Q_STEP") != args.options.end()) {
        parameters.simulation.enableGflReactivePowerStep =
            ScenarioASynGenVBRSynGenVBRGflSiemens::boolFromOption(
                args.getOptionReal("ENABLE_Q_STEP"), "ENABLE_Q_STEP");
      }
      if (args.options.find("ENABLE_BREAKER") != args.options.end()) {
        parameters.simulation.enableBreakerClosing =
            ScenarioASynGenVBRSynGenVBRGflSiemens::boolFromOption(
                args.getOptionReal("ENABLE_BREAKER"), "ENABLE_BREAKER");
      }

      if (args.options.find("ENABLE_GAS_GOVERNOR") != args.options.end()) {
        parameters.diagnostics.enableGasGovernor =
            ScenarioASynGenVBRSynGenVBRGflSiemens::boolFromOption(
                args.getOptionReal("ENABLE_GAS_GOVERNOR"),
                "ENABLE_GAS_GOVERNOR");
      }
      if (args.options.find("ENABLE_GAS_EXCITER") != args.options.end()) {
        parameters.diagnostics.enableGasExciter =
            ScenarioASynGenVBRSynGenVBRGflSiemens::boolFromOption(
                args.getOptionReal("ENABLE_GAS_EXCITER"), "ENABLE_GAS_EXCITER");
      }
      if (args.options.find("ENABLE_PSH_GOVERNOR") != args.options.end()) {
        parameters.diagnostics.enablePshGovernor =
            ScenarioASynGenVBRSynGenVBRGflSiemens::boolFromOption(
                args.getOptionReal("ENABLE_PSH_GOVERNOR"),
                "ENABLE_PSH_GOVERNOR");
      }
      if (args.options.find("ENABLE_PSH_EXCITER") != args.options.end()) {
        parameters.diagnostics.enablePshExciter =
            ScenarioASynGenVBRSynGenVBRGflSiemens::boolFromOption(
                args.getOptionReal("ENABLE_PSH_EXCITER"), "ENABLE_PSH_EXCITER");
      }
      if (args.options.find("ENABLE_GFL_CONTROL") != args.options.end()) {
        parameters.diagnostics.enableGflControl =
            ScenarioASynGenVBRSynGenVBRGflSiemens::boolFromOption(
                args.getOptionReal("ENABLE_GFL_CONTROL"), "ENABLE_GFL_CONTROL");
      }

      if (args.options.find("GFL_KFP") != args.options.end())
        parameters.gfl.frequencyToActivePowerGainPu =
            args.getOptionReal("GFL_KFP");
      if (args.options.find("GFL_KVQ") != args.options.end())
        parameters.gfl.voltageToReactivePowerGainPu =
            args.getOptionReal("GFL_KVQ");
      if (args.options.find("GFL_PLL_GAIN_SCALE") != args.options.end())
        parameters.gfl.pllGainScale = args.getOptionReal("GFL_PLL_GAIN_SCALE");
      if (args.options.find("GFL_CURRENT_GAIN_SCALE") != args.options.end())
        parameters.gfl.currentControllerGainScale =
            args.getOptionReal("GFL_CURRENT_GAIN_SCALE");
      if (args.options.find("GFL_CURRENT_FILTER_TAU") != args.options.end())
        parameters.gfl.currentMeasurementTimeConstant =
            args.getOptionReal("GFL_CURRENT_FILTER_TAU");

      if (args.options.find("GAS_XL") != args.options.end())
        parameters.gasGenerator.leakageInductance =
            args.getOptionReal("GAS_XL");
      if (args.options.find("GAS_XDPP") != args.options.end())
        parameters.gasGenerator.ldSubtransient = args.getOptionReal("GAS_XDPP");
      if (args.options.find("GAS_XQPP") != args.options.end())
        parameters.gasGenerator.lqSubtransient = args.getOptionReal("GAS_XQPP");

      if (args.options.find("LOG_DOWNSAMPLING") != args.options.end()) {
        const Real value = args.getOptionReal("LOG_DOWNSAMPLING");
        if (!(value >= 1.0) || std::floor(value) != value ||
            value > static_cast<Real>(std::numeric_limits<UInt>::max())) {
          throw std::invalid_argument(
              "LOG_DOWNSAMPLING must be an integer >= 1");
        }
        parameters.simulation.logDownsampling = static_cast<UInt>(value);
      }

      if (args.options.find("PSH_INITIAL_ANGLE_DEG") != args.options.end()) {
        parameters.pshInitialAngleDegrees =
            args.getOptionReal("PSH_INITIAL_ANGLE_DEG");
      } else if (args.options.find("PSHA_INITIAL_ANGLE_DEG") !=
                 args.options.end()) {
        parameters.pshInitialAngleDegrees =
            args.getOptionReal("PSHA_INITIAL_ANGLE_DEG");
      }
    }

    ScenarioASynGenVBRSynGenVBRGflSiemens::validateParameters(parameters);
    ScenarioASynGenVBRSynGenVBRGflSiemens::logDiagnosticSummary(parameters);

    const auto powerFlow =
        ScenarioASynGenVBRSynGenVBRGflSiemens::buildAndRunPowerFlow(parameters);

    ScenarioASynGenVBRSynGenVBRGflSiemens::runEmt(parameters, powerFlow);

    return 0;
  } catch (const std::exception &exception) {
    SPDLOG_ERROR(
        "Scenario A SynGenVBR/GFL/SynGenVBR diagnostic example failed:\n{}",
        exception.what());
    return 1;
  }
}
