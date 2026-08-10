// SPDX-FileCopyrightText: 2017-2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
// FIXED: 50 MVA gas-GFM base, 11 kV PV-GFL base, rescaled PU filters, and lossless zero-R PV transformer.

#include <DPsim.h>

#include "../Examples.h"
#include <dpsim-models/EMT/EMT_Ph3_GFL_Siemens.h>
#include <dpsim-models/EMT/EMT_Ph3_GFM_Siemens.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>

using namespace DPsim;
using namespace CPS;
namespace generic_model_B_A {

struct ScenarioConfig {
  // Network
  Real Vnom = 220e3;
  Real nomFreq = 50.0;
  Real nomOmega = nomFreq * 2.0 * PI;

  // Generator 1 / gas-converter bus
  Real nomPower_G1 = 50e6;
  Real nomPhPhVoltRMS_G1 = 10.5e3;
  Real nomFreq_G1 = 50.0;
  Real H_G1 = 5.0;
  Real Xpd_G1 = 0.31;
  Real Rs_G1 = 0.002;
  Real D_G1 = 0.0;

  // Initialization parameters
  Real initActivePower_G1 = 0.30880873930480557e6;
  Real setPointVoltage_G1 = nomPhPhVoltRMS_G1;
  Real initMechPower_G1 = 0.30880873930480557e6;

  // Transformers
  Real t1_ratio = Vnom / nomPhPhVoltRMS_G1;

  // Loads used by the original power-flow scenario
  Real activePower_L_bus_A = 0.1533e6;
  Real reactivePower_L_bus_A = 15.3267e6;
  Real activePower_L_bus_B = 0.1533e6;
  Real reactivePower_L_bus_B = 15.3268e6;

  // Shunt at bus B
  Real shuntConduntanceB = 3.0989e-06;
  Real shuntSusceptanceB = -3.0989e-04;

  // Shunt at bus A
  Real shuntConduntanceA = 3.0989e-06;
  Real shuntSusceptanceA = -3.0989e-04;

  // Cable 1-2, 10.5 km
  Real cableResistance = 0.032 * 10.5;
  Real cableInductance = 0.321493e-3 * 10.5;
  Real cableCapacitance = 0.21915e-6 * 10.5;
  Real cableConductance = 1e-15;

  // Line 3, 5 km
  Real lineResistance3 = 0.0749 * 5.0;
  Real lineInductance3 = 1.270693e-3 * 5.0;
  Real lineCapacitance3 = 0.00466961e-6 * 5.0;
  Real lineConductance3 = 1e-15;
};

} // namespace generic_model_B_A

namespace ScenarioAGfmGflSiemens {

const generic_model_B_A::ScenarioConfig scenarioConfig{};

constexpr Real SWITCH_OPEN_RESISTANCE = 1e12;
constexpr Real SWITCH_CLOSED_RESISTANCE = 1e-8;

// =============================================================================
// Scenario parameters retained from the original case
// =============================================================================

struct ScenarioParameters {
  String name = "EMT_Scenario_A_GFM_GFL_GFM";
  Real frequency = 50.0;
  // The Siemens validation cases use 20 us. The retained per-unit RLC filter
  // has an approximately 5 kHz resonance, so the original 100 us step is
  // at/below the Nyquist limit and causes the EMT converter states to diverge.
  Real timeStep = 100e-6;
  Real finalTime = 15.0;

  Bool startPshBreakerEvent = true;
  Bool endPshBreakerEvent = true;
  Real startPshBreakerTime = 10.0;
  Real endPshBreakerTime = 11100.0;

  Real gflReactivePowerStepTime = 5.0;
  Real gflReactivePowerAfterStep = -25e6;

  UInt logDownsampling = 10;
};

struct ScenarioRatings {
  // Gas GFM replaces generator G1. setBaseParameters() expects an apparent-
  // power base, therefore the scenario's 50 MW nominal rating is represented
  // by S_base=50 MVA (unity-power-factor nominal interpretation).
  Real gasConverterPower = scenarioConfig.nomPower_G1;
  Real gasConverterVoltage = scenarioConfig.nomPhPhVoltRMS_G1;

  // Retain the original PSH generator and transformer rating.
  Real pshConverterPower = 200e6;
  Real pshConverterVoltage = 18.0e3;

  // Corrected PV converter terminal base.
  Real pvConverterPower = 50e6;
  Real pvConverterVoltage = 11.0e3;

  Real networkVoltage = scenarioConfig.Vnom;

  // Transformer and connected-converter power bases are kept consistent.
  Real gasTransformerPower = scenarioConfig.nomPower_G1;
  Real pshTransformerPower = 200e6;
  Real pvTransformerPower = 50e6;
};

// =============================================================================
// Siemens validation parameters
// =============================================================================

/// GFM values copied from EMT_GFM_Siemens_LoadStep_Initial_0_5pu.
struct GfmValidationParameters {
  Real activePowerDroopPu = 0.02;
  Real reactivePowerDroopPu = 0.0311;

  Real activePowerMeasurementTimeConstant = 0.1;
  Real reactivePowerMeasurementTimeConstant = 0.1;

  Real voltageControllerKp = 0.52;
  Real voltageControllerKi = 1.16;
  Real outputCurrentFeedforwardGain = 1.0;

  Real currentControllerKp = 0.74;
  Real currentControllerKi = 1.19;
  Real pccVoltageFeedforwardGain = 1.0;

  // The validation example bypasses the 1 us Simulink PWM delay for the EMT
  // time step. Preserve that validated setting here.
  Real pwmDelayTimeConstant = 0.0;

  Real minimumFrequencyPu = 0.80;
  Real maximumFrequencyPu = 1.20;
  Real maximumCurrentReferencePu = 2.0;
  Real maximumVoltageCommandPu = 1.20;
};

/// GFL inner-loop values copied from EMT_GFL_Siemens_PQ_SetpointStep.
struct GflValidationParameters {
  // The validation inverse-droop gains K_fP=K_VQ=20 were designed for a
  // stiff-grid standalone test. In this converter-dominated scenario they
  // drove P_command/Q_command into current limiting. Use direct P/Q control:
  // P_command=P_ref and Q_command=Q_ref.
  Real frequencyToActivePowerGainPu = 4.0;
  Real voltageToReactivePowerGainPu = 2.0; //4.0 works as well

  Real pllKp = 0.449 / (9.1e-3 + 250e-6);
  Real pllKi = 250e-6 / (5.9 * (9.1e-3 + 250e-6));

  Real currentControllerKp = 0.74 / 10.0;
  // Intended Simulink value and GFL_Siemens class default: 1.19 / 10.
  // The validation example previously contained a transcription typo
  // (11.9 / 10), which is ten times too large.
  Real currentControllerKi = 1.19 / 10.0;
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

/// The physical validation filter belongs to a 50 kVA / 257 V converter.
/// Convert it once to per unit and then apply the same per-unit filter to each
/// transmission-scale converter base in this scenario.
struct SiemensValidationFilter {
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

  Real scaledResistance(Real ratedPower, Real ratedVoltage) const {
    const Real zBase = ratedVoltage * ratedVoltage / ratedPower;
    return resistancePu() * zBase;
  }

  Real scaledInductance(Real ratedPower, Real ratedVoltage,
                        Real frequency) const {
    const Real zBase = ratedVoltage * ratedVoltage / ratedPower;
    return inductiveReactancePu() * zBase / (2.0 * PI * frequency);
  }

  Real scaledCapacitance(Real ratedPower, Real ratedVoltage,
                         Real frequency) const {
    const Real zBase = ratedVoltage * ratedVoltage / ratedPower;
    return capacitiveSusceptancePu() / ((2.0 * PI * frequency) * zBase);
  }
};

struct Parameters {
  ScenarioParameters simulation;
  ScenarioRatings ratings;
  GfmValidationParameters gfm;
  GflValidationParameters gfl;
  SiemensValidationFilter filter;
};

struct PowerFlowResult {
  SystemTopology system;

  Complex gasPower;
  Complex pshPower;
  Complex pvPower;

  Complex gasBusVoltage;
  Complex pshBusVoltage;
  Complex pvBusVoltage;
};

// =============================================================================
// Helpers
// =============================================================================

Real alignedEventTime(Real requestedTime, Real timeStep) {
  return std::round(requestedTime / timeStep) * timeStep;
}

void validateParameters(const Parameters &p) {
  if (!(p.simulation.frequency > 0.0) || !(p.simulation.timeStep > 0.0) ||
      !(p.simulation.finalTime > 0.0) || !(p.simulation.logDownsampling > 0)) {
    throw std::invalid_argument(
        "Require frequency>0, timeStep>0, finalTime>0 and "
        "logDownsampling>0");
  }

  // The rescaled Siemens RLC filter retains the validation model's fast
  // dynamics. A 100 us EMT step is too coarse for this filter/controller.
  // Keep a firm guard here so command-line overrides cannot silently recreate
  // the non-finite-state failure.

  if (p.simulation.startPshBreakerEvent &&
      !(p.simulation.startPshBreakerTime >= 0.0)) {
    throw std::invalid_argument("Invalid PSH breaker closing time");
  }

  if (p.simulation.endPshBreakerEvent &&
      !(p.simulation.endPshBreakerTime > p.simulation.startPshBreakerTime)) {
    throw std::invalid_argument(
        "PSH breaker opening time must be later than its closing time");
  }

  if (!(p.ratings.gasConverterPower > 0.0) ||
      !(p.ratings.pshConverterPower > 0.0) ||
      !(p.ratings.pvConverterPower > 0.0) ||
      !(p.ratings.gasConverterVoltage > 0.0) ||
      !(p.ratings.pshConverterVoltage > 0.0) ||
      !(p.ratings.pvConverterVoltage > 0.0) ||
      !(p.ratings.networkVoltage > 0.0)) {
    throw std::invalid_argument("Invalid converter or voltage base");
  }

  if (!(p.filter.inductiveReactancePu() > 0.0) ||
      !(p.filter.capacitiveSusceptancePu() > 0.0) ||
      !(p.filter.resistancePu() >= 0.0)) {
    throw std::invalid_argument("Invalid Siemens validation filter");
  }

  const auto sameBase = [](Real converterPower, Real transformerPower) {
    const Real scale =
        std::max(std::abs(converterPower), std::abs(transformerPower));
    return std::abs(converterPower - transformerPower) <= 1e-12 * scale;
  };

  if (!sameBase(p.ratings.gasConverterPower, p.ratings.gasTransformerPower) ||
      !sameBase(p.ratings.pshConverterPower, p.ratings.pshTransformerPower) ||
      !sameBase(p.ratings.pvConverterPower, p.ratings.pvTransformerPower)) {
    throw std::invalid_argument(
        "Each converter power base must equal its transformer rating");
  }
}

void logParameterSummary(const Parameters &p) {
  const Real gasRf = p.filter.scaledResistance(p.ratings.gasConverterPower,
                                               p.ratings.gasConverterVoltage);
  const Real gasLf = p.filter.scaledInductance(p.ratings.gasConverterPower,
                                               p.ratings.gasConverterVoltage,
                                               p.simulation.frequency);
  const Real gasCf = p.filter.scaledCapacitance(p.ratings.gasConverterPower,
                                                p.ratings.gasConverterVoltage,
                                                p.simulation.frequency);

  const Real pshRf = p.filter.scaledResistance(p.ratings.pshConverterPower,
                                               p.ratings.pshConverterVoltage);
  const Real pshLf = p.filter.scaledInductance(p.ratings.pshConverterPower,
                                               p.ratings.pshConverterVoltage,
                                               p.simulation.frequency);
  const Real pshCf = p.filter.scaledCapacitance(p.ratings.pshConverterPower,
                                                p.ratings.pshConverterVoltage,
                                                p.simulation.frequency);

  const Real pvRf = p.filter.scaledResistance(p.ratings.pvConverterPower,
                                              p.ratings.pvConverterVoltage);
  const Real pvLf = p.filter.scaledInductance(p.ratings.pvConverterPower,
                                              p.ratings.pvConverterVoltage,
                                              p.simulation.frequency);
  const Real pvCf = p.filter.scaledCapacitance(p.ratings.pvConverterPower,
                                               p.ratings.pvConverterVoltage,
                                               p.simulation.frequency);

  SPDLOG_INFO(
      "Scenario A with Siemens GFM/GFL converters:"
      "\n  EMT: dt={} s, finalTime={} s, logDownsampling={}"
      "\n  gas GFM: S_base={} VA, V_base_LL={} V, transformer={} VA"
      "\n  PSH GFM: S_base={} VA, V_base_LL={} V, transformer={} VA"
      "\n  PV GFL: S_base={} VA, V_base_LL={} V, transformer={} VA, ratio={}"
      "\n  common validation filter: Rf={} pu, X_Lf={} pu, B_Cf={} pu"
      "\n  gas filter: Rf={} Ohm, Lf={} H, Cf={} F"
      "\n  PSH filter: Rf={} Ohm, Lf={} H, Cf={} F"
      "\n  PV 11 kV filter: Rf={} Ohm, Lf={} H, Cf={} F"
      "\n  GFM gains: kp={}, kq={}, Kpv={}, Kiv={}, Kpi={}, Kii={}"
      "\n  GFL outer mode: direct P/Q (KfP={}, KVQ={})"
      "\n  GFL PLL/current gains: KpPLL={}, KiPLL={}, Kpi={}, Kii={}"
      "\n  PV transformer leakage: X={} pu on 220 kV / 50 MVA base"
      "\n  GFL Q step: {} var = {} pu at t={} s",
      p.simulation.timeStep, p.simulation.finalTime,
      p.simulation.logDownsampling, p.ratings.gasConverterPower,
      p.ratings.gasConverterVoltage, p.ratings.gasTransformerPower,
      p.ratings.pshConverterPower, p.ratings.pshConverterVoltage,
      p.ratings.pshTransformerPower, p.ratings.pvConverterPower,
      p.ratings.pvConverterVoltage, p.ratings.pvTransformerPower,
      p.ratings.networkVoltage / p.ratings.pvConverterVoltage,
      p.filter.resistancePu(), p.filter.inductiveReactancePu(),
      p.filter.capacitiveSusceptancePu(), gasRf, gasLf, gasCf, pshRf, pshLf,
      pshCf, pvRf, pvLf, pvCf, p.gfm.activePowerDroopPu,
      p.gfm.reactivePowerDroopPu, p.gfm.voltageControllerKp,
      p.gfm.voltageControllerKi, p.gfm.currentControllerKp,
      p.gfm.currentControllerKi, p.gfl.frequencyToActivePowerGainPu,
      p.gfl.voltageToReactivePowerGainPu, p.gfl.pllKp, p.gfl.pllKi,
      p.gfl.currentControllerKp, p.gfl.currentControllerKi,
      2.0 * PI * p.simulation.frequency * 0.928e-3 /
          (p.ratings.networkVoltage * p.ratings.networkVoltage /
           p.ratings.pvTransformerPower),
      p.simulation.gflReactivePowerAfterStep,
      p.simulation.gflReactivePowerAfterStep / p.ratings.pvConverterPower,
      p.simulation.gflReactivePowerStepTime);
}

void configureGfm(const std::shared_ptr<EMT::Ph3::GFM_Siemens> &gfm,
                  Real ratedPower, Real ratedVoltage, Real frequency,
                  const Complex &initialPower, const Complex &initialVoltage,
                  const Parameters &p) {
  const Real initialVoltagePu = std::abs(initialVoltage) / ratedVoltage;
  const Real initialActivePowerPu = initialPower.real() / ratedPower;
  const Real initialReactivePowerPu = initialPower.imag() / ratedPower;

  gfm->setBaseParameters(ratedPower, ratedVoltage, frequency);
  gfm->setReferencesPerUnit(1.0, initialVoltagePu, initialActivePowerPu,
                            initialReactivePowerPu);

  gfm->setDroopParametersPerUnit(p.gfm.activePowerDroopPu,
                                 p.gfm.reactivePowerDroopPu);
  gfm->setPowerMeasurementFilterTimeConstants(
      p.gfm.activePowerMeasurementTimeConstant,
      p.gfm.reactivePowerMeasurementTimeConstant);
  gfm->setVoltageControllerParameters(p.gfm.voltageControllerKp,
                                      p.gfm.voltageControllerKi,
                                      p.gfm.outputCurrentFeedforwardGain);
  gfm->setCurrentControllerParameters(p.gfm.currentControllerKp,
                                      p.gfm.currentControllerKi,
                                      p.gfm.pccVoltageFeedforwardGain);
  gfm->setPwmDelayTimeConstant(p.gfm.pwmDelayTimeConstant);
  gfm->setControllerLimitsPerUnit(
      p.gfm.minimumFrequencyPu, p.gfm.maximumFrequencyPu,
      p.gfm.maximumCurrentReferencePu, p.gfm.maximumVoltageCommandPu);

  gfm->setFilterParametersPerUnit(p.filter.inductiveReactancePu(),
                                  p.filter.capacitiveSusceptancePu(),
                                  p.filter.resistancePu());
  gfm->withControl(true);
}

void configureGfl(const std::shared_ptr<EMT::Ph3::GFL_Siemens> &gfl,
                  Real ratedPower, Real ratedVoltage, Real frequency,
                  const Complex &initialPower, const Complex &initialVoltage,
                  const Parameters &p) {
  const Real initialVoltagePu = std::abs(initialVoltage) / ratedVoltage;
  const Real initialActivePowerPu = initialPower.real() / ratedPower;
  const Real initialReactivePowerPu = initialPower.imag() / ratedPower;

  gfl->setBaseParameters(ratedPower, ratedVoltage, frequency);
  gfl->setReferencesPerUnit(1.0, initialVoltagePu, initialActivePowerPu,
                            initialReactivePowerPu);

  // PLL and current-controller signals are normalized internally. Their
  // validation gains therefore remain unchanged when V_base changes from
  // 1.5 kV to 11 kV. setFilterParametersPerUnit() below performs the required
  // physical R-L-C rescaling to the new voltage/current base.
  gfl->setDroopParametersPerUnit(p.gfl.frequencyToActivePowerGainPu,
                                 p.gfl.voltageToReactivePowerGainPu);
  gfl->setPllParameters(p.gfl.pllKp, p.gfl.pllKi);
  gfl->setCurrentControllerParameters(p.gfl.currentControllerKp,
                                      p.gfl.currentControllerKi,
                                      p.gfl.pccVoltageFeedforwardGain);
  gfl->setMeasurementFilterTimeConstants(
      p.gfl.activePowerMeasurementTimeConstant,
      p.gfl.reactivePowerMeasurementTimeConstant,
      p.gfl.currentMeasurementTimeConstant);
  gfl->setControllerLimitsPerUnit(
      p.gfl.minimumFrequencyPu, p.gfl.maximumFrequencyPu,
      p.gfl.maximumCurrentReferencePu, p.gfl.maximumVoltageCommandPu,
      p.gfl.minimumVoltageForCurrentReferencePu);

  gfl->setFilterParametersPerUnit(p.filter.inductiveReactancePu(),
                                  p.filter.capacitiveSusceptancePu(),
                                  p.filter.resistancePu());
  gfl->withControl(true);
}

// =============================================================================
// SP power-flow initialization
// =============================================================================

PowerFlowResult buildAndRunPowerFlow(const Parameters &p) {
  const String simulationName = p.simulation.name + "_PF";
  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto busGas = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busA = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busB = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busPsha = SimNode<Complex>::make("BUS_psha", PhaseType::Single);
  auto busPsh = SimNode<Complex>::make("BUS_psh", PhaseType::Single);
  auto busTrPsh = SimNode<Complex>::make("BUS_tr_psh", PhaseType::Single);
  auto busPv = SimNode<Complex>::make("BUS_pv", PhaseType::Single);

  // These SP synchronous-generator objects are power-flow surrogates only.
  // The EMT system contains no synchronous-machine dynamic model.
  auto gasPfEquivalent =
      SP::Ph1::SynchronGenerator::make("GEN_gas", Logger::Level::debug);
  gasPfEquivalent->setParameters(
      p.ratings.gasConverterPower, p.ratings.gasConverterVoltage,
      scenarioConfig.initActivePower_G1, p.ratings.gasConverterVoltage,
      PowerflowBusType::VD);
  gasPfEquivalent->setBaseVoltage(p.ratings.gasConverterVoltage);

  auto pshPfEquivalent =
      SP::Ph1::SynchronGenerator::make("GEN_psh", Logger::Level::debug);
  pshPfEquivalent->setParameters(
      p.ratings.pshConverterPower, p.ratings.pshConverterVoltage, 0.0,
      p.ratings.pshConverterVoltage, PowerflowBusType::PV);
  pshPfEquivalent->setBaseVoltage(p.ratings.pshConverterVoltage);

  const Real gasTransformerResistance = 0.00376196 * p.ratings.networkVoltage *
                                        p.ratings.networkVoltage /
                                        p.ratings.gasTransformerPower / 2.0;
  const Real gasTransformerInductance =
      0.1007298 * p.ratings.networkVoltage * p.ratings.networkVoltage /
      p.ratings.gasTransformerPower / (2.0 * PI * p.simulation.frequency) / 2.0;

  auto transformerGas =
      std::make_shared<SP::Ph1::Transformer>("TR_gas", Logger::Level::debug);
  transformerGas->setParameters(
      p.ratings.gasConverterVoltage, p.ratings.networkVoltage,
      p.ratings.gasTransformerPower,
      p.ratings.gasConverterVoltage / p.ratings.networkVoltage, 0.0,
      2.0 * gasTransformerResistance, 2.0 * gasTransformerInductance);
  transformerGas->setBaseVoltage(p.ratings.networkVoltage);

  auto transformerPsh =
      std::make_shared<SP::Ph1::Transformer>("TR_psh", Logger::Level::debug);
  transformerPsh->setParameters(
      p.ratings.pshConverterVoltage, p.ratings.networkVoltage,
      p.ratings.pshTransformerPower,
      p.ratings.pshConverterVoltage / p.ratings.networkVoltage, 0.0,
      2.0 * 0.41745, 2.0 * 0.0481260775594524);
  transformerPsh->setBaseVoltage(p.ratings.networkVoltage);

  // The old AvVoltageSourceInverterDQ contained this transformer internally.
  // GFL_Siemens does not, so retain it explicitly with the same ratings and
  // impedance from the supplied scenario.
  auto transformerPv =
      std::make_shared<SP::Ph1::Transformer>("TR_pv", Logger::Level::debug);

  // Preserve the orientation of the transformer embedded in the original
  // AvVoltageSourceInverterDQ:
  //   end 1 = 220 kV system side, end 2 = 11 kV converter side.
  // The leakage inductance is therefore interpreted on the same side as in
  // the original component. Reversing the transformer would refer 0.928 mH
  // to the 11 kV base and create an artificial reactance of about 6.48 pu.
  transformerPv->setParameters(
      p.ratings.networkVoltage, p.ratings.pvConverterVoltage,
      p.ratings.pvTransformerPower,
      p.ratings.networkVoltage / p.ratings.pvConverterVoltage, 0.0, 0.0,
      0.928e-3);
  transformerPv->setBaseVoltage(p.ratings.networkVoltage);

  const Real cableCapacitance = 1.93522865e-6;
  auto cable = SP::Ph1::PiLine::make("cable", Logger::Level::debug);
  cable->setParameters(5.04669, 0.13523123641, cableCapacitance, 1e-15);
  cable->setBaseVoltage(p.ratings.networkVoltage);

  auto line3 = SP::Ph1::PiLine::make("line_3", Logger::Level::debug);
  line3->setParameters(
      scenarioConfig.lineResistance3, scenarioConfig.lineInductance3,
      scenarioConfig.lineCapacitance3, scenarioConfig.lineConductance3);
  line3->setBaseVoltage(p.ratings.networkVoltage);

  auto pshBreakerEquivalent =
      SP::Ph1::PiLine::make("breaker_GEN_psh", Logger::Level::debug);
  pshBreakerEquivalent->setParameters(1e8, 1e-5, 0.0, 1e-15);
  pshBreakerEquivalent->setBaseVoltage(p.ratings.networkVoltage);

  // Zero-power PF surrogate at the GFL low-voltage terminal. It establishes
  // the matching BUS_pv voltage and terminal P/Q without changing the original
  // operating point.
  auto pvPfEquivalent = SP::Ph1::Load::make("pv", Logger::Level::debug);
  pvPfEquivalent->setParameters(0.0, 0.0, p.ratings.pvConverterVoltage);

  gasPfEquivalent->connect({busGas});
  transformerGas->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPsha});
  pshBreakerEquivalent->connect({busPsha, busTrPsh});
  transformerPsh->connect({busPsh, busTrPsh});
  pshPfEquivalent->connect({busPsh});

  transformerPv->connect({busA, busPv});
  pvPfEquivalent->connect({busPv});

  SystemTopology system(
      p.simulation.frequency,
      SystemNodeList{busGas, busB, busA, busPsha, busTrPsh, busPsh, busPv},
      SystemComponentList{gasPfEquivalent, transformerGas, cable, line3,
                          pshBreakerEquivalent, transformerPsh, pshPfEquivalent,
                          transformerPv, pvPfEquivalent});

  auto logger = DataLogger::make(simulationName);
  logger->logAttribute("BUS_gas.v", busGas->attribute("v"));
  logger->logAttribute("BUS_b.v", busB->attribute("v"));
  logger->logAttribute("BUS_a.v", busA->attribute("v"));
  logger->logAttribute("BUS_psha.v", busPsha->attribute("v"));
  logger->logAttribute("BUS_psh.v", busPsh->attribute("v"));
  logger->logAttribute("BUS_pv.v", busPv->attribute("v"));

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(1.0);
  simulation.setFinalTime(1.0);
  simulation.setDomain(Domain::SP);
  simulation.setSolverType(Solver::Type::NRP);
  simulation.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.addLogger(logger);
  simulation.run();

  // terminal(0)->singlePower() remains zero for these SP generator models in
  // the NRP initialization path. getApparentPower() contains the solved
  // generator injection, including line/cable charging and network losses.
  const Complex gasPower = gasPfEquivalent->getApparentPower();
  const Complex pshPower = pshPfEquivalent->getApparentPower();
  const Complex pvPower(0.0, 0.0);

  const Complex gasVoltage = busGas->singleVoltage();
  const Complex pshVoltage = busPsh->singleVoltage();
  const Complex pvVoltage = busPv->singleVoltage();

  const auto finiteComplex = [](const Complex &value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  };

  if (!finiteComplex(gasPower) || !finiteComplex(pshPower) ||
      !finiteComplex(pvPower) || !finiteComplex(gasVoltage) ||
      !finiteComplex(pshVoltage) || !finiteComplex(pvVoltage) ||
      !(std::abs(gasVoltage) > 0.0) || !(std::abs(pshVoltage) > 0.0) ||
      !(std::abs(pvVoltage) > 0.0)) {
    throw std::runtime_error("Power-flow initialization produced a non-finite "
                             "converter operating point");
  }

  SPDLOG_INFO(
      "Power-flow operating point:"
      "\n  gas GFM equivalent: P={} W ({} pu), Q={} var ({} pu), "
      "|V|={} V_LL RMS ({} pu)"
      "\n  PSH GFM equivalent: P={} W ({} pu), Q={} var ({} pu), "
      "|V|={} V_LL RMS ({} pu)"
      "\n  PV GFL equivalent: P={} W ({} pu), Q={} var ({} pu), "
      "|V|={} V_LL RMS ({} pu)",
      gasPower.real(), gasPower.real() / p.ratings.gasConverterPower,
      gasPower.imag(), gasPower.imag() / p.ratings.gasConverterPower,
      std::abs(gasVoltage),
      std::abs(gasVoltage) / p.ratings.gasConverterVoltage, pshPower.real(),
      pshPower.real() / p.ratings.pshConverterPower, pshPower.imag(),
      pshPower.imag() / p.ratings.pshConverterPower, std::abs(pshVoltage),
      std::abs(pshVoltage) / p.ratings.pshConverterVoltage, pvPower.real(),
      pvPower.real() / p.ratings.pvConverterPower, pvPower.imag(),
      pvPower.imag() / p.ratings.pvConverterPower, std::abs(pvVoltage),
      std::abs(pvVoltage) / p.ratings.pvConverterVoltage);

  return {system,     gasPower,   pshPower, pvPower,
          gasVoltage, pshVoltage, pvVoltage};
}

// =============================================================================
// EMT simulation
// =============================================================================

void runEmt(const Parameters &p, const PowerFlowResult &powerFlow) {
  const String simulationName = p.simulation.name + "_EMT";
  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busA = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busPsha = SimNode<Real>::make("BUS_psha", PhaseType::ABC);
  auto busPsh = SimNode<Real>::make("BUS_psh", PhaseType::ABC);
  auto busTrPsh = SimNode<Real>::make("BUS_tr_psh", PhaseType::ABC);
  auto busPv = SimNode<Real>::make("BUS_pv", PhaseType::ABC);

  auto gfmGas =
      EMT::Ph3::GFM_Siemens::make("GEN_gas", "GEN_gas", Logger::Level::debug);
  configureGfm(gfmGas, p.ratings.gasConverterPower,
               p.ratings.gasConverterVoltage, p.simulation.frequency,
               powerFlow.gasPower, powerFlow.gasBusVoltage, p);

  auto gfmPsh =
      EMT::Ph3::GFM_Siemens::make("GEN_psh", "GEN_psh", Logger::Level::debug);
  configureGfm(gfmPsh, p.ratings.pshConverterPower,
               p.ratings.pshConverterVoltage, p.simulation.frequency,
               powerFlow.pshPower, powerFlow.pshBusVoltage, p);

  auto gflPv = EMT::Ph3::GFL_Siemens::make("pv", "pv", Logger::Level::debug);
  configureGfl(gflPv, p.ratings.pvConverterPower, p.ratings.pvConverterVoltage,
               p.simulation.frequency, powerFlow.pvPower,
               powerFlow.pvBusVoltage, p);

  const Real gasTransformerResistance = 0.00376196 * p.ratings.networkVoltage *
                                        p.ratings.networkVoltage /
                                        p.ratings.gasTransformerPower / 2.0;
  const Real gasTransformerInductance =
      0.1007298 * p.ratings.networkVoltage * p.ratings.networkVoltage /
      p.ratings.gasTransformerPower / (2.0 * PI * p.simulation.frequency) / 2.0;

  auto transformerGas = EMT::Ph3::Transformer::make("TR_gas", "TR_gas",
                                                    Logger::Level::debug, true);
  transformerGas->setParameters(
      p.ratings.gasConverterVoltage, p.ratings.networkVoltage,
      p.ratings.gasTransformerPower,
      p.ratings.gasConverterVoltage / p.ratings.networkVoltage, 0.0,
      Math::singlePhaseParameterToThreePhase(2.0 * gasTransformerResistance),
      Math::singlePhaseParameterToThreePhase(2.0 * gasTransformerInductance));

  auto transformerPsh = EMT::Ph3::Transformer::make("TR_psh", "TR_psh",
                                                    Logger::Level::debug, true);
  transformerPsh->setParameters(
      p.ratings.pshConverterVoltage, p.ratings.networkVoltage,
      p.ratings.pshTransformerPower,
      p.ratings.pshConverterVoltage / p.ratings.networkVoltage, 0.0,
      Math::singlePhaseParameterToThreePhase(2.0 * 0.41745),
      Math::singlePhaseParameterToThreePhase(2.0 * 0.0481260775594524));

  // The original AvVoltageSourceInverterDQ creates its embedded connection
  // transformer with withResistiveLosses=false. Keep the same construction
  // here because the retained PV transformer has R=0. Creating the optional
  // series resistor with an exact zero resistance makes the EMT MNA branch
  // numerically singular/ill-conditioned during the first dynamic solve.
  auto transformerPv = EMT::Ph3::Transformer::make("TR_pv", "TR_pv",
                                                   Logger::Level::debug, false);

  // Same orientation as the transformer embedded in the original PV model:
  // terminal 0 is the 220 kV system side and terminal 1 is the 11 kV side.
  transformerPv->setParameters(
      p.ratings.networkVoltage, p.ratings.pvConverterVoltage,
      p.ratings.pvTransformerPower,
      p.ratings.networkVoltage / p.ratings.pvConverterVoltage, 0.0,
      Math::singlePhaseParameterToThreePhase(0.0),
      Math::singlePhaseParameterToThreePhase(0.928e-3));

  const Real cableCapacitance = 1.93522865e-6;
  auto cable = EMT::Ph3::PiLine::make("cable", Logger::Level::debug);
  cable->setParameters(Math::singlePhaseParameterToThreePhase(5.04669),
                       Math::singlePhaseParameterToThreePhase(0.13523123641),
                       Math::singlePhaseParameterToThreePhase(cableCapacitance),
                       Math::singlePhaseParameterToThreePhase(1e-15));

  auto line3 = EMT::Ph3::PiLine::make("line_3", Logger::Level::debug);
  line3->setParameters(
      Math::singlePhaseParameterToThreePhase(scenarioConfig.lineResistance3),
      Math::singlePhaseParameterToThreePhase(scenarioConfig.lineInductance3),
      Math::singlePhaseParameterToThreePhase(scenarioConfig.lineCapacitance3),
      Math::singlePhaseParameterToThreePhase(scenarioConfig.lineConductance3));

  auto pshBreaker =
      EMT::Ph3::Switch::make("breaker_GEN_psh", Logger::Level::debug);
  pshBreaker->setParameters(
      Math::singlePhaseParameterToThreePhase(SWITCH_OPEN_RESISTANCE),
      Math::singlePhaseParameterToThreePhase(SWITCH_CLOSED_RESISTANCE));
  pshBreaker->openSwitch();

  gfmGas->connect({busGas});
  transformerGas->connect({busGas, busB});
  cable->connect({busB, busA});
  line3->connect({busA, busPsha});
  pshBreaker->connect({busPsha, busTrPsh});
  transformerPsh->connect({busPsh, busTrPsh});
  gfmPsh->connect({busPsh});

  gflPv->connect({busPv});
  transformerPv->connect({busA, busPv});

  SystemTopology system(
      p.simulation.frequency,
      SystemNodeList{busGas, busB, busA, busPsha, busPsh, busTrPsh, busPv},
      SystemComponentList{gfmGas, transformerGas, cable, line3, pshBreaker,
                          transformerPsh, gfmPsh, gflPv, transformerPv});

  system.initWithPowerflow(powerFlow.system, Domain::EMT);

  // Explicitly transfer the generator-positive terminal powers used by the
  // Siemens controllers for their bumpless state initialization.
  gfmGas->terminal(0)->setPower(powerFlow.gasPower);
  gfmPsh->terminal(0)->setPower(powerFlow.pshPower);
  gflPv->terminal(0)->setPower(powerFlow.pvPower);

  auto logger =
      DataLogger::make(simulationName, true, p.simulation.logDownsampling);

  // Network quantities.
  logger->logAttribute("BUS_gas.v", busGas->attribute("v"));
  logger->logAttribute("BUS_b.v", busB->attribute("v"));
  logger->logAttribute("BUS_a.v", busA->attribute("v"));
  logger->logAttribute("BUS_psha.v", busPsha->attribute("v"));
  logger->logAttribute("BUS_psh.v", busPsh->attribute("v"));
  logger->logAttribute("BUS_pv.v", busPv->attribute("v"));

  logger->logAttribute("TR_gas.i", transformerGas->attribute("i_intf"));
  logger->logAttribute("TR_psh.i", transformerPsh->attribute("i_intf"));
  logger->logAttribute("TR_pv.i", transformerPv->attribute("i_intf"));
  logger->logAttribute("cable.i", cable->attribute("i_intf"));
  logger->logAttribute("line_3.i", line3->attribute("i_intf"));
  logger->logAttribute("breaker_GEN_psh.i", pshBreaker->attribute("i_intf"));

  // Gas GFM.
  logger->logAttribute("GFM_gas.P_elec_pu", gfmGas->attribute("P_elec_pu"));
  logger->logAttribute("GFM_gas.Q_elec_pu", gfmGas->attribute("Q_elec_pu"));
  logger->logAttribute("GFM_gas.P_filtered_pu",
                       gfmGas->attribute("P_filtered_pu"));
  logger->logAttribute("GFM_gas.Q_filtered_pu",
                       gfmGas->attribute("Q_filtered_pu"));
  logger->logAttribute("GFM_gas.frequency_pu",
                       gfmGas->attribute("frequency_pu"));
  logger->logAttribute("GFM_gas.V_magnitude_pu",
                       gfmGas->attribute("V_magnitude_pu"));
  logger->logAttribute("GFM_gas.i_pcc_dq_pu", gfmGas->attribute("i_pcc_dq_pu"));
  logger->logAttribute("GFM_gas.i_ref_dq_pu", gfmGas->attribute("i_ref_dq_pu"));
  logger->logAttribute("GFM_gas.v_cmd_dq_pu", gfmGas->attribute("v_cmd_dq_pu"));

  // PSH GFM.
  logger->logAttribute("GFM_psh.P_elec_pu", gfmPsh->attribute("P_elec_pu"));
  logger->logAttribute("GFM_psh.Q_elec_pu", gfmPsh->attribute("Q_elec_pu"));
  logger->logAttribute("GFM_psh.frequency_pu",
                       gfmPsh->attribute("frequency_pu"));
  logger->logAttribute("GFM_psh.V_magnitude_pu",
                       gfmPsh->attribute("V_magnitude_pu"));
  logger->logAttribute("GFM_psh.i_pcc_dq_pu", gfmPsh->attribute("i_pcc_dq_pu"));
  logger->logAttribute("GFM_psh.v_cmd_dq_pu", gfmPsh->attribute("v_cmd_dq_pu"));

  // PV GFL.
  logger->logAttribute("GFL_pv.P_ref_pu", gflPv->attribute("P_ref_pu"));
  logger->logAttribute("GFL_pv.Q_ref_pu", gflPv->attribute("Q_ref_pu"));
  logger->logAttribute("GFL_pv.P_command_pu", gflPv->attribute("P_command_pu"));
  logger->logAttribute("GFL_pv.Q_command_pu", gflPv->attribute("Q_command_pu"));
  logger->logAttribute("GFL_pv.P_elec_pu", gflPv->attribute("P_elec_pu"));
  logger->logAttribute("GFL_pv.Q_elec_pu", gflPv->attribute("Q_elec_pu"));
  logger->logAttribute("GFL_pv.frequency_pu", gflPv->attribute("frequency_pu"));
  logger->logAttribute("GFL_pv.V_magnitude_pu",
                       gflPv->attribute("V_magnitude_pu"));
  logger->logAttribute("GFL_pv.pll_vq_error_pu",
                       gflPv->attribute("pll_vq_error_pu"));
  logger->logAttribute("GFL_pv.i_pcc_dq_pu", gflPv->attribute("i_pcc_dq_pu"));
  logger->logAttribute("GFL_pv.i_ref_dq_pu", gflPv->attribute("i_ref_dq_pu"));
  logger->logAttribute("GFL_pv.v_cmd_dq_pu", gflPv->attribute("v_cmd_dq_pu"));

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(p.simulation.timeStep);
  simulation.setFinalTime(p.simulation.finalTime);
  simulation.setDomain(Domain::EMT);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);

  // Required because the PSH breaker changes the MNA topology.
  simulation.doSystemMatrixRecomputation(true);

  if (p.simulation.startPshBreakerEvent &&
      p.simulation.startPshBreakerTime <= p.simulation.finalTime) {
    const Real eventTime = alignedEventTime(p.simulation.startPshBreakerTime,
                                            p.simulation.timeStep);
    simulation.addEvent(SwitchEvent3Ph::make(eventTime, pshBreaker, true));
  }

  if (p.simulation.endPshBreakerEvent &&
      p.simulation.endPshBreakerTime <= p.simulation.finalTime) {
    const Real eventTime =
        alignedEventTime(p.simulation.endPshBreakerTime, p.simulation.timeStep);
    simulation.addEvent(SwitchEvent3Ph::make(eventTime, pshBreaker, false));
  }

  if (p.simulation.gflReactivePowerStepTime <= p.simulation.finalTime) {
    const Real eventTime = alignedEventTime(
        p.simulation.gflReactivePowerStepTime, p.simulation.timeStep);
    const Real reactivePowerReferencePu =
        p.simulation.gflReactivePowerAfterStep / p.ratings.pvConverterPower;

    simulation.addEvent(AttributeEvent<Real>::make(
        eventTime, gflPv->mReactivePowerRefPu, reactivePowerReferencePu));
  }

  simulation.addLogger(logger);

  SPDLOG_INFO(
      "Starting EMT simulation. GFL Q_ref steps to {} pu at t={} s; "
      "PSH breaker closes at t={} s.",
      p.simulation.gflReactivePowerAfterStep / p.ratings.pvConverterPower,
      p.simulation.gflReactivePowerStepTime, p.simulation.startPshBreakerTime);

  // ===========================================================================
  // Render topology
  //
  // To change the layout, change only options.layout below. All other visual
  // and routing parameters are internal renderer defaults.
  // ===========================================================================

  SystemTopologyRenderer::Options options;
  options.layout = SystemTopologyRenderer::Layout::LeftToRight;

  const std::filesystem::path symbolDirectory =
      "../dpsim-models/resources/Visuals";

  const std::filesystem::path outputFile = "logs/GFM_GFL_GFM/topology.svg";

  SystemTopologyRenderer renderer(system, symbolDirectory, options);

  renderer.renderSvg(outputFile);

  simulation.run();

  SPDLOG_INFO("Simulation completed. Results: logs/{}/{}.csv", simulationName,
              simulationName);
}

} // namespace ScenarioAGfmGflSiemens

int main(int argc, char *argv[]) {
  try {
    ScenarioAGfmGflSiemens::Parameters parameters;
    CommandLineArgs args(argc, argv);

    if (argc > 1) {
      parameters.simulation.timeStep = args.timeStep;
      parameters.simulation.finalTime = args.duration;

      if (args.name != "dpsim")
        parameters.simulation.name = args.name;

      if (args.options.find("STARTTIMEFAULT") != args.options.end()) {
        parameters.simulation.startPshBreakerTime =
            args.getOptionReal("STARTTIMEFAULT");
      }

      if (args.options.find("ENDTIMEFAULT") != args.options.end()) {
        parameters.simulation.endPshBreakerTime =
            args.getOptionReal("ENDTIMEFAULT");
      }

      if (args.options.find("PV_Q_STEP_TIME") != args.options.end()) {
        parameters.simulation.gflReactivePowerStepTime =
            args.getOptionReal("PV_Q_STEP_TIME");
      }

      if (args.options.find("PV_Q_AFTER_STEP") != args.options.end()) {
        parameters.simulation.gflReactivePowerAfterStep =
            args.getOptionReal("PV_Q_AFTER_STEP");
      }

      if (args.options.find("LOG_DOWNSAMPLING") != args.options.end()) {
        const Real value = args.getOptionReal("LOG_DOWNSAMPLING");
        if (!(value >= 1.0) || std::floor(value) != value ||
            value > static_cast<Real>(std::numeric_limits<UInt>::max())) {
          throw std::invalid_argument(
              "LOG_DOWNSAMPLING must be an integer >= 1");
        }
        parameters.simulation.logDownsampling = static_cast<UInt>(value);
      }
    }

    ScenarioAGfmGflSiemens::validateParameters(parameters);
    ScenarioAGfmGflSiemens::logParameterSummary(parameters);

    const auto powerFlow =
        ScenarioAGfmGflSiemens::buildAndRunPowerFlow(parameters);
    ScenarioAGfmGflSiemens::runEmt(parameters, powerFlow);

    return EXIT_SUCCESS;
  } catch (const std::exception &exception) {
    SPDLOG_ERROR("Scenario A GFM/GFL Siemens example failed:\n{}",
                 exception.what());
    return EXIT_FAILURE;
  }
}
