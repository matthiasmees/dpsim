// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <cmath>
#include <stdexcept>

#include <DPsim.h>

using namespace DPsim;
using namespace CPS;

namespace {

struct SeriesRLParameters {
  Real resistance = 0.0; // Ohm per phase
  Real reactance = 0.0;  // Ohm per phase
  Real inductance = 0.0; // H per phase
};

SeriesRLParameters deriveSeriesRL(Real activePower, Real reactivePower,
                                  Real lineToLineVoltageRms, Real frequency) {
  const Real powerSquared =
      activePower * activePower + reactivePower * reactivePower;

  if (!(activePower > 0.0) || !(reactivePower > 0.0) ||
      !(lineToLineVoltageRms > 0.0) || !(frequency > 0.0) ||
      !(powerSquared > 0.0) || !std::isfinite(activePower) ||
      !std::isfinite(reactivePower) || !std::isfinite(lineToLineVoltageRms) ||
      !std::isfinite(frequency)) {
    throw std::invalid_argument(
        "deriveSeriesRL requires finite P>0, Q>0, V_LL>0, and f>0");
  }

  SeriesRLParameters result;
  result.resistance =
      lineToLineVoltageRms * lineToLineVoltageRms * activePower / powerSquared;
  result.reactance = lineToLineVoltageRms * lineToLineVoltageRms *
                     reactivePower / powerSquared;
  result.inductance = result.reactance / (2.0 * PI * frequency);
  return result;
}

/// Islanded GFM droop test with:
///   - an Rf-Lf output filter and series Rd-Cf damping branch,
///   - an energized base series-RL load,
///   - a switched series-RL load,
///   - an SP steady-state initialization whose physical topology mirrors EMT.
struct GfmDroopLoadStepParameters {
  // -----------------------------------------------------------------------
  // Simulation
  // -----------------------------------------------------------------------
  Real timeStep = 50e-6;
  Real finalTime = 0.6;

  Real loadStepStartTime = 0.2;
  Real loadStepEndTime = 0.4;

  // -----------------------------------------------------------------------
  // Electrical base quantities
  // -----------------------------------------------------------------------
  Real nominalFrequency = 50.0; // Hz

  Real phaseVoltageRms = 220.0; // V phase-to-neutral RMS
  Real lineToLineVoltageRms = std::sqrt(3.0) * phaseVoltageRms;
  Real phaseVoltagePeak = std::sqrt(2.0) * phaseVoltageRms;

  Real ratedApparentPower = 20e3; // VA

  // -----------------------------------------------------------------------
  // GFM output filter:
  //
  //   controlled source -- Rf -- Lf -- PCC
  //                                      |
  //                                      Rd
  //                                      |
  //                                      Cf
  //                                      |
  //                                     GND
  // -----------------------------------------------------------------------
  Real filterInductance = 3.0e-3;        // H per phase
  Real filterCapacitance = 20.0e-6;      // F per phase
  Real filterResistance = 0.05;          // Ohm per phase
  Real capacitorDampingResistance = 4.0; // Ohm per phase

  // -----------------------------------------------------------------------
  // Droop controller
  // -----------------------------------------------------------------------
  Real activePowerDroopFraction = 0.05;
  Real reactivePowerDroopFraction = 0.05;

  Real activePowerDroop =
      activePowerDroopFraction * nominalFrequency / ratedApparentPower;
  Real reactivePowerDroop =
      reactivePowerDroopFraction * phaseVoltagePeak / ratedApparentPower;

  Real voltageIntegralGain = 40.0; // 1/s

  Real powerMeasurementFilterCutoffFrequency = 10.0; // Hz
  Real powerMeasurementFilterTimeConstant =
      1.0 / (2.0 * PI * powerMeasurementFilterCutoffFrequency);

  Real minimumFrequency = 0.90 * nominalFrequency;
  Real maximumFrequency = 1.10 * nominalFrequency;
  Real minimumVoltagePeak = 0.70 * phaseVoltagePeak;
  Real maximumVoltagePeak = 1.30 * phaseVoltagePeak;

  // -----------------------------------------------------------------------
  // Loads
  //
  // Both loads are represented as physically damped series R-L branches.
  // Their impedances are derived from the requested total balanced
  // three-phase P+jQ at the nominal line-to-line RMS voltage:
  //
  //   R = V_LL^2 P / (P^2 + Q^2)
  //   X = V_LL^2 Q / (P^2 + Q^2)
  //   L = X / (2*pi*f)
  // -----------------------------------------------------------------------
  Real baseLoadActivePower = 12e3;  // W, total three-phase
  Real baseLoadReactivePower = 2e3; // var, total three-phase

  Real stepLoadActivePower = 3e3;   // W, total three-phase
  Real stepLoadReactivePower = 1e3; // var, total three-phase

  Real switchOpenResistance = 1e9;    // Ohm per phase
  Real switchClosedResistance = 1e-6; // Ohm per phase
};

} // namespace

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  const GfmDroopLoadStepParameters p;
  const String simName = "EMT_GFM_DROOP_PF_MIRRORED_TWO_SERIES_RL_LOADS";

  const SeriesRLParameters baseLoadRL =
      deriveSeriesRL(p.baseLoadActivePower, p.baseLoadReactivePower,
                     p.lineToLineVoltageRms, p.nominalFrequency);

  const SeriesRLParameters stepLoadRL =
      deriveSeriesRL(p.stepLoadActivePower, p.stepLoadReactivePower,
                     p.lineToLineVoltageRms, p.nominalFrequency);

  SPDLOG_INFO("Base series R-L load: R={} Ohm/phase, X={} Ohm/phase, "
              "L={} H/phase, tau={} s",
              baseLoadRL.resistance, baseLoadRL.reactance,
              baseLoadRL.inductance,
              baseLoadRL.inductance / baseLoadRL.resistance);

  SPDLOG_INFO("Step series R-L load: R={} Ohm/phase, X={} Ohm/phase, "
              "L={} H/phase, tau={} s",
              stepLoadRL.resistance, stepLoadRL.reactance,
              stepLoadRL.inductance,
              stepLoadRL.inductance / stepLoadRL.resistance);

  const Real filterResonanceOmega =
      1.0 / std::sqrt(p.filterInductance * p.filterCapacitance);
  const Real filterResonanceFrequency = filterResonanceOmega / (2.0 * PI);
  const Real dampingResistanceRule =
      1.0 / (3.0 * filterResonanceOmega * p.filterCapacitance);

  SPDLOG_INFO("GFM filter design: f_res={} Hz, selected Rd={} Ohm/phase, "
              "rule-of-thumb Rd={} Ohm/phase",
              filterResonanceFrequency, p.capacitorDampingResistance,
              dampingResistanceRule);

  // =========================================================================
  // SP steady-state initialization
  //
  // This is the same energized physical topology as the EMT case:
  //
  //   NetworkInjection/GFM -- Rbase -- base internal node -- Lbase -- GND
  //
  // There is no artificial PiLine and no constant-PQ Load component.
  // The switched branch is excluded because its switch is open at t=0.
  //
  // SP::Ph1::Resistor and SP::Ph1::Inductor are solved directly by the
  // steady-state MNA solver. This avoids forcing passive R-L elements into
  // the Newton-Raphson bus-power model.
  // =========================================================================
  const String simNamePF = simName + "_PF";
  Logger::setLogDir("logs/" + simNamePF);

  auto nGfmPF = SimNode<Complex>::make("nGfm", PhaseType::Single);
  auto nBaseLoadSeriesPF =
      SimNode<Complex>::make("nBaseLoadSeries", PhaseType::Single);

  auto gfmInjectionPF =
      SP::Ph1::NetworkInjection::make("GFMInjectionPF", Logger::Level::debug);
  gfmInjectionPF->setParameters(Math::polar(p.lineToLineVoltageRms, 0.0),
                                p.nominalFrequency);

  auto baseLoadResistorPF =
      SP::Ph1::Resistor::make("BaseLoadResistorPF", Logger::Level::debug);
  baseLoadResistorPF->setParameters(baseLoadRL.resistance);

  auto baseLoadInductorPF =
      SP::Ph1::Inductor::make("BaseLoadInductorPF", Logger::Level::debug);
  baseLoadInductorPF->setParameters(baseLoadRL.inductance);

  gfmInjectionPF->connect({nGfmPF});
  baseLoadResistorPF->connect({nGfmPF, nBaseLoadSeriesPF});
  baseLoadInductorPF->connect({nBaseLoadSeriesPF, SimNode<Complex>::GND});

  auto systemPF = SystemTopology(
      p.nominalFrequency, SystemNodeList{nGfmPF, nBaseLoadSeriesPF},
      SystemComponentList{gfmInjectionPF, baseLoadResistorPF,
                          baseLoadInductorPF});

  auto loggerPF = DataLogger::make(simNamePF);
  loggerPF->logAttribute("Voltage_GFM", nGfmPF->attribute("v"));
  loggerPF->logAttribute("Voltage_BaseLoad_RL_Mid",
                         nBaseLoadSeriesPF->attribute("v"));
  loggerPF->logAttribute("Current_GFMInjection",
                         gfmInjectionPF->attribute("i_intf"));
  loggerPF->logAttribute("Current_BaseLoad_Resistor",
                         baseLoadResistorPF->attribute("i_intf"));
  loggerPF->logAttribute("Current_BaseLoad_Inductor",
                         baseLoadInductorPF->attribute("i_intf"));

  Simulation simPF(simNamePF, Logger::Level::debug);
  simPF.setSystem(systemPF);
  simPF.setTimeStep(1.0);
  simPF.setFinalTime(1.0);
  simPF.setDomain(Domain::SP);
  simPF.doInitFromNodesAndTerminals(true);
  simPF.addLogger(loggerPF);
  simPF.run();

  const Complex initialGfmVoltage = nGfmPF->singleVoltage();
  const Complex initialBaseSeriesVoltage = nBaseLoadSeriesPF->singleVoltage();

  // Current is explicitly defined from the GFM node through Rbase towards
  // the base-load internal node. This avoids relying on component current
  // sign conventions when deriving the generator-positive terminal power.
  const Complex initialBaseLoadCurrent =
      (initialGfmVoltage - initialBaseSeriesVoltage) / baseLoadRL.resistance;

  const Complex initialGeneratedPower =
      initialGfmVoltage * std::conj(initialBaseLoadCurrent);

  const Real initialGfmVoltageRms = std::abs(initialGfmVoltage);
  const Real initialBaseSeriesVoltageRms = std::abs(initialBaseSeriesVoltage);

  if (!(initialGfmVoltageRms > 0.0) || !(initialBaseSeriesVoltageRms > 0.0) ||
      !std::isfinite(initialGfmVoltageRms) ||
      !std::isfinite(initialBaseSeriesVoltageRms) ||
      !std::isfinite(initialGeneratedPower.real()) ||
      !std::isfinite(initialGeneratedPower.imag()) ||
      !(initialGeneratedPower.real() > 0.0)) {
    throw std::runtime_error("Invalid SP steady-state operating point.");
  }

  const Real powerTolerance = 1.0e-6 * p.ratedApparentPower;

  if (std::abs(initialGeneratedPower.real() - p.baseLoadActivePower) >
          powerTolerance ||
      std::abs(initialGeneratedPower.imag() - p.baseLoadReactivePower) >
          powerTolerance) {
    throw std::runtime_error(
        "SP R-L operating point does not match the requested base-load P/Q.");
  }

  SPDLOG_INFO("Solved SP operating point: P_gfm={} W, Q_gfm={} var, "
              "|V_gfm|={} V_LL RMS, |V_RL_mid|={} V_LL RMS",
              initialGeneratedPower.real(), initialGeneratedPower.imag(),
              initialGfmVoltageRms, initialBaseSeriesVoltageRms);

  // =========================================================================
  // EMT simulation
  //
  // Topology mirrors the SP initialization exactly for the energized branch:
  //
  //   GFM/PCC -- Rbase -- base internal node -- Lbase -- GND
  //       |
  //       switch -- Rstep -- step internal node -- Lstep -- GND
  // =========================================================================
  const String simNameEMT = simName + "_EMT";
  Logger::setLogDir("logs/" + simNameEMT);

  auto nGfmEMT = SimNode<Real>::make("nGfm", PhaseType::ABC);
  auto nBaseLoadSeriesEMT =
      SimNode<Real>::make("nBaseLoadSeries", PhaseType::ABC);
  auto nStepLoadEMT = SimNode<Real>::make("nStepLoad", PhaseType::ABC);
  auto nStepLoadSeriesEMT =
      SimNode<Real>::make("nStepLoadSeries", PhaseType::ABC);

  auto gfm =
      EMT::Ph3::GFM_Droop::make("GFM", "GFM", Logger::Level::debug, false);

  gfm->setParameters(p.nominalFrequency, p.phaseVoltagePeak,
                     initialGeneratedPower.real(),
                     initialGeneratedPower.imag());

  gfm->setDroopParameters(p.activePowerDroop, p.reactivePowerDroop,
                          p.voltageIntegralGain);

  gfm->setPowerFilterTimeConstant(p.powerMeasurementFilterTimeConstant);

  gfm->setControllerLimits(p.minimumFrequency, p.maximumFrequency,
                           p.minimumVoltagePeak, p.maximumVoltagePeak);

  gfm->setFilterParameters(p.filterInductance, p.filterCapacitance,
                           p.filterResistance, p.capacitorDampingResistance);

  gfm->withControl(true);

  auto baseLoadResistor =
      EMT::Ph3::Resistor::make("BaseLoadResistor", Logger::Level::debug);
  baseLoadResistor->setParameters(
      Math::singlePhaseParameterToThreePhase(baseLoadRL.resistance));

  auto baseLoadInductor =
      EMT::Ph3::Inductor::make("BaseLoadInductor", Logger::Level::debug);
  baseLoadInductor->setParameters(
      Math::singlePhaseParameterToThreePhase(baseLoadRL.inductance));

  auto stepLoadResistor =
      EMT::Ph3::Resistor::make("StepLoadResistor", Logger::Level::debug);
  stepLoadResistor->setParameters(
      Math::singlePhaseParameterToThreePhase(stepLoadRL.resistance));

  auto stepLoadInductor =
      EMT::Ph3::Inductor::make("StepLoadInductor", Logger::Level::debug);
  stepLoadInductor->setParameters(
      Math::singlePhaseParameterToThreePhase(stepLoadRL.inductance));

  auto loadSwitch = EMT::Ph3::Switch::make("LoadSwitch", Logger::Level::debug);
  loadSwitch->setParameters(
      Math::singlePhaseParameterToThreePhase(p.switchOpenResistance),
      Math::singlePhaseParameterToThreePhase(p.switchClosedResistance));
  loadSwitch->openSwitch();

  gfm->connect({nGfmEMT});

  baseLoadResistor->connect({nGfmEMT, nBaseLoadSeriesEMT});
  baseLoadInductor->connect({nBaseLoadSeriesEMT, SimNode<Real>::GND});

  loadSwitch->connect({nGfmEMT, nStepLoadEMT});
  stepLoadResistor->connect({nStepLoadEMT, nStepLoadSeriesEMT});
  stepLoadInductor->connect({nStepLoadSeriesEMT, SimNode<Real>::GND});

  auto systemEMT = SystemTopology(
      p.nominalFrequency,
      SystemNodeList{nGfmEMT, nBaseLoadSeriesEMT, nStepLoadEMT,
                     nStepLoadSeriesEMT},
      SystemComponentList{gfm, baseLoadResistor, baseLoadInductor, loadSwitch,
                          stepLoadResistor, stepLoadInductor});

  // Copy the solved SP phasors for the two identically named energized nodes:
  //   nGfm and nBaseLoadSeries.
  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  // The switched-load nodes have no SP counterparts and remain at zero,
  // consistent with the switch being open at t=0.

  // Generator-positive terminal metadata used by
  // GFM_Droop::initializeParentFromNodesAndTerminals().
  gfm->terminal(0)->setPower(initialGeneratedPower);

  // -------------------------------------------------------------------------
  // Logging
  // -------------------------------------------------------------------------
  auto loggerEMT = DataLogger::make(simNameEMT);

  loggerEMT->logAttribute("Voltage_PCC", nGfmEMT->attribute("v"));
  loggerEMT->logAttribute("Voltage_GFM", nGfmEMT->attribute("v"));
  loggerEMT->logAttribute("Voltage_BaseLoad_RL_Mid",
                          nBaseLoadSeriesEMT->attribute("v"));
  loggerEMT->logAttribute("Voltage_StepLoad", nStepLoadEMT->attribute("v"));
  loggerEMT->logAttribute("Voltage_StepLoad_RL_Mid",
                          nStepLoadSeriesEMT->attribute("v"));

  loggerEMT->logAttribute("Current_BaseLoad_Resistor",
                          baseLoadResistor->attribute("i_intf"));
  loggerEMT->logAttribute("Current_BaseLoad_Inductor",
                          baseLoadInductor->attribute("i_intf"));
  loggerEMT->logAttribute("Current_StepLoad_Resistor",
                          stepLoadResistor->attribute("i_intf"));
  loggerEMT->logAttribute("Current_StepLoad_Inductor",
                          stepLoadInductor->attribute("i_intf"));

  loggerEMT->logAttribute("Voltage_Source", gfm->attribute("Vs"));
  loggerEMT->logAttribute("Voltage_Command", gfm->attribute("Vsref"));
  loggerEMT->logAttribute("Current_GFM_InterfaceConsumer",
                          gfm->attribute("i_intf"));
  loggerEMT->logAttribute("Current_GFM_PCC_GeneratorPositive",
                          gfm->attribute("i_pcc"));

  loggerEMT->logAttribute("P_elec", gfm->attribute("P_elec"));
  loggerEMT->logAttribute("Q_elec", gfm->attribute("Q_elec"));
  loggerEMT->logAttribute("P_filtered", gfm->attribute("P_filtered"));
  loggerEMT->logAttribute("Q_filtered", gfm->attribute("Q_filtered"));
  loggerEMT->logAttribute("P_ref", gfm->attribute("P_ref"));
  loggerEMT->logAttribute("Q_ref", gfm->attribute("Q_ref"));

  loggerEMT->logAttribute("Frequency", gfm->attribute("frequency"));
  loggerEMT->logAttribute("Omega", gfm->attribute("omega"));
  loggerEMT->logAttribute("Theta", gfm->attribute("theta"));

  loggerEMT->logAttribute("Voltage_Magnitude", gfm->attribute("V_magnitude"));
  loggerEMT->logAttribute("Voltage_Droop_Output", gfm->attribute("V0"));
  loggerEMT->logAttribute("Voltage_Integrator",
                          gfm->attribute("voltage_integrator"));
  loggerEMT->logAttribute("Voltage_Command_Amplitude", gfm->attribute("V1"));

  // -------------------------------------------------------------------------
  // Switched load disturbance
  // -------------------------------------------------------------------------
  Simulation sim(simNameEMT, Logger::Level::debug);

  sim.addEvent(SwitchEvent3Ph::make(p.loadStepStartTime, loadSwitch, true));

  sim.addEvent(SwitchEvent3Ph::make(p.loadStepEndTime, loadSwitch, false));

  sim.setSystem(systemEMT);
  sim.doInitFromNodesAndTerminals(true);
  sim.setTimeStep(p.timeStep);
  sim.setFinalTime(p.finalTime + p.timeStep);
  sim.setDomain(Domain::EMT);
  sim.addLogger(loggerEMT);
  sim.run();

  return 0;
}
