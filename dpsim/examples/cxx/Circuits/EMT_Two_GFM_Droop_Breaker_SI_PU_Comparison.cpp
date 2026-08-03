// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace DPsim;
using namespace CPS;

namespace GfmSiPuEquivalence {

// DataLogger writes with the stream default precision. For an equivalence test
// we need enough digits to detect conversion differences rather than comparing
// values rounded to six decimal places.
class HighPrecisionDataLogger : public DataLogger {
public:
  HighPrecisionDataLogger(String name, Bool enabled = true,
                          UInt downsampling = 1)
      : DataLogger(name, enabled, downsampling) {}

  void log(Real time, Int timeStepCount) override {
    if (!mEnabled || !(timeStepCount % mDownsampling == 0))
      return;

    if (mLogFile.tellp() == std::ofstream::pos_type(0)) {
      mLogFile << std::right << std::setw(24) << "time";
      for (const auto &entry : mAttributes)
        mLogFile << ", " << std::right << std::setw(23) << entry.first;
      mLogFile << '\n';
    }

    mLogFile << std::scientific << std::setprecision(17) << std::right
             << std::setw(24) << time;

    for (const auto &entry : mAttributes) {
      mLogFile << ", " << std::right << std::setw(23);

      if (auto realAttribute = std::dynamic_pointer_cast<CPS::Attribute<Real>>(
              entry.second.getPtr())) {
        mLogFile << realAttribute->get();
      } else if (auto intAttribute =
                     std::dynamic_pointer_cast<CPS::Attribute<Int>>(
                         entry.second.getPtr())) {
        mLogFile << intAttribute->get();
      } else {
        throw std::runtime_error(
            "HighPrecisionDataLogger received a non-scalar attribute: " +
            entry.first);
      }
    }

    mLogFile << '\n';
  }
};

// =============================================================================
// Parameters
// =============================================================================

enum class InputMode {
  SI,
  PerUnit,
};

const char *modeName(InputMode mode) {
  return mode == InputMode::SI ? "SI" : "PU";
}

struct SimulationParameters {
  String name = "EMT_Two_GFM_Droop_Breaker_SI_PU_Equivalence";
  Real frequency = 50.0;
  Real timeStep = 100e-6;
  Real finalTime = 10.0;
  Real breakerCloseTime = 2.0;

  Bool recomputeSystemMatrix = true;
  UInt comparisonDownsampling = 1;

  // Combined comparison criterion:
  //
  //   |a - b| <= absoluteTolerance
  //             + relativeTolerance * max(|a|, |b|)
  //
  // The two runs use the same physical model and differ only in which setter
  // interface is used. The default tolerances are therefore intentionally
  // strict.
  Real absoluteTolerance = 1e-7;
  Real relativeTolerance = 1e-9;
};

struct GfmParameters {
  Real ratedPower = 50e6;     // S_base [VA], total three-phase
  Real ratedVoltage = 10.5e3; // V_base [V line-to-line RMS]

  // Canonical scale-independent controller and filter parameters.
  Real filterInductiveReactancePu = 0.12981787824751215;
  Real filterCapacitiveSusceptancePu = 0.04561592533012379;
  Real filterResistancePu = 0.006887052341597797;
  Real capacitorDampingResistancePu = 0.5509641873278238;

  Real activePowerDroopPu = 0.05;
  Real reactivePowerDroopPu = 0.05;
  Real voltageIntegralGain = 40.0;                   // 1/s
  Real powerMeasurementFilterCutoffFrequency = 10.0; // Hz

  Real minimumFrequencyPu = 0.90;
  Real maximumFrequencyPu = 1.10;
  Real minimumVoltagePu = 0.70;
  Real maximumVoltagePu = 1.30;

  Real baseVoltagePhasePeak() const { return RMS3PH_TO_PEAK1PH * ratedVoltage; }

  Real baseImpedance() const {
    return ratedVoltage * ratedVoltage / ratedPower;
  }

  Real filterInductance(Real frequency) const {
    return filterInductiveReactancePu * baseImpedance() /
           (2.0 * PI * frequency);
  }

  Real filterCapacitance(Real frequency) const {
    return filterCapacitiveSusceptancePu /
           (2.0 * PI * frequency * baseImpedance());
  }

  Real filterResistance() const { return filterResistancePu * baseImpedance(); }

  Real capacitorDampingResistance() const {
    return capacitorDampingResistancePu * baseImpedance();
  }

  Real activePowerDroopSi(Real frequency) const {
    return activePowerDroopPu * frequency / ratedPower;
  }

  Real reactivePowerDroopSi() const {
    return reactivePowerDroopPu * baseVoltagePhasePeak() / ratedPower;
  }

  Real powerFilterTimeConstant() const {
    return 1.0 / (2.0 * PI * powerMeasurementFilterCutoffFrequency);
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

  Real inductance(Real /*frequency*/) const { return 0.0481260775594524 * 2.0; }
};

struct Line1Parameters {
  Real resistance = 5.04669;
  Real inductance = 0.13523123641;
  Real capacitance = 1.93522865e-6;
  Real conductance = 1e-15;
};

struct Line3Parameters {
  Real length = 5.0;

  Real resistance() const { return 0.0749 * length; }
  Real inductance() const { return 1.270693e-3 * length; }
  Real capacitance() const { return 0.00466961e-6 * length; }

  Real conductance = 1e-15;
};

struct BreakerParameters {
  Real openResistance = 1e12;
  Real closedResistance = 1e-3;
};

struct Parameters {
  SimulationParameters simulation;

  GfmParameters gasGfm;
  GfmParameters pshaGfm;

  GasTransformerParameters gasTransformer;
  PshTransformerParameters pshaTransformer;

  Line1Parameters line1;
  Line3Parameters line3;
  BreakerParameters breaker;

  Real pshaInitialAngleDegrees = 0.0;

  Parameters() {
    gasGfm.ratedPower = 50e6;
    gasGfm.ratedVoltage = 10.5e3;

    pshaGfm.ratedPower = 200e6;
    pshaGfm.ratedVoltage = 18.0e3;
  }
};

// =============================================================================
// Power-flow initialization
// =============================================================================

struct PowerFlowResult {
  SystemTopology system;

  Complex gasPower;
  Complex pshaPower;

  Complex gasBusVoltage;
  Complex pshaBusVoltage;

  Complex breakerGridSideVoltage;
  Complex breakerPshaSideVoltage;
};

Complex
generatorPositivePower(const std::shared_ptr<SP::Ph1::NetworkInjection> &source,
                       const SimNode<Complex>::Ptr &sourceBus) {
  const Complex voltage = sourceBus->singleVoltage();
  const Complex consumerCurrent = (**source->mIntfCurrent)(0, 0);
  return -voltage * std::conj(consumerCurrent);
}

PowerFlowResult buildAndRunPowerFlow(const Parameters &p) {
  const String simulationName = p.simulation.name + "_PF";

  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto busGas = SimNode<Complex>::make("BUS_gas", PhaseType::Single);
  auto busB = SimNode<Complex>::make("BUS_b", PhaseType::Single);
  auto busA = SimNode<Complex>::make("BUS_a", PhaseType::Single);
  auto busPshaGrid = SimNode<Complex>::make("BUS_psha_grid", PhaseType::Single);
  auto busPshaHv = SimNode<Complex>::make("BUS_psha_hv", PhaseType::Single);
  auto busPsha = SimNode<Complex>::make("BUS_psha", PhaseType::Single);

  auto gasInjection =
      SP::Ph1::NetworkInjection::make("GEN_gas_PF", Logger::Level::debug);
  gasInjection->setParameters(Math::polar(p.gasGfm.ratedVoltage, 0.0),
                              p.simulation.frequency);

  auto pshaInjection =
      SP::Ph1::NetworkInjection::make("GEN_psha_PF", Logger::Level::debug);
  pshaInjection->setParameters(
      Math::polar(p.pshaGfm.ratedVoltage,
                  p.pshaInitialAngleDegrees * PI / 180.0),
      p.simulation.frequency);

  auto gasTransformer = SP::Ph1::Transformer::make("TR_gas_PF", "TR_gas_PF",
                                                   Logger::Level::debug, true);
  gasTransformer->setParameters(
      p.gasTransformer.nominalVoltageLow, p.gasTransformer.nominalVoltageHigh,
      p.gasTransformer.ratedPower, p.gasTransformer.ratioMagnitude,
      p.gasTransformer.ratioPhase, p.gasTransformer.resistance(),
      p.gasTransformer.inductance(p.simulation.frequency));

  auto line1 = SP::Ph1::PiLine::make("LINE_1_PF", Logger::Level::debug);
  line1->setParameters(p.line1.resistance, p.line1.inductance,
                       p.line1.capacitance, p.line1.conductance);

  auto line3 = SP::Ph1::PiLine::make("LINE_3_PF", Logger::Level::debug);
  line3->setParameters(p.line3.resistance(), p.line3.inductance(),
                       p.line3.capacitance(), p.line3.conductance);

  auto breaker = SP::Ph1::Switch::make("BREAKER_psha_PF", Logger::Level::debug);
  breaker->setParameters(p.breaker.openResistance, p.breaker.closedResistance,
                         false);

  auto pshaTransformer = SP::Ph1::Transformer::make("TR_psha_PF", "TR_psha_PF",
                                                    Logger::Level::debug, true);
  pshaTransformer->setParameters(
      p.pshaTransformer.nominalVoltageLow, p.pshaTransformer.nominalVoltageHigh,
      p.pshaTransformer.ratedPower, p.pshaTransformer.ratioMagnitude,
      p.pshaTransformer.ratioPhase, p.pshaTransformer.resistance(),
      p.pshaTransformer.inductance(p.simulation.frequency));

  gasInjection->connect({busGas});
  gasTransformer->connect({busGas, busB});
  line1->connect({busB, busA});
  line3->connect({busA, busPshaGrid});
  breaker->connect({busPshaGrid, busPshaHv});
  pshaTransformer->connect({busPsha, busPshaHv});
  pshaInjection->connect({busPsha});

  SystemTopology system(p.simulation.frequency,
                        SystemNodeList{
                            busGas,
                            busB,
                            busA,
                            busPshaGrid,
                            busPshaHv,
                            busPsha,
                        },
                        SystemComponentList{
                            gasInjection,
                            gasTransformer,
                            line1,
                            line3,
                            breaker,
                            pshaTransformer,
                            pshaInjection,
                        });

  Simulation simulation(simulationName, Logger::Level::info);
  simulation.setSystem(system);
  simulation.setTimeStep(1.0);
  simulation.setFinalTime(1.0);
  simulation.setDomain(Domain::SP);
  simulation.setSolverType(Solver::Type::MNA);
  simulation.doInitFromNodesAndTerminals(true);
  simulation.run();

  const Complex gasPower = generatorPositivePower(gasInjection, busGas);
  const Complex pshaPower = generatorPositivePower(pshaInjection, busPsha);

  const Complex gasVoltage = busGas->singleVoltage();
  const Complex pshaVoltage = busPsha->singleVoltage();
  const Complex breakerGridVoltage = busPshaGrid->singleVoltage();
  const Complex breakerPshaVoltage = busPshaHv->singleVoltage();

  const auto finiteComplex = [](const Complex &value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  };

  if (!(std::abs(gasVoltage) > 0.0) || !(std::abs(pshaVoltage) > 0.0) ||
      !finiteComplex(gasVoltage) || !finiteComplex(pshaVoltage) ||
      !finiteComplex(gasPower) || !finiteComplex(pshaPower) ||
      !finiteComplex(breakerGridVoltage) ||
      !finiteComplex(breakerPshaVoltage)) {
    throw std::runtime_error(
        "Invalid SP steady-state result for SI/PU equivalence test");
  }

  SPDLOG_INFO("Solved common open-breaker SP initialization:"
              "\n  GEN_gas : P={} W, Q={} var, |V|={} V_LL RMS"
              "\n  GEN_psha: P={} W, Q={} var, |V|={} V_LL RMS"
              "\n  breaker mismatch: |dV|={} V_LL RMS, dAngle={} deg",
              gasPower.real(), gasPower.imag(), std::abs(gasVoltage),
              pshaPower.real(), pshaPower.imag(), std::abs(pshaVoltage),
              std::abs(breakerGridVoltage - breakerPshaVoltage),
              std::arg(breakerGridVoltage / breakerPshaVoltage) * 180.0 / PI);

  return {
      system,      gasPower,           pshaPower,          gasVoltage,
      pshaVoltage, breakerGridVoltage, breakerPshaVoltage,
  };
}

// =============================================================================
// GFM configuration through SI or PU setter interfaces
// =============================================================================

void configureGfm(const std::shared_ptr<EMT::Ph3::GFM_Droop> &gfm,
                  const GfmParameters &parameters, Real frequency,
                  const Complex &initialBusVoltage,
                  const Complex &initialGeneratorPower, InputMode mode) {
  const Real voltagePeakPhase = RMS3PH_TO_PEAK1PH * std::abs(initialBusVoltage);

  gfm->setBaseParameters(parameters.ratedPower, parameters.ratedVoltage,
                         frequency);

  if (mode == InputMode::SI) {
    gfm->setParameters(frequency, voltagePeakPhase,
                       initialGeneratorPower.real(),
                       initialGeneratorPower.imag());

    gfm->setDroopParameters(parameters.activePowerDroopSi(frequency),
                            parameters.reactivePowerDroopSi(),
                            parameters.voltageIntegralGain);

    gfm->setControllerLimits(
        parameters.minimumFrequencyPu * frequency,
        parameters.maximumFrequencyPu * frequency,
        parameters.minimumVoltagePu * parameters.baseVoltagePhasePeak(),
        parameters.maximumVoltagePu * parameters.baseVoltagePhasePeak());

    gfm->setFilterParameters(parameters.filterInductance(frequency),
                             parameters.filterCapacitance(frequency),
                             parameters.filterResistance(),
                             parameters.capacitorDampingResistance());
  } else {
    gfm->setParametersPerUnit(
        1.0, voltagePeakPhase / parameters.baseVoltagePhasePeak(),
        initialGeneratorPower.real() / parameters.ratedPower,
        initialGeneratorPower.imag() / parameters.ratedPower);

    gfm->setDroopParametersPerUnit(parameters.activePowerDroopPu,
                                   parameters.reactivePowerDroopPu,
                                   parameters.voltageIntegralGain);

    gfm->setControllerLimitsPerUnit(
        parameters.minimumFrequencyPu, parameters.maximumFrequencyPu,
        parameters.minimumVoltagePu, parameters.maximumVoltagePu);

    gfm->setFilterParametersPerUnit(parameters.filterInductiveReactancePu,
                                    parameters.filterCapacitiveSusceptancePu,
                                    parameters.filterResistancePu,
                                    parameters.capacitorDampingResistancePu);
  }

  gfm->setPowerFilterTimeConstant(parameters.powerFilterTimeConstant());
  gfm->withControl(true);

  SPDLOG_INFO(
      "Configured {} through {} setters:"
      "\n  P_ref={} pu, Q_ref={} pu"
      "\n  V_ref={} pu, f_ref=1 pu"
      "\n  k_p={} pu, k_q={} pu, k_iv={} 1/s"
      "\n  X_Lf={} pu, B_Cf={} pu, Rf={} pu, Rd={} pu",
      gfm->name(), modeName(mode),
      initialGeneratorPower.real() / parameters.ratedPower,
      initialGeneratorPower.imag() / parameters.ratedPower,
      voltagePeakPhase / parameters.baseVoltagePhasePeak(),
      parameters.activePowerDroopPu, parameters.reactivePowerDroopPu,
      parameters.voltageIntegralGain, parameters.filterInductiveReactancePu,
      parameters.filterCapacitiveSusceptancePu, parameters.filterResistancePu,
      parameters.capacitorDampingResistancePu);
}

// =============================================================================
// EMT run
// =============================================================================

std::filesystem::path runEmtVariant(const Parameters &p,
                                    const PowerFlowResult &powerFlow,
                                    InputMode mode) {
  const String simulationName =
      p.simulation.name + "_" + String(modeName(mode));

  std::filesystem::create_directories("logs/" + simulationName);
  Logger::setLogDir("logs/" + simulationName);

  auto busGas = SimNode<Real>::make("BUS_gas", PhaseType::ABC);
  auto busB = SimNode<Real>::make("BUS_b", PhaseType::ABC);
  auto busA = SimNode<Real>::make("BUS_a", PhaseType::ABC);
  auto busPshaGrid = SimNode<Real>::make("BUS_psha_grid", PhaseType::ABC);
  auto busPshaHv = SimNode<Real>::make("BUS_psha_hv", PhaseType::ABC);
  auto busPsha = SimNode<Real>::make("BUS_psha", PhaseType::ABC);

  auto gasGfm = EMT::Ph3::GFM_Droop::make("GEN_gas", "GEN_gas",
                                          Logger::Level::debug, false);
  auto pshaGfm = EMT::Ph3::GFM_Droop::make("GEN_psha", "GEN_psha",
                                           Logger::Level::debug, false);

  configureGfm(gasGfm, p.gasGfm, p.simulation.frequency,
               powerFlow.gasBusVoltage, powerFlow.gasPower, mode);

  configureGfm(pshaGfm, p.pshaGfm, p.simulation.frequency,
               powerFlow.pshaBusVoltage, powerFlow.pshaPower, mode);

  auto gasTransformer = EMT::Ph3::Transformer::make("TR_gas", "TR_gas",
                                                    Logger::Level::debug, true);
  gasTransformer->setParameters(
      p.gasTransformer.nominalVoltageLow, p.gasTransformer.nominalVoltageHigh,
      p.gasTransformer.ratedPower, p.gasTransformer.ratioMagnitude,
      p.gasTransformer.ratioPhase,
      Math::singlePhaseParameterToThreePhase(p.gasTransformer.resistance()),
      Math::singlePhaseParameterToThreePhase(
          p.gasTransformer.inductance(p.simulation.frequency)));

  auto line1 = EMT::Ph3::PiLine::make("LINE_1", Logger::Level::debug);
  line1->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line1.resistance),
      Math::singlePhaseParameterToThreePhase(p.line1.inductance),
      Math::singlePhaseParameterToThreePhase(p.line1.capacitance),
      Math::singlePhaseParameterToThreePhase(p.line1.conductance));

  auto line3 = EMT::Ph3::PiLine::make("LINE_3", Logger::Level::debug);
  line3->setParameters(
      Math::singlePhaseParameterToThreePhase(p.line3.resistance()),
      Math::singlePhaseParameterToThreePhase(p.line3.inductance()),
      Math::singlePhaseParameterToThreePhase(p.line3.capacitance()),
      Math::singlePhaseParameterToThreePhase(p.line3.conductance));

  auto breaker = EMT::Ph3::Switch::make("BREAKER_psha", Logger::Level::debug);
  breaker->setParameters(
      Math::singlePhaseParameterToThreePhase(p.breaker.openResistance),
      Math::singlePhaseParameterToThreePhase(p.breaker.closedResistance));
  breaker->openSwitch();

  auto pshaTransformer = EMT::Ph3::Transformer::make(
      "TR_psha", "TR_psha", Logger::Level::debug, true);
  pshaTransformer->setParameters(
      p.pshaTransformer.nominalVoltageLow, p.pshaTransformer.nominalVoltageHigh,
      p.pshaTransformer.ratedPower, p.pshaTransformer.ratioMagnitude,
      p.pshaTransformer.ratioPhase,
      Math::singlePhaseParameterToThreePhase(p.pshaTransformer.resistance()),
      Math::singlePhaseParameterToThreePhase(
          p.pshaTransformer.inductance(p.simulation.frequency)));

  gasGfm->connect({busGas});
  gasTransformer->connect({busGas, busB});
  line1->connect({busB, busA});
  line3->connect({busA, busPshaGrid});
  breaker->connect({busPshaGrid, busPshaHv});
  pshaTransformer->connect({busPsha, busPshaHv});
  pshaGfm->connect({busPsha});

  SystemTopology system(p.simulation.frequency,
                        SystemNodeList{
                            busGas,
                            busB,
                            busA,
                            busPshaGrid,
                            busPshaHv,
                            busPsha,
                        },
                        SystemComponentList{
                            gasGfm,
                            gasTransformer,
                            line1,
                            line3,
                            breaker,
                            pshaTransformer,
                            pshaGfm,
                        });

  system.initWithPowerflow(powerFlow.system, Domain::EMT);

  gasGfm->terminal(0)->setPower(powerFlow.gasPower);
  pshaGfm->terminal(0)->setPower(powerFlow.pshaPower);

  auto logger = std::make_shared<HighPrecisionDataLogger>(
      simulationName, true, p.simulation.comparisonDownsampling);

  // Physical EMT network quantities.
  logger->logAttribute("BUS_gas_v", busGas->attribute("v"));
  logger->logAttribute("BUS_b_v", busB->attribute("v"));
  logger->logAttribute("BUS_a_v", busA->attribute("v"));
  logger->logAttribute("BUS_psha_grid_v", busPshaGrid->attribute("v"));
  logger->logAttribute("BUS_psha_hv_v", busPshaHv->attribute("v"));
  logger->logAttribute("BUS_psha_v", busPsha->attribute("v"));

  logger->logAttribute("LINE_1_i", line1->attribute("i_intf"));
  logger->logAttribute("LINE_3_i", line3->attribute("i_intf"));
  logger->logAttribute("BREAKER_psha_i", breaker->attribute("i_intf"));
  logger->logAttribute("BREAKER_psha_v", breaker->attribute("v_intf"));

  // Gas GFM: canonical PU states and SI mirrors.
  logger->logAttribute("GEN_gas_P_elec_pu", gasGfm->attribute("P_elec_pu"));
  logger->logAttribute("GEN_gas_Q_elec_pu", gasGfm->attribute("Q_elec_pu"));
  logger->logAttribute("GEN_gas_P_filtered_pu",
                       gasGfm->attribute("P_filtered_pu"));
  logger->logAttribute("GEN_gas_Q_filtered_pu",
                       gasGfm->attribute("Q_filtered_pu"));
  logger->logAttribute("GEN_gas_frequency_pu",
                       gasGfm->attribute("frequency_pu"));
  logger->logAttribute("GEN_gas_voltage_magnitude_pu",
                       gasGfm->attribute("V_magnitude_pu"));
  logger->logAttribute("GEN_gas_voltage_command_pu",
                       gasGfm->attribute("V1_pu"));
  logger->logAttribute("GEN_gas_i_pcc_pu", gasGfm->attribute("i_pcc_pu"));
  logger->logAttribute("GEN_gas_v_source_ref_pu",
                       gasGfm->attribute("Vsref_pu"));
  logger->logAttribute("GEN_gas_theta", gasGfm->attribute("theta"));
  logger->logAttribute("GEN_gas_P_elec", gasGfm->attribute("P_elec"));
  logger->logAttribute("GEN_gas_Q_elec", gasGfm->attribute("Q_elec"));
  logger->logAttribute("GEN_gas_frequency", gasGfm->attribute("frequency"));
  logger->logAttribute("GEN_gas_voltage_magnitude",
                       gasGfm->attribute("V_magnitude"));

  // PSHA GFM: canonical PU states and SI mirrors.
  logger->logAttribute("GEN_psha_P_elec_pu", pshaGfm->attribute("P_elec_pu"));
  logger->logAttribute("GEN_psha_Q_elec_pu", pshaGfm->attribute("Q_elec_pu"));
  logger->logAttribute("GEN_psha_P_filtered_pu",
                       pshaGfm->attribute("P_filtered_pu"));
  logger->logAttribute("GEN_psha_Q_filtered_pu",
                       pshaGfm->attribute("Q_filtered_pu"));
  logger->logAttribute("GEN_psha_frequency_pu",
                       pshaGfm->attribute("frequency_pu"));
  logger->logAttribute("GEN_psha_voltage_magnitude_pu",
                       pshaGfm->attribute("V_magnitude_pu"));
  logger->logAttribute("GEN_psha_voltage_command_pu",
                       pshaGfm->attribute("V1_pu"));
  logger->logAttribute("GEN_psha_i_pcc_pu", pshaGfm->attribute("i_pcc_pu"));
  logger->logAttribute("GEN_psha_v_source_ref_pu",
                       pshaGfm->attribute("Vsref_pu"));
  logger->logAttribute("GEN_psha_theta", pshaGfm->attribute("theta"));
  logger->logAttribute("GEN_psha_P_elec", pshaGfm->attribute("P_elec"));
  logger->logAttribute("GEN_psha_Q_elec", pshaGfm->attribute("Q_elec"));
  logger->logAttribute("GEN_psha_frequency", pshaGfm->attribute("frequency"));
  logger->logAttribute("GEN_psha_voltage_magnitude",
                       pshaGfm->attribute("V_magnitude"));

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

  return std::filesystem::path("logs") / simulationName /
         (simulationName + ".csv");
}

// =============================================================================
// CSV comparison and assertion test
// =============================================================================

String trim(String value) {
  const auto notSpace = [](unsigned char c) { return !std::isspace(c); };

  value.erase(value.begin(),
              std::find_if(value.begin(), value.end(), notSpace));
  value.erase(std::find_if(value.rbegin(), value.rend(), notSpace).base(),
              value.end());
  return value;
}

std::vector<String> splitCsvLine(const String &line) {
  std::vector<String> fields;
  std::stringstream stream(line);
  String field;

  while (std::getline(stream, field, ','))
    fields.push_back(trim(field));

  return fields;
}

struct ComparisonSummary {
  std::size_t comparedRows = 0;
  std::size_t comparedValues = 0;

  Real maximumAbsoluteError = 0.0;
  Real maximumRelativeError = 0.0;

  Real worstTime = 0.0;
  String worstColumn;
  Real worstSiValue = 0.0;
  Real worstPuValue = 0.0;
};

ComparisonSummary compareCsvFiles(const std::filesystem::path &siFile,
                                  const std::filesystem::path &puFile,
                                  Real absoluteTolerance,
                                  Real relativeTolerance) {
  std::ifstream siStream(siFile);
  std::ifstream puStream(puFile);

  if (!siStream.is_open())
    throw std::runtime_error("Cannot open SI result file: " + siFile.string());
  if (!puStream.is_open())
    throw std::runtime_error("Cannot open PU result file: " + puFile.string());

  String siLine;
  String puLine;

  if (!std::getline(siStream, siLine) || !std::getline(puStream, puLine))
    throw std::runtime_error("Missing CSV header in SI/PU result files");

  const auto siHeader = splitCsvLine(siLine);
  const auto puHeader = splitCsvLine(puLine);

  if (siHeader != puHeader) {
    throw std::runtime_error("SI and PU CSV headers differ; the two runs did "
                             "not log the same signals");
  }

  if (siHeader.empty() || siHeader.front() != "time")
    throw std::runtime_error("Unexpected comparison CSV header");

  ComparisonSummary summary;
  std::size_t lineNumber = 1;

  while (true) {
    const Bool hasSiLine = static_cast<Bool>(std::getline(siStream, siLine));
    const Bool hasPuLine = static_cast<Bool>(std::getline(puStream, puLine));

    if (!hasSiLine && !hasPuLine)
      break;

    ++lineNumber;

    if (hasSiLine != hasPuLine) {
      throw std::runtime_error(
          "SI and PU CSV files contain different numbers of rows");
    }

    const auto siFields = splitCsvLine(siLine);
    const auto puFields = splitCsvLine(puLine);

    if (siFields.size() != siHeader.size() ||
        puFields.size() != puHeader.size()) {
      throw std::runtime_error("Malformed CSV row at line " +
                               std::to_string(lineNumber));
    }

    const Real siTime = std::stod(siFields[0]);
    const Real puTime = std::stod(puFields[0]);
    const Real timeError = std::abs(siTime - puTime);

    if (timeError > 1e-12) {
      throw std::runtime_error("SI and PU time columns differ at CSV line " +
                               std::to_string(lineNumber));
    }

    for (std::size_t column = 1; column < siHeader.size(); ++column) {
      const Real siValue = std::stod(siFields[column]);
      const Real puValue = std::stod(puFields[column]);

      if (!std::isfinite(siValue) || !std::isfinite(puValue)) {
        throw std::runtime_error("Non-finite comparison value in column " +
                                 siHeader[column] +
                                 " at t=" + std::to_string(siTime));
      }

      const Real absoluteError = std::abs(siValue - puValue);
      const Real scale = std::max(std::abs(siValue), std::abs(puValue));
      const Real relativeError = scale > std::numeric_limits<Real>::min()
                                     ? absoluteError / scale
                                     : absoluteError;

      if (absoluteError > summary.maximumAbsoluteError) {
        summary.maximumAbsoluteError = absoluteError;
        summary.worstTime = siTime;
        summary.worstColumn = siHeader[column];
        summary.worstSiValue = siValue;
        summary.worstPuValue = puValue;
      }

      summary.maximumRelativeError =
          std::max(summary.maximumRelativeError, relativeError);

      const Real allowedError = absoluteTolerance + relativeTolerance * scale;

      if (absoluteError > allowedError) {
        std::ostringstream message;
        message << std::setprecision(17) << "SI/PU equivalence assertion failed"
                << "\n  CSV line: " << lineNumber << "\n  time: " << siTime
                << " s"
                << "\n  column: " << siHeader[column]
                << "\n  SI value: " << siValue << "\n  PU value: " << puValue
                << "\n  absolute error: " << absoluteError
                << "\n  allowed error: " << allowedError
                << "\n  absolute tolerance: " << absoluteTolerance
                << "\n  relative tolerance: " << relativeTolerance;
        throw std::runtime_error(message.str());
      }

      ++summary.comparedValues;
    }

    ++summary.comparedRows;
  }

  if (summary.comparedRows == 0)
    throw std::runtime_error("No SI/PU result rows were compared");

  return summary;
}

void validateParameters(const Parameters &parameters) {
  if (!(parameters.simulation.timeStep > 0.0) ||
      !(parameters.simulation.finalTime > 0.0) ||
      !(parameters.simulation.breakerCloseTime >= 0.0) ||
      !(parameters.simulation.breakerCloseTime <
        parameters.simulation.finalTime) ||
      !(parameters.simulation.comparisonDownsampling > 0) ||
      !(parameters.simulation.absoluteTolerance >= 0.0) ||
      !(parameters.simulation.relativeTolerance >= 0.0)) {
    throw std::invalid_argument(
        "Require dt>0, finalTime>0, 0<=breakerCloseTime<finalTime, "
        "downsampling>0, and non-negative tolerances");
  }
}

} // namespace GfmSiPuEquivalence

int main(int argc, char *argv[]) {
  try {
    GfmSiPuEquivalence::Parameters parameters;
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

      if (args.options.find("PSHA_INITIAL_ANGLE_DEG") != args.options.end()) {
        parameters.pshaInitialAngleDegrees =
            args.getOptionReal("PSHA_INITIAL_ANGLE_DEG");
      }

      if (args.options.find("ABS_TOL") != args.options.end()) {
        parameters.simulation.absoluteTolerance = args.getOptionReal("ABS_TOL");
      }

      if (args.options.find("REL_TOL") != args.options.end()) {
        parameters.simulation.relativeTolerance = args.getOptionReal("REL_TOL");
      }

      if (args.options.find("COMPARE_DOWNSAMPLING") != args.options.end()) {
        const Real downsampling = args.getOptionReal("COMPARE_DOWNSAMPLING");
        if (!(downsampling >= 1.0) ||
            std::floor(downsampling) != downsampling ||
            downsampling >
                static_cast<Real>(std::numeric_limits<UInt>::max())) {
          throw std::invalid_argument(
              "COMPARE_DOWNSAMPLING must be an integer >= 1");
        }
        parameters.simulation.comparisonDownsampling =
            static_cast<UInt>(downsampling);
      }
    }

    GfmSiPuEquivalence::validateParameters(parameters);

    const auto powerFlow = GfmSiPuEquivalence::buildAndRunPowerFlow(parameters);

    const auto siFile = GfmSiPuEquivalence::runEmtVariant(
        parameters, powerFlow, GfmSiPuEquivalence::InputMode::SI);

    const auto puFile = GfmSiPuEquivalence::runEmtVariant(
        parameters, powerFlow, GfmSiPuEquivalence::InputMode::PerUnit);

    const auto summary = GfmSiPuEquivalence::compareCsvFiles(
        siFile, puFile, parameters.simulation.absoluteTolerance,
        parameters.simulation.relativeTolerance);

    SPDLOG_INFO("SI/PU equivalence test PASSED:"
                "\n  SI file: {}"
                "\n  PU file: {}"
                "\n  compared rows: {}"
                "\n  compared values: {}"
                "\n  maximum absolute error: {}"
                "\n  maximum relative error: {}"
                "\n  worst column: {}"
                "\n  worst time: {} s"
                "\n  SI value: {}"
                "\n  PU value: {}",
                siFile.string(), puFile.string(), summary.comparedRows,
                summary.comparedValues, summary.maximumAbsoluteError,
                summary.maximumRelativeError, summary.worstColumn,
                summary.worstTime, summary.worstSiValue, summary.worstPuValue);

    return EXIT_SUCCESS;
  } catch (const std::exception &exception) {
    SPDLOG_ERROR("SI/PU equivalence test FAILED:\n{}", exception.what());
    return EXIT_FAILURE;
  }
}
