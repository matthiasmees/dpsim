// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <cmath>
#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_Switch.h>
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim-models/SP/SP_Ph1_Load.h>
#include <dpsim-models/SP/SP_Ph1_NetworkInjection.h>
#include <dpsim-models/SP/SP_Ph1_PiLine.h>
#include <dpsim-models/SimNode.h>
#include <dpsim-models/SystemTopology.h>
#include <dpsim/Simulation.h>
#include <memory>
#include <vector>

namespace SmallGrid {
using namespace DPsim;
using namespace CPS;

constexpr Real kFrequencyHz = 60.0;
constexpr Real kVoltageRmsLineLineV = 362215.5;
constexpr Real kVoltagePhaseAngleRad = -0.22235954952;
constexpr Real kActivePowerW = 125e6;
constexpr Real kReactivePowerVar = 80.881e6;
constexpr Real kLineResistanceOhm = 0.30;
constexpr Real kLineInductanceH = 0.02;
constexpr Real kLineCapacitanceF = 1e-6;
constexpr Real kLoadStepMW = 50.0;
constexpr Real kLoadStepSeconds = 2.5;

inline Matrix parameter3(Real value) {
  return Math::singlePhaseParameterToThreePhase(value);
}

inline std::vector<Complex> pccPhasors() {
  std::vector<Complex> result(3);
  for (int phase = 0; phase < 3; ++phase)
    result[phase] = std::polar(static_cast<double>(kVoltageRmsLineLineV),
                               kVoltagePhaseAngleRad - phase * 2.0 * PI / 3.0);
  return result;
}

inline std::vector<Complex> sourcePhasors() {
  const auto pcc = pccPhasors();
  const double scale = RMS3PH_TO_PEAK1PH;
  const Complex pccSpace = scale * pcc[0];
  const double phasePeak = std::abs(pccSpace);
  const Complex passiveSpace = -(2.0 / (3.0 * phasePeak * phasePeak)) *
                               Complex(kActivePowerW, -kReactivePowerVar) *
                               pccSpace;
  const Complex lineCurrent =
      passiveSpace +
      Complex(0.0, 2.0 * PI * kFrequencyHz * kLineCapacitanceF / 2.0) *
          pccSpace;
  const Complex impedance(kLineResistanceOhm,
                          2.0 * PI * kFrequencyHz * kLineInductanceH);
  std::vector<Complex> result(3);
  for (int phase = 0; phase < 3; ++phase) {
    const Complex pccPhase = scale * pcc[phase];
    const Complex current =
        lineCurrent * std::polar(1.0, -phase * 2.0 * PI / 3.0);
    // This is the series current from source to PCC. It combines the FMU's
    // passive terminal current with the PCC half-line capacitor current so
    // that the reset voltage is a steady-state PI-line operating point.
    result[phase] = (pccPhase + impedance * current) / scale;
  }
  return result;
}

inline std::vector<Complex> threePhasePhasors(Complex positiveSequence) {
  std::vector<Complex> result(3);
  for (int phase = 0; phase < 3; ++phase)
    result[phase] = positiveSequence * std::polar(1.0, -phase * 2.0 * PI / 3.0);
  return result;
}

inline Complex periodicInverterPowerInjection() {
  // The FMU's 5 ms response filter operates on abc current states. Its
  // 60 Hz periodic output therefore has a different P/Q than its instantaneous
  // current reference; use that periodic output in the PF equivalent.
  const double omegaTau = 2.0 * PI * kFrequencyHz * 5e-3;
  const double denominator = 1.0 + omegaTau * omegaTau;
  return Complex((kActivePowerW - omegaTau * kReactivePowerVar) / denominator,
                 (kReactivePowerVar + omegaTau * kActivePowerW) / denominator);
}

struct PowerflowParts {
  SP::SimNode::Ptr sourceNode;
  SP::SimNode::Ptr pccNode;
  std::unique_ptr<SystemTopology> system;
};

inline PowerflowParts makePowerflowParts() {
  PowerflowParts p;
  p.sourceNode = SP::SimNode::make("grid_source_bus", PhaseType::Single);
  p.pccNode = SP::SimNode::make("inverter_pcc", PhaseType::Single);
  auto source = SP::Ph1::NetworkInjection::make("grid_slack_pf");
  source->setParameters(kVoltageRmsLineLineV);
  source->setBaseVoltage(kVoltageRmsLineLineV);
  source->modifyPowerFlowBusType(PowerflowBusType::VD);
  auto line = SP::Ph1::PiLine::make("feeder_pi_line_pf");
  line->setParameters(kLineResistanceOhm, kLineInductanceH, kLineCapacitanceF);
  line->setBaseVoltage(kVoltageRmsLineLineV);
  auto inverter = SP::Ph1::Load::make("inverter_powerflow_injection");
  const Complex periodicPower = periodicInverterPowerInjection();
  inverter->setParameters(-periodicPower.real(), -periodicPower.imag(),
                          kVoltageRmsLineLineV);
  inverter->modifyPowerFlowBusType(PowerflowBusType::PQ);
  source->connect({p.sourceNode});
  line->connect({p.sourceNode, p.pccNode});
  inverter->connect({p.pccNode});
  p.system = std::make_unique<SystemTopology>(
      kFrequencyHz, TopologicalNode::List{p.sourceNode, p.pccNode},
      IdentifiedObject::List{source, line, inverter});
  return p;
}

inline Matrix initialLineCurrent() {
  const Complex pccSpace = RMS3PH_TO_PEAK1PH * pccPhasors()[0];
  const double phasePeak = std::abs(pccSpace);
  const Complex passiveSpace = -(2.0 / (3.0 * phasePeak * phasePeak)) *
                               Complex(kActivePowerW, -kReactivePowerVar) *
                               pccSpace;
  const Complex currentSpace =
      passiveSpace +
      Complex(0.0, 2.0 * PI * kFrequencyHz * kLineCapacitanceF / 2.0) *
          pccSpace;
  Matrix current(3, 1);
  for (int phase = 0; phase < 3; ++phase)
    current(phase, 0) =
        (currentSpace * std::polar(1.0, -phase * 2.0 * PI / 3.0)).real();
  return current;
}

struct NetworkParts {
  EMT::SimNode::Ptr sourceNode;
  EMT::SimNode::Ptr pccNode;
  EMT::SimNode::Ptr loadNode;
  std::shared_ptr<EMT::Ph3::VoltageSource> source;
  std::shared_ptr<EMT::Ph3::PiLine> line;
  std::shared_ptr<EMT::Ph3::Resistor> load;
  std::shared_ptr<EMT::Ph3::Switch> loadSwitch;
};

inline NetworkParts
makeNetworkParts(const std::vector<Complex> &sourceStart = sourcePhasors(),
                 const std::vector<Complex> &pccStart = pccPhasors()) {
  NetworkParts p;
  p.sourceNode =
      EMT::SimNode::make("grid_source_bus", PhaseType::ABC, sourceStart);
  p.pccNode = EMT::SimNode::make("inverter_pcc", PhaseType::ABC, pccStart);
  p.loadNode = EMT::SimNode::make("load_step_node", PhaseType::ABC, pccStart);

  p.source = EMT::Ph3::VoltageSource::make("grid_voltage_source");
  MatrixComp sourceReference(3, 1);
  const auto vSource = sourceStart;
  for (int phase = 0; phase < 3; ++phase)
    sourceReference(phase, 0) = vSource[phase];
  p.source->setParameters(sourceReference, kFrequencyHz);
  p.source->connect({EMT::SimNode::GND, p.sourceNode});

  p.line = EMT::Ph3::PiLine::make("feeder_pi_line");
  p.line->setParameters(parameter3(kLineResistanceOhm),
                        parameter3(kLineInductanceH),
                        parameter3(kLineCapacitanceF), parameter3(0.0));
  p.line->connect({p.sourceNode, p.pccNode});

  p.load = EMT::Ph3::Resistor::make("pcc_step_load");
  const Real loadResistance =
      kVoltageRmsLineLineV * kVoltageRmsLineLineV / (kLoadStepMW * 1e6);
  p.load->setParameters(parameter3(loadResistance));
  p.load->connect({EMT::SimNode::GND, p.loadNode});
  p.loadSwitch = EMT::Ph3::Switch::make("pcc_load_switch");
  p.loadSwitch->setParameters(parameter3(1e12), parameter3(1e-3), false);
  p.loadSwitch->connect({p.loadNode, p.pccNode});
  return p;
}

inline SystemTopology makeSystem(NetworkParts &p,
                                 const SimPowerComp<Real>::Ptr &inverter) {
  TopologicalNode::List nodes{p.sourceNode, p.pccNode, p.loadNode};
  IdentifiedObject::List components;
  components.push_back(p.source);
  components.push_back(p.line);
  components.push_back(inverter);
  components.push_back(p.load);
  components.push_back(p.loadSwitch);
  return SystemTopology(kFrequencyHz, nodes, components);
}
} // namespace SmallGrid
