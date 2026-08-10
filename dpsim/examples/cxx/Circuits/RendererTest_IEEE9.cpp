// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <DPsim.h>

#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorIdeal.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace DPsim;
using namespace CPS;

namespace {

Real equivalentResistiveLoad(const Real lineToLineRmsVoltage,
                             const Real totalThreePhaseActivePower) {

  if (lineToLineRmsVoltage <= 0.0 || totalThreePhaseActivePower <= 0.0) {

    return 1.0e9;
  }

  return lineToLineRmsVoltage * lineToLineRmsVoltage /
         totalThreePhaseActivePower;
}

} // namespace

int main() {

  // ===========================================================================
  // IEEE-9 data
  // ===========================================================================

  constexpr Real requestedFrequency = 60.0;

  CPS::CIM::Examples::Grids::IEEE9::ScenarioConfig ieee9(requestedFrequency);

  // ===========================================================================
  // EMT nodes
  // ===========================================================================

  auto n1 = EMT::SimNode::make("BUS1", PhaseType::ABC);

  auto n2 = EMT::SimNode::make("BUS2", PhaseType::ABC);

  auto n3 = EMT::SimNode::make("BUS3", PhaseType::ABC);

  auto n4 = EMT::SimNode::make("BUS4", PhaseType::ABC);

  auto n5 = EMT::SimNode::make("BUS5", PhaseType::ABC);

  auto n6 = EMT::SimNode::make("BUS6", PhaseType::ABC);

  auto n7 = EMT::SimNode::make("BUS7", PhaseType::ABC);

  auto n8 = EMT::SimNode::make("BUS8", PhaseType::ABC);

  auto n9 = EMT::SimNode::make("BUS9", PhaseType::ABC);

  // ===========================================================================
  // Generators
  //
  // Ideal synchronous generators are sufficient here because no simulation is
  // executed. Their class name is still mapped to generator_h/v.svg by the
  // topology renderer.
  // ===========================================================================

  auto gen1 = EMT::Ph3::SynchronGeneratorIdeal::make(ieee9.gen1.Name,
                                                     Logger::Level::off);

  auto gen2 = EMT::Ph3::SynchronGeneratorIdeal::make(ieee9.gen2.Name,
                                                     Logger::Level::off);

  auto gen3 = EMT::Ph3::SynchronGeneratorIdeal::make(ieee9.gen3.Name,
                                                     Logger::Level::off);

  // ===========================================================================
  // Transformers
  // ===========================================================================

  auto transf14 =
      EMT::Ph3::Transformer::make(ieee9.transf14.Name, Logger::Level::off);

  transf14->setParameters(
      ieee9.transf14.VoltageLVSide, ieee9.transf14.VoltageHVSide,
      ieee9.transf14.RatedPower, ieee9.transf14.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf14.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf14.Inductance));

  auto transf27 =
      EMT::Ph3::Transformer::make(ieee9.transf27.Name, Logger::Level::off);

  transf27->setParameters(
      ieee9.transf27.VoltageLVSide, ieee9.transf27.VoltageHVSide,
      ieee9.transf27.RatedPower, ieee9.transf27.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf27.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf27.Inductance));

  auto transf39 =
      EMT::Ph3::Transformer::make(ieee9.transf39.Name, Logger::Level::off);

  transf39->setParameters(
      ieee9.transf39.VoltageLVSide, ieee9.transf39.VoltageHVSide,
      ieee9.transf39.RatedPower, ieee9.transf39.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf39.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf39.Inductance));

  // ===========================================================================
  // Transmission lines
  // ===========================================================================

  auto line54 = EMT::Ph3::PiLine::make(ieee9.line54.Name, Logger::Level::off);

  line54->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Conductance));

  auto line64 = EMT::Ph3::PiLine::make(ieee9.line64.Name, Logger::Level::off);

  line64->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Conductance));

  auto line75 = EMT::Ph3::PiLine::make(ieee9.line75.Name, Logger::Level::off);

  line75->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Conductance));

  auto line96 = EMT::Ph3::PiLine::make(ieee9.line96.Name, Logger::Level::off);

  line96->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Conductance));

  auto line78 = EMT::Ph3::PiLine::make(ieee9.line78.Name, Logger::Level::off);

  line78->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Conductance));

  auto line89 = EMT::Ph3::PiLine::make(ieee9.line89.Name, Logger::Level::off);

  line89->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Conductance));

  // ===========================================================================
  // Loads
  //
  // Renderer-only resistor equivalents.
  // ===========================================================================

  auto load5 = EMT::Ph3::Resistor::make("LOAD5", Logger::Level::off);

  load5->setParameters(Math::singlePhaseParameterToThreePhase(
      equivalentResistiveLoad(ieee9.load5.BaseVoltage, ieee9.load5.RealPower)));

  auto load6 = EMT::Ph3::Resistor::make("LOAD6", Logger::Level::off);

  load6->setParameters(Math::singlePhaseParameterToThreePhase(
      equivalentResistiveLoad(ieee9.load6.BaseVoltage, ieee9.load6.RealPower)));

  auto load8 = EMT::Ph3::Resistor::make("LOAD8", Logger::Level::off);

  load8->setParameters(Math::singlePhaseParameterToThreePhase(
      equivalentResistiveLoad(ieee9.load8.BaseVoltage, ieee9.load8.RealPower)));

  // ===========================================================================
  // Connections
  //
  // These are the same IEEE-9 connections used in the existing DPsim examples.
  // ===========================================================================

  gen1->connect({n1});

  gen2->connect({n2});

  gen3->connect({n3});

  transf14->connect({n1, n4});

  transf27->connect({n2, n7});

  transf39->connect({n3, n9});

  line54->connect({n5, n4});

  line64->connect({n6, n4});

  line75->connect({n7, n5});

  line96->connect({n9, n6});

  line78->connect({n7, n8});

  line89->connect({n8, n9});

  load5->connect({n5, EMT::SimNode::GND});

  load6->connect({n6, EMT::SimNode::GND});

  load8->connect({n8, EMT::SimNode::GND});

  // ===========================================================================
  // System topology
  // ===========================================================================

  SystemTopology systemEMT(
      ieee9.nomFreq, SystemNodeList{n1, n2, n3, n4, n5, n6, n7, n8, n9},
      SystemComponentList{gen1, gen2, gen3,

                          transf14, transf27, transf39,

                          line54, line64, line75, line96, line78, line89,

                          load5, load6, load8});

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

  const std::filesystem::path outputFile =
      "logs/RendererTest_IEEE9/topology.svg";

  SystemTopologyRenderer renderer(systemEMT, symbolDirectory, options);

  renderer.renderSvg(outputFile);

  return 0;
}
