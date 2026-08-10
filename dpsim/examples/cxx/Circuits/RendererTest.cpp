// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems,
// EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <DPsim.h>

#include <dpsim-models/EMT/EMT_Ph3_Capacitor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_SynchronGeneratorIdeal.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace CPS;

// =============================================================================
// Renderer test topology (all layout variants)
//
// GEN_1 -- BUS_GEN -- TR_1 -- BUS_HV -- LINE_1 -- BUS_REMOTE
//                              |                       |
//                           C_SHUNT                 R_SHUNT
//                              |                       |
//                             GND                     GND
//                                                   |
//                                                R_SERIES
//                                                   |
//                                                 BUS_R
//                                                   |
//                                                L_SERIES
//                                                   |
//                                               BUS_LOAD
//                                                   |
//                                                 R_LOAD
//                                                   |
//                                                  GND
//
// Tests:
//   Generator.svg
//   Transformer.svg
//   PiLine.svg
//   Resistor.svg
//   Inductor.svg
//   Capacitor.svg
//   BusBar.svg
//
// No simulation is executed.
// =============================================================================

int main() {

  constexpr Real frequency = 50.0;

  // ===========================================================================
  // Nodes
  //
  // Use the EMT SimNode alias explicitly.
  // CPS::EMT::SimNode == CPS::SimNode<Real>
  // ===========================================================================

  auto busGen = EMT::SimNode::make("BUS_GEN", PhaseType::ABC);

  auto busHv = EMT::SimNode::make("BUS_HV", PhaseType::ABC);

  auto busRemote = EMT::SimNode::make("BUS_REMOTE", PhaseType::ABC);

  auto busR = EMT::SimNode::make("BUS_R", PhaseType::ABC);

  auto busLoad = EMT::SimNode::make("BUS_LOAD", PhaseType::ABC);

  // ===========================================================================
  // Generator
  // ===========================================================================

  auto generator =
      EMT::Ph3::SynchronGeneratorIdeal::make("GEN_1", Logger::Level::off);

  // ===========================================================================
  // Transformer
  // ===========================================================================

  auto transformer =
      EMT::Ph3::Transformer::make("TR_1", "TR_1", Logger::Level::off, true);

  transformer->setParameters(20e3,  // nominal voltage LV [V]
                             110e3, // nominal voltage HV [V]
                             100e6, // rated apparent power [VA]
                             20e3 / 110e3,
                             0.0, // phase shift [rad]
                             Math::singlePhaseParameterToThreePhase(0.5),
                             Math::singlePhaseParameterToThreePhase(40e-3));

  // ===========================================================================
  // Pi line
  // ===========================================================================

  auto line = EMT::Ph3::PiLine::make("LINE_1", Logger::Level::off);

  line->setParameters(Math::singlePhaseParameterToThreePhase(1.0),    // R [Ohm]
                      Math::singlePhaseParameterToThreePhase(20e-3),  // L [H]
                      Math::singlePhaseParameterToThreePhase(2.0e-6), // C [F]
                      Math::singlePhaseParameterToThreePhase(1.0e-9)); // G [S]

  // ===========================================================================
  // Series resistor
  // ===========================================================================

  auto seriesResistor =
      EMT::Ph3::Resistor::make("R_SERIES", Logger::Level::off);

  seriesResistor->setParameters(Math::singlePhaseParameterToThreePhase(2.0));

  // ===========================================================================
  // Series inductor
  // ===========================================================================

  auto seriesInductor =
      EMT::Ph3::Inductor::make("L_SERIES", Logger::Level::off);

  seriesInductor->setParameters(Math::singlePhaseParameterToThreePhase(15e-3));

  // ===========================================================================
  // Shunt capacitor
  // ===========================================================================

  auto shuntCapacitor =
      EMT::Ph3::Capacitor::make("C_SHUNT", Logger::Level::off);

  shuntCapacitor->setParameters(Math::singlePhaseParameterToThreePhase(20e-6));

  // ===========================================================================
  // Shunt resistor
  // ===========================================================================

  auto shuntResistor = EMT::Ph3::Resistor::make("R_SHUNT", Logger::Level::off);

  shuntResistor->setParameters(Math::singlePhaseParameterToThreePhase(300.0));

  // ===========================================================================
  // Resistive load
  //
  // Deliberately represented by a resistor because there is currently no
  // dedicated Load.svg in the symbol set.
  // ===========================================================================

  auto load = EMT::Ph3::Resistor::make("R_LOAD", Logger::Level::off);

  load->setParameters(Math::singlePhaseParameterToThreePhase(100.0));

  // ===========================================================================
  // Connections
  // ===========================================================================

  generator->connect({busGen});

  transformer->connect({busGen, busHv});

  line->connect({busHv, busRemote});

  seriesResistor->connect({busRemote, busR});

  seriesInductor->connect({busR, busLoad});

  // Shunt capacitor
  shuntCapacitor->connect({busHv, EMT::SimNode::GND});

  // Shunt resistor
  shuntResistor->connect({busRemote, EMT::SimNode::GND});

  // Load
  load->connect({busLoad, EMT::SimNode::GND});

  // ===========================================================================
  // System topology
  // ===========================================================================

  SystemTopology systemEMT(
      frequency, DPsim::SystemNodeList{busGen, busHv, busRemote, busR, busLoad},
      DPsim::SystemComponentList{generator, transformer, line, seriesResistor,
                                 seriesInductor, shuntCapacitor, shuntResistor,
                                 load});

  // ===========================================================================
  // Renderer
  // ===========================================================================

  SystemTopologyRenderer::Options options;
  options.layout = SystemTopologyRenderer::Layout::LeftToRight;

  const std::filesystem::path symbolDirectory =
      "../dpsim-models/resources/Visuals";

  const std::filesystem::path outputFile = "logs/RendererTest/topology.svg";

  SystemTopologyRenderer renderer(systemEMT, symbolDirectory, options);

  renderer.renderSvg(outputFile);

  return 0;
}
