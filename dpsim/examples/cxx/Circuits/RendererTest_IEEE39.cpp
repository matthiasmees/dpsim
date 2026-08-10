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
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace DPsim;
using namespace CPS;

namespace {

constexpr Real kFrequency = 60.0;
constexpr Real kBasePower = 100.0e6;
constexpr Real kBaseVoltage = 345.0e3;
constexpr Real kOmega = 2.0 * PI * kFrequency;

struct BranchData {
  UInt from;
  UInt to;
  Real rPu;
  Real xPu;
  Real bPu;
  Real tapRatio;
};

struct LoadData {
  UInt bus;
  Real activePowerMW;
};

Real equivalentResistiveLoad(const Real lineToLineRmsVoltage,
                             const Real totalThreePhaseActivePower) {
  if (lineToLineRmsVoltage <= 0.0 || totalThreePhaseActivePower <= 0.0) {
    return 1.0e9;
  }
  return lineToLineRmsVoltage * lineToLineRmsVoltage /
         totalThreePhaseActivePower;
}

String branchName(const String &prefix, const UInt from, const UInt to) {
  return prefix + "_" + std::to_string(from) + "_" + std::to_string(to);
}

// Standard MATPOWER case39 branch topology.
// tapRatio == 0.0 -> Pi line
// tapRatio != 0.0 -> transformer
const std::array<BranchData, 46> kBranches = {
    {{1, 2, 0.0035, 0.0411, 0.6987, 0},   {1, 39, 0.001, 0.025, 0.75, 0},
     {2, 3, 0.0013, 0.0151, 0.2572, 0},   {2, 25, 0.007, 0.0086, 0.146, 0},
     {2, 30, 0, 0.0181, 0, 1.025},        {3, 4, 0.0013, 0.0213, 0.2214, 0},
     {3, 18, 0.0011, 0.0133, 0.2138, 0},  {4, 5, 0.0008, 0.0128, 0.1342, 0},
     {4, 14, 0.0008, 0.0129, 0.1382, 0},  {5, 6, 0.0002, 0.0026, 0.0434, 0},
     {5, 8, 0.0008, 0.0112, 0.1476, 0},   {6, 7, 0.0006, 0.0092, 0.113, 0},
     {6, 11, 0.0007, 0.0082, 0.1389, 0},  {6, 31, 0, 0.025, 0, 1.07},
     {7, 8, 0.0004, 0.0046, 0.078, 0},    {8, 9, 0.0023, 0.0363, 0.3804, 0},
     {9, 39, 0.001, 0.025, 1.2, 0},       {10, 11, 0.0004, 0.0043, 0.0729, 0},
     {10, 13, 0.0004, 0.0043, 0.0729, 0}, {10, 32, 0, 0.02, 0, 1.07},
     {12, 11, 0.0016, 0.0435, 0, 1.006},  {12, 13, 0.0016, 0.0435, 0, 1.006},
     {13, 14, 0.0009, 0.0101, 0.1723, 0}, {14, 15, 0.0018, 0.0217, 0.366, 0},
     {15, 16, 0.0009, 0.0094, 0.171, 0},  {16, 17, 0.0007, 0.0089, 0.1342, 0},
     {16, 19, 0.0016, 0.0195, 0.304, 0},  {16, 21, 0.0008, 0.0135, 0.2548, 0},
     {16, 24, 0.0003, 0.0059, 0.068, 0},  {17, 18, 0.0007, 0.0082, 0.1319, 0},
     {17, 27, 0.0013, 0.0173, 0.3216, 0}, {19, 20, 0.0007, 0.0138, 0, 1.06},
     {19, 33, 0.0007, 0.0142, 0, 1.07},   {20, 34, 0.0009, 0.018, 0, 1.009},
     {21, 22, 0.0008, 0.014, 0.2565, 0},  {22, 23, 0.0006, 0.0096, 0.1846, 0},
     {22, 35, 0, 0.0143, 0, 1.025},       {23, 24, 0.0022, 0.035, 0.361, 0},
     {23, 36, 0.0005, 0.0272, 0, 1},      {25, 26, 0.0032, 0.0323, 0.531, 0},
     {25, 37, 0.0006, 0.0232, 0, 1.025},  {26, 27, 0.0014, 0.0147, 0.2396, 0},
     {26, 28, 0.0043, 0.0474, 0.7802, 0}, {26, 29, 0.0057, 0.0625, 1.029, 0},
     {28, 29, 0.0014, 0.0151, 0.249, 0},  {29, 38, 0.0008, 0.0156, 0, 1.025}}};

const std::array<LoadData, 21> kLoads = {
    {{1, 97.6},   {3, 322},    {4, 500},  {7, 233.8}, {8, 522},  {9, 6.5},
     {12, 8.53},  {15, 320},   {16, 329}, {18, 158},  {20, 680}, {21, 274},
     {23, 247.5}, {24, 308.6}, {25, 224}, {26, 139},  {27, 281}, {28, 206},
     {29, 283.5}, {31, 9.2},   {39, 1104}}};

const std::array<UInt, 10> kGeneratorBuses = {
    {30, 31, 32, 33, 34, 35, 36, 37, 38, 39}};

} // namespace

int main() {

  // ===========================================================================
  // 39 EMT buses
  // ===========================================================================

  std::array<CPS::SimNode<Real>::Ptr, 40> buses{};
  DPsim::SystemNodeList systemNodes;
  systemNodes.reserve(39);

  for (UInt index = 1; index <= 39; ++index) {
    buses[index] =
        EMT::SimNode::make("BUS" + std::to_string(index), PhaseType::ABC);
    systemNodes.push_back(buses[index]);
  }

  // ===========================================================================
  // Components
  // ===========================================================================

  DPsim::SystemComponentList systemComponents;
  systemComponents.reserve(kBranches.size() + kLoads.size() +
                           kGeneratorBuses.size());

  const Real zBase = kBaseVoltage * kBaseVoltage / kBasePower;
  const Real yBase = 1.0 / zBase;

  UInt numberOfLines = 0;
  UInt numberOfTransformers = 0;

  for (const auto &branch : kBranches) {
    const Real resistance = branch.rPu * zBase;
    const Real inductance = branch.xPu * zBase / kOmega;

    if (branch.tapRatio == 0.0) {
      const Real capacitance = branch.bPu * yBase / kOmega;

      auto line = EMT::Ph3::PiLine::make(
          branchName("LINE", branch.from, branch.to), Logger::Level::off);

      line->setParameters(Math::singlePhaseParameterToThreePhase(resistance),
                          Math::singlePhaseParameterToThreePhase(inductance),
                          Math::singlePhaseParameterToThreePhase(capacitance),
                          Math::singlePhaseParameterToThreePhase(0.0));

      line->connect({buses[branch.from], buses[branch.to]});

      systemComponents.push_back(line);
      ++numberOfLines;

    } else {
      auto transformer = EMT::Ph3::Transformer::make(
          branchName("TR", branch.from, branch.to), Logger::Level::off);

      transformer->setParameters(
          kBaseVoltage, kBaseVoltage, kBasePower, branch.tapRatio, 0.0,
          Math::singlePhaseParameterToThreePhase(resistance),
          Math::singlePhaseParameterToThreePhase(inductance));

      transformer->connect({buses[branch.from], buses[branch.to]});

      systemComponents.push_back(transformer);
      ++numberOfTransformers;
    }
  }

  // ===========================================================================
  // Generators: standard case39 generator buses 30 ... 39
  // ===========================================================================

  for (const UInt bus : kGeneratorBuses) {
    auto generator = EMT::Ph3::SynchronGeneratorIdeal::make(
        "GEN" + std::to_string(bus), Logger::Level::off);

    generator->connect({buses[bus]});
    systemComponents.push_back(generator);
  }

  // ===========================================================================
  // Loads
  //
  // Topology-renderer-only constant-R equivalents, analogous to the IEEE-9
  // renderer test. Reactive power is intentionally not represented.
  // ===========================================================================

  for (const auto &loadData : kLoads) {
    auto load = EMT::Ph3::Resistor::make("LOAD" + std::to_string(loadData.bus),
                                         Logger::Level::off);

    load->setParameters(Math::singlePhaseParameterToThreePhase(
        equivalentResistiveLoad(kBaseVoltage, loadData.activePowerMW * 1.0e6)));

    load->connect({buses[loadData.bus], EMT::SimNode::GND});

    systemComponents.push_back(load);
  }

  // ===========================================================================
  // System topology
  // ===========================================================================

  SystemTopology systemEMT(kFrequency, systemNodes, systemComponents);

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
      "logs/RendererTest_IEEE39/topology.svg";

  SystemTopologyRenderer renderer(systemEMT, symbolDirectory, options);

  renderer.renderSvg(outputFile);

  return 0;
}
