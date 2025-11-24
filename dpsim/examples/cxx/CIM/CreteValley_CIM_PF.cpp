/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <DPsim.h>
#include <dpsim-models/CIM/Reader.h>

using namespace std;
using namespace DPsim;
using namespace CPS;
using namespace CPS::CIM;

/*
 * This example runs the powerflow for the CIGRE MV benchmark system (neglecting the tap changers of the transformers)
 */
int main(int argc, char **argv) {

  // Find CIM files
  std::list<fs::path> filenames;

  filenames = DPsim::Utils::findFiles(
      {"20240909T1759Z___DL_.xml", "20240909T1759Z__EQ_.xml",
       "20240909T1759Z___SV_.xml", "20240909T1759Z___TP_.xml",
       "20240909T1759Z___SSH_.xml", "20240909T1759Z___GL_.xml"},
      // "dpsim/dpsim/examples/cxx/CIM/CIM-Model","CIMPATH");
      "dpsim/Network_CIM_Data/Crete_equivalent_min_loading_activeOnly","CIMPATH");

  //  filenames = DPsim::Utils::findFiles(
  //    {"20250429T1349Z___DL_.xml", "20250429T1349Z___GL_.xml", "20250429T1349Z___SSH_.xml", "20250429T1349Z___SV_.xml", "20250429T1349Z___TP_.xml", "20250429T1349Z__EQ_.xml"},
  //    "dpsim/Network_CIM_Data/Crete2030","CIMPATH");


  String simName = "CreteValley_minActiveCPP";
  CPS::Real system_freq = 50;

  CIM::Reader reader(simName, Logger::Level::info, Logger::Level::info);
  SystemTopology system =
      reader.loadCIM(system_freq, filenames, CPS::Domain::SP, CPS::PhaseType::Single, CPS::GeneratorType::PVNode);

  auto logger = DPsim::DataLogger::make(simName);
  for (auto node : system.mNodes) {
    logger->logAttribute(node->name() + ".V", node->attribute("v"));
  }

  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(system);
  sim.setTimeStep(1);
  sim.setFinalTime(1440);
  sim.setDomain(Domain::SP);
  sim.setSolverType(Solver::Type::NRP);
  sim.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  sim.doInitFromNodesAndTerminals(true);
  sim.addLogger(logger);

  sim.run();

  return 0;
}
