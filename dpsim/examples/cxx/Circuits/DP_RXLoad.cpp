/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <DPsim.h>

using namespace DPsim;
using namespace CPS::DP;
using namespace CPS::DP::Ph3;

int main(int argc, char *argv[]) {

  Real timeStep = 0.0001;
  Real finalTime = 5;
  String simName = "DP_RXLoad";
  Logger::setLogDir("logs/" + simName);


  // Nodes
  auto n1 = SimNode::make("n1", PhaseType::ABC);
  auto n2 = SimNode::make("n2", PhaseType::ABC);

  // Components
  auto vs = Ph3::VoltageSource::make("vs");
  vs->setParameters(Complex(100, 0));
  auto r1 = Ph3::Resistor::make("r_1");
  Matrix r1_param = Matrix::Zero(3, 3);
  r1_param << 1., 1., 1., 1., 1., 1., 1., 1., 1.;
  r1->setParameters(r1_param);
  auto rx = Ph3::RXLoad::make("RXload");
  Matrix rx_param = Matrix::Zero(3, 3);
  rx_param << 10, 10, 10, 10, 10, 10, 10, 10, 10;
  rx->setParameters(rx_param, rx_param, 100, true);

  // Topology
  vs->connect({SimNode::GND, n1});
  r1->connect({n1, n2});
  rx->connect({n2, SimNode::GND});

  auto sys =
    SystemTopology(50,
    SystemNodeList{n1, n2},
    SystemComponentList{vs, r1, rx});

  // Logging
  auto logger = DataLogger::make(simName);
  logger->logAttribute("v1", n1->attribute("v"));
  logger->logAttribute("v2", n2->attribute("v"));
  logger->logAttribute("i12", r1->attribute("i_intf"));


  Simulation sim(simName);
  sim.setSystem(sys);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setDomain(CPS::Domain::DP);
  sim.addLogger(logger);

  sim.run();
}


