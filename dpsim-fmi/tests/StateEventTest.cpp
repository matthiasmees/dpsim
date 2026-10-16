// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/Simulation.h>
#include <iostream>
#include <stdexcept>

using namespace CPS;
using namespace DPsim;
namespace {
void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}
} // namespace

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  try {
    const std::vector<FMI::ScalarInput> inputs{
        {"u_a", 100}, {"u_b", -50}, {"u_c", 25}};
    bool rejectedByDefault = false;
    try {
      FMI::ModelExchange rejected(argv[1], inputs);
    } catch (const std::exception &e) {
      rejectedByDefault =
          std::string(e.what()).find("explicitly enabled") != std::string::npos;
    }
    require(rejectedByDefault, "state-event FMU was accepted without opt-in");

    bool rejectedReset = false;
    try {
      FMI::ModelExchange noReset(argv[1], inputs,
                                 FMI::ModelExchangeOptions{true, false});
      Matrix crossing = noReset.description().initialState;
      crossing(0, 0) += 0.01;
      noReset.prepareStep(crossing, noReset.description().initialInputs, 1e-3);
      noReset.completeStep();
    } catch (const std::exception &e) {
      rejectedReset =
          std::string(e.what()).find("reinitialization") != std::string::npos;
    }
    require(rejectedReset,
            "continuous-state reset was accepted without reset opt-in");

    auto runtime = std::make_shared<FMI::ModelExchange>(
        argv[1], inputs, FMI::ModelExchangeOptions{true, true});
    auto branch = FMI::NonlinearSSN::make(
        "branch", runtime,
        FMI::ElectricalPort3Ph{
            {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, 1});
    auto source = VoltageBoundary3Ph::make("source", [](double) -> Matrix {
      Matrix values(3, 1);
      values << 100, -50, 25;
      return values;
    });
    auto bus = SimNode<Real>::make("bus", PhaseType::ABC);
    branch->connect({SimNode<Real>::GND, bus});
    source->connect({SimNode<Real>::GND, bus});
    SystemTopology system(50, {bus}, {source, branch});
    Simulation simulation("state_event_test", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(1e-3);
    simulation.doSteadyStateInit(false);
    simulation.start();
    simulation.step();

    const double expected = 6.28 + 0.01 - 2 * std::acos(-1.0);
    require(std::abs(branch->attributeTyped<Real>("fmi.state.theta")->get() -
                     expected) < 1e-11,
            "post-event state reset was not committed to the component");
    require(branch->attributeTyped<Real>("fmi.output.event_count")->get() == 1,
            "state event was not processed exactly once");
    for (int p = 0; p < 3; ++p) {
      const std::string phase(1, "abc"[p]);
      const double voltage = p == 0 ? 100 : (p == 1 ? -50 : 25);
      const double expectedCurrent =
          0.01 * voltage * (1 + 0.2 * std::sin(expected));
      require(std::abs(
                  branch->attributeTyped<Real>("fmi.output.i_" + phase)->get() -
                  expectedCurrent) < 1e-10,
              "post-event terminal current is inconsistent");
    }
    simulation.step();
    require(branch->attributeTyped<Real>("fmi.output.event_count")->get() == 1,
            "state event was processed more than once");
    require(std::abs(branch->attributeTyped<Real>("fmi.state.theta")->get() -
                     (expected + 0.01)) < 1e-11,
            "integration did not continue from the reset state");
    simulation.stop();
    std::cout << "PASS: endpoint state event, reset publication, continuous "
                 "terminal current and no duplicate event\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
