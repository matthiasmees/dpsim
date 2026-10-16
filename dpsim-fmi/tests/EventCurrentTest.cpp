// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/Simulation.h>
#include <iostream>
int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  using namespace CPS;
  using namespace DPsim;
  auto runtime = std::make_shared<FMI::ModelExchange>(
      argv[1],
      std::vector<FMI::ScalarInput>{{"u_a", 7}, {"u_b", -11}, {"u_c", 13}});
  auto component = FMI::NonlinearSSN::make(
      "branch", runtime,
      FMI::ElectricalPort3Ph{{"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, 1});
  const Matrix initial = runtime->description().initialState;
  auto source = VoltageBoundary3Ph::make("source", [](double) -> Matrix {
    Matrix v(3, 1);
    v << 7, -11, 13;
    return v;
  });
  auto bus = SimNode<Real>::make("bus", PhaseType::ABC);
  component->connect({SimNode<Real>::GND, bus});
  source->connect({SimNode<Real>::GND, bus});
  SystemTopology system(50, {bus}, {source, component});
  Simulation simulation("event_current_rejection", Logger::Level::off);
  simulation.setDomain(Domain::EMT);
  simulation.setSystem(system);
  simulation.setTimeStep(2e-4);
  simulation.doSteadyStateInit(false);
  simulation.start();
  bool rejected = false;
  try {
    simulation.step();
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    rejected =
        std::string(e.what()).find("instantaneous algebraic event iteration") !=
        std::string::npos;
  }
  simulation.stop();
  return rejected &&
                 component->attributeTyped<Matrix>("x")->get().isApprox(initial)
             ? 0
             : 1;
}
