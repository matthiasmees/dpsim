// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>
#include <stdexcept>

using namespace CPS;
using namespace DPsim;
void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}
Matrix derivative(const Matrix &x, const Matrix &u, double gain) {
  Matrix f(6, 1);
  for (int p = 0; p < 3; ++p) {
    f(p, 0) = -(100 + 50 * p) * x(p, 0) - (10 + p) * x(3 + p, 0) +
              (20 + p + gain - 1) * u(p, 0) + 2 * u((p + 1) % 3, 0) +
              .3 * x((p + 2) % 3, 0) - .1 * std::pow(x(p, 0), 3) -
              .02 * x(p, 0) * x((p + 1) % 3, 0);
    f(3 + p, 0) = (1000 + 100 * p) * x(p, 0) + .4 * u((p + 2) % 3, 0);
  }
  return f;
}
Matrix jacobian(const Matrix &x) {
  Matrix a = Matrix::Zero(6, 6);
  for (int p = 0; p < 3; ++p) {
    a(p, p) =
        -(100 + 50 * p) - .3 * x(p, 0) * x(p, 0) - .02 * x((p + 1) % 3, 0);
    a(p, 3 + p) = -(10 + p);
    a(p, (p + 2) % 3) = .3;
    a(p, (p + 1) % 3) = -.02 * x(p, 0);
    a(3 + p, p) = 1000 + 100 * p;
  }
  return a;
}
int main(int argc, char **argv) {
  if (argc < 2 || argc > 3)
    return 2;
  try {
    const bool misaligned = argc == 3;
    Matrix u(3, 1);
    u << 7, -11, 13;
    auto runtime = std::make_shared<FMI::ModelExchange>(
        argv[1], std::vector<FMI::ScalarInput>{
                     {"u_a", 7}, {"u_b", -11}, {"u_c", 13}, {"command", 1}});
    require(runtime->description().hasDiscreteVariables,
            "discrete metadata was lost");
    auto component = FMI::NonlinearSSN::make(
        "sampled", runtime,
        FMI::ElectricalPort3Ph{
            {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, 1});
    if (!misaligned)
      component->setInputFunction(
          "command", [](double t) { return t >= .002 ? 2.0 : 1.0; }, {.002});
    auto source =
        VoltageBoundary3Ph::make("source", [=](double) -> Matrix { return u; });
    auto bus = SimNode<Real>::make("bus", PhaseType::ABC);
    component->connect({SimNode<Real>::GND, bus});
    source->connect({SimNode<Real>::GND, bus});
    SystemTopology system(50, {bus}, {source, component});
    Simulation simulation("sampled_test", Logger::Level::off);
    const double h = misaligned ? .0003 : .0002;
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(h);
    simulation.doSteadyStateInit(false);
    auto logger = DataLogger::make(String("sampled_test"));
    logger->logAttribute("source_current", source->attribute("i_intf"));
    logger->logAttribute("branch_current", component->attribute("i_intf"));
    simulation.addLogger(logger);
    simulation.start();
    Matrix expected = runtime->description().initialState;
    for (int k = 1; k <= 20; ++k) {
      const Matrix previous = component->attributeTyped<Matrix>("x")->get();
      if (misaligned) {
        try {
          simulation.step();
        } catch (const std::exception &e) {
          require(std::string(e.what()).find("align") != std::string::npos,
                  "wrong time-event rejection");
          require(
              component->attributeTyped<Matrix>("x")->get().isApprox(previous),
              "misaligned step committed state");
          simulation.stop();
          return 0;
        }
        continue;
      }
      const double gain = 1 + .1 * ((k - 1) / 5) * (k > 10 ? 2 : 1);
      const Matrix old = expected, f0 = derivative(old, u, gain);
      for (int it = 0; it < 15; ++it) {
        const Matrix r =
            expected - old - h * (f0 + derivative(expected, u, gain)) / 2;
        if (r.norm() < 1e-13)
          break;
        expected -= (Matrix::Identity(6, 6) - h * jacobian(expected) / 2)
                        .fullPivLu()
                        .solve(r);
      }
      simulation.step();
      require(
          (component->attributeTyped<Matrix>("x")->get() - expected).norm() <
              2e-9,
          "sampled controller differs from independent trapezoidal recurrence");
      require((component->attributeTyped<Matrix>("i_intf")->get() +
               source->attributeTyped<Matrix>("i_intf")->get())
                      .norm() < 1e-8,
              "sampled model violates MNA KCL");
      require(std::abs(runtime->nextTimeEvent() - (k / 5 + 1) * .001) < 1e-12,
              "wrong next controller tick");
    }
    simulation.stop();
    require(!misaligned, "misaligned event was skipped");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
