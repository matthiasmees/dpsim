// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <algorithm>
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim/Simulation.h>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace CPS;
using namespace DPsim;

namespace {
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

struct Reference {
  Matrix state;
  double time = 0;
  double period;
  int tick = 1;

  Reference(Matrix initialState, double samplePeriod)
      : state(std::move(initialState)), period(samplePeriod) {}

  int advance(double endTime, const Matrix &voltage) {
    int steps = 0;
    while (time < endTime - 1e-12) {
      const double eventTime = tick * period;
      const double endpoint = std::min(endTime, eventTime);
      const double h = endpoint - time;
      const double gain = tick == 1 ? 1.0 : 1 + .1 * (tick - 1);
      const Matrix old = state, f0 = derivative(old, voltage, gain);
      Matrix next = old + h * f0;
      for (int it = 0; it < 20; ++it) {
        const Matrix residual =
            next - old - h * (f0 + derivative(next, voltage, gain)) / 2;
        if (residual.norm() < 1e-13)
          break;
        next -= (Matrix::Identity(6, 6) - h * jacobian(next) / 2)
                    .fullPivLu()
                    .solve(residual);
      }
      state = next;
      time = endpoint;
      ++steps;
      if (std::abs(time - eventTime) < 1e-12)
        ++tick;
    }
    return steps;
  }
};
} // namespace

int main(int argc, char **argv) {
  if (argc != 3)
    return 2;
  try {
    Matrix voltage(3, 1);
    voltage << 7, -11, 13;
    auto slow = FMI::NonlinearSSN::importFMU("slow_fmu", argv[1]);
    auto fast = FMI::NonlinearSSN::importFMU("fast_fmu", argv[2]);
    auto source = VoltageBoundary3Ph::make(
        "source", [voltage](double) -> Matrix { return voltage; });
    auto bus = SimNode<Real>::make(
        "bus", PhaseType::ABC,
        std::vector<Complex>{voltage(0, 0), voltage(1, 0), voltage(2, 0)});
    slow->connect({SimNode<Real>::GND, bus});
    fast->connect({SimNode<Real>::GND, bus});
    source->connect({SimNode<Real>::GND, bus});
    SystemTopology system(50, {bus}, {source, slow, fast});
    Simulation simulation("fmi_multirate_distinct_clocks", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    constexpr double macroStep = .002;
    simulation.setTimeStep(macroStep);
    simulation.setFinalTime(5 * macroStep);
    simulation.doSteadyStateInit(false);
    simulation.start();

    Reference slowReference(slow->modelDescription().initialState, .001);
    Reference fastReference(fast->modelDescription().initialState, .0005);
    double maxSlowError = 0;
    double maxFastError = 0;
    for (int macro = 1; macro <= 5; ++macro) {
      const int slowSteps = slowReference.advance(macro * macroStep, voltage);
      const int fastSteps = fastReference.advance(macro * macroStep, voltage);
      simulation.step();
      const double slowError =
          (slow->attributeTyped<Matrix>("x")->get() - slowReference.state)
              .norm();
      const double fastError =
          (fast->attributeTyped<Matrix>("x")->get() - fastReference.state)
              .norm();
      maxSlowError = std::max(maxSlowError, slowError);
      maxFastError = std::max(maxFastError, fastError);
      require(slowError < 2e-9,
              "slow FMU diverged from its single-rate recurrence");
      require(fastError < 2e-9,
              "fast FMU diverged from its single-rate recurrence");
      require(slow->attributeTyped<Real>("fmi_internal_steps")->get() ==
                      slowSteps &&
                  fast->attributeTyped<Real>("fmi_internal_steps")->get() ==
                      fastSteps,
              "FMUs did not retain independent event-aligned schedules");
    }
    simulation.stop();
    std::cout << "max state error: 1 ms clock=" << maxSlowError
              << ", 0.5 ms clock=" << maxFastError << " (tolerance 2e-9)\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
