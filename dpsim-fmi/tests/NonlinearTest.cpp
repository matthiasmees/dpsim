// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <algorithm>
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
Matrix derivative(const Matrix &x, const Matrix &u) {
  Matrix f(6, 1);
  for (int p = 0; p < 3; ++p) {
    f(p, 0) = -(100 + 50 * p) * x(p, 0) - (10 + p) * x(3 + p, 0) +
              (20 + p) * u(p, 0) + 2 * u((p + 1) % 3, 0) +
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

std::pair<Matrix, Matrix> analyticalNorton(const Matrix &x, const Matrix &u,
                                           const Matrix &correction, double h) {
  Matrix b = Matrix::Zero(6, 3), c = Matrix::Zero(3, 6), d = Matrix::Zero(3, 3);
  for (int p = 0; p < 3; ++p) {
    b(p, p) = 20 + p;
    b(p, (p + 1) % 3) = 2;
    b(3 + p, (p + 2) % 3) = .4;
    c(p, p) = -1; // NonlinearSSN uses currentSignToPassive = -1 below.
    d(p, p) = -.02;
    d(p, (p + 1) % 3) = -.003;
  }
  const Matrix sensitivity = (Matrix::Identity(6, 6) - h * jacobian(x) / 2)
                                 .fullPivLu()
                                 .solve(h * b / 2);
  const Matrix w = d + c * sensitivity;
  Matrix current(3, 1);
  for (int p = 0; p < 3; ++p)
    current(p, 0) = -(x(p, 0) + .02 * u(p, 0) + .003 * u((p + 1) % 3, 0));
  return {w, current - correction.topRows(3) - w * u};
}

void requireClose(const Matrix &actual, const Matrix &expected,
                  double tolerance, const char *message) {
  if (actual.rows() != expected.rows() || actual.cols() != expected.cols() ||
      (actual - expected).norm() > tolerance * (1 + expected.norm()))
    throw std::runtime_error(message);
}
int main(int argc, char **argv) {
  if (argc < 2 || argc > 3)
    return 2;
  try {
    const bool reuse = argc == 3 && std::string(argv[2]) == "reuse";
    const bool tight = argc == 3 && std::string(argv[2]) == "tight";
    const bool failure = argc == 3 && !reuse && !tight;
    std::vector<FMI::ScalarInput> inputs{{"u_c", 13}, {"u_a", 7}, {"u_b", -11}};
    auto model = std::make_shared<FMI::ModelExchange>(argv[1], inputs);
    auto autoBranch = FMI::importFMU(argv[1]);
    require(autoBranch->modelDescription().inputNames.size() == 3,
            "automatic abc input mapping failed");
    FMI::NewtonSettings settings;
    settings.maxJacobianReuseSteps = reuse ? 16 : 0;
    if (tight) {
      settings.absoluteStateTolerance *= 0.1;
      settings.relativeStateTolerance *= 0.1;
      settings.absoluteCurrentTolerance *= 0.1;
      settings.relativeCurrentTolerance *= 0.1;
      settings.relativeJacobianPerturbation *= 0.5;
    }
    if (failure)
      settings.maxIterations = 1;
    auto branch = FMI::NonlinearSSN::make(
        "nonlinear", model,
        FMI::ElectricalPort3Ph{
            {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, -1},
        settings);
    auto second = FMI::NonlinearSSN::make(
        "second", std::make_shared<FMI::ModelExchange>(argv[1], inputs),
        FMI::ElectricalPort3Ph{
            {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, -1},
        settings);
    auto voltage = [](double t) -> Matrix {
      Matrix u(3, 1);
      u << 17 + 80 * std::sin(314 * t), -23 + 60 * std::sin(517 * t),
          9 - 40 * std::sin(231 * t);
      return u;
    };
    const Matrix initialNetworkVoltage = voltage(0);
    std::vector<Complex> initialNodeVoltage(3);
    for (int phase = 0; phase < 3; ++phase)
      initialNodeVoltage[phase] =
          initialNetworkVoltage(phase, 0) / RMS3PH_TO_PEAK1PH;
    auto node = SimNode<Real>::make("bus", PhaseType::ABC, initialNodeVoltage);
    auto apiNode =
        SimNode<Real>::make("api_bus", PhaseType::ABC, initialNodeVoltage);
    autoBranch->connect(apiNode);
    auto source = VoltageBoundary3Ph::make("source", voltage);
    for (auto comp : {std::static_pointer_cast<SimPowerComp<Real>>(branch),
                      std::static_pointer_cast<SimPowerComp<Real>>(second),
                      std::static_pointer_cast<SimPowerComp<Real>>(source)})
      comp->connect({SimNode<Real>::GND, node});
    SystemTopology system(50, {node}, {source, branch, second});
    Simulation simulation("nonlinear_test", Logger::Level::off);
    const double h = failure ? .01 : 2e-4;
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(h);
    simulation.setFinalTime(failure ? .01 : .02);
    simulation.doSteadyStateInit(false);
    auto logger = DataLogger::make(
        String(failure ? "nonlinear_rejection" : "nonlinear_test"));
    logger->logAttribute("source_current", source->attribute("i_intf"));
    logger->logAttribute("output", branch->attribute("fmi.output.i_a"));
    simulation.addLogger(logger);
    simulation.start();
    Matrix actualInitialState =
        branch->attributeTyped<Matrix>("norton_linearization_state")->get();
    Matrix stateCorrection =
        branch->attributeTyped<Matrix>("norton_state_correction")->get();
    require(actualInitialState.rows() == 6 && stateCorrection.rows() == 6,
            "missing nonlinear Norton linearization diagnostics");
    Matrix actualInitialVoltage =
        branch->attributeTyped<Matrix>("norton_linearization_voltage")->get();
    const auto analytic = analyticalNorton(
        actualInitialState, actualInitialVoltage, stateCorrection, h);
    requireClose(branch->attributeTyped<Matrix>("norton_W")->get(),
                 analytic.first, 2e-5,
                 "importer Norton W differs from analytical tangent");
    requireClose(branch->attributeTyped<Matrix>("norton_history")->get(),
                 analytic.second, 2e-5,
                 "importer Norton history differs from analytical intercept");
    require(std::abs(branch->attributeTyped<Real>("fmi.input.u_a")->get() -
                     initialNetworkVoltage(0, 0)) < 1e-12,
            "FMU did not adopt the connected initial terminal voltage");
    require(std::abs(branch->attributeTyped<Real>("fmi.input.u_b")->get() -
                     initialNetworkVoltage(1, 0)) < 1e-12,
            "FMU did not adopt the connected initial terminal voltage");
    require(std::abs(branch->attributeTyped<Real>("fmi.input.u_c")->get() -
                     initialNetworkVoltage(2, 0)) < 1e-12,
            "FMU did not adopt the connected initial terminal voltage");
    Matrix expected = model->description().initialState;
    if (failure) {
      bool rejected = false;
      try {
        simulation.step();
      } catch (const std::exception &e) {
        rejected =
            std::string(e.what()).find("convergence") != std::string::npos;
      }
      require(rejected, "unconverged nonlinear step was accepted");
      require(branch->attributeTyped<Matrix>("x")->get().isApprox(expected),
              "failed step committed states");
      simulation.stop();
      return 0;
    }
    Matrix previous = voltage(0);
    bool iterated = false;
    double maximumNewtonIterations = 0, maximumScaledResidual = 0;
    double maxNortonWError = 0;
    double maxNortonHistoryError = 0;
    for (int k = 1; k <= 100; ++k) {
      const Matrix u = voltage(k * h), old = expected,
                   f0 = derivative(old, previous);
      expected = old + h * f0;
      for (int it = 0; it < 15; ++it) {
        const Matrix r =
            expected - old - h * (f0 + derivative(expected, u)) / 2;
        if (r.norm() < 1e-13)
          break;
        const Matrix left = Matrix::Identity(6, 6) - h * jacobian(expected) / 2;
        expected -= left.fullPivLu().solve(r);
      }
      simulation.step();
      const Matrix stampState =
          branch->attributeTyped<Matrix>("norton_linearization_state")->get();
      const Matrix stampVoltage =
          branch->attributeTyped<Matrix>("norton_linearization_voltage")->get();
      const Matrix stampCorrection =
          branch->attributeTyped<Matrix>("norton_state_correction")->get();
      const auto stampExpected =
          analyticalNorton(stampState, stampVoltage, stampCorrection, h);
      const Matrix stampedW = branch->attributeTyped<Matrix>("norton_W")->get();
      const Matrix stampedHistory =
          branch->attributeTyped<Matrix>("norton_history")->get();
      maxNortonWError =
          std::max(maxNortonWError, (stampedW - stampExpected.first).norm() /
                                        (1 + stampExpected.first.norm()));
      maxNortonHistoryError =
          std::max(maxNortonHistoryError,
                   (stampedHistory - stampExpected.second).norm() /
                       (1 + stampExpected.second.norm()));
      require(maxNortonWError < 2e-5,
              "Norton W differs from analytic value at operating point");
      require(maxNortonHistoryError < 2e-5,
              "Norton history differs from analytic value at operating point");
      const auto &names = model->description().outputNames;
      if (std::find(names.begin(), names.end(), "evaluation_count") !=
          names.end())
        require(
            branch->attributeTyped<Real>("fmi.output.evaluation_count")
                    ->get() == k,
            "provisional FMI evaluations changed committed internal history");
      require(
          (branch->attributeTyped<Matrix>("x")->get() - expected).norm() < 2e-9,
          "FMI/MNA differs from independent nonlinear trapezoidal recurrence");
      require((second->attributeTyped<Matrix>("x")->get() - expected).norm() <
                  2e-9,
              "independent FMU instances interfered");
      Matrix current(3, 1);
      for (int p = 0; p < 3; ++p) {
        current(p, 0) =
            expected(p, 0) + .02 * u(p, 0) + .003 * u((p + 1) % 3, 0);
        const std::string phase(1, "abc"[p]);
        require(
            std::abs(
                branch->attributeTyped<Real>("fmi.output.i_" + phase)->get() -
                current(p, 0)) < 2e-9,
            "nonlinear output mapping failed");
      }
      require((source->attributeTyped<Matrix>("i_intf")->get() - 2 * current)
                      .norm() < 2e-8,
              "nonlinear injection violates network KCL");
      require(branch->attributeTyped<Real>("newton_residual")->get() <= 1,
              "nonlinear residual not converged");
      iterated |= branch->attributeTyped<Real>("newton_iterations")->get() > 1;
      maximumNewtonIterations =
          std::max(maximumNewtonIterations,
                   branch->attributeTyped<Real>("newton_iterations")->get());
      maximumScaledResidual =
          std::max(maximumScaledResidual,
                   branch->attributeTyped<Real>("newton_residual")->get());
      previous = u;
    }
    require(iterated, "fixture did not exercise repeated Newton iterations");
    require(branch->attributeTyped<Real>("fmi_evaluation_requests")->get() > 0,
            "FMI evaluation request counter did not advance");
    require(branch->attributeTyped<Real>("fmi_evaluation_wall_time_s")->get() >
                0,
            "FMI evaluation wall-time counter did not advance");
    simulation.stop();
    std::cout
        << "PASS: nonlinear FMI; maximum relative Norton W error="
        << maxNortonWError << ", history error=" << maxNortonHistoryError
        << "; FMI evaluation requests="
        << branch->attributeTyped<Real>("fmi_evaluation_requests")->get()
        << ", evaluation wall time="
        << branch->attributeTyped<Real>("fmi_evaluation_wall_time_s")->get()
        << " s; max iterations=" << maximumNewtonIterations
        << ", max scaled residual=" << maximumScaledResidual
        << "; state/Newton residual, shuffled ports, "
           "polarity, two instances and MNA KCL\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
