// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary.h"
#include <cmath>
#include <dpsim-fmi/LinearSSN.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Simulation.h>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace CPS;
using namespace DPsim;
void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}
int main(int argc, char **argv) {
  if (argc != 3)
    return 2;
  try {
    const std::string mode = argv[2];
    if (mode.rfind("reject:", 0) == 0) {
      try {
        FMI::LinearModel::load(argv[1], {"u", "i"});
      } catch (const std::exception &e) {
        require(std::string(e.what()).find(mode.substr(7)) != std::string::npos,
                "wrong malformed-FMU diagnostic");
        return 0;
      }
      throw std::runtime_error("malformed FMU was accepted");
    }
    const std::map<std::string, std::string> failures{
        {"nonlinear", "LTI validation"},
        {"affine", "affine offsets"},
        {"timevarying", "LTI validation"},
        {"event", "state events"},
        {"initerror", "enter initialization mode failed"},
        {"missingfunction", "missing function"}};
    if (failures.count(mode)) {
      try {
        FMI::LinearModel::load(argv[1], {"u", "i"});
      } catch (const std::exception &e) {
        require(std::string(e.what()).find(failures.at(mode)) !=
                    std::string::npos,
                "wrong rejection diagnostic");
        return 0;
      }
      throw std::runtime_error("unsupported model was accepted");
    }
    const auto model = FMI::LinearModel::load(argv[1], {"u", "i"});
    if (mode == "singleton") {
      try {
        FMI::LinearModel::load(argv[1], {"u", "i"});
      } catch (const std::exception &e) {
        require(std::string(e.what()).find("one instantiation") !=
                    std::string::npos,
                "wrong singleton diagnostic");
        return 0;
      }
      throw std::runtime_error("singleton FMU was instantiated twice");
    }
    Matrix a(2, 2), b(2, 1);
    a << -500, -50, 10000, 0;
    b << 50, 0;
    require(model->A.isApprox(a, 1e-12) && model->B.isApprox(b, 1e-12),
            "incorrect recovered state matrices");
    require(model->C.isApprox(Matrix::Identity(2, 2)) && model->D.isZero(),
            "incorrect recovered output matrices");
    require(model->initialState(0, 0) == 0.25 &&
                model->initialState(1, 0) == -0.5,
            "FMU initial states were lost");
    auto branch = FMI::LinearSSN::make("imported", model);
    auto clone =
        std::dynamic_pointer_cast<FMI::LinearSSN>(branch->clone("copy"));
    require(clone && clone->getStateCount() == 2, "SSN clone failed");
    const auto diagnostic = branch->attribute("fmi.output.v_C");
    require(diagnostic->getDependencies().count(branch->attribute("x")) == 1,
            "output attribute lacks scheduler state dependency");
    auto node = SimNode<Real>::make("terminal", PhaseType::Single);
    auto ground = SimNode<Real>::GND;
    branch->connect({ground, node});
    clone->connect({ground, node});
    auto source = VoltageBoundary::make(
        "boundary", [](Real t) { return 325.0 * std::sin(2 * PI * 50 * t); });
    source->connect({ground, node});
    SystemTopology system(50, {node}, {source, branch, clone});
    Simulation simulation("fmi_test");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    const double h = 1e-5;
    simulation.setTimeStep(h);
    simulation.setFinalTime(0.001);
    simulation.doSteadyStateInit(false);
    // A consumer retains the boundary post-step in the dependency scheduler.
    // Reading an attribute from C++ after step() alone does not retain its task.
    auto logger = DataLogger::make(String("fmi_test"));
    logger->logAttribute("source_current", source->attribute("i_intf"));
    logger->logAttribute("diagnostic", diagnostic);
    simulation.addLogger(logger);
    simulation.start();
    Matrix expected = model->initialState;
    require(branch->attributeTyped<Matrix>("x")->get().isApprox(expected),
            "network initialization replaced FMU initial states");
    const Matrix left = Matrix::Identity(2, 2) - h * a / 2;
    const Matrix right = Matrix::Identity(2, 2) + h * a / 2;
    double previousVoltage = 0;
    for (int k = 1; k <= 100; ++k) {
      const double voltage = 325 * std::sin(2 * PI * 50 * k * h);
      expected = left.fullPivLu().solve(
          right * expected + h * b * (voltage + previousVoltage) / 2);
      simulation.step();
      const auto state = branch->attributeTyped<Matrix>("x")->get();
      require(state.isApprox(expected, 1e-9),
              "MNA/SSN integration differs from trapezoidal reference");
      require(std::abs(branch->attributeTyped<Real>("fmi.output.i")->get() -
                       expected(0, 0)) < 1e-9,
              "current output mapping is incorrect");
      require(std::abs(branch->attributeTyped<Real>("fmi.output.v_C")->get() -
                       expected(1, 0)) < 1e-9,
              "diagnostic output mapping is incorrect");
      require(std::abs(source->attributeTyped<Matrix>("i_intf")->get()(0, 0) +
                       2 * expected(0, 0)) < 1e-9,
              "imported branch currents do not participate in network KCL");
      previousVoltage = voltage;
    }
    simulation.stop();
    std::cout << "PASS: FMI matrix recovery, initial states, attributes and "
                 "MNA integration\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}
