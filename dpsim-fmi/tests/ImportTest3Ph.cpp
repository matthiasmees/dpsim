// SPDX-License-Identifier: MPL-2.0
#include "VoltageBoundary3Ph.h"
#include <dpsim-fmi/LinearSSN.h>
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
int main(int argc, char **argv) {
  if (argc != 2 && argc != 3)
    return 2;
  try {
    const bool floating = argc == 3;
    const FMI::VoltagePorts ports{
        {"u_a", "u_b", "u_c"}, {"i_a", "i_b", "i_c"}, {7, -11, 13}};
    if (argc == 3 && std::string(argv[2]) == "reject_sampled") {
      try {
        FMI::LinearModel::loadMulti(argv[1], {{"u_a", "u_b", "u_c", "command"},
                                              {"i_a", "i_b", "i_c", "v_C_a"},
                                              {7, -11, 13, 1}});
      } catch (const std::exception &e) {
        require(std::string(e.what()).find("sampled/discrete") !=
                    std::string::npos,
                "wrong sampled LTI rejection");
        return 0;
      }
      throw std::runtime_error("sampled model was accepted as LTI");
    }
    // These known equations belong to the test fixture only. Production import
    // must recover them despite shuffled XML inputs and output declarations.
    const auto model = FMI::LinearModel::loadMulti(argv[1], ports);
    Matrix a = Matrix::Zero(6, 6), b = Matrix::Zero(6, 3);
    Matrix c = Matrix::Zero(3, 6), d = Matrix::Zero(3, 3);
    for (int p = 0; p < 3; ++p) {
      a(p, p) = -(100 + 50 * p);
      a(p, 3 + p) = -(10 + p);
      a(p, (p + 2) % 3) = .3;
      a(3 + p, p) = 1000 + 100 * p;
      b(p, p) = 20 + p;
      b(p, (p + 1) % 3) = 2;
      b(3 + p, (p + 2) % 3) = .4;
      c(p, p) = 1;
      d(p, p) = .02;
      d(p, (p + 1) % 3) = .003;
      require(model->C.row(model->currentOutputs[p]).isApprox(c.row(p)),
              "current output order is incorrect");
      require(model->D.row(model->currentOutputs[p]).isApprox(d.row(p)),
              "coupled feedthrough is incorrect");
    }
    require(model->A.isApprox(a) && model->B.isApprox(b),
            "coupled matrix recovery is incorrect");
    Matrix expected(6, 1);
    expected << .32, -.61, 1.13, 2, -3, 4;
    require(model->initialState.isApprox(expected),
            "input-dependent FMI initialization is incorrect");
    const auto waveform = [](Real t) -> Matrix {
      Matrix u(3, 1);
      u << 17 + 200 * std::sin(314 * t), -23 + 83 * std::sin(517 * t),
          9 - 137 * std::sin(231 * t);
      return u;
    };
    const auto offset = [](Real t) -> Matrix {
      Matrix u(3, 1);
      u << 17 + 5 * std::sin(100 * t), -23 + 11 * std::sin(70 * t),
          9 - 3 * std::sin(80 * t);
      return u;
    };
    const Matrix nodeVoltage =
        waveform(0) + (floating ? offset(0) : Matrix::Zero(3, 1));
    std::vector<Complex> initialNodeVoltage(3);
    std::vector<Complex> initialNegativeVoltage(3);
    for (int phase = 0; phase < 3; ++phase) {
      initialNodeVoltage[phase] = nodeVoltage(phase, 0);
      initialNegativeVoltage[phase] =
          floating ? offset(0)(phase, 0) : Complex(0, 0);
    }
    auto node =
        SimNode<Real>::make("terminal", PhaseType::ABC, initialNodeVoltage);
    auto negative = floating ? SimNode<Real>::make("negative", PhaseType::ABC,
                                                   initialNegativeVoltage)
                             : SimNode<Real>::GND;
    auto branch = FMI::LinearSSN::make("imported", model);
    branch->connect({negative, node});
    auto source = VoltageBoundary3Ph::make("boundary", [=](Real t) -> Matrix {
      Matrix u = waveform(t);
      if (floating)
        u += offset(t);
      return u;
    });
    source->connect({SimNode<Real>::GND, node});
    auto clone =
        std::dynamic_pointer_cast<FMI::LinearSSN>(branch->clone("copy"));
    require(clone && clone->getStateCount() == 6, "ABC clone failed");
    clone->connect({negative, node});
    SystemTopology system(50, {node}, {source, branch, clone});
    VoltageBoundary3Ph::Ptr negativeSource;
    if (floating) {
      negativeSource = VoltageBoundary3Ph::make("negative_boundary", offset);
      negativeSource->connect({SimNode<Real>::GND, negative});
      system.addNode(negative);
      system.addComponent(negativeSource);
    }
    auto logger = DataLogger::make(String("fmi_test_3ph"));
    logger->logAttribute("source_current", source->attribute("i_intf"));
    for (const auto &output : model->outputNames)
      logger->logAttribute(output, branch->attribute("fmi.output." + output));
    if (floating)
      logger->logAttribute("negative_current",
                           negativeSource->attribute("i_intf"));
    Simulation simulation("fmi_test_3ph");
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    const Real h = 1e-5;
    simulation.setTimeStep(h);
    simulation.setFinalTime(.001);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);
    simulation.start();
    require(branch->attributeTyped<Matrix>("x")->get().isApprox(expected),
            "network initialization replaced imported states");
    require(
        branch->attributeTyped<Matrix>("v_intf")->get().isApprox(waveform(0)),
        "initial SSN voltage history is incorrect");
    const Matrix left = Matrix::Identity(6, 6) - h * a / 2;
    const Matrix right = Matrix::Identity(6, 6) + h * a / 2;
    Matrix previous = waveform(0);
    for (int k = 1; k <= 100; ++k) {
      const Matrix u = waveform(k * h);
      expected =
          left.fullPivLu().solve(right * expected + h * b * (u + previous) / 2);
      simulation.step();
      require(
          branch->attributeTyped<Matrix>("x")->get().isApprox(expected, 1e-10),
          "three-phase SSN state differs from independent reference");
      const Matrix current = c * expected + d * u;
      require((source->attributeTyped<Matrix>("i_intf")->get() + 2 * current)
                      .norm() < 1e-9,
              "coupled currents violate MNA KCL");
      if (floating)
        require((negativeSource->attributeTyped<Matrix>("i_intf")->get() -
                 2 * current)
                        .norm() < 1e-9,
                "negative terminal currents violate MNA KCL");
      for (int p = 0; p < 3; ++p) {
        const auto suffix = std::string(1, "abc"[p]);
        require(
            std::abs(
                branch->attributeTyped<Real>("fmi.output.i_" + suffix)->get() -
                current(p, 0)) < 1e-10,
            "current attribute mapping is incorrect");
        require(
            std::abs(branch->attributeTyped<Real>("fmi.output.v_C_" + suffix)
                         ->get() -
                     expected(3 + p, 0)) < 1e-10,
            "diagnostic attribute mapping is incorrect");
      }
      previous = u;
    }
    simulation.stop();
    auto invalid = FMI::LinearSSN::make("invalid", model);
    auto scalarNode = SimNode<Real>::make("scalar", PhaseType::Single);
    invalid->connect({SimNode<Real>::GND, scalarNode});
    bool phaseRejected = false;
    try {
      invalid->initializeFromNodesAndTerminals(50);
    } catch (const std::exception &) {
      phaseRejected = true;
    }
    require(phaseRejected, "three-phase branch accepted a scalar terminal");
    auto bad = ports;
    bad.voltageInputs[1] = bad.voltageInputs[0];
    bool rejected = false;
    try {
      FMI::LinearModel::loadMulti(argv[1], bad);
    } catch (const std::exception &) {
      rejected = true;
    }
    require(rejected, "duplicate voltage mapping accepted");
    std::cout << "PASS: coupled ABC recovery, shuffled mappings, "
                 "initialization, attributes and MNA KCL\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
