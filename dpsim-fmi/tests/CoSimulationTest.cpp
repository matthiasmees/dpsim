// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <dpsim-fmi/CoSimulation.h>
#include <iostream>
#include <stdexcept>

using namespace CPS;
using namespace DPsim::FMI;

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
    CoSimulation model(argv[1]);
    require(model.description().inputNames.size() == 1,
            "expected one Co-Simulation input");
    require(model.description().outputNames.size() == 1,
            "expected one Co-Simulation output");
    model.setInput("u", 2.0);
    const auto &first = model.doStep(0.0, 0.1);
    require(std::abs(first(0, 0) - 0.2) < 1e-12,
            "first communication step output mismatch");
    const auto &second = model.doStep(0.1, 0.1);
    require(std::abs(second(0, 0) - 0.4) < 1e-12,
            "second communication step output mismatch");
    require(std::abs(model.time() - 0.2) < 1e-12,
            "FMU communication time mismatch");
    bool rejected = false;
    try {
      model.doStep(0.15, 0.1);
    } catch (const std::exception &) {
      rejected = true;
    }
    require(rejected, "nonsequential communication point was accepted");
    std::cout << "PASS: FMI Co-Simulation initialization, input mapping, "
                 "doStep sequence, and output retrieval\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
