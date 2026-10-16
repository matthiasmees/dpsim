// SPDX-License-Identifier: MPL-2.0
// Standalone lifecycle diagnostic for closed, zero-continuous-state ME FMUs.
// This does not demonstrate electrical coupling or DPsim state integration.
#include <cmath>
#include <dpsim-fmi/ModelExchange.h>
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: FMI_DiscreteProbe model.fmu stop_s step_s\n";
    return 2;
  }
  try {
    using namespace DPsim::FMI;
    ModelExchange model(argv[1], {});
    const auto &d = model.description();
    if (d.initialState.rows() != 0)
      throw std::runtime_error(
          "This lifecycle probe requires zero continuous states");
    const double stop = std::stod(argv[2]), h = std::stod(argv[3]);
    if (!std::isfinite(stop) || !std::isfinite(h) || h <= 0 || stop <= 0 ||
        std::abs(stop / h - std::round(stop / h)) > 1e-7)
      throw std::runtime_error("Invalid time grid");
    std::ofstream out("discrete_probe.csv");
    out << "time";
    for (const auto &name : d.outputNames)
      out << ",\"" << name << "\"";
    out << '\n' << std::setprecision(17);
    for (long long k = 0; k <= std::llround(stop / h); ++k) {
      const double t = k * h;
      if (k) {
        model.prepareStep(d.initialState, d.initialInputs, t);
        model.completeStep();
      }
      const auto values = model.acceptedValues().second;
      out << t;
      for (Eigen::Index j = 0; j < values.rows(); ++j)
        out << ',' << values(j, 0);
      out << '\n';
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
