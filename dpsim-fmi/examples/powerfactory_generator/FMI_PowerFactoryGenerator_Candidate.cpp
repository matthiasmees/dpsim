// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/ModelExchange.h>
#include <iostream>

using namespace DPsim::FMI;

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: FMI_PowerFactoryGenerator_Candidate model.fmu\n";
    return 2;
  }

  try {
    ModelExchange model(argv[1]);
    const auto &description = model.description();
    std::cout << "Imported PowerFactory generator candidate: "
              << description.modelName
              << "\nmodel identifier: " << description.modelIdentifier
              << "\ncontinuous states: " << description.stateNames.size()
              << "\ninputs:";
    for (const auto &name : description.inputNames)
      std::cout << ' ' << name;
    std::cout << "\noutputs:";
    for (const auto &name : description.outputNames)
      std::cout << ' ' << name;
    std::cout << "\nrollback: "
              << (model.supportsRollback() ? "available" : "unavailable")
              << "\nstatus: imported; requires an explicit P/Q or dq adapter "
                 "before native abc EMT coupling\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
