// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/ModelExchange.h>
#include <iostream>
#include <string>

using namespace DPsim::FMI;

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: FMI_MicrogridNetwork_Candidate model.fmu\n";
    return 2;
  }

  try {
    ModelExchangeOptions options;
    options.endpointStateEvents = true;
    ModelExchange model(argv[1], options);
    const auto &description = model.description();

    std::cout << "Imported microgrid candidate: " << description.modelName
              << "\nmodel identifier: " << description.modelIdentifier
              << "\ncontinuous states: " << description.stateNames.size()
              << "\nFMI inputs: " << description.inputNames.size()
              << "\nFMI outputs: " << description.outputNames.size() << '\n';

    std::cout << "inputs:";
    for (const auto &name : description.inputNames)
      std::cout << ' ' << name;
    std::cout << "\n";

    if (description.outputNames.empty()) {
      std::cout
          << "status: imported, but no causality=output variables are exposed\n"
          << "status_detail: terminal voltages are local Modelica variables; "
             "a current-driven electrical adapter is required before this FMU "
             "can be connected to a native DPsim MNA network\n";
    } else {
      std::cout << "outputs:";
      for (const auto &name : description.outputNames)
        std::cout << ' ' << name;
      std::cout << "\nstatus: output interface available for further mapping\n";
    }

    std::cout << "rollback: "
              << (model.supportsRollback() ? "available" : "unavailable")
              << "\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
