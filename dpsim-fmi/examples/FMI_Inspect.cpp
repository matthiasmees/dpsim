// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <dpsim-fmi/ModelExchange.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
int main(int argc, char **argv) {
  if (argc < 2 || argc > 4) {
    std::cerr << "Usage: FMI_Inspect model.fmu [input_values.json] "
                 "[--endpoint-events]\n";
    return 2;
  }
  try {
    std::vector<DPsim::FMI::ScalarInput> inputs;
    DPsim::FMI::ModelExchangeOptions options;
    for (int argument = 2; argument < argc; ++argument) {
      if (std::string(argv[argument]) == "--endpoint-events") {
        options.endpointStateEvents = true;
        options.allowContinuousStateResets = true;
        continue;
      }
      nlohmann::json values;
      std::ifstream file(argv[argument]);
      if (!file)
        throw std::runtime_error("Cannot open input values JSON");
      file >> values;
      for (auto it = values.begin(); it != values.end(); ++it)
        inputs.push_back({it.key(), it.value().get<double>()});
    }
    DPsim::FMI::ModelExchange model(argv[1], inputs, options);
    const auto &d = model.description();
    nlohmann::json result = {{"model", d.modelName},
                             {"states", d.stateNames},
                             {"inputs", d.inputNames},
                             {"outputs", d.outputNames},
                             {"discrete_variables", d.hasDiscreteVariables}};
    result["next_time_event_s"] = std::isfinite(model.nextTimeEvent())
                                      ? nlohmann::json(model.nextTimeEvent())
                                      : nlohmann::json(nullptr);
    std::cout << result.dump(2) << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
