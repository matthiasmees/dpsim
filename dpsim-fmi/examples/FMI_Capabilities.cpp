// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <dpsim-fmi/ModelExchange.h>
#include <iomanip>
#include <iostream>

using namespace DPsim::FMI;

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: FMI_Capabilities model.fmu\n";
    return 2;
  }
  try {
    ModelExchangeOptions options;
    options.endpointStateEvents = true;
    auto model = std::make_shared<ModelExchange>(argv[1], options);
    const auto &description = model->description();
    std::cout << "model: " << description.modelName << '\n'
              << "model_identifier: " << description.modelIdentifier << '\n'
              << "states: " << description.stateNames.size() << '\n'
              << "inputs: " << description.inputNames.size() << '\n'
              << "outputs: " << description.outputNames.size() << '\n'
              << "event_indicators: "
              << (description.hasEventIndicators ? "yes" : "no") << '\n'
              << "discrete_variables: "
              << (description.hasDiscreteVariables ? "yes" : "no") << '\n'
              << "needs_completed_integrator_step: "
              << (description.needsCompletedIntegratorStep ? "yes" : "no")
              << '\n'
              << "can_get_set_fmu_state: "
              << (description.canGetAndSetFMUState ? "yes" : "no") << '\n'
              << "can_serialize_fmu_state: "
              << (description.canSerializeFMUState ? "yes" : "no") << '\n'
              << "single_instance_process_limit: "
              << (description.canBeInstantiatedOnlyOncePerProcess ? "yes"
                                                                  : "no")
              << '\n'
              << "next_time_event: ";
    if (std::isfinite(model->nextTimeEvent()))
      std::cout << std::setprecision(17) << model->nextTimeEvent() << '\n';
    else
      std::cout << "none\n";

    std::cout << "input_names:";
    for (const auto &name : description.inputNames)
      std::cout << ' ' << name;
    std::cout << '\n' << "output_names:";
    for (const auto &name : description.outputNames)
      std::cout << ' ' << name;
    std::cout << '\n';

    std::string rollbackProbe = "not_available";
    if (model->supportsRollback()) {
      try {
        const auto replayMacroStep = [&]() {
          auto state = description.initialState;
          for (int step = 1; step <= 50; ++step) {
            double time = step * 1e-6;
            if (std::isfinite(model->nextTimeEvent()) &&
                model->nextTimeEvent() < time)
              time = model->nextTimeEvent();
            model->evaluate(state, description.initialInputs, time);
            model->completeStep();
            state = model->acceptedState();
          }
          return std::make_pair(state, model->acceptedValues());
        };
        model->saveCheckpoint();
        const auto first = replayMacroStep();
        model->restoreCheckpoint();
        const auto second = replayMacroStep();
        model->discardCheckpoint();
        rollbackProbe =
            (first.first.isApprox(second.first, 1e-10) &&
             first.second.first.isApprox(second.second.first, 1e-10) &&
             first.second.second.isApprox(second.second.second, 1e-10))
                ? "deterministic"
                : "nondeterministic";
      } catch (const std::exception &error) {
        rollbackProbe = std::string("probe_failed: ") + error.what();
      }
    }
    std::cout << "rollback_replay_probe: " << rollbackProbe << '\n'
              << "recommended_multirate_mode: "
              << (rollbackProbe == "deterministic" ? "rollback_newton"
                                                   : "single_pass_explicit")
              << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
