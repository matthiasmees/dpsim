// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fmi3Functions.h>

namespace {
struct Instance {
  double input = 0;
  double output = 0;
  double time = 0;
};
} // namespace

extern "C" {
fmi3Instance fmi3InstantiateCoSimulation(
    fmi3String, fmi3String, fmi3String, fmi3Boolean, fmi3Boolean, fmi3Boolean,
    fmi3Boolean, const fmi3ValueReference[], size_t, fmi3InstanceEnvironment,
    fmi3LogMessageCallback, fmi3IntermediateUpdateCallback) {
  return new Instance();
}
const char *fmi3GetVersion(void) { return "3.0"; }
void fmi3FreeInstance(fmi3Instance instance) {
  delete static_cast<Instance *>(instance);
}
fmi3Status fmi3EnterInitializationMode(fmi3Instance, fmi3Boolean, fmi3Float64,
                                       fmi3Float64, fmi3Boolean, fmi3Float64) {
  return fmi3OK;
}
fmi3Status fmi3ExitInitializationMode(fmi3Instance) { return fmi3OK; }
fmi3Status fmi3Terminate(fmi3Instance) { return fmi3OK; }
fmi3Status fmi3Reset(fmi3Instance instance) {
  if (!instance)
    return fmi3Error;
  *static_cast<Instance *>(instance) = Instance{};
  return fmi3OK;
}
fmi3Status fmi3GetFMUState(fmi3Instance instance, fmi3FMUState *state) {
  if (!instance || !state)
    return fmi3Error;
  *state = new Instance(*static_cast<Instance *>(instance));
  return fmi3OK;
}
fmi3Status fmi3SetFMUState(fmi3Instance instance, fmi3FMUState state) {
  if (!instance || !state)
    return fmi3Error;
  *static_cast<Instance *>(instance) = *static_cast<Instance *>(state);
  return fmi3OK;
}
fmi3Status fmi3FreeFMUState(fmi3Instance, fmi3FMUState *state) {
  if (!state || !*state)
    return fmi3Error;
  delete static_cast<Instance *>(*state);
  *state = nullptr;
  return fmi3OK;
}
fmi3Status fmi3SetFloat64(fmi3Instance instance,
                          const fmi3ValueReference references[], size_t count,
                          const fmi3Float64 values[], size_t valueCount) {
  if (!instance || count != valueCount)
    return fmi3Error;
  auto *model = static_cast<Instance *>(instance);
  for (size_t k = 0; k < count; ++k) {
    if (references[k] != 0)
      return fmi3Error;
    model->input = values[k];
  }
  return fmi3OK;
}
fmi3Status fmi3GetFloat64(fmi3Instance instance,
                          const fmi3ValueReference references[], size_t count,
                          fmi3Float64 values[], size_t valueCount) {
  if (!instance || count != valueCount)
    return fmi3Error;
  auto *model = static_cast<Instance *>(instance);
  for (size_t k = 0; k < count; ++k) {
    if (references[k] != 1)
      return fmi3Error;
    values[k] = model->output;
  }
  return fmi3OK;
}
fmi3Status fmi3DoStep(fmi3Instance instance, fmi3Float64 currentTime,
                      fmi3Float64 stepSize, fmi3Boolean,
                      fmi3Boolean *eventHandlingNeeded,
                      fmi3Boolean *terminateSimulation,
                      fmi3Boolean *earlyReturn,
                      fmi3Float64 *lastSuccessfulTime) {
  if (!instance || !eventHandlingNeeded || !terminateSimulation ||
      !earlyReturn || !lastSuccessfulTime)
    return fmi3Error;
  auto *model = static_cast<Instance *>(instance);
  if (std::abs(model->time - currentTime) > 1e-12)
    return fmi3Error;
  model->output += stepSize * model->input;
  model->time += stepSize;
  *eventHandlingNeeded = fmi3False;
  *terminateSimulation = fmi3False;
  *earlyReturn = fmi3False;
  *lastSuccessfulTime = model->time;
  return fmi3OK;
}
}
