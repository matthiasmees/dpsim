// SPDX-License-Identifier: MPL-2.0
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fmi3Functions.h>

namespace {
struct Instance {
  std::array<double, 3> inputs{};
  std::array<double, 3> outputs{};
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
fmi3Status fmi3SetFloat64(fmi3Instance instance,
                          const fmi3ValueReference references[], size_t count,
                          const fmi3Float64 values[], size_t valueCount) {
  if (!instance || count != valueCount)
    return fmi3Error;
  auto *model = static_cast<Instance *>(instance);
  for (size_t k = 0; k < count; ++k) {
    if (references[k] > 2)
      return fmi3Error;
    model->inputs[references[k]] = values[k];
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
    if (references[k] < 3 || references[k] > 5)
      return fmi3Error;
    values[k] = model->outputs[references[k] - 3];
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
  if (std::abs(currentTime - model->time) >
      1e-10 * std::max(1.0, std::abs(model->time)))
    return fmi3Error;
  for (size_t phase = 0; phase < 3; ++phase)
    model->outputs[phase] = model->inputs[phase] / 100.0;
  model->time += stepSize;
  *eventHandlingNeeded = fmi3False;
  *terminateSimulation = fmi3False;
  *earlyReturn = fmi3False;
  *lastSuccessfulTime = model->time;
  return fmi3OK;
}
}
