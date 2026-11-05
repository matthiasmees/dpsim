// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <cstdlib>
#include <fmi2FunctionTypes.h>

namespace {
struct Instance {
  double input = 0;
  double output = 0;
  double time = 0;
};
} // namespace

extern "C" {
fmi2Component fmi2Instantiate(fmi2String, int type, fmi2String, fmi2String,
                              const fmi2CallbackFunctions *, fmi2Boolean,
                              fmi2Boolean) {
  return type == 1 ? static_cast<fmi2Component>(new Instance()) : nullptr;
}
const char *fmi2GetVersion(void) { return "2.0"; }
void fmi2FreeInstance(fmi2Component component) {
  delete static_cast<Instance *>(component);
}
fmi2Status fmi2SetupExperiment(fmi2Component, fmi2Boolean, fmi2Real, fmi2Real,
                               fmi2Boolean, fmi2Real) {
  return fmi2OK;
}
fmi2Status fmi2EnterInitializationMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2ExitInitializationMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2Terminate(fmi2Component) { return fmi2OK; }
fmi2Status fmi2SetReal(fmi2Component component,
                       const fmi2ValueReference references[], size_t count,
                       const fmi2Real values[]) {
  if (!component)
    return fmi2Error;
  auto *model = static_cast<Instance *>(component);
  for (size_t k = 0; k < count; ++k) {
    if (references[k] != 0)
      return fmi2Error;
    model->input = values[k];
  }
  return fmi2OK;
}
fmi2Status fmi2GetReal(fmi2Component component,
                       const fmi2ValueReference references[], size_t count,
                       fmi2Real values[]) {
  if (!component)
    return fmi2Error;
  auto *model = static_cast<Instance *>(component);
  for (size_t k = 0; k < count; ++k) {
    if (references[k] != 1)
      return fmi2Error;
    values[k] = model->output;
  }
  return fmi2OK;
}
fmi2Status fmi2DoStep(fmi2Component component, fmi2Real currentTime,
                      fmi2Real stepSize, fmi2Boolean) {
  if (!component)
    return fmi2Error;
  auto *model = static_cast<Instance *>(component);
  if (std::abs(model->time - currentTime) > 1e-12)
    return fmi2Error;
  model->output += stepSize * model->input;
  model->time += stepSize;
  return fmi2OK;
}
}
