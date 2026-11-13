// SPDX-License-Identifier: MPL-2.0
// Minimal FMI 2.0 implementation of L*di/dt = u - R*i.
#include <fmi2FunctionTypes.h>

#include <algorithm>
#include <cmath>
#include <new>

namespace {
constexpr double kResistance = 0.1;
constexpr double kInductance = 20e-6;
struct Instance {
  double time = 0.0;
  double voltage = 0.0;
  double current = 0.0;
};
bool closeTime(double left, double right) {
  return std::abs(left - right) <= 1e-10 * std::max(1.0, std::abs(right));
}
} // namespace

extern "C" {
const char *fmi2GetTypesPlatform(void) { return "default"; }
const char *fmi2GetVersion(void) { return "2.0"; }
fmi2Status fmi2SetDebugLogging(fmi2Component, fmi2Boolean, size_t,
                               const fmi2String[]) {
  return fmi2OK;
}

fmi2Component fmi2Instantiate(fmi2String, int type, fmi2String, fmi2String,
                              const fmi2CallbackFunctions *, fmi2Boolean,
                              fmi2Boolean) {
#if defined(FMI_RL_COSIM)
  if (type != 1)
    return nullptr;
#else
  if (type != 0)
    return nullptr;
#endif
  return new (std::nothrow) Instance{};
}

void fmi2FreeInstance(fmi2Component component) {
  delete static_cast<Instance *>(component);
}
fmi2Status fmi2SetupExperiment(fmi2Component component, fmi2Boolean, fmi2Real,
                               fmi2Real startTime, fmi2Boolean, fmi2Real) {
  if (!component || !std::isfinite(startTime))
    return fmi2Error;
  static_cast<Instance *>(component)->time = startTime;
  return fmi2OK;
}
fmi2Status fmi2EnterInitializationMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2ExitInitializationMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2Terminate(fmi2Component) { return fmi2OK; }
fmi2Status fmi2Reset(fmi2Component component) {
  if (!component)
    return fmi2Error;
  *static_cast<Instance *>(component) = Instance{};
  return fmi2OK;
}
fmi2Status fmi2SetReal(fmi2Component component,
                       const fmi2ValueReference references[], size_t count,
                       const fmi2Real values[]) {
  if (!component || (count && (!references || !values)))
    return fmi2Error;
  auto &model = *static_cast<Instance *>(component);
  for (size_t k = 0; k < count; ++k) {
    if (references[k] != 0 || !std::isfinite(values[k]))
      return fmi2Error;
    model.voltage = values[k];
  }
  return fmi2OK;
}
fmi2Status fmi2GetInteger(fmi2Component, const fmi2ValueReference[], size_t n,
                          fmi2Integer[]) {
  return n == 0 ? fmi2OK : fmi2Error;
}
fmi2Status fmi2SetInteger(fmi2Component, const fmi2ValueReference[], size_t n,
                          const fmi2Integer[]) {
  return n == 0 ? fmi2OK : fmi2Error;
}
fmi2Status fmi2GetBoolean(fmi2Component, const fmi2ValueReference[], size_t n,
                          fmi2Boolean[]) {
  return n == 0 ? fmi2OK : fmi2Error;
}
fmi2Status fmi2SetBoolean(fmi2Component, const fmi2ValueReference[], size_t n,
                          const fmi2Boolean[]) {
  return n == 0 ? fmi2OK : fmi2Error;
}
fmi2Status fmi2GetString(fmi2Component, const fmi2ValueReference[], size_t n,
                         fmi2String[]) {
  return n == 0 ? fmi2OK : fmi2Error;
}
fmi2Status fmi2SetString(fmi2Component, const fmi2ValueReference[], size_t n,
                         const fmi2String[]) {
  return n == 0 ? fmi2OK : fmi2Error;
}
fmi2Status fmi2GetReal(fmi2Component component,
                       const fmi2ValueReference references[], size_t count,
                       fmi2Real values[]) {
  if (!component || (count && (!references || !values)))
    return fmi2Error;
  const auto &model = *static_cast<Instance *>(component);
  for (size_t k = 0; k < count; ++k) {
    if (references[k] != 1)
      return fmi2Error;
    values[k] = model.current;
  }
  return fmi2OK;
}

fmi2Status fmi2SetTime(fmi2Component component, fmi2Real time) {
  if (!component || !std::isfinite(time))
    return fmi2Error;
  static_cast<Instance *>(component)->time = time;
  return fmi2OK;
}
fmi2Status fmi2GetNumberOfContinuousStates(fmi2Component, size_t *count) {
  if (!count)
    return fmi2Error;
  *count = 1;
  return fmi2OK;
}
fmi2Status fmi2GetNumberOfEventIndicators(fmi2Component, size_t *count) {
  if (!count)
    return fmi2Error;
  *count = 0;
  return fmi2OK;
}
fmi2Status fmi2GetContinuousStates(fmi2Component component, fmi2Real values[],
                                   size_t count) {
  if (!component || count != 1 || !values)
    return fmi2Error;
  values[0] = static_cast<Instance *>(component)->current;
  return fmi2OK;
}
fmi2Status fmi2SetContinuousStates(fmi2Component component,
                                   const fmi2Real values[], size_t count) {
  if (!component || count != 1 || !values || !std::isfinite(values[0]))
    return fmi2Error;
  static_cast<Instance *>(component)->current = values[0];
  return fmi2OK;
}
fmi2Status fmi2GetDerivatives(fmi2Component component, fmi2Real values[],
                              size_t count) {
  if (!component || count != 1 || !values)
    return fmi2Error;
  const auto &model = *static_cast<Instance *>(component);
  values[0] = (model.voltage - kResistance * model.current) / kInductance;
  return fmi2OK;
}
fmi2Status fmi2GetEventIndicators(fmi2Component, fmi2Real[], size_t count) {
  return count == 0 ? fmi2OK : fmi2Error;
}
fmi2Status fmi2GetNominalsOfContinuousStates(fmi2Component, fmi2Real values[],
                                             size_t count) {
  if (count != 1 || !values)
    return fmi2Error;
  values[0] = 1.0;
  return fmi2OK;
}
fmi2Status fmi2CompletedIntegratorStep(fmi2Component, fmi2Boolean,
                                       fmi2Boolean *event,
                                       fmi2Boolean *terminate) {
  if (!event || !terminate)
    return fmi2Error;
  *event = fmi2False;
  *terminate = fmi2False;
  return fmi2OK;
}
fmi2Status fmi2EnterEventMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2NewDiscreteStates(fmi2Component, fmi2EventInfo *info) {
  if (!info)
    return fmi2Error;
  *info =
      fmi2EventInfo{fmi2False, fmi2False, fmi2False, fmi2False, fmi2False, 0.0};
  return fmi2OK;
}
fmi2Status fmi2EnterContinuousTimeMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2EnterStepMode(fmi2Component) { return fmi2OK; }

fmi2Status fmi2DoStep(fmi2Component component, fmi2Real currentTime,
                      fmi2Real stepSize, fmi2Boolean) {
  if (!component || !std::isfinite(currentTime) || !std::isfinite(stepSize) ||
      stepSize <= 0)
    return fmi2Error;
  auto &model = *static_cast<Instance *>(component);
  if (!closeTime(model.time, currentTime))
    return fmi2Error;
  const double decay = std::exp(-kResistance * stepSize / kInductance);
  model.current =
      decay * model.current + (1.0 - decay) * model.voltage / kResistance;
  model.time = currentTime + stepSize;
  return fmi2OK;
}
fmi2Status fmi2CancelStep(fmi2Component) { return fmi2Error; }
}
