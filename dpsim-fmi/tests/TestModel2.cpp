// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <cstdlib>
#include <fmi2FunctionTypes.h>

struct State {
  double time = 0;
  double input = 0;
  double x[2] = {0.25, -0.5};
};

extern "C" {
const char *fmi2GetVersion() { return "2.0"; }
fmi2Component fmi2Instantiate(fmi2String, int, fmi2String, fmi2String,
                              const fmi2CallbackFunctions *, fmi2Boolean,
                              fmi2Boolean) {
  return new State;
}
void fmi2FreeInstance(fmi2Component instance) {
  delete static_cast<State *>(instance);
}
fmi2Status fmi2SetupExperiment(fmi2Component, fmi2Boolean, fmi2Real, fmi2Real,
                               fmi2Boolean, fmi2Real) {
  return fmi2OK;
}
fmi2Status fmi2EnterInitializationMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2ExitInitializationMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2Terminate(fmi2Component) { return fmi2OK; }
fmi2Status fmi2EnterContinuousTimeMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2SetTime(fmi2Component instance, fmi2Real time) {
  static_cast<State *>(instance)->time = time;
  return fmi2OK;
}
fmi2Status fmi2SetReal(fmi2Component instance, const fmi2ValueReference refs[],
                       size_t count, const fmi2Real values[]) {
  auto *state = static_cast<State *>(instance);
  for (size_t k = 0; k < count; ++k)
    if (refs[k] == 0)
      state->input = values[k];
  return fmi2OK;
}
fmi2Status fmi2GetReal(fmi2Component instance, const fmi2ValueReference refs[],
                       size_t count, fmi2Real values[]) {
  auto *state = static_cast<State *>(instance);
  for (size_t k = 0; k < count; ++k) {
    if (refs[k] == 1)
      values[k] = state->x[0];
    else if (refs[k] == 2)
      values[k] = state->x[1];
    else
      values[k] = 0;
  }
  return fmi2OK;
}
fmi2Status fmi2GetNumberOfContinuousStates(fmi2Component, size_t *count) {
  *count = 2;
  return fmi2OK;
}
fmi2Status fmi2GetNumberOfEventIndicators(fmi2Component, size_t *count) {
  *count = 0;
  return fmi2OK;
}
fmi2Status fmi2GetEventIndicators(fmi2Component, fmi2Real[], size_t) {
  return fmi2OK;
}
fmi2Status fmi2GetContinuousStates(fmi2Component instance, fmi2Real values[],
                                   size_t count) {
  auto *state = static_cast<State *>(instance);
  for (size_t k = 0; k < count; ++k)
    values[k] = state->x[k];
  return fmi2OK;
}
fmi2Status fmi2SetContinuousStates(fmi2Component instance,
                                   const fmi2Real values[], size_t count) {
  auto *state = static_cast<State *>(instance);
  for (size_t k = 0; k < count; ++k)
    state->x[k] = values[k];
  return fmi2OK;
}
fmi2Status fmi2GetDerivatives(fmi2Component instance, fmi2Real values[],
                              size_t count) {
  auto *state = static_cast<State *>(instance);
  if (count > 0)
    values[0] = -500 * state->x[0] - 50 * state->x[1] + 50 * state->input;
  if (count > 1)
    values[1] = 10000 * state->x[0];
  return fmi2OK;
}
fmi2Status fmi2GetNominalsOfContinuousStates(fmi2Component, fmi2Real values[],
                                             size_t count) {
  for (size_t k = 0; k < count; ++k)
    values[k] = 1;
  return fmi2OK;
}
fmi2Status fmi2CompletedIntegratorStep(fmi2Component, fmi2Boolean,
                                       fmi2Boolean *event,
                                       fmi2Boolean *terminate) {
  *event = fmi2False;
  *terminate = fmi2False;
  return fmi2OK;
}
fmi2Status fmi2EnterEventMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2NewDiscreteStates(fmi2Component, fmi2EventInfo *info) {
  info->newDiscreteStatesNeeded = fmi2False;
  info->terminateSimulation = fmi2False;
  info->nominalsOfContinuousStatesChanged = fmi2False;
  info->valuesOfContinuousStatesChanged = fmi2False;
  info->nextEventTimeDefined = fmi2False;
  info->nextEventTime = 0;
  return fmi2OK;
}
}
