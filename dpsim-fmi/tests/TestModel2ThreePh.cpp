// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <cstring>
#include <fmi2FunctionTypes.h>

struct State {
  double x[6]{.25, -.5, 1, 2, -3, 4};
  double u[3]{7, -11, 13};
  double time = 0;
};
static State *state(fmi2Component c) { return static_cast<State *>(c); }
static void derivative(const State *s, double *d) {
  for (int p = 0; p < 3; ++p) {
    d[p] = -(100 + 50 * p) * s->x[p] - (10 + p) * s->x[3 + p] +
           (20 + p) * s->u[p] + 2 * s->u[(p + 1) % 3] + .3 * s->x[(p + 2) % 3] -
           .1 * s->x[p] * s->x[p] * s->x[p] - .02 * s->x[p] * s->x[(p + 1) % 3];
    d[3 + p] = (1000 + 100 * p) * s->x[p] + .4 * s->u[(p + 2) % 3];
  }
}
extern "C" {
const char *fmi2GetVersion() { return "2.0"; }
fmi2Component fmi2Instantiate(fmi2String, int, fmi2String, fmi2String,
                              const fmi2CallbackFunctions *, fmi2Boolean,
                              fmi2Boolean) {
  return new State;
}
void fmi2FreeInstance(fmi2Component c) { delete state(c); }
fmi2Status fmi2SetupExperiment(fmi2Component, fmi2Boolean, fmi2Real, fmi2Real,
                               fmi2Boolean, fmi2Real) {
  return fmi2OK;
}
fmi2Status fmi2EnterInitializationMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2ExitInitializationMode(fmi2Component c) {
  auto *s = state(c);
  for (int p = 0; p < 3; ++p)
    s->x[p] += .01 * s->u[p];
  return fmi2OK;
}
fmi2Status fmi2Terminate(fmi2Component) { return fmi2OK; }
fmi2Status fmi2EnterContinuousTimeMode(fmi2Component) { return fmi2OK; }
fmi2Status fmi2SetTime(fmi2Component c, fmi2Real t) {
  state(c)->time = t;
  return fmi2OK;
}
fmi2Status fmi2SetReal(fmi2Component c, const fmi2ValueReference refs[],
                       size_t n, const fmi2Real values[]) {
  auto *s = state(c);
  for (size_t k = 0; k < n; ++k)
    if (refs[k] < 3)
      s->u[refs[k]] = values[k];
  return fmi2OK;
}
fmi2Status fmi2GetReal(fmi2Component c, const fmi2ValueReference refs[],
                       size_t n, fmi2Real values[]) {
  auto *s = state(c);
  double d[6];
  derivative(s, d);
  for (size_t k = 0; k < n; ++k) {
    const auto ref = refs[k];
    if (ref < 3)
      values[k] = s->u[ref];
    else if (ref < 6)
      values[k] =
          s->x[ref - 3] + .02 * s->u[ref - 3] + .003 * s->u[(ref - 2) % 3];
    else if (ref < 9)
      values[k] = s->x[ref - 6];
    else
      values[k] = d[ref - 9];
  }
  return fmi2OK;
}
fmi2Status fmi2GetNumberOfContinuousStates(fmi2Component, size_t *n) {
  *n = 6;
  return fmi2OK;
}
fmi2Status fmi2GetNumberOfEventIndicators(fmi2Component, size_t *n) {
  *n = 0;
  return fmi2OK;
}
fmi2Status fmi2GetEventIndicators(fmi2Component, fmi2Real[], size_t) {
  return fmi2OK;
}
fmi2Status fmi2GetContinuousStates(fmi2Component c, fmi2Real x[], size_t n) {
  std::memcpy(x, state(c)->x, n * sizeof(double));
  return fmi2OK;
}
fmi2Status fmi2SetContinuousStates(fmi2Component c, const fmi2Real x[],
                                   size_t n) {
  std::memcpy(state(c)->x, x, n * sizeof(double));
  return fmi2OK;
}
fmi2Status fmi2GetDerivatives(fmi2Component c, fmi2Real d[], size_t) {
  derivative(state(c), d);
  return fmi2OK;
}
fmi2Status fmi2GetNominalsOfContinuousStates(fmi2Component, fmi2Real x[],
                                             size_t n) {
  for (size_t k = 0; k < n; ++k)
    x[k] = 1;
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
