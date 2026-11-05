// SPDX-License-Identifier: MPL-2.0
// Minimal FMI ABI fixture for importer tests, not a certified general-purpose FMU.
#include <cstring>
#include <fmi3Functions.h>
#include <new>

namespace {
struct State {
  double x[2]{0.25, -0.5};
  double u = 0, t = 0;
  int mode = 0;
};
State *state(fmi3Instance i) { return static_cast<State *>(i); }
void derivative(State *s, double *d) {
  d[0] = -500 * s->x[0] - 50 * s->x[1] + 50 * s->u;
  d[1] = 10000 * s->x[0];
#if FIXTURE_KIND == 1
  d[0] += 0.1 * s->x[0] * s->x[0];
#elif FIXTURE_KIND == 2
  d[0] += 1;
#elif FIXTURE_KIND == 3
  d[0] += s->t * s->x[0];
#endif
}
} // namespace
const char *fmi3GetVersion() { return "3.0"; }
fmi3Instance fmi3InstantiateModelExchange(fmi3String, fmi3String, fmi3String,
                                          fmi3Boolean, fmi3Boolean,
                                          fmi3InstanceEnvironment,
                                          fmi3LogMessageCallback) {
  return new (std::nothrow) State;
}
void fmi3FreeInstance(fmi3Instance instance) { delete state(instance); }
fmi3Status fmi3Terminate(fmi3Instance instance) {
  state(instance)->mode = 4;
  return fmi3OK;
}
fmi3Status fmi3EnterInitializationMode(fmi3Instance instance, fmi3Boolean,
                                       fmi3Float64, fmi3Float64, fmi3Boolean,
                                       fmi3Float64) {
  if (FIXTURE_KIND == 5)
    return fmi3Error;
  state(instance)->mode = 1;
  return fmi3OK;
}
fmi3Status fmi3ExitInitializationMode(fmi3Instance instance) {
  if (state(instance)->mode != 1)
    return fmi3Error;
  state(instance)->mode = 2;
  return fmi3OK;
}
fmi3Status fmi3UpdateDiscreteStates(fmi3Instance instance, fmi3Boolean *again,
                                    fmi3Boolean *terminate,
                                    fmi3Boolean *nominals, fmi3Boolean *values,
                                    fmi3Boolean *timeEvent,
                                    fmi3Float64 *nextTime) {
  if (state(instance)->mode != 2)
    return fmi3Error;
  *again = *terminate = *nominals = *values = *timeEvent = fmi3False;
  *nextTime = 0;
  return fmi3OK;
}
fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance instance) {
  if (state(instance)->mode != 2)
    return fmi3Error;
  state(instance)->mode = 3;
  return fmi3OK;
}
fmi3Status fmi3GetNumberOfContinuousStates(fmi3Instance, size_t *n) {
  *n = 2;
  return fmi3OK;
}
fmi3Status fmi3GetNumberOfEventIndicators(fmi3Instance, size_t *n) {
  *n = FIXTURE_KIND == 4 ? 1 : 0;
  return fmi3OK;
}
fmi3Status fmi3SetTime(fmi3Instance instance, fmi3Float64 time) {
  if (state(instance)->mode != 3)
    return fmi3Error;
  state(instance)->t = time;
  return fmi3OK;
}
fmi3Status fmi3GetContinuousStates(fmi3Instance instance, fmi3Float64 x[],
                                   size_t n) {
  if (n != 2)
    return fmi3Error;
  std::memcpy(x, state(instance)->x, sizeof(state(instance)->x));
  return fmi3OK;
}
fmi3Status fmi3SetContinuousStates(fmi3Instance instance, const fmi3Float64 x[],
                                   size_t n) {
  if (n != 2 || state(instance)->mode != 3)
    return fmi3Error;
  std::memcpy(state(instance)->x, x, sizeof(state(instance)->x));
  return fmi3OK;
}
#if FIXTURE_KIND != 7
fmi3Status fmi3GetContinuousStateDerivatives(fmi3Instance instance,
                                             fmi3Float64 d[], size_t n) {
  if (n != 2 || state(instance)->mode != 3)
    return fmi3Error;
  derivative(state(instance), d);
  return fmi3OK;
}
#endif
fmi3Status fmi3SetFloat64(fmi3Instance instance,
                          const fmi3ValueReference refs[], size_t nRefs,
                          const fmi3Float64 values[], size_t nValues) {
  if (nRefs != 1 || nValues != 1 || refs[0] != 0)
    return fmi3Error;
  state(instance)->u = values[0];
  return fmi3OK;
}
fmi3Status fmi3GetFloat64(fmi3Instance instance,
                          const fmi3ValueReference refs[], size_t nRefs,
                          fmi3Float64 values[], size_t nValues) {
  if (nRefs != nValues || state(instance)->mode != 3)
    return fmi3Error;
  double d[2];
  derivative(state(instance), d);
  for (size_t k = 0; k < nRefs; ++k) {
    switch (refs[k]) {
    case 0:
      values[k] = state(instance)->u;
      break;
    case 1:
    case 4:
      values[k] = state(instance)->x[0];
      break;
    case 2:
    case 5:
      values[k] = state(instance)->x[1];
      break;
    case 3:
      values[k] = state(instance)->t;
      break;
    case 6:
      values[k] = d[0];
      break;
    case 7:
      values[k] = d[1];
      break;
    default:
      return fmi3Error;
    }
  }
  return fmi3OK;
}
