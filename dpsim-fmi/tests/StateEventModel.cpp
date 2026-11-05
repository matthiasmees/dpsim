// SPDX-License-Identifier: MPL-2.0
// FMI 3 fixture for an endpoint-detected state event with a periodic reset.
#include <cmath>
#include <cstring>
#include <fmi3Functions.h>
#include <new>

namespace {
constexpr double twoPi = 6.283185307179586476925286766559;
struct State {
  double x[6]{6.28, 0, 0, 0, 0, 0};
  double u[3]{0, 0, 0};
  double time = 0;
  int mode = 0;
  int eventCount = 0;
};
State *state(fmi3Instance instance) { return static_cast<State *>(instance); }
void derivatives(double *values) {
  values[0] = 10;
  for (int k = 1; k < 6; ++k)
    values[k] = 0;
}
double current(const State *s, int phase) {
  return 0.01 * s->u[phase] * (1 + 0.2 * std::sin(s->x[0]));
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
                                    fmi3Boolean *nominalsChanged,
                                    fmi3Boolean *statesChanged,
                                    fmi3Boolean *nextEventDefined,
                                    fmi3Float64 *nextEventTime) {
  auto *s = state(instance);
  if (s->mode != 2)
    return fmi3Error;
  *again = *terminate = *nominalsChanged = *nextEventDefined = fmi3False;
  *statesChanged = fmi3False;
  *nextEventTime = 0;
  if (s->x[0] >= twoPi) {
    s->x[0] -= twoPi;
    ++s->eventCount;
    *statesChanged = fmi3True;
  }
  return fmi3OK;
}
fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance instance) {
  if (state(instance)->mode != 2)
    return fmi3Error;
  state(instance)->mode = 3;
  return fmi3OK;
}
fmi3Status fmi3EnterEventMode(fmi3Instance instance) {
  if (state(instance)->mode != 3)
    return fmi3Error;
  state(instance)->mode = 2;
  return fmi3OK;
}
fmi3Status fmi3GetNumberOfContinuousStates(fmi3Instance, size_t *count) {
  *count = 6;
  return fmi3OK;
}
fmi3Status fmi3GetNumberOfEventIndicators(fmi3Instance, size_t *count) {
  *count = 1;
  return fmi3OK;
}
fmi3Status fmi3SetTime(fmi3Instance instance, fmi3Float64 time) {
  if (state(instance)->mode != 3)
    return fmi3Error;
  state(instance)->time = time;
  return fmi3OK;
}
fmi3Status fmi3GetContinuousStates(fmi3Instance instance, fmi3Float64 values[],
                                   size_t count) {
  if (count != 6)
    return fmi3Error;
  std::memcpy(values, state(instance)->x, sizeof(state(instance)->x));
  return fmi3OK;
}
fmi3Status fmi3SetContinuousStates(fmi3Instance instance,
                                   const fmi3Float64 values[], size_t count) {
  if (count != 6 || state(instance)->mode != 3)
    return fmi3Error;
  std::memcpy(state(instance)->x, values, sizeof(state(instance)->x));
  return fmi3OK;
}
fmi3Status fmi3GetContinuousStateDerivatives(fmi3Instance instance,
                                             fmi3Float64 values[],
                                             size_t count) {
  if (count != 6 || state(instance)->mode != 3)
    return fmi3Error;
  derivatives(values);
  return fmi3OK;
}
fmi3Status fmi3GetEventIndicators(fmi3Instance instance, fmi3Float64 values[],
                                  size_t count) {
  if (count != 1 || state(instance)->mode != 3)
    return fmi3Error;
  values[0] = twoPi - state(instance)->x[0];
  return fmi3OK;
}
fmi3Status fmi3SetFloat64(fmi3Instance instance,
                          const fmi3ValueReference refs[], size_t nRefs,
                          const fmi3Float64 values[], size_t nValues) {
  if (nRefs != nValues)
    return fmi3Error;
  for (size_t k = 0; k < nRefs; ++k) {
    if (refs[k] >= 3)
      return fmi3Error;
    state(instance)->u[refs[k]] = values[k];
  }
  return fmi3OK;
}
fmi3Status fmi3GetFloat64(fmi3Instance instance,
                          const fmi3ValueReference refs[], size_t nRefs,
                          fmi3Float64 values[], size_t nValues) {
  if (nRefs != nValues || state(instance)->mode != 3)
    return fmi3Error;
  double dx[6];
  derivatives(dx);
  const auto *s = state(instance);
  for (size_t k = 0; k < nRefs; ++k) {
    const auto vr = refs[k];
    if (vr < 3)
      values[k] = s->u[vr];
    else if (vr < 6)
      values[k] = current(s, static_cast<int>(vr - 3));
    else if (vr == 9)
      values[k] = s->time;
    else if (vr >= 10 && vr < 16)
      values[k] = s->x[vr - 10];
    else if (vr >= 16 && vr < 22)
      values[k] = dx[vr - 16];
    else if (vr == 30)
      values[k] = twoPi - s->x[0];
    else if (vr == 31)
      values[k] = s->eventCount;
    else
      return fmi3Error;
  }
  return fmi3OK;
}
fmi3Status fmi3CompletedIntegratorStep(fmi3Instance instance, fmi3Boolean,
                                       fmi3Boolean *event,
                                       fmi3Boolean *terminate) {
  if (state(instance)->mode != 3)
    return fmi3Error;
  *event = *terminate = fmi3False;
  return fmi3OK;
}
fmi3Status fmi3GetNominalsOfContinuousStates(fmi3Instance, fmi3Float64 values[],
                                             size_t count) {
  if (count != 6)
    return fmi3Error;
  for (size_t k = 0; k < count; ++k)
    values[k] = 1;
  return fmi3OK;
}
