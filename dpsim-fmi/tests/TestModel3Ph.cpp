// SPDX-License-Identifier: MPL-2.0
// Minimal FMI ABI fixture for importer tests, not a certified general-purpose FMU.
#define FIXTURE_KIND 0
#ifndef FIXTURE_NONLINEAR
#define FIXTURE_NONLINEAR 0
#endif
#include <cstring>
#include <fmi3Functions.h>
#include <new>

namespace {
struct State {
  double x[6]{0.25, -0.5, 1, 2, -3, 4};
  double u[3]{0, 0, 0};
  double t = 0;
  int mode = 0;
  bool currentJump = false;
  double sampledGain = 1;
  double command = 1;
  int tick = 1;
  int evaluationCount = 0;
};
State *state(fmi3Instance i) { return static_cast<State *>(i); }
void derivative(State *s, double *d) {
  for (int p = 0; p < 3; ++p) {
    d[p] = -(100 + 50 * p) * s->x[p] - (10 + p) * s->x[3 + p] +
           (20 + p) * s->u[p] + 2 * s->u[(p + 1) % 3] + .3 * s->x[(p + 2) % 3];
#if FIXTURE_NONLINEAR
    d[p] -=
        .1 * s->x[p] * s->x[p] * s->x[p] + .02 * s->x[p] * s->x[(p + 1) % 3];
#endif
#ifdef FIXTURE_SAMPLED
    d[p] += (s->sampledGain - 1) * s->u[p];
#endif
    d[3 + p] = (1000 + 100 * p) * s->x[p] + .4 * s->u[(p + 2) % 3];
  }
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
  for (int p = 0; p < 3; ++p)
    state(instance)->x[p] += .01 * state(instance)->u[p];
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
#ifdef FIXTURE_EVENT_JUMP
  if (state(instance)->t > 0)
    state(instance)->currentJump = true;
#endif
  *again = *terminate = *nominals = *values = *timeEvent = fmi3False;
  *nextTime = 0;
#ifdef FIXTURE_SAMPLED
  if (state(instance)->t >= state(instance)->tick * .001 - 1e-12) {
    state(instance)->sampledGain =
        1 + .1 * state(instance)->tick * state(instance)->command;
    ++state(instance)->tick;
  }
  *timeEvent = fmi3True;
  *nextTime = state(instance)->tick * .001;
#endif
  return fmi3OK;
}
fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance instance) {
  if (state(instance)->mode != 2)
    return fmi3Error;
  state(instance)->mode = 3;
  return fmi3OK;
}
fmi3Status fmi3GetNumberOfContinuousStates(fmi3Instance, size_t *n) {
  *n = 6;
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
  if (n != 6)
    return fmi3Error;
  std::memcpy(x, state(instance)->x, sizeof(state(instance)->x));
  return fmi3OK;
}
fmi3Status fmi3SetContinuousStates(fmi3Instance instance, const fmi3Float64 x[],
                                   size_t n) {
  if (n != 6 || state(instance)->mode != 3)
    return fmi3Error;
  std::memcpy(state(instance)->x, x, sizeof(state(instance)->x));
  return fmi3OK;
}
#if FIXTURE_KIND != 7
fmi3Status fmi3GetContinuousStateDerivatives(fmi3Instance instance,
                                             fmi3Float64 d[], size_t n) {
  if (n != 6 || state(instance)->mode != 3)
    return fmi3Error;
  derivative(state(instance), d);
  return fmi3OK;
}
#endif
fmi3Status fmi3SetFloat64(fmi3Instance instance,
                          const fmi3ValueReference refs[], size_t nRefs,
                          const fmi3Float64 values[], size_t nValues) {
  if (nRefs != nValues)
    return fmi3Error;
  for (size_t k = 0; k < nRefs; ++k) {
#ifdef FIXTURE_SAMPLED
    if (refs[k] == 23) {
      state(instance)->command = values[k];
      continue;
    }
#endif
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
  double d[6];
  const auto s = state(instance);
#ifdef FIXTURE_SNAPSHOTS
  if (s->t > 0)
    ++s->evaluationCount;
#endif
  derivative(s, d);
  for (size_t k = 0; k < nRefs; ++k) {
    const auto vr = refs[k];
    if (vr < 3)
      values[k] = s->u[vr];
    else if (vr < 6) {
      const auto p = vr - 3;
      values[k] = s->x[p] + .02 * s->u[p] + .003 * s->u[(p + 1) % 3] +
                  (s->currentJump ? 1 : 0);
    } else if (vr < 9)
      values[k] = s->x[vr - 3];
    else if (vr == 9)
      values[k] = s->t;
    else if (vr < 16)
      values[k] = s->x[vr - 10];
    else if (vr < 22)
      values[k] = d[vr - 16];
    else if (vr == 22)
      values[k] = s->sampledGain;
    else if (vr == 24)
      values[k] = s->evaluationCount;
    else
      return fmi3Error;
  }
  return fmi3OK;
}

fmi3Status fmi3CompletedIntegratorStep(fmi3Instance instance, fmi3Boolean,
                                       fmi3Boolean *event,
                                       fmi3Boolean *terminate) {
#ifdef FIXTURE_FORBID_COMPLETED
  return fmi3Error;
#endif
  if (state(instance)->mode != 3)
    return fmi3Error;
  *event = *terminate = fmi3False;
#ifdef FIXTURE_EVENT_JUMP
  *event = !state(instance)->currentJump;
#endif
  return fmi3OK;
}
fmi3Status fmi3GetNominalsOfContinuousStates(fmi3Instance, fmi3Float64 values[],
                                             size_t n) {
  if (n != 6)
    return fmi3Error;
  for (size_t k = 0; k < n; ++k)
    values[k] = 1;
  return fmi3OK;
}
fmi3Status fmi3EnterEventMode(fmi3Instance instance) {
  if (state(instance)->mode != 3)
    return fmi3Error;
  state(instance)->mode = 2;
  return fmi3OK;
}

#ifdef FIXTURE_SNAPSHOTS
fmi3Status fmi3GetFMUState(fmi3Instance instance, fmi3FMUState *snapshot) {
  if (!*snapshot)
    *snapshot = new (std::nothrow) State;
  if (!*snapshot)
    return fmi3Error;
  *static_cast<State *>(*snapshot) = *state(instance);
  return fmi3OK;
}
fmi3Status fmi3SetFMUState(fmi3Instance instance, fmi3FMUState snapshot) {
  if (!snapshot)
    return fmi3Error;
  *state(instance) = *static_cast<State *>(snapshot);
  return fmi3OK;
}
fmi3Status fmi3FreeFMUState(fmi3Instance, fmi3FMUState *snapshot) {
  delete static_cast<State *>(*snapshot);
  *snapshot = nullptr;
  return fmi3OK;
}
#endif
