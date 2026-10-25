// SPDX-License-Identifier: MPL-2.0
/*
 * FMI 3.0 Model Exchange example: three-phase RL grid load controlled by
 * one-step finite-control-set MPC. This is an educational FMU fixture, not a
 * protection-grade or production inverter controller.
 */
#include <fmi3FunctionTypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum {
  VR_VA = 0,
  VR_VB,
  VR_VC,
  VR_IA,
  VR_IB,
  VR_IC,
  VR_REFA,
  VR_REFB,
  VR_REFC,
  VR_UA,
  VR_UB,
  VR_UC,
  VR_COST,
  VR_TICKS,
  VR_IA_STATE,
  VR_IB_STATE,
  VR_IC_STATE,
  VR_DIA,
  VR_DIB,
  VR_DIC
};
static const double Ts = 50e-6;
static const double Lf = 25e-3;
static const double Rf = 0.4;
static const double VdcHalf = 350.0;
static const double IpkBeforeStep = 4.0;
static const double IpkAfterStep = -4.0;
static const double setpointStepTime = 5.0;
static const double omega = 2.0 * 3.14159265358979323846 * 50.0;
static const double switchPenalty = 0.05;

typedef struct {
  double i[3], v[3], uInv[3], t, nextTick, cost;
  unsigned long ticks;
  int mode;
} MpcState;

static MpcState *getState(fmi3Instance instance) {
  return (MpcState *)instance;
}
static double currentRef(double time, int phase) {
  /* At 5 s, reverse active-power flow: consume before the step, export after. */
  const double amplitude =
      time < setpointStepTime ? IpkBeforeStep : IpkAfterStep;
  return amplitude *
         cos(omega * time - phase * 2.0 * 3.14159265358979323846 / 3.0);
}
static void chooseSwitchState(MpcState *s) {
  double best = INFINITY;
  double bestU[3] = {0, 0, 0};
  const double targetTime = s->t + Ts;
  for (int state = 0; state < 8; ++state) {
    double cost = 0.0;
    double candidate[3];
    for (int p = 0; p < 3; ++p) {
      candidate[p] = (state & (1 << p)) ? VdcHalf : -VdcHalf;
      const double predicted =
          s->i[p] + Ts / Lf * (s->v[p] - candidate[p] - Rf * s->i[p]);
      const double error = predicted - currentRef(targetTime, p);
      cost += error * error;
      if ((candidate[p] > 0) != (s->uInv[p] > 0))
        cost += switchPenalty;
    }
    if (cost < best) {
      best = cost;
      memcpy(bestU, candidate, sizeof(bestU));
    }
  }
  memcpy(s->uInv, bestU, sizeof(bestU));
  s->cost = best;
  ++s->ticks;
}

const char *fmi3GetVersion(void) { return "3.0"; }
fmi3Instance fmi3InstantiateModelExchange(
    fmi3String instanceName, fmi3String token, fmi3String resourcePath,
    fmi3Boolean visible, fmi3Boolean loggingOn,
    fmi3InstanceEnvironment environment, fmi3LogMessageCallback logger) {
  (void)instanceName;
  (void)token;
  (void)resourcePath;
  (void)visible;
  (void)loggingOn;
  (void)environment;
  (void)logger;
  MpcState *s = (MpcState *)calloc(1, sizeof(MpcState));
  if (s) {
    s->nextTick = 0.0;
    s->mode = 0;
  }
  return (fmi3Instance)s;
}
void fmi3FreeInstance(fmi3Instance instance) { free(instance); }
fmi3Status fmi3Terminate(fmi3Instance instance) {
  if (!instance)
    return fmi3Error;
  getState(instance)->mode = 4;
  return fmi3OK;
}
fmi3Status
fmi3EnterInitializationMode(fmi3Instance instance, fmi3Boolean toleranceDefined,
                            fmi3Float64 tolerance, fmi3Float64 startTime,
                            fmi3Boolean stopTimeDefined, fmi3Float64 stopTime) {
  (void)toleranceDefined;
  (void)tolerance;
  (void)stopTimeDefined;
  (void)stopTime;
  if (!instance)
    return fmi3Error;
  MpcState *s = getState(instance);
  s->t = startTime;
  s->mode = 1;
  return fmi3OK;
}
fmi3Status fmi3ExitInitializationMode(fmi3Instance instance) {
  if (!instance || getState(instance)->mode != 1)
    return fmi3Error;
  getState(instance)->mode = 2; /* FMI event mode */
  return fmi3OK;
}
fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance instance) {
  if (!instance || getState(instance)->mode != 2)
    return fmi3Error;
  getState(instance)->mode = 3;
  return fmi3OK;
}
fmi3Status fmi3EnterEventMode(fmi3Instance instance) {
  if (!instance || getState(instance)->mode != 3)
    return fmi3Error;
  getState(instance)->mode = 2;
  return fmi3OK;
}
fmi3Status fmi3UpdateDiscreteStates(
    fmi3Instance instance, fmi3Boolean *discreteStatesNeedUpdate,
    fmi3Boolean *terminateSimulation,
    fmi3Boolean *nominalsOfContinuousStatesChanged,
    fmi3Boolean *valuesOfContinuousStatesChanged,
    fmi3Boolean *nextEventTimeDefined, fmi3Float64 *nextEventTime) {
  if (!instance || getState(instance)->mode != 2 || !discreteStatesNeedUpdate ||
      !terminateSimulation || !nominalsOfContinuousStatesChanged ||
      !valuesOfContinuousStatesChanged || !nextEventTimeDefined ||
      !nextEventTime)
    return fmi3Error;
  MpcState *s = getState(instance);
  *discreteStatesNeedUpdate = fmi3False;
  *terminateSimulation = fmi3False;
  *nominalsOfContinuousStatesChanged = fmi3False;
  *valuesOfContinuousStatesChanged = fmi3False;
  if (s->t + 1e-12 >= s->nextTick) {
    chooseSwitchState(s);
    s->nextTick += Ts;
    *valuesOfContinuousStatesChanged = fmi3False;
  }
  *nextEventTimeDefined = fmi3True;
  *nextEventTime = s->nextTick;
  return fmi3OK;
}
fmi3Status fmi3SetTime(fmi3Instance instance, fmi3Float64 time) {
  if (!instance || getState(instance)->mode != 3 || !isfinite(time))
    return fmi3Error;
  getState(instance)->t = time;
  return fmi3OK;
}
fmi3Status fmi3GetNumberOfContinuousStates(fmi3Instance instance, size_t *n) {
  if (!instance || !n)
    return fmi3Error;
  *n = 3;
  return fmi3OK;
}
fmi3Status fmi3GetNumberOfEventIndicators(fmi3Instance instance, size_t *n) {
  if (!instance || !n)
    return fmi3Error;
  *n = 0;
  return fmi3OK;
}
fmi3Status fmi3GetContinuousStates(fmi3Instance instance, fmi3Float64 x[],
                                   size_t n) {
  if (!instance || !x || n != 3)
    return fmi3Error;
  memcpy(x, getState(instance)->i, 3 * sizeof(double));
  return fmi3OK;
}
fmi3Status fmi3SetContinuousStates(fmi3Instance instance, const fmi3Float64 x[],
                                   size_t n) {
  if (!instance || !x || n != 3 || getState(instance)->mode != 3)
    return fmi3Error;
  memcpy(getState(instance)->i, x, 3 * sizeof(double));
  return fmi3OK;
}
fmi3Status fmi3GetContinuousStateDerivatives(fmi3Instance instance,
                                             fmi3Float64 d[], size_t n) {
  if (!instance || !d || n != 3 || getState(instance)->mode != 3)
    return fmi3Error;
  MpcState *s = getState(instance);
  for (int p = 0; p < 3; ++p)
    d[p] = (s->v[p] - s->uInv[p] - Rf * s->i[p]) / Lf;
  return fmi3OK;
}
fmi3Status fmi3GetNominalsOfContinuousStates(fmi3Instance instance,
                                             fmi3Float64 x[], size_t n) {
  if (!instance || !x || n != 3)
    return fmi3Error;
  for (size_t k = 0; k < n; ++k)
    x[k] = 4.0;
  return fmi3OK;
}
fmi3Status fmi3SetFloat64(fmi3Instance instance, const fmi3ValueReference vr[],
                          size_t nvr, const fmi3Float64 values[],
                          size_t nvalues) {
  if (!instance || !vr || !values || nvr != nvalues)
    return fmi3Error;
  MpcState *s = getState(instance);
  for (size_t k = 0; k < nvr; ++k) {
    if (vr[k] < VR_VA || vr[k] > VR_VC || !isfinite(values[k]))
      return fmi3Error;
    s->v[vr[k] - VR_VA] = values[k];
  }
  return fmi3OK;
}
fmi3Status fmi3GetFloat64(fmi3Instance instance, const fmi3ValueReference vr[],
                          size_t nvr, fmi3Float64 values[], size_t nvalues) {
  if (!instance || !vr || !values || nvr != nvalues ||
      getState(instance)->mode != 3)
    return fmi3Error;
  MpcState *s = getState(instance);
  for (size_t k = 0; k < nvr; ++k) {
    if (vr[k] >= VR_IA && vr[k] <= VR_IC)
      values[k] = s->i[vr[k] - VR_IA];
    else if (vr[k] >= VR_REFA && vr[k] <= VR_REFC)
      values[k] = currentRef(s->t, (int)(vr[k] - VR_REFA));
    else if (vr[k] >= VR_UA && vr[k] <= VR_UC)
      values[k] = s->uInv[vr[k] - VR_UA];
    else if (vr[k] == VR_COST)
      values[k] = s->cost;
    else if (vr[k] == VR_TICKS)
      values[k] = (double)s->ticks;
    else
      return fmi3Error;
  }
  return fmi3OK;
}
fmi3Status fmi3GetFMUState(fmi3Instance instance, fmi3FMUState *fmustate) {
  if (!instance || !fmustate)
    return fmi3Error;
  MpcState *copy = (MpcState *)malloc(sizeof(MpcState));
  if (!copy)
    return fmi3Error;
  *copy = *getState(instance);
  *fmustate = (fmi3FMUState)copy;
  return fmi3OK;
}
fmi3Status fmi3SetFMUState(fmi3Instance instance, fmi3FMUState fmustate) {
  if (!instance || !fmustate)
    return fmi3Error;
  *getState(instance) = *(MpcState *)fmustate;
  return fmi3OK;
}
fmi3Status fmi3FreeFMUState(fmi3Instance instance, fmi3FMUState *fmustate) {
  (void)instance;
  if (!fmustate || !*fmustate)
    return fmi3Error;
  free(*fmustate);
  *fmustate = NULL;
  return fmi3OK;
}
