#define S_FUNCTION_NAME ieee39_network_sfun
#define S_FUNCTION_LEVEL 2
#include "simstruc.h"
#include <cstring>

extern "C" {
void *ieee39_network_create(double, double, double, double *, double *, char *,
                            int);
int ieee39_network_step(void *, const double *, double *, double *, char *,
                        int);
void ieee39_network_destroy(void *);
}

namespace {
constexpr int kBusVoltageCount = 39 * 3;
constexpr int kOutputCount = kBusVoltageCount + 2;
} // namespace

static void mdlInitializeSizes(SimStruct *S) {
  ssSetNumSFcnParams(S, 3);
  if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
    return;
  ssSetSFcnParamTunable(S, 0, 0);
  ssSetSFcnParamTunable(S, 1, 0);
  ssSetSFcnParamTunable(S, 2, 0);
  if (!ssSetNumInputPorts(S, 1))
    return;
  ssSetInputPortWidth(S, 0, 3);
  ssSetInputPortDataType(S, 0, SS_DOUBLE);
  ssSetInputPortDirectFeedThrough(S, 0, 0);
  ssSetInputPortRequiredContiguous(S, 0, 1);
  if (!ssSetNumOutputPorts(S, 1))
    return;
  ssSetOutputPortWidth(S, 0, kOutputCount);
  ssSetOutputPortDataType(S, 0, SS_DOUBLE);
  ssSetNumSampleTimes(S, 1);
  ssSetNumPWork(S, 1);
  ssSetNumDWork(S, 1);
  ssSetDWorkWidth(S, 0, kOutputCount);
  ssSetDWorkDataType(S, 0, SS_DOUBLE);
  ssSetNumContStates(S, 0);
  ssSetNumDiscStates(S, 0);
  ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE |
                      SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME);
}

static void mdlInitializeSampleTimes(SimStruct *S) {
  ssSetSampleTime(S, 0, mxGetScalar(ssGetSFcnParam(S, 0)));
  ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct *S) {
  const double dt = mxGetScalar(ssGetSFcnParam(S, 0));
  const double stop = mxGetScalar(ssGetSFcnParam(S, 1));
  const double loadStep = mxGetScalar(ssGetSFcnParam(S, 2));
  double initial[kBusVoltageCount] = {};
  double initSeconds = 0.0;
  char error[1024] = {};
  void *network = ieee39_network_create(dt, stop, loadStep, initial,
                                        &initSeconds, error, sizeof(error));
  if (!network) {
    ssSetErrorStatus(S, error[0] ? error
                                 : "Failed to initialize IEEE39 DPsim network");
    return;
  }
  ssSetPWorkValue(S, 0, network);
  real_T *state = static_cast<real_T *>(ssGetDWork(S, 0));
  std::memcpy(state, initial, sizeof(initial));
  state[kBusVoltageCount] = 0.0;
  state[kBusVoltageCount + 1] = initSeconds;
}

static void mdlOutputs(SimStruct *S, int_T) {
  const real_T *state = static_cast<const real_T *>(ssGetDWork(S, 0));
  real_T *output = ssGetOutputPortRealSignal(S, 0);
  std::memcpy(output, state, kOutputCount * sizeof(real_T));
}

#define MDL_UPDATE
static void mdlUpdate(SimStruct *S, int_T) {
  void *network = ssGetPWorkValue(S, 0);
  const real_T *current =
      static_cast<const real_T *>(ssGetInputPortSignal(S, 0));
  real_T *state = static_cast<real_T *>(ssGetDWork(S, 0));
  double nextVoltage[kBusVoltageCount] = {};
  double stepSeconds = 0.0;
  char error[1024] = {};
  if (ieee39_network_step(network, current, nextVoltage, &stepSeconds, error,
                          sizeof(error)) != 0) {
    ssSetErrorStatus(S, error[0] ? error : "DPsim IEEE39 network step failed");
    return;
  }
  std::memcpy(state, nextVoltage, sizeof(nextVoltage));
  state[kBusVoltageCount] = stepSeconds;
}

static void mdlTerminate(SimStruct *S) {
  void *network = ssGetPWorkValue(S, 0);
  if (network)
    ieee39_network_destroy(network);
}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
