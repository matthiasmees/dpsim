/*
 * Minimal FMI 2.0 Model Exchange function declarations used by DPsim.
 * The declarations follow the FMI 2.0 specification and are released under
 * the same 2-Clause BSD license as the official FMI headers.
 */
#pragma once

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  fmi2OK,
  fmi2Warning,
  fmi2Discard,
  fmi2Error,
  fmi2Fatal,
  fmi2Pending
} fmi2Status;
typedef void *fmi2Component;
typedef void *fmi2ComponentEnvironment;
typedef void *fmi2FMUstate;
typedef unsigned int fmi2ValueReference;
typedef double fmi2Real;
typedef int fmi2Integer;
typedef int fmi2Boolean;
typedef const char *fmi2String;
typedef unsigned char fmi2Byte;

#define fmi2True 1
#define fmi2False 0

typedef void (*fmi2CallbackLogger)(fmi2ComponentEnvironment, fmi2String,
                                   fmi2Status, fmi2String, fmi2String, ...);
typedef void *(*fmi2CallbackAllocateMemory)(size_t, size_t);
typedef void (*fmi2CallbackFreeMemory)(void *);
typedef void (*fmi2StepFinished)(fmi2ComponentEnvironment, fmi2Status);

typedef struct {
  fmi2CallbackLogger logger;
  fmi2CallbackAllocateMemory allocateMemory;
  fmi2CallbackFreeMemory freeMemory;
  fmi2StepFinished stepFinished;
  fmi2ComponentEnvironment componentEnvironment;
} fmi2CallbackFunctions;

typedef struct {
  fmi2Boolean newDiscreteStatesNeeded;
  fmi2Boolean terminateSimulation;
  fmi2Boolean nominalsOfContinuousStatesChanged;
  fmi2Boolean valuesOfContinuousStatesChanged;
  fmi2Boolean nextEventTimeDefined;
  fmi2Real nextEventTime;
} fmi2EventInfo;

typedef const char *fmi2GetVersionTYPE(void);
typedef fmi2Component fmi2InstantiateTYPE(fmi2String, int, fmi2String,
                                          fmi2String,
                                          const fmi2CallbackFunctions *,
                                          fmi2Boolean, fmi2Boolean);
typedef void fmi2FreeInstanceTYPE(fmi2Component);
typedef fmi2Status fmi2SetupExperimentTYPE(fmi2Component, fmi2Boolean, fmi2Real,
                                           fmi2Real, fmi2Boolean, fmi2Real);
typedef fmi2Status fmi2EnterInitializationModeTYPE(fmi2Component);
typedef fmi2Status fmi2ExitInitializationModeTYPE(fmi2Component);
typedef fmi2Status fmi2TerminateTYPE(fmi2Component);
typedef fmi2Status fmi2ResetTYPE(fmi2Component);
typedef fmi2Status fmi2SetTimeTYPE(fmi2Component, fmi2Real);
typedef fmi2Status fmi2SetRealTYPE(fmi2Component, const fmi2ValueReference[],
                                   size_t, const fmi2Real[]);
typedef fmi2Status fmi2GetRealTYPE(fmi2Component, const fmi2ValueReference[],
                                   size_t, fmi2Real[]);
typedef fmi2Status fmi2GetContinuousStatesTYPE(fmi2Component, fmi2Real[],
                                               size_t);
typedef fmi2Status fmi2SetContinuousStatesTYPE(fmi2Component, const fmi2Real[],
                                               size_t);
typedef fmi2Status fmi2GetDerivativesTYPE(fmi2Component, fmi2Real[], size_t);
typedef fmi2Status fmi2GetEventIndicatorsTYPE(fmi2Component, fmi2Real[],
                                              size_t);
typedef fmi2Status fmi2GetNominalsOfContinuousStatesTYPE(fmi2Component,
                                                         fmi2Real[], size_t);
typedef fmi2Status fmi2GetNumberOfContinuousStatesTYPE(fmi2Component, size_t *);
typedef fmi2Status fmi2GetNumberOfEventIndicatorsTYPE(fmi2Component, size_t *);
typedef fmi2Status fmi2CompletedIntegratorStepTYPE(fmi2Component, fmi2Boolean,
                                                   fmi2Boolean *,
                                                   fmi2Boolean *);
typedef fmi2Status fmi2EnterEventModeTYPE(fmi2Component);
typedef fmi2Status fmi2NewDiscreteStatesTYPE(fmi2Component, fmi2EventInfo *);
typedef fmi2Status fmi2EnterContinuousTimeModeTYPE(fmi2Component);
typedef fmi2Status fmi2GetFMUstateTYPE(fmi2Component, fmi2FMUstate *);
typedef fmi2Status fmi2SetFMUstateTYPE(fmi2Component, fmi2FMUstate);
typedef fmi2Status fmi2FreeFMUstateTYPE(fmi2Component, fmi2FMUstate *);

#ifdef __cplusplus
}
#endif
