// SPDX-License-Identifier: MPL-2.0
// Small FMI 3 dual-interface test FMU for an independent three-phase R-L load.
#include <algorithm>
#include <array>
#include <cmath>
#include <fmi3Functions.h>
#include <new>

namespace {
struct Instance {
  std::array<double, 3> voltage{};
  std::array<double, 3> current{};
  double time = 0;
  bool initialized = false;
};

Instance *instance(fmi3Instance value) {
  return static_cast<Instance *>(value);
}

bool validReferences(const fmi3ValueReference references[], size_t count,
                     fmi3ValueReference maxValue) {
  for (size_t index = 0; index < count; ++index)
    if (references[index] > maxValue)
      return false;
  return true;
}

void derivative(const Instance *model, double values[3]) {
  constexpr double resistance = 8.0;
  constexpr double inductance = 0.02;
  for (size_t phase = 0; phase < 3; ++phase)
    values[phase] =
        (model->voltage[phase] - resistance * model->current[phase]) /
        inductance;
}
} // namespace

extern "C" {
const char *fmi3GetVersion(void) { return "3.0"; }

fmi3Instance fmi3InstantiateModelExchange(
    fmi3String, fmi3String, fmi3String, fmi3Boolean, fmi3Boolean,
    fmi3InstanceEnvironment, fmi3LogMessageCallback) {
  return new (std::nothrow) Instance;
}

fmi3Instance fmi3InstantiateCoSimulation(
    fmi3String, fmi3String, fmi3String, fmi3Boolean, fmi3Boolean, fmi3Boolean,
    fmi3Boolean, const fmi3ValueReference[], size_t,
    fmi3InstanceEnvironment, fmi3LogMessageCallback,
    fmi3IntermediateUpdateCallback) {
  return new (std::nothrow) Instance;
}

void fmi3FreeInstance(fmi3Instance value) { delete instance(value); }

fmi3Status fmi3EnterInitializationMode(fmi3Instance value, fmi3Boolean,
                                       fmi3Float64, fmi3Float64, fmi3Boolean,
                                       fmi3Float64) {
  if (!value)
    return fmi3Error;
  auto *model = instance(value);
  model->time = 0;
  model->initialized = false;
  return fmi3OK;
}

fmi3Status fmi3ExitInitializationMode(fmi3Instance value) {
  if (!value)
    return fmi3Error;
  instance(value)->initialized = true;
  return fmi3OK;
}

fmi3Status fmi3UpdateDiscreteStates(fmi3Instance value, fmi3Boolean *again,
                                    fmi3Boolean *terminate,
                                    fmi3Boolean *nominalsChanged,
                                    fmi3Boolean *valuesChanged,
                                    fmi3Boolean *timeEvent,
                                    fmi3Float64 *nextTime) {
  if (!value || !again || !terminate || !nominalsChanged || !valuesChanged ||
      !timeEvent || !nextTime)
    return fmi3Error;
  *again = *terminate = *nominalsChanged = *valuesChanged = *timeEvent =
      fmi3False;
  *nextTime = 0;
  return fmi3OK;
}

fmi3Status fmi3EnterContinuousTimeMode(fmi3Instance value) {
  return value ? fmi3OK : fmi3Error;
}

fmi3Status fmi3EnterEventMode(fmi3Instance value) {
  return value ? fmi3OK : fmi3Error;
}

fmi3Status fmi3Terminate(fmi3Instance value) {
  return value ? fmi3OK : fmi3Error;
}

fmi3Status fmi3SetFloat64(fmi3Instance value,
                          const fmi3ValueReference references[], size_t count,
                          const fmi3Float64 values[], size_t valueCount) {
  if (!value || count != valueCount || !validReferences(references, count, 2))
    return fmi3Error;
  auto *model = instance(value);
  for (size_t index = 0; index < count; ++index)
    model->voltage[references[index]] = values[index];
  return fmi3OK;
}

fmi3Status fmi3GetFloat64(fmi3Instance value,
                          const fmi3ValueReference references[], size_t count,
                          fmi3Float64 values[], size_t valueCount) {
  if (!value || !instance(value)->initialized || count != valueCount ||
      !validReferences(references, count, 12))
    return fmi3Error;
  auto *model = instance(value);
  double rates[3];
  derivative(model, rates);
  for (size_t index = 0; index < count; ++index) {
    const auto reference = references[index];
    if (reference < 3)
      values[index] = model->voltage[reference];
    else if (reference < 6)
      values[index] = model->current[reference - 3];
    else if (reference < 9)
      values[index] = model->current[reference - 6];
    else if (reference == 9)
      values[index] = model->time;
    else
      values[index] = rates[reference - 10];
  }
  return fmi3OK;
}

fmi3Status fmi3SetTime(fmi3Instance value, fmi3Float64 time) {
  if (!value || !instance(value)->initialized || !std::isfinite(time))
    return fmi3Error;
  instance(value)->time = time;
  return fmi3OK;
}

fmi3Status fmi3GetNumberOfContinuousStates(fmi3Instance, size_t *count) {
  if (!count)
    return fmi3Error;
  *count = 3;
  return fmi3OK;
}

fmi3Status fmi3GetNumberOfEventIndicators(fmi3Instance, size_t *count) {
  if (!count)
    return fmi3Error;
  *count = 0;
  return fmi3OK;
}

fmi3Status fmi3GetEventIndicators(fmi3Instance, fmi3Float64[], size_t count) {
  return count == 0 ? fmi3OK : fmi3Error;
}

fmi3Status fmi3GetContinuousStates(fmi3Instance value, fmi3Float64 states[],
                                   size_t count) {
  if (!value || count != 3)
    return fmi3Error;
  std::copy(instance(value)->current.begin(), instance(value)->current.end(),
            states);
  return fmi3OK;
}

fmi3Status fmi3SetContinuousStates(fmi3Instance value,
                                   const fmi3Float64 states[], size_t count) {
  if (!value || count != 3)
    return fmi3Error;
  std::copy(states, states + 3, instance(value)->current.begin());
  return fmi3OK;
}

fmi3Status fmi3GetContinuousStateDerivatives(fmi3Instance value,
                                             fmi3Float64 states[],
                                             size_t count) {
  if (!value || count != 3)
    return fmi3Error;
  derivative(instance(value), states);
  return fmi3OK;
}

fmi3Status fmi3DoStep(fmi3Instance value, fmi3Float64 currentTime,
                      fmi3Float64 stepSize, fmi3Boolean,
                      fmi3Boolean *eventHandlingNeeded,
                      fmi3Boolean *terminateSimulation,
                      fmi3Boolean *earlyReturn,
                      fmi3Float64 *lastSuccessfulTime) {
  if (!value || !instance(value)->initialized || !eventHandlingNeeded ||
      !terminateSimulation || !earlyReturn || !lastSuccessfulTime ||
      !std::isfinite(currentTime) || !std::isfinite(stepSize) || stepSize <= 0)
    return fmi3Error;
  auto *model = instance(value);
  if (std::abs(currentTime - model->time) >
      1e-10 * std::max(1.0, std::abs(model->time)))
    return fmi3Error;
  constexpr double resistance = 8.0;
  constexpr double inductance = 0.02;
  const double decay = std::exp(-resistance * stepSize / inductance);
  for (size_t phase = 0; phase < 3; ++phase)
    model->current[phase] = model->current[phase] * decay +
                            model->voltage[phase] / resistance * (1 - decay);
  model->time += stepSize;
  *eventHandlingNeeded = fmi3False;
  *terminateSimulation = fmi3False;
  *earlyReturn = fmi3False;
  *lastSuccessfulTime = model->time;
  return fmi3OK;
}

fmi3Status fmi3CompletedIntegratorStep(fmi3Instance value, fmi3Boolean,
                                       fmi3Boolean *eventHandlingNeeded,
                                       fmi3Boolean *terminateSimulation) {
  if (!value || !eventHandlingNeeded || !terminateSimulation)
    return fmi3Error;
  *eventHandlingNeeded = fmi3False;
  *terminateSimulation = fmi3False;
  return fmi3OK;
}

fmi3Status fmi3GetNominalsOfContinuousStates(fmi3Instance, fmi3Float64 values[],
                                             size_t count) {
  if (!values || count != 3)
    return fmi3Error;
  std::fill(values, values + count, 1.0);
  return fmi3OK;
}
} // extern "C"
