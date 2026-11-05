// SPDX-License-Identifier: MPL-2.0
// FMI 3 dual-interface averaged grid-following converter used as the IEEE39
// bus-30 generator replacement. Current outputs use passive sign convention.
#include <algorithm>
#include <array>
#include <cmath>
#include <fmi3Functions.h>
#include <new>

namespace {
struct Instance {
  std::array<double, 3> voltage{};
  std::array<double, 3> current{};
  double time = 0.0;
  bool initialized = false;
};

Instance *instance(fmi3Instance value) { return static_cast<Instance *>(value); }

bool validReferences(const fmi3ValueReference references[], size_t count,
                     fmi3ValueReference maximum) {
  for (size_t i = 0; i < count; ++i)
    if (references[i] > maximum)
      return false;
  return true;
}

std::array<double, 3> passiveCurrentReference(const Instance *model) {
  // Clarke transform with amplitude-invariant phase conventions.
  constexpr double root3over2 = 0.86602540378443864676;
  const double alpha = (2.0 / 3.0) *
                       (model->voltage[0] - 0.5 * model->voltage[1] -
                        0.5 * model->voltage[2]);
  const double beta = (2.0 / 3.0) * root3over2 *
                      (model->voltage[1] - model->voltage[2]);
  const double voltageSquared = alpha * alpha + beta * beta;
  if (!std::isfinite(voltageSquared) || voltageSquared < 1.0)
    return {0.0, 0.0, 0.0};

  // IEEE39 generator 30 operating point, interpreted as positive power
  // injection. The FMU reports the opposite sign (current into the device).
  constexpr double activePowerW = 125e6;
  constexpr double reactivePowerVar = 80.881e6;
  constexpr double nominalLineVoltageV = 345e3;
  constexpr double currentLimitPu = 1.20;
  constexpr double responseTimeS = 5e-3;
  (void)responseTimeS;
  const double factor = 2.0 / (3.0 * voltageSquared);
  double injectionAlpha = factor * (activePowerW * alpha +
                                    reactivePowerVar * beta);
  double injectionBeta = factor * (activePowerW * beta -
                                   reactivePowerVar * alpha);
  const double ratedPeak = std::sqrt(2.0) *
                           std::hypot(activePowerW, reactivePowerVar) /
                           (std::sqrt(3.0) * nominalLineVoltageV);
  const double magnitude = std::hypot(injectionAlpha, injectionBeta);
  const double maximum = currentLimitPu * ratedPeak;
  if (magnitude > maximum) {
    const double scale = maximum / magnitude;
    injectionAlpha *= scale;
    injectionBeta *= scale;
  }
  const double ia = injectionAlpha;
  const double ib = -0.5 * injectionAlpha + root3over2 * injectionBeta;
  const double ic = -0.5 * injectionAlpha - root3over2 * injectionBeta;
  return {-ia, -ib, -ic};
}

std::array<double, 3> derivatives(const Instance *model) {
  constexpr double responseTimeS = 5e-3;
  const auto target = passiveCurrentReference(model);
  std::array<double, 3> result{};
  for (size_t phase = 0; phase < 3; ++phase)
    result[phase] = (target[phase] - model->current[phase]) / responseTimeS;
  return result;
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
                                       fmi3Float64, fmi3Float64 startTime,
                                       fmi3Boolean, fmi3Float64) {
  if (!value || !std::isfinite(startTime))
    return fmi3Error;
  auto *model = instance(value);
  model->time = startTime;
  model->current.fill(0.0);
  model->initialized = false;
  return fmi3OK;
}

fmi3Status fmi3ExitInitializationMode(fmi3Instance value) {
  if (!value)
    return fmi3Error;
  auto *model = instance(value);
  model->current = passiveCurrentReference(model);
  model->initialized = true;
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
  *nextTime = 0.0;
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
  if (!value || !references || !values || count != valueCount ||
      !validReferences(references, count, 2))
    return fmi3Error;
  auto *model = instance(value);
  for (size_t i = 0; i < count; ++i) {
    if (!std::isfinite(values[i]))
      return fmi3Error;
    model->voltage[references[i]] = values[i];
  }
  return fmi3OK;
}

fmi3Status fmi3GetFloat64(fmi3Instance value,
                          const fmi3ValueReference references[], size_t count,
                          fmi3Float64 values[], size_t valueCount) {
  if (!value || !instance(value)->initialized || !references || !values ||
      count != valueCount || !validReferences(references, count, 12))
    return fmi3Error;
  const auto *model = instance(value);
  const auto rates = derivatives(model);
  for (size_t i = 0; i < count; ++i) {
    const auto reference = references[i];
    if (reference < 3)
      values[i] = model->voltage[reference];
    else if (reference < 6)
      values[i] = model->current[reference - 3];
    else if (reference < 9)
      values[i] = model->current[reference - 6];
    else if (reference == 9)
      values[i] = model->time;
    else
      values[i] = rates[reference - 10];
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
  if (!value || !states || count != 3)
    return fmi3Error;
  std::copy(instance(value)->current.begin(), instance(value)->current.end(),
            states);
  return fmi3OK;
}
fmi3Status fmi3SetContinuousStates(fmi3Instance value,
                                   const fmi3Float64 states[], size_t count) {
  if (!value || !states || count != 3)
    return fmi3Error;
  std::copy(states, states + 3, instance(value)->current.begin());
  return fmi3OK;
}
fmi3Status fmi3GetContinuousStateDerivatives(fmi3Instance value,
                                             fmi3Float64 states[],
                                             size_t count) {
  if (!value || !states || count != 3)
    return fmi3Error;
  const auto rates = derivatives(instance(value));
  std::copy(rates.begin(), rates.end(), states);
  return fmi3OK;
}
fmi3Status fmi3GetNominalsOfContinuousStates(fmi3Instance, fmi3Float64 values[],
                                             size_t count) {
  if (!values || count != 3)
    return fmi3Error;
  std::fill(values, values + count, 1.0);
  return fmi3OK;
}
fmi3Status fmi3CompletedIntegratorStep(fmi3Instance value, fmi3Boolean,
                                       fmi3Boolean *eventHandlingNeeded,
                                       fmi3Boolean *terminateSimulation) {
  if (!value || !eventHandlingNeeded || !terminateSimulation)
    return fmi3Error;
  *eventHandlingNeeded = *terminateSimulation = fmi3False;
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
  constexpr double responseTimeS = 5e-3;
  const auto target = passiveCurrentReference(model);
  const double decay = std::exp(-stepSize / responseTimeS);
  for (size_t phase = 0; phase < 3; ++phase)
    model->current[phase] = target[phase] +
                            (model->current[phase] - target[phase]) * decay;
  model->time += stepSize;
  *eventHandlingNeeded = *terminateSimulation = *earlyReturn = fmi3False;
  *lastSuccessfulTime = model->time;
  return fmi3OK;
}
} // extern "C"
