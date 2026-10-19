// SPDX-License-Identifier: MPL-2.0
#include <algorithm>
#include <cmath>
#include <dpsim-fmi/LinearModel.h>
#include <dpsim-fmi/ModelExchange.h>
#include <limits>
#include <set>
#include <stdexcept>
namespace DPsim::FMI {
namespace {
using CPS::Matrix;
[[noreturn]] void fail(const std::string &s) {
  throw std::runtime_error("FMI import: " + s);
}
bool approximately(const Matrix &actual, const Matrix &expected) {
  if (actual.size() == 0)
    return true;
  return (actual - expected).cwiseAbs().maxCoeff() <=
         1e-8 * (1.0 + actual.cwiseAbs().maxCoeff());
}
} // namespace

LinearModel::Ptr LinearModel::load(const std::string &path,
                                   const VoltagePort &port) {
  return loadMulti(path, {{port.voltageInput}, {port.currentOutput}, {0.0}});
}

LinearModel::Ptr LinearModel::loadMulti(const std::string &path,
                                        const VoltagePorts &port) {
  const auto channelsConfigured = port.voltageInputs.size();
  if (!channelsConfigured || channelsConfigured != port.currentOutputs.size() ||
      (!port.initialVoltages.empty() &&
       port.initialVoltages.size() != channelsConfigured))
    fail("voltage/current mappings must have equal nonzero channel counts and "
         "matching initial voltages");
  if (std::set<std::string>(port.currentOutputs.begin(),
                            port.currentOutputs.end())
          .size() != channelsConfigured)
    fail("duplicate current output mapping");
  std::vector<ScalarInput> inputs;
  for (std::size_t k = 0; k < channelsConfigured; ++k)
    inputs.push_back({port.voltageInputs[k], port.initialVoltages.empty()
                                                 ? 0.0
                                                 : port.initialVoltages[k]});
  ModelExchange runtime(path, inputs);
  const auto &description = runtime.description();
  // MATLAB R2025a uses DBL_MAX as its FMI 3 "no scheduled event" sentinel.
  const bool hasScheduledEvent =
      std::isfinite(runtime.nextTimeEvent()) &&
      runtime.nextTimeEvent() < std::numeric_limits<double>::max() / 2.0;
  if (description.hasDiscreteVariables || hasScheduledEvent)
    fail("sampled/discrete dynamics require the persistent nonlinear importer");
  LinearModel model;
  model.modelName = description.modelName;
  model.initialState = description.initialState;
  model.initialVoltage = description.initialInputs;
  model.stateNames = description.stateNames;
  model.outputNames = description.outputNames;
  model.voltageInputs = description.inputNames;
  model.voltageInput = model.voltageInputs.front();
  for (const auto &name : port.currentOutputs) {
    const auto it =
        std::find(model.outputNames.begin(), model.outputNames.end(), name);
    if (it == model.outputNames.end())
      fail("mapped current output does not exist: " + name);
    model.currentOutputs.push_back(it - model.outputNames.begin());
  }
  model.currentOutput = model.currentOutputs.front();
  const auto n = model.initialState.rows();
  const auto outputs = model.outputNames.size();
  const auto channels = model.voltageInputs.size();
  model.A.resize(n, n);
  model.B.resize(n, channels);
  model.C.resize(outputs, n);
  model.D.resize(outputs, channels);
  Matrix zero = Matrix::Zero(n, 1);
  const Matrix zeroInput = Matrix::Zero(channels, 1);
  const auto origin = runtime.evaluate(zero, zeroInput);
  if (!approximately(origin.first, zero) ||
      !approximately(origin.second, Matrix::Zero(outputs, 1)))
    fail("affine offsets are unsupported by the homogeneous linear SSN "
         "importer");
  for (Eigen::Index column = 0; column < n; ++column) {
    Matrix state = zero;
    state(column, 0) = 1;
    const auto positive = runtime.evaluate(state, zeroInput);
    const auto negative = runtime.evaluate(-state, zeroInput);
    model.A.col(column) = (positive.first - negative.first) / 2;
    model.C.col(column) = (positive.second - negative.second) / 2;
  }
  for (std::size_t column = 0; column < channels; ++column) {
    Matrix input = zeroInput;
    input(column, 0) = 1;
    const auto positive = runtime.evaluate(zero, input);
    const auto negative = runtime.evaluate(zero, -input);
    model.B.col(column) = (positive.first - negative.first) / 2;
    model.D.col(column) = (positive.second - negative.second) / 2;
  }
  // Mixed states, signs, amplitudes, initialized states and times guard against
  // accidentally treating a local linearization as a reusable LTI model.
  for (int probe = 0; probe < 10; ++probe) {
    Matrix state(n, 1);
    for (Eigen::Index row = 0; row < n; ++row)
      state(row, 0) = std::sin((row + 1) * (probe + 0.37)) * (probe + 0.5);
    if (probe == 0)
      state = model.initialState;
    Matrix input(channels, 1);
    for (std::size_t channel = 0; channel < channels; ++channel)
      input(channel, 0) =
          std::cos(probe + 0.29 + channel * 1.13) * (probe + 1) * 10;
    if (probe == 0)
      input = model.initialVoltage;
    const auto actual = runtime.evaluate(state, input, probe * 0.013);
    if (!approximately(actual.first, model.A * state + model.B * input) ||
        !approximately(actual.second, model.C * state + model.D * input))
      fail("FMU failed LTI validation; nonlinear, time-varying or hidden-state "
           "behavior is unsupported");
  }
  return std::make_shared<const LinearModel>(std::move(model));
}
} // namespace DPsim::FMI
