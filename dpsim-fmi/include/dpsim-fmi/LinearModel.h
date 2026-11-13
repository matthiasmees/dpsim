// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-models/Definitions.h>
#include <memory>
#include <string>
#include <vector>

namespace DPsim::FMI {

/// Explicit electrical interpretation; FMI variable names alone do not define ports.
struct VoltagePort {
  std::string voltageInput;
  std::string currentOutput;
};

/// Ordered voltage/current channels, e.g. {a,b,c}. Input names are mapped by
/// name, independently of XML declaration order. Initial voltages are
/// instantaneous SI values at t=0; an empty vector means all zero.
struct VoltagePorts {
  std::vector<std::string> voltageInputs;
  std::vector<std::string> currentOutputs;
  std::vector<double> initialVoltages;
};

/// A continuous, homogeneous LTI model recovered from FMI function evaluations.
/// Loading is an assertion by the caller that the FMU is LTI. Numerical probes
/// catch incompatible models but cannot prove global linearity.
struct LinearModel {
  using Ptr = std::shared_ptr<const LinearModel>;
  CPS::Matrix A, B, C, D;
  CPS::Matrix initialState;
  CPS::Matrix initialVoltage;
  std::string modelName;
  std::string voltageInput;
  std::vector<std::string> stateNames;
  std::vector<std::string> outputNames;
  std::size_t currentOutput = 0;
  std::vector<std::string> voltageInputs;
  std::vector<std::size_t> currentOutputs;

  /// Import scalar Float64 voltage input and scalar Float64 outputs. Continuous
  /// states are discovered at runtime. Eventful, affine and nonlinear models
  /// are unsupported. No FMI doStep or external integration solver is used.
  static Ptr load(const std::string &fmuPath, const VoltagePort &port);
  /// Recover a multi-input model. The LinearSSN electrical adapter supports
  /// one channel or three ordered abc channels, including cross-phase coupling.
  static Ptr loadMulti(const std::string &fmuPath, const VoltagePorts &ports);
};
} // namespace DPsim::FMI
