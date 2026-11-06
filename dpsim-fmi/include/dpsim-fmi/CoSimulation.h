// SPDX-License-Identifier: MPL-2.0
#pragma once

#include <dpsim-fmi/ModelExchange.h>
#include <memory>
#include <string>
#include <vector>

namespace DPsim::FMI {

/// Options for importing an FMI 2.0 or FMI 3.0 Co-Simulation FMU.
struct CoSimulationOptions {
  bool loggingOn = false;
  /// Reinitialize the FMU with terminal voltages transferred from the
  /// connected network's initialized nodes (for example after
  /// SystemTopology::initWithPowerflow). The FMU's own initialization
  /// equations determine its internal state; additional operating-point
  /// inputs must be supplied through the normal FMU input configuration.
  /// Requires the FMU to export Reset.
  bool initializeFromPowerflow = false;
};

/// Persistent FMI Co-Simulation instance. Inputs are set at a communication
/// point and `doStep` advances the FMU exactly once over the requested
/// interval. The class intentionally exposes communication-step semantics;
/// it does not pretend that an FMU's internal continuous states are available
/// to DPsim's MNA Newton solver.
class CoSimulation final {
public:
  class Backend;
  using Ptr = std::shared_ptr<CoSimulation>;

private:
  std::unique_ptr<Backend> mImpl;
  bool mInitializeFromPowerflow = false;

public:
  explicit CoSimulation(const std::string &path,
                        CoSimulationOptions options = {});
  CoSimulation(const std::string &path, const std::vector<ScalarInput> &inputs,
               CoSimulationOptions options = {});
  ~CoSimulation();
  CoSimulation(const CoSimulation &) = delete;
  CoSimulation &operator=(const CoSimulation &) = delete;

  const ModelDescription &description() const;
  bool initializeFromPowerflow() const { return mInitializeFromPowerflow; }
  const CPS::Matrix &inputs() const;
  const CPS::Matrix &outputs() const;
  double time() const;
  void setInputs(const CPS::Matrix &inputs);
  void setInput(const std::string &name, double value);
  /// Reset and re-enter FMI initialization mode with a new vector of inputs.
  /// This is intended for network-derived initial conditions before stepping.
  void reinitializeFromInputs(const CPS::Matrix &inputs);
  /// Advance from currentCommunicationPoint by communicationStepSize. The
  /// next output vector is returned in modelDescription output order.
  const CPS::Matrix &doStep(double currentCommunicationPoint,
                            double communicationStepSize);
};

} // namespace DPsim::FMI
