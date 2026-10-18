// SPDX-License-Identifier: MPL-2.0
#pragma once
#include <dpsim-models/Definitions.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace DPsim::FMI {
struct ScalarInput {
  std::string name;
  double initialValue;
};
struct ModelExchangeOptions {
  /// Detect FMI event-indicator sign changes at accepted fixed-step endpoints.
  /// Event times are therefore quantized to the DPsim timestep.
  bool endpointStateEvents = false;
  /// Accept continuous-state values reinitialized by FMI event iteration.
  bool allowContinuousStateResets = false;
};
struct ModelDescription {
  std::string modelName;
  std::string modelIdentifier, instantiationToken;
  bool needsCompletedIntegratorStep = false;
  bool canGetAndSetFMUState = false;
  bool canSerializeFMUState = false;
  bool canBeInstantiatedOnlyOncePerProcess = false;
  bool hasEventIndicators = false;
  bool hasDiscreteVariables = false;
  std::vector<std::string> stateNames, inputNames, outputNames;
  CPS::Matrix initialState, initialInputs, stateNominals;
};

/// Persistent in-process FMI 3 Model Exchange instance. Evaluations set all
/// continuous states and inputs; accepted-step/event methods commit history.
/// Each instance must be owned by one component (no concurrent evaluations).
class ModelExchange final {
public:
  class Backend;

private:
  std::unique_ptr<Backend> mImpl;

public:
  using Ptr = std::shared_ptr<ModelExchange>;
  ModelExchange(const std::string &path, const std::vector<ScalarInput> &inputs,
                ModelExchangeOptions options = {});
  /// Import every continuous scalar Float64 input using its declared start.
  explicit ModelExchange(const std::string &path,
                         ModelExchangeOptions options = {});
  ~ModelExchange();
  ModelExchange(const ModelExchange &) = delete;
  ModelExchange &operator=(const ModelExchange &) = delete;
  const ModelDescription &description() const;
  std::pair<CPS::Matrix, CPS::Matrix> evaluate(const CPS::Matrix &state,
                                               const CPS::Matrix &inputs,
                                               double time = 0);
  /// Absolute next internal time event, or positive infinity when none.
  double nextTimeEvent() const;
  /// Cached right-side values at the last committed point (no new FMI probes).
  const std::pair<CPS::Matrix, CPS::Matrix> &acceptedValues() const;
  /// Actual continuous state after the last accepted FMI event iteration.
  const CPS::Matrix &acceptedState() const;
  /// Whether the exporter supplies rollback snapshots. Multirate integration
  /// needs this to replay a macro step for each network Newton trial.
  bool supportsRollback() const;
  /// Save, repeatedly restore, and finally discard one accepted FMI point.
  /// Checkpoints include discrete history, continuous state, time-event state,
  /// and the importer's accepted-value caches.
  void saveCheckpoint();
  void restoreCheckpoint();
  void discardCheckpoint();
  /// Select an accepted endpoint without querying provisional outputs.
  void prepareStep(const CPS::Matrix &state, const CPS::Matrix &inputs,
                   double time);
  void completeStep();
  /// Apply right-side inputs in the same event iteration as a coincident time event.
  void completeStep(const CPS::Matrix &eventInputs);
  /// Scheduled external input discontinuity at an accepted endpoint. State
  /// reinitialization and state events remain unsupported.
  void inputEvent(const CPS::Matrix &state, const CPS::Matrix &inputs,
                  double time);
};
} // namespace DPsim::FMI
