// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/ModelExchange.h>
#include <fmi2FunctionTypes.h>
#include <fmi3FunctionTypes.h>
#include <libxml/parser.h>
#include <miniz.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace DPsim::FMI {
namespace {
using CPS::Matrix;
namespace fs = std::filesystem;

[[noreturn]] void fail(const std::string &message) {
  throw std::runtime_error("FMI import: " + message);
}

struct Archive {
  fs::path directory;
  Archive() {
    std::random_device random;
    for (int attempt = 0; attempt < 100; ++attempt) {
      auto path =
          fs::temp_directory_path() / ("dpsim-fmi-" + std::to_string(random()) +
                                       "-" + std::to_string(random()));
      if (fs::create_directory(path)) {
        directory = path;
        fs::permissions(directory, fs::perms::owner_all);
        return;
      }
    }
    fail("cannot create temporary extraction directory");
  }
  ~Archive() {
    std::error_code error;
    fs::remove_all(directory, error);
  }
  void extract(const std::string &path) {
    mz_zip_archive zip{};
    if (!mz_zip_reader_init_file(&zip, path.c_str(), 0))
      fail("cannot read FMU ZIP: " + path);
    struct Close {
      mz_zip_archive *zip;
      ~Close() { mz_zip_reader_end(zip); }
    } close{&zip};
    const auto count = mz_zip_reader_get_num_files(&zip);
    if (count > 10000)
      fail("archive has too many entries");
    std::set<std::string> paths;
    mz_uint64 total = 0;
    for (mz_uint index = 0; index < count; ++index) {
      mz_zip_archive_file_stat stat{};
      if (!mz_zip_reader_file_stat(&zip, index, &stat))
        fail("cannot inspect ZIP entry");
      std::string name = stat.m_filename;
      if (name.size() == sizeof(stat.m_filename) - 1)
        fail("archive entry name is too long");
      fs::path relative(name);
      if (name.empty() || name.find('\\') != std::string::npos ||
          name.find(':') != std::string::npos || relative.is_absolute())
        fail("unsafe archive path: " + name);
      for (const auto &part : relative)
        if (part == ".." || part == ".")
          fail("unsafe archive path: " + name);
      if (!paths.insert(name).second)
        fail("duplicate archive entry: " + name);
      const auto kind = (stat.m_external_attr >> 16) & 0170000;
      if (kind != 0 && kind != 0100000 && kind != 0040000)
        fail("archive links and special files are unsupported");
      total += stat.m_uncomp_size;
      if (stat.m_uncomp_size > 512ULL * 1024 * 1024 ||
          total > 1024ULL * 1024 * 1024)
        fail("archive exceeds extraction size limit");
      const auto destination = directory / relative;
      if (mz_zip_reader_is_file_a_directory(&zip, index)) {
        fs::create_directories(destination);
      } else {
        fs::create_directories(destination.parent_path());
        if (!mz_zip_reader_extract_to_file(&zip, index,
                                           destination.string().c_str(), 0))
          fail("cannot extract ZIP entry: " + name);
      }
    }
  }
};

std::string property(xmlNode *node, const char *name,
                     const std::string &fallback = "") {
  xmlChar *value = xmlGetProp(node, BAD_CAST name);
  if (!value)
    return fallback;
  std::string result(reinterpret_cast<char *>(value));
  xmlFree(value);
  return result;
}
bool tag(xmlNode *node, const char *name) {
  return node->type == XML_ELEMENT_NODE &&
         xmlStrEqual(node->name, BAD_CAST name);
}
bool flag(xmlNode *node, const char *name) {
  const auto value = property(node, name, "false");
  if (value == "true" || value == "1")
    return true;
  if (value == "false" || value == "0")
    return false;
  fail(std::string("invalid Boolean capability: ") + name);
}
xmlNode *child(xmlNode *node, const char *name) {
  for (auto entry = node->children; entry; entry = entry->next)
    if (tag(entry, name))
      return entry;
  fail(std::string("missing XML element: ") + name);
}
fmi3ValueReference reference(const std::string &text) {
  std::size_t consumed = 0;
  const auto value = std::stoull(text, &consumed);
  if (consumed != text.size() || value > UINT32_MAX)
    fail("invalid value reference: " + text);
  return static_cast<fmi3ValueReference>(value);
}

struct Variable {
  std::string name, type, causality, derivative, variability, start;
  fmi3ValueReference ref;
  bool scalar;
};

} // namespace

class ModelExchange::Backend {
public:
  ModelDescription model;
  virtual ~Backend() = default;
  virtual std::pair<Matrix, Matrix> evaluate(const Matrix &, const Matrix &,
                                             double) = 0;
  virtual double nextTimeEvent() const = 0;
  virtual const std::pair<Matrix, Matrix> &acceptedValues() const = 0;
  virtual const Matrix &acceptedState() const = 0;
  virtual bool supportsRollback() const = 0;
  virtual void saveCheckpoint() = 0;
  virtual void restoreCheckpoint() = 0;
  virtual void discardCheckpoint() = 0;
  virtual void prepareStep(const Matrix &, const Matrix &, double) = 0;
  virtual void completeStep(const Matrix *eventInputs = nullptr) = 0;
  virtual void inputEvent(const Matrix &, const Matrix &, double) = 0;
};

class Fmi3Backend final : public ModelExchange::Backend {
  Archive mArchive;
#ifdef _WIN32
  HMODULE mLibrary = nullptr;
#else
  void *mLibrary = nullptr;
#endif
  fmi3Instance mInstance = nullptr;
  bool mContinuous = false;
  bool mValid = true;
  bool mNeedsCompletedIntegratorStep = false;
  bool mUseTrialSnapshots = false;
  ModelExchangeOptions mOptions;
  double mTime = 0, mAcceptedTime = 0;
  Matrix mLastState, mLastInputs, mAcceptedState, mAcceptedInputs;
  Matrix mAcceptedIndicators;
  std::pair<Matrix, Matrix> mAcceptedValues;
  double mNextTimeEvent = std::numeric_limits<double>::infinity();
  std::string mLastMessage;
  std::vector<fmi3ValueReference> mInputs;
  std::vector<fmi3ValueReference> mOutputs;
  std::size_t mEventIndicatorMetadata = 0;
  struct Checkpoint {
    fmi3FMUState runtime = nullptr;
    bool continuous = false;
    bool useTrialSnapshots = false;
    double time = 0;
    double acceptedTime = 0;
    double nextTimeEvent = std::numeric_limits<double>::infinity();
    Matrix lastState, lastInputs, acceptedState, acceptedInputs;
    Matrix acceptedIndicators, stateNominals;
    std::pair<Matrix, Matrix> acceptedValues;
  };
  std::unique_ptr<Checkpoint> mCheckpoint;

#define FMI_FUNCTION(name) fmi3##name##TYPE *name = nullptr
  FMI_FUNCTION(InstantiateModelExchange);
  FMI_FUNCTION(GetVersion);
  FMI_FUNCTION(FreeInstance);
  FMI_FUNCTION(Terminate);
  FMI_FUNCTION(EnterInitializationMode);
  FMI_FUNCTION(ExitInitializationMode);
  FMI_FUNCTION(UpdateDiscreteStates);
  FMI_FUNCTION(EnterContinuousTimeMode);
  FMI_FUNCTION(SetTime);
  FMI_FUNCTION(SetFloat64);
  FMI_FUNCTION(GetFloat64);
  FMI_FUNCTION(GetContinuousStates);
  FMI_FUNCTION(SetContinuousStates);
  FMI_FUNCTION(GetContinuousStateDerivatives);
  FMI_FUNCTION(GetNumberOfContinuousStates);
  FMI_FUNCTION(GetNumberOfEventIndicators);
  FMI_FUNCTION(GetEventIndicators);
  FMI_FUNCTION(CompletedIntegratorStep);
  FMI_FUNCTION(EnterEventMode);
  FMI_FUNCTION(GetNominalsOfContinuousStates);
  FMI_FUNCTION(GetFMUState);
  FMI_FUNCTION(SetFMUState);
  FMI_FUNCTION(FreeFMUState);
#undef FMI_FUNCTION

  template <class T> T symbol(const char *name, bool required = true) {
#ifdef _WIN32
    auto address = GetProcAddress(mLibrary, name);
#else
    auto address = dlsym(mLibrary, name);
#endif
    if (!address && required)
      fail(std::string("missing function: ") + name);
    return reinterpret_cast<T>(address);
  }
  void check(fmi3Status status, const char *operation) {
    if (status > fmi3Warning) {
      mValid = false;
      fail(std::string(operation) + " failed (status " +
           std::to_string(status) + "): " + mLastMessage);
    }
  }
  static void log(fmi3InstanceEnvironment environment, fmi3Status, fmi3String,
                  fmi3String message) noexcept {
    try {
      static_cast<Fmi3Backend *>(environment)->mLastMessage =
          message ? message : "";
    } catch (...) {
      // Never allow a C++ exception across the FMI C callback boundary.
    }
  }
  void close() noexcept {
    if (mInstance) {
      if (mCheckpoint && mCheckpoint->runtime && FreeFMUState)
        FreeFMUState(mInstance, &mCheckpoint->runtime);
      mCheckpoint.reset();
      if (mContinuous && mValid && Terminate)
        Terminate(mInstance);
      if (FreeInstance)
        FreeInstance(mInstance);
      mInstance = nullptr;
    }
    if (mLibrary) {
#ifdef _WIN32
      FreeLibrary(mLibrary);
#else
      dlclose(mLibrary);
#endif
      mLibrary = nullptr;
    }
  }
  static std::string binaryDirectory() {
#if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64)
    const std::string architecture = "aarch64";
#elif defined(__x86_64__) || defined(_M_X64)
    const std::string architecture = "x86_64";
#else
    fail("unsupported host architecture");
#endif
#ifdef __APPLE__
    return architecture + "-darwin";
#elif defined(_WIN32)
    return architecture + "-windows";
#elif defined(__linux__)
    return architecture + "-linux";
#else
    fail("unsupported host operating system");
#endif
  }

public:
  Fmi3Backend(const std::string &path,
              const std::vector<ScalarInput> &configuration,
              ModelExchangeOptions options)
      : mOptions(options) {
    try {
      initialize(path, configuration);
    } catch (...) {
      close();
      throw;
    }
  }
  ~Fmi3Backend() override { close(); }
  Fmi3Backend(const Fmi3Backend &) = delete;
  Fmi3Backend &operator=(const Fmi3Backend &) = delete;

  void initialize(const std::string &path,
                  const std::vector<ScalarInput> &configuration) {
    const bool automaticInputs = configuration.empty();
    const auto channels = configuration.size();
    if (channels > 64)
      fail("input count exceeds maximum 64");
    model.initialInputs = Matrix::Zero(channels, 1);
    std::set<std::string> configuredNames;
    for (std::size_t i = 0; i < channels; ++i) {
      if (configuration[i].name.empty() ||
          !configuredNames.insert(configuration[i].name).second)
        fail("empty or duplicate input mapping");
      model.inputNames.push_back(configuration[i].name);
      model.initialInputs(i, 0) = configuration[i].initialValue;
    }
    if (!model.initialInputs.allFinite())
      fail("initial inputs must be finite");
    mInputs.resize(channels);
    std::vector<bool> foundInputs(channels, false);
    std::vector<double> discoveredInputStarts;
    mArchive.extract(path);
    const auto xmlPath = mArchive.directory / "modelDescription.xml";
    std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> document(
        xmlReadFile(xmlPath.string().c_str(), nullptr, XML_PARSE_NONET),
        &xmlFreeDoc);
    if (!document || document->intSubset || document->extSubset)
      fail("invalid modelDescription.xml (DTDs are unsupported)");
    auto root = xmlDocGetRootElement(document.get());
    if (!root || !tag(root, "fmiModelDescription") ||
        property(root, "fmiVersion") != "3.0")
      fail("requires an FMI 3.0 model description");
    auto me = child(root, "ModelExchange");
    mNeedsCompletedIntegratorStep = flag(me, "needsCompletedIntegratorStep");
    model.needsCompletedIntegratorStep = mNeedsCompletedIntegratorStep;
    model.canGetAndSetFMUState = flag(me, "canGetAndSetFMUState");
    model.canSerializeFMUState = flag(me, "canSerializeFMUState");
    model.canBeInstantiatedOnlyOncePerProcess =
        flag(me, "canBeInstantiatedOnlyOncePerProcess");
    if (flag(me, "needsExecutionTool"))
      fail("FMUs requiring an external execution tool are unsupported");
    const auto identifier = property(me, "modelIdentifier");
    const auto token = property(root, "instantiationToken");
    if (identifier.empty() ||
        identifier.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNO"
                                     "PQRSTUVWXYZ0123456789_") !=
            std::string::npos ||
        token.empty())
      fail("invalid modelIdentifier or instantiationToken");
    model.modelName = property(root, "modelName");
    model.modelIdentifier = identifier;
    model.instantiationToken = token;
    if (model.modelName.empty())
      fail("modelName is missing");
    std::map<fmi3ValueReference, Variable> variables;
    std::set<std::string> names;
    auto modelVariables = child(root, "ModelVariables");
    for (auto v = modelVariables->children; v; v = v->next) {
      if (v->type != XML_ELEMENT_NODE)
        continue;
      Variable variable{property(v, "name"),
                        reinterpret_cast<const char *>(v->name),
                        property(v, "causality", "local"),
                        property(v, "derivative"),
                        property(v, "variability", "continuous"),
                        property(v, "start", "0"),
                        reference(property(v, "valueReference")),
                        true};
      for (auto dimension = v->children; dimension; dimension = dimension->next)
        if (tag(dimension, "Dimension"))
          variable.scalar = false;
      if (!names.insert(variable.name).second)
        fail("duplicate FMI variable name");
      // FMI permits aliases: several variables may intentionally share one
      // value reference. Keep the first descriptor for derivative metadata;
      // names and the input/output reference lists retain every variable.
      variables.emplace(variable.ref, variable);
      if (variable.type == "Clock")
        fail("clocked models are unsupported");
      if (variable.variability == "discrete")
        model.hasDiscreteVariables = true;
      if (variable.causality == "structuralParameter")
        fail("structural parameters are unsupported");
      // Parameters retain exported defaults, including parameters marked tunable.
      // Internal discrete variables remain owned and updated by FMI event mode.
      if (variable.causality == "input" || variable.causality == "output") {
        if (variable.type != "Float64" || !variable.scalar ||
            variable.variability != "continuous")
          fail("requires continuous scalar Float64 inputs and outputs");
        if (variable.causality == "input") {
          if (automaticInputs) {
            std::size_t consumed = 0;
            const double start = std::stod(variable.start, &consumed);
            if (consumed != variable.start.size() || !std::isfinite(start))
              fail("invalid Float64 input start for " + variable.name);
            model.inputNames.push_back(variable.name);
            mInputs.push_back(variable.ref);
            discoveredInputStarts.push_back(start);
            continue;
          }
          const auto found = std::find(model.inputNames.begin(),
                                       model.inputNames.end(), variable.name);
          if (found == model.inputNames.end())
            fail("all inputs must be mapped; unexpected input " +
                 variable.name);
          const auto index =
              static_cast<std::size_t>(found - model.inputNames.begin());
          mInputs[index] = variable.ref;
          foundInputs[index] = true;
        } else {
          mOutputs.push_back(variable.ref);
          model.outputNames.push_back(variable.name);
        }
      }
    }
    if (std::find(foundInputs.begin(), foundInputs.end(), false) !=
        foundInputs.end())
      fail("mapped input does not exist");
    if (automaticInputs) {
      if (discoveredInputStarts.size() > 64)
        fail("input count exceeds maximum 64");
      model.initialInputs.resize(discoveredInputStarts.size(), 1);
      for (std::size_t index = 0; index < discoveredInputStarts.size(); ++index)
        model.initialInputs(index, 0) = discoveredInputStarts[index];
    }
    auto structure = child(root, "ModelStructure");
    for (auto entry = structure->children; entry; entry = entry->next) {
      if (tag(entry, "EventIndicator")) {
        ++mEventIndicatorMetadata;
        model.hasEventIndicators = true;
        if (!mOptions.endpointStateEvents)
          fail("state events are unsupported unless endpoint handling is "
               "explicitly enabled");
      }
      if (!tag(entry, "ContinuousStateDerivative"))
        continue;
      const auto derivative =
          variables.find(reference(property(entry, "valueReference")));
      if (derivative == variables.end() ||
          derivative->second.derivative.empty())
        fail("invalid continuous-state derivative metadata");
      const auto state =
          variables.find(reference(derivative->second.derivative));
      if (state == variables.end() || state->second.type != "Float64" ||
          !state->second.scalar || derivative->second.type != "Float64" ||
          !derivative->second.scalar)
        fail("requires scalar Float64 continuous states and derivatives");
      model.stateNames.push_back(state->second.name);
    }
    auto binary =
        mArchive.directory / "binaries" / binaryDirectory() / identifier;
#ifdef __APPLE__
    binary += ".dylib";
#elif defined(_WIN32)
    binary += ".dll";
#else
    binary += ".so";
#endif
    if (!fs::is_regular_file(binary))
      fail("missing host binary: " +
           binary.lexically_relative(mArchive.directory).string());
#ifdef _WIN32
    mLibrary = LoadLibraryW(binary.wstring().c_str());
#else
    mLibrary = dlopen(binary.string().c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
    if (!mLibrary) {
#ifdef _WIN32
      fail("cannot load FMU binary (error " + std::to_string(GetLastError()) +
           ")");
#else
      fail(std::string("cannot load FMU binary: ") + dlerror());
#endif
    }
#define LOAD(name) name = symbol<decltype(name)>("fmi3" #name)
    LOAD(InstantiateModelExchange);
    LOAD(GetVersion);
    LOAD(FreeInstance);
    LOAD(Terminate);
    LOAD(EnterInitializationMode);
    LOAD(ExitInitializationMode);
    LOAD(UpdateDiscreteStates);
    LOAD(EnterContinuousTimeMode);
    LOAD(SetTime);
    LOAD(SetFloat64);
    LOAD(GetFloat64);
    LOAD(GetContinuousStates);
    LOAD(SetContinuousStates);
    LOAD(GetContinuousStateDerivatives);
    LOAD(GetNumberOfContinuousStates);
    LOAD(GetNumberOfEventIndicators);
#undef LOAD
    GetEventIndicators =
        symbol<decltype(GetEventIndicators)>("fmi3GetEventIndicators", false);
    CompletedIntegratorStep = symbol<decltype(CompletedIntegratorStep)>(
        "fmi3CompletedIntegratorStep", false);
    EnterEventMode =
        symbol<decltype(EnterEventMode)>("fmi3EnterEventMode", false);
    GetNominalsOfContinuousStates =
        symbol<decltype(GetNominalsOfContinuousStates)>(
            "fmi3GetNominalsOfContinuousStates", false);
    if (flag(me, "canGetAndSetFMUState")) {
      GetFMUState = symbol<decltype(GetFMUState)>("fmi3GetFMUState");
      SetFMUState = symbol<decltype(SetFMUState)>("fmi3SetFMUState");
      FreeFMUState = symbol<decltype(FreeFMUState)>("fmi3FreeFMUState");
    }
    const auto version = GetVersion();
    if (!version || std::string(version) != "3.0")
      fail("binary does not implement FMI 3.0");
    // Respect this capability even after an instance has been freed.
    if (flag(me, "canBeInstantiatedOnlyOncePerProcess")) {
      static std::mutex mutex;
      static std::set<std::string> instantiated;
      std::lock_guard<std::mutex> guard(mutex);
      if (!instantiated.insert(identifier + ":" + token).second)
        fail("FMU permits only one instantiation per process; reuse the "
             "imported LinearModel; persistent nonlinear components cannot "
             "clone a singleton");
    }
    const auto resources = mArchive.directory / "resources";
    const std::string resourcePath =
        resources.string() + fs::path::preferred_separator;
    mInstance = InstantiateModelExchange(
        model.modelName.c_str(), token.c_str(),
        fs::is_directory(resources) ? resourcePath.c_str() : nullptr, fmi3False,
        fmi3False, this, log);
    if (!mInstance)
      fail("instantiation failed: " + mLastMessage);
    check(EnterInitializationMode(mInstance, fmi3False, 0, 0, fmi3False, 0),
          "enter initialization mode");
    check(SetFloat64(mInstance, mInputs.data(), mInputs.size(),
                     model.initialInputs.data(), mInputs.size()),
          "set initial voltage");
    check(ExitInitializationMode(mInstance), "exit initialization mode");
    bool settled = false;
    for (int iteration = 0; iteration < 100; ++iteration) {
      fmi3Boolean again = fmi3False, terminate = fmi3False;
      fmi3Boolean nominalsChanged = fmi3False, statesChanged = fmi3False;
      fmi3Boolean timeEvent = fmi3False;
      double nextTime = std::numeric_limits<double>::infinity();
      check(UpdateDiscreteStates(mInstance, &again, &terminate,
                                 &nominalsChanged, &statesChanged, &timeEvent,
                                 &nextTime),
            "initial event iteration");
      if (terminate)
        fail("initialization requested termination");
      if (!again) {
        updateNextTimeEvent(timeEvent, nextTime);
        settled = true;
        break;
      }
    }
    if (!settled)
      fail("initial event iteration did not converge");
    check(EnterContinuousTimeMode(mInstance), "enter continuous time mode");
    mContinuous = true;
    size_t states = 0, indicators = 0;
    check(GetNumberOfContinuousStates(mInstance, &states), "get state count");
    check(GetNumberOfEventIndicators(mInstance, &indicators),
          "get event indicator count");
    if (indicators != mEventIndicatorMetadata)
      fail("event indicator count is inconsistent with modelDescription.xml");
    if (indicators && (!mOptions.endpointStateEvents || !GetEventIndicators))
      fail("state events require explicitly enabled endpoint handling and "
           "fmi3GetEventIndicators");
    if (states > 512 || states != model.stateNames.size())
      fail(
          "state count is excessive or inconsistent with modelDescription.xml");
    model.initialState.resize(states, 1);
    if (states)
      check(GetContinuousStates(mInstance, model.initialState.data(), states),
            "get initial states");
    model.stateNominals = Matrix::Ones(states, 1);
    if (states && GetNominalsOfContinuousStates)
      check(GetNominalsOfContinuousStates(mInstance, model.stateNominals.data(),
                                          states),
            "get state nominals");
    if (!model.stateNominals.allFinite() ||
        (model.stateNominals.array() <= 0).any())
      fail("invalid continuous state nominals");
    if (!model.initialState.allFinite())
      fail("non-finite initial states");
    mLastState = mAcceptedState = model.initialState;
    mLastInputs = mAcceptedInputs = model.initialInputs;
    mAcceptedValues = readValues();
    mAcceptedIndicators = readEventIndicators();
    mUseTrialSnapshots = mOptions.snapshotTrialEvaluations && GetFMUState &&
                         (model.hasDiscreteVariables ||
                          std::isfinite(mNextTimeEvent) || indicators != 0);
  }

  Matrix readEventIndicators() {
    Matrix indicators(mEventIndicatorMetadata, 1);
    if (indicators.rows())
      check(GetEventIndicators(mInstance, indicators.data(), indicators.rows()),
            "evaluate event indicators");
    if (!indicators.allFinite())
      fail("non-finite FMI event indicator");
    return indicators;
  }
  static bool crossedEventIndicator(const Matrix &before, const Matrix &after) {
    for (Eigen::Index k = 0; k < before.rows(); ++k)
      if ((before(k, 0) > 0) != (after(k, 0) > 0))
        return true;
    return false;
  }

  std::pair<Matrix, Matrix> readValues() {
    Matrix derivatives(model.initialState.rows(), 1),
        outputs(mOutputs.size(), 1);
    check(GetFloat64(mInstance, mOutputs.data(), mOutputs.size(),
                     outputs.data(), mOutputs.size()),
          "evaluate outputs");
    if (derivatives.rows())
      check(GetContinuousStateDerivatives(mInstance, derivatives.data(),
                                          derivatives.rows()),
            "evaluate derivatives");
    if (!derivatives.allFinite() || !outputs.allFinite())
      fail("non-finite derivatives or outputs during model recovery");
    return {derivatives, outputs};
  }
  void applyPoint(const Matrix &state, const Matrix &inputs, double time) {
    check(SetTime(mInstance, time), "set evaluation time");
    if (state.rows())
      check(SetContinuousStates(mInstance, state.data(), state.rows()),
            "set trial states");
    if (mInputs.size())
      check(SetFloat64(mInstance, mInputs.data(), mInputs.size(), inputs.data(),
                       mInputs.size()),
            "set trial inputs");
  }
  void cacheAccepted(const Matrix *indicators = nullptr) {
    mAcceptedTime = mTime;
    mAcceptedState = mLastState;
    mAcceptedInputs = mLastInputs;
    mAcceptedValues = readValues();
    mAcceptedIndicators = indicators ? *indicators : readEventIndicators();
  }
  const std::pair<Matrix, Matrix> &acceptedValues() const override {
    return mAcceptedValues;
  }
  const Matrix &acceptedState() const override { return mAcceptedState; }
  bool supportsRollback() const override {
    return GetFMUState && SetFMUState && FreeFMUState;
  }
  void saveCheckpoint() override {
    if (!supportsRollback())
      fail("multirate integration requires canGetAndSetFMUState");
    if (mCheckpoint)
      fail("an FMI checkpoint is already active");
    auto checkpoint = std::make_unique<Checkpoint>();
    check(GetFMUState(mInstance, &checkpoint->runtime),
          "save multirate checkpoint");
    checkpoint->continuous = mContinuous;
    checkpoint->useTrialSnapshots = mUseTrialSnapshots;
    checkpoint->time = mTime;
    checkpoint->acceptedTime = mAcceptedTime;
    checkpoint->nextTimeEvent = mNextTimeEvent;
    checkpoint->lastState = mLastState;
    checkpoint->lastInputs = mLastInputs;
    checkpoint->acceptedState = mAcceptedState;
    checkpoint->acceptedInputs = mAcceptedInputs;
    checkpoint->acceptedIndicators = mAcceptedIndicators;
    checkpoint->stateNominals = model.stateNominals;
    checkpoint->acceptedValues = mAcceptedValues;
    mCheckpoint = std::move(checkpoint);
  }
  void restoreCheckpoint() override {
    if (!mCheckpoint || !mCheckpoint->runtime)
      fail("no active FMI checkpoint to restore");
    check(SetFMUState(mInstance, mCheckpoint->runtime),
          "restore multirate checkpoint");
    mContinuous = mCheckpoint->continuous;
    mUseTrialSnapshots = mCheckpoint->useTrialSnapshots;
    mTime = mCheckpoint->time;
    mAcceptedTime = mCheckpoint->acceptedTime;
    mNextTimeEvent = mCheckpoint->nextTimeEvent;
    mLastState = mCheckpoint->lastState;
    mLastInputs = mCheckpoint->lastInputs;
    mAcceptedState = mCheckpoint->acceptedState;
    mAcceptedInputs = mCheckpoint->acceptedInputs;
    mAcceptedIndicators = mCheckpoint->acceptedIndicators;
    model.stateNominals = mCheckpoint->stateNominals;
    mAcceptedValues = mCheckpoint->acceptedValues;
  }
  void discardCheckpoint() override {
    if (!mCheckpoint)
      return;
    if (mCheckpoint->runtime)
      check(FreeFMUState(mInstance, &mCheckpoint->runtime),
            "free multirate checkpoint");
    mCheckpoint.reset();
  }

  std::pair<Matrix, Matrix> evaluate(const Matrix &state, const Matrix &voltage,
                                     double time = 0) override {
    prepareStep(state, voltage, time);
    if (std::abs(mTime - mAcceptedTime) <= 1e-12 && state == mAcceptedState &&
        voltage == mAcceptedInputs)
      return mAcceptedValues;
    // Preserve discrete history and lazy output caches across provisional Newton
    // evaluations when the exporter supplies FMI state snapshots.
    fmi3FMUState saved = nullptr;
    if (mUseTrialSnapshots)
      check(GetFMUState(mInstance, &saved), "snapshot trial evaluation");
    try {
      applyPoint(state, voltage, mTime);
      auto result = readValues();
      if (saved) {
        check(SetFMUState(mInstance, saved), "restore trial evaluation");
        check(FreeFMUState(mInstance, &saved), "free trial snapshot");
      }
      return result;
    } catch (...) {
      if (saved) {
        if (mValid)
          SetFMUState(mInstance, saved);
        FreeFMUState(mInstance, &saved);
      }
      throw;
    }
  }
  void prepareStep(const Matrix &state, const Matrix &voltage,
                   double time) override {
    if (state.rows() != model.initialState.rows() || state.cols() != 1 ||
        voltage.rows() != model.initialInputs.rows() || voltage.cols() != 1 ||
        !state.allFinite() || !voltage.allFinite() || !std::isfinite(time))
      fail("invalid state/input/time in evaluation");
    const double tolerance = 1e-10 * std::max(1.0, std::abs(time));
    if (time > mNextTimeEvent + tolerance)
      fail("evaluation crossed a scheduled FMI time event");
    mTime =
        std::abs(time - mNextTimeEvent) <= tolerance ? mNextTimeEvent : time;
    mLastState = state;
    mLastInputs = voltage;
  }
  void updateNextTimeEvent(bool defined, double next) {
    // Some exported FMUs use the defined flag together with +inf or DBL_MAX
    // to mean that no further scheduled event exists. Treat those sentinels
    // like an undefined event, while retaining strict validation for NaN and
    // stale finite event times.
    if (defined && ((std::isinf(next) && next > 0) ||
                    next >= std::numeric_limits<double>::max() / 2.0)) {
      mNextTimeEvent = std::numeric_limits<double>::infinity();
      if (mOptions.snapshotTrialEvaluations && GetFMUState &&
          model.hasDiscreteVariables)
        mUseTrialSnapshots = true;
      return;
    }
    if (defined && (!std::isfinite(next) || next <= mTime)) {
      std::ostringstream message;
      message << std::setprecision(17)
              << "FMI next time event must advance time (current=" << mTime
              << ", next=" << next << ")";
      fail(message.str());
    }
    mNextTimeEvent = defined ? next : std::numeric_limits<double>::infinity();
    if (mOptions.snapshotTrialEvaluations && GetFMUState &&
        (defined || model.hasDiscreteVariables))
      mUseTrialSnapshots = true;
  }
  double nextTimeEvent() const override { return mNextTimeEvent; }
  void settleEvent(const char *operation, const Matrix *eventInputs = nullptr) {
    if (!EnterEventMode)
      fail("missing function: fmi3EnterEventMode required for events");
    Matrix before(model.initialState.rows(), 1);
    if (before.rows())
      check(GetContinuousStates(mInstance, before.data(), before.rows()),
            "snapshot event states");
    check(EnterEventMode(mInstance), "enter event mode");
    mContinuous = false;
    if (eventInputs) {
      if (eventInputs->rows() != model.initialInputs.rows() ||
          eventInputs->cols() != 1 || !eventInputs->allFinite())
        fail("invalid event inputs");
      check(SetFloat64(mInstance, mInputs.data(), mInputs.size(),
                       eventInputs->data(), mInputs.size()),
            "apply event inputs");
    }
    bool settled = false;
    for (int k = 0; k < 100; ++k) {
      fmi3Boolean again = fmi3False, terminate = fmi3False;
      fmi3Boolean nominals = fmi3False, values = fmi3False;
      fmi3Boolean timeEvent = fmi3False;
      double next = std::numeric_limits<double>::infinity();
      check(UpdateDiscreteStates(mInstance, &again, &terminate, &nominals,
                                 &values, &timeEvent, &next),
            operation);
      if (terminate)
        fail("event requested termination");
      if (nominals && GetNominalsOfContinuousStates && before.rows()) {
        check(GetNominalsOfContinuousStates(
                  mInstance, model.stateNominals.data(), before.rows()),
              "update state nominals");
        if (!model.stateNominals.allFinite() ||
            (model.stateNominals.array() <= 0).any())
          fail("invalid updated state nominals");
      }
      if (!again) {
        updateNextTimeEvent(timeEvent, next);
        settled = true;
        break;
      }
    }
    if (!settled)
      fail("event iteration did not converge");
    check(EnterContinuousTimeMode(mInstance), "leave event mode");
    mContinuous = true;
    if (before.rows()) {
      Matrix after(before.rows(), 1);
      check(GetContinuousStates(mInstance, after.data(), after.rows()),
            "get post-event states");
      if (!after.allFinite())
        fail("event produced non-finite continuous states");
      const bool changed = (after - before).cwiseAbs().maxCoeff() >
                           1e-12 * (1 + before.cwiseAbs().maxCoeff());
      if (changed && !mOptions.allowContinuousStateResets)
        fail("event requested unsupported continuous-state reinitialization");
      mLastState = after;
    }
  }
  void completeStep(const Matrix *eventInputs = nullptr) override {
    // Apply the accepted point without evaluating provisional outputs first.
    applyPoint(mLastState, mLastInputs, mTime);
    fmi3Boolean event = fmi3False, terminate = fmi3False;
    if (mNeedsCompletedIntegratorStep) {
      if (!CompletedIntegratorStep)
        fail("missing function: fmi3CompletedIntegratorStep required by "
             "capability flag");
      check(CompletedIntegratorStep(mInstance, fmi3True, &event, &terminate),
            "complete integrator step");
      if (terminate)
        fail("integrator step requested termination");
    }
    const bool due = std::isfinite(mNextTimeEvent) &&
                     std::abs(mTime - mNextTimeEvent) <=
                         1e-10 * std::max(1.0, std::abs(mTime));
    const Matrix endpointIndicators = readEventIndicators();
    const bool stateEvent =
        crossedEventIndicator(mAcceptedIndicators, endpointIndicators);
    if (event || due || eventInputs || stateEvent)
      settleEvent("accepted-step event", eventInputs);
    if (eventInputs)
      mLastInputs = *eventInputs;
    if (event || due || eventInputs || stateEvent)
      cacheAccepted();
    else
      cacheAccepted(&endpointIndicators);
  }
  void inputEvent(const Matrix &state, const Matrix &inputs,
                  double time) override {
    prepareStep(state, inputs, time);
    applyPoint(mLastState, mLastInputs, mTime);
    settleEvent("external input event");
    cacheAccepted();
  }
};

class Fmi2Backend final : public ModelExchange::Backend {
  Archive mArchive;
#ifdef _WIN32
  HMODULE mLibrary = nullptr;
#else
  void *mLibrary = nullptr;
#endif
  fmi2Component mInstance = nullptr;
  bool mContinuous = false;
  bool mValid = true;
  bool mCanGetAndSetFMUState = false;
  bool mUseTrialSnapshots = false;
  ModelExchangeOptions mOptions;
  double mTime = 0, mAcceptedTime = 0;
  Matrix mLastState, mLastInputs, mAcceptedState, mAcceptedInputs;
  Matrix mAcceptedIndicators;
  std::pair<Matrix, Matrix> mAcceptedValues;
  double mNextTimeEvent = std::numeric_limits<double>::infinity();
  std::string mLastMessage;
  std::vector<fmi2ValueReference> mInputs, mOutputs;
  std::size_t mEventIndicatorMetadata = 0;
  std::unique_ptr<fmi2FMUstate> mCheckpoint;

#define FMI2_FUNCTION(name) fmi2##name##TYPE *name = nullptr
  FMI2_FUNCTION(Instantiate);
  FMI2_FUNCTION(GetVersion);
  FMI2_FUNCTION(FreeInstance);
  FMI2_FUNCTION(Terminate);
  FMI2_FUNCTION(SetupExperiment);
  FMI2_FUNCTION(EnterInitializationMode);
  FMI2_FUNCTION(ExitInitializationMode);
  FMI2_FUNCTION(EnterEventMode);
  FMI2_FUNCTION(NewDiscreteStates);
  FMI2_FUNCTION(EnterContinuousTimeMode);
  FMI2_FUNCTION(SetTime);
  FMI2_FUNCTION(SetReal);
  FMI2_FUNCTION(GetReal);
  FMI2_FUNCTION(GetContinuousStates);
  FMI2_FUNCTION(SetContinuousStates);
  FMI2_FUNCTION(GetDerivatives);
  FMI2_FUNCTION(GetNumberOfContinuousStates);
  FMI2_FUNCTION(GetNumberOfEventIndicators);
  FMI2_FUNCTION(GetEventIndicators);
  FMI2_FUNCTION(CompletedIntegratorStep);
  FMI2_FUNCTION(GetNominalsOfContinuousStates);
  FMI2_FUNCTION(GetFMUstate);
  FMI2_FUNCTION(SetFMUstate);
  FMI2_FUNCTION(FreeFMUstate);
#undef FMI2_FUNCTION

  template <class T> T symbol(const char *name, bool required = true) {
#ifdef _WIN32
    auto address = GetProcAddress(mLibrary, name);
#else
    auto address = dlsym(mLibrary, name);
#endif
    if (!address && required)
      fail(std::string("missing function: ") + name);
    return reinterpret_cast<T>(address);
  }
  void check(fmi2Status status, const char *operation) {
    if (status > fmi2Warning) {
      mValid = false;
      fail(std::string(operation) + " failed (status " +
           std::to_string(status) + "): " + mLastMessage);
    }
  }
  static void log(fmi2ComponentEnvironment environment, fmi2String, fmi2Status,
                  fmi2String, fmi2String message, ...) noexcept {
    try {
      static_cast<Fmi2Backend *>(environment)->mLastMessage =
          message ? message : "";
    } catch (...) {
    }
  }
  static void *allocate(size_t n, size_t size) { return std::calloc(n, size); }
  static void release(void *pointer) { std::free(pointer); }
  void close() noexcept {
    if (mInstance) {
      if (mCheckpoint && *mCheckpoint && FreeFMUstate)
        FreeFMUstate(mInstance, mCheckpoint.get());
      mCheckpoint.reset();
      if (mContinuous && mValid && Terminate)
        Terminate(mInstance);
      if (FreeInstance)
        FreeInstance(mInstance);
      mInstance = nullptr;
    }
    if (mLibrary) {
#ifdef _WIN32
      FreeLibrary(mLibrary);
#else
      dlclose(mLibrary);
#endif
      mLibrary = nullptr;
    }
  }
  static std::string binaryDirectory() {
#ifdef __APPLE__
    return "darwin64";
#elif defined(_WIN32)
    return "win64";
#elif defined(__linux__)
    return "linux64";
#else
    fail("unsupported host operating system");
#endif
  }
  struct Variable {
    std::string name, causality, variability, start;
    fmi2ValueReference ref = 0;
    bool real = false;
  };
  static fmi2ValueReference reference(const std::string &text) {
    std::size_t consumed = 0;
    auto value = std::stoul(text, &consumed);
    if (consumed != text.size())
      fail("invalid FMI 2.0 value reference: " + text);
    return static_cast<fmi2ValueReference>(value);
  }
  std::pair<Matrix, Matrix> readValues() {
    Matrix derivatives(model.initialState.rows(), 1),
        outputs(mOutputs.size(), 1);
    check(GetReal(mInstance, mOutputs.data(), mOutputs.size(), outputs.data()),
          "evaluate outputs");
    if (derivatives.rows())
      check(GetDerivatives(mInstance, derivatives.data(), derivatives.rows()),
            "evaluate derivatives");
    if (!derivatives.allFinite() || !outputs.allFinite())
      fail("non-finite derivatives or outputs during model recovery");
    return {derivatives, outputs};
  }
  Matrix readEventIndicators() {
    Matrix indicators(mEventIndicatorMetadata, 1);
    if (indicators.rows())
      check(GetEventIndicators(mInstance, indicators.data(), indicators.rows()),
            "evaluate event indicators");
    return indicators;
  }
  void applyPoint(const Matrix &state, const Matrix &inputs, double time) {
    check(SetTime(mInstance, time), "set evaluation time");
    if (state.rows())
      check(SetContinuousStates(mInstance, state.data(), state.rows()),
            "set trial states");
    if (mInputs.size())
      check(SetReal(mInstance, mInputs.data(), mInputs.size(), inputs.data()),
            "set trial inputs");
  }
  static bool crossed(const Matrix &before, const Matrix &after) {
    for (Eigen::Index k = 0; k < before.rows(); ++k)
      if ((before(k, 0) > 0) != (after(k, 0) > 0))
        return true;
    return false;
  }
  void updateEventInfo(const fmi2EventInfo &info) {
    if (info.nextEventTimeDefined) {
      if (!std::isfinite(info.nextEventTime) || info.nextEventTime <= mTime)
        fail("FMI 2.0 next event time must advance time");
      mNextTimeEvent = info.nextEventTime;
    } else {
      mNextTimeEvent = std::numeric_limits<double>::infinity();
    }
  }
  void cacheAccepted(const Matrix *indicators = nullptr) {
    mAcceptedTime = mTime;
    mAcceptedState = mLastState;
    mAcceptedInputs = mLastInputs;
    mAcceptedValues = readValues();
    mAcceptedIndicators = indicators ? *indicators : readEventIndicators();
  }
  void settleEvent(const Matrix *eventInputs = nullptr) {
    if (!EnterEventMode || !NewDiscreteStates || !EnterContinuousTimeMode)
      fail("FMI 2.0 event functions are unavailable");
    Matrix before = mLastState;
    check(EnterEventMode(mInstance), "enter FMI 2.0 event mode");
    mContinuous = false;
    if (eventInputs)
      check(SetReal(mInstance, mInputs.data(), mInputs.size(),
                    eventInputs->data()),
            "apply event inputs");
    fmi2EventInfo info{};
    bool settled = false;
    for (int k = 0; k < 100; ++k) {
      check(NewDiscreteStates(mInstance, &info), "FMI 2.0 event iteration");
      if (info.terminateSimulation)
        fail("initialization/event requested termination");
      if (!info.newDiscreteStatesNeeded) {
        updateEventInfo(info);
        settled = true;
        break;
      }
    }
    if (!settled)
      fail("FMI 2.0 event iteration did not converge");
    check(EnterContinuousTimeMode(mInstance), "leave FMI 2.0 event mode");
    mContinuous = true;
    if (before.rows()) {
      Matrix after(before.rows(), 1);
      check(GetContinuousStates(mInstance, after.data(), after.rows()),
            "get post-event states");
      if ((after - before).cwiseAbs().maxCoeff() >
              1e-12 * (1 + before.cwiseAbs().maxCoeff()) &&
          !mOptions.allowContinuousStateResets)
        fail("FMI 2.0 event requested unsupported state reinitialization");
      mLastState = after;
    }
  }

public:
  Fmi2Backend(const std::string &path,
              const std::vector<ScalarInput> &configuration,
              ModelExchangeOptions options)
      : mOptions(options) {
    try {
      initialize(path, configuration);
    } catch (...) {
      close();
      throw;
    }
  }
  ~Fmi2Backend() override { close(); }
  void initialize(const std::string &path,
                  const std::vector<ScalarInput> &configuration) {
    const bool automaticInputs = configuration.empty();
    if (configuration.size() > 64)
      fail("input count exceeds maximum 64");
    if (!automaticInputs) {
      model.initialInputs = Matrix::Zero(configuration.size(), 1);
      for (std::size_t k = 0; k < configuration.size(); ++k) {
        if (configuration[k].name.empty() ||
            !std::isfinite(configuration[k].initialValue))
          fail("invalid FMI 2.0 input mapping");
        model.inputNames.push_back(configuration[k].name);
        model.initialInputs(k, 0) = configuration[k].initialValue;
      }
    }
    mArchive.extract(path);
    const auto xmlPath = mArchive.directory / "modelDescription.xml";
    std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> document(
        xmlReadFile(xmlPath.string().c_str(), nullptr, XML_PARSE_NONET),
        &xmlFreeDoc);
    if (!document || document->intSubset || document->extSubset)
      fail("invalid modelDescription.xml (DTDs are unsupported)");
    auto root = xmlDocGetRootElement(document.get());
    if (!root || !tag(root, "fmiModelDescription") ||
        property(root, "fmiVersion") != "2.0")
      fail("requires an FMI 2.0 model description");
    auto me = child(root, "ModelExchange");
    const auto identifier = property(me, "modelIdentifier");
    const auto guid = property(root, "guid");
    if (identifier.empty() || guid.empty())
      fail("invalid FMI 2.0 modelIdentifier or guid");
    model.modelName = property(root, "modelName", identifier);
    model.modelIdentifier = identifier;
    model.instantiationToken = guid;
    mCanGetAndSetFMUState = flag(me, "canGetAndSetFMUstate");
    model.canGetAndSetFMUState = mCanGetAndSetFMUState;
    model.canSerializeFMUState = flag(me, "canSerializeFMUstate");
    model.canBeInstantiatedOnlyOncePerProcess =
        flag(me, "canBeInstantiatedOnlyOncePerProcess");
    std::map<fmi2ValueReference, Variable> variables;
    std::set<std::string> names;
    std::vector<std::string> variableNames;
    std::vector<bool> foundInputs(configuration.size(), false);
    std::vector<double> discoveredStarts;
    auto modelVariables = child(root, "ModelVariables");
    for (auto v = modelVariables->children; v; v = v->next) {
      if (!tag(v, "ScalarVariable"))
        continue;
      variableNames.push_back(property(v, "name"));
      xmlNode *real = nullptr;
      for (auto type = v->children; type; type = type->next)
        if (tag(type, "Real")) {
          real = type;
          break;
        }
      if (!real)
        continue;
      Variable variable{property(v, "name"),
                        property(v, "causality", "local"),
                        property(v, "variability", "continuous"),
                        property(real, "start", "0"),
                        reference(property(v, "valueReference")),
                        true};
      if (variable.variability == "discrete")
        model.hasDiscreteVariables = true;
      if (!names.insert(variable.name).second)
        fail("duplicate FMI variable name");
      // FMI permits aliases: several variables may intentionally share one
      // value reference. Keep the first descriptor for derivative metadata;
      // names and the input/output reference lists retain every variable.
      variables.emplace(variable.ref, variable);
      if (variable.causality == "input") {
        if (automaticInputs) {
          std::size_t consumed = 0;
          double start = std::stod(variable.start, &consumed);
          if (consumed != variable.start.size() || !std::isfinite(start))
            fail("invalid FMI 2.0 input start for " + variable.name);
          model.inputNames.push_back(variable.name);
          mInputs.push_back(variable.ref);
          discoveredStarts.push_back(start);
        } else {
          auto found = std::find(model.inputNames.begin(),
                                 model.inputNames.end(), variable.name);
          if (found == model.inputNames.end())
            continue;
          foundInputs[found - model.inputNames.begin()] = true;
          mInputs.resize(configuration.size());
          mInputs[found - model.inputNames.begin()] = variable.ref;
        }
      } else if (variable.causality == "output") {
        model.outputNames.push_back(variable.name);
        mOutputs.push_back(variable.ref);
      }
    }
    if (!automaticInputs) {
      if (std::find(foundInputs.begin(), foundInputs.end(), false) !=
          foundInputs.end())
        fail("mapped input does not exist");
    } else {
      model.initialInputs = Matrix::Zero(discoveredStarts.size(), 1);
      for (std::size_t k = 0; k < discoveredStarts.size(); ++k)
        model.initialInputs(k, 0) = discoveredStarts[k];
    }
    auto structure = child(root, "ModelStructure");
    for (auto entry = structure->children; entry; entry = entry->next) {
      if (tag(entry, "Derivatives")) {
        for (auto unknown = entry->children; unknown; unknown = unknown->next)
          if (tag(unknown, "Unknown")) {
            const auto derivative = property(unknown, "derivative");
            std::size_t index = 0, consumed = 0;
            if (!derivative.empty())
              index = std::stoul(derivative, &consumed);
            if (consumed == derivative.size() && index > 0 &&
                index <= variableNames.size())
              model.stateNames.push_back(variableNames[index - 1]);
            else
              model.stateNames.push_back(
                  "state" + std::to_string(model.stateNames.size()));
          }
      }
      if (tag(entry, "EventIndicators")) {
        for (auto unknown = entry->children; unknown; unknown = unknown->next)
          if (tag(unknown, "Unknown"))
            ++mEventIndicatorMetadata;
        model.hasEventIndicators = mEventIndicatorMetadata != 0;
      }
    }
    if (model.hasEventIndicators && !mOptions.endpointStateEvents)
      fail("state events are unsupported unless endpoint handling is "
           "explicitly enabled");
    auto binary =
        mArchive.directory / "binaries" / binaryDirectory() / identifier;
#ifdef __APPLE__
    binary += ".dylib";
#elif defined(_WIN32)
    binary += ".dll";
#else
    binary += ".so";
#endif
    if (!fs::is_regular_file(binary))
      fail("missing FMI 2.0 host binary: " +
           binary.lexically_relative(mArchive.directory).string());
#ifdef _WIN32
    mLibrary = LoadLibraryW(binary.wstring().c_str());
#else
    mLibrary = dlopen(binary.string().c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
    if (!mLibrary)
      fail("cannot load FMI 2.0 FMU binary");
#define LOAD2(name) name = symbol<decltype(name)>("fmi2" #name)
    LOAD2(Instantiate);
    LOAD2(GetVersion);
    LOAD2(FreeInstance);
    LOAD2(Terminate);
    LOAD2(SetupExperiment);
    LOAD2(EnterInitializationMode);
    LOAD2(ExitInitializationMode);
    LOAD2(EnterEventMode);
    LOAD2(NewDiscreteStates);
    LOAD2(EnterContinuousTimeMode);
    LOAD2(SetTime);
    LOAD2(SetReal);
    LOAD2(GetReal);
    LOAD2(GetContinuousStates);
    LOAD2(SetContinuousStates);
    LOAD2(GetDerivatives);
    LOAD2(GetNumberOfContinuousStates);
    LOAD2(GetNumberOfEventIndicators);
    LOAD2(GetEventIndicators);
    LOAD2(CompletedIntegratorStep);
    LOAD2(GetNominalsOfContinuousStates);
    GetFMUstate = symbol<decltype(GetFMUstate)>("fmi2GetFMUstate", false);
    SetFMUstate = symbol<decltype(SetFMUstate)>("fmi2SetFMUstate", false);
    FreeFMUstate = symbol<decltype(FreeFMUstate)>("fmi2FreeFMUstate", false);
#undef LOAD2
    if (!GetVersion() || std::string(GetVersion()) != "2.0")
      fail("binary does not implement FMI 2.0");
    if (model.canBeInstantiatedOnlyOncePerProcess) {
      static std::mutex mutex;
      static std::set<std::string> instantiated;
      std::lock_guard<std::mutex> guard(mutex);
      if (!instantiated.insert(identifier + ":" + guid).second)
        fail("FMU permits only one instantiation per process; load a separate "
             "export or reuse the imported model");
    }
    fmi2CallbackFunctions callbacks{log, allocate, release, nullptr, this};
    mInstance = Instantiate(model.modelName.c_str(), 0, guid.c_str(), nullptr,
                            &callbacks, fmi2False, fmi2False);
    if (!mInstance)
      fail("FMI 2.0 instantiation failed: " + mLastMessage);
    check(SetupExperiment(mInstance, fmi2False, 0, 0, fmi2False, 0),
          "setup FMI 2.0 experiment");
    check(EnterInitializationMode(mInstance), "enter FMI 2.0 initialization");
    check(SetReal(mInstance, mInputs.data(), mInputs.size(),
                  model.initialInputs.data()),
          "set initial inputs");
    check(ExitInitializationMode(mInstance), "exit FMI 2.0 initialization");
    check(EnterContinuousTimeMode(mInstance), "enter FMI 2.0 continuous mode");
    mContinuous = true;
    size_t states = 0, indicators = 0;
    check(GetNumberOfContinuousStates(mInstance, &states),
          "get FMI 2.0 state count");
    check(GetNumberOfEventIndicators(mInstance, &indicators),
          "get FMI 2.0 event indicator count");
    if (indicators != mEventIndicatorMetadata)
      fail("FMI 2.0 event indicator count is inconsistent with XML");
    if (states > 512)
      fail("FMI 2.0 state count is excessive");
    model.stateNames.resize(states);
    for (std::size_t k = 0; k < states; ++k)
      if (model.stateNames[k].empty())
        model.stateNames[k] = "state" + std::to_string(k);
    model.initialState.resize(states, 1);
    if (states)
      check(GetContinuousStates(mInstance, model.initialState.data(), states),
            "get FMI 2.0 initial states");
    model.stateNominals = Matrix::Ones(states, 1);
    if (states)
      check(GetNominalsOfContinuousStates(mInstance, model.stateNominals.data(),
                                          states),
            "get FMI 2.0 state nominals");
    if (!model.initialState.allFinite() || !model.stateNominals.allFinite() ||
        (model.stateNominals.array() <= 0).any())
      fail("invalid FMI 2.0 initial state data");
    mLastState = mAcceptedState = model.initialState;
    mLastInputs = mAcceptedInputs = model.initialInputs;
    mAcceptedValues = readValues();
    mAcceptedIndicators = readEventIndicators();
    mUseTrialSnapshots = mOptions.snapshotTrialEvaluations &&
                         mCanGetAndSetFMUState && GetFMUstate && SetFMUstate &&
                         FreeFMUstate &&
                         (model.hasDiscreteVariables || indicators != 0);
  }
  std::pair<Matrix, Matrix> evaluate(const Matrix &state, const Matrix &inputs,
                                     double time) override {
    prepareStep(state, inputs, time);
    if (std::abs(mTime - mAcceptedTime) <= 1e-12 && state == mAcceptedState &&
        inputs == mAcceptedInputs)
      return mAcceptedValues;
    fmi2FMUstate saved = nullptr;
    if (mUseTrialSnapshots)
      check(GetFMUstate(mInstance, &saved), "snapshot FMI 2.0 trial");
    try {
      applyPoint(state, inputs, mTime);
      auto result = readValues();
      if (saved) {
        check(SetFMUstate(mInstance, saved), "restore FMI 2.0 trial");
        check(FreeFMUstate(mInstance, &saved), "free FMI 2.0 trial");
      }
      return result;
    } catch (...) {
      if (saved) {
        if (mValid)
          SetFMUstate(mInstance, saved);
        FreeFMUstate(mInstance, &saved);
      }
      throw;
    }
  }
  double nextTimeEvent() const override { return mNextTimeEvent; }
  const std::pair<Matrix, Matrix> &acceptedValues() const override {
    return mAcceptedValues;
  }
  const Matrix &acceptedState() const override { return mAcceptedState; }
  bool supportsRollback() const override {
    return mCanGetAndSetFMUState && GetFMUstate && SetFMUstate && FreeFMUstate;
  }
  void saveCheckpoint() override {
    if (!supportsRollback())
      fail("multirate integration requires FMI 2.0 FMU state support");
    if (mCheckpoint)
      fail("an FMI checkpoint is already active");
    auto state = std::make_unique<fmi2FMUstate>(nullptr);
    check(GetFMUstate(mInstance, state.get()), "save FMI 2.0 checkpoint");
    mCheckpoint = std::move(state);
  }
  void restoreCheckpoint() override {
    if (!mCheckpoint || !*mCheckpoint)
      fail("no active FMI 2.0 checkpoint");
    check(SetFMUstate(mInstance, *mCheckpoint), "restore FMI 2.0 checkpoint");
  }
  void discardCheckpoint() override {
    if (mCheckpoint && *mCheckpoint)
      check(FreeFMUstate(mInstance, mCheckpoint.get()),
            "free FMI 2.0 checkpoint");
    mCheckpoint.reset();
  }
  void prepareStep(const Matrix &state, const Matrix &inputs,
                   double time) override {
    if (state.rows() != model.initialState.rows() || state.cols() != 1 ||
        inputs.rows() != model.initialInputs.rows() || inputs.cols() != 1 ||
        !state.allFinite() || !inputs.allFinite() || !std::isfinite(time))
      fail("invalid state/input/time in FMI 2.0 evaluation");
    const double tolerance = 1e-10 * std::max(1.0, std::abs(time));
    if (time > mNextTimeEvent + tolerance)
      fail("evaluation crossed a scheduled FMI 2.0 time event");
    mTime =
        std::abs(time - mNextTimeEvent) <= tolerance ? mNextTimeEvent : time;
    mLastState = state;
    mLastInputs = inputs;
  }
  void completeStep(const Matrix *eventInputs = nullptr) override {
    applyPoint(mLastState, mLastInputs, mTime);
    fmi2Boolean event = fmi2False, terminate = fmi2False;
    check(CompletedIntegratorStep(mInstance, fmi2True, &event, &terminate),
          "complete FMI 2.0 integrator step");
    if (terminate)
      fail("FMI 2.0 integrator step requested termination");
    const Matrix endpoint = readEventIndicators();
    const bool due = std::isfinite(mNextTimeEvent) &&
                     std::abs(mTime - mNextTimeEvent) <=
                         1e-10 * std::max(1.0, std::abs(mTime));
    if (event || due || eventInputs || crossed(mAcceptedIndicators, endpoint))
      settleEvent(eventInputs);
    if (eventInputs)
      mLastInputs = *eventInputs;
    if (event || due || eventInputs || crossed(mAcceptedIndicators, endpoint))
      cacheAccepted();
    else
      cacheAccepted(&endpoint);
  }
  void inputEvent(const Matrix &state, const Matrix &inputs,
                  double time) override {
    prepareStep(state, inputs, time);
    applyPoint(mLastState, mLastInputs, mTime);
    settleEvent();
    cacheAccepted();
  }
};

std::string detectFmiVersion(const std::string &path) {
  Archive archive;
  archive.extract(path);
  const auto xmlPath = archive.directory / "modelDescription.xml";
  std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> document(
      xmlReadFile(xmlPath.string().c_str(), nullptr, XML_PARSE_NONET),
      &xmlFreeDoc);
  if (!document)
    fail("invalid modelDescription.xml");
  auto root = xmlDocGetRootElement(document.get());
  if (!root || !tag(root, "fmiModelDescription"))
    fail("invalid fmiModelDescription root element");
  return property(root, "fmiVersion");
}

ModelExchange::ModelExchange(const std::string &path,
                             const std::vector<ScalarInput> &inputs,
                             ModelExchangeOptions options)
    : mImpl(detectFmiVersion(path) == "2.0"
                ? std::unique_ptr<Backend>(
                      std::make_unique<Fmi2Backend>(path, inputs, options))
                : std::unique_ptr<Backend>(
                      std::make_unique<Fmi3Backend>(path, inputs, options))) {}
ModelExchange::ModelExchange(const std::string &path,
                             ModelExchangeOptions options)
    : ModelExchange(path, std::vector<ScalarInput>{}, options) {}
ModelExchange::~ModelExchange() = default;
const ModelDescription &ModelExchange::description() const {
  return mImpl->model;
}
std::pair<Matrix, Matrix> ModelExchange::evaluate(const Matrix &state,
                                                  const Matrix &inputs,
                                                  double time) {
  return mImpl->evaluate(state, inputs, time);
}
const std::pair<Matrix, Matrix> &ModelExchange::acceptedValues() const {
  return mImpl->acceptedValues();
}
const Matrix &ModelExchange::acceptedState() const {
  return mImpl->acceptedState();
}
bool ModelExchange::supportsRollback() const {
  return mImpl->supportsRollback();
}
void ModelExchange::saveCheckpoint() { mImpl->saveCheckpoint(); }
void ModelExchange::restoreCheckpoint() { mImpl->restoreCheckpoint(); }
void ModelExchange::discardCheckpoint() { mImpl->discardCheckpoint(); }
double ModelExchange::nextTimeEvent() const { return mImpl->nextTimeEvent(); }
void ModelExchange::prepareStep(const Matrix &state, const Matrix &inputs,
                                double time) {
  mImpl->prepareStep(state, inputs, time);
}
void ModelExchange::completeStep() { mImpl->completeStep(); }
void ModelExchange::completeStep(const Matrix &eventInputs) {
  mImpl->completeStep(&eventInputs);
}
void ModelExchange::inputEvent(const Matrix &state, const Matrix &inputs,
                               double time) {
  mImpl->inputEvent(state, inputs, time);
}
} // namespace DPsim::FMI
