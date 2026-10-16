// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/ModelExchange.h>
#include <fmi3FunctionTypes.h>
#include <libxml/parser.h>
#include <miniz.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <set>
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
  std::string name, type, causality, derivative, variability;
  fmi3ValueReference ref;
  bool scalar;
};

} // namespace

class ModelExchange::Impl {
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
      static_cast<Impl *>(environment)->mLastMessage = message ? message : "";
    } catch (...) {
      // Never allow a C++ exception across the FMI C callback boundary.
    }
  }
  void close() noexcept {
    if (mInstance) {
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
  ModelDescription model;
  Impl(const std::string &path, const std::vector<ScalarInput> &configuration,
       ModelExchangeOptions options)
      : mOptions(options) {
    try {
      initialize(path, configuration);
    } catch (...) {
      close();
      throw;
    }
  }
  ~Impl() { close(); }
  Impl(const Impl &) = delete;
  Impl &operator=(const Impl &) = delete;

  void initialize(const std::string &path,
                  const std::vector<ScalarInput> &configuration) {
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
                        reference(property(v, "valueReference")),
                        true};
      for (auto dimension = v->children; dimension; dimension = dimension->next)
        if (tag(dimension, "Dimension"))
          variable.scalar = false;
      if (!names.insert(variable.name).second ||
          !variables.emplace(variable.ref, variable).second)
        fail("duplicate variable name or value reference");
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
    auto structure = child(root, "ModelStructure");
    for (auto entry = structure->children; entry; entry = entry->next) {
      if (tag(entry, "EventIndicator")) {
        ++mEventIndicatorMetadata;
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
      fmi3Boolean again, terminate, nominalsChanged, statesChanged, timeEvent;
      double nextTime;
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
    mUseTrialSnapshots =
        GetFMUState && (model.hasDiscreteVariables ||
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
  const std::pair<Matrix, Matrix> &acceptedValues() const {
    return mAcceptedValues;
  }
  const Matrix &acceptedState() const { return mAcceptedState; }

  std::pair<Matrix, Matrix> evaluate(const Matrix &state, const Matrix &voltage,
                                     double time = 0) {
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
  void prepareStep(const Matrix &state, const Matrix &voltage, double time) {
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
    if (defined && (!std::isfinite(next) || next <= mTime))
      fail("FMI next time event must advance time");
    mNextTimeEvent = defined ? next : std::numeric_limits<double>::infinity();
    if (GetFMUState && (defined || model.hasDiscreteVariables))
      mUseTrialSnapshots = true;
  }
  double nextTimeEvent() const { return mNextTimeEvent; }
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
      fmi3Boolean again, terminate, nominals, values, timeEvent;
      double next;
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
  void completeStep(const Matrix *eventInputs = nullptr) {
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
  void inputEvent(const Matrix &state, const Matrix &inputs, double time) {
    prepareStep(state, inputs, time);
    applyPoint(mLastState, mLastInputs, mTime);
    settleEvent("external input event");
    cacheAccepted();
  }
};

ModelExchange::ModelExchange(const std::string &path,
                             const std::vector<ScalarInput> &inputs,
                             ModelExchangeOptions options)
    : mImpl(std::make_unique<Impl>(path, inputs, options)) {}
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
