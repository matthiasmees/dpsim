// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/CoSimulation.h>
#include <fmi2FunctionTypes.h>
#include <fmi3FunctionTypes.h>
#include <libxml/parser.h>
#include <miniz.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
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
  throw std::runtime_error("FMI Co-Simulation import: " + message);
}

class Archive {
public:
  fs::path directory;
  Archive() {
    std::random_device random;
    for (int attempt = 0; attempt < 100; ++attempt) {
      const auto path = fs::temp_directory_path() /
                        ("dpsim-fmi-cs-" + std::to_string(random()) + "-" +
                         std::to_string(random()));
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
      const std::string name = stat.m_filename;
      if (name.size() == sizeof(stat.m_filename) - 1)
        fail("archive entry name is too long");
      const fs::path relative(name);
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

bool tag(xmlNode *node, const char *name) {
  return node->type == XML_ELEMENT_NODE &&
         xmlStrEqual(node->name, BAD_CAST name);
}
std::string property(xmlNode *node, const char *name,
                     const std::string &fallback = "") {
  xmlChar *value = xmlGetProp(node, BAD_CAST name);
  if (!value)
    return fallback;
  std::string result(reinterpret_cast<const char *>(value));
  xmlFree(value);
  return result;
}
bool flag(xmlNode *node, const char *name, bool fallback = false) {
  const auto value = property(node, name, fallback ? "true" : "false");
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
std::uint32_t reference(const std::string &text) {
  std::size_t consumed = 0;
  const auto value = std::stoull(text, &consumed);
  if (consumed != text.size() || value > UINT32_MAX)
    fail("invalid value reference: " + text);
  return static_cast<std::uint32_t>(value);
}
std::string binaryDirectory() {
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
fs::path libraryPath(const fs::path &root, const std::string &identifier) {
  auto result = root / "binaries" / binaryDirectory() / identifier;
#ifdef __APPLE__
  result += ".dylib";
#elif defined(_WIN32)
  result += ".dll";
#else
  result += ".so";
#endif
  return result;
}
fs::path fmi2LibraryPath(const fs::path &root, const std::string &identifier) {
#ifdef __APPLE__
  const auto directory = "darwin64";
  const auto extension = ".dylib";
#elif defined(_WIN32)
  const auto directory = "win64";
  const auto extension = ".dll";
#elif defined(__linux__)
  const auto directory = "linux64";
  const auto extension = ".so";
#else
  fail("unsupported FMI 2.0 host operating system");
#endif
  return root / "binaries" / directory / (identifier + extension);
}
void enforceSingleInstantiation(bool once, const std::string &key) {
  if (!once)
    return;
  static std::mutex mutex;
  static std::set<std::string> instantiated;
  std::lock_guard<std::mutex> guard(mutex);
  if (!instantiated.insert(key).second)
    fail("FMU permits only one instantiation per process");
}
} // namespace

class CoSimulation::Backend {
public:
  ModelDescription model;
  Matrix inputs, outputs;
  std::vector<std::string> inputNames;
  std::vector<std::uint32_t> inputReferences, outputReferences;
  virtual ~Backend() = default;
  virtual void setInputs(const Matrix &) = 0;
  virtual void doStep(double, double) = 0;
  virtual double time() const = 0;
};

namespace {
class Fmi3CoSimulation final : public CoSimulation::Backend {
  Archive mArchive;
#ifdef _WIN32
  HMODULE mLibrary = nullptr;
#else
  void *mLibrary = nullptr;
#endif
  fmi3Instance mInstance = nullptr;
  bool mValid = true;
  bool mTerminated = false;
  bool mEventModeUsed = false;
  bool mLoggingOn = false;
  double mTime = 0;
  std::string mLastMessage;

#define FMI3_FUNCTION(name) fmi3##name##TYPE *name = nullptr
  FMI3_FUNCTION(InstantiateCoSimulation);
  FMI3_FUNCTION(GetVersion);
  FMI3_FUNCTION(FreeInstance);
  FMI3_FUNCTION(Terminate);
  FMI3_FUNCTION(EnterInitializationMode);
  FMI3_FUNCTION(ExitInitializationMode);
  FMI3_FUNCTION(SetFloat64);
  FMI3_FUNCTION(GetFloat64);
  FMI3_FUNCTION(DoStep);
  FMI3_FUNCTION(EnterStepMode);
  FMI3_FUNCTION(EnterEventMode);
  FMI3_FUNCTION(UpdateDiscreteStates);
#undef FMI3_FUNCTION

  template <class T> T symbol(const char *name, bool required = true) {
#ifdef _WIN32
    auto address = GetProcAddress(mLibrary, name);
#else
    auto address = dlsym(mLibrary, name);
#endif
    if (!address && required)
      fail(std::string("missing FMI 3 function: ") + name);
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
      static_cast<Fmi3CoSimulation *>(environment)->mLastMessage =
          message ? message : "";
    } catch (...) {
    }
  }
  void close() noexcept {
    if (mInstance) {
      if (!mTerminated && mValid && Terminate)
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
  void settleEventMode(const char *operation, bool enterEventMode) {
    if (!mEventModeUsed || !EnterEventMode || !UpdateDiscreteStates ||
        !EnterStepMode)
      fail(std::string(operation) + " requested event mode but required FMI 3 "
                                    "event functions are absent");
    if (enterEventMode)
      check(EnterEventMode(mInstance), "enter FMI 3 event mode");
    bool settled = false;
    for (int iteration = 0; iteration < 100; ++iteration) {
      fmi3Boolean again = fmi3False, terminate = fmi3False;
      fmi3Boolean nominalsChanged = fmi3False, statesChanged = fmi3False;
      fmi3Boolean nextTimeDefined = fmi3False;
      fmi3Float64 nextTime = 0;
      check(UpdateDiscreteStates(mInstance, &again, &terminate,
                                 &nominalsChanged, &statesChanged,
                                 &nextTimeDefined, &nextTime),
            "update FMI 3 discrete states");
      if (terminate)
        fail("FMU requested termination during event iteration");
      if (!again) {
        settled = true;
        break;
      }
    }
    if (!settled)
      fail("FMI 3 event iteration did not converge");
    check(EnterStepMode(mInstance), "enter FMI 3 step mode");
  }

public:
  Fmi3CoSimulation(const std::string &path,
                   const std::vector<ScalarInput> &configuration,
                   CoSimulationOptions options) {
    mLoggingOn = options.loggingOn;
    try {
      initialize(path, configuration);
    } catch (...) {
      close();
      throw;
    }
  }
  ~Fmi3CoSimulation() override { close(); }

  void initialize(const std::string &path,
                  const std::vector<ScalarInput> &configuration) {
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
    auto cs = child(root, "CoSimulation");
    mEventModeUsed = flag(cs, "eventModeUsed");
    model.canBeInstantiatedOnlyOncePerProcess =
        flag(cs, "canBeInstantiatedOnlyOncePerProcess");
    if (flag(cs, "needsExecutionTool"))
      fail("FMUs requiring an external execution tool are unsupported");
    const auto identifier = property(cs, "modelIdentifier");
    const auto token = property(root, "instantiationToken");
    if (identifier.empty() || token.empty())
      fail("invalid FMI 3 modelIdentifier or instantiationToken");
    model.modelName = property(root, "modelName");
    model.modelIdentifier = identifier;
    model.instantiationToken = token;
    if (model.modelName.empty())
      fail("modelName is missing");
    const auto count = configuration.size();
    if (count > 64)
      fail("input count exceeds maximum 64");
    const bool automaticInputs = configuration.empty();
    if (!automaticInputs) {
      std::set<std::string> unique;
      model.initialInputs = Matrix::Zero(count, 1);
      for (std::size_t i = 0; i < count; ++i) {
        if (configuration[i].name.empty() ||
            !unique.insert(configuration[i].name).second ||
            !std::isfinite(configuration[i].initialValue))
          fail("invalid or duplicate input mapping");
        inputNames.push_back(configuration[i].name);
        model.initialInputs(i, 0) = configuration[i].initialValue;
      }
    }
    std::map<std::string, std::size_t> mapped;
    for (std::size_t i = 0; i < inputNames.size(); ++i)
      mapped.emplace(inputNames[i], i);
    std::vector<bool> found(inputNames.size(), false);
    std::vector<double> starts;
    std::set<std::string> variableNames;
    auto variables = child(root, "ModelVariables");
    for (auto variable = variables->children; variable;
         variable = variable->next) {
      if (variable->type != XML_ELEMENT_NODE)
        continue;
      const auto name = property(variable, "name");
      const auto causality = property(variable, "causality", "local");
      if (!variableNames.insert(name).second)
        fail("duplicate FMI variable name");
      if (causality != "input" && causality != "output")
        continue;
      bool scalar = true;
      for (auto dimension = variable->children; dimension;
           dimension = dimension->next)
        if (tag(dimension, "Dimension"))
          scalar = false;
      if (!tag(variable, "Float64") || !scalar)
        fail("requires scalar Float64 Co-Simulation inputs and outputs");
      const auto ref = reference(property(variable, "valueReference"));
      if (causality == "input") {
        if (automaticInputs) {
          inputNames.push_back(name);
          inputReferences.push_back(ref);
          std::size_t consumed = 0;
          const auto text = property(variable, "start", "0");
          const double start = std::stod(text, &consumed);
          if (consumed != text.size() || !std::isfinite(start))
            fail("invalid Float64 input start for " + name);
          starts.push_back(start);
        } else {
          const auto foundName = mapped.find(name);
          if (foundName != mapped.end()) {
            const auto index = foundName->second;
            if (found[index])
              fail("duplicate mapped input variable " + name);
            if (inputReferences.empty())
              inputReferences.resize(inputNames.size());
            inputReferences[index] = ref;
            found[index] = true;
          }
        }
      } else {
        model.outputNames.push_back(name);
        outputReferences.push_back(ref);
      }
    }
    if (inputReferences.size() != inputNames.size())
      fail("all mapped input variables must exist exactly once");
    if (!automaticInputs &&
        std::find(found.begin(), found.end(), false) != found.end())
      fail("mapped input does not exist");
    if (automaticInputs) {
      model.initialInputs.resize(starts.size(), 1);
      for (std::size_t i = 0; i < starts.size(); ++i)
        model.initialInputs(i, 0) = starts[i];
    }
    model.inputNames = inputNames;
    model.initialState.resize(0, 1);
    model.stateNominals.resize(0, 1);
    if (model.outputNames.empty())
      fail("FMU has no scalar Float64 outputs");
    auto binary = libraryPath(mArchive.directory, identifier);
    if (!fs::is_regular_file(binary))
      fail("missing host binary: " +
           binary.lexically_relative(mArchive.directory).string());
#ifdef _WIN32
    mLibrary = LoadLibraryW(binary.wstring().c_str());
#else
    mLibrary = dlopen(binary.string().c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
    if (!mLibrary)
      fail("cannot load FMU binary");
#define LOAD3(name) name = symbol<decltype(name)>("fmi3" #name)
    LOAD3(InstantiateCoSimulation);
    LOAD3(GetVersion);
    LOAD3(FreeInstance);
    LOAD3(Terminate);
    LOAD3(EnterInitializationMode);
    LOAD3(ExitInitializationMode);
    LOAD3(SetFloat64);
    LOAD3(GetFloat64);
    LOAD3(DoStep);
#undef LOAD3
    EnterStepMode = symbol<decltype(EnterStepMode)>("fmi3EnterStepMode", false);
    EnterEventMode =
        symbol<decltype(EnterEventMode)>("fmi3EnterEventMode", false);
    UpdateDiscreteStates = symbol<decltype(UpdateDiscreteStates)>(
        "fmi3UpdateDiscreteStates", false);
    if (!GetVersion() || std::string(GetVersion()) != "3.0")
      fail("binary does not implement FMI 3.0");
    enforceSingleInstantiation(model.canBeInstantiatedOnlyOncePerProcess,
                               "fmi3:" + identifier + ":" + token);
    const auto resources = mArchive.directory / "resources";
    const auto resourcePath =
        resources.string() + fs::path::preferred_separator;
    mInstance = InstantiateCoSimulation(
        model.modelName.c_str(), token.c_str(),
        fs::is_directory(resources) ? resourcePath.c_str() : nullptr, fmi3False,
        mLoggingOn ? fmi3True : fmi3False,
        mEventModeUsed ? fmi3True : fmi3False, fmi3False, nullptr, 0, this, log,
        nullptr);
    if (!mInstance)
      fail("instantiation failed: " + mLastMessage);
    check(EnterInitializationMode(mInstance, fmi3False, 0, 0, fmi3False, 0),
          "enter initialization mode");
    setInputs(model.initialInputs);
    check(ExitInitializationMode(mInstance), "exit initialization mode");
    if (mEventModeUsed)
      settleEventMode("initialization", false);
    inputs = model.initialInputs;
    outputs.resize(model.outputNames.size(), 1);
    if (!outputReferences.empty())
      check(GetFloat64(mInstance, outputReferences.data(),
                       outputReferences.size(), outputs.data(), outputs.size()),
            "read initial outputs");
    if (!outputs.allFinite())
      fail("initial outputs contain non-finite values");
  }
  void setInputs(const Matrix &values) override {
    if (values.rows() != static_cast<Eigen::Index>(inputReferences.size()) ||
        values.cols() != 1 || !values.allFinite())
      fail("input vector has invalid dimensions or values");
    if (!inputReferences.empty())
      check(SetFloat64(mInstance, inputReferences.data(),
                       inputReferences.size(), values.data(), values.size()),
            "set Co-Simulation inputs");
    inputs = values;
  }
  void doStep(double currentTime, double stepSize) override {
    if (!std::isfinite(currentTime) || !std::isfinite(stepSize) ||
        stepSize <= 0 ||
        std::abs(currentTime - mTime) > 1e-10 * std::max(1.0, std::abs(mTime)))
      fail("communication point must match current FMU time and step size be "
           "positive");
    fmi3Boolean eventHandlingNeeded = fmi3False;
    fmi3Boolean terminateSimulation = fmi3False;
    fmi3Boolean earlyReturn = fmi3False;
    fmi3Float64 lastSuccessfulTime = currentTime;
    check(DoStep(mInstance, currentTime, stepSize, fmi3True,
                 &eventHandlingNeeded, &terminateSimulation, &earlyReturn,
                 &lastSuccessfulTime),
          "doStep");
    if (terminateSimulation)
      fail("FMU requested termination during doStep");
    if (earlyReturn ||
        std::abs(lastSuccessfulTime - (currentTime + stepSize)) >
            1e-10 * std::max(1.0, std::abs(currentTime + stepSize)))
      fail("early return from doStep is unsupported");
    mTime = currentTime + stepSize;
    if (eventHandlingNeeded)
      settleEventMode("doStep", true);
    check(GetFloat64(mInstance, outputReferences.data(),
                     outputReferences.size(), outputs.data(), outputs.size()),
          "read outputs after doStep");
    if (!outputs.allFinite())
      fail("output vector contains non-finite values");
  }
  double time() const override { return mTime; }
};

class Fmi2CoSimulation final : public CoSimulation::Backend {
  Archive mArchive;
#ifdef _WIN32
  HMODULE mLibrary = nullptr;
#else
  void *mLibrary = nullptr;
#endif
  fmi2Component mInstance = nullptr;
  bool mValid = true;
  bool mTerminated = false;
  double mTime = 0;
  std::string mLastMessage;

#define FMI2_FUNCTION(name) fmi2##name##TYPE *name = nullptr
  FMI2_FUNCTION(Instantiate);
  FMI2_FUNCTION(GetVersion);
  FMI2_FUNCTION(FreeInstance);
  FMI2_FUNCTION(Terminate);
  FMI2_FUNCTION(SetupExperiment);
  FMI2_FUNCTION(EnterInitializationMode);
  FMI2_FUNCTION(ExitInitializationMode);
  FMI2_FUNCTION(SetReal);
  FMI2_FUNCTION(GetReal);
  FMI2_FUNCTION(DoStep);
  FMI2_FUNCTION(CancelStep);
#undef FMI2_FUNCTION

  template <class T> T symbol(const char *name, bool required = true) {
#ifdef _WIN32
    auto address = GetProcAddress(mLibrary, name);
#else
    auto address = dlsym(mLibrary, name);
#endif
    if (!address && required)
      fail(std::string("missing FMI 2 function: ") + name);
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
      static_cast<Fmi2CoSimulation *>(environment)->mLastMessage =
          message ? message : "";
    } catch (...) {
    }
  }
  static void *allocate(size_t n, size_t size) { return std::calloc(n, size); }
  static void release(void *pointer) { std::free(pointer); }
  void close() noexcept {
    if (mInstance) {
      if (!mTerminated && mValid && Terminate)
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

public:
  Fmi2CoSimulation(const std::string &path,
                   const std::vector<ScalarInput> &configuration,
                   CoSimulationOptions options) {
    try {
      initialize(path, configuration, options);
    } catch (...) {
      close();
      throw;
    }
  }
  ~Fmi2CoSimulation() override { close(); }
  void initialize(const std::string &path,
                  const std::vector<ScalarInput> &configuration,
                  CoSimulationOptions options) {
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
    auto cs = child(root, "CoSimulation");
    const auto identifier = property(cs, "modelIdentifier");
    const auto guid = property(root, "guid");
    if (identifier.empty() || guid.empty())
      fail("invalid FMI 2 modelIdentifier or guid");
    model.modelName = property(root, "modelName", identifier);
    model.modelIdentifier = identifier;
    model.instantiationToken = guid;
    model.canGetAndSetFMUState = flag(cs, "canGetAndSetFMUstate");
    model.canSerializeFMUState = flag(cs, "canSerializeFMUstate");
    model.canBeInstantiatedOnlyOncePerProcess =
        flag(cs, "canBeInstantiatedOnlyOncePerProcess");
    if (flag(cs, "canRunAsynchronuously"))
      fail("asynchronous FMI 2 Co-Simulation FMUs are unsupported");
    if (configuration.size() > 64)
      fail("input count exceeds maximum 64");
    const bool automaticInputs = configuration.empty();
    std::map<std::string, std::size_t> mapped;
    std::vector<bool> found(configuration.size(), false);
    if (!automaticInputs) {
      model.initialInputs = Matrix::Zero(configuration.size(), 1);
      for (std::size_t i = 0; i < configuration.size(); ++i) {
        if (configuration[i].name.empty() ||
            !std::isfinite(configuration[i].initialValue) ||
            !mapped.emplace(configuration[i].name, i).second)
          fail("invalid or duplicate input mapping");
        inputNames.push_back(configuration[i].name);
        model.initialInputs(i, 0) = configuration[i].initialValue;
      }
    }
    std::set<std::string> names;
    std::vector<double> starts;
    auto variables = child(root, "ModelVariables");
    for (auto variable = variables->children; variable;
         variable = variable->next) {
      if (!tag(variable, "ScalarVariable"))
        continue;
      const auto name = property(variable, "name");
      if (!names.insert(name).second)
        fail("duplicate FMI variable name");
      const auto causality = property(variable, "causality", "local");
      xmlNode *real = nullptr;
      for (auto type = variable->children; type; type = type->next)
        if (tag(type, "Real")) {
          real = type;
          break;
        }
      if (causality != "input" && causality != "output")
        continue;
      if (!real)
        fail("requires scalar Real Co-Simulation inputs and outputs");
      for (auto dimension = variable->children; dimension;
           dimension = dimension->next)
        if (tag(dimension, "Dimension"))
          fail("array-valued Co-Simulation inputs and outputs are unsupported");
      const auto ref = reference(property(variable, "valueReference"));
      if (causality == "input") {
        if (automaticInputs) {
          inputNames.push_back(name);
          inputReferences.push_back(ref);
          const auto text = property(real, "start", "0");
          std::size_t consumed = 0;
          const double value = std::stod(text, &consumed);
          if (consumed != text.size() || !std::isfinite(value))
            fail("invalid Real input start for " + name);
          starts.push_back(value);
        } else {
          const auto item = mapped.find(name);
          if (item != mapped.end()) {
            const auto index = item->second;
            if (found[index])
              fail("duplicate mapped input variable " + name);
            if (inputReferences.empty())
              inputReferences.resize(inputNames.size());
            inputReferences[index] = ref;
            found[index] = true;
          }
        }
      } else {
        outputReferences.push_back(ref);
        model.outputNames.push_back(name);
      }
    }
    if (inputReferences.size() != inputNames.size())
      fail("all mapped input variables must exist exactly once");
    if (!automaticInputs &&
        std::find(found.begin(), found.end(), false) != found.end())
      fail("mapped input does not exist");
    if (automaticInputs) {
      model.initialInputs.resize(starts.size(), 1);
      for (std::size_t i = 0; i < starts.size(); ++i)
        model.initialInputs(i, 0) = starts[i];
    }
    model.inputNames = inputNames;
    model.initialState.resize(0, 1);
    model.stateNominals.resize(0, 1);
    if (model.outputNames.empty())
      fail("FMU has no Real outputs");
    // FMI 2.0 specifies platform folders such as linux64; FMI 3.0 uses
    // architecture-OS triplets handled by libraryPath above.
    const auto binary = fmi2LibraryPath(mArchive.directory, identifier);
    if (!fs::is_regular_file(binary))
      fail("missing host binary: " +
           binary.lexically_relative(mArchive.directory).string());
#ifdef _WIN32
    mLibrary = LoadLibraryW(binary.wstring().c_str());
#else
    mLibrary = dlopen(binary.string().c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
    if (!mLibrary)
      fail("cannot load FMU binary");
#define LOAD2(name) name = symbol<decltype(name)>("fmi2" #name)
    LOAD2(Instantiate);
    LOAD2(GetVersion);
    LOAD2(FreeInstance);
    LOAD2(Terminate);
    LOAD2(SetupExperiment);
    LOAD2(EnterInitializationMode);
    LOAD2(ExitInitializationMode);
    LOAD2(SetReal);
    LOAD2(GetReal);
    LOAD2(DoStep);
#undef LOAD2
    CancelStep = symbol<decltype(CancelStep)>("fmi2CancelStep", false);
    if (!GetVersion() || std::string(GetVersion()) != "2.0")
      fail("binary does not implement FMI 2.0");
    enforceSingleInstantiation(model.canBeInstantiatedOnlyOncePerProcess,
                               "fmi2:" + identifier + ":" + guid);
    fmi2CallbackFunctions callbacks{log, allocate, release, nullptr, this};
    const bool visible = false;
    mInstance = Instantiate(model.modelName.c_str(), 1, guid.c_str(), nullptr,
                            &callbacks, visible ? fmi2True : fmi2False,
                            options.loggingOn ? fmi2True : fmi2False);
    if (!mInstance)
      fail("instantiation failed: " + mLastMessage);
    check(SetupExperiment(mInstance, fmi2False, 0, 0, fmi2False, 0),
          "setup experiment");
    check(EnterInitializationMode(mInstance), "enter initialization mode");
    setInputs(model.initialInputs);
    check(ExitInitializationMode(mInstance), "exit initialization mode");
    inputs = model.initialInputs;
    outputs.resize(model.outputNames.size(), 1);
    check(GetReal(mInstance, outputReferences.data(), outputReferences.size(),
                  outputs.data()),
          "read initial outputs");
    if (!outputs.allFinite())
      fail("initial outputs contain non-finite values");
  }
  void setInputs(const Matrix &values) override {
    if (values.rows() != static_cast<Eigen::Index>(inputReferences.size()) ||
        values.cols() != 1 || !values.allFinite())
      fail("input vector has invalid dimensions or values");
    if (!inputReferences.empty())
      check(SetReal(mInstance, inputReferences.data(), inputReferences.size(),
                    values.data()),
            "set Co-Simulation inputs");
    inputs = values;
  }
  void doStep(double currentTime, double stepSize) override {
    if (!std::isfinite(currentTime) || !std::isfinite(stepSize) ||
        stepSize <= 0 ||
        std::abs(currentTime - mTime) > 1e-10 * std::max(1.0, std::abs(mTime)))
      fail("communication point must match current FMU time and step size be "
           "positive");
    const auto status = DoStep(mInstance, currentTime, stepSize, fmi2True);
    if (status == fmi2Pending) {
      if (CancelStep)
        CancelStep(mInstance);
      mValid = false;
      fail("asynchronous doStep returned Pending; asynchronous stepping is "
           "unsupported");
    }
    check(status, "doStep");
    if (status == fmi2Discard) {
      mValid = false;
      fail("FMU discarded doStep; partial communication steps are unsupported");
    }
    mTime = currentTime + stepSize;
    check(GetReal(mInstance, outputReferences.data(), outputReferences.size(),
                  outputs.data()),
          "read outputs after doStep");
    if (!outputs.allFinite())
      fail("output vector contains non-finite values");
  }
  double time() const override { return mTime; }
};
} // namespace

CoSimulation::CoSimulation(const std::string &path, CoSimulationOptions options)
    : CoSimulation(path, std::vector<ScalarInput>{}, options) {}
CoSimulation::CoSimulation(const std::string &path,
                           const std::vector<ScalarInput> &inputs,
                           CoSimulationOptions options) {
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
    fail("invalid fmiModelDescription root");
  const auto version = property(root, "fmiVersion");
  if (version == "3.0")
    mImpl = std::make_unique<Fmi3CoSimulation>(path, inputs, options);
  else if (version == "2.0")
    mImpl = std::make_unique<Fmi2CoSimulation>(path, inputs, options);
  else
    fail("supports FMI 2.0 and FMI 3.0 Co-Simulation only");
}
CoSimulation::~CoSimulation() = default;
const ModelDescription &CoSimulation::description() const {
  return mImpl->model;
}
const Matrix &CoSimulation::inputs() const { return mImpl->inputs; }
const Matrix &CoSimulation::outputs() const { return mImpl->outputs; }
double CoSimulation::time() const { return mImpl->time(); }
void CoSimulation::setInputs(const Matrix &inputs) { mImpl->setInputs(inputs); }
void CoSimulation::setInput(const std::string &name, double value) {
  if (!std::isfinite(value))
    fail("input value must be finite");
  const auto found = std::find(mImpl->model.inputNames.begin(),
                               mImpl->model.inputNames.end(), name);
  if (found == mImpl->model.inputNames.end())
    fail("input does not exist: " + name);
  Matrix next = mImpl->inputs;
  next(static_cast<Eigen::Index>(found - mImpl->model.inputNames.begin()), 0) =
      value;
  mImpl->setInputs(next);
}
const Matrix &CoSimulation::doStep(double currentCommunicationPoint,
                                   double communicationStepSize) {
  mImpl->doStep(currentCommunicationPoint, communicationStepSize);
  return mImpl->outputs;
}
} // namespace DPsim::FMI
