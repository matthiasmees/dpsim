// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/NonlinearSSN.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

void addFMIBindings(py::module_ root) {
  auto fmi = root.def_submodule("fmi", "FMI Model Exchange components");

  py::enum_<DPsim::FMI::MultirateMode>(fmi, "MultirateMode")
      .value("RollbackNewton", DPsim::FMI::MultirateMode::RollbackNewton)
      .value("SinglePassExplicit",
             DPsim::FMI::MultirateMode::SinglePassExplicit);

  py::class_<DPsim::FMI::ModelExchangeOptions>(fmi, "ModelExchangeOptions")
      .def(py::init<>())
      .def_readwrite("endpoint_state_events",
                     &DPsim::FMI::ModelExchangeOptions::endpointStateEvents)
      .def_readwrite(
          "allow_continuous_state_resets",
          &DPsim::FMI::ModelExchangeOptions::allowContinuousStateResets)
      .def_readwrite(
          "snapshot_trial_evaluations",
          &DPsim::FMI::ModelExchangeOptions::snapshotTrialEvaluations);

  py::class_<DPsim::FMI::NewtonSettings>(fmi, "NewtonSettings")
      .def(py::init<>())
      .def_readwrite("absolute_state_tolerance",
                     &DPsim::FMI::NewtonSettings::absoluteStateTolerance)
      .def_readwrite("relative_state_tolerance",
                     &DPsim::FMI::NewtonSettings::relativeStateTolerance)
      .def_readwrite("absolute_current_tolerance",
                     &DPsim::FMI::NewtonSettings::absoluteCurrentTolerance)
      .def_readwrite("relative_current_tolerance",
                     &DPsim::FMI::NewtonSettings::relativeCurrentTolerance)
      .def_readwrite("relative_jacobian_perturbation",
                     &DPsim::FMI::NewtonSettings::relativeJacobianPerturbation)
      .def_readwrite("max_iterations",
                     &DPsim::FMI::NewtonSettings::maxIterations)
      .def_readwrite("automatic_event_substepping",
                     &DPsim::FMI::NewtonSettings::automaticEventSubstepping)
      .def_readwrite("multirate_mode",
                     &DPsim::FMI::NewtonSettings::multirateMode)
      .def_readwrite("max_internal_steps",
                     &DPsim::FMI::NewtonSettings::maxInternalSteps);

  py::class_<DPsim::FMI::ScalarInput>(fmi, "ScalarInput")
      .def(py::init<std::string, double>(), "name"_a, "initial_value"_a)
      .def_readwrite("name", &DPsim::FMI::ScalarInput::name)
      .def_readwrite("initial_value", &DPsim::FMI::ScalarInput::initialValue);

  py::class_<DPsim::FMI::ModelDescription>(fmi, "ModelDescription")
      .def_readonly("model_name", &DPsim::FMI::ModelDescription::modelName)
      .def_readonly("state_names", &DPsim::FMI::ModelDescription::stateNames)
      .def_readonly("input_names", &DPsim::FMI::ModelDescription::inputNames)
      .def_readonly("output_names", &DPsim::FMI::ModelDescription::outputNames)
      .def_readonly("can_get_and_set_fmu_state",
                    &DPsim::FMI::ModelDescription::canGetAndSetFMUState);

  py::class_<DPsim::FMI::ElectricalPort3Ph>(fmi, "ElectricalPort3Ph")
      .def(py::init<>())
      .def_readwrite("voltage_inputs",
                     &DPsim::FMI::ElectricalPort3Ph::voltageInputs)
      .def_readwrite("current_outputs",
                     &DPsim::FMI::ElectricalPort3Ph::currentOutputs)
      .def_readwrite("current_sign_to_passive",
                     &DPsim::FMI::ElectricalPort3Ph::currentSignToPassive);

  py::class_<DPsim::FMI::ModelExchange, DPsim::FMI::ModelExchange::Ptr>(
      fmi, "ModelExchange")
      .def(py::init<const std::string &, DPsim::FMI::ModelExchangeOptions>(),
           "path"_a, "options"_a = DPsim::FMI::ModelExchangeOptions{})
      .def(py::init<const std::string &,
                    const std::vector<DPsim::FMI::ScalarInput> &,
                    DPsim::FMI::ModelExchangeOptions>(),
           "path"_a, "inputs"_a,
           "options"_a = DPsim::FMI::ModelExchangeOptions{})
      .def_property_readonly("description",
                             &DPsim::FMI::ModelExchange::description,
                             py::return_value_policy::reference_internal)
      .def("supports_rollback", &DPsim::FMI::ModelExchange::supportsRollback);

  py::class_<DPsim::FMI::NonlinearSSN,
             std::shared_ptr<DPsim::FMI::NonlinearSSN>,
             CPS::SimPowerComp<CPS::Real>>(fmi, "NonlinearFMU",
                                           py::multiple_inheritance())
      .def("connect",
           static_cast<void (DPsim::FMI::NonlinearSSN::*)(
               CPS::SimNode<CPS::Real>::Ptr)>(
               &DPsim::FMI::NonlinearSSN::connect),
           "node"_a,
           "Connect the FMU between ground and this three-phase node.")
      .def(
          "connect_terminals",
          [](DPsim::FMI::NonlinearSSN &component,
             CPS::SimNode<CPS::Real>::List nodes) {
            component.CPS::SimPowerComp<CPS::Real>::connect(std::move(nodes));
          },
          "nodes"_a, "Connect explicit negative and positive terminal nodes.")
      .def_property_readonly("description",
                             &DPsim::FMI::NonlinearSSN::modelDescription,
                             py::return_value_policy::reference_internal)
      .def("set_input_function", &DPsim::FMI::NonlinearSSN::setInputFunction,
           "name"_a, "function"_a,
           "discontinuity_times"_a = std::vector<double>{});

  fmi.def(
      "import_fmu",
      [](const std::string &path, const std::string &name,
         const py::object &optionsObject,
         const DPsim::FMI::NewtonSettings &settings) {
        auto options =
            optionsObject.is_none()
                ? DPsim::FMI::ModelExchangeOptions{}
                : optionsObject.cast<DPsim::FMI::ModelExchangeOptions>();
        if (optionsObject.is_none())
          options.endpointStateEvents = true;
        if (name.empty())
          return DPsim::FMI::importFMU(path, options, settings);
        return DPsim::FMI::NonlinearSSN::importFMU(name, path, options,
                                                   settings);
      },
      "path"_a, "name"_a = "", "model_options"_a = py::none(),
      "settings"_a = DPsim::FMI::NewtonSettings{},
      "Import a conventional abc Model Exchange FMU as an EMT component.");

  fmi.def(
      "import_fmu_with_mapping",
      [](const std::string &path, const std::vector<std::string> &voltageInputs,
         const std::vector<std::string> &currentOutputs, double currentSign,
         const std::string &name,
         const std::vector<DPsim::FMI::ScalarInput> &initialInputs,
         const DPsim::FMI::ModelExchangeOptions &options,
         const DPsim::FMI::NewtonSettings &settings) {
        if (voltageInputs.size() != 3 || currentOutputs.size() != 3)
          throw std::invalid_argument(
              "three voltage inputs and three current outputs are required");
        auto model =
            initialInputs.empty()
                ? std::make_shared<DPsim::FMI::ModelExchange>(path, options)
                : std::make_shared<DPsim::FMI::ModelExchange>(
                      path, initialInputs, options);
        DPsim::FMI::ElectricalPort3Ph port{
            {voltageInputs[0], voltageInputs[1], voltageInputs[2]},
            {currentOutputs[0], currentOutputs[1], currentOutputs[2]},
            currentSign};
        const auto componentName =
            name.empty() ? model->description().modelName : name;
        return DPsim::FMI::NonlinearSSN::make(componentName, std::move(model),
                                              port, settings);
      },
      "path"_a, "voltage_inputs"_a, "current_outputs"_a,
      "current_sign_to_passive"_a = 1.0, "name"_a = "",
      "initial_inputs"_a = std::vector<DPsim::FMI::ScalarInput>{},
      "model_options"_a = DPsim::FMI::ModelExchangeOptions{},
      "settings"_a = DPsim::FMI::NewtonSettings{},
      "Import an FMU with explicit terminal variable names and polarity.");
}
