// SPDX-License-Identifier: MPL-2.0
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_Switch.h>
#include <dpsim-models/EMT/EMT_Ph3_VoltageSource.h>
#include <dpsim/DataLogger.h>
#include <dpsim/Event.h>
#include <dpsim/Simulation.h>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

// Real MNA connection: voltage source -- native R/L -- PCC -- FMU and native
// load.
using namespace DPsim;
using namespace CPS;

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "Usage: FMI_EMT_Ph3_Mixed model.fmu config.json "
                 "[stop_time_s [step_s]]\n";
    return 2;
  }
  try {

    // ----- PARAMETERS AND FMU INTERFACE -----

    nlohmann::json config;
    std::ifstream file(argv[2]);
    file >> config;
    const double frequency = config.at("frequency_Hz");
    const double phaseVoltagePeak = config.at("voltage_peak_V");
    const double sourcePhase = config.value("phase_rad", 0.0);
    const double timeStep =
        argc == 5 ? std::stod(argv[4]) : config.at("step_s").get<double>();
    const double finalTime =
        argc >= 4 ? std::stod(argv[3]) : config.at("stop_time_s").get<double>();
    const double lineResistance = config.at("series_R_ohm");
    const double lineInductance = config.value("series_L_H", 0.0);
    const double loadResistance = config.at("load_R_ohm");
    const double initialVoltageScale =
        config.value("initial_voltage_scale", 1.0);
    if (timeStep <= 0 || finalTime <= 0 || lineResistance <= 0 ||
        lineInductance < 0 || loadResistance <= 0 ||
        std::abs(finalTime / timeStep - std::round(finalTime / timeStep)) >
            1e-7) {
      throw std::invalid_argument("Invalid network parameters or time grid");
    }

    // ----- GRID DISTURBANCES -----

    const Real gridEventTime = config.value("grid_event_time_s", 0.17);
    const Real sagEndTime = config.value("grid_sag_end_s", 0.19);
    const Real sagScale = config.value("grid_sag_scale", 1.0);
    std::array<Real, 3> sagScaleABC{sagScale, sagScale, sagScale};
    if (config.contains("grid_sag_scale_abc"))
      sagScaleABC = config.at("grid_sag_scale_abc").get<std::array<Real, 3>>();
    for (const auto scale : sagScaleABC)
      if (!std::isfinite(scale) || scale < 0)
        throw std::invalid_argument(
            "Grid sag scales must be finite and nonnegative");
    const Real phaseJump = config.value("grid_phase_jump_rad", 0.0);
    const Real frequencyChange = config.value("grid_frequency_change_Hz", 0.0);
    const auto voltageScale = [=](Real time, int phase) {
      return time >= gridEventTime - 1e-12 && time < sagEndTime - 1e-12
                 ? sagScaleABC[phase]
                 : 1.0;
    };
    const auto phaseShift = [=](Real time) {
      return time >= gridEventTime - 1e-12
                 ? phaseJump + 2 * PI * frequencyChange * (time - gridEventTime)
                 : 0.0;
    };

    // ----- INITIAL CONDITIONS -----

    const auto gridVoltage = [=](double time) -> Matrix {
      Matrix phaseVoltages(3, 1);
      for (int phase = 0; phase < 3; ++phase) {
        phaseVoltages(phase, 0) =
            voltageScale(time, phase) * phaseVoltagePeak *
            std::cos(2 * PI * frequency * time + sourcePhase +
                     phaseShift(time) - phase * 2 * PI / 3);
      }
      return phaseVoltages;
    };
    const bool initializeLineCurrentFromLoad = config.value(
        "initialize_line_current_from_load", std::abs(initialVoltageScale) > 0);
    const Matrix initialGridVoltage = gridVoltage(0) * initialVoltageScale;
    Matrix initialVoltage = initialGridVoltage;
    Matrix initialLineCurrent = Matrix::Zero(3, 1);
    Matrix initialResistorVoltage = Matrix::Zero(3, 1);
    if (initializeLineCurrentFromLoad) {
      const double omega = 2 * PI * frequency;
      const Complex lineImpedance(lineResistance, omega * lineInductance);
      for (int phase = 0; phase < 3; ++phase) {
        const Complex gridPhasor =
            std::polar(phaseVoltagePeak * initialVoltageScale,
                       sourcePhase - phase * 2 * PI / 3);
        const Complex pccPhasor =
            gridPhasor * loadResistance / (loadResistance + lineImpedance);
        const Complex lineCurrentPhasor = -pccPhasor / loadResistance;
        initialVoltage(phase, 0) = pccPhasor.real();
        initialLineCurrent(phase, 0) = lineCurrentPhasor.real();
        initialResistorVoltage(phase, 0) =
            (lineResistance * lineCurrentPhasor).real();
      }
    }
    // ----- IMPORTED CONVERTER -----

    // The common v_a/v_b/v_c and i_a/i_b/i_c interface is inferred. An FMU
    // with nonstandard names can provide the optional fmi_mapping object.
    FMI::ModelExchangeOptions options;
    options.endpointStateEvents = config.value("endpoint_state_events", false);
    options.allowContinuousStateResets =
        config.value("allow_continuous_state_resets", false);
    FMI::NewtonSettings newtonSettings;
    const std::string multirateMode =
        config.value("multirate_mode", std::string("rollback_newton"));
    if (multirateMode == "rollback_newton") {
      newtonSettings.multirateMode = FMI::MultirateMode::RollbackNewton;
    } else if (multirateMode == "single_pass_explicit") {
      newtonSettings.multirateMode = FMI::MultirateMode::SinglePassExplicit;
    } else {
      throw std::invalid_argument(
          "multirate_mode must be rollback_newton or single_pass_explicit");
    }
    // The x86-64 Linux FMU exported by MathWorks R2025a retains substantial
    // memory across repeated FMI state snapshot/free cycles. Single-pass
    // coupling does not need provisional rollback, so leave it disabled unless
    // a configuration explicitly asks for snapshot-backed trial evaluations.
    options.snapshotTrialEvaluations = config.value(
        "snapshot_trial_evaluations",
        newtonSettings.multirateMode != FMI::MultirateMode::SinglePassExplicit);
    newtonSettings.allowInstantaneousCurrentEvents =
        config.value("allow_instantaneous_current_events", false);
    newtonSettings.maxJacobianReuseSteps =
        config.value("jacobian_reuse_steps", 0);
    const double toleranceScale = config.value("tolerance_scale", 1.0);
    newtonSettings.absoluteStateTolerance *= toleranceScale;
    newtonSettings.relativeStateTolerance *= toleranceScale;
    newtonSettings.absoluteCurrentTolerance *= toleranceScale;
    newtonSettings.relativeCurrentTolerance *= toleranceScale;
    newtonSettings.multirateAbsoluteCurrentTolerance *= toleranceScale;
    newtonSettings.multirateRelativeCurrentTolerance *= toleranceScale;
    FMI::ModelExchange::Ptr model;
    std::shared_ptr<FMI::NonlinearSSN> converter;
    if (config.contains("fmi_mapping")) {
      const auto &mapping = config.at("fmi_mapping");
      const auto voltageInputs =
          mapping.at("voltage_inputs").get<std::array<std::string, 3>>();
      const auto currentOutputs =
          mapping.at("current_outputs").get<std::array<std::string, 3>>();
      std::vector<FMI::ScalarInput> inputs;
      for (int phase = 0; phase < 3; ++phase)
        inputs.push_back({voltageInputs[phase], initialVoltage(phase, 0)});
      model = std::make_shared<FMI::ModelExchange>(argv[1], inputs, options);
      converter = FMI::NonlinearSSN::make(
          "converter", model,
          FMI::ElectricalPort3Ph{voltageInputs, currentOutputs,
                                 mapping.value("current_sign_to_passive", 1.0)},
          newtonSettings);
    } else {
      converter = FMI::importFMU(argv[1], options, newtonSettings);
    }

    // ----- NODES -----

    const auto initialNodeVoltage = [](const Matrix &voltage) {
      std::vector<Complex> result(3);
      for (int phase = 0; phase < 3; ++phase)
        // SimNode<Real> converts its three-phase initial phasors to phase
        // voltages with RMS3PH_TO_PEAK1PH during initialization.
        result[phase] = voltage(phase, 0) / RMS3PH_TO_PEAK1PH;
      return result;
    };
    auto gridNode = SimNode<Real>::make("grid", PhaseType::ABC,
                                        initialNodeVoltage(initialGridVoltage));
    auto pccNode = SimNode<Real>::make("pcc", PhaseType::ABC,
                                       initialNodeVoltage(initialVoltage));
    auto lineNode = SimNode<Real>::make(
        "line_rl", PhaseType::ABC,
        initialNodeVoltage(initialGridVoltage + initialResistorVoltage));

    // ----- NATIVE EMT COMPONENTS -----

    auto gridSource = EMT::Ph3::VoltageSource::make("grid_source");
    MatrixComp phasors(3, 1);
    for (int phase = 0; phase < 3; ++phase) {
      phasors(phase, 0) = std::polar(phaseVoltagePeak / RMS3PH_TO_PEAK1PH,
                                     sourcePhase - phase * 2 * PI / 3);
    }
    gridSource->setParameters(phasors, frequency);
    auto lineResistor = EMT::Ph3::Resistor::make("line_resistor");
    lineResistor->setParameters(lineResistance * Matrix::Identity(3, 3));
    auto nativeLoad = EMT::Ph3::Resistor::make("native_load");
    nativeLoad->setParameters(loadResistance * Matrix::Identity(3, 3));
    std::shared_ptr<EMT::Ph3::Inductor> lineInductor;
    if (lineInductance > 0) {
      lineInductor = EMT::Ph3::Inductor::make("line_inductor");
      lineInductor->setParameters(lineInductance * Matrix::Identity(3, 3));
    }

    // ----- TOPOLOGY -----

    gridSource->connect({SimNode<Real>::GND, gridNode});
    lineResistor->connect({gridNode, lineInductance > 0 ? lineNode : pccNode});
    nativeLoad->connect({SimNode<Real>::GND, pccNode});
    converter->connect(pccNode);

    SystemTopology system(frequency, {gridNode, pccNode},
                          {gridSource, lineResistor, nativeLoad, converter});
    if (lineInductor) {
      lineInductor->connect({lineNode, pccNode});
      system.addNode(lineNode);
      system.addComponent(lineInductor);
    }

    // ----- OPTIONAL SWITCHED LOAD -----

    std::shared_ptr<EMT::Ph3::Switch> loadSwitch;
    std::shared_ptr<EMT::Ph3::Resistor> additionalLoad;
    if (config.contains("step_load_R_ohm_abc")) {
      const auto phaseResistances =
          config.at("step_load_R_ohm_abc").get<std::array<double, 3>>();
      Matrix stepLoadResistance = Matrix::Zero(3, 3);
      for (int phase = 0; phase < 3; ++phase) {
        if (!std::isfinite(phaseResistances[phase]) ||
            phaseResistances[phase] <= 0) {
          throw std::invalid_argument("Step-load resistance must be positive");
        }
        stepLoadResistance(phase, phase) = phaseResistances[phase];
      }
      auto stepLoadNode = SimNode<Real>::make("step_load_node", PhaseType::ABC);
      additionalLoad = EMT::Ph3::Resistor::make("step_load_resistor");
      additionalLoad->setParameters(stepLoadResistance);
      additionalLoad->connect({SimNode<Real>::GND, stepLoadNode});
      loadSwitch = EMT::Ph3::Switch::make("native_load_switch");
      loadSwitch->setParameters(1e12 * Matrix::Identity(3, 3),
                                1e-3 * Matrix::Identity(3, 3), false);
      loadSwitch->connect({stepLoadNode, pccNode});
      system.addNode(stepLoadNode);
      system.addComponent(additionalLoad);
      system.addComponent(loadSwitch);
    }

    // ----- LOGGING -----

    const UInt logDownsampling = config.value("log_downsampling", 1u);
    auto logger = DataLogger::make(String("fmi_mixed"), true, logDownsampling);
    if (loadSwitch) {
      logger->logAttribute("step_load_current",
                           additionalLoad->attribute("i_intf"));
    }
    for (const auto &name : converter->modelDescription().outputNames) {
      logger->logAttribute(name, converter->attribute("fmi.output." + name));
    }
    logger->logAttribute("grid_voltage", gridSource->attribute("v_intf"));
    logger->logAttribute("pcc_voltage", pccNode->attribute("v"));
    logger->logAttribute("converter_current", converter->attribute("i_intf"));
    logger->logAttribute("line_current", lineResistor->attribute("i_intf"));
    logger->logAttribute("load_current", nativeLoad->attribute("i_intf"));
    logger->logAttribute("newton_iterations",
                         converter->attribute("newton_iterations"));
    logger->logAttribute("newton_residual",
                         converter->attribute("newton_residual"));
    logger->logAttribute("fmi_internal_steps",
                         converter->attribute("fmi_internal_steps"));
    logger->logAttribute("fmi_coupling_error_A",
                         converter->attribute("fmi_coupling_error_A"));
    logger->logAttribute("fmi_evaluation_requests",
                         converter->attribute("fmi_evaluation_requests"));
    logger->logAttribute("fmi_evaluation_wall_time_s",
                         converter->attribute("fmi_evaluation_wall_time_s"));

    // ----- SIMULATION -----

    Simulation simulation("fmi_mixed", Logger::Level::off);
    simulation.setDomain(Domain::EMT);
    simulation.setSystem(system);
    simulation.setTimeStep(timeStep);
    simulation.setFinalTime(finalTime);
    simulation.doSteadyStateInit(false);
    simulation.addLogger(logger);

    // ----- LOAD EVENT -----

    if (loadSwitch) {
      const double loadStepTime = config.at("load_step_time_s");
      if (loadStepTime <= 0 || loadStepTime >= finalTime ||
          std::abs(loadStepTime / timeStep -
                   std::round(loadStepTime / timeStep)) > 1e-7) {
        throw std::invalid_argument(
            "Load event must align with the integration grid");
      }
      simulation.addEvent(SwitchEvent3Ph::make(loadStepTime, loadSwitch, true));
    }

    // ----- CONSISTENTLY INITIALIZE THE LOADED NETWORK AND RUN -----

    // Initialize components first so that the FMU's initial terminal current
    // is available. Then seed the line, resistor, and load histories before
    // Simulation::start logs the t=0 sample. The line current satisfies PCC
    // KCL for the configured load and the FMU's initialized current.
    simulation.initialize();
    **gridSource->attributeTyped<Matrix>("v_intf") = initialGridVoltage;
    if (lineInductor) {
      if (initializeLineCurrentFromLoad) {
        // Interface currents are oriented from terminal 0 to terminal 1.
        // The line therefore carries the opposite of the load and converter
        // currents computed at the initialized PCC voltage.
        const Matrix initialLoadCurrent = initialVoltage / loadResistance;
        initialLineCurrent = -initialLoadCurrent;
        initialLineCurrent -= **converter->attributeTyped<Matrix>("i_intf");
        **nativeLoad->attributeTyped<Matrix>("v_intf") = initialVoltage;
        **nativeLoad->attributeTyped<Matrix>("i_intf") = initialLoadCurrent;
        const Matrix kclResidual =
            initialLineCurrent + initialLoadCurrent +
            **converter->attributeTyped<Matrix>("i_intf");
        if (!kclResidual.allFinite() ||
            kclResidual.lpNorm<Eigen::Infinity>() > 1e-10) {
          throw std::runtime_error(
              "Loaded-start initialization does not satisfy PCC KCL");
        }
      }
      initialResistorVoltage = lineResistance * initialLineCurrent;
      **gridNode->attributeTyped<Matrix>("v") = initialGridVoltage;
      **lineNode->attributeTyped<Matrix>("v") =
          initialGridVoltage + initialResistorVoltage;
      **pccNode->attributeTyped<Matrix>("v") = initialVoltage;
      **lineResistor->attributeTyped<Matrix>("v_intf") = initialResistorVoltage;
      **lineResistor->attributeTyped<Matrix>("i_intf") = initialLineCurrent;
      **lineInductor->attributeTyped<Matrix>("v_intf") =
          initialVoltage - initialGridVoltage - initialResistorVoltage;
      **lineInductor->attributeTyped<Matrix>("i_intf") = initialLineCurrent;
    }
    simulation.start();
    const auto numberOfSteps = std::llround(finalTime / timeStep);
    for (long long step = 0; step < numberOfSteps; ++step) {
      // Native source keeps its continuous base-frequency oscillator. Updating
      // its phasor adds the prescribed amplitude and continuous phase offset.
      const Real time = (step + 1) * timeStep;
      MatrixComp disturbedPhasors = phasors;
      for (int phase = 0; phase < 3; ++phase) {
        disturbedPhasors(phase, 0) *=
            std::polar(voltageScale(time, phase), phaseShift(time));
      }
      **gridSource->mVoltageRef = disturbedPhasors;
      simulation.step();
    }
    simulation.stop();
    std::cout << "Mixed native EMT ABC / FMI completed: "
              << converter->modelDescription().stateNames.size()
              << " FMU states, " << std::llround(finalTime / timeStep)
              << " steps\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
