// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"
#include "../GeneratorFactory.h"

#include <DPsim.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace DPsim;
using namespace CPS;

namespace {

// Controller parameters and their equation-based derivation live in Examples.h.
using GfmParams = CPS::CIM::Examples::Components::GFM::Ieee9SsnGridForming;
using GflParams = CPS::CIM::Examples::Components::GFL::Ieee9AvVsi;

using SynchronousGenerator = CPS::EMT::Ph3::SynchronGenerator4OrderVBR;
using GridFormingInverter = CPS::EMT::Ph3::SSN_GFM;
using GridFollowingInverter = CPS::EMT::Ph3::SSN_GFL;

struct ScenarioArtifacts {
  SystemTopology system;
  std::vector<SimNode<Real>::Ptr> buses;

  std::shared_ptr<SynchronousGenerator> synchronousGenerator;
  std::shared_ptr<GridFormingInverter> gridFormingInverter;
  std::shared_ptr<GridFollowingInverter> gridFollowingInverter;

  // Line-to-line RMS base voltages of BUS1 ... BUS9.
  std::array<Real, 9> busBaseVoltages;

  Real nominalOmega;

  Real synchronousGeneratorRatedPower;
  Real synchronousGeneratorRatedVoltage;

  Real gridFormingRatedPower;
  Real gridFormingRatedVoltage;

  Real gridFollowingRatedPower;
  Real gridFollowingRatedVoltage;

  Complex gridFormingPowerFlowPower;
  Complex gridFollowingPccPower;
};

struct DiagnosticSettings {
  Int csvEverySteps = 10;
  Int consoleEverySteps = 1000;

  // These limits are deliberately broad. Their purpose is to identify the
  // first diverging subsystem, not to implement protection.
  Real maximumVoltagePu = 3.0;
  Real maximumCurrentPu = 20.0;
  Real minimumFrequencyPu = 0.5;
  Real maximumFrequencyPu = 1.5;
  Real maximumStateGrowthFactor = 1e8;

  Bool abortOnLimit = true;
};

Real maxAbs(const Matrix &values) {
  if (values.size() == 0)
    return 0.0;

  return values.cwiseAbs().maxCoeff();
}

Real peakPhaseVoltage(Real lineToLineRmsVoltage) {
  return RMS3PH_TO_PEAK1PH * lineToLineRmsVoltage;
}

Real peakPhaseCurrent(Real apparentPower, Real lineToLineRmsVoltage) {
  if (apparentPower <= 0.0 || lineToLineRmsVoltage <= 0.0)
    return 1.0;

  return std::sqrt(2.0) * apparentPower /
         (std::sqrt(3.0) * lineToLineRmsVoltage);
}

struct LocalModeSummary {
  Bool valid = false;
  Real spectralRadius = std::numeric_limits<Real>::quiet_NaN();
  Real maximumEquivalentContinuousGrowth =
      std::numeric_limits<Real>::quiet_NaN();
  Complex dominantEigenvalue = Complex(std::numeric_limits<Real>::quiet_NaN(),
                                       std::numeric_limits<Real>::quiet_NaN());
};

LocalModeSummary summarizeLocalDiscreteModel(const Matrix &discreteA,
                                             Real timeStep) {
  LocalModeSummary result;

  if (discreteA.rows() == 0 || discreteA.rows() != discreteA.cols() ||
      !discreteA.allFinite() || timeStep <= 0.0)
    return result;

  Eigen::EigenSolver<Matrix> eigenSolver(discreteA, false);

  if (eigenSolver.info() != Eigen::Success)
    return result;

  const auto eigenvalues = eigenSolver.eigenvalues();

  Real spectralRadius = 0.0;
  Real maximumGrowth = -std::numeric_limits<Real>::infinity();
  Complex dominant(0.0, 0.0);

  for (Eigen::Index index = 0; index < eigenvalues.size(); ++index) {
    const Complex eigenvalue = eigenvalues(index);
    const Real magnitude = std::abs(eigenvalue);

    if (magnitude > spectralRadius) {
      spectralRadius = magnitude;
      dominant = eigenvalue;
    }

    const Real equivalentGrowth = magnitude > 0.0
                                      ? std::log(magnitude) / timeStep
                                      : -std::numeric_limits<Real>::infinity();

    maximumGrowth = std::max(maximumGrowth, equivalentGrowth);
  }

  result.valid = true;
  result.spectralRadius = spectralRadius;
  result.maximumEquivalentContinuousGrowth = maximumGrowth;
  result.dominantEigenvalue = dominant;

  return result;
}

class RuntimeDiagnostics {
public:
  RuntimeDiagnostics(const String &simulationName,
                     const ScenarioArtifacts &scenario,
                     const DiagnosticSettings &settings, Real timeStep)
      : mScenario(scenario), mSettings(settings), mTimeStep(timeStep),
        mLog(CPS::Logger::get(simulationName + "_diagnostics",
                              CPS::Logger::Level::info,
                              CPS::Logger::Level::info)) {

    const std::filesystem::path directory =
        std::filesystem::path("logs") / simulationName;
    std::filesystem::create_directories(directory);

    const std::filesystem::path csvPath =
        directory / (simulationName + "_diagnostics.csv");

    mCsv.open(csvPath);

    if (!mCsv.is_open())
      throw std::runtime_error("Could not open diagnostics CSV: " +
                               csvPath.string());

    mCsv << std::setprecision(17);
    writeHeader();

    mLog->info("Runtime diagnostics CSV: {}", csvPath.string());
    mLog->info("Diagnostic limits: Vmax={:.3f} pu, Imax={:.3f} pu, "
               "frequency=[{:.3f},{:.3f}] pu, state-growth abort={:.3e}",
               mSettings.maximumVoltagePu, mSettings.maximumCurrentPu,
               mSettings.minimumFrequencyPu, mSettings.maximumFrequencyPu,
               mSettings.maximumStateGrowthFactor);
  }

  ~RuntimeDiagnostics() {
    if (mCsv.is_open())
      mCsv.flush();
  }

  void logInitialStateDetails() const {
    logNamedState("GEN2 SSN_GFM",
                  mScenario.gridFormingInverter->getLocalStateNames(),
                  mScenario.gridFormingInverter->getState(),
                  mScenario.gridFormingInverter->getStateDerivative());

    logNamedState("GEN3 SSN_GFL",
                  mScenario.gridFollowingInverter->getLocalStateNames(),
                  mScenario.gridFollowingInverter->getState(),
                  mScenario.gridFollowingInverter->getStateDerivative());

    const Complex pfPower = mScenario.gridFormingPowerFlowPower;

    mLog->info("GEN2 PF operating point used by SSN_GFM: "
               "P={:.9e} W, Q={:.9e} var",
               pfPower.real(), pfPower.imag());

    mLog->info("GEN3 PCC-side PF image: P={:.9e} W, Q={:.9e} var",
               mScenario.gridFollowingPccPower.real(),
               mScenario.gridFollowingPccPower.imag());
  }

  void sample(Real time, Int step, Bool force = false,
              Bool enforceLimits = true) {
    const Bool writeCsv =
        force || step % std::max<Int>(1, mSettings.csvEverySteps) == 0;
    const Bool writeConsole =
        force || step % std::max<Int>(1, mSettings.consoleEverySteps) == 0;

    const Matrix gfmState = mScenario.gridFormingInverter->getState();
    const Matrix gfmDerivative =
        mScenario.gridFormingInverter->getStateDerivative();
    const Matrix gfmVoltage =
        mScenario.gridFormingInverter->getInterfaceVoltage();
    const Matrix gfmCurrent =
        mScenario.gridFormingInverter->getInterfaceCurrent();

    const Matrix gflState = mScenario.gridFollowingInverter->getState();
    const Matrix gflDerivative =
        mScenario.gridFollowingInverter->getStateDerivative();
    const Matrix gflVoltage =
        mScenario.gridFollowingInverter->getInterfaceVoltage();
    const Matrix gflCurrent =
        mScenario.gridFollowingInverter->getInterfaceCurrent();

    const Matrix sgVoltage =
        **mScenario.synchronousGenerator->attributeTyped<Matrix>("v_intf");
    const Matrix sgCurrent =
        **mScenario.synchronousGenerator->attributeTyped<Matrix>("i_intf");

    const Real sgOmegaPu =
        **mScenario.synchronousGenerator->attributeTyped<Real>("w_r");
    const Real sgDelta =
        **mScenario.synchronousGenerator->attributeTyped<Real>("delta");

    const Real gfmP =
        **mScenario.gridFormingInverter->attributeTyped<Real>("p_inst");
    const Real gfmQ =
        **mScenario.gridFormingInverter->attributeTyped<Real>("q_inst");
    const Real gfmOmega =
        **mScenario.gridFormingInverter->attributeTyped<Real>("omega_gfm");

    const Real gflP =
        **mScenario.gridFollowingInverter->attributeTyped<Real>("p_inst");
    const Real gflQ =
        **mScenario.gridFollowingInverter->attributeTyped<Real>("q_inst");
    const Real gflOmega =
        **mScenario.gridFollowingInverter->attributeTyped<Real>("omega_pll");

    std::array<Real, 9> busMaximumVoltagePu{};
    Real systemMaximumVoltagePu = 0.0;

    for (std::size_t busIndex = 0; busIndex < mScenario.buses.size();
         ++busIndex) {
      const Matrix busVoltage = mScenario.buses[busIndex]->voltage();
      const Real voltagePu =
          maxAbs(busVoltage) /
          peakPhaseVoltage(mScenario.busBaseVoltages[busIndex]);

      busMaximumVoltagePu[busIndex] = voltagePu;
      systemMaximumVoltagePu = std::max(systemMaximumVoltagePu, voltagePu);
    }

    const Real sgVoltagePu =
        maxAbs(sgVoltage) /
        peakPhaseVoltage(mScenario.synchronousGeneratorRatedVoltage);
    const Real sgCurrentPu =
        maxAbs(sgCurrent) /
        peakPhaseCurrent(mScenario.synchronousGeneratorRatedPower,
                         mScenario.synchronousGeneratorRatedVoltage);

    const Real gfmVoltagePu =
        maxAbs(gfmVoltage) /
        peakPhaseVoltage(mScenario.gridFormingRatedVoltage);
    const Real gfmCurrentPu =
        maxAbs(gfmCurrent) /
        peakPhaseCurrent(mScenario.gridFormingRatedPower,
                         mScenario.gridFormingRatedVoltage);

    const Real gflVoltagePu =
        maxAbs(gflVoltage) /
        peakPhaseVoltage(mScenario.gridFollowingRatedVoltage);
    const Real gflCurrentPu =
        maxAbs(gflCurrent) /
        peakPhaseCurrent(mScenario.gridFollowingRatedPower,
                         mScenario.gridFollowingRatedVoltage);

    const Real gfmOmegaPu = gfmOmega / mScenario.nominalOmega;
    const Real gflOmegaPu = gflOmega / mScenario.nominalOmega;

    const Real gfmStateNorm = gfmState.norm();
    const Real gflStateNorm = gflState.norm();

    if (!mInitialStateCaptured) {
      mInitialGfmStateNorm = std::max<Real>(1.0, gfmStateNorm);
      mInitialGflStateNorm = std::max<Real>(1.0, gflStateNorm);
      mInitialStateCaptured = true;
    }

    const Real gfmStateGrowth = gfmStateNorm / mInitialGfmStateNorm;
    const Real gflStateGrowth = gflStateNorm / mInitialGflStateNorm;

    LocalModeSummary gfmModes;
    LocalModeSummary gflModes;

    if (writeCsv || writeConsole || force) {
      // These are local fixed-input SSN modes. They are not the eigenvalues of
      // the complete coupled IEEE-9 system, but they identify whether one
      // inverter's local discretized model has already become unstable.
      gfmModes = summarizeLocalDiscreteModel(
          mScenario.gridFormingInverter->getDiscreteA(), mTimeStep);
      gflModes = summarizeLocalDiscreteModel(
          mScenario.gridFollowingInverter->getDiscreteA(), mTimeStep);
    }

    const Bool allFinite = gfmState.allFinite() && gfmDerivative.allFinite() &&
                           gfmVoltage.allFinite() && gfmCurrent.allFinite() &&
                           gflState.allFinite() && gflDerivative.allFinite() &&
                           gflVoltage.allFinite() && gflCurrent.allFinite() &&
                           sgVoltage.allFinite() && sgCurrent.allFinite() &&
                           std::isfinite(sgOmegaPu) && std::isfinite(sgDelta) &&
                           std::isfinite(gfmP) && std::isfinite(gfmQ) &&
                           std::isfinite(gfmOmegaPu) && std::isfinite(gflP) &&
                           std::isfinite(gflQ) && std::isfinite(gflOmegaPu);

    const Real maximumCurrentPu =
        std::max({sgCurrentPu, gfmCurrentPu, gflCurrentPu});

    if (writeCsv) {
      writeRow(time, step, busMaximumVoltagePu, sgOmegaPu, sgDelta, sgVoltagePu,
               sgCurrentPu, gfmP, gfmQ, gfmOmegaPu, gfmVoltagePu, gfmCurrentPu,
               gfmState, gfmDerivative, gfmStateGrowth, gfmModes, gflP, gflQ,
               gflOmegaPu, gflVoltagePu, gflCurrentPu, gflState, gflDerivative,
               gflStateGrowth, gflModes);
    }

    if (writeConsole) {
      mLog->info("t={:.6f} s, step={}: Vmax={:.4f} pu, Imax={:.4f} pu, "
                 "SG omega={:.6f} pu, GFM omega={:.6f} pu, "
                 "GFL omega={:.6f} pu, GFM |x|={:.4e}, GFL |x|={:.4e}, "
                 "rho(Ad_GFM)={:.8f}, rho(Ad_GFL)={:.8f}",
                 time, step, systemMaximumVoltagePu, maximumCurrentPu,
                 sgOmegaPu, gfmOmegaPu, gflOmegaPu, gfmStateNorm, gflStateNorm,
                 gfmModes.spectralRadius, gflModes.spectralRadius);
    }

    if (!enforceLimits)
      return;

    if (!allFinite)
      fail("A non-finite value was detected", time, step, gfmState,
           gfmDerivative, gflState, gflDerivative);

    if (systemMaximumVoltagePu > mSettings.maximumVoltagePu)
      fail("Bus voltage exceeded the diagnostic limit", time, step, gfmState,
           gfmDerivative, gflState, gflDerivative);

    if (maximumCurrentPu > mSettings.maximumCurrentPu)
      fail("Source current exceeded the diagnostic limit", time, step, gfmState,
           gfmDerivative, gflState, gflDerivative);

    if (sgOmegaPu < mSettings.minimumFrequencyPu ||
        sgOmegaPu > mSettings.maximumFrequencyPu ||
        gfmOmegaPu < mSettings.minimumFrequencyPu ||
        gfmOmegaPu > mSettings.maximumFrequencyPu ||
        gflOmegaPu < mSettings.minimumFrequencyPu ||
        gflOmegaPu > mSettings.maximumFrequencyPu)
      fail("A source frequency left the diagnostic interval", time, step,
           gfmState, gfmDerivative, gflState, gflDerivative);

    if (gfmStateGrowth > mSettings.maximumStateGrowthFactor)
      fail("GEN2 SSN_GFM state norm exceeded its growth limit", time, step,
           gfmState, gfmDerivative, gflState, gflDerivative);

    if (gflStateGrowth > mSettings.maximumStateGrowthFactor)
      fail("GEN3 SSN_GFL state norm exceeded its growth limit", time, step,
           gfmState, gfmDerivative, gflState, gflDerivative);
  }

  void reportException(Real time, Int step, const std::exception &exception) {
    mLog->critical("Simulation exception at t={:.9e} s, step={}: {}", time,
                   step, exception.what());

    try {
      sample(time, step, true, false);
    } catch (...) {
      mLog->critical(
          "The final diagnostic snapshot could not be evaluated safely.");
    }
  }

private:
  void writeHeader() {
    mCsv << "time,step";

    for (Int bus = 1; bus <= 9; ++bus)
      mCsv << ",BUS" << bus << "_vmax_pu";

    mCsv << ",GEN1_omega_pu,GEN1_delta,GEN1_vmax_pu,GEN1_imax_pu"
         << ",GEN2_p,GEN2_q,GEN2_omega_pu,GEN2_vmax_pu,GEN2_imax_pu"
         << ",GEN2_state_norm,GEN2_state_max_abs,GEN2_derivative_norm"
         << ",GEN2_derivative_max_abs,GEN2_state_growth"
         << ",GEN2_local_spectral_radius,GEN2_local_max_growth_rate"
         << ",GEN3_p,GEN3_q,GEN3_omega_pu,GEN3_vmax_pu,GEN3_imax_pu"
         << ",GEN3_state_norm,GEN3_state_max_abs,GEN3_derivative_norm"
         << ",GEN3_derivative_max_abs,GEN3_state_growth"
         << ",GEN3_local_spectral_radius,GEN3_local_max_growth_rate\n";
  }

  void writeRow(Real time, Int step,
                const std::array<Real, 9> &busMaximumVoltagePu, Real sgOmegaPu,
                Real sgDelta, Real sgVoltagePu, Real sgCurrentPu, Real gfmP,
                Real gfmQ, Real gfmOmegaPu, Real gfmVoltagePu,
                Real gfmCurrentPu, const Matrix &gfmState,
                const Matrix &gfmDerivative, Real gfmStateGrowth,
                const LocalModeSummary &gfmModes, Real gflP, Real gflQ,
                Real gflOmegaPu, Real gflVoltagePu, Real gflCurrentPu,
                const Matrix &gflState, const Matrix &gflDerivative,
                Real gflStateGrowth, const LocalModeSummary &gflModes) {
    mCsv << time << ',' << step;

    for (const Real voltagePu : busMaximumVoltagePu)
      mCsv << ',' << voltagePu;

    mCsv << ',' << sgOmegaPu << ',' << sgDelta << ',' << sgVoltagePu << ','
         << sgCurrentPu << ',' << gfmP << ',' << gfmQ << ',' << gfmOmegaPu
         << ',' << gfmVoltagePu << ',' << gfmCurrentPu << ',' << gfmState.norm()
         << ',' << maxAbs(gfmState) << ',' << gfmDerivative.norm() << ','
         << maxAbs(gfmDerivative) << ',' << gfmStateGrowth << ','
         << gfmModes.spectralRadius << ','
         << gfmModes.maximumEquivalentContinuousGrowth << ',' << gflP << ','
         << gflQ << ',' << gflOmegaPu << ',' << gflVoltagePu << ','
         << gflCurrentPu << ',' << gflState.norm() << ',' << maxAbs(gflState)
         << ',' << gflDerivative.norm() << ',' << maxAbs(gflDerivative) << ','
         << gflStateGrowth << ',' << gflModes.spectralRadius << ','
         << gflModes.maximumEquivalentContinuousGrowth << '\n';

    mCsv.flush();
  }

  void logNamedState(const String &label, const std::vector<String> &stateNames,
                     const Matrix &state, const Matrix &stateDerivative) const {
    std::ostringstream stream;

    stream << '\n' << label << " initial states and derivatives:";

    for (Eigen::Index index = 0; index < state.rows(); ++index) {
      const String stateName =
          index < static_cast<Eigen::Index>(stateNames.size())
              ? stateNames[static_cast<std::size_t>(index)]
              : "x" + std::to_string(index);

      stream << "\n  [" << index << "] " << stateName << " = "
             << std::setprecision(12) << state(index, 0)
             << ", derivative = " << stateDerivative(index, 0);
    }

    stream << "\n  state norm = " << state.norm()
           << "\n  derivative norm = " << stateDerivative.norm();

    mLog->info("{}", stream.str());
  }

  void fail(const String &reason, Real time, Int step, const Matrix &gfmState,
            const Matrix &gfmDerivative, const Matrix &gflState,
            const Matrix &gflDerivative) const {
    mLog->critical("{} at t={:.9e} s, step={}."
                   "\nGEN2 state: {}"
                   "\nGEN2 derivative: {}"
                   "\nGEN3 state: {}"
                   "\nGEN3 derivative: {}",
                   reason, time, step, CPS::Logger::matrixToString(gfmState),
                   CPS::Logger::matrixToString(gfmDerivative),
                   CPS::Logger::matrixToString(gflState),
                   CPS::Logger::matrixToString(gflDerivative));

    if (mSettings.abortOnLimit)
      throw std::runtime_error(reason);
  }

private:
  const ScenarioArtifacts &mScenario;
  DiagnosticSettings mSettings;
  Real mTimeStep;

  CPS::Logger::Log mLog;
  std::ofstream mCsv;

  Bool mInitialStateCaptured = false;
  Real mInitialGfmStateNorm = 1.0;
  Real mInitialGflStateNorm = 1.0;
};

} // namespace

ScenarioArtifacts buildTopology(CommandLineArgs &args,
                                std::shared_ptr<DataLoggerInterface> logger) {

  String simName = args.name;

  CPS::CIM::Examples::Grids::IEEE9::ScenarioConfig ieee9(args.sysFreq);

  // POWER FLOW FOR INITIALIZATION
  CPS::Logger::get(args.name)->info("Creating power flow initialization.");

  String simNamePF = simName + "_PF";
  CPS::Logger::setLogDir("logs/" + simNamePF);

  // Nodes
  auto n1PF = SimNode<Complex>::make("BUS1", PhaseType::Single);
  auto n2PF = SimNode<Complex>::make("BUS2", PhaseType::Single);
  auto n3PF = SimNode<Complex>::make("BUS3", PhaseType::Single);
  auto n4PF = SimNode<Complex>::make("BUS4", PhaseType::Single);
  auto n5PF = SimNode<Complex>::make("BUS5", PhaseType::Single);
  auto n6PF = SimNode<Complex>::make("BUS6", PhaseType::Single);
  auto n7PF = SimNode<Complex>::make("BUS7", PhaseType::Single);
  auto n8PF = SimNode<Complex>::make("BUS8", PhaseType::Single);
  auto n9PF = SimNode<Complex>::make("BUS9", PhaseType::Single);

  auto gen1PF = SP::Ph1::SynchronGenerator::make(ieee9.gen1.Name,
                                                 CPS::Logger::Level::off);
  gen1PF->setParameters(ieee9.gen1.RatedPower, ieee9.gen1.RatedVoltage,
                        ieee9.gen1.InitialPower, ieee9.gen1.InitialVoltage,
                        ieee9.gen1.BusType);
  gen1PF->setBaseVoltage(ieee9.gen1.RatedVoltage);

  auto gen2PF = SP::Ph1::SynchronGenerator::make(ieee9.gen2.Name,
                                                 CPS::Logger::Level::off);
  gen2PF->setParameters(ieee9.gen2.RatedPower, ieee9.gen2.RatedVoltage,
                        ieee9.gen2.InitialPower, ieee9.gen2.InitialVoltage,
                        ieee9.gen2.BusType);
  gen2PF->setBaseVoltage(ieee9.gen2.RatedVoltage);

  // gen3's PF image: a negative PQ load injecting the PCC-side power (the
  // rc-corrected filter reference).
  const GflParams gfl;
  const auto [gfl3PPcc, gfl3QPcc] = Math::pccPowerFromFilterPowerReference(
      ieee9.gen3.InitialPower, ieee9.gen3.InitialPowerReactive, gfl.Rc,
      ieee9.gen3.RatedVoltage);
  auto gfl3PF = SP::Ph1::Load::make(ieee9.gen3.Name, CPS::Logger::Level::off);
  gfl3PF->setParameters(-gfl3PPcc, -gfl3QPcc, ieee9.gen3.RatedVoltage);
  gfl3PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // Loads
  auto load5PF = SP::Ph1::Load::make(ieee9.load5.Name, CPS::Logger::Level::off);
  load5PF->setParameters(ieee9.load5.RealPower, ieee9.load5.ReactivePower,
                         ieee9.load5.BaseVoltage);
  load5PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  auto load6PF = SP::Ph1::Load::make(ieee9.load6.Name, CPS::Logger::Level::off);
  load6PF->setParameters(ieee9.load6.RealPower, ieee9.load6.ReactivePower,
                         ieee9.load6.BaseVoltage);
  load6PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  auto load8PF = SP::Ph1::Load::make(ieee9.load8.Name, CPS::Logger::Level::off);
  load8PF->setParameters(ieee9.load8.RealPower, ieee9.load8.ReactivePower,
                         ieee9.load8.BaseVoltage);
  load8PF->modifyPowerFlowBusType(PowerflowBusType::PQ);

  // Transmission Lines

  auto line54PF =
      SP::Ph1::PiLine::make(ieee9.line54.Name, CPS::Logger::Level::off);
  line54PF->setParameters(ieee9.line54.Resistance, ieee9.line54.Inductance,
                          ieee9.line54.Capacitance, ieee9.line54.Conductance);
  line54PF->setBaseVoltage(ieee9.line54.BaseVoltage);

  auto line64PF =
      SP::Ph1::PiLine::make(ieee9.line64.Name, CPS::Logger::Level::off);
  line64PF->setParameters(ieee9.line64.Resistance, ieee9.line64.Inductance,
                          ieee9.line64.Capacitance, ieee9.line64.Conductance);
  line64PF->setBaseVoltage(ieee9.line64.BaseVoltage);

  auto line75PF =
      SP::Ph1::PiLine::make(ieee9.line75.Name, CPS::Logger::Level::off);
  line75PF->setParameters(ieee9.line75.Resistance, ieee9.line75.Inductance,
                          ieee9.line75.Capacitance, ieee9.line75.Conductance);
  line75PF->setBaseVoltage(ieee9.line75.BaseVoltage);

  auto line96PF =
      SP::Ph1::PiLine::make(ieee9.line96.Name, CPS::Logger::Level::off);
  line96PF->setParameters(ieee9.line96.Resistance, ieee9.line96.Inductance,
                          ieee9.line96.Capacitance, ieee9.line96.Conductance);
  line96PF->setBaseVoltage(ieee9.line96.BaseVoltage);

  auto line78PF =
      SP::Ph1::PiLine::make(ieee9.line78.Name, CPS::Logger::Level::off);
  line78PF->setParameters(ieee9.line78.Resistance, ieee9.line78.Inductance,
                          ieee9.line78.Capacitance, ieee9.line78.Conductance);
  line78PF->setBaseVoltage(ieee9.line78.BaseVoltage);

  auto line89PF =
      SP::Ph1::PiLine::make(ieee9.line89.Name, CPS::Logger::Level::off);
  line89PF->setParameters(ieee9.line89.Resistance, ieee9.line89.Inductance,
                          ieee9.line89.Capacitance, ieee9.line89.Conductance);
  line89PF->setBaseVoltage(ieee9.line89.BaseVoltage);

  // Transformers

  auto transf14PF =
      SP::Ph1::Transformer::make(ieee9.transf14.Name, CPS::Logger::Level::off);
  transf14PF->setParameters(
      ieee9.transf14.VoltageLVSide, ieee9.transf14.VoltageHVSide,
      ieee9.transf14.Ratio, 0.0, // No phase shift (ratioPhase = 0.0)
      ieee9.transf14.Resistance, ieee9.transf14.Inductance);
  transf14PF->setBaseVoltage(ieee9.transf14.VoltageHVSide);

  auto transf27PF =
      SP::Ph1::Transformer::make(ieee9.transf27.Name, CPS::Logger::Level::off);
  transf27PF->setParameters(ieee9.transf27.VoltageLVSide,
                            ieee9.transf27.VoltageHVSide, ieee9.transf27.Ratio,
                            0.0, ieee9.transf27.Resistance,
                            ieee9.transf27.Inductance);
  transf27PF->setBaseVoltage(ieee9.transf27.VoltageHVSide);

  auto transf39PF =
      SP::Ph1::Transformer::make(ieee9.transf39.Name, CPS::Logger::Level::off);
  transf39PF->setParameters(ieee9.transf39.VoltageLVSide,
                            ieee9.transf39.VoltageHVSide, ieee9.transf39.Ratio,
                            0.0, ieee9.transf39.Resistance,
                            ieee9.transf39.Inductance);
  transf39PF->setBaseVoltage(ieee9.transf39.VoltageHVSide);

  // Connect components
  gen1PF->connect({n1PF});
  gen2PF->connect({n2PF});
  gfl3PF->connect({n3PF});

  load5PF->connect({n5PF});
  load6PF->connect({n6PF});
  load8PF->connect({n8PF});

  line54PF->connect({n5PF, n4PF});
  line64PF->connect({n6PF, n4PF});
  line75PF->connect({n7PF, n5PF});
  line96PF->connect({n9PF, n6PF});
  line78PF->connect({n7PF, n8PF});
  line89PF->connect({n8PF, n9PF});

  transf14PF->connect({n1PF, n4PF});
  transf27PF->connect({n2PF, n7PF});
  transf39PF->connect({n3PF, n9PF});

  // Create system topology
  auto systemPF = SystemTopology(
      ieee9.nomFreq,
      SystemNodeList{n1PF, n2PF, n3PF, n4PF, n5PF, n6PF, n7PF, n8PF, n9PF},
      SystemComponentList{gen1PF, gen2PF, gfl3PF, load5PF, load6PF, load8PF,
                          line54PF, line64PF, line75PF, line96PF, line78PF,
                          line89PF, transf14PF, transf27PF, transf39PF});

  // Logger
  auto loggerPF = DataLogger::make(simNamePF, CPS::Logger::Level::off);
  // Log node voltages
  loggerPF->logAttribute("v_bus1", n1PF->attribute("v"));
  loggerPF->logAttribute("v_bus2", n2PF->attribute("v"));
  loggerPF->logAttribute("v_bus3", n3PF->attribute("v"));
  loggerPF->logAttribute("v_bus4", n4PF->attribute("v"));
  loggerPF->logAttribute("v_bus5", n5PF->attribute("v"));
  loggerPF->logAttribute("v_bus6", n6PF->attribute("v"));
  loggerPF->logAttribute("v_bus7", n7PF->attribute("v"));
  loggerPF->logAttribute("v_bus8", n8PF->attribute("v"));
  loggerPF->logAttribute("v_bus9", n9PF->attribute("v"));
  // Log node powers
  loggerPF->logAttribute("s_bus1", n1PF->attribute("s"));
  loggerPF->logAttribute("s_bus2", n2PF->attribute("s"));
  loggerPF->logAttribute("s_bus3", n3PF->attribute("s"));
  loggerPF->logAttribute("s_bus4", n4PF->attribute("s"));
  loggerPF->logAttribute("s_bus5", n5PF->attribute("s"));
  loggerPF->logAttribute("s_bus6", n6PF->attribute("s"));
  loggerPF->logAttribute("s_bus7", n7PF->attribute("s"));
  loggerPF->logAttribute("s_bus8", n8PF->attribute("s"));
  loggerPF->logAttribute("s_bus9", n9PF->attribute("s"));

  // Run power flow simulation
  Simulation simPF(simNamePF, CPS::Logger::Level::off);
  simPF.setSystem(systemPF);
  simPF.setTimeStep(args.timeStep);
  simPF.setFinalTime(1 * args.timeStep);
  simPF.setDomain(Domain::SP);
  simPF.setSolverType(Solver::Type::NRP);
  simPF.setSolverAndComponentBehaviour(Solver::Behaviour::Initialization);
  simPF.doInitFromNodesAndTerminals(false);
  simPF.addLogger(loggerPF);
  simPF.run();
  CPS::Logger::get(args.name)->info("Power flow simulation finished.");

  // Use the solved PV reactive power for GEN2. Initializing the GFM with the
  // configured Q value while the PF solved another Q creates a controller
  // disequilibrium at t=0.
  const Complex gen2PowerPF = gen2PF->getApparentPower();

  CPS::Logger::get(args.name)->info(
      "GEN2 solved PF injection: P={:.9e} W, Q={:.9e} var; "
      "configured values: P={:.9e} W, Q={:.9e} var",
      gen2PowerPF.real(), gen2PowerPF.imag(), ieee9.gen2.InitialPower,
      ieee9.gen2.InitialPowerReactive);

  // Seed voltages for the step-2 inverters: converged n2 (GFM PCC) and
  // n3 (GFL PCC) node voltages, magnitude [kV] and angle [deg]. Dedicated
  // console-enabled logger so the seeds survive the setLogDir churn.
  Complex v2PF = n2PF->singleVoltage();
  Complex v3PF = n3PF->singleVoltage();
  auto seedLog = CPS::Logger::get(simNamePF + "_seed", CPS::Logger::Level::info,
                                  CPS::Logger::Level::info);
  seedLog->info(
      "PF seed n2 (BUS2): |V| = {:.4f} kV ({:.4f} pu), angle = {:.4f} deg",
      std::abs(v2PF) / 1e3, std::abs(v2PF) / ieee9.gen2.RatedVoltage,
      std::arg(v2PF) * 180.0 / PI);
  seedLog->info(
      "PF seed n3 (BUS3): |V| = {:.4f} kV ({:.4f} pu), angle = {:.4f} deg",
      std::abs(v3PF) / 1e3, std::abs(v3PF) / ieee9.gen3.RatedVoltage,
      std::arg(v3PF) * 180.0 / PI);

  // DYNAMIC SIMULATION - EMT
  CPS::Logger::get(args.name)->info("Dynamic simulation initialization.");
  String simNameEMT = simName + "_EMT";
  CPS::Logger::setLogDir("logs/" + simNameEMT);

  // Nodes
  auto n1EMT = SimNode<Real>::make("BUS1", PhaseType::ABC);
  auto n2EMT = SimNode<Real>::make("BUS2", PhaseType::ABC);
  auto n3EMT = SimNode<Real>::make("BUS3", PhaseType::ABC);
  auto n4EMT = SimNode<Real>::make("BUS4", PhaseType::ABC);
  auto n5EMT = SimNode<Real>::make("BUS5", PhaseType::ABC);
  auto n6EMT = SimNode<Real>::make("BUS6", PhaseType::ABC);
  auto n7EMT = SimNode<Real>::make("BUS7", PhaseType::ABC);
  auto n8EMT = SimNode<Real>::make("BUS8", PhaseType::ABC);
  auto n9EMT = SimNode<Real>::make("BUS9", PhaseType::ABC);

  // Generators
  const Bool componentInfoLog =
      args.options.find("component_log") != args.options.end() &&
      args.getOptionBool("component_log");
  const auto componentLogLevel =
      componentInfoLog ? CPS::Logger::Level::info : CPS::Logger::Level::off;

  auto gen1EMT = EMT::Ph3::SynchronGenerator4OrderVBR::make(ieee9.gen1.Name,
                                                            componentLogLevel);

  gen1EMT->setOperationalParametersPerUnit(
      ieee9.gen1.RatedPower,   // nomPower [VA]
      ieee9.gen1.RatedVoltage, // nomVolt [V]
      ieee9.nomFreq,           // nomFreq [Hz]
      ieee9.gen1.H, ieee9.gen1.Xd, ieee9.gen1.Xq, ieee9.gen1.Xa,
      ieee9.gen1.XdPrime, ieee9.gen1.XqPrime, ieee9.gen1.TdoPrime,
      ieee9.gen1.TqoPrime);

  auto exciter1Params = std::make_shared<Signal::ExciterDC1SimpParameters>();
  exciter1Params->Ta = ieee9.exc1.TA;
  exciter1Params->Ka = ieee9.exc1.KA;
  exciter1Params->Tef = ieee9.exc1.TE;
  exciter1Params->Kef = ieee9.exc1.KE;
  exciter1Params->Tf = ieee9.exc1.TF;
  exciter1Params->Kf = ieee9.exc1.KF;
  exciter1Params->Tr = 0.01;
  exciter1Params->MaxVa = ieee9.exc1.VRmax;
  exciter1Params->MinVa = ieee9.exc1.VRmin;
  exciter1Params->Bef = std::log(ieee9.exc1.S_EX2 / ieee9.exc1.S_EX1) /
                        (ieee9.exc1.EX2 - ieee9.exc1.EX1);
  exciter1Params->Aef =
      ieee9.exc1.S_EX1 / std::exp(exciter1Params->Bef * ieee9.exc1.EX1);
  auto exciter1 =
      Signal::ExciterDC1Simp::make("Gen1_Exciter", CPS::Logger::Level::off);
  exciter1->setParameters(exciter1Params);
  gen1EMT->addExciter(exciter1);

  // Adaptation of the governor model parameters to the dpsim implementation
  CPS::Real T4 = 1.0;
  CPS::Real T5 = 1.0;

  std::shared_ptr<Signal::TurbineGovernorType1> turbineGovernor1 =
      Signal::TurbineGovernorType1::make("Gen1_TurbineGovernor",
                                         CPS::Logger::Level::off);

  turbineGovernor1->setParameters(ieee9.gov1.T2, T4, T5, ieee9.gov1.T3,
                                  ieee9.gov1.T1, ieee9.gov1.R, ieee9.gov1.Vmin,
                                  ieee9.gov1.Vmax, 1.0);

  gen1EMT->addGovernor(turbineGovernor1);

  const Real omegaN = 2.0 * PI * ieee9.nomFreq;

  // GEN2 is replaced by the grid-forming SSN inverter while retaining the
  // component identity required by initWithPowerflow().
  //
  // By default, use the solved BUS2 voltage and solved GEN2 P/Q. This removes
  // the initial mismatch between a PV power-flow operating point and a GFM
  // initialized from an unrelated reactive-power reference.
  const Bool useSolvedGfmOperatingPoint =
      args.options.find("gfm_use_pf_operating_point") == args.options.end()
          ? true
          : args.getOptionBool("gfm_use_pf_operating_point");

  const Real gfmNominalVoltage =
      RMS3PH_TO_PEAK1PH *
      (useSolvedGfmOperatingPoint ? std::abs(v2PF) : ieee9.gen2.InitialVoltage);

  const Real gfmPReference =
      useSolvedGfmOperatingPoint ? gen2PowerPF.real() : ieee9.gen2.InitialPower;
  const Real gfmQReference = useSolvedGfmOperatingPoint
                                 ? gen2PowerPF.imag()
                                 : ieee9.gen2.InitialPowerReactive;

  GfmParams gfm;
  // Optional overrides for studying the grid-forming tuning from a notebook.
  auto opt = [&](const String &key, Real def) {
    return args.options.find(key) != args.options.end()
               ? args.getOptionReal(key)
               : def;
  };
  gfm.dampingCoefficient = opt("gfm_d", gfm.dampingCoefficient);
  gfm.KpVoltage = opt("gfm_kpv", gfm.KpVoltage);
  gfm.KiVoltage = opt("gfm_kiv", gfm.KiVoltage);
  gfm.reactivePowerDroop = opt("gfm_dq", gfm.reactivePowerDroop);
  gfm.reactiveDroopCutoff = opt("gfm_dqc", gfm.reactiveDroopCutoff);
  const Real gfmFeedforward = opt("gfm_ff", gfm.gridCurrentFeedforward);
  const Real gfmVirtualResistance = opt("gfm_rv", 0.0);

  auto gen2EMT = EMT::Ph3::SSN_GFM::make(ieee9.gen2.Name, ieee9.gen2.Name,
                                         componentLogLevel);
  gen2EMT->setNumericalLinearizationParameters(1e-6, 1e-8);
  gen2EMT->setLinearizationUpdateInterval(1);
  gen2EMT->setParameters(
      gfm.Lf, gfm.Cf, gfm.Rf, gfm.Rc, gfmNominalVoltage, omegaN, gfmPReference,
      gfmQReference, gfm.virtualInertia, gfm.dampingCoefficient,
      gfm.voltageDroopGain, gfm.reactiveIntegralGain, gfm.KpVoltage,
      gfm.KiVoltage, gfm.KpCurrent, gfm.KiCurrent, gfm.activeDampingGain,
      gfm.powerFilterCutoff, gfm.delayBandwidth);
  // Grid-connected control: no grid-current feedforward, proportional Q-V droop.
  gen2EMT->setGridCurrentFeedforward(gfmFeedforward);
  gen2EMT->setVirtualImpedance(gfmVirtualResistance, 0.0);
  gen2EMT->setReactivePowerDroop(gfm.reactivePowerDroop,
                                 gfm.reactiveDroopCutoff);

  // gen3 replaced by a grid-following averaged VSI (SSN), keeping the GEN3
  // identity so the topology wiring is unchanged.
  auto gen3EMT = EMT::Ph3::SSN_GFL::make(ieee9.gen3.Name, ieee9.gen3.Name,
                                         componentLogLevel);
  // Optional overrides for studying the grid-following tuning from a notebook.
  gen3EMT->setParameters(
      gfl.Lf, gfl.Cf, gfl.Rf, gfl.Rc, omegaN, opt("gfl_kppll", gfl.KpPLL),
      opt("gfl_kipll", gfl.KiPLL), omegaN, ieee9.gen3.InitialPower,
      ieee9.gen3.InitialPowerReactive, opt("gfl_kpp", gfl.KpPowerCtrl),
      opt("gfl_kip", gfl.KiPowerCtrl), opt("gfl_kpi", gfl.KpCurrCtrl),
      opt("gfl_kii", gfl.KiCurrCtrl));

  // Loads
  auto load5EMT =
      EMT::Ph3::RXLoad::make(ieee9.load5.Name, CPS::Logger::Level::off);
  load5EMT->setParameters(
      Math::singlePhasePowerToThreePhase(ieee9.load5.RealPower),
      Math::singlePhasePowerToThreePhase(ieee9.load5.ReactivePower),
      ieee9.load5.BaseVoltage);

  auto load6EMT =
      EMT::Ph3::RXLoad::make(ieee9.load6.Name, CPS::Logger::Level::off);
  load6EMT->setParameters(
      Math::singlePhasePowerToThreePhase(ieee9.load6.RealPower),
      Math::singlePhasePowerToThreePhase(ieee9.load6.ReactivePower),
      ieee9.load6.BaseVoltage);

  auto load8EMT =
      EMT::Ph3::RXLoad::make(ieee9.load8.Name, CPS::Logger::Level::off);
  load8EMT->setParameters(
      Math::singlePhasePowerToThreePhase(ieee9.load8.RealPower),
      Math::singlePhasePowerToThreePhase(ieee9.load8.ReactivePower),
      ieee9.load8.BaseVoltage);

  // Lines
  auto line54EMT =
      EMT::Ph3::PiLine::make(ieee9.line54.Name, CPS::Logger::Level::off);
  line54EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line54.Conductance));

  auto line64EMT =
      EMT::Ph3::PiLine::make(ieee9.line64.Name, CPS::Logger::Level::off);
  line64EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line64.Conductance));

  auto line75EMT =
      EMT::Ph3::PiLine::make(ieee9.line75.Name, CPS::Logger::Level::off);
  line75EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line75.Conductance));

  auto line96EMT =
      EMT::Ph3::PiLine::make(ieee9.line96.Name, CPS::Logger::Level::off);
  line96EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line96.Conductance));

  auto line78EMT =
      EMT::Ph3::PiLine::make(ieee9.line78.Name, CPS::Logger::Level::off);
  line78EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line78.Conductance));

  auto line89EMT =
      EMT::Ph3::PiLine::make(ieee9.line89.Name, CPS::Logger::Level::off);
  line89EMT->setParameters(
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Inductance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Capacitance),
      Math::singlePhaseParameterToThreePhase(ieee9.line89.Conductance));

  // Transformers
  auto transf14EMT =
      EMT::Ph3::Transformer::make(ieee9.transf14.Name, CPS::Logger::Level::off);
  transf14EMT->setParameters(
      ieee9.transf14.VoltageLVSide, ieee9.transf14.VoltageHVSide,
      ieee9.transf14.RatedPower, ieee9.transf14.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf14.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf14.Inductance));

  auto transf27EMT =
      EMT::Ph3::Transformer::make(ieee9.transf27.Name, CPS::Logger::Level::off);
  transf27EMT->setParameters(
      ieee9.transf27.VoltageLVSide, ieee9.transf27.VoltageHVSide,
      ieee9.transf27.RatedPower, ieee9.transf27.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf27.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf27.Inductance));

  auto transf39EMT =
      EMT::Ph3::Transformer::make(ieee9.transf39.Name, CPS::Logger::Level::off);
  transf39EMT->setParameters(
      ieee9.transf39.VoltageLVSide, ieee9.transf39.VoltageHVSide,
      ieee9.transf39.RatedPower, ieee9.transf39.Ratio, 0.0,
      Math::singlePhaseParameterToThreePhase(ieee9.transf39.Resistance),
      Math::singlePhaseParameterToThreePhase(ieee9.transf39.Inductance));

  // Connect components to nodes
  gen1EMT->connect({n1EMT});
  // Inverter terminals: 0 = GND, 1 = PCC.
  gen2EMT->connect({SimNode<Real>::GND, n2EMT});
  gen3EMT->connect({SimNode<Real>::GND, n3EMT});

  load5EMT->connect({n5EMT});
  load6EMT->connect({n6EMT});
  load8EMT->connect({n8EMT});

  line54EMT->connect({n5EMT, n4EMT});
  line64EMT->connect({n6EMT, n4EMT});
  line75EMT->connect({n7EMT, n5EMT});
  line96EMT->connect({n9EMT, n6EMT});
  line78EMT->connect({n7EMT, n8EMT});
  line89EMT->connect({n8EMT, n9EMT});

  transf14EMT->connect({n1EMT, n4EMT});
  transf27EMT->connect({n2EMT, n7EMT});
  transf39EMT->connect({n3EMT, n9EMT});

  // Create system topology
  auto systemEMT = SystemTopology(
      ieee9.nomFreq,
      SystemNodeList{n1EMT, n2EMT, n3EMT, n4EMT, n5EMT, n6EMT, n7EMT, n8EMT,
                     n9EMT},
      SystemComponentList{gen1EMT, gen2EMT, gen3EMT, load5EMT, load6EMT,
                          load8EMT, line54EMT, line64EMT, line75EMT, line96EMT,
                          line78EMT, line89EMT, transf14EMT, transf27EMT,
                          transf39EMT});

  systemEMT.initWithPowerflow(systemPF, Domain::EMT);

  // Logger
  if (logger) {
    // Logging
    logger->logAttribute("BUS1", n1EMT->attribute("v"));
    logger->logAttribute("BUS2", n2EMT->attribute("v"));
    logger->logAttribute("BUS3", n3EMT->attribute("v"));
    logger->logAttribute("BUS4", n4EMT->attribute("v"));
    logger->logAttribute("BUS5", n5EMT->attribute("v"));
    logger->logAttribute("BUS6", n6EMT->attribute("v"));
    logger->logAttribute("BUS7", n7EMT->attribute("v"));
    logger->logAttribute("BUS8", n8EMT->attribute("v"));
    logger->logAttribute("BUS9", n9EMT->attribute("v"));

    // GFM inverter (gen2) signals
    logger->logAttribute("GEN2.I", gen2EMT->attribute("i_intf"));
    logger->logAttribute("GEN2.V", gen2EMT->attribute("v_intf"));
    logger->logAttribute("GEN2.p_inst", gen2EMT->attribute("p_inst"));
    logger->logAttribute("GEN2.q_inst", gen2EMT->attribute("q_inst"));
    logger->logAttribute("GEN2.omega", gen2EMT->attribute("omega_gfm"));
    logger->logAttribute("GEN2.vc_d", gen2EMT->attribute("vc_d"));
    logger->logAttribute("GEN2.vc_q", gen2EMT->attribute("vc_q"));
    logger->logAttribute("GEN2.theta", gen2EMT->attribute("theta_gfm"));
    logger->logAttribute("GEN2.voltage_magnitude",
                         gen2EMT->attribute("voltage_magnitude_gfm"));
    logger->logAttribute("GEN2.i_grid_d", gen2EMT->attribute("i_grid_d"));
    logger->logAttribute("GEN2.i_grid_q", gen2EMT->attribute("i_grid_q"));
    logger->logAttribute("GEN2.if_d", gen2EMT->attribute("if_d"));
    logger->logAttribute("GEN2.if_q", gen2EMT->attribute("if_q"));
    logger->logAttribute("GEN2.v_ref_d", gen2EMT->attribute("v_ref_d"));
    logger->logAttribute("GEN2.v_ref_q", gen2EMT->attribute("v_ref_q"));
    logger->logAttribute("GEN2.x", gen2EMT->attribute("x"));

    // GFL inverter (gen3) signals
    logger->logAttribute("GEN3.I", gen3EMT->attribute("i_intf"));
    logger->logAttribute("GEN3.V", gen3EMT->attribute("v_intf"));
    logger->logAttribute("GEN3.p_inst", gen3EMT->attribute("p_inst"));
    logger->logAttribute("GEN3.q_inst", gen3EMT->attribute("q_inst"));
    logger->logAttribute("GEN3.omega_pll", gen3EMT->attribute("omega_pll"));
    logger->logAttribute("GEN3.vc_d", gen3EMT->attribute("vc_d"));
    logger->logAttribute("GEN3.vc_q", gen3EMT->attribute("vc_q"));
    logger->logAttribute("GEN3.irc_d", gen3EMT->attribute("irc_d"));
    logger->logAttribute("GEN3.irc_q", gen3EMT->attribute("irc_q"));
    logger->logAttribute("GEN3.x", gen3EMT->attribute("x"));

    // Synchronous-generator internal quantities.
    logger->logAttribute("GEN1.Vdq0", gen1EMT->attribute("Vdq0"));
    logger->logAttribute("GEN1.Idq0", gen1EMT->attribute("Idq0"));
    logger->logAttribute("GEN1.Te", gen1EMT->attribute("Te"));
    logger->logAttribute("GEN1.Tm", gen1EMT->attribute("Tm"));
    logger->logAttribute("GEN1.Ef", gen1EMT->attribute("Ef"));
    logger->logAttribute("GEN1.Edq0_t", gen1EMT->attribute("Edq0_t"));

    // Loads are useful for distinguishing source divergence from a network
    // initialization problem.
    logger->logAttribute("LOAD5.I", load5EMT->attribute("i_intf"));
    logger->logAttribute("LOAD5.V", load5EMT->attribute("v_intf"));
    logger->logAttribute("LOAD6.I", load6EMT->attribute("i_intf"));
    logger->logAttribute("LOAD6.V", load6EMT->attribute("v_intf"));
    logger->logAttribute("LOAD8.I", load8EMT->attribute("i_intf"));
    logger->logAttribute("LOAD8.V", load8EMT->attribute("v_intf"));

    // log generator's current
    for (auto comp : systemEMT.mComponents) {
      if (std::dynamic_pointer_cast<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>(
              comp)) {
        logger->logAttribute(comp->name() + ".I", comp->attribute("i_intf"));
        logger->logAttribute(comp->name() + ".V", comp->attribute("v_intf"));
        logger->logAttribute(comp->name() + ".omega", comp->attribute("w_r"));
        logger->logAttribute(comp->name() + ".delta", comp->attribute("delta"));
      }
    }

    // log transfomers voltages & currents
    for (auto comp : systemEMT.mComponents) {
      if (std::dynamic_pointer_cast<CPS::EMT::Ph3::Transformer>(comp)) {
        logger->logAttribute(comp->name() + ".I", comp->attribute("i_intf"));
        logger->logAttribute(comp->name() + ".V", comp->attribute("v_intf"));
      }
    }

    // log Lines voltages & currents
    for (auto comp : systemEMT.mComponents) {
      if (std::dynamic_pointer_cast<CPS::EMT::Ph3::PiLine>(comp)) {
        logger->logAttribute(comp->name() + ".I", comp->attribute("i_intf"));
        logger->logAttribute(comp->name() + ".V", comp->attribute("v_intf"));
      }
    }
  }

  const std::array<Real, 9> busBaseVoltages = {
      ieee9.gen1.RatedVoltage,      ieee9.gen2.RatedVoltage,
      ieee9.gen3.RatedVoltage,      ieee9.transf14.VoltageHVSide,
      ieee9.load5.BaseVoltage,      ieee9.load6.BaseVoltage,
      ieee9.transf27.VoltageHVSide, ieee9.load8.BaseVoltage,
      ieee9.transf39.VoltageHVSide,
  };

  return {
      systemEMT,
      {n1EMT, n2EMT, n3EMT, n4EMT, n5EMT, n6EMT, n7EMT, n8EMT, n9EMT},
      gen1EMT,
      gen2EMT,
      gen3EMT,
      busBaseVoltages,
      omegaN,
      ieee9.gen1.RatedPower,
      ieee9.gen1.RatedVoltage,
      ieee9.gen2.RatedPower,
      ieee9.gen2.RatedVoltage,
      ieee9.gen3.RatedPower,
      ieee9.gen3.RatedVoltage,
      gen2PowerPF,
      Complex(gfl3PPcc, gfl3QPcc),
  };
}

int main(int argc, char *argv[]) {

  CommandLineArgs args(argc, argv, "EMT_Ph3_IEEE9_SSN_newGFM",
                       100e-6, // EMT time step
                       5.0,    // duration
                       60.0,   // IEEE-9 frequency
                       -1, CPS::Logger::Level::info, CPS::Logger::Level::off,
                       false, false, false, CPS::Domain::EMT);

  CPS::Logger::setLogDir("./logs/" + args.name);

  // Detailed signal logging is enabled by default for this diagnostic case.
  const Bool log = args.options.find("log") == args.options.end()
                       ? true
                       : args.getOptionBool("log");

  const auto optionInt = [&](const String &name, Int defaultValue) {
    return args.options.find(name) == args.options.end()
               ? defaultValue
               : args.getOptionInt(name);
  };

  const auto optionReal = [&](const String &name, Real defaultValue) {
    return args.options.find(name) == args.options.end()
               ? defaultValue
               : args.getOptionReal(name);
  };

  const auto optionBool = [&](const String &name, Bool defaultValue) {
    return args.options.find(name) == args.options.end()
               ? defaultValue
               : args.getOptionBool(name);
  };

  DiagnosticSettings diagnosticSettings;
  diagnosticSettings.csvEverySteps =
      optionInt("diag_every", diagnosticSettings.csvEverySteps);
  diagnosticSettings.consoleEverySteps =
      optionInt("diag_console_every", diagnosticSettings.consoleEverySteps);
  diagnosticSettings.maximumVoltagePu =
      optionReal("diag_vmax_pu", diagnosticSettings.maximumVoltagePu);
  diagnosticSettings.maximumCurrentPu =
      optionReal("diag_imax_pu", diagnosticSettings.maximumCurrentPu);
  diagnosticSettings.minimumFrequencyPu =
      optionReal("diag_fmin_pu", diagnosticSettings.minimumFrequencyPu);
  diagnosticSettings.maximumFrequencyPu =
      optionReal("diag_fmax_pu", diagnosticSettings.maximumFrequencyPu);
  diagnosticSettings.maximumStateGrowthFactor = optionReal(
      "diag_state_growth", diagnosticSettings.maximumStateGrowthFactor);
  diagnosticSettings.abortOnLimit =
      optionBool("diag_abort", diagnosticSettings.abortOnLimit);

  std::shared_ptr<DataLoggerInterface> logger = nullptr;

  if (log) {
    logger = DataLogger::make(args.name, CPS::Logger::Level::off);
  }

  auto scenario = buildTopology(args, logger);

  // buildTopology() uses separate PF/EMT log directories while constructing
  // the systems. Restore the main diagnostic directory before creating the
  // runtime loggers.

  RuntimeDiagnostics diagnostics(args.name, scenario, diagnosticSettings,
                                 args.timeStep);

  Simulation sim(args.name, args);
  sim.setSystem(scenario.system);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);

  // Global MNA state-space extraction is deliberately not enabled here:
  // SynchronGenerator4OrderVBR is not currently an extraction contributor.
  // The diagnostics instead report the live local discrete SSN modes of both
  // inverters together with all network/source trajectories.
  sim.doStateSpaceExtraction(false);

  const Bool collectStepTimes = optionBool("step_times", false);
  sim.setLogStepTimes(collectStepTimes);

  if (log)
    sim.addLogger(logger);

  Bool simulationStarted = false;

  try {
    sim.start();
    simulationStarted = true;

    // start() logs t=0 and then advances the internal simulation time to dt.
    // The component values are nevertheless still the initialized t=0 values.
    diagnostics.logInitialStateDetails();
    diagnostics.sample(0.0, 0, true);

    while (sim.time() < sim.finalTime() + DOUBLE_EPSILON) {
      sim.step();
      diagnostics.sample(sim.time(), sim.timeStepCount());
    }

    sim.stop();
    simulationStarted = false;

  } catch (const std::exception &exception) {
    diagnostics.reportException(sim.time(), sim.timeStepCount(), exception);

    if (simulationStarted) {
      try {
        sim.stop();
      } catch (...) {
        CPS::Logger::get(args.name)->critical(
            "Simulation cleanup failed after the diagnostic exception.");
      }
    }

    CPS::Logger::get(args.name)->critical(
        "Diagnostic IEEE-9 simulation aborted: {}", exception.what());

    return 1;
  }

  if (collectStepTimes)
    sim.logStepTimes(args.name + "_step_times");

  CPS::Logger::get(args.name)->info(
      "Simulation finished without crossing a diagnostic abort limit.");

  return 0;
}
