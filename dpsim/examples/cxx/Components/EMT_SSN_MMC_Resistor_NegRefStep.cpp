// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <cmath>
#include <complex>
#include <stdexcept>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Inductor.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Resistor.h>
#include <dpsim-models/EMT/EMT_DC_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>

using namespace DPsim;
using namespace CPS;

namespace {

/// Harmony's V_m is the phase-voltage amplitude. DPsim's three-phase EMT
/// voltage source stores a phase-to-phase RMS phasor and internally multiplies
/// it by RMS3PH_TO_PEAK1PH when producing the instantaneous waveform.
MatrixComp
balancedVoltageReferenceFromHarmonyAmplitude(Real phaseVoltageAmplitude) {
  MatrixComp voltageReference = MatrixComp::Zero(3, 1);

  const Real lineToLineRmsVoltage = phaseVoltageAmplitude / RMS3PH_TO_PEAK1PH;
  const Complex phaseA(lineToLineRmsVoltage, 0.0);
  voltageReference(0, 0) = phaseA;
  voltageReference(1, 0) = phaseA * SHIFT_TO_PHASE_B;
  voltageReference(2, 0) = phaseA * SHIFT_TO_PHASE_C;

  return voltageReference;
}

struct AcInitialOperatingPoint {
  Complex sourceVoltagePeak;
  Complex afterSourceResistanceVoltagePeak;
  Complex afterBranchResistanceVoltagePeak;
  Complex mmcVoltagePeak;
  Complex currentPeak;
};

/// Calculates the exact high-voltage solution for the balanced zero-Q
/// operating point
///
///   source -- R_source -- R_branch -- jX_branch -- MMC(P,Q=0).
///
/// All phasors use phase-to-neutral peak values. The equation is
///
///   S = 3/2 * V_mmc * conj(I)
///   V_source = V_mmc + (R_total + jX_branch) * I.
///
/// No additional physical parameter is introduced; every quantity is derived
/// from the supplied Harmony source, branch and power data.
AcInitialOperatingPoint
calculateZeroQAcOperatingPoint(Real sourceVoltageAmplitude, Real activePower,
                               Real sourceResistance, Real branchResistance,
                               Real branchReactance) {

  if (sourceVoltageAmplitude <= 0.0)
    throw std::invalid_argument("sourceVoltageAmplitude must be positive.");
  if (!std::isfinite(activePower))
    throw std::invalid_argument("activePower must be finite.");

  const Real totalResistance = sourceResistance + branchResistance;
  const Real powerFactor = activePower / 1.5;

  // Let y = |V_mmc|^2. For Q=0:
  //
  //   |V_source|^2 y =
  //       (y + R_total * powerFactor)^2
  //       + (X_branch * powerFactor)^2.
  //
  // The larger root is the physically relevant high-voltage solution.
  const Real quadraticB = 2.0 * totalResistance * powerFactor -
                          sourceVoltageAmplitude * sourceVoltageAmplitude;
  const Real quadraticC =
      (totalResistance * totalResistance + branchReactance * branchReactance) *
      powerFactor * powerFactor;

  const Real discriminant = quadraticB * quadraticB - 4.0 * quadraticC;
  if (discriminant <= 0.0)
    throw std::runtime_error(
        "Harmony AC operating point has no high-voltage solution.");

  const Real terminalVoltageSquared =
      0.5 * (-quadraticB + std::sqrt(discriminant));
  const Real terminalVoltageMagnitude = std::sqrt(terminalVoltageSquared);
  // Signed current amplitude. A negative value represents power flow
  // from the DC side through the MMC into the AC grid.
  const Real signedCurrentAmplitude = powerFactor / terminalVoltageMagnitude;

  const Real terminalAngle = -std::atan2(
      branchReactance * signedCurrentAmplitude,
      terminalVoltageMagnitude + totalResistance * signedCurrentAmplitude);

  const Complex sourceVoltagePeak(sourceVoltageAmplitude, 0.0);
  const Complex mmcVoltagePeak =
      std::polar(terminalVoltageMagnitude, terminalAngle);
  const Complex currentPeak =
      signedCurrentAmplitude * std::polar(1.0, terminalAngle);

  AcInitialOperatingPoint result;
  result.sourceVoltagePeak = sourceVoltagePeak;
  result.currentPeak = currentPeak;
  result.afterSourceResistanceVoltagePeak =
      sourceVoltagePeak - sourceResistance * currentPeak;
  result.afterBranchResistanceVoltagePeak =
      sourceVoltagePeak - (sourceResistance + branchResistance) * currentPeak;
  result.mmcVoltagePeak = mmcVoltagePeak;

  return result;
}

Complex phasePeakToDpsimInitialPhasor(const Complex &phasePeakVoltage) {
  return phasePeakVoltage / RMS3PH_TO_PEAK1PH;
}

} // namespace

int main(int argc, char *argv[]) {
  // -------------------------------------------------------------------------
  // DPsim-only numerical settings.
  // Harmony's point-to-point example is an OPF example and does not specify
  // an EMT time step or transient simulation duration. These are retained
  // from the preceding DPsim MMC test rather than presented as Harmony data.
  // -------------------------------------------------------------------------
  const String simName = "EMT_SSN_MMC_NegStep";
  const Real timeStep = 100e-6;
  const Real finalTime = 0.80;

  // Bidirectional validation sequence:
  //
  //   +50 MW -> 0 MW -> -20 MW -> 0 MW -> +50 MW
  //
  // The intermediate zero-power intervals make the current and power signs
  // unambiguous and avoid hiding a sign error inside one very large reversal.
  // The negative reference is a DPsim validation disturbance; it is not an
  // additional Harmony hardware parameter.
  const Real zeroPowerTime1 = 0.20;
  const Real negativePowerTime = 0.35;
  const Real zeroPowerTime2 = 0.50;
  const Real restorePositivePowerTime = 0.65;

  const Real zeroActivePowerReference = 0.0;
  const Real negativeActivePowerReference = -20.0e6;

  // -------------------------------------------------------------------------
  // Exact MMC1 converter parameters from the supplied Harmony point-to-point
  // example:
  //   omega, Pac, Qac, theta, V_m, Pdc, Vdc,
  //   Larm, Rarm, Carm, N, Lreactor, Rreactor, delay
  // -------------------------------------------------------------------------
  const Real nominalFrequency = 50.0;
  const Real nominalOmega = 2.0 * PI * nominalFrequency;

  const Real activePowerReference = 50.0e6;
  const Real reactivePowerReference = 0.0;
  const Real acVoltageAmplitude = 345.0e3;
  const Real dcPower = 50.0e6;
  const Real nominalDcVoltage = 440.0e3;
  const Real dcLoadResistance = 44.0e3;
  const Real dcLoadInductance = 1.0;
  const Real dcLoadInitialCurrent = nominalDcVoltage / dcLoadResistance;

  const Real armInductance = 0.05;
  const Real armResistance = 1.07;
  const Real submoduleCapacitance = 0.01;
  const UInt numberOfSubmodules = 400;
  const Real reactorInductance = 0.0005;
  const Real reactorResistance = 0.0001;
  const Real modulationDelay = 0.0;
  (void)dcPower;
  (void)modulationDelay;

  // -------------------------------------------------------------------------
  // Exact AC-source and AC-branch data from Harmony.
  // AC_source: Zsrc = 5 ohm.
  // Branch br1_ac: Z = 5 + j140 ohm.
  // DPsim's EMT inductor takes L, therefore L = X / omega is a direct unit
  // conversion of Harmony's branch reactance, not an additional parameter.
  // -------------------------------------------------------------------------
  const Real sourceResistance = 5.0;
  const Real branchResistance = 5.0;
  const Real branchReactance = 140.0;
  const Real branchInductance = branchReactance / nominalOmega;

  // -------------------------------------------------------------------------
  // Exact MMC1 controller gains and references from the supplied Harmony
  // controller_params1 vector.
  // -------------------------------------------------------------------------
  const Real kpPll = 0.001103374;
  const Real kiPll = 0.00073;

  const Real kpActivePower = 6.6667e-7;
  const Real kiActivePower = 3.3333e-4;

  const Real kpReactivePower = 6.6667e-7;
  const Real kiReactivePower = 3.3333e-4;

  const Real kpEnergy = 120.0;
  const Real kiEnergy = 400.0;

  const Real kpZeroSequenceCurrent = 19.93;
  const Real kiZeroSequenceCurrent = 4500.0;

  const Real kpOutputCurrent = 117.93;
  const Real kiOutputCurrent = 8.5e4;

  const Real kpCirculatingCurrent = 19.93;
  const Real kiCirculatingCurrent = 4500.0;

  const Real iSigmaDReference = 0.0;
  const Real iSigmaQReference = 0.0;

  // Exact ZCC reference contained in Harmony controller_params1.
  // A Harmony-equivalent energy-controller implementation must replace this
  // reference at runtime when the energy controller is active; the test does
  // not alter the published constructor value to compensate for model code.
  const Real iSigmaZReference = 0.0;

  Logger::setLogDir("logs/" + simName);

  // -------------------------------------------------------------------------
  // Nodes
  // -------------------------------------------------------------------------
  auto nAcSource = SimNode<Real>::make("nAcSource", PhaseType::ABC);
  auto nAcAfterSourceResistance =
      SimNode<Real>::make("nAcAfterSourceResistance", PhaseType::ABC);
  auto nAcAfterBranchResistance =
      SimNode<Real>::make("nAcAfterBranchResistance", PhaseType::ABC);
  auto nMmcAc = SimNode<Real>::make("nMmcAc", PhaseType::ABC);
  auto nDcPositive = SimNode<Real>::make("nDcPositive", PhaseType::DC);
  auto nDcLoadMid = SimNode<Real>::make("nDcLoadMid", PhaseType::DC);

  // Build one consistent loaded AC operating point using the unchanged total
  // external impedance:
  //
  //   Ztotal = sourceResistance + branchResistance + j*branchReactance.
  //
  // The branch implementation changes from separate EMT R and L components to
  // one SSN RL state model, but the physical Harmony data remain unchanged.
  const AcInitialOperatingPoint acOperatingPoint =
      calculateZeroQAcOperatingPoint(acVoltageAmplitude, activePowerReference,
                                     sourceResistance, branchResistance,
                                     branchReactance);

  nAcSource->setInitialVoltage(
      phasePeakToDpsimInitialPhasor(acOperatingPoint.sourceVoltagePeak));
  nAcAfterSourceResistance->setInitialVoltage(phasePeakToDpsimInitialPhasor(
      acOperatingPoint.afterSourceResistanceVoltagePeak));
  nAcAfterBranchResistance->setInitialVoltage(phasePeakToDpsimInitialPhasor(
      acOperatingPoint.afterBranchResistanceVoltagePeak));
  nMmcAc->setInitialVoltage(
      phasePeakToDpsimInitialPhasor(acOperatingPoint.mmcVoltagePeak));
  nDcPositive->setInitialVoltage(Complex(nominalDcVoltage, 0.0));
  nDcLoadMid->setInitialVoltage(
      Complex(dcLoadResistance * dcLoadInitialCurrent, 0.0));

  // -------------------------------------------------------------------------
  // Harmony AC source and network impedance
  //
  // Source -- 5 ohm -- [5 ohm + j140 ohm] -- MMC
  //
  // The branch R and L use the existing EMT three-phase primitives.
  // -------------------------------------------------------------------------
  auto acSource =
      EMT::Ph3::NetworkInjection::make("AcSource", Logger::Level::info);
  acSource->setParameters(
      balancedVoltageReferenceFromHarmonyAmplitude(acVoltageAmplitude),
      nominalFrequency);

  auto sourceResistor =
      EMT::Ph3::Resistor::make("SourceResistance", Logger::Level::info);
  sourceResistor->setParameters(
      Math::singlePhaseParameterToThreePhase(sourceResistance));

  auto branchResistor =
      EMT::Ph3::Resistor::make("BranchResistance", Logger::Level::info);
  branchResistor->setParameters(
      Math::singlePhaseParameterToThreePhase(branchResistance));

  auto branchInductor =
      EMT::Ph3::Inductor::make("BranchInductance", Logger::Level::info);
  branchInductor->setParameters(
      Math::singlePhaseParameterToThreePhase(branchInductance));

  acSource->connect({nAcSource});
  sourceResistor->connect({nAcSource, nAcAfterSourceResistance});

  branchResistor->connect({nAcAfterSourceResistance, nAcAfterBranchResistance});
  branchInductor->connect({nAcAfterBranchResistance, nMmcAc});

  // -------------------------------------------------------------------------
  // Stiff 440 kV DC terminal
  //
  // EMT::DC::VoltageSource enforces v(terminal 1) - v(terminal 0) = Vref.
  // Connecting terminal 0 to ground therefore fixes nDcPositive at +440 kV.
  // -------------------------------------------------------------------------
  auto dcSource = EMT::DC::VoltageSource::make("DcSource", Logger::Level::info);
  dcSource->setParameters(nominalDcVoltage);
  dcSource->connect({SimNode<Real>::GND, nDcPositive});

  // Scalar-DC series R-L regression. Positive branch current flows
  // nDcPositive -> inductor -> resistor -> ground, matching each component's
  // terminal-1-to-terminal-0 passive convention.
  auto dcLoadInductor =
      EMT::DC::SSN::Inductor::make("DcLoadInductor", Logger::Level::info);
  dcLoadInductor->setParameters(dcLoadInductance, dcLoadInitialCurrent);
  dcLoadInductor->connect({nDcLoadMid, nDcPositive});

  auto dcLoadResistor =
      EMT::DC::SSN::Resistor::make("DcLoadResistor", Logger::Level::info);
  dcLoadResistor->setParameters(dcLoadResistance);
  dcLoadResistor->connect({SimNode<Real>::GND, nDcLoadMid});

  // -------------------------------------------------------------------------
  // Harmony MMC1
  // -------------------------------------------------------------------------
  auto mmc = EMT::Ph3::SSN_MMC::make("MMC1", "MMC1", Logger::Level::debug);

  mmc->setParameters(nominalFrequency, acVoltageAmplitude, nominalDcVoltage,
                     armInductance, armResistance, submoduleCapacitance,
                     numberOfSubmodules, reactorInductance, reactorResistance);

  // Harmony's standalone converter parameter contains theta = 0. The
  // network operating point introduces a terminal-bus angle through the
  // specified 10 + j140 ohm source/branch impedance. Use that derived bus
  // angle for the EMT dq clock so the PLL starts with v_q = 0.
  mmc->setInitialAngle(std::arg(acOperatingPoint.mmcVoltagePeak));
  mmc->setPLL(kpPll, kiPll, true);

  mmc->setActivePowerControl(activePowerReference, kpActivePower,
                             kiActivePower);

  auto activePowerReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(activePowerReference));
  mmc->setReactivePowerControl(reactivePowerReference, kpReactivePower,
                               kiReactivePower);

  mmc->setEnergyController(kpEnergy, kiEnergy, true);

  mmc->setOutputCurrentController(kpOutputCurrent, kiOutputCurrent);
  mmc->setCirculatingCurrentController(kpCirculatingCurrent,
                                       kiCirculatingCurrent);
  mmc->setZeroSequenceCurrentController(kpZeroSequenceCurrent,
                                        kiZeroSequenceCurrent);

  mmc->setCirculatingCurrentReferences(iSigmaDReference, iSigmaQReference,
                                       iSigmaZReference);

  // Harmony specifies delay = 0.0 for this converter and supplies no filter
  // parameter vector. No delay or measurement-filter values are invented here.

  // DPsim-only numerical diagnostics. They do not alter the physical model.
  // The class also performs its nonlinear operating-point initialization by
  // default before building the first SSN equivalent.
  mmc->setEigenvalueDiagnostics(true, 100);
  mmc->setDiagnosticTimeStep(timeStep);

  // SSN_MMC terminal order:
  //   0 -> AC abc
  //   1 -> DC+
  //   2 -> DC-
  mmc->connect({nMmcAc, nDcPositive, SimNode<Real>::GND});

  // -------------------------------------------------------------------------
  // System topology
  // -------------------------------------------------------------------------
  auto system = SystemTopology(nominalFrequency,
                               SystemNodeList{
                                   nAcSource,
                                   nAcAfterSourceResistance,
                                   nAcAfterBranchResistance,
                                   nMmcAc,
                                   nDcPositive,
                                   nDcLoadMid,
                               },
                               SystemComponentList{
                                   acSource,
                                   sourceResistor,
                                   branchResistor,
                                   branchInductor,
                                   dcSource,
                                   dcLoadInductor,
                                   dcLoadResistor,
                                   mmc,
                               });

  // -------------------------------------------------------------------------
  // Logging
  // -------------------------------------------------------------------------
  auto logger = DataLogger::make(simName);

  logger->logAttribute("Voltage_AC_Source", nAcSource->attribute("v"));
  logger->logAttribute("Voltage_AC_MMC", nMmcAc->attribute("v"));
  logger->logAttribute("Voltage_DC_Positive", nDcPositive->attribute("v"));
  logger->logAttribute("Voltage_DC_LoadMid", nDcLoadMid->attribute("v"));

  logger->logAttribute("Current_SourceResistance",
                       sourceResistor->attribute("i_intf"));
  logger->logAttribute("Current_BranchSeriesRL",
                       branchInductor->attribute("i_intf"));
  logger->logAttribute("Current_DC_Source", dcSource->attribute("i_intf"));
  logger->logAttribute("Current_DC_LoadInductor",
                       dcLoadInductor->attribute("i_intf"));
  logger->logAttribute("Current_DC_LoadResistor",
                       dcLoadResistor->attribute("i_intf"));

  logger->logAttribute("MMC_InterfaceVoltage", mmc->attribute("v_intf"));
  logger->logAttribute("MMC_InterfaceCurrent", mmc->attribute("i_intf"));

  logger->logAttribute("MMC_Vdc", mmc->attribute("vdc"));
  logger->logAttribute("MMC_Vdcp", mmc->attribute("vdcp"));
  logger->logAttribute("MMC_Vdcn", mmc->attribute("vdcn"));
  logger->logAttribute("MMC_Idc", mmc->attribute("idc"));
  logger->logAttribute("MMC_ActivePowerReference", activePowerReferenceSignal);
  logger->logAttribute("MMC_Pac", mmc->attribute("p_ac"));
  logger->logAttribute("MMC_Qac", mmc->attribute("q_ac"));
  logger->logAttribute("MMC_Pdc", mmc->attribute("p_dc"));
  logger->logAttribute("MMC_PowerBalanceError",
                       mmc->attribute("power_balance_error"));

  logger->logAttribute("MMC_IDeltaD", mmc->attribute("i_delta_d"));
  logger->logAttribute("MMC_IDeltaDReference", mmc->attribute("i_delta_d_ref"));
  logger->logAttribute("MMC_IDeltaQ", mmc->attribute("i_delta_q"));
  logger->logAttribute("MMC_IDeltaQReference", mmc->attribute("i_delta_q_ref"));
  logger->logAttribute("MMC_ISigmaZ", mmc->attribute("i_sigma_z"));
  logger->logAttribute("MMC_ISigmaZReference", mmc->attribute("i_sigma_z_ref"));
  logger->logAttribute("MMC_VacMagnitude", mmc->attribute("v_ac_mag"));
  logger->logAttribute("MMC_StoredEnergy", mmc->attribute("stored_energy"));
  logger->logAttribute("MMC_Theta", mmc->attribute("theta"));
  logger->logAttribute("MMC_PLLFrequency", mmc->attribute("pll_frequency"));

  logger->logAttribute("MMC_PLLError", mmc->attribute("pll_error"));
  logger->logAttribute("MMC_PLLAngleDeviation",
                       mmc->attribute("pll_angle_deviation"));
  logger->logAttribute("MMC_ControlVoltageD", mmc->attribute("v_control_d"));
  logger->logAttribute("MMC_ControlVoltageQ", mmc->attribute("v_control_q"));

  logger->logAttribute("MMC_StateNorm", mmc->attribute("state_norm"));
  logger->logAttribute("MMC_StateDerivativeNorm",
                       mmc->attribute("state_derivative_norm"));
  logger->logAttribute("MMC_EquilibriumResidualNorm",
                       mmc->attribute("equilibrium_residual_norm"));
  logger->logAttribute("MMC_JacobianMaxRealEigenvalue",
                       mmc->attribute("jacobian_max_real_eigenvalue"));
  logger->logAttribute("MMC_JacobianMaxAbsEigenvalue",
                       mmc->attribute("jacobian_max_abs_eigenvalue"));
  logger->logAttribute("MMC_JacobianDominantFrequency",
                       mmc->attribute("jacobian_dominant_frequency"));
  logger->logAttribute("MMC_DiagnosticsValid",
                       mmc->attribute("diagnostics_valid"));

  logger->logAttribute("MMC_JacobianMaxDiscreteMagnitude",
                       mmc->attribute("jacobian_max_discrete_magnitude"));
  logger->logAttribute("MMC_JacobianDiscreteDominantFrequency",
                       mmc->attribute("jacobian_discrete_dominant_frequency"));

  // -------------------------------------------------------------------------
  // Simulation
  // -------------------------------------------------------------------------
  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);

  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);

  sim.addLogger(logger);

  Bool setFirstZeroReference = true;
  Bool setNegativeReference = true;
  Bool setSecondZeroReference = true;
  Bool restorePositiveReference = true;

  auto applyActivePowerReference = [&](Real reference,
                                       const char *description) {
    SPDLOG_INFO("{}: setting MMC active-power reference to {:.3f} MW "
                "at t = {:.6f} s.",
                description, reference / 1.0e6, sim.time());

    mmc->setActivePowerControl(reference, kpActivePower, kiActivePower);

    activePowerReferenceSignal->set(reference);
  };

  sim.initialize();
  sim.start();

  while (sim.time() < sim.finalTime()) {
    if (setFirstZeroReference && sim.time() >= zeroPowerTime1) {
      applyActivePowerReference(zeroActivePowerReference,
                                "First zero-power transition");
      setFirstZeroReference = false;
    }

    if (setNegativeReference && sim.time() >= negativePowerTime) {
      applyActivePowerReference(negativeActivePowerReference,
                                "Negative-power transition");
      setNegativeReference = false;
    }

    if (setSecondZeroReference && sim.time() >= zeroPowerTime2) {
      applyActivePowerReference(zeroActivePowerReference,
                                "Second zero-power transition");
      setSecondZeroReference = false;
    }

    if (restorePositiveReference && sim.time() >= restorePositivePowerTime) {
      applyActivePowerReference(activePowerReference,
                                "Positive-power restoration");
      restorePositiveReference = false;
    }

    sim.step();
  }

  sim.stop();

  if (!mmc->getState().allFinite() || !mmc->getInterfaceVoltage().allFinite() ||
      !mmc->getInterfaceCurrent().allFinite())
    throw std::runtime_error("MMC state or interface contains NaN or Inf.");

  const Real inductorCurrent = dcLoadInductor->intfCurrent()(0, 0);
  const Real resistorCurrent = dcLoadResistor->intfCurrent()(0, 0);
  const Real dcNodeKclResidual = inductorCurrent +
                                 dcSource->intfCurrent()(0, 0) +
                                 mmc->getInterfaceCurrent()(3, 0);
  const Real finalPac = mmc->attributeTyped<Real>("p_ac")->get();
  const Real finalPdc = mmc->attributeTyped<Real>("p_dc")->get();

  if (!std::isfinite(dcNodeKclResidual) || std::abs(dcNodeKclResidual) > 1.0e-6)
    throw std::runtime_error(
        fmt::format("scalar DC-node KCL residual is {} A", dcNodeKclResidual));
  if (std::abs(inductorCurrent - resistorCurrent) > 1.0e-6)
    throw std::runtime_error(
        fmt::format("scalar DC series R-L current discontinuity is {} A",
                    inductorCurrent - resistorCurrent));

  SPDLOG_INFO("Scalar-DC MMC R-L validation: Vdc={} V, IL={} A, IR={} A, "
              "KCL residual={} A, Pac={} W, Pdc={} W, Pac-Pdc={} W",
              nDcPositive->voltage()(0, 0), inductorCurrent, resistorCurrent,
              dcNodeKclResidual, finalPac, finalPdc, finalPac - finalPdc);

  return 0;
}
