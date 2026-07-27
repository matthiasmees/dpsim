// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <cmath>
#include <complex>
#include <stdexcept>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_CurrentSource.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Capacitor.h>
#include <dpsim-models/EMT/EMT_DC_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_DC_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>

using namespace DPsim;
using namespace CPS;

namespace {

MatrixComp balancedVoltageReference(Real phaseVoltageAmplitude) {
  MatrixComp reference = MatrixComp::Zero(3, 1);
  const Real lineToLineRms = phaseVoltageAmplitude / RMS3PH_TO_PEAK1PH;
  const Complex phaseA(lineToLineRms, 0.0);
  reference(0, 0) = phaseA;
  reference(1, 0) = phaseA * SHIFT_TO_PHASE_B;
  reference(2, 0) = phaseA * SHIFT_TO_PHASE_C;
  return reference;
}

void requireFinite(const Matrix &value, const String &name) {
  if (!value.allFinite())
    throw std::runtime_error(name + " contains NaN or Inf.");
}

void requireNear(Real actual, Real expected, Real tolerance,
                 const String &name) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
    throw std::runtime_error(
        fmt::format("{}: actual={} expected={} tolerance={}", name, actual,
                    expected, tolerance));
}

} // namespace

int main(int argc, char *argv[]) {
  const String simName = "EMT_SSN_MMC_PiLine";

  const Real timeStep = 50.0e-6;
  const Real finalTime = 0.10;

  const Real nominalFrequency = 50.0;
  const Real nominalOmega = 2.0 * PI * nominalFrequency;
  const Real acVoltageAmplitude = 345.0e3;
  const Real nominalDcVoltage = 440.0e3;
  const Real nominalPoleVoltage = nominalDcVoltage / 2.0;

  const Real acSeriesResistance = 10.0;
  const Real acSeriesReactance = 140.0;
  const Real acSeriesInductance = acSeriesReactance / nominalOmega;

  const Real armInductance = 0.05;
  const Real armResistance = 1.07;
  const Real submoduleCapacitance = 0.01;
  const UInt numberOfSubmodules = 400;
  const Real reactorInductance = 0.0005;
  const Real reactorResistance = 0.0001;

  const Real effectiveDcCapacitance =
      6.0 * submoduleCapacitance / static_cast<Real>(numberOfSubmodules);

  // The Delft script gives R_cable=0.5 ohm and L_cable=15 mH but does not
  // state textually whether these are per conductor or loop totals. This
  // Phase-1 study explicitly treats them as pole-to-pole loop totals and
  // divides them equally between two symmetric scalar pole conductors.
  const Real dcCableLoopResistance = 0.5;
  const Real dcCableLoopInductance = 15.0e-3;
  const Real dcLineSeriesResistance = dcCableLoopResistance / 2.0;
  const Real dcLineSeriesInductance = dcCableLoopInductance / 2.0;

  // The reference script provides no cable capacitance or conductance. Keep
  // these explicit example parameters at zero rather than inventing values.
  const Real dcLineParallelCapacitance = 0.0;
  const Real dcLineParallelConductance = 0.0;
  const Real dcLineInitialCurrent = 0.0;

  const Real kpPll = 0.001103374;
  const Real kiPll = 0.00073;
  const Real kpOutputCurrent = 117.93;
  const Real kiOutputCurrent = 8.5e4;
  const Real kpCirculatingCurrent = 19.93;
  const Real kiCirculatingCurrent = 4500.0;
  const Real kpZeroSequenceCurrent = 19.93;
  const Real kiZeroSequenceCurrent = 4500.0;
  const Real kpEnergy = 120.0;
  const Real kiEnergy = 400.0;

  const Real dcVoltageLoopBandwidth = 2.0 * PI * 3.0;
  const Real dcVoltageLoopDamping = 1.0;
  const Real acCurrentToDcCurrentGain =
      1.5 * acVoltageAmplitude / nominalDcVoltage;
  const Real kpDcVoltage = 2.0 * dcVoltageLoopDamping * dcVoltageLoopBandwidth *
                           effectiveDcCapacitance / acCurrentToDcCurrentGain;
  const Real kiDcVoltage = dcVoltageLoopBandwidth * dcVoltageLoopBandwidth *
                           effectiveDcCapacitance / acCurrentToDcCurrentGain;

  Logger::setLogDir("logs/" + simName);

  auto nAcSource = SimNode<Real>::make("nAcSource", PhaseType::ABC);
  auto nAcAfterResistance =
      SimNode<Real>::make("nAcAfterResistance", PhaseType::ABC);
  auto nMmcAc = SimNode<Real>::make("nMmcAc", PhaseType::ABC);
  auto nMmcDcPositive = SimNode<Real>::make("nMmcDcPositive", PhaseType::DC);
  auto nMmcDcNegative = SimNode<Real>::make("nMmcDcNegative", PhaseType::DC);
  auto nRemoteDcPositive =
      SimNode<Real>::make("nRemoteDcPositive", PhaseType::DC);
  auto nRemoteDcNegative =
      SimNode<Real>::make("nRemoteDcNegative", PhaseType::DC);

  const Complex initialAcPhasor(acVoltageAmplitude / RMS3PH_TO_PEAK1PH, 0.0);
  nAcSource->setInitialVoltage(initialAcPhasor);
  nAcAfterResistance->setInitialVoltage(initialAcPhasor);
  nMmcAc->setInitialVoltage(initialAcPhasor);
  nMmcDcPositive->setInitialVoltage(Complex(+nominalPoleVoltage, 0.0));
  nMmcDcNegative->setInitialVoltage(Complex(-nominalPoleVoltage, 0.0));
  nRemoteDcPositive->setInitialVoltage(Complex(+nominalPoleVoltage, 0.0));
  nRemoteDcNegative->setInitialVoltage(Complex(-nominalPoleVoltage, 0.0));

  auto acSource =
      EMT::Ph3::NetworkInjection::make("AcSource", Logger::Level::info);
  acSource->setParameters(balancedVoltageReference(acVoltageAmplitude),
                          nominalFrequency);
  acSource->connect({nAcSource});

  auto acSeriesResistor =
      EMT::Ph3::Resistor::make("AcSeriesResistance", Logger::Level::info);
  acSeriesResistor->setParameters(
      Math::singlePhaseParameterToThreePhase(acSeriesResistance));
  acSeriesResistor->connect({nAcSource, nAcAfterResistance});

  auto acSeriesInductor =
      EMT::Ph3::Inductor::make("AcSeriesInductance", Logger::Level::info);
  acSeriesInductor->setParameters(
      Math::singlePhaseParameterToThreePhase(acSeriesInductance));
  acSeriesInductor->connect({nAcAfterResistance, nMmcAc});

  auto positivePoleLine =
      EMT::DC::SSN::PiLine::make("DcPositivePolePiLine", Logger::Level::info);
  positivePoleLine->setParameters(
      dcLineSeriesResistance, dcLineSeriesInductance, dcLineParallelCapacitance,
      dcLineParallelConductance, dcLineInitialCurrent);

  // PiLine positive current is terminal 1 -> terminal 0. This orientation
  // makes positive current flow from the MMC positive pole to the remote
  // positive pole.
  positivePoleLine->connect({nRemoteDcPositive, nMmcDcPositive});

  auto negativePoleLine =
      EMT::DC::SSN::PiLine::make("DcNegativePolePiLine", Logger::Level::info);
  negativePoleLine->setParameters(
      dcLineSeriesResistance, dcLineSeriesInductance, dcLineParallelCapacitance,
      dcLineParallelConductance, dcLineInitialCurrent);

  // Positive return current flows from the remote negative pole to the MMC
  // negative pole, again terminal 1 -> terminal 0.
  negativePoleLine->connect({nMmcDcNegative, nRemoteDcNegative});

  // This is a single-station line test, not the complete P2P system. Two
  // scalar remote sources represent the second DC system and establish the
  // otherwise-unconstrained bipolar common mode at +/-220 kV.
  auto remotePositiveVoltageSource = EMT::DC::VoltageSource::make(
      "RemotePositiveVoltageSource", Logger::Level::info);
  remotePositiveVoltageSource->setParameters(+nominalPoleVoltage);
  remotePositiveVoltageSource->connect({SimNode<Real>::GND, nRemoteDcPositive});

  auto remoteNegativeVoltageSource = EMT::DC::VoltageSource::make(
      "RemoteNegativeVoltageSource", Logger::Level::info);
  remoteNegativeVoltageSource->setParameters(+nominalPoleVoltage);
  remoteNegativeVoltageSource->connect({nRemoteDcNegative, SimNode<Real>::GND});

  auto dcCurrentSource = EMT::DC::CurrentSource::make("DcRemoteCurrentSource",
                                                      Logger::Level::info);
  dcCurrentSource->setParameters(0.0);

  // Positive source current flows terminal 1 -> terminal 0, from the remote
  // negative pole to the remote positive pole.
  dcCurrentSource->connect({nRemoteDcPositive, nRemoteDcNegative});

  auto effectiveDcCapacitor = EMT::DC::SSN::Capacitor::make(
      "EffectiveDcCapacitor", Logger::Level::info);
  effectiveDcCapacitor->setParameters(effectiveDcCapacitance);
  effectiveDcCapacitor->connect({nMmcDcNegative, nMmcDcPositive});

  auto mmc = EMT::Ph3::SSN_MMC::make("MMC1", "MMC1", Logger::Level::debug);
  mmc->setParameters(nominalFrequency, acVoltageAmplitude, nominalDcVoltage,
                     armInductance, armResistance, submoduleCapacitance,
                     numberOfSubmodules, reactorInductance, reactorResistance);
  mmc->setInitialAngle(0.0);
  mmc->setPLL(kpPll, kiPll, true);
  mmc->setReactiveControlOpenLoop(0.0);
  mmc->setDcVoltageControl(nominalDcVoltage, kpDcVoltage, kiDcVoltage);
  mmc->setEnergyController(kpEnergy, kiEnergy, true);
  mmc->setOutputCurrentController(kpOutputCurrent, kiOutputCurrent);
  mmc->setCirculatingCurrentController(kpCirculatingCurrent,
                                       kiCirculatingCurrent);
  mmc->setZeroSequenceCurrentController(kpZeroSequenceCurrent,
                                        kiZeroSequenceCurrent);
  mmc->setCirculatingCurrentReferences(0.0, 0.0, 0.0);
  mmc->setLimits(100.0, 100.0, 2.0);
  mmc->setOperatingPointInitialization(true, 50, 1.0e-8);
  mmc->setEigenvalueDiagnostics(true, 200);
  mmc->setDiagnosticTimeStep(timeStep);
  mmc->connect({nMmcAc, nMmcDcPositive, nMmcDcNegative});

  auto dcCurrentReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(0.0));
  auto dcVoltageReferenceSignal =
      Attribute<Real>::Ptr(AttributeStatic<Real>::make(nominalDcVoltage));
  // The present averaged MMC is initialized deblocked and has no blocked
  // plant mode. Log that explicit Phase-1 limitation as a constant state.
  auto deblockedSignal = Attribute<Real>::Ptr(AttributeStatic<Real>::make(1.0));

  auto system = SystemTopology(
      nominalFrequency,
      SystemNodeList{nAcSource, nAcAfterResistance, nMmcAc, nMmcDcPositive,
                     nMmcDcNegative, nRemoteDcPositive, nRemoteDcNegative},
      SystemComponentList{acSource, acSeriesResistor, acSeriesInductor,
                          positivePoleLine, negativePoleLine,
                          remotePositiveVoltageSource,
                          remoteNegativeVoltageSource, dcCurrentSource,
                          effectiveDcCapacitor, mmc});

  auto logger = DataLogger::make(simName);
  logger->logAttribute("Voltage_AC_Source", nAcSource->attribute("v"));
  logger->logAttribute("Voltage_AC_MMC", nMmcAc->attribute("v"));
  logger->logAttribute("Current_AC_SeriesRL",
                       acSeriesInductor->attribute("i_intf"));
  logger->logAttribute("MMC_InterfaceCurrent", mmc->attribute("i_intf"));
  logger->logAttribute("Voltage_DC_MMC_Positive",
                       nMmcDcPositive->attribute("v"));
  logger->logAttribute("Voltage_DC_MMC_Negative",
                       nMmcDcNegative->attribute("v"));
  logger->logAttribute("Voltage_DC_Remote_Positive",
                       nRemoteDcPositive->attribute("v"));
  logger->logAttribute("Voltage_DC_Remote_Negative",
                       nRemoteDcNegative->attribute("v"));
  logger->logAttribute("Voltage_DC_PositiveLine",
                       positivePoleLine->attribute("v_intf"));
  logger->logAttribute("Voltage_DC_NegativeLine",
                       negativePoleLine->attribute("v_intf"));
  logger->logAttribute("Current_DC_PositiveLine",
                       positivePoleLine->attribute("i_intf"));
  logger->logAttribute("Current_DC_NegativeLine",
                       negativePoleLine->attribute("i_intf"));
  logger->logAttribute("DcCurrentReference", dcCurrentReferenceSignal);
  logger->logAttribute("Current_DC_RemoteSource",
                       dcCurrentSource->attribute("i_intf"));
  logger->logAttribute("Current_DC_RemotePositiveVoltageSource",
                       remotePositiveVoltageSource->attribute("i_intf"));
  logger->logAttribute("Current_DC_RemoteNegativeVoltageSource",
                       remoteNegativeVoltageSource->attribute("i_intf"));
  logger->logAttribute("Current_DC_EffectiveCapacitor",
                       effectiveDcCapacitor->attribute("i_intf"));
  logger->logAttribute("VdcReference", dcVoltageReferenceSignal);
  logger->logAttribute("MMC_Vdc", mmc->attribute("vdc"));
  logger->logAttribute("MMC_Vdcp", mmc->attribute("vdcp"));
  logger->logAttribute("MMC_Vdcn", mmc->attribute("vdcn"));
  logger->logAttribute("MMC_Idc", mmc->attribute("idc"));
  logger->logAttribute("MMC_Pac", mmc->attribute("p_ac"));
  logger->logAttribute("MMC_Qac", mmc->attribute("q_ac"));
  logger->logAttribute("MMC_Pdc", mmc->attribute("p_dc"));
  logger->logAttribute("MMC_PowerBalanceError",
                       mmc->attribute("power_balance_error"));
  logger->logAttribute("MMC_StoredEnergy", mmc->attribute("stored_energy"));
  logger->logAttribute("MMC_StateNorm", mmc->attribute("state_norm"));
  logger->logAttribute("MMC_StateDerivativeNorm",
                       mmc->attribute("state_derivative_norm"));
  logger->logAttribute("MMC_DiagnosticsValid",
                       mmc->attribute("diagnostics_valid"));
  logger->logAttribute("MMC_Deblocked", deblockedSignal);

  Simulation sim(simName, Logger::Level::info);
  sim.setSystem(system);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setDomain(Domain::EMT);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.addLogger(logger);

  Real positiveCurrentBeforeStep = dcLineInitialCurrent;
  Real negativeCurrentBeforeStep = dcLineInitialCurrent;
  Real localVdcBeforeStep = nominalDcVoltage;
  Bool continuityChecked = false;

  sim.initialize();
  positiveCurrentBeforeStep = positivePoleLine->intfCurrent()(0, 0);
  negativeCurrentBeforeStep = negativePoleLine->intfCurrent()(0, 0);
  localVdcBeforeStep =
      nMmcDcPositive->voltage()(0, 0) - nMmcDcNegative->voltage()(0, 0);
  sim.start();

  while (sim.time() < sim.finalTime()) {
    sim.step();

    if (!continuityChecked) {
      requireNear(positivePoleLine->intfCurrent()(0, 0),
                  positiveCurrentBeforeStep, 0.05,
                  "positive-pole PiLine current continuity");
      requireNear(negativePoleLine->intfCurrent()(0, 0),
                  negativeCurrentBeforeStep, 0.05,
                  "negative-pole PiLine current continuity");
      const Real localVdc =
          nMmcDcPositive->voltage()(0, 0) - nMmcDcNegative->voltage()(0, 0);
      requireNear(localVdc, localVdcBeforeStep, 5.0,
                  "effective DC capacitor voltage continuity");
      continuityChecked = true;
    }
  }

  sim.stop();

  requireFinite(mmc->getState(), "MMC state");
  requireFinite(mmc->getInterfaceVoltage(), "MMC interface voltage");
  requireFinite(mmc->getInterfaceCurrent(), "MMC interface current");
  requireFinite(positivePoleLine->intfVoltage(), "positive-pole line voltage");
  requireFinite(positivePoleLine->intfCurrent(), "positive-pole line current");
  requireFinite(negativePoleLine->intfVoltage(), "negative-pole line voltage");
  requireFinite(negativePoleLine->intfCurrent(), "negative-pole line current");

  const Real finalVdcp = nMmcDcPositive->voltage()(0, 0);
  const Real finalVdcn = nMmcDcNegative->voltage()(0, 0);
  const Real finalVdc = finalVdcp - finalVdcn;
  const Real finalPositiveLineCurrent = positivePoleLine->intfCurrent()(0, 0);
  const Real finalNegativeLineCurrent = negativePoleLine->intfCurrent()(0, 0);
  const Real finalCapacitorCurrent = effectiveDcCapacitor->intfCurrent()(0, 0);
  const Real positiveNodeKclResidual = finalPositiveLineCurrent +
                                       finalCapacitorCurrent +
                                       mmc->getInterfaceCurrent()(3, 0);
  const Real negativeNodeKclResidual = -finalNegativeLineCurrent -
                                       finalCapacitorCurrent +
                                       mmc->getInterfaceCurrent()(4, 0);
  const Real finalPac = mmc->attributeTyped<Real>("p_ac")->get();
  const Real finalPdc = mmc->attributeTyped<Real>("p_dc")->get();

  if (!(finalVdcp > 0.0 && finalVdcn < 0.0))
    throw std::runtime_error("MMC DC pole-voltage polarity is invalid.");
  requireNear(finalVdc, nominalDcVoltage, 2000.0,
              "final MMC pole-to-pole voltage");
  requireNear(finalPositiveLineCurrent, finalNegativeLineCurrent, 0.05,
              "bipolar line-current continuity");
  requireNear(positiveNodeKclResidual, 0.0, 1.0e-6, "local dc+ node KCL");
  requireNear(negativeNodeKclResidual, 0.0, 1.0e-6, "local dc- node KCL");

  SPDLOG_INFO(
      "Scalar-DC MMC PiLine validation: Vdcp={} V, Vdcn={} V, Vdc={} V, "
      "Ipositive={} A, Inegative={} A, dc+ KCL={} A, dc- KCL={} A, "
      "Pac={} W, Pdc={} W, Pac-Pdc={} W",
      finalVdcp, finalVdcn, finalVdc, finalPositiveLineCurrent,
      finalNegativeLineCurrent, positiveNodeKclResidual,
      negativeNodeKclResidual, finalPac, finalPdc, finalPac - finalPdc);

  return 0;
}
