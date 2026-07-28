// SPDX-License-Identifier: MPL-2.0
#include "../Examples.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Capacitor.h>
#include <dpsim-models/EMT/EMT_DC_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_DC_VoltageSource.h>
#include <dpsim-models/EMT/EMT_Ph1_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMCStation.h>

using namespace CPS;
using namespace DPsim;

namespace {
MatrixComp balancedVoltageReference(Real amplitude) {
  MatrixComp reference(3, 1);
  const Complex phaseA(amplitude / RMS3PH_TO_PEAK1PH, 0.0);
  reference << phaseA, phaseA * SHIFT_TO_PHASE_B, phaseA * SHIFT_TO_PHASE_C;
  return reference;
}
void require(Bool condition, const String &message) {
  if (!condition)
    throw std::runtime_error(message);
}

enum class DiagnosticMode {
  PlantOnly,
  InnerLoop,
  GainSign,
  ProportionalVdc,
  FullPi,
  BoundaryAudit,
  InternalPlant,
  InternalOpenLoop,
  NoEnergy,
  NoCirculatingDq,
  CirculatingDqProportionalOnly,
  CirculatingDqIntegralOnly,
  NoZeroSequence,
  NoPll
};

DiagnosticMode parseMode(int argc, char **argv) {
  if (argc == 1 || String(argv[1]) == "full-pi")
    return DiagnosticMode::FullPi;
  if (String(argv[1]) == "plant-only")
    return DiagnosticMode::PlantOnly;
  if (String(argv[1]) == "inner-loop")
    return DiagnosticMode::InnerLoop;
  if (String(argv[1]) == "gain-sign")
    return DiagnosticMode::GainSign;
  if (String(argv[1]) == "proportional-vdc")
    return DiagnosticMode::ProportionalVdc;
  if (String(argv[1]) == "boundary-audit")
    return DiagnosticMode::BoundaryAudit;
  if (String(argv[1]) == "internal-plant")
    return DiagnosticMode::InternalPlant;
  if (String(argv[1]) == "internal-open-loop")
    return DiagnosticMode::InternalOpenLoop;
  if (String(argv[1]) == "no-energy")
    return DiagnosticMode::NoEnergy;
  if (String(argv[1]) == "no-circulating-dq")
    return DiagnosticMode::NoCirculatingDq;
  if (String(argv[1]) == "circulating-dq-kp-only")
    return DiagnosticMode::CirculatingDqProportionalOnly;
  if (String(argv[1]) == "circulating-dq-ki-only")
    return DiagnosticMode::CirculatingDqIntegralOnly;
  if (String(argv[1]) == "no-zero-sequence")
    return DiagnosticMode::NoZeroSequence;
  if (String(argv[1]) == "no-pll")
    return DiagnosticMode::NoPll;
  throw std::invalid_argument(
      "usage: EMT_SSN_MMCStation_DCPiLine "
      "[plant-only|inner-loop|gain-sign|proportional-vdc|full-pi|"
      "boundary-audit|internal-plant|internal-open-loop|no-energy|"
      "no-circulating-dq|circulating-dq-kp-only|"
      "circulating-dq-ki-only|no-zero-sequence|no-pll]");
}

struct BoundarySnapshot {
  Matrix differential;
  Matrix common;
  Matrix modulation;
  Matrix realized;
  Matrix derivative;
  Matrix A;
  Matrix B;
  Matrix C;
  Matrix D;
};

BoundarySnapshot snapshot(const EMT::Ph3::SSN_MMC::Ptr &mmc) {
  BoundarySnapshot result;
  result.derivative = mmc->getStateDerivative();
  result.differential = **mmc->appliedDifferentialVoltageAttribute();
  result.common = **mmc->appliedCommonModeVoltageAttribute();
  result.modulation = **mmc->appliedModulationAttribute();
  result.realized = **mmc->realizedConverterVoltageAttribute();
  mmc->getLocalLinearization(result.A, result.B, result.C, result.D);
  // The finite differences above leave diagnostic attributes at the final
  // perturbation. Re-evaluate once at the unperturbed state.
  result.derivative = mmc->getStateDerivative();
  result.differential = **mmc->appliedDifferentialVoltageAttribute();
  result.common = **mmc->appliedCommonModeVoltageAttribute();
  result.modulation = **mmc->appliedModulationAttribute();
  result.realized = **mmc->realizedConverterVoltageAttribute();
  return result;
}

Real relativeDifference(const Matrix &a, const Matrix &b) {
  return (a - b).cwiseAbs().maxCoeff() /
         std::max({1.0, a.cwiseAbs().maxCoeff(), b.cwiseAbs().maxCoeff()});
}

Real stateValueByName(const EMT::Ph3::SSN_MMC::Ptr &mmc,
                      const String &stateName) {
  const auto names = mmc->getLocalStateNames();
  const Matrix state = mmc->getState();
  for (std::size_t index = 0; index < names.size(); ++index) {
    if (names[index] == stateName)
      return state(static_cast<Eigen::Index>(index), 0);
  }
  throw std::logic_error("MMC state not found: " + stateName);
}

Real calculateMmcInductorEnergy(const Matrix &state, Real armInductance,
                                Real reactorInductance) {
  const Real iDeltaSquared =
      state(0, 0) * state(0, 0) + state(1, 0) * state(1, 0);
  const Real iSigmaSquared =
      state(3, 0) * state(3, 0) + state(4, 0) * state(4, 0);
  const Real iSigmaZSquared = state(2, 0) * state(2, 0);
  const Real acEquivalentInductance = armInductance / 2.0 + reactorInductance;

  // Energy associated with the dq current states represented explicitly by
  // the averaged MMC model. The capacitor-energy attribute alone does not
  // include these inductive terms.
  return 0.75 * acEquivalentInductance * iDeltaSquared +
         1.5 * armInductance * iSigmaSquared +
         3.0 * armInductance * iSigmaZSquared;
}

Real calculateMmcTotalEnergy(const EMT::Ph3::SSN_MMC::Ptr &mmc,
                             Real armInductance, Real reactorInductance) {
  return **mmc->storedEnergyAttribute() +
         calculateMmcInductorEnergy(mmc->getState(), armInductance,
                                    reactorInductance);
}

struct TraceSample {
  Real time;
  Real localVdc;
  Real remoteVdc;
  Real lineCurrent;
  Real stationVdcError;
  Real stationIdReference;
  Real stationId;
  Real stationVdReference;
  Real stationVqReference;
  Real internalFilteredVdc;
  Real internalIdReference;
  Real internalActiveIntegrator;
  Real mmcDcCurrent;
  Real positiveRcCurrent;
  Real negativeRcCurrent;
  Real mmcCapacitorEnergy;
  Real mmcTotalEnergy;
  Real stationOuterProportional;
  Real stationOuterIntegral;
  Real plantModulationD;
  Real plantModulationQ;
  Real rcCapEnergy;
  Real lineEnergy;
  Real sourcePower;
  Real externalPowerEnergyResidual;
  Real mmcPowerEnergyResidual;
};
} // namespace

int main(int argc, char **argv) {
  const DiagnosticMode mode = parseMode(argc, argv);
  const Real theta = argc > 2 ? std::stod(argv[2]) : 0.5;
  require(theta == 0.5 || theta == 0.55 || theta == 0.6,
          "Diagnostic theta must be exactly 0.5, 0.55, or 0.6.");
  const Real timeStep = argc > 3 ? std::stod(argv[3]) : 40e-6;
  require(timeStep == 40e-6 || timeStep == 20e-6 || timeStep == 10e-6 ||
              timeStep == 5e-6 || timeStep == 2e-6,
          "Diagnostic timestep must be 40, 20, 10, 5, or 2 us.");
  const Bool vdcLoopMode =
      mode == DiagnosticMode::ProportionalVdc || mode == DiagnosticMode::FullPi;
  const Real finalTime = mode == DiagnosticMode::GainSign ? 0.08 : 0.30;
  const Real frequency = 50.0;
  const Real omega = 2.0 * PI * frequency;
  const Real acVoltage = 345e3;
  const Real nominalDcVoltage = 440e3;
  const Real nominalPower = 1e9;

  // MMC plant parameters used both by the component and by the diagnostics.
  const Real armInductance = 0.05;
  const Real armResistance = 1.07;
  const Real submoduleCapacitance = 0.01;
  const UInt numberOfSubmodules = 400;
  const Real reactorInductance = 0.0005;
  const Real reactorResistance = 0.0001;

  // Target DC network from the reference model:
  // two identical RL cable sections per pole and one series RC-to-ground
  // branch at each converter pole.
  const Real cableSectionResistance = 0.25;
  const Real cableSectionInductance = 7.5e-3;
  const Real groundingResistance = 100.0;
  const Real groundingCapacitance = 50e-9;
  // const Real differentialGroundCapacitance = groundingCapacitance / 2.0;

  auto acNode = SimNode<Real>::make("ac", PhaseType::ABC);
  auto localPositive = SimNode<Real>::make("localPositive", PhaseType::DC);
  auto localNegative = SimNode<Real>::make("localNegative", PhaseType::DC);
  auto positiveCableMid =
      SimNode<Real>::make("positiveCableMid", PhaseType::DC);
  auto negativeCableMid =
      SimNode<Real>::make("negativeCableMid", PhaseType::DC);
  auto remotePositive = SimNode<Real>::make("remotePositive", PhaseType::DC);
  auto remoteNegative = SimNode<Real>::make("remoteNegative", PhaseType::DC);
  auto positiveRcNode = SimNode<Real>::make("positiveRcNode", PhaseType::DC);
  auto negativeRcNode = SimNode<Real>::make("negativeRcNode", PhaseType::DC);

  const Real initialPositivePoleVoltage = nominalDcVoltage / 2.0;
  const Real initialNegativePoleVoltage = -nominalDcVoltage / 2.0;
  acNode->setInitialVoltage(Complex(acVoltage / RMS3PH_TO_PEAK1PH, 0.0));
  localPositive->setInitialVoltage(Complex(initialPositivePoleVoltage, 0.0));
  localNegative->setInitialVoltage(Complex(initialNegativePoleVoltage, 0.0));
  positiveCableMid->setInitialVoltage(Complex(initialPositivePoleVoltage, 0.0));
  negativeCableMid->setInitialVoltage(Complex(initialNegativePoleVoltage, 0.0));
  remotePositive->setInitialVoltage(Complex(initialPositivePoleVoltage, 0.0));
  remoteNegative->setInitialVoltage(Complex(initialNegativePoleVoltage, 0.0));
  positiveRcNode->setInitialVoltage(Complex(initialPositivePoleVoltage, 0.0));
  negativeRcNode->setInitialVoltage(Complex(initialNegativePoleVoltage, 0.0));

  auto acSource = EMT::Ph3::NetworkInjection::make("acSource");
  acSource->setParameters(balancedVoltageReference(acVoltage), frequency);
  acSource->connect({acNode});

  auto positiveLineLocal = EMT::DC::SSN::PiLine::make("positiveLineLocal");
  auto positiveLineRemote = EMT::DC::SSN::PiLine::make("positiveLineRemote");
  auto negativeLineLocal = EMT::DC::SSN::PiLine::make("negativeLineLocal");
  auto negativeLineRemote = EMT::DC::SSN::PiLine::make("negativeLineRemote");

  for (const auto &line : {positiveLineLocal, positiveLineRemote,
                           negativeLineLocal, negativeLineRemote}) {
    line->setParameters(cableSectionResistance, cableSectionInductance, 0.0,
                        0.0, 0.0);
    // Explicit validation-study setting. The component default remains the
    // mathematical trapezoidal reference theta=0.5.
    line->setTheta(theta);
  }

  // Positive-pole current orientation: local converter -> remote source.
  positiveLineLocal->connect({positiveCableMid, localPositive});
  positiveLineRemote->connect({remotePositive, positiveCableMid});

  // Negative-pole current orientation: remote source -> local converter.
  negativeLineLocal->connect({localNegative, negativeCableMid});
  negativeLineRemote->connect({negativeCableMid, remoteNegative});

  // Reference-model pole-to-ground branches:
  //
  //   localPositive -- Rg -- positiveRcNode -- Cg -- GND
  //   localNegative -- Rg -- negativeRcNode -- Cg -- GND
  //
  // There is deliberately no artificial pole-to-pole DC-link capacitor.
  auto positiveGroundResistor =
      EMT::Ph1::Resistor::make("positiveGroundResistor");
  auto negativeGroundResistor =
      EMT::Ph1::Resistor::make("negativeGroundResistor");
  positiveGroundResistor->setParameters(groundingResistance);
  negativeGroundResistor->setParameters(groundingResistance);
  positiveGroundResistor->connect({positiveRcNode, localPositive});
  negativeGroundResistor->connect({localNegative, negativeRcNode});

  auto positiveGroundCapacitor =
      EMT::DC::SSN::Capacitor::make("positiveGroundCapacitor");
  auto negativeGroundCapacitor =
      EMT::DC::SSN::Capacitor::make("negativeGroundCapacitor");
  positiveGroundCapacitor->setParameters(groundingCapacitance);
  negativeGroundCapacitor->setParameters(groundingCapacitance);
  positiveGroundCapacitor->setTheta(theta);
  negativeGroundCapacitor->setTheta(theta);
  positiveGroundCapacitor->connect({SimNode<Real>::GND, positiveRcNode});
  negativeGroundCapacitor->connect({negativeRcNode, SimNode<Real>::GND});
  // This single-station validation follows the validated Phase-1 boundary:
  // ideal remote pole sources represent the second DC system. Their symmetric
  // small reference steps exercise both differential-voltage directions
  // without introducing an unvalidated second converter.
  auto remotePositiveSource =
      EMT::DC::VoltageSource::make("remotePositiveSource");
  auto remoteNegativeSource =
      EMT::DC::VoltageSource::make("remoteNegativeSource");
  remotePositiveSource->setParameters(nominalDcVoltage / 2.0);
  remoteNegativeSource->setParameters(nominalDcVoltage / 2.0);
  remotePositiveSource->connect({SimNode<Real>::GND, remotePositive});
  remoteNegativeSource->connect({remoteNegative, SimNode<Real>::GND});

  auto mmc = EMT::Ph3::SSN_MMC::make("mmc");
  mmc->setParameters(frequency, acVoltage, nominalDcVoltage, armInductance,
                     armResistance, submoduleCapacitance, numberOfSubmodules,
                     reactorInductance, reactorResistance);
  mmc->setInitialAngle(0.0);
  // Retain the validated MMC plant PLL for its internal dq-dependent
  // electrical controls, except in the explicit PLL-isolation mode.
  mmc->setPLL(0.001103374, 0.00073, mode != DiagnosticMode::NoPll);

  const Bool useInternalDcVoltageControl =
      mode == DiagnosticMode::BoundaryAudit ||
      mode == DiagnosticMode::InternalPlant ||
      mode == DiagnosticMode::NoEnergy ||
      mode == DiagnosticMode::NoCirculatingDq ||
      mode == DiagnosticMode::NoZeroSequence || mode == DiagnosticMode::NoPll;

  if (useInternalDcVoltageControl) {
    // Diagnostic gains retained from the previous DC-voltage-control test.
    // These are not derived from the 50 nF cable-grounding capacitors.
    constexpr Real kpDc = 4.80800266984177e-3;
    constexpr Real kiDc = 4.53143575980451e-2;
    mmc->setDcVoltageControl(nominalDcVoltage, kpDc, kiDc);
    SPDLOG_INFO("Internal MMC DC-voltage controller: Kp={}, Ki={}", kpDc, kiDc);
  } else {
    mmc->setActiveControlOpenLoop(0.0);
  }

  mmc->setReactiveControlOpenLoop(0.0);
  mmc->setEnergyController(120.0, 400.0, mode != DiagnosticMode::NoEnergy);
  mmc->setOutputCurrentController(117.93, 8.5e4);

  switch (mode) {
  case DiagnosticMode::NoCirculatingDq:
    mmc->setCirculatingCurrentController(0.0, 0.0);
    break;

  case DiagnosticMode::CirculatingDqProportionalOnly:
    mmc->setCirculatingCurrentController(19.93, 0.0);
    break;

  case DiagnosticMode::CirculatingDqIntegralOnly:
    mmc->setCirculatingCurrentController(0.0, 4500.0);
    break;

  default:
    mmc->setCirculatingCurrentController(19.93, 4500.0);
    break;
  }
  mmc->setCirculatingCurrentReferences(0.0, 0.0, 0.0);
  mmc->setLimits(1000.0, 1000.0, 2.0);
  mmc->setOperatingPointInitialization(true, 50, 1e-8);
  mmc->setEigenvalueDiagnostics(true, 5000);
  mmc->setDiagnosticTimeStep(timeStep);
  mmc->setTheta(theta);
  mmc->connect({acNode, localPositive, localNegative});

  auto station = EMT::Ph3::SSN_MMCStation::make("station", mmc);
  EMT::Ph3::SSN_MMCStationParameters parameters;
  parameters.nominalPower = nominalPower;
  parameters.nominalAcLineLineRms = acVoltage;
  parameters.nominalDcVoltage = nominalDcVoltage;
  parameters.nominalFrequencyHz = frequency;
  parameters.controllerTimeStep = timeStep;

  // Match the station current-controller feedforward to the actual MMC
  // differential-current plant impedance. SSN_MMCStation divides these
  // values by two internally.
  const Real vBase = std::sqrt(2.0 / 3.0) * acVoltage;
  const Real iBase = (2.0 / 3.0) * nominalPower / vBase;
  const Real zBase = vBase / iBase;
  const Real rEqAc = armResistance / 2.0 + reactorResistance;
  const Real lEqAc = armInductance / 2.0 + reactorInductance;
  parameters.armResistancePu = 2.0 * rEqAc / zBase;
  parameters.armInductancePu = 2.0 * omega * lEqAc / zBase;
  if (mode == DiagnosticMode::ProportionalVdc)
    parameters.dcVoltageIntegralGain = 0.0;
  station->setParameters(parameters);
  station->setControlMode(
      EMT::Ph3::SSN_MMCStation::ControlMode::DCVoltageReactivePower);
  if (mode == DiagnosticMode::PlantOnly || mode == DiagnosticMode::InnerLoop ||
      mode == DiagnosticMode::GainSign)
    station->setOuterLoopsEnabled(false);
  station->mAngle->set(0.0);
  station->mAngularFrequency->set(omega);

  SystemTopology system(
      frequency,
      SystemNodeList{acNode, localPositive, localNegative, positiveCableMid,
                     negativeCableMid, remotePositive, remoteNegative,
                     positiveRcNode, negativeRcNode},
      SystemComponentList{
          acSource, positiveLineLocal, positiveLineRemote, negativeLineLocal,
          negativeLineRemote, positiveGroundResistor, negativeGroundResistor,
          positiveGroundCapacitor, negativeGroundCapacitor,
          remotePositiveSource, remoteNegativeSource, mmc, station});
  Simulation sim("EMT_SSN_MMCStation_DCPiLine", Logger::Level::off);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.initialize();
  if (mode == DiagnosticMode::BoundaryAudit) {
    const auto stateNames = mmc->getLocalStateNames();
    mmc->setControlSource(
        EMT::Ph3::SSN_MMC::ControlSource::InternalControllers);
    const BoundarySnapshot internal = snapshot(mmc);

    // Transparent shadow: copy the selected final command into the external
    // attribute but deliberately leave InternalControllers selected.
    mmc->setExternalDifferentialVoltageCommand(internal.differential(0, 0),
                                               internal.differential(1, 0));
    const BoundarySnapshot shadow = snapshot(mmc);

    mmc->setControlSource(
        EMT::Ph3::SSN_MMC::ControlSource::ExternalDifferentialVoltage);
    const BoundarySnapshot external = snapshot(mmc);

    SPDLOG_INFO(
        "Boundary internal-shadow differences: differential abs={} rel={}, "
        "common abs={} rel={}, modulation abs={} rel={}, realized abs={} "
        "rel={}, derivative abs={} rel={}, A abs={} rel={}, B abs={} rel={}, "
        "C abs={} rel={}, D abs={} rel={}",
        (internal.differential - shadow.differential).cwiseAbs().maxCoeff(),
        relativeDifference(internal.differential, shadow.differential),
        (internal.common - shadow.common).cwiseAbs().maxCoeff(),
        relativeDifference(internal.common, shadow.common),
        (internal.modulation - shadow.modulation).cwiseAbs().maxCoeff(),
        relativeDifference(internal.modulation, shadow.modulation),
        (internal.realized - shadow.realized).cwiseAbs().maxCoeff(),
        relativeDifference(internal.realized, shadow.realized),
        (internal.derivative - shadow.derivative).cwiseAbs().maxCoeff(),
        relativeDifference(internal.derivative, shadow.derivative),
        (internal.A - shadow.A).cwiseAbs().maxCoeff(),
        relativeDifference(internal.A, shadow.A),
        (internal.B - shadow.B).cwiseAbs().maxCoeff(),
        relativeDifference(internal.B, shadow.B),
        (internal.C - shadow.C).cwiseAbs().maxCoeff(),
        relativeDifference(internal.C, shadow.C),
        (internal.D - shadow.D).cwiseAbs().maxCoeff(),
        relativeDifference(internal.D, shadow.D));

    SPDLOG_INFO(
        "Boundary internal-external differences: differential abs={} rel={}, "
        "common abs={} rel={}, modulation abs={} rel={}, realized abs={} "
        "rel={}, all_derivative abs={} rel={}, electrical_derivative abs={} "
        "rel={}, A abs={} rel={}, B abs={} rel={}, C abs={} rel={}, D abs={} "
        "rel={}",
        (internal.differential - external.differential).cwiseAbs().maxCoeff(),
        relativeDifference(internal.differential, external.differential),
        (internal.common - external.common).cwiseAbs().maxCoeff(),
        relativeDifference(internal.common, external.common),
        (internal.modulation - external.modulation).cwiseAbs().maxCoeff(),
        relativeDifference(internal.modulation, external.modulation),
        (internal.realized - external.realized).cwiseAbs().maxCoeff(),
        relativeDifference(internal.realized, external.realized),
        (internal.derivative - external.derivative).cwiseAbs().maxCoeff(),
        relativeDifference(internal.derivative, external.derivative),
        (internal.derivative.topRows(12) - external.derivative.topRows(12))
            .cwiseAbs()
            .maxCoeff(),
        relativeDifference(internal.derivative.topRows(12),
                           external.derivative.topRows(12)),
        (internal.A - external.A).cwiseAbs().maxCoeff(),
        relativeDifference(internal.A, external.A),
        (internal.B - external.B).cwiseAbs().maxCoeff(),
        relativeDifference(internal.B, external.B),
        (internal.C - external.C).cwiseAbs().maxCoeff(),
        relativeDifference(internal.C, external.C),
        (internal.D - external.D).cwiseAbs().maxCoeff(),
        relativeDifference(internal.D, external.D));

    const Matrix x = mmc->getState();
    const Matrix f = external.derivative;
    Real maximumNormalizedResidual = 0.0;
    String maximumResidualState = "none";
    for (Eigen::Index index = 0; index < x.rows(); ++index) {
      const String &name = stateNames[static_cast<std::size_t>(index)];
      if (name == "grid_angle") {
        SPDLOG_INFO("Cyclic state={} x={} derivative={} rad/s; excluded from "
                    "equilibrium residual.",
                    name, x(index, 0), f(index, 0));
        continue;
      }

      Real base = 1.0;
      if (name.rfind("i", 0) == 0)
        base = 1000.0;
      else if (name.rfind("vC", 0) == 0 || name.rfind("filter_v", 0) == 0)
        base = nominalDcVoltage;
      else if (name.rfind("filter_p", 0) == 0 || name.rfind("filter_q", 0) == 0)
        base = nominalPower;
      const Real normalized =
          std::abs(f(index, 0)) / std::max(base, std::abs(x(index, 0)));
      if (normalized > maximumNormalizedResidual) {
        maximumNormalizedResidual = normalized;
        maximumResidualState = name;
      }
      SPDLOG_INFO(
          "Equilibrium state={} x={} derivative={} base={} abs_residual={} "
          "relative_rate={}",
          name, x(index, 0), f(index, 0), base, std::abs(f(index, 0)),
          normalized);
    }

    const Real iPositiveLocal = positiveLineLocal->intfCurrent()(0, 0);
    const Real iPositiveRemote = positiveLineRemote->intfCurrent()(0, 0);
    const Real iNegativeLocal = negativeLineLocal->intfCurrent()(0, 0);
    const Real iNegativeRemote = negativeLineRemote->intfCurrent()(0, 0);
    const Real iPositiveRc =
        (localPositive->voltage()(0, 0) - positiveRcNode->voltage()(0, 0)) /
        groundingResistance;
    const Real iNegativeRc =
        (negativeRcNode->voltage()(0, 0) - localNegative->voltage()(0, 0)) /
        groundingResistance;
    const Real iPositiveCap = positiveGroundCapacitor->intfCurrent()(0, 0);
    const Real iNegativeCap = negativeGroundCapacitor->intfCurrent()(0, 0);
    const Real iRemotePositiveSource =
        remotePositiveSource->intfCurrent()(0, 0);
    const Real iRemoteNegativeSource =
        remoteNegativeSource->intfCurrent()(0, 0);

    const Real localPositiveKcl =
        iPositiveLocal + iPositiveRc + mmc->getInterfaceCurrent()(3, 0);
    const Real localNegativeKcl =
        -iNegativeLocal - iNegativeRc + mmc->getInterfaceCurrent()(4, 0);
    const Real positiveMidKcl = -iPositiveLocal + iPositiveRemote;
    const Real negativeMidKcl = iNegativeLocal - iNegativeRemote;
    const Real positiveRcKcl = -iPositiveRc + iPositiveCap;
    const Real negativeRcKcl = iNegativeRc - iNegativeCap;
    const Real remotePositiveKcl = -iPositiveRemote + iRemotePositiveSource;
    const Real remoteNegativeKcl = iNegativeRemote - iRemoteNegativeSource;

    const Real dIPositiveLocal =
        (localPositive->voltage()(0, 0) - positiveCableMid->voltage()(0, 0) -
         cableSectionResistance * iPositiveLocal) /
        cableSectionInductance;
    const Real dIPositiveRemote =
        (positiveCableMid->voltage()(0, 0) - remotePositive->voltage()(0, 0) -
         cableSectionResistance * iPositiveRemote) /
        cableSectionInductance;
    const Real dINegativeLocal =
        (negativeCableMid->voltage()(0, 0) - localNegative->voltage()(0, 0) -
         cableSectionResistance * iNegativeLocal) /
        cableSectionInductance;
    const Real dINegativeRemote =
        (remoteNegative->voltage()(0, 0) - negativeCableMid->voltage()(0, 0) -
         cableSectionResistance * iNegativeRemote) /
        cableSectionInductance;

    const Real dVPositiveRc = iPositiveCap / groundingCapacitance;
    // The negative capacitor is oriented GND -> negativeRcNode, so its
    // interface current is -C*d(v_negativeRc)/dt.
    const Real dVNegativeRc = -iNegativeCap / groundingCapacitance;

    const auto includeNetworkResidual = [&](const String &name, Real derivative,
                                            Real base) {
      const Real normalized = std::abs(derivative) / base;
      SPDLOG_INFO("Equilibrium state={} derivative={} base={} relative_rate={}",
                  name, derivative, base, normalized);
      if (normalized > maximumNormalizedResidual) {
        maximumNormalizedResidual = normalized;
        maximumResidualState = name;
      }
    };

    includeNetworkResidual("positive_line_local_current", dIPositiveLocal,
                           1000.0);
    includeNetworkResidual("positive_line_remote_current", dIPositiveRemote,
                           1000.0);
    includeNetworkResidual("negative_line_local_current", dINegativeLocal,
                           1000.0);
    includeNetworkResidual("negative_line_remote_current", dINegativeRemote,
                           1000.0);
    includeNetworkResidual("positive_ground_capacitor_voltage", dVPositiveRc,
                           nominalDcVoltage);
    includeNetworkResidual("negative_ground_capacitor_voltage", dVNegativeRc,
                           nominalDcVoltage);

    SPDLOG_INFO(
        "Boundary RC-network KCL [local+,local-,mid+,mid-,rc+,rc-,remote+,"
        "remote-]=[{},{},{},{},{},{},{},{}] A",
        localPositiveKcl, localNegativeKcl, positiveMidKcl, negativeMidKcl,
        positiveRcKcl, negativeRcKcl, remotePositiveKcl, remoteNegativeKcl);
    SPDLOG_INFO("Maximum complete normalized equilibrium residual={} 1/s at {}",
                maximumNormalizedResidual, maximumResidualState);

    // The previous two-state analytical augmentation assumed a direct
    // pole-to-pole capacitor. It is intentionally not reused for this RC
    // topology because the terminal voltage is algebraically coupled through
    // both 100-ohm series branches. MMC-local A/B/C/D diagnostics above remain
    // valid; a complete RC-network eigenanalysis must augment the MNA
    // companion states explicitly.
    SPDLOG_INFO("RC-network boundary audit: old pole-to-pole-capacitor coupled "
                "linearization intentionally disabled.");

    const Real dLineEnergy =
        cableSectionInductance *
        (iPositiveLocal * dIPositiveLocal + iPositiveRemote * dIPositiveRemote +
         iNegativeLocal * dINegativeLocal + iNegativeRemote * dINegativeRemote);
    const Real dCapEnergy = positiveRcNode->voltage()(0, 0) * iPositiveCap -
                            negativeRcNode->voltage()(0, 0) * iNegativeCap;
    const Real lineLoss =
        cableSectionResistance *
        (iPositiveLocal * iPositiveLocal + iPositiveRemote * iPositiveRemote +
         iNegativeLocal * iNegativeLocal + iNegativeRemote * iNegativeRemote);
    const Real groundingLoss =
        groundingResistance *
        (iPositiveRc * iPositiveRc + iNegativeRc * iNegativeRc);
    const Real remoteVdc =
        remotePositive->voltage()(0, 0) - remoteNegative->voltage()(0, 0);
    const Real sourcePower = remotePositive->voltage()(0, 0) * iPositiveRemote -
                             remoteNegative->voltage()(0, 0) * iNegativeRemote;
    const Real initialVdc = **mmc->dcVoltageAttribute();
    const Real initialIdc = **mmc->dcCurrentAttribute();
    const Real mmcDcPower = initialVdc * initialIdc;
    const Real dExternalEnergy = dLineEnergy + dCapEnergy;
    const Real externalResidual =
        dExternalEnergy + lineLoss + groundingLoss - mmcDcPower + sourcePower;

    const Real energyScale = 3.0 * submoduleCapacitance /
                             (2.0 * static_cast<Real>(numberOfSubmodules));
    const Real mmcCapacitorEnergyDerivative =
        2.0 * energyScale *
        (x(5, 0) * f(5, 0) + x(6, 0) * f(6, 0) + x(7, 0) * f(7, 0) +
         x(8, 0) * f(8, 0) + x(9, 0) * f(9, 0) + x(10, 0) * f(10, 0) +
         2.0 * x(11, 0) * f(11, 0));
    const Real acEquivalentInductance = armInductance / 2.0 + reactorInductance;
    const Real mmcInductorEnergyDerivative =
        1.5 * acEquivalentInductance * (x(0, 0) * f(0, 0) + x(1, 0) * f(1, 0)) +
        3.0 * armInductance * (x(3, 0) * f(3, 0) + x(4, 0) * f(4, 0)) +
        6.0 * armInductance * x(2, 0) * f(2, 0);
    const Real mmcTotalEnergyDerivative =
        mmcCapacitorEnergyDerivative + mmcInductorEnergyDerivative;
    const Real iDeltaSquared = x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0);
    const Real iSigmaSquared = x(3, 0) * x(3, 0) + x(4, 0) * x(4, 0);
    const Real representedLoss =
        armResistance * (6.0 * x(2, 0) * x(2, 0) + 3.0 * iSigmaSquared +
                         0.75 * iDeltaSquared) +
        reactorResistance * 1.5 * iDeltaSquared;
    const Real acPower = **mmc->activePowerAttribute();
    const Real mmcResidual =
        mmcTotalEnergyDerivative + acPower + representedLoss - mmcDcPower;
    SPDLOG_INFO(
        "Continuous energy residuals: external={} W (dE={}, line_loss={}, "
        "ground_loss={}, Pmmc_dc={}, Psource={}, Vremote={}), MMC={} W "
        "(dEcap={}, dEind={}, Pac={}, Pdc={}, represented_loss={})",
        externalResidual, dExternalEnergy, lineLoss, groundingLoss, mmcDcPower,
        sourcePower, remoteVdc, mmcResidual, mmcCapacitorEnergyDerivative,
        mmcInductorEnergyDerivative, acPower, mmcDcPower, representedLoss);
    return 0;
  }
  if (mode == DiagnosticMode::PlantOnly) {
    // Station initialization has calculated and installed the held
    // operating-point command. Select it explicitly while all station
    // controller dynamics remain disabled in Ready.
    const Matrix internalCommand = **mmc->appliedDifferentialVoltageAttribute();
    mmc->setExternalDifferentialVoltageCommand(internalCommand(0, 0),
                                               internalCommand(1, 0));
    mmc->setControlSource(
        EMT::Ph3::SSN_MMC::ControlSource::ExternalDifferentialVoltage);
  } else if (mode == DiagnosticMode::InternalPlant ||
             mode == DiagnosticMode::InternalOpenLoop ||
             mode == DiagnosticMode::NoEnergy ||
             mode == DiagnosticMode::NoCirculatingDq ||
             mode == DiagnosticMode::CirculatingDqProportionalOnly ||
             mode == DiagnosticMode::CirculatingDqIntegralOnly ||
             mode == DiagnosticMode::NoZeroSequence ||
             mode == DiagnosticMode::NoPll) {
    require(mmc->controlSource() ==
                EMT::Ph3::SSN_MMC::ControlSource::InternalControllers,
            "Internal-controller diagnostic did not retain internal control.");
  } else if (mode == DiagnosticMode::InnerLoop ||
             mode == DiagnosticMode::GainSign) {
    station->setCurrentReferences(**station->mIdPu, **station->mIqPu);
    require(station->requestEnable(2e-5, 1.0),
            "RC cable station enable rejected; diagnostic=" +
                std::to_string(**station->mEnableDiagnostic));
  } else {
    require(station->requestEnable(2e-5, 1.0),
            "RC cable station enable rejected; diagnostic=" +
                std::to_string(**station->mEnableDiagnostic));
  }

  SPDLOG_INFO(
      "RC cable diagnostic enable: mode={}, local_Vdc_error={} V, "
      "remote_Vdc_error={} V, line_current={} A, Pac={} W, Pdc={} W, "
      "required_id_ref={} pu, first_id_ref={} pu, command=[{},{}] V",
      static_cast<Int>(mode), **mmc->dcVoltageAttribute() - nominalDcVoltage,
      remotePositive->voltage()(0, 0) - remoteNegative->voltage()(0, 0) -
          nominalDcVoltage,
      positiveLineLocal->intfCurrent()(0, 0), **mmc->activePowerAttribute(),
      **station->mDcPower, **station->mIdPu, **station->mIdReferencePu,
      (**station->mPlantDifferentialVoltageCommand)(0, 0),
      (**station->mPlantDifferentialVoltageCommand)(1, 0));
  SPDLOG_INFO("RC cable diagnostic theta={}", theta);
  require(std::abs(**station->mIdReferencePu - **station->mIdPu) < 2e-5,
          "Initialized Vdc PI output does not reproduce equilibrium id_ref.");
  SPDLOG_INFO("RC cable Vdc PI initialization: error={} pu, P={} pu, I={} pu, "
              "raw_output={} pu, mapped_id_ref={} pu",
              **station->mDcOuterErrorPu,
              **station->mDcOuterProportionalContribution,
              **station->mDcOuterIntegralContribution,
              **station->mDcOuterUnsaturatedOutput, **station->mIdReferencePu);

  const auto internalIdReferenceAttribute =
      mmc->attributeTyped<Real>("i_delta_d_ref");
  const auto internalFilteredVdcAttribute =
      mmc->attributeTyped<Real>("vdc_filtered");

  Real maxLocalPositiveKcl = 0.0;
  Real maxLocalNegativeKcl = 0.0;
  Real maxPositiveMidKcl = 0.0;
  Real maxNegativeMidKcl = 0.0;
  Real maxPositiveRcKcl = 0.0;
  Real maxNegativeRcKcl = 0.0;
  Real maxRemotePositiveKcl = 0.0;
  Real maxRemoteNegativeKcl = 0.0;
  Real minLocalVdc = std::numeric_limits<Real>::infinity();
  Real maxLocalVdc = -std::numeric_limits<Real>::infinity();
  Real minRemoteVdc = std::numeric_limits<Real>::infinity();
  Real maxRemoteVdc = -std::numeric_limits<Real>::infinity();
  Real minEnergy = std::numeric_limits<Real>::infinity();
  Real maxEnergy = -std::numeric_limits<Real>::infinity();
  Real finalLocalVdc = nominalDcVoltage;
  Real finalRemoteVdc = nominalDcVoltage;
  Real positiveQ = -std::numeric_limits<Real>::infinity();
  Real finalQ = 0.0;
  std::deque<TraceSample> trace;
  Bool earlyAbort = false;
  Real positiveStartVdc = nominalDcVoltage;
  Real positiveEndVdc = nominalDcVoltage;
  Real negativeStartVdc = nominalDcVoltage;
  Real negativeEndVdc = nominalDcVoltage;
  Real previousDcStoredEnergy = groundingCapacitance *
                                (nominalDcVoltage / 2.0) *
                                (nominalDcVoltage / 2.0);
  Real maximumPowerEnergyResidual = 0.0;
  Real previousMmcStoredEnergy =
      calculateMmcTotalEnergy(mmc, armInductance, reactorInductance);
  Real previousAcPower = **mmc->activePowerAttribute();
  Real previousMmcDcPower = nominalDcVoltage * **mmc->dcCurrentAttribute();
  Real previousSourcePower = 0.0;
  Real previousExternalLoss = 0.0;
  Real previousMmcLoss = 0.0;
  Real maximumExternalContinuousResidual = 0.0;
  Real maximumMmcContinuousResidual = 0.0;
  Real maximumTotalContinuousResidual = 0.0;
  Real maximumExternalDiscreteResidual = 0.0;
  Real maximumMmcDiscreteResidual = 0.0;
  Real maximumTotalDiscreteResidual = 0.0;
  Real earlyExternalDiscreteResidual = 0.0;
  Real earlyMmcDiscreteResidual = 0.0;
  Real earlyAggregateDiscreteResidual = 0.0;

  sim.start();
  while (sim.time() < sim.finalTime()) {
    station->mAngle->set(omega * sim.time());
    if (mode == DiagnosticMode::GainSign) {
      if (sim.time() >= 0.02 && sim.time() < 0.024)
        station->mIdReferencePu->set(0.001);
      else if (sim.time() >= 0.06 && sim.time() < 0.064)
        station->mIdReferencePu->set(-0.001);
      else
        station->mIdReferencePu->set(0.0);
    }
    Real remotePoleVoltage = nominalDcVoltage / 2.0;
    if (vdcLoopMode && sim.time() >= 0.03 && sim.time() < 0.05)
      remotePoleVoltage += 50.0;
    else if (vdcLoopMode && sim.time() >= 0.08 && sim.time() < 0.10)
      remotePoleVoltage -= 50.0;
    remotePositiveSource->mVoltageRef->set(remotePoleVoltage);
    remoteNegativeSource->mVoltageRef->set(remotePoleVoltage);
    if (mode == DiagnosticMode::FullPi && sim.time() >= 0.23 &&
        sim.time() < 0.25)
      station->mReactivePowerReferencePu->set(0.01);
    else
      station->mReactivePowerReferencePu->set(0.0);
    sim.step();

    require(mmc->getState().allFinite(),
            "RC cable MMC state became non-finite.");
    const Real positiveLocalCurrent = positiveLineLocal->intfCurrent()(0, 0);
    const Real positiveRemoteCurrent = positiveLineRemote->intfCurrent()(0, 0);
    const Real negativeLocalCurrent = negativeLineLocal->intfCurrent()(0, 0);
    const Real negativeRemoteCurrent = negativeLineRemote->intfCurrent()(0, 0);
    const Real positiveRcCurrent =
        (localPositive->voltage()(0, 0) - positiveRcNode->voltage()(0, 0)) /
        groundingResistance;

    const Real negativeRcCurrent =
        (negativeRcNode->voltage()(0, 0) - localNegative->voltage()(0, 0)) /
        groundingResistance;
    const Real positiveCapCurrent =
        positiveGroundCapacitor->intfCurrent()(0, 0);
    const Real negativeCapCurrent =
        negativeGroundCapacitor->intfCurrent()(0, 0);
    const Real remotePositiveSourceCurrent =
        remotePositiveSource->intfCurrent()(0, 0);
    const Real remoteNegativeSourceCurrent =
        remoteNegativeSource->intfCurrent()(0, 0);

    const Real localPositiveKcl = positiveLocalCurrent + positiveRcCurrent +
                                  mmc->getInterfaceCurrent()(3, 0);
    const Real localNegativeKcl = -negativeLocalCurrent - negativeRcCurrent +
                                  mmc->getInterfaceCurrent()(4, 0);
    const Real positiveMidKcl = -positiveLocalCurrent + positiveRemoteCurrent;
    const Real negativeMidKcl = negativeLocalCurrent - negativeRemoteCurrent;
    const Real positiveRcKcl = -positiveRcCurrent + positiveCapCurrent;
    const Real negativeRcKcl = negativeRcCurrent - negativeCapCurrent;
    const Real remotePositiveKcl =
        -positiveRemoteCurrent + remotePositiveSourceCurrent;
    const Real remoteNegativeKcl =
        negativeRemoteCurrent - remoteNegativeSourceCurrent;

    maxLocalPositiveKcl =
        std::max(maxLocalPositiveKcl, std::abs(localPositiveKcl));
    maxLocalNegativeKcl =
        std::max(maxLocalNegativeKcl, std::abs(localNegativeKcl));
    maxPositiveMidKcl = std::max(maxPositiveMidKcl, std::abs(positiveMidKcl));
    maxNegativeMidKcl = std::max(maxNegativeMidKcl, std::abs(negativeMidKcl));
    maxPositiveRcKcl = std::max(maxPositiveRcKcl, std::abs(positiveRcKcl));
    maxNegativeRcKcl = std::max(maxNegativeRcKcl, std::abs(negativeRcKcl));
    maxRemotePositiveKcl =
        std::max(maxRemotePositiveKcl, std::abs(remotePositiveKcl));
    maxRemoteNegativeKcl =
        std::max(maxRemoteNegativeKcl, std::abs(remoteNegativeKcl));

    const Real localVdc = **mmc->dcVoltageAttribute();
    const Real remoteVdc =
        remotePositive->voltage()(0, 0) - remoteNegative->voltage()(0, 0);
    minLocalVdc = std::min(minLocalVdc, localVdc);
    maxLocalVdc = std::max(maxLocalVdc, localVdc);
    minRemoteVdc = std::min(minRemoteVdc, remoteVdc);
    maxRemoteVdc = std::max(maxRemoteVdc, remoteVdc);

    const Real mmcCapacitorEnergy = **mmc->storedEnergyAttribute();
    const Real mmcTotalEnergy =
        calculateMmcTotalEnergy(mmc, armInductance, reactorInductance);
    const Real positiveRcVoltage = positiveRcNode->voltage()(0, 0);
    const Real negativeRcVoltage = negativeRcNode->voltage()(0, 0);
    const Real rcCapEnergy = 0.5 * groundingCapacitance *
                             (positiveRcVoltage * positiveRcVoltage +
                              negativeRcVoltage * negativeRcVoltage);
    const Real lineEnergy = 0.5 * cableSectionInductance *
                            (positiveLocalCurrent * positiveLocalCurrent +
                             positiveRemoteCurrent * positiveRemoteCurrent +
                             negativeLocalCurrent * negativeLocalCurrent +
                             negativeRemoteCurrent * negativeRemoteCurrent);
    const Real dcStoredEnergy = rcCapEnergy + lineEnergy;
    const Real storedEnergyDerivative =
        (dcStoredEnergy - previousDcStoredEnergy) / timeStep;
    previousDcStoredEnergy = dcStoredEnergy;

    const Real lineLoss = cableSectionResistance *
                          (positiveLocalCurrent * positiveLocalCurrent +
                           positiveRemoteCurrent * positiveRemoteCurrent +
                           negativeLocalCurrent * negativeLocalCurrent +
                           negativeRemoteCurrent * negativeRemoteCurrent);
    const Real groundingLoss =
        groundingResistance * (positiveRcCurrent * positiveRcCurrent +
                               negativeRcCurrent * negativeRcCurrent);
    const Real externalLoss = lineLoss + groundingLoss;
    const Real sourcePower =
        remotePositive->voltage()(0, 0) * positiveRemoteCurrent -
        remoteNegative->voltage()(0, 0) * negativeRemoteCurrent;
    const Real converterDcPower = localVdc * **mmc->dcCurrentAttribute();
    const Real externalPowerEnergyResidual =
        storedEnergyDerivative + externalLoss - converterDcPower + sourcePower;
    maximumPowerEnergyResidual = std::max(
        maximumPowerEnergyResidual, std::abs(externalPowerEnergyResidual));
    const Matrix mmcState = mmc->getState();
    const Real iDeltaSquared =
        mmcState(0, 0) * mmcState(0, 0) + mmcState(1, 0) * mmcState(1, 0);
    const Real iSigmaSquared =
        mmcState(3, 0) * mmcState(3, 0) + mmcState(4, 0) * mmcState(4, 0);
    const Real mmcLoss =
        armResistance * (6.0 * mmcState(2, 0) * mmcState(2, 0) +
                         3.0 * iSigmaSquared + 0.75 * iDeltaSquared) +
        reactorResistance * 1.5 * iDeltaSquared;
    const Real acPower = **mmc->activePowerAttribute();
    const Real mmcStoredEnergyDerivative =
        (mmcTotalEnergy - previousMmcStoredEnergy) / timeStep;
    const Real externalContinuousResidual = externalPowerEnergyResidual;
    const Real mmcContinuousResidual =
        mmcStoredEnergyDerivative + acPower + mmcLoss - converterDcPower;
    const Real totalContinuousResidual =
        externalContinuousResidual + mmcContinuousResidual;
    const Real weightedExternalLoss =
        (1.0 - theta) * previousExternalLoss + theta * externalLoss;
    const Real weightedSourcePower =
        (1.0 - theta) * previousSourcePower + theta * sourcePower;
    const Real weightedMmcDcPower =
        (1.0 - theta) * previousMmcDcPower + theta * converterDcPower;
    const Real weightedAcPower =
        (1.0 - theta) * previousAcPower + theta * acPower;
    const Real weightedMmcLoss =
        (1.0 - theta) * previousMmcLoss + theta * mmcLoss;
    const Real externalDiscreteResidual =
        storedEnergyDerivative + weightedExternalLoss - weightedMmcDcPower +
        weightedSourcePower;
    const Real mmcDiscreteResidual = mmcStoredEnergyDerivative +
                                     weightedAcPower + weightedMmcLoss -
                                     weightedMmcDcPower;
    const Real totalDiscreteResidual =
        externalDiscreteResidual + mmcDiscreteResidual;
    // Skip only the initialization step, whose previous algebraic powers are
    // node guesses rather than a solved MNA sample.
    if (sim.time() > 2.0 * timeStep) {
      maximumExternalContinuousResidual =
          std::max(maximumExternalContinuousResidual,
                   std::abs(externalContinuousResidual));
      maximumMmcContinuousResidual = std::max(maximumMmcContinuousResidual,
                                              std::abs(mmcContinuousResidual));
      maximumTotalContinuousResidual = std::max(
          maximumTotalContinuousResidual, std::abs(totalContinuousResidual));
      maximumExternalDiscreteResidual = std::max(
          maximumExternalDiscreteResidual, std::abs(externalDiscreteResidual));
      maximumMmcDiscreteResidual =
          std::max(maximumMmcDiscreteResidual, std::abs(mmcDiscreteResidual));
      maximumTotalDiscreteResidual = std::max(maximumTotalDiscreteResidual,
                                              std::abs(totalDiscreteResidual));
      if (sim.time() <= 0.05) {
        earlyExternalDiscreteResidual = std::max(
            earlyExternalDiscreteResidual, std::abs(externalDiscreteResidual));
        earlyMmcDiscreteResidual =
            std::max(earlyMmcDiscreteResidual, std::abs(mmcDiscreteResidual));
        earlyAggregateDiscreteResidual = std::max(
            earlyAggregateDiscreteResidual, std::abs(totalDiscreteResidual));
      }
    }
    previousMmcStoredEnergy = mmcTotalEnergy;
    previousAcPower = acPower;
    previousMmcDcPower = converterDcPower;
    previousSourcePower = sourcePower;
    previousExternalLoss = externalLoss;
    previousMmcLoss = mmcLoss;
    minEnergy = std::min(minEnergy, mmcCapacitorEnergy);
    maxEnergy = std::max(maxEnergy, mmcCapacitorEnergy);
    if (sim.time() >= 0.23 && sim.time() < 0.25)
      positiveQ = std::max(positiveQ, **station->mReactivePowerPu);
    finalLocalVdc = localVdc;
    finalRemoteVdc = remoteVdc;
    finalQ = **station->mReactivePowerPu;
    if (mode == DiagnosticMode::GainSign) {
      if (sim.time() >= 0.02 && sim.time() < 0.02005)
        positiveStartVdc = localVdc;
      if (sim.time() >= 0.02395 && sim.time() < 0.02405)
        positiveEndVdc = localVdc;
      if (sim.time() >= 0.06 && sim.time() < 0.06005)
        negativeStartVdc = localVdc;
      if (sim.time() >= 0.06395 && sim.time() < 0.06405)
        negativeEndVdc = localVdc;
    }

    const Matrix plantModulation = **mmc->appliedModulationAttribute();
    trace.push_back(
        {sim.time(),
         localVdc,
         remoteVdc,
         positiveLocalCurrent,
         **station->mFilteredDcVoltage - **station->mDcVoltageReference,
         **station->mIdReferencePu,
         **station->mIdPu,
         **station->mVdReferencePu,
         **station->mVqReferencePu,
         **internalFilteredVdcAttribute,
         **internalIdReferenceAttribute,
         stateValueByName(mmc, "xi_active"),
         **mmc->dcCurrentAttribute(),
         positiveRcCurrent,
         negativeRcCurrent,
         mmcCapacitorEnergy,
         mmcTotalEnergy,
         **station->mDcOuterProportionalContribution,
         **station->mDcOuterIntegralContribution,
         plantModulation(0, 0),
         plantModulation(1, 0),
         rcCapEnergy,
         lineEnergy,
         sourcePower,
         externalPowerEnergyResidual,
         mmcContinuousResidual});
    if (trace.size() > 50)
      trace.pop_front();
    const Real nominalEnergy = 14.52e6;
    if (std::abs(localVdc) > 1.5 * nominalDcVoltage ||
        std::abs(**mmc->dcCurrentAttribute()) >
            1.5 * nominalPower / nominalDcVoltage ||
        mmcCapacitorEnergy < 0.25 * nominalEnergy ||
        mmcCapacitorEnergy > 2.0 * nominalEnergy ||
        mmc->getState().cwiseAbs().maxCoeff() > 10.0 * nominalDcVoltage) {
      earlyAbort = true;
      break;
    }
  }
  sim.stop();

  if (earlyAbort) {
    std::cerr << "EARLY ABORT: physically meaningful bound exceeded\n";
    std::cerr << "t,local_vdc,remote_vdc,line_i,station_vdc_error,"
                 "station_id_ref,station_id,station_vd_ref,station_vq_ref,"
                 "internal_vdc_filtered,internal_id_ref,internal_xi_active,"
                 "mmc_idc,positive_rc_i,negative_rc_i,mmc_cap_energy,"
                 "mmc_total_energy,station_vdc_p,station_vdc_i,plant_mod_d,"
                 "plant_mod_q,rc_cap_energy,line_energy,source_power,"
                 "external_energy_residual,mmc_energy_residual\n";
    std::cerr << std::setprecision(12);
    for (const auto &sample : trace)
      std::cerr << sample.time << ',' << sample.localVdc << ','
                << sample.remoteVdc << ',' << sample.lineCurrent << ','
                << sample.stationVdcError << ',' << sample.stationIdReference
                << ',' << sample.stationId << ',' << sample.stationVdReference
                << ',' << sample.stationVqReference << ','
                << sample.internalFilteredVdc << ','
                << sample.internalIdReference << ','
                << sample.internalActiveIntegrator << ',' << sample.mmcDcCurrent
                << ',' << sample.positiveRcCurrent << ','
                << sample.negativeRcCurrent << ',' << sample.mmcCapacitorEnergy
                << ',' << sample.mmcTotalEnergy << ','
                << sample.stationOuterProportional << ','
                << sample.stationOuterIntegral << ',' << sample.plantModulationD
                << ',' << sample.plantModulationQ << ',' << sample.rcCapEnergy
                << ',' << sample.lineEnergy << ',' << sample.sourcePower << ','
                << sample.externalPowerEnergyResidual << ','
                << sample.mmcPowerEnergyResidual << '\n';
  }

  SPDLOG_INFO(
      "Stage3F RC cable: local_Vdc=[{},{}] final={}, remote_Vdc=[{},{}] "
      "final={}, KCL=[{},{},{},{},{},{},{},{}] A, energy=[{},{}] J, "
      "Q_peak={} pu, Q_final={} pu, line_currents=[{},{},{},{}] A, "
      "RC_currents=[{},{}] A",
      minLocalVdc, maxLocalVdc, finalLocalVdc, minRemoteVdc, maxRemoteVdc,
      finalRemoteVdc, maxLocalPositiveKcl, maxLocalNegativeKcl,
      maxPositiveMidKcl, maxNegativeMidKcl, maxPositiveRcKcl, maxNegativeRcKcl,
      maxRemotePositiveKcl, maxRemoteNegativeKcl, minEnergy, maxEnergy,
      positiveQ, finalQ, positiveLineLocal->intfCurrent()(0, 0),
      positiveLineRemote->intfCurrent()(0, 0),
      negativeLineLocal->intfCurrent()(0, 0),
      negativeLineRemote->intfCurrent()(0, 0),
      (localPositive->voltage()(0, 0) - positiveRcNode->voltage()(0, 0)) /
          groundingResistance,
      (negativeRcNode->voltage()(0, 0) - localNegative->voltage()(0, 0)) /
          groundingResistance);
  SPDLOG_INFO("RC cable maximum pre-abort DC power/energy residual={} W",
              maximumPowerEnergyResidual);
  SPDLOG_INFO("RC cable energy audit dt={} theta={}: continuous max "
              "[external,MMC,aggregate]=[{},{},{}] W; discrete-theta max "
              "[external,MMC,aggregate]=[{},{},{}] W",
              timeStep, theta, maximumExternalContinuousResidual,
              maximumMmcContinuousResidual, maximumTotalContinuousResidual,
              maximumExternalDiscreteResidual, maximumMmcDiscreteResidual,
              maximumTotalDiscreteResidual);
  SPDLOG_INFO("RC cable early (t<=50 ms) discrete energy residuals "
              "[external,MMC,aggregate]=[{},{},{}] W",
              earlyExternalDiscreteResidual, earlyMmcDiscreteResidual,
              earlyAggregateDiscreteResidual);
  SPDLOG_INFO(
      "MMC-local finite-difference Jacobian: max_real={} 1/s, "
      "max_continuous_magnitude={} 1/s, discrete_spectral_radius={}, "
      "dominant_frequency={} Hz",
      mmc->attributeTyped<Real>("jacobian_max_real_eigenvalue")->get(),
      mmc->attributeTyped<Real>("jacobian_max_abs_eigenvalue")->get(),
      mmc->attributeTyped<Real>("jacobian_max_discrete_magnitude")->get(),
      mmc->attributeTyped<Real>("jacobian_discrete_dominant_frequency")->get());
  if (mode == DiagnosticMode::GainSign) {
    const Real positiveSlope = (positiveEndVdc - positiveStartVdc) / 0.004;
    const Real negativeSlope = (negativeEndVdc - negativeStartVdc) / 0.004;
    SPDLOG_INFO("RC cable id/Vdc gain: positive_slope={} V/s, "
                "negative_slope={} V/s",
                positiveSlope, negativeSlope);
    require(positiveSlope * negativeSlope < 0.0,
            "RC cable id perturbations did not produce opposite Vdc slopes.");
  }

  require(!earlyAbort, "RC cable diagnostic exceeded physical bounds.");

  require(maxLocalPositiveKcl < 1e-6 && maxLocalNegativeKcl < 1e-6 &&
              maxPositiveMidKcl < 1e-6 && maxNegativeMidKcl < 1e-6 &&
              maxPositiveRcKcl < 1e-6 && maxNegativeRcKcl < 1e-6 &&
              maxRemotePositiveKcl < 1e-6 && maxRemoteNegativeKcl < 1e-6,
          "RC cable KCL residual exceeds 1e-6 A.");

  if (vdcLoopMode) {
    require(maxRemoteVdc > nominalDcVoltage + 1.0 &&
                minRemoteVdc < nominalDcVoltage - 1.0,
            "RC cable disturbances did not exercise both Vdc directions.");
  }

  if (vdcLoopMode) {
    const Bool dcVoltageRestored =
        std::abs(finalLocalVdc - nominalDcVoltage) < 200.0 &&
        std::abs(finalRemoteVdc - nominalDcVoltage) < 200.0;

    if (!dcVoltageRestored) {
      SPDLOG_WARN("RC cable DC voltage did not restore: "
                  "finalLocalVdc={} V, finalRemoteVdc={} V, nominalVdc={} V",
                  finalLocalVdc, finalRemoteVdc, nominalDcVoltage);
    }
  } else
    require(minLocalVdc > 0.9 * nominalDcVoltage &&
                maxLocalVdc < 1.1 * nominalDcVoltage,
            "Isolated RC cable plant did not remain bounded within 10% Vdc.");
  if (mode == DiagnosticMode::FullPi)
    require(positiveQ > 0.001 && std::abs(finalQ) < 0.002,
            "RC cable reactive-power step or recovery failed.");
  require(minEnergy > 0.0 && maxEnergy / minEnergy < 1.1,
          "RC cable MMC energy is unbounded.");
  require(std::abs(positiveLineLocal->intfCurrent()(0, 0) -
                   positiveLineRemote->intfCurrent()(0, 0)) < 0.1 &&
              std::abs(negativeLineLocal->intfCurrent()(0, 0) -
                       negativeLineRemote->intfCurrent()(0, 0)) < 0.1 &&
              std::abs(positiveLineRemote->intfCurrent()(0, 0) -
                       negativeLineRemote->intfCurrent()(0, 0)) < 0.1,
          "RC cable section or pole-current continuity failed.");
  return 0;
}
