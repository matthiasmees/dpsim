// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Eigenvalues>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_Ph1_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_RxLine.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>

using namespace CPS;
using namespace DPsim;

namespace {

void require(Bool condition, const String &message) {
  if (!condition)
    throw std::runtime_error(message);
}

Matrix diagonalThreePhase(Real value) { return value * Matrix::Identity(3, 3); }

/// DPsim's EMT source reference and SimNode positive-sequence initialization
/// use the line-line RMS phasor. The analytical network calculation below uses
/// phase-to-neutral RMS phasors.
Complex phaseRmsToLineLineRms(const Complex &phaseVoltage) {
  return std::sqrt(3.0) * phaseVoltage;
}

MatrixComp balancedLineLineRmsReference(const Complex &phaseA) {
  MatrixComp reference(3, 1);
  reference << phaseA, phaseA * SHIFT_TO_PHASE_B, phaseA * SHIFT_TO_PHASE_C;
  return reference;
}

struct HarmonyHotStartPoint {
  Complex sourceIdealPhaseRms;
  Complex sourceBusPhaseRms;
  Complex mmc1PhaseRms;
  Complex mmc2PhaseRms;
  Complex loadBusPhaseRms;
  Complex loadInternalPhaseRms;

  Complex sourceCurrentPhaseRms;
  Complex loadCurrentPhaseRms;

  Complex mmc1AcPower;
  Complex mmc2AcPower;
  Complex loadPower;

  Real dcCurrent;
  Real mmc1DcVoltage;
  Real mmc2DcVoltage;
  Real mmc1DcPower;
  Real mmc2DcPower;
  Real mmc1RepresentedLoss;
  Real mmc2RepresentedLoss;
};

/// Calculate a physically consistent equilibrium for Harmony's point-to-point
/// network after removing both converter transformers.
///
/// The unchanged receiving-side network fixes Q/P. Consequently Harmony's
/// original MMC2 reference of -20 Mvar cannot coexist with P=-50 MW once the
/// transformer impedances are removed. The reactive-power reference is derived
/// from the remaining 5+j110-ohm branch and the 2280-ohm/1.457-H load.
HarmonyHotStartPoint makeHarmonyHotStartPoint(
    Real frequency, Real receivingActivePower, Real receivingDcVoltage,
    Real sourceResistance, Real acLineResistance, Real acLineReactance,
    Real loadResistance, Real loadInductance, Real dcPoleResistance,
    Real armResistance, Real reactorResistance) {
  require(receivingActivePower > 0.0,
          "Receiving-side active power must be positive.");
  require(receivingDcVoltage > 0.0,
          "Receiving-side DC voltage must be positive.");

  const Real omega = 2.0 * PI * frequency;
  const Complex zAcLine(acLineResistance, acLineReactance);
  const Complex zLoad(loadResistance, omega * loadInductance);
  const Complex zReceiving = zAcLine + zLoad;
  const Complex yReceiving = 1.0 / zReceiving;

  require(yReceiving.real() > 0.0,
          "Receiving-side equivalent conductance must be positive.");

  // Select the transformer-free MMC2 terminal voltage such that the unchanged
  // line and load absorb exactly 50 MW.
  const Real mmc2PhaseRmsMagnitude =
      std::sqrt(receivingActivePower / (3.0 * yReceiving.real()));
  const Complex mmc2PhaseRms(mmc2PhaseRmsMagnitude, 0.0);
  const Complex loadCurrent = mmc2PhaseRms / zReceiving;
  const Complex loadBusPhaseRms = loadCurrent * zLoad;
  const Complex loadInternalPhaseRms =
      loadCurrent * Complex(0.0, omega * loadInductance);

  const Complex receivingNetworkPower =
      3.0 * mmc2PhaseRms * std::conj(loadCurrent);
  const Complex loadPower = 3.0 * loadBusPhaseRms * std::conj(loadCurrent);

  // Converter terminal powers use the SSN_MMC convention: positive Pac enters
  // the converter. MMC2 supplies the passive receiving network.
  const Complex mmc2AcPower = -receivingNetworkPower;

  // Include the losses represented explicitly by the averaged MMC model in the
  // hot-start power balance:
  //   Pdc = Pac + Ploss
  // with i_dc = 3*iSigma_z.
  const Real vPeakPhase = std::sqrt(2.0) * mmc2PhaseRmsMagnitude;
  const Real iDeltaD2 = (2.0 / 3.0) * mmc2AcPower.real() / vPeakPhase;
  const Real iDeltaQ2 = -(2.0 / 3.0) * mmc2AcPower.imag() / vPeakPhase;

  const Real differentialCurrentLossCoefficient =
      0.75 * armResistance + 1.5 * reactorResistance;
  const Real zeroSequenceLossCoefficient = (2.0 / 3.0) * armResistance;
  const Real mmc2DifferentialLoss = differentialCurrentLossCoefficient *
                                    (iDeltaD2 * iDeltaD2 + iDeltaQ2 * iDeltaQ2);

  // At MMC2:
  //   -Vdc2*Idc = -Preceiving + Ploss2.
  // This gives a scalar quadratic because the zero-sequence arm loss is
  // (2/3)*Rarm*Idc^2.
  const Real quadraticA = zeroSequenceLossCoefficient;
  const Real quadraticB = receivingDcVoltage;
  const Real quadraticC = -receivingActivePower + mmc2DifferentialLoss;
  const Real discriminant =
      quadraticB * quadraticB - 4.0 * quadraticA * quadraticC;
  require(discriminant > 0.0,
          "No positive DC-current solution for the hot-start point.");

  const Real dcCurrent =
      (-quadraticB + std::sqrt(discriminant)) / (2.0 * quadraticA);
  require(dcCurrent > 0.0, "Calculated DC current must be positive.");

  const Real dcLoopResistance = 2.0 * dcPoleResistance;
  const Real mmc1DcVoltage = receivingDcVoltage + dcLoopResistance * dcCurrent;
  const Real zeroSequenceLoss =
      zeroSequenceLossCoefficient * dcCurrent * dcCurrent;
  const Real mmc2RepresentedLoss = zeroSequenceLoss + mmc2DifferentialLoss;

  // Keep the same converter-side AC-voltage magnitude on both ends. MMC1 is
  // initialized at Q=0. Its active power follows from its DC power minus its
  // represented loss. Since the differential-current loss depends on P1^2,
  // solve the resulting scalar quadratic exactly.
  const Complex mmc1PhaseRms = mmc2PhaseRms;
  const Real pToIDelta = (2.0 / 3.0) / vPeakPhase;
  const Real mmc1PowerSquaredLossCoefficient =
      differentialCurrentLossCoefficient * pToIDelta * pToIDelta;
  const Real availableMmc1AcPower =
      mmc1DcVoltage * dcCurrent - zeroSequenceLoss;
  const Real mmc1PowerDiscriminant =
      1.0 + 4.0 * mmc1PowerSquaredLossCoefficient * availableMmc1AcPower;
  require(mmc1PowerDiscriminant > 0.0,
          "No positive MMC1 active-power solution for the hot start.");

  const Real mmc1ActivePower = (-1.0 + std::sqrt(mmc1PowerDiscriminant)) /
                               (2.0 * mmc1PowerSquaredLossCoefficient);
  const Complex mmc1AcPower(mmc1ActivePower, 0.0);
  const Real mmc1RepresentedLoss =
      zeroSequenceLoss +
      mmc1PowerSquaredLossCoefficient * mmc1ActivePower * mmc1ActivePower;

  // Adapt the ideal source voltage so that the calculated MMC1 operating point
  // is obtained through Harmony's 5-ohm source impedance and 5+j110-ohm branch.
  const Complex sourceCurrent = std::conj(mmc1AcPower / (3.0 * mmc1PhaseRms));
  const Complex sourceBusPhaseRms = mmc1PhaseRms + zAcLine * sourceCurrent;
  const Complex sourceIdealPhaseRms =
      sourceBusPhaseRms + sourceResistance * sourceCurrent;

  HarmonyHotStartPoint point;
  point.sourceIdealPhaseRms = sourceIdealPhaseRms;
  point.sourceBusPhaseRms = sourceBusPhaseRms;
  point.mmc1PhaseRms = mmc1PhaseRms;
  point.mmc2PhaseRms = mmc2PhaseRms;
  point.loadBusPhaseRms = loadBusPhaseRms;
  point.loadInternalPhaseRms = loadInternalPhaseRms;
  point.sourceCurrentPhaseRms = sourceCurrent;
  point.loadCurrentPhaseRms = loadCurrent;
  point.mmc1AcPower = mmc1AcPower;
  point.mmc2AcPower = mmc2AcPower;
  point.loadPower = loadPower;
  point.dcCurrent = dcCurrent;
  point.mmc1DcVoltage = mmc1DcVoltage;
  point.mmc2DcVoltage = receivingDcVoltage;
  point.mmc1DcPower = mmc1DcVoltage * dcCurrent;
  point.mmc2DcPower = -receivingDcVoltage * dcCurrent;
  point.mmc1RepresentedLoss = mmc1RepresentedLoss;
  point.mmc2RepresentedLoss = mmc2RepresentedLoss;
  return point;
}

void configureCommonHarmonyMmc(const EMT::Ph3::SSN_MMC::Ptr &mmc,
                               Real frequency, Real harmonyAcVoltageParameter,
                               Real nominalDcVoltage, Real armInductance,
                               Real armResistance, Real submoduleCapacitance,
                               UInt numberOfSubmodules, Real reactorInductance,
                               Real reactorResistance, Real initialActivePower,
                               Real initialReactivePower, Real initialDcCurrent,
                               Real theta, Real timeStep) {
  mmc->setParameters(frequency, harmonyAcVoltageParameter, nominalDcVoltage,
                     armInductance, armResistance, submoduleCapacitance,
                     numberOfSubmodules, reactorInductance, reactorResistance);
  mmc->setInitialAngle(0.0);
  mmc->setInitialOperatingPoint(initialActivePower, initialReactivePower);

  // Harmony point-to-point inner controllers.
  mmc->setPLL(0.001103374, 0.00073, true);
  mmc->setEnergyController(120.0, 400.0, true);
  mmc->setZeroSequenceCurrentController(19.93, 4500.0);
  mmc->setOutputCurrentController(117.93, 8.5e4);
  mmc->setCirculatingCurrentController(19.93, 4500.0);

  // At equilibrium Idc = 3*iSigma_z. The sign follows each converter's
  // DC-port current orientation.
  const Real initialSigmaZReference = initialDcCurrent / 3.0;
  mmc->setCirculatingCurrentReferences(0.0, 0.0, initialSigmaZReference);

  mmc->setLimits(1000.0, 1000.0, 2.0);
  mmc->setOperatingPointInitialization(true, 100, 1e-9);
  mmc->setDiagnosticTimeStep(timeStep);
  mmc->setTheta(theta);
}

Real maximumAbsoluteValue(const Matrix &matrix) {
  if (matrix.size() == 0)
    return 0.0;
  return matrix.cwiseAbs().maxCoeff();
}

struct MmcDiagnosticSnapshot {
  Matrix state;
  Matrix derivative;
  Matrix interfaceVoltage;
  Matrix interfaceCurrent;
  Matrix differentialVoltage;
  Matrix commonModeVoltage;
  Matrix modulation;
  Matrix realizedVoltage;

  Real vdc = 0.0;
  Real idc = 0.0;
  Real activePower = 0.0;
  Real reactivePower = 0.0;
  Real storedEnergy = 0.0;
  Real dcPower = 0.0;
  Real powerBalanceError = 0.0;
  Real stateNorm = 0.0;
  Real derivativeNorm = 0.0;
  Real equilibriumResidual = 0.0;
  Real pllError = 0.0;
  Real pllAngleDeviation = 0.0;
  Real gridAngle = 0.0;
  Real pllFrequency = 0.0;
  Real gridVoltageD = 0.0;
  Real gridVoltageQ = 0.0;
  Real controlVoltageD = 0.0;
  Real controlVoltageQ = 0.0;
  Real deltaCurrentD = 0.0;
  Real deltaCurrentQ = 0.0;
  Real sigmaCurrentZ = 0.0;
  Real deltaCurrentReferenceD = 0.0;
  Real deltaCurrentReferenceQ = 0.0;
  Real sigmaCurrentReferenceZ = 0.0;

  String largestStateName = "none";
  String largestDerivativeName = "none";
  Real largestStateMagnitude = 0.0;
  Real largestDerivativeMagnitude = 0.0;
};

Real scalarAttribute(const EMT::Ph3::SSN_MMC::Ptr &mmc, const String &name) {
  return mmc->attributeTyped<Real>(name)->get();
}

MmcDiagnosticSnapshot captureMmcSnapshot(const EMT::Ph3::SSN_MMC::Ptr &mmc,
                                         Bool evaluateDerivative) {
  MmcDiagnosticSnapshot snapshot;
  snapshot.state = mmc->getState();
  snapshot.derivative =
      evaluateDerivative
          ? mmc->getStateDerivative()
          : Matrix::Constant(snapshot.state.rows(), 1,
                             std::numeric_limits<Real>::quiet_NaN());
  snapshot.interfaceVoltage = mmc->getInterfaceVoltage();
  snapshot.interfaceCurrent = mmc->getInterfaceCurrent();
  snapshot.differentialVoltage = **mmc->appliedDifferentialVoltageAttribute();
  snapshot.commonModeVoltage = **mmc->appliedCommonModeVoltageAttribute();
  snapshot.modulation = **mmc->appliedModulationAttribute();
  snapshot.realizedVoltage = **mmc->realizedConverterVoltageAttribute();

  snapshot.vdc = **mmc->dcVoltageAttribute();
  snapshot.idc = **mmc->dcCurrentAttribute();
  snapshot.activePower = **mmc->activePowerAttribute();
  snapshot.reactivePower = **mmc->reactivePowerAttribute();
  snapshot.storedEnergy = **mmc->storedEnergyAttribute();
  snapshot.dcPower = scalarAttribute(mmc, "p_dc");
  snapshot.powerBalanceError = scalarAttribute(mmc, "power_balance_error");
  snapshot.stateNorm = scalarAttribute(mmc, "state_norm");
  snapshot.derivativeNorm = scalarAttribute(mmc, "state_derivative_norm");
  snapshot.equilibriumResidual =
      scalarAttribute(mmc, "equilibrium_residual_norm");
  snapshot.pllError = scalarAttribute(mmc, "pll_error");
  snapshot.pllAngleDeviation = scalarAttribute(mmc, "pll_angle_deviation");
  snapshot.gridAngle = scalarAttribute(mmc, "grid_angle");
  snapshot.pllFrequency = scalarAttribute(mmc, "pll_frequency");
  snapshot.gridVoltageD = scalarAttribute(mmc, "v_grid_d");
  snapshot.gridVoltageQ = scalarAttribute(mmc, "v_grid_q");
  snapshot.controlVoltageD = scalarAttribute(mmc, "v_control_d");
  snapshot.controlVoltageQ = scalarAttribute(mmc, "v_control_q");
  snapshot.deltaCurrentD = scalarAttribute(mmc, "i_delta_d");
  snapshot.deltaCurrentQ = scalarAttribute(mmc, "i_delta_q");
  snapshot.sigmaCurrentZ = scalarAttribute(mmc, "i_sigma_z");
  snapshot.deltaCurrentReferenceD = scalarAttribute(mmc, "i_delta_d_ref");
  snapshot.deltaCurrentReferenceQ = scalarAttribute(mmc, "i_delta_q_ref");
  snapshot.sigmaCurrentReferenceZ = scalarAttribute(mmc, "i_sigma_z_ref");

  const auto stateNames = mmc->getLocalStateNames();
  if (snapshot.state.rows() > 0) {
    Eigen::Index row = 0;
    Eigen::Index column = 0;
    snapshot.largestStateMagnitude =
        snapshot.state.cwiseAbs().maxCoeff(&row, &column);
    if (row >= 0 && static_cast<std::size_t>(row) < stateNames.size())
      snapshot.largestStateName = stateNames[static_cast<std::size_t>(row)];
  }

  if (snapshot.derivative.rows() > 0 && snapshot.derivative.allFinite()) {
    Eigen::Index row = 0;
    Eigen::Index column = 0;
    snapshot.largestDerivativeMagnitude =
        snapshot.derivative.cwiseAbs().maxCoeff(&row, &column);
    if (row >= 0 && static_cast<std::size_t>(row) < stateNames.size())
      snapshot.largestDerivativeName =
          stateNames[static_cast<std::size_t>(row)];
  } else {
    snapshot.largestDerivativeMagnitude = std::numeric_limits<Real>::infinity();
    snapshot.largestDerivativeName = "non-finite";
  }

  return snapshot;
}

struct P2PDiagnosticSample {
  Real time = 0.0;
  MmcDiagnosticSnapshot mmc1;
  MmcDiagnosticSnapshot mmc2;
  Real positiveDcLineCurrent = 0.0;
  Real negativeDcLineCurrent = 0.0;
  Real dcLineLoss = 0.0;
  Real dcPowerResidual = 0.0;
};

P2PDiagnosticSample
captureP2PSample(Real time, const EMT::Ph3::SSN_MMC::Ptr &mmc1,
                 const EMT::Ph3::SSN_MMC::Ptr &mmc2,
                 const EMT::Ph1::Resistor::Ptr &positiveDcLine,
                 const EMT::Ph1::Resistor::Ptr &negativeDcLine,
                 Real dcPoleResistance, Bool evaluateDerivative) {
  P2PDiagnosticSample sample;
  sample.time = time;
  sample.mmc1 = captureMmcSnapshot(mmc1, evaluateDerivative);
  sample.mmc2 = captureMmcSnapshot(mmc2, evaluateDerivative);
  sample.positiveDcLineCurrent = positiveDcLine->intfCurrent()(0, 0);
  sample.negativeDcLineCurrent = negativeDcLine->intfCurrent()(0, 0);
  sample.dcLineLoss =
      dcPoleResistance *
      (sample.positiveDcLineCurrent * sample.positiveDcLineCurrent +
       sample.negativeDcLineCurrent * sample.negativeDcLineCurrent);
  sample.dcPowerResidual =
      sample.mmc1.dcPower + sample.mmc2.dcPower - sample.dcLineLoss;
  return sample;
}

void printLocalJacobianDiagnostic(const String &name,
                                  const EMT::Ph3::SSN_MMC::Ptr &mmc) {
  Matrix A;
  Matrix B;
  Matrix C;
  Matrix D;
  mmc->getLocalLinearization(A, B, C, D);

  Eigen::EigenSolver<Matrix> solver(A, false);
  require(solver.info() == Eigen::Success,
          name + " local eigenvalue calculation failed.");

  Real maximumRealPart = -std::numeric_limits<Real>::infinity();
  Real maximumMagnitude = 0.0;
  Real dominantFrequency = 0.0;
  for (Eigen::Index index = 0; index < solver.eigenvalues().rows(); ++index) {
    const Complex eigenvalue = solver.eigenvalues()(index);
    if (eigenvalue.real() > maximumRealPart) {
      maximumRealPart = eigenvalue.real();
      dominantFrequency = std::abs(eigenvalue.imag()) / (2.0 * PI);
    }
    maximumMagnitude = std::max(maximumMagnitude, std::abs(eigenvalue));
  }

  // Numerical differentiation updates command diagnostics. Restore them at
  // the unperturbed state before the simulation starts.
  (void)mmc->getStateDerivative();

  SPDLOG_INFO(
      "{} initial local Jacobian: max_real={} 1/s, max_magnitude={} 1/s, "
      "dominant_frequency={} Hz",
      name, maximumRealPart, maximumMagnitude, dominantFrequency);
}

void printCompactDiagnostic(const P2PDiagnosticSample &sample,
                            Real initialMmc1Energy, Real initialMmc2Energy) {
  SPDLOG_INFO(
      "P2P diagnostic t={} s: "
      "MMC1[Vdc={},Idc={},P={},Q={},dE={},pll_err={},"
      "max_x={}({}),max_dx={}({}),|m|max={}] "
      "MMC2[Vdc={},Idc={},P={},Q={},dE={},pll_err={},"
      "max_x={}({}),max_dx={}({}),|m|max={}] "
      "DC[I+={},I-={},loss={},power_residual={}]",
      sample.time, sample.mmc1.vdc, sample.mmc1.idc, sample.mmc1.activePower,
      sample.mmc1.reactivePower, sample.mmc1.storedEnergy - initialMmc1Energy,
      sample.mmc1.pllError, sample.mmc1.largestStateMagnitude,
      sample.mmc1.largestStateName, sample.mmc1.largestDerivativeMagnitude,
      sample.mmc1.largestDerivativeName,
      maximumAbsoluteValue(sample.mmc1.modulation), sample.mmc2.vdc,
      sample.mmc2.idc, sample.mmc2.activePower, sample.mmc2.reactivePower,
      sample.mmc2.storedEnergy - initialMmc2Energy, sample.mmc2.pllError,
      sample.mmc2.largestStateMagnitude, sample.mmc2.largestStateName,
      sample.mmc2.largestDerivativeMagnitude, sample.mmc2.largestDerivativeName,
      maximumAbsoluteValue(sample.mmc2.modulation),
      sample.positiveDcLineCurrent, sample.negativeDcLineCurrent,
      sample.dcLineLoss, sample.dcPowerResidual);
}

void writeMmcColumns(std::ofstream &stream, const String &prefix,
                     const std::vector<String> &stateNames) {
  stream << ',' << prefix << "_vdc" << ',' << prefix << "_idc" << ',' << prefix
         << "_p" << ',' << prefix << "_q" << ',' << prefix << "_energy" << ','
         << prefix << "_pdc" << ',' << prefix << "_power_balance_error" << ','
         << prefix << "_pll_error" << ',' << prefix << "_pll_angle_deviation"
         << ',' << prefix << "_grid_angle" << ',' << prefix << "_pll_frequency"
         << ',' << prefix << "_v_grid_d" << ',' << prefix << "_v_grid_q" << ','
         << prefix << "_v_control_d" << ',' << prefix << "_v_control_q" << ','
         << prefix << "_i_delta_d" << ',' << prefix << "_i_delta_q" << ','
         << prefix << "_i_sigma_z" << ',' << prefix << "_i_delta_d_ref" << ','
         << prefix << "_i_delta_q_ref" << ',' << prefix << "_i_sigma_z_ref"
         << ',' << prefix << "_modulation_max" << ',' << prefix
         << "_differential_voltage_max" << ',' << prefix
         << "_common_voltage_max" << ',' << prefix << "_realized_voltage_max";
  for (const auto &stateName : stateNames)
    stream << ',' << prefix << "_x_" << stateName;
  for (const auto &stateName : stateNames)
    stream << ',' << prefix << "_dx_" << stateName;
}

void writeMmcValues(std::ofstream &stream,
                    const MmcDiagnosticSnapshot &snapshot) {
  stream << ',' << snapshot.vdc << ',' << snapshot.idc << ','
         << snapshot.activePower << ',' << snapshot.reactivePower << ','
         << snapshot.storedEnergy << ',' << snapshot.dcPower << ','
         << snapshot.powerBalanceError << ',' << snapshot.pllError << ','
         << snapshot.pllAngleDeviation << ',' << snapshot.gridAngle << ','
         << snapshot.pllFrequency << ',' << snapshot.gridVoltageD << ','
         << snapshot.gridVoltageQ << ',' << snapshot.controlVoltageD << ','
         << snapshot.controlVoltageQ << ',' << snapshot.deltaCurrentD << ','
         << snapshot.deltaCurrentQ << ',' << snapshot.sigmaCurrentZ << ','
         << snapshot.deltaCurrentReferenceD << ','
         << snapshot.deltaCurrentReferenceQ << ','
         << snapshot.sigmaCurrentReferenceZ << ','
         << maximumAbsoluteValue(snapshot.modulation) << ','
         << maximumAbsoluteValue(snapshot.differentialVoltage) << ','
         << maximumAbsoluteValue(snapshot.commonModeVoltage) << ','
         << maximumAbsoluteValue(snapshot.realizedVoltage);
  for (Eigen::Index row = 0; row < snapshot.state.rows(); ++row)
    stream << ',' << snapshot.state(row, 0);
  for (Eigen::Index row = 0; row < snapshot.derivative.rows(); ++row)
    stream << ',' << snapshot.derivative(row, 0);
}

void writeFailureTrace(const String &fileName,
                       const std::deque<P2PDiagnosticSample> &trace,
                       const std::vector<String> &mmc1StateNames,
                       const std::vector<String> &mmc2StateNames) {
  std::ofstream stream(fileName);
  if (!stream)
    throw std::runtime_error("Could not open diagnostic trace file: " +
                             fileName);

  stream << std::setprecision(17);
  stream << "time,positive_dc_line_current,negative_dc_line_current,"
            "dc_line_loss,dc_power_residual";
  writeMmcColumns(stream, "mmc1", mmc1StateNames);
  writeMmcColumns(stream, "mmc2", mmc2StateNames);
  stream << '\n';

  for (const auto &sample : trace) {
    stream << sample.time << ',' << sample.positiveDcLineCurrent << ','
           << sample.negativeDcLineCurrent << ',' << sample.dcLineLoss << ','
           << sample.dcPowerResidual;
    writeMmcValues(stream, sample.mmc1);
    writeMmcValues(stream, sample.mmc2);
    stream << '\n';
  }
}

void dumpMmcStateToConsole(const String &name,
                           const MmcDiagnosticSnapshot &snapshot,
                           const std::vector<String> &stateNames) {
  std::cerr << "\n" << name << " final safe state\n";
  std::cerr << "index,name,x,dx\n";
  std::cerr << std::setprecision(17);
  const Eigen::Index count =
      std::min(snapshot.state.rows(), snapshot.derivative.rows());
  for (Eigen::Index row = 0; row < count; ++row) {
    const String stateName = static_cast<std::size_t>(row) < stateNames.size()
                                 ? stateNames[static_cast<std::size_t>(row)]
                                 : "unknown";
    std::cerr << row << ',' << stateName << ',' << snapshot.state(row, 0) << ','
              << snapshot.derivative(row, 0) << '\n';
  }

  std::cerr << name
            << " interface_voltage=" << snapshot.interfaceVoltage.transpose()
            << '\n';
  std::cerr << name
            << " interface_current=" << snapshot.interfaceCurrent.transpose()
            << '\n';
  std::cerr << name << " differential_voltage="
            << snapshot.differentialVoltage.transpose() << '\n';
  std::cerr << name
            << " common_mode_voltage=" << snapshot.commonModeVoltage.transpose()
            << '\n';
  std::cerr << name << " modulation=" << snapshot.modulation.transpose()
            << '\n';
  std::cerr << name
            << " realized_voltage=" << snapshot.realizedVoltage.transpose()
            << '\n';
}

} // namespace

int main(int argc, char **argv) {
  const Real timeStep = argc > 1 ? std::stod(argv[1]) : 20e-6;
  const Real finalTime = argc > 2 ? std::stod(argv[2]) : 0.5;
  const Real theta = argc > 3 ? std::stod(argv[3]) : 0.5;
  const UInt reportIntervalSteps =
      argc > 4 ? static_cast<UInt>(std::stoul(argv[4])) : 500;

  require(timeStep > 0.0, "The EMT time step must be positive.");
  require(finalTime > 0.0, "The final time must be positive.");
  require(theta >= 0.5 && theta <= 1.0,
          "The SSN integration theta must be in [0.5, 1.0].");
  require(reportIntervalSteps > 0,
          "The diagnostic report interval must be positive.");

  const String simName = "EMT_SSN_MMC_Harmony_P2P_HotStart";
  Logger::setLogDir("logs/" + simName);

  const Real frequency = 50.0;
  const Real omega = 2.0 * PI * frequency;

  // Harmony point-to-point data. Harmony's V_m=345 kV is retained as the MMC
  // model parameter. The actual ideal-source voltage is calculated below for
  // the transformer-free loaded equilibrium.
  const Real harmonyAcVoltageParameter = 345e3;
  const Real nominalDcVoltage = 440e3;
  const Real receivingActivePower = 50e6;

  const Real sourceResistance = 5.0;
  const Real acLineResistance = 5.0;
  const Real acLineReactance = 110.0;
  const Real acLineInductance = acLineReactance / omega;
  const Real loadResistance = 2.28e3;
  const Real loadInductance = 1.457;
  const Real dcPoleResistance = 20.0;

  const Real armInductance = 0.05;
  const Real armResistance = 1.07;
  const Real submoduleCapacitance = 0.01;
  const UInt numberOfSubmodules = 400;
  const Real reactorInductance = 0.0005;
  const Real reactorResistance = 0.0001;

  // Harmony outer-loop gains.
  const Real activePowerKp = 6.6667e-7;
  const Real activePowerKi = 3.3333e-4;
  const Real reactivePowerKp = 6.6667e-7;
  const Real reactivePowerKi = 3.3333e-4;
  const Real dcVoltageKp = 2.0;
  const Real dcVoltageKi = 82.0;

  const HarmonyHotStartPoint operatingPoint = makeHarmonyHotStartPoint(
      frequency, receivingActivePower, nominalDcVoltage, sourceResistance,
      acLineResistance, acLineReactance, loadResistance, loadInductance,
      dcPoleResistance, armResistance, reactorResistance);

  SPDLOG_INFO(
      "Harmony transformer-free hot start: source_LL={} V angle={} deg, "
      "MMC1_LL={} V, MMC2_LL={} V, load_LL={} V, "
      "MMC1_PQ=[{},{}] W/var, MMC2_PQ=[{},{}] W/var, "
      "Vdc=[{},{}] V, Idc={} A, MMC_losses=[{},{}] W",
      std::abs(phaseRmsToLineLineRms(operatingPoint.sourceIdealPhaseRms)),
      std::arg(operatingPoint.sourceIdealPhaseRms) * 180.0 / PI,
      std::abs(phaseRmsToLineLineRms(operatingPoint.mmc1PhaseRms)),
      std::abs(phaseRmsToLineLineRms(operatingPoint.mmc2PhaseRms)),
      std::abs(phaseRmsToLineLineRms(operatingPoint.loadBusPhaseRms)),
      operatingPoint.mmc1AcPower.real(), operatingPoint.mmc1AcPower.imag(),
      operatingPoint.mmc2AcPower.real(), operatingPoint.mmc2AcPower.imag(),
      operatingPoint.mmc1DcVoltage, operatingPoint.mmc2DcVoltage,
      operatingPoint.dcCurrent, operatingPoint.mmc1RepresentedLoss,
      operatingPoint.mmc2RepresentedLoss);

  // AC topology from Harmony with both converter transformers removed:
  //
  // ideal source -- 5 ohm -- (5+j110 ohm) -- MMC1
  //                                            || bipolar DC 2 x 20 ohm
  //                               MMC2 -- (5+j110 ohm) -- R-L load
  auto sourceIdealNode = SimNode<Real>::make("ACBUS01_IDEAL", PhaseType::ABC);
  auto sourceBusNode = SimNode<Real>::make("ACBUS01", PhaseType::ABC);
  auto mmc1AcNode = SimNode<Real>::make("ACBUS02", PhaseType::ABC);
  auto mmc2AcNode = SimNode<Real>::make("ACBUS05", PhaseType::ABC);
  auto loadBusNode = SimNode<Real>::make("ACBUS06", PhaseType::ABC);
  auto loadInternalNode =
      SimNode<Real>::make("LOAD02_INTERNAL", PhaseType::ABC);

  sourceIdealNode->setInitialVoltage(
      phaseRmsToLineLineRms(operatingPoint.sourceIdealPhaseRms));
  sourceBusNode->setInitialVoltage(
      phaseRmsToLineLineRms(operatingPoint.sourceBusPhaseRms));
  mmc1AcNode->setInitialVoltage(
      phaseRmsToLineLineRms(operatingPoint.mmc1PhaseRms));
  mmc2AcNode->setInitialVoltage(
      phaseRmsToLineLineRms(operatingPoint.mmc2PhaseRms));
  loadBusNode->setInitialVoltage(
      phaseRmsToLineLineRms(operatingPoint.loadBusPhaseRms));
  loadInternalNode->setInitialVoltage(
      phaseRmsToLineLineRms(operatingPoint.loadInternalPhaseRms));

  auto mmc1Positive = SimNode<Real>::make("DCBUS01_POS", PhaseType::DC);
  auto mmc1Negative = SimNode<Real>::make("DCBUS01_NEG", PhaseType::DC);
  auto mmc2Positive = SimNode<Real>::make("DCBUS02_POS", PhaseType::DC);
  auto mmc2Negative = SimNode<Real>::make("DCBUS02_NEG", PhaseType::DC);

  mmc1Positive->setInitialVoltage(
      Complex(operatingPoint.mmc1DcVoltage / 2.0, 0.0));
  mmc1Negative->setInitialVoltage(
      Complex(-operatingPoint.mmc1DcVoltage / 2.0, 0.0));
  mmc2Positive->setInitialVoltage(
      Complex(operatingPoint.mmc2DcVoltage / 2.0, 0.0));
  mmc2Negative->setInitialVoltage(
      Complex(-operatingPoint.mmc2DcVoltage / 2.0, 0.0));

  auto acSource = EMT::Ph3::NetworkInjection::make("SRC01");
  acSource->setParameters(balancedLineLineRmsReference(phaseRmsToLineLineRms(
                              operatingPoint.sourceIdealPhaseRms)),
                          frequency);
  acSource->connect({sourceIdealNode});

  auto sourceImpedance = EMT::Ph3::Resistor::make("SRC01_ZSRC");
  sourceImpedance->setParameters(diagonalThreePhase(sourceResistance));
  sourceImpedance->connect({sourceIdealNode, sourceBusNode});

  auto acLine1 = EMT::Ph3::RxLine::make("br1_ac");
  acLine1->setParameters(diagonalThreePhase(acLineResistance),
                         diagonalThreePhase(acLineInductance));
  acLine1->connect({sourceBusNode, mmc1AcNode});

  auto acLine2 = EMT::Ph3::RxLine::make("br2_ac");
  acLine2->setParameters(diagonalThreePhase(acLineResistance),
                         diagonalThreePhase(acLineInductance));
  acLine2->connect({mmc2AcNode, loadBusNode});

  // Exact Harmony LOAD02: a balanced series R-L-C load with
  // R=2280 ohm, L=1.457 H and C=0.
  auto loadResistor = EMT::Ph3::Resistor::make("LOAD02_R");
  loadResistor->setParameters(diagonalThreePhase(loadResistance));
  loadResistor->connect({loadBusNode, loadInternalNode});

  auto loadInductor = EMT::Ph3::Inductor::make("LOAD02_L");
  loadInductor->setParameters(diagonalThreePhase(loadInductance));
  loadInductor->connect({loadInternalNode, SimNode<Real>::GND});

  // Harmony's two-conductor DC impedance: 20 ohm in each pole.
  auto positiveDcLine = EMT::Ph1::Resistor::make("br1_dc_positive");
  auto negativeDcLine = EMT::Ph1::Resistor::make("br1_dc_negative");
  positiveDcLine->setParameters(dcPoleResistance);
  negativeDcLine->setParameters(dcPoleResistance);
  positiveDcLine->connect({mmc1Positive, mmc2Positive});
  negativeDcLine->connect({mmc2Negative, mmc1Negative});

  auto mmc1 = EMT::Ph3::SSN_MMC::make("MMC1");
  auto mmc2 = EMT::Ph3::SSN_MMC::make("MMC2");

  configureCommonHarmonyMmc(
      mmc1, frequency, harmonyAcVoltageParameter, nominalDcVoltage,
      armInductance, armResistance, submoduleCapacitance, numberOfSubmodules,
      reactorInductance, reactorResistance, operatingPoint.mmc1AcPower.real(),
      operatingPoint.mmc1AcPower.imag(), operatingPoint.dcCurrent, theta,
      timeStep);
  configureCommonHarmonyMmc(
      mmc2, frequency, harmonyAcVoltageParameter, nominalDcVoltage,
      armInductance, armResistance, submoduleCapacitance, numberOfSubmodules,
      reactorInductance, reactorResistance, operatingPoint.mmc2AcPower.real(),
      operatingPoint.mmc2AcPower.imag(), -operatingPoint.dcCurrent, theta,
      timeStep);

  // Exact Harmony control allocation.
  mmc1->setActivePowerControl(operatingPoint.mmc1AcPower.real(), activePowerKp,
                              activePowerKi);
  mmc1->setReactivePowerControl(operatingPoint.mmc1AcPower.imag(),
                                reactivePowerKp, reactivePowerKi);

  mmc2->setDcVoltageControl(operatingPoint.mmc2DcVoltage, dcVoltageKp,
                            dcVoltageKi);
  mmc2->setReactivePowerControl(operatingPoint.mmc2AcPower.imag(),
                                reactivePowerKp, reactivePowerKi);

  mmc1->connect({mmc1AcNode, mmc1Positive, mmc1Negative});
  mmc2->connect({mmc2AcNode, mmc2Positive, mmc2Negative});

  SystemTopology system(
      frequency,
      SystemNodeList{sourceIdealNode, sourceBusNode, mmc1AcNode, mmc2AcNode,
                     loadBusNode, loadInternalNode, mmc1Positive, mmc1Negative,
                     mmc2Positive, mmc2Negative},
      SystemComponentList{acSource, sourceImpedance, acLine1, acLine2,
                          loadResistor, loadInductor, positiveDcLine,
                          negativeDcLine, mmc1, mmc2});

  auto logger = DataLogger::make(simName);
  logger->logAttribute("ACBUS01.v", sourceBusNode->attribute("v"));
  logger->logAttribute("ACBUS02.v", mmc1AcNode->attribute("v"));
  logger->logAttribute("ACBUS05.v", mmc2AcNode->attribute("v"));
  logger->logAttribute("ACBUS06.v", loadBusNode->attribute("v"));
  logger->logAttribute("br1_ac.i", acLine1->attribute("i_intf"));
  logger->logAttribute("br2_ac.i", acLine2->attribute("i_intf"));
  logger->logAttribute("br1_dc_positive.i",
                       positiveDcLine->attribute("i_intf"));
  logger->logAttribute("br1_dc_negative.i",
                       negativeDcLine->attribute("i_intf"));

  logger->logAttribute("MMC1.vdc", mmc1->dcVoltageAttribute());
  logger->logAttribute("MMC2.vdc", mmc2->dcVoltageAttribute());
  logger->logAttribute("MMC1.idc", mmc1->dcCurrentAttribute());
  logger->logAttribute("MMC2.idc", mmc2->dcCurrentAttribute());
  logger->logAttribute("MMC1.p", mmc1->activePowerAttribute());
  logger->logAttribute("MMC2.p", mmc2->activePowerAttribute());
  logger->logAttribute("MMC1.q", mmc1->reactivePowerAttribute());
  logger->logAttribute("MMC2.q", mmc2->reactivePowerAttribute());
  logger->logAttribute("MMC1.energy", mmc1->storedEnergyAttribute());
  logger->logAttribute("MMC2.energy", mmc2->storedEnergyAttribute());
  logger->logAttribute("MMC1.modulation", mmc1->appliedModulationAttribute());
  logger->logAttribute("MMC2.modulation", mmc2->appliedModulationAttribute());
  logger->logAttribute("MMC1.eq_residual",
                       mmc1->attribute("equilibrium_residual_norm"));
  logger->logAttribute("MMC2.eq_residual",
                       mmc2->attribute("equilibrium_residual_norm"));
  logger->logAttribute("MMC1.state_norm", mmc1->attribute("state_norm"));
  logger->logAttribute("MMC2.state_norm", mmc2->attribute("state_norm"));
  logger->logAttribute("MMC1.derivative_norm",
                       mmc1->attribute("state_derivative_norm"));
  logger->logAttribute("MMC2.derivative_norm",
                       mmc2->attribute("state_derivative_norm"));
  logger->logAttribute("MMC1.pdc", mmc1->attribute("p_dc"));
  logger->logAttribute("MMC2.pdc", mmc2->attribute("p_dc"));
  logger->logAttribute("MMC1.power_balance_error",
                       mmc1->attribute("power_balance_error"));
  logger->logAttribute("MMC2.power_balance_error",
                       mmc2->attribute("power_balance_error"));
  logger->logAttribute("MMC1.pll_error", mmc1->attribute("pll_error"));
  logger->logAttribute("MMC2.pll_error", mmc2->attribute("pll_error"));
  logger->logAttribute("MMC1.pll_angle_deviation",
                       mmc1->attribute("pll_angle_deviation"));
  logger->logAttribute("MMC2.pll_angle_deviation",
                       mmc2->attribute("pll_angle_deviation"));
  logger->logAttribute("MMC1.i_delta_d_ref", mmc1->attribute("i_delta_d_ref"));
  logger->logAttribute("MMC1.i_delta_q_ref", mmc1->attribute("i_delta_q_ref"));
  logger->logAttribute("MMC2.i_delta_d_ref", mmc2->attribute("i_delta_d_ref"));
  logger->logAttribute("MMC2.i_delta_q_ref", mmc2->attribute("i_delta_q_ref"));
  logger->logAttribute("MMC1.applied_differential_voltage",
                       mmc1->appliedDifferentialVoltageAttribute());
  logger->logAttribute("MMC2.applied_differential_voltage",
                       mmc2->appliedDifferentialVoltageAttribute());
  logger->logAttribute("MMC1.applied_common_mode_voltage",
                       mmc1->appliedCommonModeVoltageAttribute());
  logger->logAttribute("MMC2.applied_common_mode_voltage",
                       mmc2->appliedCommonModeVoltageAttribute());
  logger->logAttribute("MMC1.realized_converter_voltage",
                       mmc1->realizedConverterVoltageAttribute());
  logger->logAttribute("MMC2.realized_converter_voltage",
                       mmc2->realizedConverterVoltageAttribute());

  Simulation sim(simName, Logger::Level::off);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.addLogger(logger);
  sim.initialize();

  const Real initialMmc1Vdc = **mmc1->dcVoltageAttribute();
  const Real initialMmc2Vdc = **mmc2->dcVoltageAttribute();
  const Real initialMmc1Energy = **mmc1->storedEnergyAttribute();
  const Real initialMmc2Energy = **mmc2->storedEnergyAttribute();
  const Real initialMmc1Power = **mmc1->activePowerAttribute();
  const Real initialMmc2Power = **mmc2->activePowerAttribute();

  SPDLOG_INFO("Initialized Harmony P2P: MMC1 "
              "[P,Q,Vdc,Idc,E,eq_res]=[{},{},{},{},{},{}], "
              "MMC2 [P,Q,Vdc,Idc,E,eq_res]=[{},{},{},{},{},{}]",
              initialMmc1Power, **mmc1->reactivePowerAttribute(),
              initialMmc1Vdc, **mmc1->dcCurrentAttribute(), initialMmc1Energy,
              mmc1->attributeTyped<Real>("equilibrium_residual_norm")->get(),
              initialMmc2Power, **mmc2->reactivePowerAttribute(),
              initialMmc2Vdc, **mmc2->dcCurrentAttribute(), initialMmc2Energy,
              mmc2->attributeTyped<Real>("equilibrium_residual_norm")->get());

  printLocalJacobianDiagnostic("MMC1", mmc1);
  printLocalJacobianDiagnostic("MMC2", mmc2);

  Real maxMmc1VdcDeviation = 0.0;
  Real maxMmc2VdcDeviation = 0.0;
  Real maxMmc1EnergyDeviation = 0.0;
  Real maxMmc2EnergyDeviation = 0.0;

  const auto mmc1StateNames = mmc1->getLocalStateNames();
  const auto mmc2StateNames = mmc2->getLocalStateNames();
  constexpr std::size_t traceCapacity = 500;
  std::deque<P2PDiagnosticSample> trace;

  auto appendTrace = [&](P2PDiagnosticSample sample) {
    trace.push_back(std::move(sample));
    if (trace.size() > traceCapacity)
      trace.pop_front();
  };

  appendTrace(captureP2PSample(0.0, mmc1, mmc2, positiveDcLine, negativeDcLine,
                               dcPoleResistance, true));
  printCompactDiagnostic(trace.back(), initialMmc1Energy, initialMmc2Energy);

  Bool failed = false;
  String failureMessage;
  Real failureTime = 0.0;
  UInt stepCount = 0;

  // No switching, precharge, deblocking, reference steps or artificial test
  // disturbances. The initialized steady operating point is held throughout.
  sim.start();
  while (sim.time() < sim.finalTime()) {
    try {
      sim.step();
      ++stepCount;

      P2PDiagnosticSample sample =
          captureP2PSample(sim.time(), mmc1, mmc2, positiveDcLine,
                           negativeDcLine, dcPoleResistance, true);
      appendTrace(std::move(sample));

      const auto &latest = trace.back();
      if (!latest.mmc1.state.allFinite() || !latest.mmc2.state.allFinite() ||
          !latest.mmc1.derivative.allFinite() ||
          !latest.mmc2.derivative.allFinite()) {
        failed = true;
        failureMessage =
            "A state or independently evaluated derivative became non-finite.";
        failureTime = sim.time();
        break;
      }

      maxMmc1VdcDeviation = std::max(
          maxMmc1VdcDeviation, std::abs(latest.mmc1.vdc - initialMmc1Vdc));
      maxMmc2VdcDeviation = std::max(
          maxMmc2VdcDeviation, std::abs(latest.mmc2.vdc - initialMmc2Vdc));
      maxMmc1EnergyDeviation =
          std::max(maxMmc1EnergyDeviation,
                   std::abs(latest.mmc1.storedEnergy - initialMmc1Energy));
      maxMmc2EnergyDeviation =
          std::max(maxMmc2EnergyDeviation,
                   std::abs(latest.mmc2.storedEnergy - initialMmc2Energy));

      if (stepCount % reportIntervalSteps == 0)
        printCompactDiagnostic(latest, initialMmc1Energy, initialMmc2Energy);
    } catch (const std::exception &exception) {
      failed = true;
      failureMessage = exception.what();
      failureTime = sim.time();
      break;
    }
  }

  try {
    sim.stop();
  } catch (const std::exception &exception) {
    if (!failed) {
      failed = true;
      failureMessage = String("sim.stop() failed: ") + exception.what();
      failureTime = sim.time();
    }
  }

  if (failed) {
    const String traceFile = simName + "_failure_trace.csv";
    try {
      writeFailureTrace(traceFile, trace, mmc1StateNames, mmc2StateNames);
      SPDLOG_CRITICAL(
          "Harmony P2P failed at t={} s after {} successful steps: {}. "
          "The last {} safe samples were written to {}.",
          failureTime, stepCount, failureMessage, trace.size(), traceFile);
    } catch (const std::exception &traceException) {
      SPDLOG_CRITICAL(
          "Harmony P2P failed at t={} s after {} successful steps: {}. "
          "Writing the failure trace also failed: {}",
          failureTime, stepCount, failureMessage, traceException.what());
    }

    if (!trace.empty()) {
      printCompactDiagnostic(trace.back(), initialMmc1Energy,
                             initialMmc2Energy);
      dumpMmcStateToConsole("MMC1", trace.back().mmc1, mmc1StateNames);
      dumpMmcStateToConsole("MMC2", trace.back().mmc2, mmc2StateNames);
    }
    return 2;
  }

  SPDLOG_INFO(
      "Harmony P2P hot-start drift after {} s: "
      "MMC1 [max_dVdc={}, final_dVdc={}, max_dE={}, final_dE={}, final_P={}, "
      "final_Q={}], MMC2 [max_dVdc={}, final_dVdc={}, max_dE={}, final_dE={}, "
      "final_P={}, final_Q={}]",
      finalTime, maxMmc1VdcDeviation,
      **mmc1->dcVoltageAttribute() - initialMmc1Vdc, maxMmc1EnergyDeviation,
      **mmc1->storedEnergyAttribute() - initialMmc1Energy,
      **mmc1->activePowerAttribute(), **mmc1->reactivePowerAttribute(),
      maxMmc2VdcDeviation, **mmc2->dcVoltageAttribute() - initialMmc2Vdc,
      maxMmc2EnergyDeviation,
      **mmc2->storedEnergyAttribute() - initialMmc2Energy,
      **mmc2->activePowerAttribute(), **mmc2->reactivePowerAttribute());

  return 0;
}
