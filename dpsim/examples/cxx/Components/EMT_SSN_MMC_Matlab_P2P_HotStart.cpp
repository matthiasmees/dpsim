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
#include <Eigen/LU>

#include <DPsim.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Capacitor.h>
#include <dpsim-models/EMT/EMT_DC_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_DC_SSN_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>
#include <dpsim-models/EMT/EMT_Ph3_NetworkInjection.h>
#include <dpsim-models/EMT/EMT_Ph3_RXLoad.h>
#include <dpsim-models/EMT/EMT_Ph3_Resistor.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_MMC.h>
#include <dpsim-models/EMT/EMT_Ph3_Transformer.h>
#include <dpsim-models/MNASimPowerComp.h>

using namespace CPS;
using namespace DPsim;

namespace {

void require(Bool condition, const String &message) {
  if (!condition)
    throw std::runtime_error(message);
}

enum class ClosedLoopDiagnosticMode {
  InternalControllers,
  HoldInitialDifferentialVoltage,
  HoldInitialFullConverterVoltage
};

ClosedLoopDiagnosticMode parseClosedLoopDiagnosticMode(const String &value) {
  if (value == "internal")
    return ClosedLoopDiagnosticMode::InternalControllers;
  if (value == "hold-differential")
    return ClosedLoopDiagnosticMode::HoldInitialDifferentialVoltage;
  if (value == "hold-full")
    return ClosedLoopDiagnosticMode::HoldInitialFullConverterVoltage;
  throw std::invalid_argument(
      "Diagnostic mode must be 'internal', 'hold-differential', "
      "or 'hold-full'.");
}

String closedLoopDiagnosticModeName(ClosedLoopDiagnosticMode mode) {
  switch (mode) {
  case ClosedLoopDiagnosticMode::InternalControllers:
    return "internal";
  case ClosedLoopDiagnosticMode::HoldInitialDifferentialVoltage:
    return "hold-differential";
  case ClosedLoopDiagnosticMode::HoldInitialFullConverterVoltage:
    return "hold-full";
  }
  return "unknown";
}

Matrix diagonalThreePhase(Real value) { return value * Matrix::Identity(3, 3); }

Complex phaseRmsToLineLineRms(const Complex &phaseVoltage) {
  return std::sqrt(3.0) * phaseVoltage;
}

MatrixComp balancedLineLineRmsReference(const Complex &lineLinePhaseA) {
  MatrixComp reference(3, 1);
  reference << lineLinePhaseA, lineLinePhaseA * SHIFT_TO_PHASE_B,
      lineLinePhaseA * SHIFT_TO_PHASE_C;
  return reference;
}

/// Inverse amplitude-invariant Park transformation used by SSN_MMC.
/// The d/q quantities and returned abc quantities are instantaneous peak values.
/// This must match SSN_MMC::dqToAbc().
Matrix dqToAbcPeak(Real d, Real q, Real theta) {
  Matrix abc(3, 1);
  abc(0, 0) = d * std::cos(theta) - q * std::sin(theta);
  abc(1, 0) = d * std::cos(theta - 2.0 * PI / 3.0) -
              q * std::sin(theta - 2.0 * PI / 3.0);
  abc(2, 0) = d * std::cos(theta + 2.0 * PI / 3.0) -
              q * std::sin(theta + 2.0 * PI / 3.0);
  return abc;
}

struct MatlabCaseParameters {
  Real frequency = 50.0;
  Real ratedPower = 1000e6;
  Real primaryVoltage = 400e3;
  Real secondaryVoltage = 333e3;
  Real nominalDcVoltage = 640e3;

  UInt submodulesPerArm = 36;
  Real submoduleCapacitance = 1.758e-3;

  Real transformerResistancePu = 0.003;
  Real transformerReactancePu = 0.12;

  Real armInductancePu = 0.15;
  Real armResistancePu = 0.0015;

  Real primaryLoadPower = 1000e6 * 20.0 / 30.0;
  Real secondaryLoadPower = 500e3;

  Real rightGridResistance = 0.8929;
  Real rightGridInductance = 16.58e-3;

  Real cableResistancePerPole = 0.5;
  Real cableInductancePerPole = 15e-3;

  Real groundingResistance = 100.0;
  Real groundingCapacitance = 50e-9;

  Real powerReferencePu = 0.99;
  Real reactivePowerReferencePu = 0.0;
  Real controlTimeStep = 40e-6;
  Real dqMeasurementFilterFrequency = 1000.0;

  Real pllKpPu = 180.0;
  Real pllKiPu = 3200.0;
  Real dcVoltageKpPu = 4.0;
  Real dcVoltageKiPu = 100.0;
  Real reactivePowerKpPu = 0.5 / 3.0;
  Real reactivePowerKiPu = 1.0;
  Real currentKpPu = 0.6;
  Real currentKiPu = 6.0;
  Real circulatingKpPu = 1.0;
  Real circulatingKiPu = 5.0;
};

struct ControllerGainsSi {
  Real pllKp;
  Real pllKi;
  Real dcVoltageKp;
  Real dcVoltageKi;
  Real reactivePowerKp;
  Real reactivePowerKi;
  Real currentKp;
  Real currentKi;
  Real circulatingKp;
  Real circulatingKi;
};

ControllerGainsSi convertControllerGainsToSi(const MatlabCaseParameters &p) {
  const Real vBaseAc = std::sqrt(2.0 / 3.0) * p.secondaryVoltage;
  const Real iBaseAc = (2.0 / 3.0) * p.ratedPower / vBaseAc;
  const Real zBaseAc = vBaseAc / iBaseAc;

  ControllerGainsSi gains;
  gains.pllKp = p.pllKpPu / vBaseAc;
  gains.pllKi = p.pllKiPu / vBaseAc;
  gains.dcVoltageKp = p.dcVoltageKpPu * iBaseAc / p.nominalDcVoltage;
  gains.dcVoltageKi = p.dcVoltageKiPu * iBaseAc / p.nominalDcVoltage;
  gains.reactivePowerKp = p.reactivePowerKpPu * iBaseAc / p.ratedPower;
  gains.reactivePowerKi = p.reactivePowerKiPu * iBaseAc / p.ratedPower;
  gains.currentKp = p.currentKpPu * zBaseAc;
  gains.currentKi = p.currentKiPu * zBaseAc;
  gains.circulatingKp = p.circulatingKpPu * zBaseAc;
  gains.circulatingKi = p.circulatingKiPu * zBaseAc;
  return gains;
}

struct AcSideOperatingPoint {
  // External 400-kV bus on which the primary constant-Z load is connected.
  Complex primaryBusPhaseRms;

  // High-voltage transformer terminal after the explicit winding resistor.
  Complex transformerPrimaryPhaseRms;

  // 333-kV converter-side transformer terminal.
  Complex secondaryPhaseRms;
};

struct HotStartPoint {
  AcSideOperatingPoint rectifier;
  AcSideOperatingPoint inverter;
  Complex leftSourcePhaseRms;
  Complex rightSourcePhaseRms;

  Complex rectifierPower;
  Complex inverterPower;

  Real dcCurrent;
  Real rectifierDcVoltage;
  Real inverterDcVoltage;

  Real armInductance;
  Real armResistance;
  Real transformerResistance;
  Real transformerInductance;

  Real rectifierRepresentedLoss;
  Real inverterRepresentedLoss;
};

Vector complexResidual(const std::vector<Complex> &values) {
  Vector result(2 * static_cast<Eigen::Index>(values.size()));
  for (std::size_t index = 0; index < values.size(); ++index) {
    result(2 * static_cast<Eigen::Index>(index), 0) = values[index].real();
    result(2 * static_cast<Eigen::Index>(index) + 1, 0) = values[index].imag();
  }
  return result;
}

template <typename ResidualFunction>
Vector solveNewton(Vector initial, const ResidualFunction &residual,
                   UInt maximumIterations = 50, Real tolerance = 1e-8) {
  Vector x = std::move(initial);

  for (UInt iteration = 0; iteration < maximumIterations; ++iteration) {
    const Vector f = residual(x);
    if (!f.allFinite())
      throw std::runtime_error("Hot-start Newton residual became non-finite.");
    if (f.lpNorm<Eigen::Infinity>() < tolerance)
      return x;

    Matrix jacobian(f.rows(), x.rows());
    for (Eigen::Index column = 0; column < x.rows(); ++column) {
      Vector perturbed = x;
      const Real step = 1e-6 * std::max<Real>(1.0, std::abs(x(column, 0)));
      perturbed(column, 0) += step;
      jacobian.col(column) = (residual(perturbed) - f) / step;
    }

    const Vector correction = jacobian.fullPivLu().solve(-f);
    if (!correction.allFinite())
      throw std::runtime_error(
          "Hot-start Newton correction became non-finite.");

    x += correction;

    if (correction.lpNorm<Eigen::Infinity>() <
        tolerance * std::max<Real>(1.0, x.lpNorm<Eigen::Infinity>()))
      return x;
  }

  throw std::runtime_error(
      "Hot-start AC-network Newton iteration did not converge.");
}

Complex vectorComplex(const Vector &x, Eigen::Index firstIndex) {
  return Complex(x(firstIndex, 0), x(firstIndex + 1, 0));
}

AcSideOperatingPoint
solveLeftAcOperatingPoint(const MatlabCaseParameters &p, Complex converterPower,
                          Complex sourcePhaseRms, Real transformerResistance,
                          Complex transformerLeakageImpedance,
                          Complex transformerHighShuntAdmittance,
                          Complex transformerLowShuntAdmittance) {
  const Real ratio = p.primaryVoltage / p.secondaryVoltage;

  Vector initial(4, 1);
  initial << sourcePhaseRms.real(), sourcePhaseRms.imag(),
      (sourcePhaseRms / ratio).real(), (sourcePhaseRms / ratio).imag();

  const Vector solution = solveNewton(initial, [&](const Vector &x) {
    const Complex transformerPrimary = vectorComplex(x, 0);
    const Complex secondary = vectorComplex(x, 2);

    const Complex resistorCurrent =
        (sourcePhaseRms - transformerPrimary) / transformerResistance;
    const Complex transformerCurrent =
        (transformerPrimary - ratio * secondary) / transformerLeakageImpedance;
    // converterPower is generation-positive, while converterCurrent is
    // positive from the AC node into the converter.
    const Complex converterCurrent =
        -std::conj(converterPower / (3.0 * secondary));

    const Complex transformerPrimaryResidual =
        resistorCurrent - transformerHighShuntAdmittance * transformerPrimary -
        transformerCurrent;
    const Complex secondaryResidual =
        ratio * transformerCurrent - transformerLowShuntAdmittance * secondary -
        converterCurrent;

    return complexResidual({transformerPrimaryResidual, secondaryResidual});
  });

  return {sourcePhaseRms, vectorComplex(solution, 0),
          vectorComplex(solution, 2)};
}

AcSideOperatingPoint solveRightAcOperatingPoint(
    const MatlabCaseParameters &p, Complex converterPower,
    Complex sourcePhaseRms, Complex sourceImpedance, Real transformerResistance,
    Complex transformerLeakageImpedance, Complex primaryLoadAdmittance,
    Complex transformerHighShuntAdmittance,
    Complex transformerLowShuntAdmittance) {
  const Real ratio = p.primaryVoltage / p.secondaryVoltage;

  Vector initial(6, 1);
  initial << sourcePhaseRms.real(), sourcePhaseRms.imag(),
      sourcePhaseRms.real(), sourcePhaseRms.imag(),
      (sourcePhaseRms / ratio).real(), (sourcePhaseRms / ratio).imag();

  const Vector solution = solveNewton(initial, [&](const Vector &x) {
    const Complex primaryBus = vectorComplex(x, 0);
    const Complex transformerPrimary = vectorComplex(x, 2);
    const Complex secondary = vectorComplex(x, 4);

    const Complex sourceCurrent =
        (sourcePhaseRms - primaryBus) / sourceImpedance;
    const Complex resistorCurrent =
        (primaryBus - transformerPrimary) / transformerResistance;
    const Complex transformerCurrent =
        (transformerPrimary - ratio * secondary) / transformerLeakageImpedance;
    // converterPower is generation-positive, while converterCurrent is
    // positive from the AC node into the converter.
    const Complex converterCurrent =
        -std::conj(converterPower / (3.0 * secondary));

    const Complex primaryBusResidual =
        sourceCurrent - primaryLoadAdmittance * primaryBus - resistorCurrent;
    const Complex transformerPrimaryResidual =
        resistorCurrent - transformerHighShuntAdmittance * transformerPrimary -
        transformerCurrent;
    const Complex secondaryResidual =
        ratio * transformerCurrent - transformerLowShuntAdmittance * secondary -
        converterCurrent;

    return complexResidual(
        {primaryBusResidual, transformerPrimaryResidual, secondaryResidual});
  });

  return {vectorComplex(solution, 0), vectorComplex(solution, 2),
          vectorComplex(solution, 4)};
}

Real differentialCurrentLoss(Complex converterPower, Complex secondaryPhaseRms,
                             Real armResistance, Real reactorResistance) {
  const Real peakPhaseVoltage = std::sqrt(2.0) * std::abs(secondaryPhaseRms);
  const Real iDeltaD = (2.0 / 3.0) * converterPower.real() / peakPhaseVoltage;
  const Real iDeltaQ = -(2.0 / 3.0) * converterPower.imag() / peakPhaseVoltage;
  const Real coefficient = 0.75 * armResistance + 1.5 * reactorResistance;
  return coefficient * (iDeltaD * iDeltaD + iDeltaQ * iDeltaQ);
}

HotStartPoint makeHotStartPoint(const MatlabCaseParameters &p) {
  const Real omega = 2.0 * PI * p.frequency;
  const Real acBaseImpedance =
      p.secondaryVoltage * p.secondaryVoltage / p.ratedPower;
  const Real primaryBaseImpedance =
      p.primaryVoltage * p.primaryVoltage / p.ratedPower;

  const Real armInductance = p.armInductancePu * acBaseImpedance / omega;
  const Real armResistance = p.armResistancePu * acBaseImpedance;
  const Real transformerResistance =
      p.transformerResistancePu * primaryBaseImpedance;
  const Real transformerInductance =
      p.transformerReactancePu * primaryBaseImpedance / omega;

  // DPsim's EMT transformer initializes its internal inductor correctly only
  // when no internal series resistor is requested. Therefore the Simulink
  // transformer's total winding resistance is represented by an explicit EMT
  // resistor immediately in front of the transformer. The transformer itself
  // contains the complete leakage inductance and ideal voltage ratio.
  const Complex transformerLeakageImpedance(0.0, omega * transformerInductance);
  const Complex rightSourceImpedance(p.rightGridResistance,
                                     omega * p.rightGridInductance);

  // DPsim's standard EMT transformer adds terminal snubbers. They are not part
  // of the MATLAB network, but they cannot be disabled through the current
  // Transformer API. Include them in the analytical hot start so they do not
  // create a first-step current discontinuity.
  const Real snubberActivePower = P_SNUB_TRANSFORMER * p.ratedPower;
  const Real snubberReactivePower = Q_SNUB_TRANSFORMER * p.ratedPower;

  const Complex primaryLoadAdmittance(
      p.primaryLoadPower / (p.primaryVoltage * p.primaryVoltage), 0.0);
  const Complex transformerHighShuntAdmittance(
      snubberActivePower / (p.primaryVoltage * p.primaryVoltage), 0.0);
  const Complex transformerLowShuntAdmittance(
      (p.secondaryLoadPower + snubberActivePower) /
          (p.secondaryVoltage * p.secondaryVoltage),
      snubberReactivePower / (p.secondaryVoltage * p.secondaryVoltage));

  const Complex leftSourcePhaseRms(p.primaryVoltage / std::sqrt(3.0), 0.0);
  const Complex rightSourcePhaseRms(p.primaryVoltage / std::sqrt(3.0), 0.0);

  // MATLAB power convention: positive P/Q is generated by the converter.
  // Controllers3 forms Id_ref = Pref/Vd, so Pref=+0.99 pu is an inverter
  // operating point, not -0.99 pu.
  const Complex inverterPower(p.powerReferencePu * p.ratedPower,
                              p.reactivePowerReferencePu * p.ratedPower);

  const AcSideOperatingPoint inverter = solveRightAcOperatingPoint(
      p, inverterPower, rightSourcePhaseRms, rightSourceImpedance,
      transformerResistance, transformerLeakageImpedance, primaryLoadAdmittance,
      transformerHighShuntAdmittance, transformerLowShuntAdmittance);

  const Real reactorResistance = 0.0;
  const Real inverterDifferentialLoss =
      differentialCurrentLoss(inverterPower, inverter.secondaryPhaseRms,
                              armResistance, reactorResistance);

  const Real zeroSequenceLossCoefficient = (2.0 / 3.0) * armResistance;
  const Real loopResistance = 2.0 * p.cableResistancePerPole;

  // Positive pole current flows from the Vdc-controlled rectifier to the
  // P-controlled inverter. With generation-positive converter power:
  //
  //   P_inv + P_loss,inv
  //       = Vdc_inv * Idc
  //       = (Vdc_rec - R_loop*Idc)*Idc
  //
  // and P_loss,inv contains kz*Idc^2. Therefore
  //
  //   (R_loop + kz)*Idc^2 - Vdc_rec*Idc
  //       + P_inv + P_loss,diff,inv = 0.
  const Real quadraticA = loopResistance + zeroSequenceLossCoefficient;
  const Real quadraticB = -p.nominalDcVoltage;
  const Real quadraticC = inverterPower.real() + inverterDifferentialLoss;
  const Real discriminant =
      quadraticB * quadraticB - 4.0 * quadraticA * quadraticC;
  require(discriminant > 0.0,
          "No physical DC-current root for the MATLAB hot start.");

  const Real root1 =
      (-quadraticB - std::sqrt(discriminant)) / (2.0 * quadraticA);
  const Real root2 =
      (-quadraticB + std::sqrt(discriminant)) / (2.0 * quadraticA);
  const Real dcCurrent = (root1 > 0.0 && root1 < root2) ? root1 : root2;
  require(dcCurrent > 0.0, "Calculated hot-start DC current is not positive.");

  const Real rectifierDcVoltage = p.nominalDcVoltage;
  const Real inverterDcVoltage =
      rectifierDcVoltage - loopResistance * dcCurrent;
  const Real zeroSequenceLoss =
      zeroSequenceLossCoefficient * dcCurrent * dcCurrent;

  // The sending-end rectifier absorbs AC power and exports DC power:
  //   P_rec + Vdc_rec*Idc + P_loss,rec = 0.
  Real rectifierPower = -rectifierDcVoltage * dcCurrent - zeroSequenceLoss;

  AcSideOperatingPoint rectifier;
  Real rectifierDifferentialLoss = 0.0;
  for (UInt iteration = 0; iteration < 30; ++iteration) {
    rectifier = solveLeftAcOperatingPoint(
        p, Complex(rectifierPower, 0.0), leftSourcePhaseRms,
        transformerResistance, transformerLeakageImpedance,
        transformerHighShuntAdmittance, transformerLowShuntAdmittance);
    rectifierDifferentialLoss = differentialCurrentLoss(
        Complex(rectifierPower, 0.0), rectifier.secondaryPhaseRms,
        armResistance, reactorResistance);
    const Real updatedPower = -rectifierDcVoltage * dcCurrent -
                              zeroSequenceLoss - rectifierDifferentialLoss;
    if (std::abs(updatedPower - rectifierPower) < 1.0) {
      rectifierPower = updatedPower;
      break;
    }
    rectifierPower = updatedPower;
  }

  HotStartPoint point;
  point.rectifier = rectifier;
  point.inverter = inverter;
  point.leftSourcePhaseRms = leftSourcePhaseRms;
  point.rightSourcePhaseRms = rightSourcePhaseRms;
  point.rectifierPower = Complex(rectifierPower, 0.0);
  point.inverterPower = inverterPower;
  point.dcCurrent = dcCurrent;
  point.rectifierDcVoltage = rectifierDcVoltage;
  point.inverterDcVoltage = inverterDcVoltage;
  point.armInductance = armInductance;
  point.armResistance = armResistance;
  point.transformerResistance = transformerResistance;
  point.transformerInductance = transformerInductance;
  point.rectifierRepresentedLoss = zeroSequenceLoss + rectifierDifferentialLoss;
  point.inverterRepresentedLoss = zeroSequenceLoss + inverterDifferentialLoss;
  return point;
}

Real maximumAbsoluteValue(const Matrix &matrix) {
  if (matrix.size() == 0)
    return 0.0;
  return matrix.cwiseAbs().maxCoeff();
}

Real scalarAttribute(const EMT::Ph3::SSN_MMC::Ptr &mmc, const String &name) {
  return mmc->attributeTyped<Real>(name)->get();
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
  Real equilibriumResidual = 0.0;
  Real pllError = 0.0;
  Real pllAngleDeviation = 0.0;
  Real pllFrequency = 0.0;
  Real gridVoltageD = 0.0;
  Real gridVoltageQ = 0.0;
  Real deltaCurrentD = 0.0;
  Real deltaCurrentQ = 0.0;
  Real sigmaCurrentZ = 0.0;
  Real deltaCurrentReferenceD = 0.0;
  Real deltaCurrentReferenceQ = 0.0;
  Real sigmaCurrentReferenceZ = 0.0;
  Real filteredDaxisVoltage = 0.0;
  Real heldActiveCurrentReference = 0.0;

  String largestStateName = "none";
  String largestDerivativeName = "none";
  Real largestStateMagnitude = 0.0;
  Real largestDerivativeMagnitude = 0.0;
};

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
  snapshot.equilibriumResidual =
      scalarAttribute(mmc, "equilibrium_residual_norm");
  snapshot.pllError = scalarAttribute(mmc, "pll_error");
  snapshot.pllAngleDeviation = scalarAttribute(mmc, "pll_angle_deviation");
  snapshot.pllFrequency = scalarAttribute(mmc, "pll_frequency");
  snapshot.gridVoltageD = scalarAttribute(mmc, "v_grid_d");
  snapshot.gridVoltageQ = scalarAttribute(mmc, "v_grid_q");
  snapshot.deltaCurrentD = scalarAttribute(mmc, "i_delta_d");
  snapshot.deltaCurrentQ = scalarAttribute(mmc, "i_delta_q");
  snapshot.sigmaCurrentZ = scalarAttribute(mmc, "i_sigma_z");
  snapshot.deltaCurrentReferenceD = scalarAttribute(mmc, "i_delta_d_ref");
  snapshot.deltaCurrentReferenceQ = scalarAttribute(mmc, "i_delta_q_ref");
  snapshot.sigmaCurrentReferenceZ = scalarAttribute(mmc, "i_sigma_z_ref");
  snapshot.filteredDaxisVoltage =
      scalarAttribute(mmc, "v_d_feedforward_filtered");
  snapshot.heldActiveCurrentReference =
      scalarAttribute(mmc, "i_delta_d_feedforward_held");

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
  MmcDiagnosticSnapshot rectifier;
  MmcDiagnosticSnapshot inverter;

  Real positiveLocalCurrent = 0.0;
  Real positiveRemoteCurrent = 0.0;
  Real negativeLocalCurrent = 0.0;
  Real negativeRemoteCurrent = 0.0;
  Real positiveGroundResistorCurrent = 0.0;
  Real negativeGroundResistorCurrent = 0.0;
  Real positiveGroundCapacitorCurrent = 0.0;
  Real negativeGroundCapacitorCurrent = 0.0;

  Real rectifierPositivePole = 0.0;
  Real rectifierNegativePole = 0.0;
  Real inverterPositivePole = 0.0;
  Real inverterNegativePole = 0.0;
  Real rectifierCommonMode = 0.0;
  Real inverterCommonMode = 0.0;

  Real cableLoss = 0.0;
  Real groundingLoss = 0.0;
  Real dcPowerResidual = 0.0;
  Real maximumKclResidual = 0.0;
};

P2PDiagnosticSample
captureP2PSample(Real time, const EMT::Ph3::SSN_MMC::Ptr &rectifier,
                 const EMT::Ph3::SSN_MMC::Ptr &inverter,
                 const EMT::DC::SSN::PiLine::Ptr &positiveLineLocal,
                 const EMT::DC::SSN::PiLine::Ptr &positiveLineRemote,
                 const EMT::DC::SSN::PiLine::Ptr &negativeLineLocal,
                 const EMT::DC::SSN::PiLine::Ptr &negativeLineRemote,
                 const EMT::DC::SSN::Resistor::Ptr &positiveGroundResistor,
                 const EMT::DC::SSN::Resistor::Ptr &negativeGroundResistor,
                 const EMT::DC::SSN::Capacitor::Ptr &positiveGroundCapacitor,
                 const EMT::DC::SSN::Capacitor::Ptr &negativeGroundCapacitor,
                 const SimNode<Real>::Ptr &rectifierPositive,
                 const SimNode<Real>::Ptr &rectifierNegative,
                 const SimNode<Real>::Ptr &positiveCableMid,
                 const SimNode<Real>::Ptr &negativeCableMid,
                 const SimNode<Real>::Ptr &inverterPositive,
                 const SimNode<Real>::Ptr &inverterNegative,
                 Real cableSectionResistance, Real groundingResistance,
                 Bool evaluateDerivative) {
  P2PDiagnosticSample sample;
  sample.time = time;
  sample.rectifier = captureMmcSnapshot(rectifier, evaluateDerivative);
  sample.inverter = captureMmcSnapshot(inverter, evaluateDerivative);

  sample.positiveLocalCurrent = positiveLineLocal->intfCurrent()(0, 0);
  sample.positiveRemoteCurrent = positiveLineRemote->intfCurrent()(0, 0);
  sample.negativeLocalCurrent = negativeLineLocal->intfCurrent()(0, 0);
  sample.negativeRemoteCurrent = negativeLineRemote->intfCurrent()(0, 0);
  sample.positiveGroundResistorCurrent =
      positiveGroundResistor->intfCurrent()(0, 0);
  sample.negativeGroundResistorCurrent =
      negativeGroundResistor->intfCurrent()(0, 0);
  sample.positiveGroundCapacitorCurrent =
      positiveGroundCapacitor->intfCurrent()(0, 0);
  sample.negativeGroundCapacitorCurrent =
      negativeGroundCapacitor->intfCurrent()(0, 0);

  sample.rectifierPositivePole = rectifierPositive->voltage()(0, 0);
  sample.rectifierNegativePole = rectifierNegative->voltage()(0, 0);
  sample.inverterPositivePole = inverterPositive->voltage()(0, 0);
  sample.inverterNegativePole = inverterNegative->voltage()(0, 0);
  sample.rectifierCommonMode =
      0.5 * (sample.rectifierPositivePole + sample.rectifierNegativePole);
  sample.inverterCommonMode =
      0.5 * (sample.inverterPositivePole + sample.inverterNegativePole);

  sample.cableLoss =
      cableSectionResistance *
      (sample.positiveLocalCurrent * sample.positiveLocalCurrent +
       sample.positiveRemoteCurrent * sample.positiveRemoteCurrent +
       sample.negativeLocalCurrent * sample.negativeLocalCurrent +
       sample.negativeRemoteCurrent * sample.negativeRemoteCurrent);
  sample.groundingLoss =
      groundingResistance * (sample.positiveGroundResistorCurrent *
                                 sample.positiveGroundResistorCurrent +
                             sample.negativeGroundResistorCurrent *
                                 sample.negativeGroundResistorCurrent);
  // Converter p_dc is positive when DC power is absorbed by a converter.
  // Cable and grounding losses are also absorbed powers, so their sum must
  // vanish for the complete isolated DC network.
  sample.dcPowerResidual = sample.rectifier.dcPower + sample.inverter.dcPower +
                           sample.cableLoss + sample.groundingLoss;

  const Real localPositiveKcl = sample.positiveLocalCurrent +
                                sample.positiveGroundResistorCurrent +
                                sample.rectifier.interfaceCurrent(3, 0);
  const Real localNegativeKcl = -sample.negativeLocalCurrent -
                                sample.negativeGroundResistorCurrent +
                                sample.rectifier.interfaceCurrent(4, 0);
  const Real positiveMidKcl =
      -sample.positiveLocalCurrent + sample.positiveRemoteCurrent;
  const Real negativeMidKcl =
      sample.negativeLocalCurrent - sample.negativeRemoteCurrent;
  const Real remotePositiveKcl =
      -sample.positiveRemoteCurrent + sample.inverter.interfaceCurrent(3, 0);
  const Real remoteNegativeKcl =
      sample.negativeRemoteCurrent + sample.inverter.interfaceCurrent(4, 0);

  sample.maximumKclResidual =
      std::max({std::abs(localPositiveKcl), std::abs(localNegativeKcl),
                std::abs(positiveMidKcl), std::abs(negativeMidKcl),
                std::abs(remotePositiveKcl), std::abs(remoteNegativeKcl)});

  (void)positiveCableMid;
  (void)negativeCableMid;
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

  (void)mmc->getStateDerivative();

  SPDLOG_INFO("{} initial local Jacobian: max_real={} 1/s, "
              "max_magnitude={} 1/s, dominant_frequency={} Hz",
              name, maximumRealPart, maximumMagnitude, dominantFrequency);
}

void printCompactDiagnostic(const P2PDiagnosticSample &sample,
                            Real initialRectifierEnergy,
                            Real initialInverterEnergy) {
  SPDLOG_INFO(
      "MATLAB P2P diagnostic t={} s: "
      "REC[Vdc={},Idc={},P={},Q={},dE={},pll_err={},Vd_f={},Id_ref_h={},"
      "max_x={}({}),max_dx={}({}),|m|max={}] "
      "INV[Vdc={},Idc={},P={},Q={},dE={},pll_err={},Vd_f={},Id_ref_h={},"
      "max_x={}({}),max_dx={}({}),|m|max={}] "
      "DC[I+=[{},{}],I-=[{},{}],Vcm=[{},{}],"
      "loss={},power_residual={},max_KCL={}]",
      sample.time, sample.rectifier.vdc, sample.rectifier.idc,
      sample.rectifier.activePower, sample.rectifier.reactivePower,
      sample.rectifier.storedEnergy - initialRectifierEnergy,
      sample.rectifier.pllError, sample.rectifier.filteredDaxisVoltage,
      sample.rectifier.heldActiveCurrentReference,
      sample.rectifier.largestStateMagnitude, sample.rectifier.largestStateName,
      sample.rectifier.largestDerivativeMagnitude,
      sample.rectifier.largestDerivativeName,
      maximumAbsoluteValue(sample.rectifier.modulation), sample.inverter.vdc,
      sample.inverter.idc, sample.inverter.activePower,
      sample.inverter.reactivePower,
      sample.inverter.storedEnergy - initialInverterEnergy,
      sample.inverter.pllError, sample.inverter.filteredDaxisVoltage,
      sample.inverter.heldActiveCurrentReference,
      sample.inverter.largestStateMagnitude, sample.inverter.largestStateName,
      sample.inverter.largestDerivativeMagnitude,
      sample.inverter.largestDerivativeName,
      maximumAbsoluteValue(sample.inverter.modulation),
      sample.positiveLocalCurrent, sample.positiveRemoteCurrent,
      sample.negativeLocalCurrent, sample.negativeRemoteCurrent,
      sample.rectifierCommonMode, sample.inverterCommonMode,
      sample.cableLoss + sample.groundingLoss, sample.dcPowerResidual,
      sample.maximumKclResidual);
}

void writeMmcColumns(std::ofstream &stream, const String &prefix,
                     const std::vector<String> &stateNames) {
  stream << ',' << prefix << "_vdc" << ',' << prefix << "_idc" << ',' << prefix
         << "_p" << ',' << prefix << "_q" << ',' << prefix << "_energy" << ','
         << prefix << "_pdc" << ',' << prefix << "_power_balance_error" << ','
         << prefix << "_pll_error" << ',' << prefix << "_pll_angle_deviation"
         << ',' << prefix << "_pll_frequency" << ',' << prefix << "_v_grid_d"
         << ',' << prefix << "_v_grid_q" << ',' << prefix << "_i_delta_d" << ','
         << prefix << "_i_delta_q" << ',' << prefix << "_i_sigma_z" << ','
         << prefix << "_i_delta_d_ref" << ',' << prefix << "_i_delta_q_ref"
         << ',' << prefix << "_i_sigma_z_ref" << ',' << prefix
         << "_v_d_feedforward_filtered" << ',' << prefix
         << "_i_delta_d_feedforward_held" << ',' << prefix << "_modulation_max"
         << ',' << prefix << "_differential_voltage_max" << ',' << prefix
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
         << snapshot.pllAngleDeviation << ',' << snapshot.pllFrequency << ','
         << snapshot.gridVoltageD << ',' << snapshot.gridVoltageQ << ','
         << snapshot.deltaCurrentD << ',' << snapshot.deltaCurrentQ << ','
         << snapshot.sigmaCurrentZ << ',' << snapshot.deltaCurrentReferenceD
         << ',' << snapshot.deltaCurrentReferenceQ << ','
         << snapshot.sigmaCurrentReferenceZ << ','
         << snapshot.filteredDaxisVoltage << ','
         << snapshot.heldActiveCurrentReference << ','
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
                       const std::vector<String> &rectifierStateNames,
                       const std::vector<String> &inverterStateNames) {
  std::ofstream stream(fileName);
  if (!stream)
    throw std::runtime_error("Could not open diagnostic trace file: " +
                             fileName);

  stream << std::setprecision(17);
  stream << "time,"
            "positive_local_current,"
            "positive_remote_current,"
            "negative_local_current,"
            "negative_remote_current,"
            "positive_ground_resistor_current,"
            "negative_ground_resistor_current,"
            "positive_ground_capacitor_current,"
            "negative_ground_capacitor_current,"
            "rectifier_positive_pole,"
            "rectifier_negative_pole,"
            "inverter_positive_pole,"
            "inverter_negative_pole,"
            "rectifier_common_mode,"
            "inverter_common_mode,"
            "cable_loss,"
            "grounding_loss,"
            "dc_power_residual,"
            "maximum_kcl_residual";
  writeMmcColumns(stream, "rectifier", rectifierStateNames);
  writeMmcColumns(stream, "inverter", inverterStateNames);
  stream << '\n';

  for (const auto &sample : trace) {
    stream << sample.time << ',' << sample.positiveLocalCurrent << ','
           << sample.positiveRemoteCurrent << ',' << sample.negativeLocalCurrent
           << ',' << sample.negativeRemoteCurrent << ','
           << sample.positiveGroundResistorCurrent << ','
           << sample.negativeGroundResistorCurrent << ','
           << sample.positiveGroundCapacitorCurrent << ','
           << sample.negativeGroundCapacitorCurrent << ','
           << sample.rectifierPositivePole << ','
           << sample.rectifierNegativePole << ',' << sample.inverterPositivePole
           << ',' << sample.inverterNegativePole << ','
           << sample.rectifierCommonMode << ',' << sample.inverterCommonMode
           << ',' << sample.cableLoss << ',' << sample.groundingLoss << ','
           << sample.dcPowerResidual << ',' << sample.maximumKclResidual;
    writeMmcValues(stream, sample.rectifier);
    writeMmcValues(stream, sample.inverter);
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

void configureMmcCommon(const EMT::Ph3::SSN_MMC::Ptr &mmc,
                        const MatlabCaseParameters &p,
                        const ControllerGainsSi &gains, const HotStartPoint &op,
                        Real initialActivePower, Real initialReactivePower,
                        Real initialDcCurrent, Real initialAngle, Real theta,
                        Real timeStep) {
  const Real vBaseAc = std::sqrt(2.0 / 3.0) * p.secondaryVoltage;
  const Real iBaseAc = (2.0 / 3.0) * p.ratedPower / vBaseAc;
  const Real iBaseDc = p.ratedPower / p.nominalDcVoltage;

  mmc->setParameters(p.frequency, p.secondaryVoltage, p.nominalDcVoltage,
                     op.armInductance, op.armResistance, p.submoduleCapacitance,
                     p.submodulesPerArm, 0.0, 0.0);
  mmc->setInitialAngle(initialAngle);
  mmc->setInitialOperatingPoint(initialActivePower, initialReactivePower);

  mmc->setPLL(gains.pllKp, gains.pllKi, true);
  mmc->setOutputCurrentController(gains.currentKp, gains.currentKi);
  mmc->setCirculatingCurrentController(gains.circulatingKp,
                                       gains.circulatingKi);

  // The supplied Simulink controller has circulating-current suppression,
  // but no separate total-energy or zero-sequence-current PI.
  mmc->setEnergyController(0.0, 0.0, false);
  mmc->setZeroSequenceCurrentController(0.0, 0.0);

  // Positive iSigma_z is converter DC absorption. Therefore the sending
  // rectifier uses a negative local DC-current reference and the receiving
  // inverter uses a positive one.
  mmc->setCirculatingCurrentReferences(0.0, 0.0, initialDcCurrent / 3.0);
  mmc->setLimits(2.0 * iBaseAc, 2.0 * iBaseDc, 2.0);
  mmc->setOperatingPointInitialization(true, 100, 1e-9);
  mmc->setDiagnosticTimeStep(timeStep);
  mmc->setTheta(theta);
}

} // namespace

int main(int argc, char **argv) {
  const Real timeStep = argc > 1 ? std::stod(argv[1]) : 20e-6;
  const Real finalTime = argc > 2 ? std::stod(argv[2]) : 0.5;
  const Real theta = argc > 3 ? std::stod(argv[3]) : 0.5;
  const UInt reportIntervalSteps =
      argc > 4 ? static_cast<UInt>(std::stoul(argv[4])) : 500;
  const ClosedLoopDiagnosticMode diagnosticMode =
      argc > 5 ? parseClosedLoopDiagnosticMode(argv[5])
               : ClosedLoopDiagnosticMode::InternalControllers;

  require(timeStep > 0.0, "The EMT time step must be positive.");
  require(finalTime > 0.0, "The final time must be positive.");
  require(theta >= 0.5 && theta <= 1.0, "The SSN theta must be in [0.5,1.0].");
  require(reportIntervalSteps > 0, "The report interval must be positive.");

  const String simName = "EMT_SSN_MMC_Matlab_P2P_HotStart";
  Logger::setLogDir("logs/" + simName);

  const MatlabCaseParameters p;
  const ControllerGainsSi gains = convertControllerGainsToSi(p);
  const HotStartPoint op = makeHotStartPoint(p);

  const Real omega = 2.0 * PI * p.frequency;
  const Real cableSectionResistance = p.cableResistancePerPole / 2.0;
  const Real cableSectionInductance = p.cableInductancePerPole / 2.0;

  SPDLOG_INFO(
      "MATLAB P2P hot start: mode={}, "
      "Pnom={} VA, Vac=[{},{}] V, Vdc_nom={} V, "
      "MMC P/Q=[{},{},{},{}] W/var, "
      "MMC Vac_LL=[{},{}] V, "
      "TR_HV_LL=[{},{}] V, "
      "Vdc=[{},{}] V, Idc={} A, "
      "MMC losses=[{},{}] W, "
      "Larm={} H, Rarm={} ohm, "
      "transformer R/L=[{},{}], "
      "right_grid_I_rms={} A angle={} deg",
      closedLoopDiagnosticModeName(diagnosticMode), p.ratedPower,
      p.primaryVoltage, p.secondaryVoltage, p.nominalDcVoltage,
      op.rectifierPower.real(), op.rectifierPower.imag(),
      op.inverterPower.real(), op.inverterPower.imag(),
      std::abs(phaseRmsToLineLineRms(op.rectifier.secondaryPhaseRms)),
      std::abs(phaseRmsToLineLineRms(op.inverter.secondaryPhaseRms)),
      std::abs(phaseRmsToLineLineRms(op.rectifier.transformerPrimaryPhaseRms)),
      std::abs(phaseRmsToLineLineRms(op.inverter.transformerPrimaryPhaseRms)),
      op.rectifierDcVoltage, op.inverterDcVoltage, op.dcCurrent,
      op.rectifierRepresentedLoss, op.inverterRepresentedLoss, op.armInductance,
      op.armResistance, op.transformerResistance, op.transformerInductance,
      std::abs((op.rightSourcePhaseRms - op.inverter.primaryBusPhaseRms) /
               Complex(p.rightGridResistance, omega * p.rightGridInductance)),
      std::arg((op.rightSourcePhaseRms - op.inverter.primaryBusPhaseRms) /
               Complex(p.rightGridResistance, omega * p.rightGridInductance)) *
          180.0 / PI);

  SPDLOG_INFO("Converted MATLAB controller gains: "
              "PLL=[{},{}], Vdc=[{},{}], Q=[{},{}], "
              "I=[{},{}], CCC=[{},{}]",
              gains.pllKp, gains.pllKi, gains.dcVoltageKp, gains.dcVoltageKi,
              gains.reactivePowerKp, gains.reactivePowerKi, gains.currentKp,
              gains.currentKi, gains.circulatingKp, gains.circulatingKi);

  // --------------------------------------------------------------------------
  // AC topology of P2PHVDCMMC.slx in its post-startup conducting state:
  //
  // ideal 400-kV source
  //      |---- 666.667-MW primary constant-Z load
  // explicit transformer winding R + 400/333-kV transformer leakage L/ratio
  //      |---- 0.5-MW secondary constant-Z load
  // rectifier MMC
  //
  // inverter MMC
  //      |---- 0.5-MW secondary constant-Z load
  // 333/400-kV transformer leakage L/ratio + explicit winding R
  //      |---- 666.667-MW primary constant-Z load
  // 0.8929-ohm resistor + 16.58-mH inductor
  // ideal 400-kV source
  //
  // Startup resistors and breakers are absent because the case starts after
  // their bypass state has already been reached.
  // --------------------------------------------------------------------------
  auto leftPrimary = SimNode<Real>::make("LEFT_400KV", PhaseType::ABC);
  auto leftTransformerPrimary =
      SimNode<Real>::make("LEFT_TRANSFORMER_400KV", PhaseType::ABC);
  auto leftSecondary = SimNode<Real>::make("LEFT_333KV", PhaseType::ABC);
  auto rightSecondary = SimNode<Real>::make("RIGHT_333KV", PhaseType::ABC);
  auto rightTransformerPrimary =
      SimNode<Real>::make("RIGHT_TRANSFORMER_400KV", PhaseType::ABC);
  auto rightPrimary = SimNode<Real>::make("RIGHT_400KV", PhaseType::ABC);
  auto rightGridAfterResistance =
      SimNode<Real>::make("RIGHT_GRID_AFTER_R", PhaseType::ABC);
  auto rightSourceNode = SimNode<Real>::make("RIGHT_SOURCE", PhaseType::ABC);

  leftPrimary->setInitialVoltage(
      phaseRmsToLineLineRms(op.rectifier.primaryBusPhaseRms));
  leftTransformerPrimary->setInitialVoltage(
      phaseRmsToLineLineRms(op.rectifier.transformerPrimaryPhaseRms));
  leftSecondary->setInitialVoltage(
      phaseRmsToLineLineRms(op.rectifier.secondaryPhaseRms));
  rightSecondary->setInitialVoltage(
      phaseRmsToLineLineRms(op.inverter.secondaryPhaseRms));
  rightTransformerPrimary->setInitialVoltage(
      phaseRmsToLineLineRms(op.inverter.transformerPrimaryPhaseRms));
  rightPrimary->setInitialVoltage(
      phaseRmsToLineLineRms(op.inverter.primaryBusPhaseRms));

  const Complex rightGridImpedance(p.rightGridResistance,
                                   omega * p.rightGridInductance);
  const Complex rightGridCurrentPhaseRms =
      (op.rightSourcePhaseRms - op.inverter.primaryBusPhaseRms) /
      rightGridImpedance;
  const Complex rightGridAfterResistancePhaseRms =
      op.rightSourcePhaseRms - p.rightGridResistance * rightGridCurrentPhaseRms;

  rightGridAfterResistance->setInitialVoltage(
      phaseRmsToLineLineRms(rightGridAfterResistancePhaseRms));
  rightSourceNode->setInitialVoltage(
      phaseRmsToLineLineRms(op.rightSourcePhaseRms));

  auto leftSource = EMT::Ph3::NetworkInjection::make("LEFT_SOURCE");
  leftSource->setParameters(balancedLineLineRmsReference(
                                phaseRmsToLineLineRms(op.leftSourcePhaseRms)),
                            p.frequency);
  leftSource->connect({leftPrimary});

  auto rightSource = EMT::Ph3::NetworkInjection::make("RIGHT_SOURCE_IDEAL");
  rightSource->setParameters(balancedLineLineRmsReference(
                                 phaseRmsToLineLineRms(op.rightSourcePhaseRms)),
                             p.frequency);
  rightSource->connect({rightSourceNode});

  // Represent the Simulink right-grid series R-L branch using separate
  // primitive EMT components. EMT::Ph3::RxLine is a composite R-L element
  // whose internal inductor does not receive the loaded branch current at
  // initialization. Splitting the branch is electrically identical but
  // permits an exact state hot start.
  auto rightGridResistance = EMT::Ph3::Resistor::make("RIGHT_GRID_R");
  auto rightGridInductance = EMT::Ph3::Inductor::make("RIGHT_GRID_L");

  rightGridResistance->setParameters(diagonalThreePhase(p.rightGridResistance));
  rightGridInductance->setParameters(diagonalThreePhase(p.rightGridInductance));

  // Both component current attributes use terminal-1 minus terminal-0
  // orientation. These connection orders therefore give the same current
  // direction in the series R and L components.
  rightGridResistance->connect({rightSourceNode, rightGridAfterResistance});
  rightGridInductance->connect({rightGridAfterResistance, rightPrimary});

  // The total Simulink transformer winding resistance is placed in an
  // explicit algebraic EMT resistor. This is electrically equivalent in the
  // balanced positive-sequence case, while allowing every dynamic inductor to
  // receive an exact steady-state current at initialization.
  auto leftTransformerResistance = EMT::Ph3::Resistor::make("TR_RECTIFIER_R");
  auto rightTransformerResistance = EMT::Ph3::Resistor::make("TR_INVERTER_R");
  leftTransformerResistance->setParameters(
      diagonalThreePhase(op.transformerResistance));
  rightTransformerResistance->setParameters(
      diagonalThreePhase(op.transformerResistance));
  leftTransformerResistance->connect({leftTransformerPrimary, leftPrimary});
  rightTransformerResistance->connect({rightTransformerPrimary, rightPrimary});

  auto leftTransformer = std::make_shared<EMT::Ph3::Transformer>(
      "TR_RECTIFIER", "TR_RECTIFIER", Logger::Level::off, false);
  auto rightTransformer = std::make_shared<EMT::Ph3::Transformer>(
      "TR_INVERTER", "TR_INVERTER", Logger::Level::off, false);

  const Real transformerRatio = p.primaryVoltage / p.secondaryVoltage;
  const Matrix transformerInductance =
      diagonalThreePhase(op.transformerInductance);

  leftTransformer->setParameters(p.primaryVoltage, p.secondaryVoltage,
                                 p.ratedPower, transformerRatio, 0.0,
                                 Matrix::Zero(3, 3), transformerInductance);
  rightTransformer->setParameters(p.primaryVoltage, p.secondaryVoltage,
                                  p.ratedPower, transformerRatio, 0.0,
                                  Matrix::Zero(3, 3), transformerInductance);
  leftTransformer->connect({leftTransformerPrimary, leftSecondary});
  rightTransformer->connect({rightTransformerPrimary, rightSecondary});

  auto leftPrimaryLoad = EMT::Ph3::RXLoad::make("LEFT_PRIMARY_LOAD");
  auto rightPrimaryLoad = EMT::Ph3::RXLoad::make("RIGHT_PRIMARY_LOAD");
  auto leftSecondaryLoad = EMT::Ph3::RXLoad::make("LEFT_SECONDARY_LOAD");
  auto rightSecondaryLoad = EMT::Ph3::RXLoad::make("RIGHT_SECONDARY_LOAD");

  leftPrimaryLoad->setParameters(
      Math::singlePhasePowerToThreePhase(p.primaryLoadPower),
      Matrix::Zero(3, 3), p.primaryVoltage);
  rightPrimaryLoad->setParameters(
      Math::singlePhasePowerToThreePhase(p.primaryLoadPower),
      Matrix::Zero(3, 3), p.primaryVoltage);
  leftSecondaryLoad->setParameters(
      Math::singlePhasePowerToThreePhase(p.secondaryLoadPower),
      Matrix::Zero(3, 3), p.secondaryVoltage);
  rightSecondaryLoad->setParameters(
      Math::singlePhasePowerToThreePhase(p.secondaryLoadPower),
      Matrix::Zero(3, 3), p.secondaryVoltage);

  leftPrimaryLoad->connect({leftPrimary});
  rightPrimaryLoad->connect({rightPrimary});
  leftSecondaryLoad->connect({leftSecondary});
  rightSecondaryLoad->connect({rightSecondary});

  // --------------------------------------------------------------------------
  // DC topology of the Simulink subsystem:
  //
  // rectifier + -- RL/2 -- midpoint -- RL/2 -- inverter +
  // rectifier - -- RL/2 -- midpoint -- RL/2 -- inverter -
  //
  // rectifier + -- 100 ohm -- C=50 nF -- GND
  // rectifier - -- 100 ohm -- C=50 nF -- GND
  //
  // All cable, resistor and capacitor elements use native EMT::DC components.
  // --------------------------------------------------------------------------
  auto rectifierPositive =
      SimNode<Real>::make("RECTIFIER_DC_POS", PhaseType::DC);
  auto rectifierNegative =
      SimNode<Real>::make("RECTIFIER_DC_NEG", PhaseType::DC);
  auto positiveCableMid = SimNode<Real>::make("DC_POS_MID", PhaseType::DC);
  auto negativeCableMid = SimNode<Real>::make("DC_NEG_MID", PhaseType::DC);
  auto inverterPositive = SimNode<Real>::make("INVERTER_DC_POS", PhaseType::DC);
  auto inverterNegative = SimNode<Real>::make("INVERTER_DC_NEG", PhaseType::DC);
  auto positiveGroundRcNode =
      SimNode<Real>::make("DC_POS_GROUND_RC", PhaseType::DC);
  auto negativeGroundRcNode =
      SimNode<Real>::make("DC_NEG_GROUND_RC", PhaseType::DC);

  const Real rectifierPositiveVoltage = op.rectifierDcVoltage / 2.0;
  const Real rectifierNegativeVoltage = -op.rectifierDcVoltage / 2.0;
  const Real inverterPositiveVoltage = op.inverterDcVoltage / 2.0;
  const Real inverterNegativeVoltage = -op.inverterDcVoltage / 2.0;
  const Real positiveMidVoltage =
      0.5 * (rectifierPositiveVoltage + inverterPositiveVoltage);
  const Real negativeMidVoltage =
      0.5 * (rectifierNegativeVoltage + inverterNegativeVoltage);

  rectifierPositive->setInitialVoltage(Complex(rectifierPositiveVoltage, 0.0));
  rectifierNegative->setInitialVoltage(Complex(rectifierNegativeVoltage, 0.0));
  positiveCableMid->setInitialVoltage(Complex(positiveMidVoltage, 0.0));
  negativeCableMid->setInitialVoltage(Complex(negativeMidVoltage, 0.0));
  inverterPositive->setInitialVoltage(Complex(inverterPositiveVoltage, 0.0));
  inverterNegative->setInitialVoltage(Complex(inverterNegativeVoltage, 0.0));
  positiveGroundRcNode->setInitialVoltage(
      Complex(rectifierPositiveVoltage, 0.0));
  negativeGroundRcNode->setInitialVoltage(
      Complex(rectifierNegativeVoltage, 0.0));

  auto positiveLineLocal = EMT::DC::SSN::PiLine::make("DC_POS_LOCAL");
  auto positiveLineRemote = EMT::DC::SSN::PiLine::make("DC_POS_REMOTE");
  auto negativeLineLocal = EMT::DC::SSN::PiLine::make("DC_NEG_LOCAL");
  auto negativeLineRemote = EMT::DC::SSN::PiLine::make("DC_NEG_REMOTE");

  // Hot-start the four series inductors with the physical pole current.
  //
  // PiLine's fifth parameter is the initial series current. Leaving this at
  // zero would create a full-load pole-current discontinuity at the first
  // MNA step, despite the MMCs being initialized at the loaded equilibrium.
  positiveLineLocal->setParameters(
      cableSectionResistance, cableSectionInductance, 0.0, 0.0, op.dcCurrent);
  positiveLineRemote->setParameters(
      cableSectionResistance, cableSectionInductance, 0.0, 0.0, op.dcCurrent);
  negativeLineLocal->setParameters(
      cableSectionResistance, cableSectionInductance, 0.0, 0.0, op.dcCurrent);
  negativeLineRemote->setParameters(
      cableSectionResistance, cableSectionInductance, 0.0, 0.0, op.dcCurrent);

  for (const auto &line : {positiveLineLocal, positiveLineRemote,
                           negativeLineLocal, negativeLineRemote})
    line->setTheta(theta);

  // Same current orientations as the validated DC-cable diagnostic:
  // positive pole: rectifier -> inverter
  // negative pole: inverter -> rectifier
  positiveLineLocal->connect({positiveCableMid, rectifierPositive});
  positiveLineRemote->connect({inverterPositive, positiveCableMid});
  negativeLineLocal->connect({rectifierNegative, negativeCableMid});
  negativeLineRemote->connect({negativeCableMid, inverterNegative});

  auto positiveGroundResistor = EMT::DC::SSN::Resistor::make("DC_POS_GROUND_R");
  auto negativeGroundResistor = EMT::DC::SSN::Resistor::make("DC_NEG_GROUND_R");
  positiveGroundResistor->setParameters(p.groundingResistance);
  negativeGroundResistor->setParameters(p.groundingResistance);
  positiveGroundResistor->connect({positiveGroundRcNode, rectifierPositive});
  negativeGroundResistor->connect({rectifierNegative, negativeGroundRcNode});

  auto positiveGroundCapacitor =
      EMT::DC::SSN::Capacitor::make("DC_POS_GROUND_C");
  auto negativeGroundCapacitor =
      EMT::DC::SSN::Capacitor::make("DC_NEG_GROUND_C");
  positiveGroundCapacitor->setParameters(p.groundingCapacitance);
  negativeGroundCapacitor->setParameters(p.groundingCapacitance);
  positiveGroundCapacitor->setTheta(theta);
  negativeGroundCapacitor->setTheta(theta);
  positiveGroundCapacitor->connect({SimNode<Real>::GND, positiveGroundRcNode});
  negativeGroundCapacitor->connect({negativeGroundRcNode, SimNode<Real>::GND});

  auto rectifier = EMT::Ph3::SSN_MMC::make("MMC_RECTIFIER");
  auto inverter = EMT::Ph3::SSN_MMC::make("MMC_INVERTER");

  configureMmcCommon(rectifier, p, gains, op, op.rectifierPower.real(),
                     op.rectifierPower.imag(), -op.dcCurrent,
                     std::arg(op.rectifier.secondaryPhaseRms), theta, timeStep);
  configureMmcCommon(inverter, p, gains, op, op.inverterPower.real(),
                     op.inverterPower.imag(), op.dcCurrent,
                     std::arg(op.inverter.secondaryPhaseRms), theta, timeStep);

  // Sending-end rectifier: Vdc/Q control. MATLAB positive AC power is
  // generation, so this operating point has negative P and negative iDelta_d.
  rectifier->setDcVoltageControl(p.nominalDcVoltage, gains.dcVoltageKp,
                                 gains.dcVoltageKi);
  rectifier->setReactivePowerControl(0.0, gains.reactivePowerKp,
                                     gains.reactivePowerKi);

  // Receiving-end inverter: permanent MATLAB Controllers3 path.
  //
  //   Id_ref = (2/3) * Pref / Vd_filtered
  //
  // The critically damped 1-kHz second-order voltage filter, 40-us sampled
  // execution and zero-order-held current reference are implemented inside
  // SSN_MMC. The example only changes Pref at the MATLAB event times.
  const Real initialInverterIdReference =
      (2.0 / 3.0) * op.inverterPower.real() /
      (std::sqrt(2.0) * std::abs(op.inverter.secondaryPhaseRms));
  inverter->setActivePowerFeedforwardControl(op.inverterPower.real(),
                                             p.dqMeasurementFilterFrequency,
                                             p.controlTimeStep);
  inverter->setReactivePowerControl(0.0, gains.reactivePowerKp,
                                    gains.reactivePowerKi);

  rectifier->connect({leftSecondary, rectifierPositive, rectifierNegative});
  inverter->connect({rightSecondary, inverterPositive, inverterNegative});

  SystemTopology system(
      p.frequency,
      SystemNodeList{leftPrimary, leftTransformerPrimary, leftSecondary,
                     rightSecondary, rightTransformerPrimary, rightPrimary,
                     rightGridAfterResistance, rightSourceNode,
                     rectifierPositive, rectifierNegative, positiveCableMid,
                     negativeCableMid, inverterPositive, inverterNegative,
                     positiveGroundRcNode, negativeGroundRcNode},
      SystemComponentList{leftSource,
                          rightSource,
                          rightGridResistance,
                          rightGridInductance,
                          leftTransformerResistance,
                          rightTransformerResistance,
                          leftTransformer,
                          rightTransformer,
                          leftPrimaryLoad,
                          rightPrimaryLoad,
                          leftSecondaryLoad,
                          rightSecondaryLoad,
                          positiveLineLocal,
                          positiveLineRemote,
                          negativeLineLocal,
                          negativeLineRemote,
                          positiveGroundResistor,
                          negativeGroundResistor,
                          positiveGroundCapacitor,
                          negativeGroundCapacitor,
                          rectifier,
                          inverter});

  auto logger = DataLogger::make(simName);
  logger->logAttribute("left_400kV.v", leftPrimary->attribute("v"));
  logger->logAttribute("left_transformer_400kV.v",
                       leftTransformerPrimary->attribute("v"));
  logger->logAttribute("left_333kV.v", leftSecondary->attribute("v"));
  logger->logAttribute("right_333kV.v", rightSecondary->attribute("v"));
  logger->logAttribute("right_transformer_400kV.v",
                       rightTransformerPrimary->attribute("v"));
  logger->logAttribute("right_400kV.v", rightPrimary->attribute("v"));
  logger->logAttribute("right_grid_after_R.v",
                       rightGridAfterResistance->attribute("v"));
  logger->logAttribute("right_grid_R.i",
                       rightGridResistance->attribute("i_intf"));
  logger->logAttribute("right_grid_L.i",
                       rightGridInductance->attribute("i_intf"));
  logger->logAttribute("transformer_rectifier_R.i",
                       leftTransformerResistance->attribute("i_intf"));
  logger->logAttribute("transformer_inverter_R.i",
                       rightTransformerResistance->attribute("i_intf"));
  logger->logAttribute("transformer_rectifier.i",
                       leftTransformer->attribute("i_intf"));
  logger->logAttribute("transformer_inverter.i",
                       rightTransformer->attribute("i_intf"));

  logger->logAttribute("dc_rectifier_pos.v", rectifierPositive->attribute("v"));
  logger->logAttribute("dc_rectifier_neg.v", rectifierNegative->attribute("v"));
  logger->logAttribute("dc_inverter_pos.v", inverterPositive->attribute("v"));
  logger->logAttribute("dc_inverter_neg.v", inverterNegative->attribute("v"));

  logger->logAttribute("dc_pos_local.i",
                       positiveLineLocal->attribute("i_intf"));
  logger->logAttribute("dc_pos_remote.i",
                       positiveLineRemote->attribute("i_intf"));
  logger->logAttribute("dc_neg_local.i",
                       negativeLineLocal->attribute("i_intf"));
  logger->logAttribute("dc_neg_remote.i",
                       negativeLineRemote->attribute("i_intf"));
  logger->logAttribute("dc_pos_ground_r.i",
                       positiveGroundResistor->attribute("i_intf"));
  logger->logAttribute("dc_neg_ground_r.i",
                       negativeGroundResistor->attribute("i_intf"));
  logger->logAttribute("dc_pos_ground_c.i",
                       positiveGroundCapacitor->attribute("i_intf"));
  logger->logAttribute("dc_neg_ground_c.i",
                       negativeGroundCapacitor->attribute("i_intf"));

  for (const auto &[prefix, mmc] :
       std::vector<std::pair<String, EMT::Ph3::SSN_MMC::Ptr>>{
           {"rectifier", rectifier}, {"inverter", inverter}}) {
    logger->logAttribute(prefix + ".vdc", mmc->dcVoltageAttribute());
    logger->logAttribute(prefix + ".idc", mmc->dcCurrentAttribute());
    logger->logAttribute(prefix + ".p", mmc->activePowerAttribute());
    logger->logAttribute(prefix + ".q", mmc->reactivePowerAttribute());
    logger->logAttribute(prefix + ".energy", mmc->storedEnergyAttribute());
    logger->logAttribute(prefix + ".modulation",
                         mmc->appliedModulationAttribute());
    logger->logAttribute(prefix + ".eq_residual",
                         mmc->attribute("equilibrium_residual_norm"));
    logger->logAttribute(prefix + ".state_norm", mmc->attribute("state_norm"));
    logger->logAttribute(prefix + ".derivative_norm",
                         mmc->attribute("state_derivative_norm"));
    logger->logAttribute(prefix + ".pdc", mmc->attribute("p_dc"));
    logger->logAttribute(prefix + ".power_balance_error",
                         mmc->attribute("power_balance_error"));
    logger->logAttribute(prefix + ".pll_error", mmc->attribute("pll_error"));
    logger->logAttribute(prefix + ".pll_angle_deviation",
                         mmc->attribute("pll_angle_deviation"));
    logger->logAttribute(prefix + ".i_delta_d_ref",
                         mmc->attribute("i_delta_d_ref"));
    logger->logAttribute(prefix + ".i_delta_q_ref",
                         mmc->attribute("i_delta_q_ref"));
    logger->logAttribute(prefix + ".applied_differential_voltage",
                         mmc->appliedDifferentialVoltageAttribute());
    logger->logAttribute(prefix + ".applied_common_mode_voltage",
                         mmc->appliedCommonModeVoltageAttribute());
    logger->logAttribute(prefix + ".external_common_mode_voltage",
                         mmc->externalCommonModeVoltageAttribute());
    logger->logAttribute(prefix + ".realized_converter_voltage",
                         mmc->realizedConverterVoltageAttribute());
  }

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

  // Verify that the explicit right-grid series R and L start with the same
  // instantaneous phase currents before any controller or MNA step is taken.
  const Matrix initialRightGridResistanceCurrent =
      rightGridResistance->intfCurrent();
  const Matrix initialRightGridInductanceCurrent =
      rightGridInductance->intfCurrent();
  const Real initialRightGridCurrentMismatch =
      (initialRightGridResistanceCurrent - initialRightGridInductanceCurrent)
          .cwiseAbs()
          .maxCoeff();

  SPDLOG_INFO("Initial right-grid branch: R_i={}, L_i={}, max_mismatch={} A",
              Logger::matrixToString(initialRightGridResistanceCurrent),
              Logger::matrixToString(initialRightGridInductanceCurrent),
              initialRightGridCurrentMismatch);

  require(initialRightGridCurrentMismatch < 1e-6,
          "Right-grid resistor and inductor hot-start currents disagree.");

  Matrix expectedRightGridCurrent = Matrix::Zero(3, 1);
  const Complex expectedRightGridPeakA =
      std::sqrt(2.0) * rightGridCurrentPhaseRms;
  expectedRightGridCurrent(0, 0) = -expectedRightGridPeakA.real();
  expectedRightGridCurrent(1, 0) =
      -(expectedRightGridPeakA * SHIFT_TO_PHASE_B).real();
  expectedRightGridCurrent(2, 0) =
      -(expectedRightGridPeakA * SHIFT_TO_PHASE_C).real();

  const Real initialRightGridOperatingPointError =
      (initialRightGridInductanceCurrent - expectedRightGridCurrent)
          .cwiseAbs()
          .maxCoeff();
  SPDLOG_INFO(
      "Initial right-grid expected current={}, max_operating_point_error={} A",
      Logger::matrixToString(expectedRightGridCurrent),
      initialRightGridOperatingPointError);
  require(initialRightGridOperatingPointError < 1e-3,
          "Right-grid inductor is not initialized at the analytical "
          "operating current.");

  const Real initialRectifierVdc = **rectifier->dcVoltageAttribute();
  const Real initialInverterVdc = **inverter->dcVoltageAttribute();
  const Real initialRectifierEnergy = **rectifier->storedEnergyAttribute();
  const Real initialInverterEnergy = **inverter->storedEnergyAttribute();

  SPDLOG_INFO("Initialized MATLAB P2P: "
              "REC[P,Q,Vdc,Idc,E,eq_res]=[{},{},{},{},{},{}], "
              "INV[P,Q,Vdc,Idc,E,eq_res]=[{},{},{},{},{},{}]",
              **rectifier->activePowerAttribute(),
              **rectifier->reactivePowerAttribute(), initialRectifierVdc,
              **rectifier->dcCurrentAttribute(), initialRectifierEnergy,
              scalarAttribute(rectifier, "equilibrium_residual_norm"),
              **inverter->activePowerAttribute(),
              **inverter->reactivePowerAttribute(), initialInverterVdc,
              **inverter->dcCurrentAttribute(), initialInverterEnergy,
              scalarAttribute(inverter, "equilibrium_residual_norm"));

  // Validate the complete t=0 sign convention before taking an EMT step.
  // DPsim interface current is positive into the component. MATLAB iDelta and
  // p_ac are positive for AC export; p_dc is positive for DC absorption.
  const auto validateMmcPortSigns = [&](const String &name,
                                        const EMT::Ph3::SSN_MMC::Ptr &mmc) {
    const Matrix interfaceVoltage = **mmc->interfaceVoltageAttribute();
    const Matrix interfaceCurrent = **mmc->interfaceCurrentAttribute();

    const Real iDeltaD = scalarAttribute(mmc, "i_delta_d");
    const Real iDeltaQ = scalarAttribute(mmc, "i_delta_q");
    const Real gridAngle = scalarAttribute(mmc, "grid_angle");
    const Matrix expectedAcInterfaceCurrent =
        -dqToAbcPeak(iDeltaD, iDeltaQ, gridAngle);
    const Real acCurrentError =
        (interfaceCurrent.block(0, 0, 3, 1) - expectedAcInterfaceCurrent)
            .cwiseAbs()
            .maxCoeff();

    const Real pAcFromMna =
        -interfaceVoltage.block(0, 0, 3, 1)
             .cwiseProduct(interfaceCurrent.block(0, 0, 3, 1))
             .sum();
    const Real pAcAttribute = **mmc->activePowerAttribute();

    const Real pDcFromMna = interfaceVoltage(3, 0) * interfaceCurrent(3, 0) +
                            interfaceVoltage(4, 0) * interfaceCurrent(4, 0);
    const Real pDcAttribute = scalarAttribute(mmc, "p_dc");

    SPDLOG_INFO("{} initial port signs: max_ac_current_error={} A, "
                "P_ac[MNA,attribute]=[{},{}] W, "
                "P_dc[MNA,attribute]=[{},{}] W",
                name, acCurrentError, pAcFromMna, pAcAttribute, pDcFromMna,
                pDcAttribute);

    require(acCurrentError < 1e-6,
            name + " AC interface-current sign is inconsistent.");
    require(std::abs(pAcFromMna - pAcAttribute) <
                std::max<Real>(1.0, 1e-9 * std::abs(pAcAttribute)),
            name + " AC port power sign is inconsistent.");
    require(std::abs(pDcFromMna - pDcAttribute) <
                std::max<Real>(1.0, 1e-9 * std::abs(pDcAttribute)),
            name + " DC port power sign is inconsistent.");
  };

  validateMmcPortSigns("MMC_RECTIFIER", rectifier);
  validateMmcPortSigns("MMC_INVERTER", inverter);

  printLocalJacobianDiagnostic("MMC_RECTIFIER", rectifier);
  printLocalJacobianDiagnostic("MMC_INVERTER", inverter);

  Matrix heldRectifierDifferentialVoltage =
      **rectifier->appliedDifferentialVoltageAttribute();
  Matrix heldInverterDifferentialVoltage =
      **inverter->appliedDifferentialVoltageAttribute();
  Matrix heldRectifierCommonModeVoltage =
      **rectifier->appliedCommonModeVoltageAttribute();
  Matrix heldInverterCommonModeVoltage =
      **inverter->appliedCommonModeVoltageAttribute();

  if (diagnosticMode != ClosedLoopDiagnosticMode::InternalControllers) {
    rectifier->setExternalDifferentialVoltageCommand(
        heldRectifierDifferentialVoltage(0, 0),
        heldRectifierDifferentialVoltage(1, 0));
    inverter->setExternalDifferentialVoltageCommand(
        heldInverterDifferentialVoltage(0, 0),
        heldInverterDifferentialVoltage(1, 0));

    if (diagnosticMode ==
        ClosedLoopDiagnosticMode::HoldInitialFullConverterVoltage) {
      rectifier->setExternalCommonModeVoltageCommand(
          heldRectifierCommonModeVoltage(0, 0),
          heldRectifierCommonModeVoltage(1, 0),
          heldRectifierCommonModeVoltage(2, 0));
      inverter->setExternalCommonModeVoltageCommand(
          heldInverterCommonModeVoltage(0, 0),
          heldInverterCommonModeVoltage(1, 0),
          heldInverterCommonModeVoltage(2, 0));

      rectifier->setControlSource(
          EMT::Ph3::SSN_MMC::ControlSource::ExternalFullConverterVoltage);
      inverter->setControlSource(
          EMT::Ph3::SSN_MMC::ControlSource::ExternalFullConverterVoltage);
    } else {
      rectifier->setControlSource(
          EMT::Ph3::SSN_MMC::ControlSource::ExternalDifferentialVoltage);
      inverter->setControlSource(
          EMT::Ph3::SSN_MMC::ControlSource::ExternalDifferentialVoltage);
    }

    const Matrix rectifierDerivativeAfterHold = rectifier->getStateDerivative();
    const Matrix inverterDerivativeAfterHold = inverter->getStateDerivative();

    SPDLOG_INFO("Held converter-voltage diagnostic enabled: mode={}, "
                "REC_vDelta=[{},{}] V, REC_vSigma=[{},{},{}] V, "
                "INV_vDelta=[{},{}] V, INV_vSigma=[{},{},{}] V, "
                "REC_max_electrical_dx={} 1/s-equivalent, "
                "INV_max_electrical_dx={} 1/s-equivalent",
                closedLoopDiagnosticModeName(diagnosticMode),
                heldRectifierDifferentialVoltage(0, 0),
                heldRectifierDifferentialVoltage(1, 0),
                heldRectifierCommonModeVoltage(0, 0),
                heldRectifierCommonModeVoltage(1, 0),
                heldRectifierCommonModeVoltage(2, 0),
                heldInverterDifferentialVoltage(0, 0),
                heldInverterDifferentialVoltage(1, 0),
                heldInverterCommonModeVoltage(0, 0),
                heldInverterCommonModeVoltage(1, 0),
                heldInverterCommonModeVoltage(2, 0),
                rectifierDerivativeAfterHold.topRows(12).cwiseAbs().maxCoeff(),
                inverterDerivativeAfterHold.topRows(12).cwiseAbs().maxCoeff());
  }

  const auto rectifierStateNames = rectifier->getLocalStateNames();
  const auto inverterStateNames = inverter->getLocalStateNames();

  constexpr std::size_t traceCapacity = 1000;
  std::deque<P2PDiagnosticSample> trace;
  auto appendTrace = [&](P2PDiagnosticSample sample) {
    trace.push_back(std::move(sample));
    if (trace.size() > traceCapacity)
      trace.pop_front();
  };

  auto capture = [&](Real time, Bool evaluateDerivative) {
    return captureP2PSample(
        time, rectifier, inverter, positiveLineLocal, positiveLineRemote,
        negativeLineLocal, negativeLineRemote, positiveGroundResistor,
        negativeGroundResistor, positiveGroundCapacitor,
        negativeGroundCapacitor, rectifierPositive, rectifierNegative,
        positiveCableMid, negativeCableMid, inverterPositive, inverterNegative,
        cableSectionResistance, p.groundingResistance, evaluateDerivative);
  };

  appendTrace(capture(0.0, true));
  printCompactDiagnostic(trace.back(), initialRectifierEnergy,
                         initialInverterEnergy);

  // A true hot start requires the dynamic cable currents and the converter
  // DC-port currents to agree before the first simulation step.
  const auto &initialSample = trace.back();
  const Real initialCurrentTolerance =
      std::max<Real>(1e-6, 1e-9 * std::abs(op.dcCurrent));
  require(std::abs(initialSample.positiveLocalCurrent - op.dcCurrent) <
                  initialCurrentTolerance &&
              std::abs(initialSample.positiveRemoteCurrent - op.dcCurrent) <
                  initialCurrentTolerance &&
              std::abs(initialSample.negativeLocalCurrent - op.dcCurrent) <
                  initialCurrentTolerance &&
              std::abs(initialSample.negativeRemoteCurrent - op.dcCurrent) <
                  initialCurrentTolerance,
          "DC cable series-current hot start is inconsistent with the MMC "
          "operating point.");
  require(initialSample.maximumKclResidual < 1e-6,
          "Initial DC-node KCL residual exceeds 1e-6 A.");
  require(std::abs(initialSample.dcPowerResidual) < 10.0,
          "Initial DC-network power residual exceeds 10 W.");

  SPDLOG_INFO("Permanent inverter sampled Pref/Vd control enabled: "
              "Pref={} W, filter=[fn={} Hz,zeta=1], TsControl={} s, "
              "initial_Id_ref={} A",
              op.inverterPower.real(), p.dqMeasurementFilterFrequency,
              p.controlTimeStep, initialInverterIdReference);

  Bool matlabStepAt4SecondsApplied = false;
  Bool matlabStepAt6SecondsApplied = false;
  Bool matlabStepAt8SecondsApplied = false;

  Bool failed = false;
  String failureMessage;
  Real failureTime = 0.0;
  UInt stepCount = 0;

  Real maximumRectifierVdcDeviation = 0.0;
  Real maximumInverterVdcDeviation = 0.0;
  Real maximumRectifierEnergyDeviation = 0.0;
  Real maximumInverterEnergyDeviation = 0.0;
  Real maximumKclResidual = 0.0;
  Real maximumDcPowerResidual = 0.0;

  sim.start();
  while (sim.time() < sim.finalTime()) {
    try {
      if (diagnosticMode == ClosedLoopDiagnosticMode::InternalControllers) {
        const Real currentTime = sim.time();

        if (!matlabStepAt4SecondsApplied &&
            currentTime >= 4.0 - 0.5 * timeStep) {
          inverter->setActivePowerFeedforwardReference(0.20 * p.ratedPower);
          matlabStepAt4SecondsApplied = true;
          SPDLOG_INFO("MATLAB Pref step at t={} s: Pref={} W", currentTime,
                      0.20 * p.ratedPower);
        }

        if (!matlabStepAt6SecondsApplied &&
            currentTime >= 6.0 - 0.5 * timeStep) {
          inverter->setActivePowerFeedforwardReference(-0.35 * p.ratedPower);
          matlabStepAt6SecondsApplied = true;
          SPDLOG_INFO("MATLAB Pref step at t={} s: Pref={} W", currentTime,
                      -0.35 * p.ratedPower);
        }

        if (!matlabStepAt8SecondsApplied &&
            currentTime >= 8.0 - 0.5 * timeStep) {
          inverter->setActivePowerFeedforwardReference(0.50 * p.ratedPower);
          matlabStepAt8SecondsApplied = true;
          SPDLOG_INFO("MATLAB Pref step at t={} s: Pref={} W", currentTime,
                      0.50 * p.ratedPower);
        }
      }

      sim.step();
      ++stepCount;

      P2PDiagnosticSample sample = capture(sim.time(), true);
      appendTrace(std::move(sample));
      const auto &latest = trace.back();

      if (!latest.rectifier.state.allFinite() ||
          !latest.inverter.state.allFinite() ||
          !latest.rectifier.derivative.allFinite() ||
          !latest.inverter.derivative.allFinite()) {
        failed = true;
        failureMessage = "A state or independently evaluated "
                         "derivative became non-finite.";
        failureTime = sim.time();
        break;
      }

      maximumRectifierVdcDeviation =
          std::max(maximumRectifierVdcDeviation,
                   std::abs(latest.rectifier.vdc - initialRectifierVdc));
      maximumInverterVdcDeviation =
          std::max(maximumInverterVdcDeviation,
                   std::abs(latest.inverter.vdc - initialInverterVdc));
      maximumRectifierEnergyDeviation = std::max(
          maximumRectifierEnergyDeviation,
          std::abs(latest.rectifier.storedEnergy - initialRectifierEnergy));
      maximumInverterEnergyDeviation = std::max(
          maximumInverterEnergyDeviation,
          std::abs(latest.inverter.storedEnergy - initialInverterEnergy));
      maximumKclResidual =
          std::max(maximumKclResidual, latest.maximumKclResidual);
      maximumDcPowerResidual =
          std::max(maximumDcPowerResidual, std::abs(latest.dcPowerResidual));

      const Real currentBound = 4.0 * p.ratedPower / p.nominalDcVoltage;
      const Bool rectifierVdcBound =
          std::abs(latest.rectifier.vdc) > 1.5 * p.nominalDcVoltage;
      const Bool inverterVdcBound =
          std::abs(latest.inverter.vdc) > 1.5 * p.nominalDcVoltage;
      const Bool rectifierCurrentBound =
          std::abs(latest.rectifier.idc) > currentBound;
      const Bool inverterCurrentBound =
          std::abs(latest.inverter.idc) > currentBound;
      const Bool rectifierEnergyBound =
          latest.rectifier.storedEnergy < 0.5 * initialRectifierEnergy ||
          latest.rectifier.storedEnergy > 1.5 * initialRectifierEnergy;
      const Bool inverterEnergyBound =
          latest.inverter.storedEnergy < 0.5 * initialInverterEnergy ||
          latest.inverter.storedEnergy > 1.5 * initialInverterEnergy;
      const Bool rectifierCommonModeBound =
          std::abs(latest.rectifierCommonMode) > 0.25 * p.nominalDcVoltage;
      const Bool inverterCommonModeBound =
          std::abs(latest.inverterCommonMode) > 0.25 * p.nominalDcVoltage;

      const Bool physicalBoundExceeded =
          rectifierVdcBound || inverterVdcBound || rectifierCurrentBound ||
          inverterCurrentBound || rectifierEnergyBound || inverterEnergyBound ||
          rectifierCommonModeBound || inverterCommonModeBound;

      if (physicalBoundExceeded) {
        std::ostringstream reason;
        reason << "Physical diagnostic bound exceeded:";
        if (rectifierVdcBound)
          reason << " rectifier_Vdc=" << latest.rectifier.vdc;
        if (inverterVdcBound)
          reason << " inverter_Vdc=" << latest.inverter.vdc;
        if (rectifierCurrentBound)
          reason << " rectifier_Idc=" << latest.rectifier.idc;
        if (inverterCurrentBound)
          reason << " inverter_Idc=" << latest.inverter.idc;
        if (rectifierEnergyBound)
          reason << " rectifier_E=" << latest.rectifier.storedEnergy;
        if (inverterEnergyBound)
          reason << " inverter_E=" << latest.inverter.storedEnergy;
        if (rectifierCommonModeBound)
          reason << " rectifier_Vcm=" << latest.rectifierCommonMode;
        if (inverterCommonModeBound)
          reason << " inverter_Vcm=" << latest.inverterCommonMode;

        failed = true;
        failureMessage = reason.str();
        failureTime = sim.time();
        break;
      }

      if (stepCount % reportIntervalSteps == 0)
        printCompactDiagnostic(latest, initialRectifierEnergy,
                               initialInverterEnergy);
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
      writeFailureTrace(traceFile, trace, rectifierStateNames,
                        inverterStateNames);
      SPDLOG_CRITICAL("MATLAB P2P mode={} failed at t={} s after {} "
                      "successful steps: {}. "
                      "The last {} samples were written to {}.",
                      closedLoopDiagnosticModeName(diagnosticMode), failureTime,
                      stepCount, failureMessage, trace.size(), traceFile);
    } catch (const std::exception &traceException) {
      SPDLOG_CRITICAL("MATLAB P2P mode={} failed at t={} s after {} "
                      "successful steps: {}. "
                      "Writing the trace also failed: {}",
                      closedLoopDiagnosticModeName(diagnosticMode), failureTime,
                      stepCount, failureMessage, traceException.what());
    }

    if (!trace.empty()) {
      printCompactDiagnostic(trace.back(), initialRectifierEnergy,
                             initialInverterEnergy);
      dumpMmcStateToConsole("MMC_RECTIFIER", trace.back().rectifier,
                            rectifierStateNames);
      dumpMmcStateToConsole("MMC_INVERTER", trace.back().inverter,
                            inverterStateNames);
    }
    return 2;
  }

  const auto finalSample = capture(sim.time(), true);
  SPDLOG_INFO("MATLAB P2P mode={} hot-start drift after {} s: "
              "REC[max_dVdc={},final_dVdc={},max_dE={},final_dE={},"
              "final_P={},final_Q={}], "
              "INV[max_dVdc={},final_dVdc={},max_dE={},final_dE={},"
              "final_P={},final_Q={}], "
              "DC[max_KCL={},max_power_residual={},final_Vcm=[{},{}]]",
              closedLoopDiagnosticModeName(diagnosticMode), finalTime,
              maximumRectifierVdcDeviation,
              finalSample.rectifier.vdc - initialRectifierVdc,
              maximumRectifierEnergyDeviation,
              finalSample.rectifier.storedEnergy - initialRectifierEnergy,
              finalSample.rectifier.activePower,
              finalSample.rectifier.reactivePower, maximumInverterVdcDeviation,
              finalSample.inverter.vdc - initialInverterVdc,
              maximumInverterEnergyDeviation,
              finalSample.inverter.storedEnergy - initialInverterEnergy,
              finalSample.inverter.activePower,
              finalSample.inverter.reactivePower, maximumKclResidual,
              maximumDcPowerResidual, finalSample.rectifierCommonMode,
              finalSample.inverterCommonMode);

  return 0;
}
