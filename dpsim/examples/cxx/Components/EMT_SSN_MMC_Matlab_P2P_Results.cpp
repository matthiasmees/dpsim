// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include "../Examples.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

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

void configureMmcCommon(const EMT::Ph3::SSN_MMC::Ptr &mmc,
                        const MatlabCaseParameters &p,
                        const ControllerGainsSi &gains, const HotStartPoint &op,
                        Real initialActivePower, Real initialReactivePower,
                        Real initialDcCurrent, Real initialAngle, Real theta) {
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
  mmc->setTheta(theta);
}

} // namespace

int main(int argc, char **argv) {
  const Real timeStep = argc > 1 ? std::stod(argv[1]) : 20e-6;
  const Real finalTime = argc > 2 ? std::stod(argv[2]) : 10.0;
  const Real theta = argc > 3 ? std::stod(argv[3]) : 0.5;
  const UInt logEveryNthStep =
      argc > 4 ? static_cast<UInt>(std::stoul(argv[4])) : 10;

  require(timeStep > 0.0, "The EMT time step must be positive.");
  require(finalTime > 0.0, "The final time must be positive.");
  require(theta >= 0.5 && theta <= 1.0, "The SSN theta must be in [0.5,1.0].");
  require(logEveryNthStep > 0,
          "The logging interval must be at least one EMT step.");

  const String simName = "EMT_SSN_MMC_Matlab_P2P_Results";
  Logger::setLogDir("logs/" + simName);

  const MatlabCaseParameters p;
  const ControllerGainsSi gains = convertControllerGainsToSi(p);
  const HotStartPoint op = makeHotStartPoint(p);

  const Real omega = 2.0 * PI * p.frequency;
  const Real cableSectionResistance = p.cableResistancePerPole / 2.0;
  const Real cableSectionInductance = p.cableInductancePerPole / 2.0;

  SPDLOG_INFO("MATLAB P2P result run: dt={} s, finalTime={} s, theta={}, "
              "logEvery={} steps ({} s), Pref_initial={} W, "
              "Vdc=[{},{}] V, Idc={} A",
              timeStep, finalTime, theta, logEveryNthStep,
              logEveryNthStep * timeStep, op.inverterPower.real(),
              op.rectifierDcVoltage, op.inverterDcVoltage, op.dcCurrent);

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

  // DC cable current orientations:
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
                     std::arg(op.rectifier.secondaryPhaseRms), theta);
  configureMmcCommon(inverter, p, gains, op, op.inverterPower.real(),
                     op.inverterPower.imag(), op.dcCurrent,
                     std::arg(op.inverter.secondaryPhaseRms), theta);

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
  //   const Real initialInverterIdReference =
  //       (2.0 / 3.0) *
  //       op.inverterPower.real() /
  //       (std::sqrt(2.0) *
  //        std::abs(
  //            op.inverter.secondaryPhaseRms));
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

  // The reference is logged together with the measured quantities so the
  // complete MATLAB step sequence can be plotted from one CSV file.
  auto inverterActivePowerReference =
      AttributeStatic<Real>::make(op.inverterPower.real());

  auto logger = DataLogger::make(simName);
  logger->logAttribute("reference.p_inverter", inverterActivePowerReference);

  // AC terminal waveforms.
  logger->logAttribute("ac.rectifier_terminal.v_abc",
                       leftSecondary->attribute("v"));
  logger->logAttribute("ac.inverter_terminal.v_abc",
                       rightSecondary->attribute("v"));
  logger->logAttribute("ac.right_grid.i_abc",
                       rightGridInductance->attribute("i_intf"));

  // DC pole quantities.
  logger->logAttribute("dc.rectifier_positive.v",
                       rectifierPositive->attribute("v"));
  logger->logAttribute("dc.rectifier_negative.v",
                       rectifierNegative->attribute("v"));
  logger->logAttribute("dc.inverter_positive.v",
                       inverterPositive->attribute("v"));
  logger->logAttribute("dc.inverter_negative.v",
                       inverterNegative->attribute("v"));
  logger->logAttribute("dc.positive_pole.i",
                       positiveLineLocal->attribute("i_intf"));
  logger->logAttribute("dc.negative_pole.i",
                       negativeLineLocal->attribute("i_intf"));

  for (const auto &[prefix, mmc] :
       std::vector<std::pair<String, EMT::Ph3::SSN_MMC::Ptr>>{
           {"rectifier", rectifier}, {"inverter", inverter}}) {
    // Powers, DC operating point and stored energy.
    logger->logAttribute(prefix + ".p_ac", mmc->activePowerAttribute());
    logger->logAttribute(prefix + ".q_ac", mmc->reactivePowerAttribute());
    logger->logAttribute(prefix + ".p_dc", mmc->attribute("p_dc"));
    logger->logAttribute(prefix + ".vdc", mmc->dcVoltageAttribute());
    logger->logAttribute(prefix + ".idc", mmc->dcCurrentAttribute());
    logger->logAttribute(prefix + ".energy", mmc->storedEnergyAttribute());

    // PLL-frame voltages and currents.
    logger->logAttribute(prefix + ".v_grid_d", mmc->attribute("v_grid_d"));
    logger->logAttribute(prefix + ".v_grid_q", mmc->attribute("v_grid_q"));
    logger->logAttribute(prefix + ".i_delta_d", mmc->attribute("i_delta_d"));
    logger->logAttribute(prefix + ".i_delta_q", mmc->attribute("i_delta_q"));
    logger->logAttribute(prefix + ".i_sigma_z", mmc->attribute("i_sigma_z"));

    // Controller references and filtered measurements.
    logger->logAttribute(prefix + ".i_delta_d_ref",
                         mmc->attribute("i_delta_d_ref"));
    logger->logAttribute(prefix + ".i_delta_q_ref",
                         mmc->attribute("i_delta_q_ref"));
    logger->logAttribute(prefix + ".i_sigma_z_ref",
                         mmc->attribute("i_sigma_z_ref"));
    logger->logAttribute(prefix + ".p_filtered", mmc->attribute("p_filtered"));
    logger->logAttribute(prefix + ".q_filtered", mmc->attribute("q_filtered"));
    logger->logAttribute(prefix + ".vdc_filtered",
                         mmc->attribute("vdc_filtered"));
    logger->logAttribute(prefix + ".v_d_feedforward_filtered",
                         mmc->attribute("v_d_feedforward_filtered"));
    logger->logAttribute(prefix + ".i_delta_d_feedforward_held",
                         mmc->attribute("i_delta_d_feedforward_held"));

    // Synchronization and converter commands.
    logger->logAttribute(prefix + ".pll_frequency",
                         mmc->attribute("pll_frequency"));
    logger->logAttribute(prefix + ".pll_angle_deviation",
                         mmc->attribute("pll_angle_deviation"));
    logger->logAttribute(prefix + ".modulation",
                         mmc->appliedModulationAttribute());
    logger->logAttribute(prefix + ".differential_voltage",
                         mmc->appliedDifferentialVoltageAttribute());

    const auto commonModeVoltage = mmc->appliedCommonModeVoltageAttribute();
    logger->logAttribute(prefix + ".common_mode_voltage_z",
                         commonModeVoltage->deriveCoeff<Real>(2, 0));
  }

  Simulation sim(simName, Logger::Level::off);
  sim.setSystem(system);
  sim.setDomain(Domain::EMT);
  sim.setTimeStep(timeStep);
  sim.setFinalTime(finalTime);
  sim.setSolverType(Solver::Type::MNA);
  sim.doSystemMatrixRecomputation(true);
  sim.doInitFromNodesAndTerminals(true);
  sim.initialize();

  SPDLOG_INFO(
      "Initialized operating point: "
      "REC[P={},Q={},Vdc={},Idc={}], "
      "INV[P={},Q={},Vdc={},Idc={}]. "
      "Results are logged under logs/{}/.",
      **rectifier->activePowerAttribute(),
      **rectifier->reactivePowerAttribute(), **rectifier->dcVoltageAttribute(),
      **rectifier->dcCurrentAttribute(), **inverter->activePowerAttribute(),
      **inverter->reactivePowerAttribute(), **inverter->dcVoltageAttribute(),
      **inverter->dcCurrentAttribute(), simName);

  Bool matlabStepAt4SecondsApplied = false;
  Bool matlabStepAt6SecondsApplied = false;
  Bool matlabStepAt8SecondsApplied = false;

  logger->start();
  logger->log(0.0, 0);

  UInt stepCount = 0;
  sim.start();
  while (sim.time() < sim.finalTime()) {
    const Real currentTime = sim.time();

    if (!matlabStepAt4SecondsApplied && currentTime >= 4.0 - 0.5 * timeStep) {
      const Real reference = 0.20 * p.ratedPower;
      **inverterActivePowerReference = reference;
      inverter->setActivePowerFeedforwardReference(reference);
      matlabStepAt4SecondsApplied = true;
      SPDLOG_INFO("MATLAB Pref step at t={} s: Pref={} W", currentTime,
                  reference);
    }

    if (!matlabStepAt6SecondsApplied && currentTime >= 6.0 - 0.5 * timeStep) {
      const Real reference = -0.35 * p.ratedPower;
      **inverterActivePowerReference = reference;
      inverter->setActivePowerFeedforwardReference(reference);
      matlabStepAt6SecondsApplied = true;
      SPDLOG_INFO("MATLAB Pref step at t={} s: Pref={} W", currentTime,
                  reference);
    }

    if (!matlabStepAt8SecondsApplied && currentTime >= 8.0 - 0.5 * timeStep) {
      const Real reference = 0.50 * p.ratedPower;
      **inverterActivePowerReference = reference;
      inverter->setActivePowerFeedforwardReference(reference);
      matlabStepAt8SecondsApplied = true;
      SPDLOG_INFO("MATLAB Pref step at t={} s: Pref={} W", currentTime,
                  reference);
    }

    sim.step();
    ++stepCount;

    if (stepCount % logEveryNthStep == 0 || sim.time() >= sim.finalTime())
      logger->log(sim.time(), static_cast<Int>(stepCount));
  }
  sim.stop();
  logger->stop();

  SPDLOG_INFO("Simulation completed at t={} s. "
              "Plot data are available under logs/{}/.",
              sim.time(), simName);

  return 0;
}
