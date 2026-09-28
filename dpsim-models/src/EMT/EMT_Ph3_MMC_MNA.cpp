// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <dpsim-models/EMT/EMT_Ph3_MMC_MNA.h>
#include <stdexcept>

using namespace CPS;
using CPS::EMT::Ph3::MMC_MNA;

namespace {
void finite(Real value, const char *name) {
  if (!std::isfinite(value))
    throw std::invalid_argument(String("MMC_MNA: non-finite ") + name);
}
void positive(Real value, const char *name) {
  finite(value, name);
  if (value <= 0.0)
    throw std::invalid_argument(String("MMC_MNA: expected positive ") + name);
}
Real limit(Real value, Real magnitude) {
  return std::max(-magnitude, std::min(magnitude, value));
}
Real conditionalError(Real error, Real raw, Real limited, Real sign = 1.0) {
  return (raw - limited) * sign * error > 0.0 ? 0.0 : error;
}
} // namespace

MMC_MNA::MMC_MNA(String uid, String name, Logger::Level level)
    : MNASimPowerComp<Real>(uid, name, true, true, level) {
  // Auxiliary branch currents and capacitor voltages are scalar MNA unknowns.
  mPhaseType = PhaseType::DC;
  setVirtualNodeNumber(BranchCount);
  mPhaseType = PhaseType::ABC;
  setTerminalNumber(3);
  **mIntfVoltage = Matrix::Zero(5, 1);
  **mIntfCurrent = Matrix::Zero(5, 1);
  mElectricalAttribute =
      mAttributes->create<Matrix>("electrical", Matrix::Zero(12, 1));
  mAcVoltageAttribute =
      mAttributes->create<Matrix>("ac_terminal_voltage", Matrix::Zero(3, 1));
  mAcCurrentAttribute =
      mAttributes->create<Matrix>("ac_terminal_current", Matrix::Zero(3, 1));
  mModulationAttribute =
      mAttributes->create<Matrix>("applied_modulation", Matrix::Zero(5, 1));
  mDifferentialAttribute = mAttributes->create<Matrix>(
      "applied_differential_voltage", Matrix::Zero(2, 1));
  mCommonAttribute = mAttributes->create<Matrix>("applied_common_mode_voltage",
                                                 Matrix::Zero(3, 1));
  for (const auto *key : {"vdc",
                          "vdcp",
                          "vdcn",
                          "idc",
                          "p_ac",
                          "q_ac",
                          "p_dc",
                          "stored_energy",
                          "v_grid_d",
                          "v_grid_q",
                          "i_delta_d",
                          "i_delta_q",
                          "i_sigma_z",
                          "i_sigma_d",
                          "i_sigma_q",
                          "i_delta_d_ref",
                          "i_delta_q_ref",
                          "i_sigma_z_ref",
                          "p_filtered",
                          "q_filtered",
                          "vdc_filtered",
                          "v_d_feedforward_filtered",
                          "i_delta_d_feedforward_held",
                          "pll_frequency",
                          "pll_angle_deviation",
                          "pll_error",
                          "grid_angle",
                          "power_balance_error",
                          "equilibrium_residual_norm"})
    mScalars.emplace(key, mAttributes->create<Real>(key, 0.0));
}

MMC_MNA::PIParameters MMC_MNA::gains(Real kp, Real ki) {
  finite(kp, "proportional gain");
  finite(ki, "integral gain");
  if (kp < 0.0 || ki < 0.0)
    throw std::invalid_argument("MMC_MNA: PI gains must be nonnegative");
  return {kp, ki};
}
void MMC_MNA::setParameters(Real frequency, Real acVoltage, Real dcVoltage,
                            Real armInductance, Real armResistance,
                            Real submoduleCapacitance, UInt submodules,
                            Real reactorInductance, Real reactorResistance) {
  positive(frequency, "frequency");
  positive(acVoltage, "AC voltage");
  positive(dcVoltage, "DC voltage");
  positive(armInductance, "arm inductance");
  positive(submoduleCapacitance, "submodule capacitance");
  finite(armResistance, "arm resistance");
  finite(reactorInductance, "reactor inductance");
  finite(reactorResistance, "reactor resistance");
  if (!submodules || armResistance < 0 || reactorInductance < 0 ||
      reactorResistance < 0)
    throw std::invalid_argument("MMC_MNA: invalid arm/reactor parameters");
  mOmega = 2 * PI * frequency;
  mAcVoltage = acVoltage;
  mDcVoltage = dcVoltage;
  mArmL = armInductance;
  mArmR = armResistance;
  mCapacitance = submoduleCapacitance;
  mSubmodules = submodules;
  mAcL = armInductance / 2 + reactorInductance;
  mAcR = armResistance / 2 + reactorResistance;
  mStorage.fill(submoduleCapacitance / submodules);
  mStorage[IDeltaD] = mStorage[IDeltaQ] = mAcL;
  mStorage[ISigmaZ] = mStorage[ISigmaD] = mStorage[ISigmaQ] = mArmL;
  mParametersSet = true;
}
void MMC_MNA::setInitialAngle(Real value) {
  finite(value, "angle");
  mInitialAngle = value;
}
void MMC_MNA::setInitialOperatingPoint(Real p, Real q) {
  finite(p, "initial P");
  finite(q, "initial Q");
  mInitialP = p;
  mInitialQ = q;
  mLoadedStart = true;
}
void MMC_MNA::setOutputCurrentController(Real kp, Real ki) {
  mOutput = gains(kp, ki);
}
void MMC_MNA::setCirculatingCurrentController(Real kp, Real ki) {
  mCirculating = gains(kp, ki);
}
void MMC_MNA::setZeroSequenceCurrentController(Real kp, Real ki) {
  mZero = gains(kp, ki);
}
void MMC_MNA::setEnergyController(Real kp, Real ki, Bool enabled) {
  mEnergy = gains(kp, ki);
  mEnergyEnabled = enabled;
}
void MMC_MNA::setPLL(Real kp, Real ki, Bool enabled) {
  mPll = gains(kp, ki);
  mPllEnabled = enabled;
}
void MMC_MNA::setActivePowerControl(Real reference, Real kp, Real ki) {
  finite(reference, "P reference");
  mActive = gains(kp, ki);
  mPRef = reference;
  mActiveMode = ActiveControlMode::ActivePower;
}
void MMC_MNA::setActivePowerFeedforwardControl(Real reference, Real cutoff,
                                               Real sampleTime,
                                               Real minimumVd) {
  finite(reference, "P reference");
  positive(cutoff, "filter cutoff");
  positive(sampleTime, "sample time");
  finite(minimumVd, "minimum Vd");
  if (minimumVd < 0)
    throw std::invalid_argument("MMC_MNA: negative minimum Vd");
  mPRef = reference;
  mCutoff = cutoff;
  mSampleTime = sampleTime;
  mMinimumVd = minimumVd;
  mActiveMode = ActiveControlMode::ActivePowerFeedforward;
}
void MMC_MNA::setActivePowerFeedforwardReference(Real reference) {
  finite(reference, "P reference");
  if (mActiveMode != ActiveControlMode::ActivePowerFeedforward)
    throw std::logic_error("MMC_MNA: sampled feedforward is not configured");
  mPRef = reference;
}
void MMC_MNA::setDcVoltageControl(Real reference, Real kp, Real ki) {
  positive(reference, "DC reference");
  mActive = gains(kp, ki);
  mVdcRef = reference;
  mActiveMode = ActiveControlMode::DcVoltage;
}
void MMC_MNA::setDcDroopControl(Real power, Real voltage, Real droop) {
  finite(power, "P reference");
  positive(voltage, "DC reference");
  positive(droop, "droop");
  mPRef = power;
  mVdcRef = voltage;
  mDroop = droop;
  mActiveMode = ActiveControlMode::DcDroop;
}
void MMC_MNA::setActiveControlOpenLoop(Real reference) {
  finite(reference, "Id reference");
  mIdRef = reference;
  mActiveMode = ActiveControlMode::OpenLoop;
}
void MMC_MNA::setReactivePowerControl(Real reference, Real kp, Real ki) {
  finite(reference, "Q reference");
  mReactive = gains(kp, ki);
  mQRef = reference;
  mReactiveMode = ReactiveControlMode::ReactivePower;
}
void MMC_MNA::setAcVoltageControl(Real reference, Real kp, Real ki) {
  positive(reference, "AC reference");
  mReactive = gains(kp, ki);
  mVacRef = reference;
  mReactiveMode = ReactiveControlMode::AcVoltage;
}
void MMC_MNA::setReactiveControlOpenLoop(Real reference) {
  finite(reference, "Iq reference");
  mIqRef = reference;
  mReactiveMode = ReactiveControlMode::OpenLoop;
}
void MMC_MNA::setCirculatingCurrentReferences(Real d, Real q, Real z) {
  finite(d, "Isigma d");
  finite(q, "Isigma q");
  finite(z, "Isigma z");
  mSigmaRef = {d, q, z};
}
void MMC_MNA::setMeasurementFilters(Real dq, Real p, Real q, Real dc, Real ac) {
  for (Real value : {dq, p, q, dc, ac})
    finite(value, "filter time constant");
  mFilterTime = {dq, p, q, dc, ac};
}
void MMC_MNA::setModulationDelay(Real delay, UInt order) {
  finite(delay, "modulation delay");
  if (order != 2 || delay < 0)
    throw std::invalid_argument(
        "MMC_MNA: expected nonnegative second-order delay");
  mDelay = delay;
}
void MMC_MNA::setLimits(Real ac, Real circulating, Real modulation) {
  positive(ac, "AC limit");
  positive(circulating, "circulating limit");
  positive(modulation, "modulation limit");
  mMaxAc = ac;
  mMaxSigma = circulating;
  mMaxModulation = modulation;
}
void MMC_MNA::setTheta(Real theta) {
  finite(theta, "theta");
  if (theta < 0.5 || theta > 1 || mStep > 0)
    throw std::invalid_argument(
        "MMC_MNA: set theta in [0.5,1] before initialization");
  mTheta = theta;
}
void MMC_MNA::setOperatingPointInitialization(Bool enabled, UInt iterations,
                                              Real tolerance) {
  positive(tolerance, "initialization tolerance");
  if (!iterations)
    throw std::invalid_argument("MMC_MNA: zero initialization iterations");
  mSolveStart = enabled;
  mInitialIterations = iterations;
  mInitialTolerance = tolerance;
}

MMC_MNA::Dq MMC_MNA::rotate(Dq value, Real angle) {
  const Real c = std::cos(angle), s = std::sin(angle);
  return {c * value[0] - s * value[1], s * value[0] + c * value[1]};
}
Real MMC_MNA::regularize(Real value) {
  return std::abs(value) >= 1.0 ? value : (value < 0 ? -1.0 : 1.0);
}
MMC_MNA::Measurement MMC_MNA::measure(const Matrix &voltage, Real angle) const {
  Measurement result{0.0, 0.0, voltage(3) - voltage(4)};
  for (UInt phase = 0; phase < 3; ++phase) {
    const Real a = angle - phase * 2 * PI / 3;
    result.vd += (2.0 / 3.0) * std::cos(a) * voltage(phase);
    result.vq -= (2.0 / 3.0) * std::sin(a) * voltage(phase);
  }
  return result;
}
Real MMC_MNA::energy(const Electrical &e) const {
  Real squares = e[VCSigmaZ] * e[VCSigmaZ];
  for (UInt i = VCDeltaD; i < BranchCount; ++i)
    squares += e[i] * e[i];
  return 1.5 * mCapacitance / mSubmodules * squares;
}
MMC_MNA::Memory MMC_MNA::advance(const Memory &m, const Memory &r, Real h) {
  Memory next = m;
  for (UInt i = 0; i < 8; ++i)
    next.integral[i] += h * r.integral[i];
  next.pllIntegral += h * r.pllIntegral;
  next.pllAngle += h * r.pllAngle;
  for (UInt i = 0; i < 10; ++i) {
    next.filter[i] += h * r.filter[i];
    next.delay[i] += h * r.delay[i];
  }
  return next;
}

void MMC_MNA::electricalCouplings(
    const Modulation &m,
    const std::function<void(UInt, UInt, Real)> &add) const {
  const Real d = m[0], q = m[1], sd = m[2], sq = m[3], sz = m[4];
  const Real c = mCapacitance / mSubmodules;
  // Inductor KVL: L di/dt = controlled arm voltages - R i - terminal voltage.
  add(IDeltaD, IDeltaD, -mAcR);
  add(IDeltaD, IDeltaQ, mOmega * mAcL);
  add(IDeltaQ, IDeltaQ, -mAcR);
  add(IDeltaQ, IDeltaD, -mOmega * mAcL);
  add(ISigmaD, ISigmaD, -mArmR);
  add(ISigmaD, ISigmaQ, -2 * mOmega * mArmL);
  add(ISigmaQ, ISigmaQ, -mArmR);
  add(ISigmaQ, ISigmaD, 2 * mOmega * mArmL);
  add(ISigmaZ, ISigmaZ, -mArmR);

  add(IDeltaD, VCDeltaD, -sd / 4 - sz / 2);
  add(IDeltaD, VCDeltaQ, sq / 4);
  add(IDeltaD, VCDeltaZd, -sd / 4);
  add(IDeltaD, VCDeltaZq, sq / 4);
  add(IDeltaD, VCSigmaD, -d / 4);
  add(IDeltaD, VCSigmaQ, q / 4);
  add(IDeltaD, VCSigmaZ, -d / 2);
  add(IDeltaQ, VCDeltaD, sq / 4);
  add(IDeltaQ, VCDeltaQ, sd / 4 - sz / 2);
  add(IDeltaQ, VCDeltaZd, -sq / 4);
  add(IDeltaQ, VCDeltaZq, -sd / 4);
  add(IDeltaQ, VCSigmaD, q / 4);
  add(IDeltaQ, VCSigmaQ, d / 4);
  add(IDeltaQ, VCSigmaZ, -q / 2);
  add(ISigmaD, VCDeltaD, -d / 4);
  add(ISigmaD, VCDeltaQ, q / 4);
  add(ISigmaD, VCDeltaZd, -d / 4);
  add(ISigmaD, VCDeltaZq, -q / 4);
  add(ISigmaD, VCSigmaD, -sz / 2);
  add(ISigmaD, VCSigmaZ, -sd / 2);
  add(ISigmaQ, VCDeltaD, q / 4);
  add(ISigmaQ, VCDeltaQ, d / 4);
  add(ISigmaQ, VCDeltaZd, -q / 4);
  add(ISigmaQ, VCDeltaZq, d / 4);
  add(ISigmaQ, VCSigmaQ, -sz / 2);
  add(ISigmaQ, VCSigmaZ, -sq / 2);
  add(ISigmaZ, VCDeltaD, -d / 4);
  add(ISigmaZ, VCDeltaQ, -q / 4);
  add(ISigmaZ, VCSigmaD, -sd / 4);
  add(ISigmaZ, VCSigmaQ, -sq / 4);
  add(ISigmaZ, VCSigmaZ, -sz / 2);

  // Capacitor KCL: (Csm/N) dv/dt = inserted arm currents + frame couplings.
  add(VCSigmaD, ISigmaD, sz / 2);
  add(VCSigmaD, ISigmaZ, sd / 2);
  add(VCSigmaD, IDeltaD, d / 8);
  add(VCSigmaD, IDeltaQ, -q / 8);
  add(VCSigmaD, VCSigmaQ, -2 * c * mOmega);
  add(VCSigmaQ, ISigmaQ, sz / 2);
  add(VCSigmaQ, ISigmaZ, sq / 2);
  add(VCSigmaQ, IDeltaD, -q / 8);
  add(VCSigmaQ, IDeltaQ, -d / 8);
  add(VCSigmaQ, VCSigmaD, 2 * c * mOmega);
  add(VCSigmaZ, IDeltaD, d / 8);
  add(VCSigmaZ, IDeltaQ, q / 8);
  add(VCSigmaZ, ISigmaD, sd / 4);
  add(VCSigmaZ, ISigmaQ, sq / 4);
  add(VCSigmaZ, ISigmaZ, sz / 2);
  add(VCDeltaD, ISigmaZ, d / 2);
  add(VCDeltaD, ISigmaD, d / 4);
  add(VCDeltaD, ISigmaQ, -q / 4);
  add(VCDeltaD, IDeltaD, sd / 8 + sz / 4);
  add(VCDeltaD, IDeltaQ, -sq / 8);
  add(VCDeltaD, VCDeltaQ, c * mOmega);
  add(VCDeltaQ, ISigmaZ, q / 2);
  add(VCDeltaQ, ISigmaD, -q / 4);
  add(VCDeltaQ, ISigmaQ, -d / 4);
  add(VCDeltaQ, IDeltaD, -sq / 8);
  add(VCDeltaQ, IDeltaQ, -sd / 8 + sz / 4);
  add(VCDeltaQ, VCDeltaD, -c * mOmega);
  add(VCDeltaZd, IDeltaD, sd / 8);
  add(VCDeltaZd, IDeltaQ, sq / 8);
  add(VCDeltaZd, ISigmaD, d / 4);
  add(VCDeltaZd, ISigmaQ, q / 4);
  add(VCDeltaZd, VCDeltaZq, 3 * c * mOmega);
  add(VCDeltaZq, IDeltaD, -sq / 8);
  add(VCDeltaZq, IDeltaQ, sd / 8);
  add(VCDeltaZq, ISigmaD, q / 4);
  add(VCDeltaZq, ISigmaQ, -d / 4);
  add(VCDeltaZq, VCDeltaZd, -3 * c * mOmega);
}

MMC_MNA::Electrical MMC_MNA::electricalBalance(const Electrical &e,
                                               const Modulation &m,
                                               const Measurement &v) const {
  Electrical balance{};
  electricalCouplings(m, [&](UInt row, UInt col, Real value) {
    balance[row] += value * e[col];
  });
  balance[IDeltaD] -= v.vd;
  balance[IDeltaQ] -= v.vq;
  balance[ISigmaZ] += 0.5 * v.vdc;
  return balance;
}

MMC_MNA::Command MMC_MNA::control(const Electrical &e, const Memory &memory,
                                  const Measurement &v) const {
  Command out;
  auto &rate = out.rate;
  const Dq rawControl = rotate({v.vd, v.vq}, -memory.pllAngle);
  out.pllError = rawControl[1];
  if (mPllEnabled) {
    rate.pllIntegral = out.pllError;
    out.deltaOmega = mPll.kp * out.pllError + mPll.ki * memory.pllIntegral;
    rate.pllAngle = out.deltaOmega;
  }
  const Real omega = mOmega + out.deltaOmega;
  auto firstOrder = [&](UInt i, Real input, Real tau) {
    if (tau <= 0)
      return input;
    rate.filter[i] = (input - memory.filter[i]) / tau;
    return memory.filter[i];
  };
  auto secondOrder = [&](UInt i, Real input, Real tau) {
    if (tau <= 0)
      return input;
    rate.filter[i] = (input - memory.filter[i]) / tau;
    rate.filter[i + 1] = (memory.filter[i] - memory.filter[i + 1]) / tau;
    return memory.filter[i + 1];
  };
  const Real vd = firstOrder(0, rawControl[0], mFilterTime[0]);
  const Real vq = firstOrder(1, rawControl[1], mFilterTime[0]);
  const Real p = secondOrder(2, 1.5 * (v.vd * e[IDeltaD] + v.vq * e[IDeltaQ]),
                             mFilterTime[1]);
  const Real q = secondOrder(4, 1.5 * (v.vq * e[IDeltaD] - v.vd * e[IDeltaQ]),
                             mFilterTime[2]);
  const Real dc = secondOrder(6, regularize(v.vdc), mFilterTime[3]);
  const Real vac = secondOrder(8, 1.5 * std::hypot(vd, vq), mFilterTime[4]);
  out.pFiltered = p;
  out.qFiltered = q;
  out.vdcFiltered = dc;
  Real id = mIdRef, iq = mIqRef, pError = 0.0, qError = 0.0, qSign = 1.0;
  switch (mActiveMode) {
  case ActiveControlMode::ActivePower:
    pError = mPRef - p;
    id = mActive.kp * pError + mActive.ki * memory.integral[Active];
    break;
  case ActiveControlMode::DcVoltage:
    pError = dc - mVdcRef;
    id = mActive.kp * pError + mActive.ki * memory.integral[Active];
    break;
  case ActiveControlMode::ActivePowerFeedforward:
    id = mHeldId;
    break;
  case ActiveControlMode::DcDroop:
    id = (mPRef + (mVdcRef - dc) / mDroop) / std::max(1.0, std::abs(vd));
    break;
  case ActiveControlMode::OpenLoop:
    break;
  }
  switch (mReactiveMode) {
  case ReactiveControlMode::ReactivePower:
    qError = mQRef - q;
    qSign = -1.0;
    iq = -(mReactive.kp * qError + mReactive.ki * memory.integral[Reactive]);
    break;
  case ReactiveControlMode::AcVoltage:
    qError = mVacRef - vac;
    iq = mReactive.kp * qError + mReactive.ki * memory.integral[Reactive];
    break;
  case ReactiveControlMode::OpenLoop:
    break;
  }
  const Real scale =
      std::hypot(id, iq) > mMaxAc ? mMaxAc / std::hypot(id, iq) : 1.0;
  out.idRef = scale * id;
  out.iqRef = scale * iq;
  if (mActive.ki > 0)
    rate.integral[Active] = conditionalError(pError, id, out.idRef);
  if (mReactive.ki > 0)
    rate.integral[Reactive] = conditionalError(qError, iq, out.iqRef, qSign);
  const Dq current = rotate({e[IDeltaD], e[IDeltaQ]}, -memory.pllAngle);
  const Dq reference = rotate({out.idRef, out.iqRef}, -memory.pllAngle);
  const Dq error{reference[0] - current[0], reference[1] - current[1]};
  rate.integral[OccD] = error[0];
  rate.integral[OccQ] = error[1];
  out.differential =
      rotate({vd + mOutput.kp * error[0] + mOutput.ki * memory.integral[OccD] +
                  mAcR * reference[0] - omega * mAcL * reference[1],
              vq + mOutput.kp * error[1] + mOutput.ki * memory.integral[OccQ] +
                  mAcR * reference[1] + omega * mAcL * reference[0]},
             memory.pllAngle);
  const Dq sigma = rotate({e[ISigmaD], e[ISigmaQ]}, -2 * memory.pllAngle);
  const Dq sigmaRef =
      rotate({mSigmaRef[0], mSigmaRef[1]}, -2 * memory.pllAngle);
  rate.integral[CccD] = sigmaRef[0] - sigma[0];
  rate.integral[CccQ] = sigmaRef[1] - sigma[1];
  const Dq common = rotate({-(mCirculating.kp * rate.integral[CccD] +
                              mCirculating.ki * memory.integral[CccD] +
                              2 * omega * mArmL * sigma[1]),
                            -(mCirculating.kp * rate.integral[CccQ] +
                              mCirculating.ki * memory.integral[CccQ] -
                              2 * omega * mArmL * sigma[0])},
                           2 * memory.pllAngle);
  Real izRef = mSigmaRef[2], energyError = 0.0;
  if (mEnergyEnabled) {
    energyError =
        3 * mCapacitance / mSubmodules * mDcVoltage * mDcVoltage - energy(e);
    izRef =
        (mEnergy.kp * energyError + mEnergy.ki * memory.integral[Energy] + p) /
        (3 * regularize(dc));
  }
  out.izRef = limit(izRef, mMaxSigma);
  if (mEnergyEnabled && mEnergy.ki > 0)
    rate.integral[Energy] = conditionalError(energyError, izRef, out.izRef);
  rate.integral[Zcc] = out.izRef - e[ISigmaZ];
  out.common = {common[0], common[1],
                dc / 2 - mZero.kp * rate.integral[Zcc] -
                    mZero.ki * memory.integral[Zcc]};
  const Real dcReg = regularize(dc);
  auto &m = out.modulation;
  m = {-2 * out.differential[0] / dcReg, -2 * out.differential[1] / dcReg,
       2 * out.common[0] / dcReg, 2 * out.common[1] / dcReg,
       2 * out.common[2] / dcReg};
  if (std::hypot(m[0], m[1]) > mMaxModulation) {
    const Real ratio = mMaxModulation / std::hypot(m[0], m[1]);
    const Dq drive =
        rotate({mOutput.ki * error[0], mOutput.ki * error[1]}, memory.pllAngle);
    if (out.differential[0] * drive[0] + out.differential[1] * drive[1] > 0)
      rate.integral[OccD] = rate.integral[OccQ] = 0.0;
    m[0] *= ratio;
    m[1] *= ratio;
  }
  if (std::hypot(m[2], m[3]) > mMaxModulation) {
    const Real ratio = mMaxModulation / std::hypot(m[2], m[3]);
    m[2] *= ratio;
    m[3] *= ratio;
  }
  m[4] = limit(m[4], mMaxModulation);
  if (mDelay > 0) {
    const Real a0 = 12 / (mDelay * mDelay), a1 = 6 / mDelay;
    for (UInt channel = 0; channel < 5; ++channel) {
      const UInt i = 2 * channel;
      rate.delay[i] = memory.delay[i + 1];
      rate.delay[i + 1] =
          -a0 * memory.delay[i] - a1 * memory.delay[i + 1] + m[channel];
      m[channel] -= 2 * a1 * memory.delay[i + 1];
    }
  }
  return out;
}

void MMC_MNA::sampleFeedforward() {
  if (mActiveMode != ActiveControlMode::ActivePowerFeedforward)
    return;
  const UInt stride =
      std::max<UInt>(1, static_cast<UInt>(std::llround(mSampleTime / mStep)));
  if (mSampleCounter == 0) {
    const Real h = stride * mStep, w = 2 * PI * mCutoff, a = std::exp(-w * h);
    const Real error = mFilteredVd - mMeasurement.vd;
    const Real y =
        mMeasurement.vd + a * ((1 + w * h) * error + h * mFilteredVdDerivative);
    mFilteredVdDerivative =
        a * ((1 - w * h) * mFilteredVdDerivative - w * w * h * error);
    mFilteredVd = y;
    const Real minimum =
        mMinimumVd > 0 ? mMinimumVd : 0.1 * std::sqrt(2.0 / 3.0) * mAcVoltage;
    if (!std::isfinite(y) || y <= minimum)
      throw std::runtime_error("MMC_MNA: invalid filtered feedforward voltage");
    mHeldId = (2.0 / 3.0) * mPRef / y;
  }
  mSampleCounter = (mSampleCounter + 1) % stride;
}

void MMC_MNA::validateTerminals() const {
  auto *self = const_cast<MMC_MNA *>(this);
  if (self->node(0)->isGround() || self->node(0)->phaseType() != PhaseType::ABC)
    throw std::invalid_argument(
        "MMC_MNA: terminal 0 requires a non-grounded ABC node");
  for (UInt i = 1; i < 3; ++i)
    if (!self->node(i)->isGround() &&
        self->node(i)->phaseType() != PhaseType::DC)
      throw std::invalid_argument(
          "MMC_MNA: DC terminals require DC nodes or ground");
}

void MMC_MNA::initializeFromNodesAndTerminals(Real) {
  if (!mParametersSet)
    throw std::logic_error("MMC_MNA: physical parameters not set");
  validateTerminals();
  // SimPowerComp initializes ABC attributes with three rows. The mixed AC/DC
  // interface needs five entries after that generic initialization.
  **mIntfVoltage = Matrix::Zero(5, 1);
  **mIntfCurrent = Matrix::Zero(5, 1);
  const Complex va = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
  (**mIntfVoltage)(0) = va.real();
  (**mIntfVoltage)(1) = (va * SHIFT_TO_PHASE_B).real();
  (**mIntfVoltage)(2) = (va * SHIFT_TO_PHASE_C).real();
  (**mIntfVoltage)(3) = initialSingleVoltage(1).real();
  (**mIntfVoltage)(4) = initialSingleVoltage(2).real();
  if (std::abs((**mIntfVoltage)(3) - (**mIntfVoltage)(4)) < 1.0) {
    (**mIntfVoltage)(3) =
        terminalNotGrounded(1)
            ? (terminalNotGrounded(2) ? mDcVoltage / 2 : mDcVoltage)
            : 0;
    (**mIntfVoltage)(4) = (**mIntfVoltage)(3) - mDcVoltage;
  }
  mAngle = mNextAngle = mInitialAngle;
  mStepCount = 0;
  mSampleCounter = 0;
  mMemory = {};
  mElectrical = {};
  mSourceConductance = {};
  mMeasurement = measure(**mIntfVoltage, mAngle);
  initializeOperatingPoint();
  mStepModulation = control(mElectrical, mMemory, mMeasurement).modulation;
  updateAttributes();
}

void MMC_MNA::initializeOperatingPoint() {
  const auto &v = mMeasurement;
  const Real p = mLoadedStart ? mInitialP : mPRef,
             q = mLoadedStart ? mInitialQ : mQRef;
  const Real v2 = v.vd * v.vd + v.vq * v.vq;
  mElectrical[IDeltaD] = mIdRef;
  mElectrical[IDeltaQ] = mIqRef;
  if (v2 > 0 && (mLoadedStart || mActiveMode != ActiveControlMode::OpenLoop ||
                 mReactiveMode == ReactiveControlMode::ReactivePower)) {
    mElectrical[IDeltaD] = (2.0 / 3.0) * (v.vd * p + v.vq * q) / v2;
    mElectrical[IDeltaQ] = (2.0 / 3.0) * (v.vq * p - v.vd * q) / v2;
  }
  mElectrical[ISigmaZ] =
      mEnergyEnabled ? p / (3 * regularize(v.vdc)) : mSigmaRef[2];
  mElectrical[ISigmaD] = mSigmaRef[0];
  mElectrical[ISigmaQ] = mSigmaRef[1];
  mElectrical[VCSigmaZ] = std::abs(regularize(v.vdc));
  if (mActive.ki > 0)
    mMemory.integral[Active] = mElectrical[IDeltaD] / mActive.ki;
  if (mReactive.ki > 0)
    mMemory.integral[Reactive] =
        (mReactiveMode == ReactiveControlMode::ReactivePower ? -1 : 1) *
        mElectrical[IDeltaQ] / mReactive.ki;
  if (mCirculating.ki > 0) {
    mMemory.integral[CccD] = mArmR * mElectrical[ISigmaD] / mCirculating.ki;
    mMemory.integral[CccQ] = mArmR * mElectrical[ISigmaQ] / mCirculating.ki;
  }
  if (mZero.ki > 0)
    mMemory.integral[Zcc] = mArmR * mElectrical[ISigmaZ] / mZero.ki;
  if (mEnergyEnabled && mEnergy.ki > 0)
    mMemory.integral[Energy] =
        (3 * v.vdc * mElectrical[ISigmaZ] - p) / mEnergy.ki;
  if (mActiveMode == ActiveControlMode::ActivePowerFeedforward) {
    const Real minimum =
        mMinimumVd > 0 ? mMinimumVd : 0.1 * std::sqrt(2.0 / 3.0) * mAcVoltage;
    if (v.vd <= minimum)
      throw std::runtime_error("MMC_MNA: invalid hot-start feedforward Vd");
    mFilteredVd = v.vd;
    mFilteredVdDerivative = 0;
    mHeldId = (2.0 / 3.0) * mPRef / v.vd;
  }

  // Only the static operating point uses a Newton solve. Its unknowns are
  // branch currents, capacitor voltages and active PI integrator memories.
  // This Jacobian is discarded before time stepping; it is not a dynamical model.
  std::vector<UInt> integrators;
  if ((mActiveMode == ActiveControlMode::ActivePower ||
       mActiveMode == ActiveControlMode::DcVoltage) &&
      mActive.ki > 0)
    integrators.push_back(Active);
  if (mReactiveMode != ReactiveControlMode::OpenLoop && mReactive.ki > 0)
    integrators.push_back(Reactive);
  if (mOutput.ki > 0) {
    integrators.push_back(OccD);
    integrators.push_back(OccQ);
  }
  if (mCirculating.ki > 0) {
    integrators.push_back(CccD);
    integrators.push_back(CccQ);
  }
  if (mZero.ki > 0)
    integrators.push_back(Zcc);
  if (mEnergyEnabled && mEnergy.ki > 0)
    integrators.push_back(Energy);
  const UInt size = BranchCount + integrators.size();
  auto equilibriumMemory = [&](const Electrical &e, Memory memory) {
    const Real pac = 1.5 * (v.vd * e[IDeltaD] + v.vq * e[IDeltaQ]);
    const Real qac = 1.5 * (v.vq * e[IDeltaD] - v.vd * e[IDeltaQ]);
    memory.filter = {v.vd,
                     v.vq,
                     pac,
                     pac,
                     qac,
                     qac,
                     v.vdc,
                     v.vdc,
                     1.5 * std::hypot(v.vd, v.vq),
                     1.5 * std::hypot(v.vd, v.vq)};
    // At equilibrium the delay has unity gain and zero derivative.
    if (mDelay > 0) {
      const auto command = control(e, memory, v);
      for (UInt i = 0; i < 5; ++i) {
        memory.delay[2 * i] = command.modulation[i] * mDelay * mDelay / 12;
        memory.delay[2 * i + 1] = 0;
      }
    }
    return memory;
  };
  auto unpack = [&](const Vector &values, Electrical &e, Memory &memory) {
    e = mElectrical;
    memory = mMemory;
    for (UInt i = 0; i < BranchCount; ++i)
      e[i] = values(i);
    for (UInt i = 0; i < integrators.size(); ++i)
      memory.integral[integrators[i]] = values(BranchCount + i);
    memory = equilibriumMemory(e, memory);
  };
  auto residual = [&](const Vector &values) {
    Electrical e;
    Memory memory;
    unpack(values, e, memory);
    const auto command = control(e, memory, v);
    const auto balance = electricalBalance(e, command.modulation, v);
    Vector result(size);
    for (UInt i = 0; i < BranchCount; ++i)
      result(i) = balance[i] / mStorage[i];
    for (UInt i = 0; i < integrators.size(); ++i) {
      const UInt k = integrators[i];
      result(BranchCount + i) =
          mLoadedStart && mActiveMode == ActiveControlMode::DcVoltage &&
                  k == Active
              ? 1.5 * (v.vd * e[IDeltaD] + v.vq * e[IDeltaQ]) - mInitialP
              : command.rate.integral[k];
    }
    return result;
  };
  auto scaling = [&](const Vector &values) {
    Vector scale(size);
    for (UInt i = 0; i < size; ++i)
      scale(i) = std::max(1.0, mOmega * std::max(1.0, std::abs(values(i))));
    for (UInt i = 0; i < integrators.size(); ++i)
      if (mLoadedStart && mActiveMode == ActiveControlMode::DcVoltage &&
          integrators[i] == Active)
        scale(BranchCount + i) = std::max(1.0, std::abs(mInitialP));
    return scale;
  };
  Vector values(size);
  for (UInt i = 0; i < BranchCount; ++i)
    values(i) = mElectrical[i];
  for (UInt i = 0; i < integrators.size(); ++i)
    values(BranchCount + i) = mMemory.integral[integrators[i]];
  Real norm = 0;
  for (UInt iteration = 0;; ++iteration) {
    const Vector r = residual(values), scale = scaling(values);
    norm = r.cwiseQuotient(scale).norm();
    if (!r.allFinite())
      throw std::runtime_error("MMC_MNA: non-finite initial residual");
    if (!mSolveStart || norm <= mInitialTolerance)
      break;
    if (iteration >= mInitialIterations)
      throw std::runtime_error(
          "MMC_MNA: operating-point solve did not converge");
    Matrix jacobian(size, size);
    Vector columnScale(size);
    for (UInt col = 0; col < size; ++col) {
      columnScale(col) = std::max(1.0, std::abs(values(col)));
      const Real h = 1e-8 + 1e-6 * columnScale(col);
      Vector plus = values, minus = values;
      plus(col) += h;
      minus(col) -= h;
      jacobian.col(col) =
          ((residual(plus) - residual(minus)) / (2 * h)).cwiseQuotient(scale) *
          columnScale(col);
    }
    const Vector scaledCorrection =
        jacobian.colPivHouseholderQr().solve(-r.cwiseQuotient(scale));
    const Vector correction = scaledCorrection.cwiseProduct(columnScale);
    if (!correction.allFinite())
      throw std::runtime_error("MMC_MNA: invalid initial correction");
    Bool accepted = false;
    for (Real fraction = 1; fraction >= 1.0 / 2048; fraction *= 0.5) {
      const Vector next = values + fraction * correction;
      if (residual(next).cwiseQuotient(scaling(next)).norm() < norm) {
        values = next;
        accepted = true;
        break;
      }
    }
    if (!accepted)
      throw std::runtime_error("MMC_MNA: initial line search failed");
  }
  unpack(values, mElectrical, mMemory);
  logScalar("equilibrium_residual_norm", norm);
}

void MMC_MNA::mnaCompInitialize(Real, Real step, Attribute<Matrix>::Ptr) {
  positive(step, "simulation step");
  mStep = step;
  // Explicit filter predictor must resolve the fastest configured pole.
  // The P2P example only uses the unconditionally stable sampled filter.
  for (Real tau : mFilterTime)
    if (tau > 0 && step >= tau)
      throw std::invalid_argument("MMC_MNA: use an EMT step smaller than each "
                                  "measurement-filter time constant");
  if (mDelay > 0 && step >= mDelay / 6)
    throw std::invalid_argument(
        "MMC_MNA: use an EMT step smaller than modulation delay / 6");
  validateTerminals();
  updateMatrixNodeIndices();
  mTerminalIndices.fill(-1);
  for (UInt i = 0; i < 3; ++i)
    mTerminalIndices[i] = matrixNodeIndex(0, i);
  if (terminalNotGrounded(1))
    mTerminalIndices[3] = matrixNodeIndex(1, 0);
  if (terminalNotGrounded(2))
    mTerminalIndices[4] = matrixNodeIndex(2, 0);
  for (UInt i = 0; i < BranchCount; ++i)
    mBranchIndices[i] = mVirtualNodes[i]->matrixNodeIndex(PhaseType::DC);
  const auto command = control(mElectrical, mMemory, mMeasurement);
  const auto balance =
      electricalBalance(mElectrical, command.modulation, mMeasurement);
  for (UInt i = 0; i < BranchCount; ++i)
    mHistory[i] = mStorage[i] / (mTheta * mStep) * mElectrical[i] +
                  (1 - mTheta) / mTheta * balance[i];
}

void MMC_MNA::mnaCompApplySystemMatrixStamp(SparseMatrixRow &matrix) {
  auto stamp = [&](UInt row, UInt col, Real value) {
    Math::addToMatrixElement(matrix, mBranchIndices[row], mBranchIndices[col],
                             value);
  };
  // Generalized-theta companion coefficients: L/(theta*h) and C/(theta*h).
  for (UInt i = 0; i < BranchCount; ++i)
    stamp(i, i, mStorage[i] / (mTheta * mStep));
  electricalCouplings(mStepModulation, [&](UInt row, UInt col, Real value) {
    stamp(row, col, -value);
  });
  for (UInt row = 0; row < BranchCount; ++row)
    for (UInt col = 0; col < BranchCount; ++col)
      stamp(row, col, -mSourceConductance[row][col]);
  for (UInt phase = 0; phase < 3; ++phase) {
    const Real angle = mNextAngle - phase * 2 * PI / 3;
    const Real c = std::cos(angle), s = std::sin(angle);
    const UInt nodeIndex = mTerminalIndices[phase];
    // Terminal KCL: component current = -abc(iDelta).
    Math::addToMatrixElement(matrix, nodeIndex, mBranchIndices[IDeltaD], -c);
    Math::addToMatrixElement(matrix, nodeIndex, mBranchIndices[IDeltaQ], s);
    // Branch KVL: +v_grid_d/q on the left hand side.
    Math::addToMatrixElement(matrix, mBranchIndices[IDeltaD], nodeIndex,
                             (2.0 / 3.0) * c);
    Math::addToMatrixElement(matrix, mBranchIndices[IDeltaQ], nodeIndex,
                             -(2.0 / 3.0) * s);
    for (UInt row = 0; row < BranchCount; ++row)
      Math::addToMatrixElement(
          matrix, mBranchIndices[row], nodeIndex,
          -(2.0 / 3.0) * (c * mSourceConductance[row][BranchCount] -
                          s * mSourceConductance[row][BranchCount + 1]));
  }
  for (UInt pole = 0; pole < 2; ++pole) {
    const Int nodeIndex = mTerminalIndices[3 + pole];
    if (nodeIndex < 0)
      continue;
    const Real sign = pole == 0 ? 1.0 : -1.0;
    Math::addToMatrixElement(matrix, nodeIndex, mBranchIndices[ISigmaZ],
                             3 * sign);
    Math::addToMatrixElement(matrix, mBranchIndices[ISigmaZ], nodeIndex,
                             -0.5 * sign);
    for (UInt row = 0; row < BranchCount; ++row)
      Math::addToMatrixElement(matrix, mBranchIndices[row], nodeIndex,
                               -sign *
                                   mSourceConductance[row][BranchCount + 2]);
  }
}
void MMC_MNA::mnaCompApplyRightSideVectorStamp(Matrix &rhs) {
  for (UInt i = 0; i < BranchCount; ++i)
    Math::setVectorElement(rhs, mBranchIndices[i], mHistory[i]);
}
void MMC_MNA::mnaCompPreStep(Real, Int) {
  sampleFeedforward();
  const auto old = control(mElectrical, mMemory, mMeasurement);
  mOldRate = old.rate;
  const auto balance =
      electricalBalance(mElectrical, old.modulation, mMeasurement);
  const Electrical &operatingPoint = mElectrical;
  for (UInt i = 0; i < BranchCount; ++i) {
    mHistory[i] = mStorage[i] / (mTheta * mStep) * mElectrical[i] +
                  (1 - mTheta) / mTheta * balance[i];
  }
  mPredictedMemory = advance(mMemory, mOldRate, mStep);
  // Linearize the circuit sources at the known electrical operating point.
  // An explicit electrical extrapolation can select the opposite side of a
  // modulation limiter before the network has actually reached that limit.
  // The PI memories are predicted, while electrical feedback is implicit.
  mStepModulation =
      control(operatingPoint, mPredictedMemory, mMeasurement).modulation;
  // Direct incremental conductances of the controlled arm sources. Holding
  // terminal-voltage feedforward explicitly would destabilize the inductive
  // AC/DC network. These contributions stay in the global MNA circuit solve.
  std::array<Electrical, 5> insertionCoupling;
  const auto constant = electricalBalance(operatingPoint, {}, {0, 0, 0});
  for (UInt j = 0; j < 5; ++j) {
    Modulation unit{};
    unit[j] = 1;
    insertionCoupling[j] = electricalBalance(operatingPoint, unit, {0, 0, 0});
    for (UInt row = 0; row < BranchCount; ++row)
      insertionCoupling[j][row] -= constant[row];
  }
  for (UInt col = 0; col < BranchCount + 3; ++col) {
    Electrical ep = operatingPoint, em = operatingPoint;
    Measurement vp = mMeasurement, vm = mMeasurement;
    Real value = col < BranchCount        ? operatingPoint[col]
                 : col == BranchCount     ? mMeasurement.vd
                 : col == BranchCount + 1 ? mMeasurement.vq
                                          : mMeasurement.vdc;
    const Real increment = 1e-6 * std::max(1.0, std::abs(value));
    if (col < BranchCount) {
      ep[col] += increment;
      em[col] -= increment;
    } else if (col == BranchCount) {
      vp.vd += increment;
      vm.vd -= increment;
    } else if (col == BranchCount + 1) {
      vp.vq += increment;
      vm.vq -= increment;
    } else {
      vp.vdc += increment;
      vm.vdc -= increment;
    }
    const auto plus = control(ep, mPredictedMemory, vp).modulation;
    const auto minus = control(em, mPredictedMemory, vm).modulation;
    for (UInt row = 0; row < BranchCount; ++row) {
      Real conductance = 0;
      for (UInt j = 0; j < 5; ++j)
        conductance +=
            insertionCoupling[j][row] * (plus[j] - minus[j]) / (2 * increment);
      mSourceConductance[row][col] = conductance;
      mHistory[row] -= conductance * value;
    }
  }
  mNextAngle = mInitialAngle + (mStepCount + 1) * mStep * mOmega;
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}
void MMC_MNA::mnaCompPostStep(Real, Int, Attribute<Matrix>::Ptr &leftVector) {
  for (UInt i = 0; i < 5; ++i)
    (**mIntfVoltage)(i) =
        mTerminalIndices[i] < 0 ? 0.0 : (**leftVector)(mTerminalIndices[i]);
  for (UInt i = 0; i < BranchCount; ++i) {
    mElectrical[i] = (**leftVector)(mBranchIndices[i]);
    if (!std::isfinite(mElectrical[i]))
      throw std::runtime_error("MMC_MNA: non-finite MNA solution");
  }
  ++mStepCount;
  mAngle = mNextAngle;
  mMeasurement = measure(**mIntfVoltage, mAngle);
  const auto corrected = control(mElectrical, mPredictedMemory, mMeasurement);
  mMemory = advance(advance(mMemory, mOldRate, (1 - mTheta) * mStep),
                    corrected.rate, mTheta * mStep);
  updateAttributes();
}
void MMC_MNA::mnaCompAddPreStepDependencies(AttributeBase::List &previous,
                                            AttributeBase::List &,
                                            AttributeBase::List &modified) {
  previous.push_back(mElectricalAttribute);
  previous.push_back(mIntfVoltage);
  modified.push_back(mRightVector);
}
void MMC_MNA::mnaCompAddPostStepDependencies(
    AttributeBase::List &, AttributeBase::List &dependencies,
    AttributeBase::List &modified, Attribute<Matrix>::Ptr &leftVector) {
  dependencies.push_back(leftVector);
  modified.push_back(mElectricalAttribute);
  modified.push_back(mIntfVoltage);
  modified.push_back(mIntfCurrent);
  modified.push_back(mAcVoltageAttribute);
  modified.push_back(mAcCurrentAttribute);
  modified.push_back(mModulationAttribute);
  modified.push_back(mDifferentialAttribute);
  modified.push_back(mCommonAttribute);
  for (const auto &entry : mScalars)
    modified.push_back(entry.second);
}
void MMC_MNA::logScalar(const String &name, Real value) {
  **mScalars.at(name) = value;
}
Attribute<Real>::Ptr MMC_MNA::scalar(const String &name) const {
  return mScalars.at(name);
}
Attribute<Matrix>::Ptr MMC_MNA::acTerminalVoltageAttribute() const {
  return mAcVoltageAttribute;
}
Attribute<Matrix>::Ptr MMC_MNA::acTerminalCurrentAttribute() const {
  return mAcCurrentAttribute;
}
Attribute<Real>::Ptr MMC_MNA::dcVoltageAttribute() const {
  return scalar("vdc");
}
Attribute<Real>::Ptr MMC_MNA::dcCurrentAttribute() const {
  return scalar("idc");
}
Attribute<Real>::Ptr MMC_MNA::activePowerAttribute() const {
  return scalar("p_ac");
}
Attribute<Real>::Ptr MMC_MNA::reactivePowerAttribute() const {
  return scalar("q_ac");
}
Attribute<Real>::Ptr MMC_MNA::storedEnergyAttribute() const {
  return scalar("stored_energy");
}

void MMC_MNA::updateAttributes() {
  const auto &e = mElectrical;
  const auto &v = mMeasurement;
  const auto command = control(e, mMemory, v);
  for (UInt i = 0; i < BranchCount; ++i)
    (**mElectricalAttribute)(i) = e[i];
  for (UInt phase = 0; phase < 3; ++phase) {
    const Real angle = mAngle - phase * 2 * PI / 3;
    (**mIntfCurrent)(phase) =
        -std::cos(angle) * e[IDeltaD] + std::sin(angle) * e[IDeltaQ];
  }
  (**mIntfCurrent)(3) = 3 * e[ISigmaZ];
  (**mIntfCurrent)(4) = -3 * e[ISigmaZ];
  **mAcVoltageAttribute = (**mIntfVoltage).topRows(3);
  **mAcCurrentAttribute = (**mIntfCurrent).topRows(3);
  // Report the controller command at the solved operating point, as the
  // reference example does. The circuit used its incremental source stamps.
  for (UInt i = 0; i < 5; ++i)
    (**mModulationAttribute)(i) = command.modulation[i];
  for (UInt i = 0; i < 2; ++i)
    (**mDifferentialAttribute)(i) = command.differential[i];
  for (UInt i = 0; i < 3; ++i)
    (**mCommonAttribute)(i) = command.common[i];
  const Real p = 1.5 * (v.vd * e[IDeltaD] + v.vq * e[IDeltaQ]);
  const Real q = 1.5 * (v.vq * e[IDeltaD] - v.vd * e[IDeltaQ]);
  logScalar("p_ac", p);
  logScalar("q_ac", q);
  logScalar("vdc", v.vdc);
  logScalar("vdcp", (**mIntfVoltage)(3));
  logScalar("vdcn", (**mIntfVoltage)(4));
  logScalar("idc", 3 * e[ISigmaZ]);
  logScalar("p_dc", v.vdc * 3 * e[ISigmaZ]);
  logScalar("stored_energy", energy(e));
  logScalar("v_grid_d", v.vd);
  logScalar("v_grid_q", v.vq);
  logScalar("i_delta_d", e[IDeltaD]);
  logScalar("i_delta_q", e[IDeltaQ]);
  logScalar("i_sigma_z", e[ISigmaZ]);
  logScalar("i_sigma_d", e[ISigmaD]);
  logScalar("i_sigma_q", e[ISigmaQ]);
  logScalar("i_delta_d_ref", command.idRef);
  logScalar("i_delta_q_ref", command.iqRef);
  logScalar("i_sigma_z_ref", command.izRef);
  logScalar("p_filtered", command.pFiltered);
  logScalar("q_filtered", command.qFiltered);
  logScalar("vdc_filtered", command.vdcFiltered);
  logScalar("v_d_feedforward_filtered", mFilteredVd);
  logScalar("i_delta_d_feedforward_held", mHeldId);
  logScalar("pll_frequency", (mOmega + command.deltaOmega) / (2 * PI));
  logScalar("pll_angle_deviation", std::remainder(mMemory.pllAngle, 2 * PI));
  logScalar("pll_error", command.pllError);
  logScalar("grid_angle", std::remainder(mAngle, 2 * PI));
  const Real loss =
      1.5 * mAcR * (e[IDeltaD] * e[IDeltaD] + e[IDeltaQ] * e[IDeltaQ]) +
      3 * mArmR * (e[ISigmaD] * e[ISigmaD] + e[ISigmaQ] * e[ISigmaQ]) +
      6 * mArmR * e[ISigmaZ] * e[ISigmaZ];
  logScalar("power_balance_error", p - v.vdc * 3 * e[ISigmaZ] + loss);
}
