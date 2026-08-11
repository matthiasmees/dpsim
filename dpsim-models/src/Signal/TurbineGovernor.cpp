/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/Signal/TurbineGovernor.h>

#include <algorithm>
#include <cmath>

using namespace CPS;
using namespace CPS::Signal;

namespace {

  
constexpr Real MAXIMUM_VALVE_OPENING_RATE = 0.3;
constexpr Real MAXIMUM_VALVE_CLOSING_RATE = -0.3;

Real exactFirstOrderStep(Real previousState, Real timeConstant,
                         Real steadyStateInput, Real timeStep) {
  const Real alpha = std::exp(-timeStep / timeConstant);
  return alpha * previousState + (1.0 - alpha) * steadyStateInput;
}

Real turbineSteadyStateGain(Real F1a, Real Fa, Real Fb, Real Fc) {
  return F1a + Fa + Fb + Fc;
}

} // namespace

void TurbineGovernor::setParameters(Real Ta, Real Tb, Real Tc, Real F1a, Real Fa, Real Fb,
                                    Real Fc, Real K, Real Tsr, Real Tsm) {
  const Real steadyStateGain = turbineSteadyStateGain(F1a, Fa, Fb, Fc);

  if (!(Ta > 0.0) || !(Tb > 0.0) || !(Tc > 0.0) || !(Tsr > 0.0) ||
      !(Tsm > 0.0) || !(steadyStateGain > 0.0) || !std::isfinite(Ta) ||
      !std::isfinite(Tb) || !std::isfinite(Tc) || !std::isfinite(Fa) ||
      !std::isfinite(Fb) || !std::isfinite(Fc) || !std::isfinite(K) ||
      !std::isfinite(Tsr) || !std::isfinite(Tsm)) {
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "Invalid TurbineGovernor parameters: require Ta, Tb, Tc, Tsr, Tsm "
        "> 0, finite coefficients, and positive turbine steady-state gain");
    throw CPS::InvalidArgumentException();
  }

  mTa = Ta;
  mTb = Tb;
  mTc = Tc;
  mF1a = F1a;
  mFa = Fa;
  mFb = Fb;
  mFc = Fc;

  mK = K;
  mTsr = Tsr;
  mTsm = Tsm;

  SPDLOG_LOGGER_INFO(mSLog,
                     "Turbine parameters:"
                     "\nTa: {:e}"
                     "\nTb: {:e}"
                     "\nTc: {:e}"
                     "\nF1a: {:e}"
                     "\nFa: {:e}"
                     "\nFb: {:e}"
                     "\nFc: {:e}"
                     "\nsteady-state gain: {:e}",
                     mTa, mTb, mTc, mF1a, mFa, mFb, mFc,
                     steadyStateGain);

  SPDLOG_LOGGER_INFO(mSLog,
                     "Governor parameters:"
                     "\nK: {:e}"
                     "\nTsr: {:e}"
                     "\nTsm: {:e}",
                     mK, mTsr, mTsm);
}

void TurbineGovernor::initialize(Real PmRef, Real Tm_init) {
  if (!std::isfinite(PmRef) || !std::isfinite(Tm_init)) {
    SPDLOG_LOGGER_ERROR(mSLog, "TurbineGovernor initial values must be finite");
    throw CPS::InvalidArgumentException();
  }

  const Real steadyStateGain = turbineSteadyStateGain(mF1a, mFa, mFb, mFc);

  // PmRef and Tm_init are desired turbine-output torque/power in pu. The
  // internal governor signal is a valve position, so it must be divided by
  // the complete static turbine gain.
  const Real valveReference = PmRef / steadyStateGain;
  const Real valveInitial = Tm_init / steadyStateGain;

  if (valveInitial < 0.0 || valveInitial > 1.0 || valveReference < 0.0 ||
      valveReference > 1.0) {
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "TurbineGovernor operating point is outside valve limits: "
        "valve_init={}, valve_ref={}, limits=[0,1]",
        valveInitial, valveReference);
    throw CPS::InvalidArgumentException();
  }

  // Exact equilibrium of the complete governor and three-stage turbine.
  Psr_in = valveReference;
  Psm_in = valveInitial;
  mVcv = valveInitial;
  mpVcv = 0.0;

  T1 = valveInitial;
  T2 = valveInitial;
  T3 = valveInitial;

  T1a = mFa * T1 + mF1a * mVcv;
  T2b = mFb * T2;
  T3c = mFc * T3;
  mTm = T1a + T2b + T3c;

  SPDLOG_LOGGER_INFO(mSLog,
                     "TurbineGovernor steady-state initialization:"
                     "\nPmRef: {:e}"
                     "\nTm_init: {:e}"
                     "\nvalve_ref: {:e}"
                     "\nvalve_init: {:e}"
                     "\nT1: {:e}"
                     "\nT2: {:e}"
                     "\nT3: {:e}"
                     "\nTm: {:e}",
                     PmRef, Tm_init, valveReference, valveInitial, T1, T2, T3,
                     mTm);
}

Real TurbineGovernor::step(Real Om, Real OmRef, Real PmRef, Real dt) {
  if (!(dt > 0.0) || !std::isfinite(dt) || !std::isfinite(Om) ||
      !std::isfinite(OmRef) || !std::isfinite(PmRef)) {
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "TurbineGovernor inputs and time step must be finite and dt must be "
        "positive");
    throw CPS::InvalidArgumentException();
  }

  const Real steadyStateGain = turbineSteadyStateGain(mF1a, mFa, mFb, mFc);

  // Convert the desired mechanical-torque command to the valve domain. This
  // is essential because SynchronGeneratorVBR supplies PmRef in turbine-output
  // pu, while the turbine's internal static gain is generally not one.
  Psr_in = (PmRef + (OmRef - Om) * mK) / steadyStateGain;
  Psr_in = std::clamp(Psr_in, 0.0, 1.0);

  Psm_in = exactFirstOrderStep(Psm_in, mTsr, Psr_in, dt);

  mpVcv = (Psm_in - mVcv) / mTsm;
  mpVcv =
      std::clamp(mpVcv, MAXIMUM_VALVE_CLOSING_RATE, MAXIMUM_VALVE_OPENING_RATE);

  mVcv = std::clamp(mVcv + dt * mpVcv, 0.0, 1.0);

  // Exact zero-order-hold discretization of the three turbine lags. A correct
  // steady-state initialization remains exactly stationary at the first step.
  T1 = exactFirstOrderStep(T1, mTa, mVcv, dt);
  T2 = exactFirstOrderStep(T2, mTb, T1, dt);
  T3 = exactFirstOrderStep(T3, mTc, T2, dt);

  T1a = mFa * T1 + mF1a * mVcv;
  T2b = mFb * T2;
  T3c = mFc * T3;
  mTm = T1a + T2b + T3c;

  if (!std::isfinite(mTm)) {
    SPDLOG_LOGGER_ERROR(
        mSLog, "TurbineGovernor produced non-finite mechanical torque");
    throw CPS::InvalidArgumentException();
  }

  return mTm;
}
