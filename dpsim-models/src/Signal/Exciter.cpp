/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/Signal/Exciter.h>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace CPS;
using namespace CPS::Signal;

namespace {

Real exactFirstOrderStep(Real previousState, Real timeConstant,
                         Real steadyStateInput, Real timeStep) {
  const Real alpha = std::exp(-timeStep / timeConstant);
  return alpha * previousState + (1.0 - alpha) * steadyStateInput;
}

Real exactFieldStep(Real previousFieldVoltage, Real Ke, Real Te,
                    Real regulatorMinusSaturation, Real timeStep) {
  // dEf/dt = (-Ke*Ef + regulatorMinusSaturation) / Te.
  if (std::abs(Ke) > std::numeric_limits<Real>::epsilon()) {
    const Real alpha = std::exp(-Ke * timeStep / Te);
    const Real steadyStateFieldVoltage = regulatorMinusSaturation / Ke;
    return alpha * previousFieldVoltage +
           (1.0 - alpha) * steadyStateFieldVoltage;
  }

  // Ke == 0 degenerates to a pure integrator.
  return previousFieldVoltage + timeStep * regulatorMinusSaturation / Te;
}

} // namespace

Exciter::Exciter(const String &name, CPS::Logger::Level logLevel)
    : SimSignalComp(name, name, logLevel),
      mVm(mAttributes->create<Real>("Vm", 0.0)),
      mVh(mAttributes->create<Real>("Vh", 0.0)),
      mVis(mAttributes->create<Real>("Vis", 0.0)),
      mVse(mAttributes->create<Real>("Vse", 0.0)),
      mVr(mAttributes->create<Real>("Vr", 0.0)),
      mEf(mAttributes->create<Real>("Ef", 0.0)) {}

void Exciter::setParameters(
    std::shared_ptr<Base::ExciterParameters> parameters) {
  auto typedParameters =
      std::dynamic_pointer_cast<Signal::ExciterParameters>(parameters);

  if (!typedParameters) {
    SPDLOG_LOGGER_ERROR(
        mSLog, "Exciter requires Signal::ExciterParameters, but a different "
               "parameter type was supplied");
    throw CPS::TypeException();
  }

  mParameters = typedParameters;
  setParameters(mParameters->Ta, mParameters->Ka, mParameters->Te,
                mParameters->Ke, mParameters->Tf, mParameters->Kf,
                mParameters->Tr, mParameters->maxVr, mParameters->minVr);
}

void Exciter::setParameters(Real Ta, Real Ka, Real Te, Real Ke, Real Tf,
                            Real Kf, Real Tr, Real maxVr, Real minVr) {
  if (!(Ta > 0.0) || !(Ka > 0.0) || !(Te > 0.0) || !(Ke >= 0.0) ||
      !(Tf > 0.0) || !(Kf >= 0.0) || !(Tr > 0.0) || !(maxVr > minVr) ||
      !std::isfinite(Ta) || !std::isfinite(Ka) || !std::isfinite(Te) ||
      !std::isfinite(Ke) || !std::isfinite(Tf) || !std::isfinite(Kf) ||
      !std::isfinite(Tr) || !std::isfinite(maxVr) || !std::isfinite(minVr)) {
    SPDLOG_LOGGER_ERROR(
        mSLog, "Invalid Exciter parameters: require Ta, Ka, Te, Tf, Tr > 0, "
               "Ke, Kf >= 0, and maxVr > minVr");
    throw CPS::InvalidArgumentException();
  }

  mTa = Ta;
  mKa = Ka;
  mTe = Te;
  mKe = Ke;
  mTf = Tf;
  mKf = Kf;
  mTr = Tr;
  mMaxVr = maxVr;
  mMinVr = minVr;

  SPDLOG_LOGGER_INFO(mSLog,
                     "Exciter parameters:"
                     "\nTa: {:e}"
                     "\nKa: {:e}"
                     "\nTe: {:e}"
                     "\nKe: {:e}"
                     "\nTf: {:e}"
                     "\nKf: {:e}"
                     "\nTr: {:e}"
                     "\nMaximum regulator voltage: {:e}"
                     "\nMinimum regulator voltage: {:e}",
                     mTa, mKa, mTe, mKe, mTf, mKf, mTr, mMaxVr, mMinVr);
}

Real Exciter::saturation(Real fieldVoltage) const {
  constexpr Real E1 = 3.9;
  constexpr Real Se1 = 0.0001;
  constexpr Real E2 = 5.2;
  constexpr Real Se2 = 0.001;

  const Real sq = std::sqrt((E1 * Se1) / (E2 * Se2));
  const Real Asq = (E1 - E2 * sq) / (1.0 - sq);
  const Real Bsq = (E2 * Se2) / ((E2 - Asq) * (E2 - Asq));

  if (fieldVoltage > Asq) {
    const Real difference = fieldVoltage - Asq;
    return Bsq * difference * difference;
  }

  return 0.0;
}

void Exciter::initializeStates(Real VhInit, Real EfInit) {
  initialize(VhInit, EfInit);
}

void Exciter::initialize(Real VhInit, Real EfInit) {
  if (!(VhInit >= 0.0) || !std::isfinite(VhInit) || !std::isfinite(EfInit)) {
    SPDLOG_LOGGER_ERROR(
        mSLog, "Exciter initial values must be finite and Vh_init must be "
               "non-negative");
    throw CPS::InvalidArgumentException();
  }

  // Exact steady-state initialization of all four dynamic blocks:
  //   Vm_dot  = (-Vm + Vh) / Tr
  //   Vis_dot = (-Vis + Ef) / Tf
  //   Vr_dot  = (-Vr + Ka*(Vref - Vm - Vstab)) / Ta
  //   Ef_dot  = (-Ke*Ef + Vr - Vse(Ef)) / Te
  **mVh = VhInit;
  **mVm = VhInit;
  **mEf = EfInit;
  **mVis = EfInit;
  **mVse = saturation(EfInit);

  const Real requiredRegulatorVoltage = mKe * EfInit + **mVse;
  if (requiredRegulatorVoltage < mMinVr || requiredRegulatorVoltage > mMaxVr) {
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "Exciter cannot be initialized in steady state: required Vr={} is "
        "outside limiter [{}, {}]",
        requiredRegulatorVoltage, mMinVr, mMaxVr);
    throw CPS::InvalidArgumentException();
  }

  **mVr = requiredRegulatorVoltage;

  // Vis=Ef makes the washout/stabilizing-feedback output exactly zero.
  mVref = VhInit + requiredRegulatorVoltage / mKa;

  mVmPrev = **mVm;
  mVisPrev = **mVis;
  mVrPrev = **mVr;
  mEfPrev = **mEf;

  SPDLOG_LOGGER_INFO(mSLog,
                     "Exciter steady-state initialization:"
                     "\nVref: {:e}"
                     "\nVh: {:e}"
                     "\nVm: {:e}"
                     "\nEf: {:e}"
                     "\nVse: {:e}"
                     "\nVr: {:e}"
                     "\nVis: {:e}",
                     mVref, **mVh, **mVm, **mEf, **mVse, **mVr, **mVis);
}

Real Exciter::step(Real Vd, Real Vq, Real dt, Real Vpss) {
  if (!(dt > 0.0) || !std::isfinite(dt) || !std::isfinite(Vd) ||
      !std::isfinite(Vq) || !std::isfinite(Vpss)) {
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "Exciter inputs and time step must be finite and dt must be positive");
    throw CPS::InvalidArgumentException();
  }

  **mVh = std::hypot(Vd, Vq);

  mVmPrev = **mVm;
  mVisPrev = **mVis;
  mVrPrev = **mVr;
  mEfPrev = **mEf;

  // Exact zero-order-hold discretization of each first-order block. This
  // preserves a correctly initialized equilibrium exactly and removes the
  // artificial numerical kick caused by chained explicit-Euler updates.
  **mVm = exactFirstOrderStep(mVmPrev, mTr, **mVh, dt);

  **mVse = saturation(mEfPrev);

  **mVis = exactFirstOrderStep(mVisPrev, mTf, mEfPrev, dt);
  const Real stabilizingFeedback = (mKf / mTf) * (mEfPrev - **mVis);

  const Real regulatorTarget =
      mKa * (mVref + Vpss - **mVm - stabilizingFeedback);
  **mVr = exactFirstOrderStep(mVrPrev, mTa, regulatorTarget, dt);
  **mVr = std::clamp(**mVr, mMinVr, mMaxVr);

  **mEf = exactFieldStep(mEfPrev, mKe, mTe, **mVr - **mVse, dt);

  if (!std::isfinite(**mVm) || !std::isfinite(**mVis) ||
      !std::isfinite(**mVr) || !std::isfinite(**mEf)) {
    SPDLOG_LOGGER_ERROR(mSLog, "Exciter produced a non-finite internal state");
    throw CPS::InvalidArgumentException();
  }

  return **mEf;
}
