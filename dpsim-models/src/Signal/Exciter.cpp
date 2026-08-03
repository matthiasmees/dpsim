/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/MathUtils.h>
#include <dpsim-models/Signal/Exciter.h>

#include <cmath>

using namespace CPS;
using namespace CPS::Signal;

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
  if (!(Ta > 0.0) || !(Te > 0.0) || !(Tf > 0.0) || !(Tr > 0.0) ||
      !(Ka != 0.0) || !(maxVr > minVr) || !std::isfinite(Ta) ||
      !std::isfinite(Ka) || !std::isfinite(Te) || !std::isfinite(Ke) ||
      !std::isfinite(Tf) || !std::isfinite(Kf) || !std::isfinite(Tr) ||
      !std::isfinite(maxVr) || !std::isfinite(minVr)) {
    SPDLOG_LOGGER_ERROR(
        mSLog,
        "Invalid Exciter parameters: require Ta, Te, Tf, Tr > 0, Ka != 0, "
        "and maxVr > minVr");
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
  // Exact quadratic saturation characteristic from the supplied legacy class.
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
  if (!std::isfinite(VhInit) || !std::isfinite(EfInit)) {
    SPDLOG_LOGGER_ERROR(mSLog, "Exciter initial values must be finite");
    throw CPS::InvalidArgumentException();
  }

  SPDLOG_LOGGER_INFO(mSLog,
                     "Initially set excitation-system values:"
                     "\nVh_init: {:e}"
                     "\nEf_init: {:e}",
                     VhInit, EfInit);

  **mVm = VhInit;
  **mVh = VhInit;
  **mEf = EfInit;

  **mVse = saturation(**mEf);

  // Vis is vr2 in the original PSAT-oriented implementation.
  **mVis = **mEf;

  // Vr is vr1 in the original implementation.
  **mVr = mKe * **mEf + **mVse;
  if (**mVr > mMaxVr)
    **mVr = mMaxVr;
  else if (**mVr < mMinVr)
    **mVr = mMinVr;

  mVref = **mVr / mKa + **mVm;

  // Keep explicit previous-value storage consistent before the first step.
  mVmPrev = **mVm;
  mVisPrev = **mVis;
  mVrPrev = **mVr;
  mEfPrev = **mEf;

  SPDLOG_LOGGER_INFO(mSLog,
                     "Applied excitation-system initial values:"
                     "\nVref: {:e}"
                     "\nVm: {:e}"
                     "\nEf: {:e}"
                     "\nVse: {:e}"
                     "\nVr: {:e}"
                     "\nVis: {:e}",
                     mVref, **mVm, **mEf, **mVse, **mVr, **mVis);
}

Real Exciter::step(Real Vd, Real Vq, Real dt, Real Vpss) {
  if (!(dt > 0.0) || !std::isfinite(dt)) {
    SPDLOG_LOGGER_ERROR(mSLog, "Exciter time step must be finite and positive");
    throw CPS::InvalidArgumentException();
  }

  // Voltage magnitude calculation.
  **mVh = std::hypot(Vd, Vq);

  // States at k-1.
  mVmPrev = **mVm;
  mVisPrev = **mVis;
  mVrPrev = **mVr;
  mEfPrev = **mEf;

  // Voltage transducer equation. This intentionally retains the Euler update
  // from the supplied legacy model.
  **mVm = Math::StateSpaceEuler(mVmPrev, -1.0 / mTr, 1.0 / mTr, dt, **mVh);

  // Quadratic saturation based on the previous field voltage.
  **mVse = saturation(mEfPrev);

  // Stabilizing-feedback state, retained from the legacy implementation.
  **mVis = Math::StateSpaceEuler(mVisPrev, -1.0 / mTf, 1.0 / mTf, dt, mEfPrev);
  const Real stabilizingFeedback =
      (-mKf / mTf) * **mVis + (mKf / mTf) * mEfPrev;

  // Voltage-regulator equation. Vpss is zero in the supplied PSHA scenario;
  // adding it here only connects the current Base::Exciter interface without
  // changing the no-PSS equations.
  **mVr = Math::StateSpaceEuler(mVrPrev, -1.0 / mTa, mKa / mTa, dt,
                                mVref + Vpss - **mVm - stabilizingFeedback);

  if (**mVr > mMaxVr)
    **mVr = mMaxVr;
  else if (**mVr < mMinVr)
    **mVr = mMinVr;

  // Exciter equation, retained from the supplied class.
  **mEf =
      Math::StateSpaceEuler(mEfPrev, -mKe / mTe, 1.0 / mTe, dt, **mVr - **mVse);

  return **mEf;
}
