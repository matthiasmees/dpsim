/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/MathUtils.h>
#include <dpsim-models/Signal/TurbineGovernor.h>

using namespace CPS;
using namespace CPS::Signal;

void TurbineGovernor::setParameters(Real Ta, Real Tb, Real Tc, Real Fa, Real Fb,
                                    Real Fc, Real K, Real Tsr, Real Tsm) {
  mTa = Ta;
  mTb = Tb;
  mTc = Tc;
  mFa = Fa;
  mFb = Fb;
  mFc = Fc;

  mK = K;
  mTsr = Tsr;
  mTsm = Tsm;

  SPDLOG_LOGGER_INFO(mSLog,
                     "Turbine parameters: \n"
                     "Ta: {:e}\nTb: {:e}\nTc: {:e}\n"
                     "Fa: {:e}\nFb: {:e}\nFc: {:e}\n",
                     mTa, mTb, mTc, mFa, mFb, mFc);

  SPDLOG_LOGGER_INFO(mSLog,
                     "Governor parameters: \n"
                     "K: {:e}\nTsr: {:e}\nTsm: {:e}\n",
                     mK, mTsr, mTsm);
}

void TurbineGovernor::initialize(Real PmRef, Real Tm_init) {
  mTm = Tm_init;
  //T1 = (1 - mFa) * PmRef;
  //T1 = (mTm - PmRef*mFa)/(mFb + mFc);
  //T2 = mFa * PmRef;
  T1 = PmRef;
  T2 = PmRef;
  T3 = PmRef;

  SPDLOG_LOGGER_INFO(mSLog,
                     "Turbine initial values: \n"
                     "init_Tm: {:e}\ninit_T1: {:e}\ninit_T2: {:e}\n",
                     mTm, T1, T2);

  mVcv = PmRef;
  mpVcv = 0;
  Psm_in = PmRef;

  SPDLOG_LOGGER_INFO(mSLog,
                     "Governor initial values: \n"
                     "init_Vcv: {:e}\ninit_pVcv: {:e}\ninit_Psm_in: {:e}\n",
                     mVcv, mpVcv, Psm_in);
}

Real TurbineGovernor::step(Real Om, Real OmRef, Real PmRef, Real dt) {
  // ### Governing ###
  // Modelled according to V. Sapucaia-Gunzenhauser, "Modeling and Simulation of Rotating Machines", p.45

  // Input of speed relay
  Psr_in = PmRef + (OmRef - Om) * mK;

  // Input of servor motor
  Psm_in = Math::StateSpaceEuler(Psm_in, -1 / mTsr, 1 / mTsr, dt,
                                 Psr_in); // ovo je bilo
  //Psm_in = Math::StateSpaceTrapezoidal(Psm_in, -1 / mTsr, 1 / mTsr, dt, Psr_in); // oov sam ja dodao
  // rate of change of valve
  // including limiter with upper bound 0.1pu/s and lower bound -1.0pu/s
  mpVcv = (Psm_in - mVcv) / mTsm;
  if (mpVcv >= 0.3)
    mpVcv = 0.3;
  else if (mpVcv <= -0.3)
    mpVcv = -0.3;

  // valve position
  // including limiter with upper bound 1 and lower bound 0
  mVcv = mVcv + dt * mpVcv;
  if (mVcv >= 1)
    mVcv = 1;
  else if (mVcv <= 0)
    mVcv = 0;

  // ### Turbine ###
  // Simplified Single Reheat Tandem-Compound Steam Turbine
  // Modelled according to V. Sapucaia-Gunzenhauser, "Modeling and Simulation of Rotating Machines", p.44

  // T1 = Math::StateSpaceEuler(T1, -1 / mTb, 1 / mTb, dt, mVcv); // ovo je bilo
  Real F1 = 0.3;
  T1 = Math::StateSpaceTrapezoidal(T1, -1 / mTa, 1 / mTa, dt,
                                   mVcv); // ovo sam ja dodao
  T1a = mFa * T1 + F1 * mVcv;

  T2 = Math::StateSpaceTrapezoidal(T2, -1 / mTb, 1 / mTb, dt, T1);
  T2b = mFb * T2;

  T3 = Math::StateSpaceTrapezoidal(T3, -1 / mTc, 1 / mTc, dt, T2);
  T3c = mFc * T3;

  mTm = T1a + T2b + T3c;
  //T2 = mVcv * mFa;
  //mTm = (mFb + mFc) * T1 + T2;

  return mTm;
}
