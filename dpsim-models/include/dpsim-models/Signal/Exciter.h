/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <memory>

#include <dpsim-models/Base/Base_Exciter.h>
#include <dpsim-models/Logger.h>
#include <dpsim-models/SimSignalComp.h>

namespace CPS {
namespace Signal {

/// Parameter object for the legacy DPsim excitation-system equations that were
/// previously configured through the nine-argument SynchronGenerator::addExciter
/// overload.
class ExciterParameters final : public Base::ExciterParameters,
                                public SharedFactory<ExciterParameters> {
public:
  Real Ta = 0.0;
  Real Ka = 0.0;
  Real Te = 0.0;
  Real Ke = 0.0;
  Real Tf = 0.0;
  Real Kf = 0.0;
  Real Tr = 0.0;
  Real maxVr = 0.0;
  Real minVr = 0.0;
};

/// Legacy excitation-system model used by the original PSHA scenario.
///
/// The electrical/controller equations are intentionally retained. The class
/// now implements Base::Exciter so that it can be attached through the current
/// SynchronGenerator API.
class Exciter final : public Base::Exciter,
                      public SimSignalComp,
                      public SharedFactory<Exciter> {
public:
  Exciter(const String &name, Logger::Level logLevel = Logger::Level::info);

  /// Current polymorphic API used by Base::SynchronGenerator.
  void setParameters(std::shared_ptr<Base::ExciterParameters> parameters) final;
  void initializeStates(Real VhInit, Real EfInit) final;
  Real step(Real Vd, Real Vq, Real dt, Real Vpss = 0.0) final;

  /// Legacy convenience API retained for existing callers.
  void setParameters(Real Ta, Real Ka, Real Te, Real Ke, Real Tf, Real Kf,
                     Real Tr, Real maxVr, Real minVr);
  void initialize(Real VhInit, Real EfInit);

private:
  Real saturation(Real fieldVoltage) const;

  std::shared_ptr<Signal::ExciterParameters> mParameters;

  Real mTa = 0.0;
  Real mKa = 0.0;
  Real mTe = 0.0;
  Real mKe = 0.0;
  Real mTf = 0.0;
  Real mKf = 0.0;
  Real mTr = 0.0;
  Real mMaxVr = 0.0;
  Real mMinVr = 0.0;

  Real mVref = 0.0;

  Real mVmPrev = 0.0;
  Real mVisPrev = 0.0;
  Real mVrPrev = 0.0;
  Real mEfPrev = 0.0;

  const Attribute<Real>::Ptr mVm;
  const Attribute<Real>::Ptr mVh;
  const Attribute<Real>::Ptr mVis;
  const Attribute<Real>::Ptr mVse;
  const Attribute<Real>::Ptr mVr;
  const Attribute<Real>::Ptr mEf;
};

} // namespace Signal
} // namespace CPS
