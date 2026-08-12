/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/AttributeList.h>
#include <dpsim-models/Base/Base_Exciter.h>
#include <dpsim-models/Base/Base_Governor.h>
#include <dpsim-models/Base/Base_PSS.h>
#include <dpsim-models/Base/Base_Turbine.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/Signal/TurbineGovernor.h>
#include <dpsim-models/Signal/TurbineGovernorType1.h>

namespace CPS {
namespace Base {

/// @brief Base synchronous generator model
class SynchronGenerator {
public:
  /// \brief State type of machine variables.
  enum class StateType { perUnit, statorReferred, rotorReferred };

  /// \brief Machine parameters type.
  enum class ParameterType { perUnit, statorReferred, operational };

  /// Initialization mode for the legacy Signal::TurbineGovernor.
  ///
  /// FromPowerflow means that the governor parameters are attached when
  /// addGovernor() is called, while its equilibrium state is initialized later
  /// from the synchronous-generator operating point established by
  /// initializeFromNodesAndTerminals().
  enum class LegacyGovernorInitialization { FromPowerflow };

  // Legacy TurbineGovernor with explicit initial mechanical operating point.
  // Kept unchanged for backwards compatibility.
  void addGovernor(Real Ta, Real Tb, Real Tc, Real F1a, Real Fa, Real Fb,
                   Real Fc, Real K, Real Tsr, Real Tsm, Real Tm_init,
                   Real PmRef);

  // Legacy TurbineGovernor with automatic initial operating point derived from
  // the generator's PF-initialized mechanical power.
  void addGovernor(Real Ta, Real Tb, Real Tc, Real F1a, Real Fa, Real Fb,
                   Real Fc, Real K, Real Tsr, Real Tsm,
                   LegacyGovernorInitialization initialization);

  /// Add TurbineGovernorType1 (already constructed and initialised)
  void
  addGovernor(std::shared_ptr<Signal::TurbineGovernorType1> turbineGovernor);

  /// Deprecated scalar convenience for TurbineGovernorType1.
  void addGovernor(Real T3, Real T4, Real T5, Real Tc, Real Ts, Real R,
                   Real Tmin, Real Tmax, Real OmRef, Real TmRef);

  /// Add a modular governor + turbine pair
  void
  addGovernorAndTurbine(std::shared_ptr<Base::Governor> governor,
                        std::shared_ptr<Base::GovernorParameters> govParams,
                        std::shared_ptr<Base::Turbine> turbine,
                        std::shared_ptr<Base::TurbineParameters> turbineParams);

  /// Add a pre-constructed governor + turbine pair
  void addGovernorAndTurbine(std::shared_ptr<Base::Governor> governor,
                             std::shared_ptr<Base::Turbine> turbine);

  /// Add voltage regulator and exciter
  void addExciter(std::shared_ptr<Base::Exciter> exciter,
                  std::shared_ptr<Base::ExciterParameters> params);

  /// Add already constructed regulator and exciter
  void addExciter(std::shared_ptr<Base::Exciter> exciter);

  // Deprecated method
  void addExciter(Real Ta, Real Ka, Real Te, Real Ke, Real Tf, Real Kf,
                  Real Tr);

  /// Attach a PSS (initialised separately) to this generator
  void addPSS(std::shared_ptr<Base::PSS> pss,
              std::shared_ptr<Base::PSSParameters> parameters);

  /// Attach a pre-constructed PSS
  void addPSS(std::shared_ptr<Base::PSS> pss);

protected:
  NumericalMethod
      mNumericalMethod; // not needed if sundials used; could instead determine implicit / explicit solve

  /// Simulation angular system speed
  Real mSystemOmega;

  /// Simulation time step
  Real mTimeStep;

  /// specifies if the machine parameters are transformed to per unit
  StateType mStateType = StateType::perUnit;

  ParameterType mParameterType;

  /// Flag to remember when initial values are set
  Bool mInitialValuesSet = false;

  // ### Machine parameters ###
  Real mNomPower = 0;
  Real mNomVolt = 0;
  Real mNomFreq = 0;
  Real mNomOmega = 0;
  Real mNomFieldCur = 0;
  Int mNumDampingWindings = 0;
  Int mPoleNumber = 0;

  Real mLmd = 0;
  Real mLmq = 0;
  Real mRfd = 0;
  Real mLlfd = 0;
  Real mLfd = 0;
  Real mLf = 0;
  Real mRkd = 0;
  Real mLlkd = 0;
  Real mLkd = 0;
  Real mRkq1 = 0;
  Real mLlkq1 = 0;
  Real mLkq1 = 0;
  Real mRkq2 = 0;
  Real mLlkq2 = 0;
  Real mLkq2 = 0;

public:
  const Attribute<Real>::Ptr mRs;
  const Attribute<Real>::Ptr mLl;
  const Attribute<Real>::Ptr mLd;
  const Attribute<Real>::Ptr mLq;

  const Attribute<Real>::Ptr mLd_t;
  const Attribute<Real>::Ptr mLq_t;
  const Attribute<Real>::Ptr mLd_s;
  const Attribute<Real>::Ptr mLq_s;
  const Attribute<Real>::Ptr mTd0_t;
  const Attribute<Real>::Ptr mTq0_t;
  const Attribute<Real>::Ptr mTd0_s;
  const Attribute<Real>::Ptr mTq0_s;

protected:
  // #### Initial Values ####
  Complex mInitElecPower = 0;
  Complex mInitTermVoltage = 0;
  Real mInitMechPower = 0;

  // ### Stator base values ###
  Real mBase_V = 0;
  Real mBase_V_RMS = 0;
  Real mBase_I = 0;
  Real mBase_I_RMS = 0;
  Real mBase_Z = 0;
  Real mBase_OmElec = 0;
  Real mBase_OmMech = 0;
  Real mBase_L = 0;
  Real mBase_T = 0;
  Real mBase_Psi = 0;

  Real mBase_ifd = 0;
  Real mBase_vfd = 0;
  Real mBase_Zfd = 0;
  Real mBase_Lfd = 0;

  // ### Useful Matrices ###
  Matrix mInductanceMat;
  Matrix mResistanceMat;
  Matrix mInvInductanceMat;

  // ### State variables ###
  Real mThetaMech = 0;

public:
  const Attribute<Real>::Ptr mDelta;
  const Attribute<Real>::Ptr mMechTorque;
  const Attribute<Real>::Ptr mInertia;
  const Attribute<Real>::Ptr mOmMech;
  const Attribute<Real>::Ptr mElecActivePower;
  const Attribute<Real>::Ptr mElecReactivePower;
  const Attribute<Real>::Ptr mMechPower;
  const Attribute<Real>::Ptr mElecTorque;

protected:
  Matrix mVsr;
  Matrix mIsr;
  Matrix mPsisr;

  // #### dq-frame specific variables ####
  Matrix mVdq0;
  Matrix mIdq0;
  Matrix mFluxStateSpaceMat;
  Matrix mOmegaFluxMat;
  Matrix mFluxToCurrentMat;
  Real mLad;
  Real mLaq;
  Real mDetLd;
  Real mDetLq;

  Bool mCompensationOn;
  Real mRcomp;

  /// Initializes states in per unit.
  void initPerUnitStates();

  // #### Controllers ####
  Bool mHasTurbineGovernor = false;
  Bool mHasTurbineGovernorType1 = false;
  Bool mHasGovernorAndTurbine = false;
  Bool mHasExciter = false;
  Bool mHasPSS = false;

  /// True only for the new legacy-governor PF initialization mode.
  Bool mLegacyGovernorInitFromPowerflow = false;

  // Deprecated
  Real mInitTerminalVoltage = 0;
  Real mInitVoltAngle = 0;

  /// Constructor
  explicit SynchronGenerator(CPS::AttributeList::Ptr attributeList)
      : mRs(attributeList->create<Real>("Rs", 0)),
        mLl(attributeList->create<Real>("Ll", 0)),
        mLd(attributeList->create<Real>("Ld", 0)),
        mLq(attributeList->create<Real>("Lq", 0)),
        mLd_t(attributeList->create<Real>("Ld_t", 0)),
        mLq_t(attributeList->create<Real>("Lq_t", 0)),
        mLd_s(attributeList->create<Real>("Ld_s", 0)),
        mLq_s(attributeList->create<Real>("Lq_s", 0)),
        mTd0_t(attributeList->create<Real>("Td0_t", 0)),
        mTq0_t(attributeList->create<Real>("Tq0_t", 0)),
        mTd0_s(attributeList->create<Real>("Td0_s", 0)),
        mTq0_s(attributeList->create<Real>("Tq0_s", 0)),
        mDelta(attributeList->create<Real>("delta_r", 0)),
        mMechTorque(attributeList->create<Real>("T_m", 0)),
        mInertia(attributeList->create<Real>("inertia", 0)),
        mOmMech(attributeList->create<Real>("w_r", 0)),
        mElecActivePower(attributeList->create<Real>("P_elec", 0)),
        mElecReactivePower(attributeList->create<Real>("Q_elec", 0)),
        mMechPower(attributeList->create<Real>("P_mech", 0)),
        mElecTorque(attributeList->create<Real>("T_e", 0)){};

  void setBaseParameters(Real nomPower, Real nomVolt, Real nomFreq);
  void setBaseParameters(Real nomPower, Real nomVolt, Real nomFreq,
                         Real nomFieldCur);

  void calcStateSpaceMatrixDQ();
  Real calcHfromJ(Real J, Real omegaNominal, Int polePairNumber);

public:
  virtual ~SynchronGenerator() {}

  void setBaseAndFundamentalPerUnitParameters(Real nomPower, Real nomVolt,
                                              Real nomFreq, Real nomFieldCur,
                                              Int poleNumber, Real Rs, Real Ll,
                                              Real Lmd, Real Lmq, Real Rfd,
                                              Real Llfd, Real Rkd, Real Llkd,
                                              Real Rkq1, Real Llkq1, Real Rkq2,
                                              Real Llkq2, Real inertia);

  void setBaseAndOperationalPerUnitParameters(
      Real nomPower, Real nomVolt, Real nomFreq, Int poleNumber,
      Real nomFieldCur, Real Rs, Real Ld, Real Lq, Real Ld_t, Real Lq_t,
      Real Ld_s, Real Lq_s, Real Ll, Real Td0_t, Real Tq0_t, Real Td0_s,
      Real Tq0_s, Real inertia);

  void setFundamentalPerUnitParameters(Int poleNumber, Real Rs, Real Ll,
                                       Real Lmd, Real Lmq, Real Rfd, Real Llfd,
                                       Real Rkd, Real Llkd, Real Rkq1,
                                       Real Llkq1, Real Rkq2, Real Llkq2,
                                       Real inertia);

  void applyFundamentalPerUnitParameters();

  void setAndApplyFundamentalPerUnitParameters(Int poleNumber, Real Rs, Real Ll,
                                               Real Lmd, Real Lmq, Real Rfd,
                                               Real Llfd, Real Rkd, Real Llkd,
                                               Real Rkq1, Real Llkq1, Real Rkq2,
                                               Real Llkq2, Real inertia);

  void setOperationalPerUnitParameters(Int poleNumber, Real inertia, Real Rs,
                                       Real Ld, Real Lq, Real Ll, Real Ld_t,
                                       Real Lq_t, Real Ld_s, Real Lq_s,
                                       Real Td0_t, Real Tq0_t, Real Td0_s,
                                       Real Tq0_s);

  void calculateFundamentalFromOperationalParameters();

  void setInitialValues(Real initActivePower, Real initReactivePower,
                        Real initTerminalVolt, Real initVoltAngle,
                        Real initMechPower);

  void setNumericalMethod(NumericalMethod method) { mNumericalMethod = method; }

  std::shared_ptr<Signal::TurbineGovernor> mTurbineGovernor;
  std::shared_ptr<Signal::TurbineGovernorType1> mTurbineGovernorType1;
  std::shared_ptr<Base::Governor> mGovernor;
  std::shared_ptr<Base::Turbine> mTurbine;
  std::shared_ptr<Base::Exciter> mExciter;
  std::shared_ptr<Base::PSS> mPSS;
};

} // namespace Base
} // namespace CPS
