/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/DP/DP_Ph3_Capacitor.h>
#include <dpsim-models/DP/DP_Ph3_Inductor.h>
#include <dpsim-models/DP/DP_Ph3_Resistor.h>
#include <dpsim-models/Definitions.h>
#include <dpsim-models/Logger.h>
#include <dpsim-models/MNASimPowerComp.h>
#include <dpsim-models/Solver/MNAInterface.h>

namespace CPS {
namespace DP {
namespace Ph3 {
///
class RXLoad : public MNASimPowerComp<Complex>,
                public SharedFactory<RXLoad> {


protected:
  /// Power [Watt]
  MatrixComp mPower;
  /// Resistance [Ohm]
  Matrix mResistance;
  /// Reactance [Ohm]
  Matrix mReactance;
  /// Inductance [H]
  Matrix mInductance;
  /// Capacitance [F]
  Matrix mCapacitance;

  /// If set to true, the reactance is in series with the resistor. Otherwise it is parallel to the resistor.
  Bool mReactanceInSeries;

  /// Internal inductor
  std::shared_ptr<DP::Ph3::Inductor> mSubInductor;
  /// Internal capacitor
  std::shared_ptr<DP::Ph3::Capacitor> mSubCapacitor;
  /// Internal resistance
  std::shared_ptr<DP::Ph3::Resistor> mSubResistor;
  /// Right side vectors of subcomponents
  std::vector<const Matrix *> mRightVectorStamps;

public:
/// Active power [Watt]
  const Attribute<Matrix>::Ptr mActivePower;
  /// Reactive power [VAr]
  const Attribute<Matrix>::Ptr mReactivePower;
  /// Nominal voltage [V]
  const Attribute<Real>::Ptr mNomVoltage;

  /// Defines UID, name and logging level
  RXLoad(String uid, String name, Logger::Level logLevel = Logger::Level::off);
  /// Defines name, component parameters and logging level
  RXLoad(String name, Logger::Level logLevel = Logger::Level::off);
  /// Defines name, component parameters and logging level
  RXLoad(String name, Matrix activePower, Matrix reactivePower, Real volt,
         Logger::Level logLevel = Logger::Level::off);
  virtual ~RXLoad();
  

  // #### General ####

  /// Set model specific parameters
  void setParameters(Matrix activePower, Matrix reactivePower, Real volt, bool reactanceInSeries = false);

  /// Initializes component from power flow data
  void initializeFromNodesAndTerminals(Real frequency) override;


  // #### MNA section ####
  ///
  void mnaCompInitialize(Real omega, Real timeStep,
                         Attribute<Matrix>::Ptr leftVector) override;
  /// Stamps system matrix
  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  ///
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;
  ///
  void mnaCompUpdateCurrent(const Matrix &leftVector) override;

  /// Add MNA post step dependencies
  void
  mnaCompAddPostStepDependencies(AttributeBase::List &prevStepDependencies,
                                 AttributeBase::List &attributeDependencies,
                                 AttributeBase::List &modifiedAttributes,
                                 Attribute<Matrix>::Ptr &leftVector) override;
  void mnaCompPostStep(Real time, Int timeStepCount,
                       Attribute<Matrix>::Ptr &leftVector) override;
};
} // namespace Ph3
} // namespace DP
} // namespace CPS
