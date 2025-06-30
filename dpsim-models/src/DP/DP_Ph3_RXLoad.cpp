/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/DP/DP_Ph3_RXLoad.h>

using namespace CPS;


DP::Ph3::RXLoad::RXLoad(String uid, String name, Logger::Level logLevel)
    : MNASimPowerComp<Complex>(uid, name, false, true, logLevel),
      mActivePower(mAttributes->create<Matrix>("P")),
      mReactivePower(mAttributes->create<Matrix>("Q")),
      mNomVoltage(mAttributes->create<Real>("V_nom")),
      mReactanceInSeries(false)  {


  mPhaseType = PhaseType::ABC;
  setTerminalNumber(2);

  SPDLOG_LOGGER_INFO(mSLog, "Create {} {}", this->type(), name);

  **mIntfVoltage = MatrixComp::Zero(3, 1); // Interface voltage for logging
  **mIntfCurrent = MatrixComp::Zero(3, 1); // interface current for logging
}

DP::Ph3::RXLoad::RXLoad(String name, Logger::Level logLevel)
    : RXLoad(name, name, logLevel) {}

DP::Ph3::RXLoad::RXLoad(String name, Matrix activePower, Matrix reactivePower, Real volt, Logger::Level logLevel)
    : RXLoad(name, logLevel) {
  **mActivePower = activePower;
  **mReactivePower = reactivePower;
  mPower = MatrixComp::Zero(3, 3);
  mPower << Complex((**mActivePower)(0, 0), (**mReactivePower)(0, 0)),
      Complex((**mActivePower)(0, 1), (**mReactivePower)(0, 1)),
      Complex((**mActivePower)(0, 2), (**mReactivePower)(0, 2)),
      Complex((**mActivePower)(1, 0), (**mReactivePower)(1, 0)),
      Complex((**mActivePower)(1, 1), (**mReactivePower)(1, 1)),
      Complex((**mActivePower)(1, 2), (**mReactivePower)(1, 2)),
      Complex((**mActivePower)(2, 0), (**mReactivePower)(2, 0)),
      Complex((**mActivePower)(2, 1), (**mReactivePower)(2, 1)),
      Complex((**mActivePower)(2, 2), (**mReactivePower)(2, 2));

  **mNomVoltage = volt;

}

DP::Ph3::RXLoad::~RXLoad() = default;

void DP::Ph3::RXLoad::setParameters(Matrix activePower, Matrix reactivePower, Real volt, bool reactanceInSeries) {
  **mActivePower = activePower;
  **mReactivePower = reactivePower;
  mReactanceInSeries = reactanceInSeries;

  if (mReactanceInSeries) {
    setVirtualNodeNumber(1);
  }
  // complex power
  mPower = MatrixComp::Zero(3, 3);
  mPower(0, 0) = {(**mActivePower)(0, 0), (**mReactivePower)(0, 0)};
  mPower(1, 1) = {(**mActivePower)(1, 1), (**mReactivePower)(1, 1)};
  mPower(2, 2) = {(**mActivePower)(2, 2), (**mReactivePower)(2, 2)};

  **mNomVoltage = volt;

  SPDLOG_LOGGER_INFO(mSLog,
                     "\nActive Power [W]: {}"
                     "\nReactive Power [VAr]: {}",
                     Logger::matrixToString(**mActivePower),
                     Logger::matrixToString(**mReactivePower));
  SPDLOG_LOGGER_INFO(mSLog, "Nominal Voltage={} [V]", **mNomVoltage);

  //initPowerFromTerminal = false;
}

void DP::Ph3::RXLoad::initializeFromNodesAndTerminals(Real frequency) {

  // Split Power in three phases
    **mActivePower = Matrix::Zero(3, 3);
    (**mActivePower)(0, 0) = mTerminals[0]->singleActivePower() / 3.;
    (**mActivePower)(1, 1) = mTerminals[0]->singleActivePower() / 3.;
    (**mActivePower)(2, 2) = mTerminals[0]->singleActivePower() / 3.;

    **mReactivePower = Matrix::Zero(3, 3);
    (**mReactivePower)(0, 0) = mTerminals[0]->singleReactivePower() / 3.;
    (**mReactivePower)(1, 1) = mTerminals[0]->singleReactivePower() / 3.;
    (**mReactivePower)(2, 2) = mTerminals[0]->singleReactivePower() / 3.;

    // complex power
    mPower = MatrixComp::Zero(3, 3);
    mPower(0, 0) = {(**mActivePower)(0, 0), (**mReactivePower)(0, 0)};
    mPower(1, 1) = {(**mActivePower)(1, 1), (**mReactivePower)(1, 1)};
    mPower(2, 2) = {(**mActivePower)(2, 2), (**mReactivePower)(2, 2)};

    **mNomVoltage = std::abs(mTerminals[0]->initialSingleVoltage());

    SPDLOG_LOGGER_INFO(mSLog,
                       "\nActive Power [W]: {}"
                       "\nReactive Power [VAr]: {}",
                       Logger::matrixToString(**mActivePower),
                       Logger::matrixToString(**mReactivePower));
    SPDLOG_LOGGER_INFO(mSLog, "Nominal Voltage={} [V]", **mNomVoltage);


  MatrixComp vInitABC = MatrixComp::Zero(3, 1);
  vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * mTerminals[0]->initialSingleVoltage();
  vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
  vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
  **mIntfVoltage = vInitABC.real();

  // Resistance Value
  if ((**mActivePower)(0, 0) != 0) {
    mResistance =
        std::pow(**mNomVoltage / sqrt(3), 2) * (**mActivePower).inverse();
  }

  // Inductor / Capacitor reactance
  if ((**mReactivePower)(0, 0) != 0)
    mReactance =
        std::pow(**mNomVoltage / sqrt(3), 2) * (**mReactivePower).inverse();
  else
    mReactance = Matrix::Zero(1, 1);


  // If reactance chosen create virtual node
  if (mReactanceInSeries) {
    MatrixComp impedance = MatrixComp::Zero(3, 3);
    impedance << Complex(mResistance(0, 0), mReactance(0, 0)),
        Complex(mResistance(0, 1), mReactance(0, 1)),
        Complex(mResistance(0, 2), mReactance(0, 2)),
        Complex(mResistance(1, 0), mReactance(1, 0)),
        Complex(mResistance(1, 1), mReactance(1, 1)),
        Complex(mResistance(1, 2), mReactance(1, 2)),
        Complex(mResistance(2, 0), mReactance(2, 0)),
        Complex(mResistance(2, 1), mReactance(2, 1)),
        Complex(mResistance(2, 2), mReactance(2, 2));
    **mIntfCurrent = (impedance.inverse() * vInitABC).real();

    // Initialization of virtual node
    // Initial voltage of phase B,C is set after A
    MatrixComp vInitTerm0 = MatrixComp::Zero(3, 1);
    vInitTerm0(0, 0) = initialSingleVoltage(0);
    vInitTerm0(1, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_B;
    vInitTerm0(2, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_C;
    mVirtualNodes[0]->setInitialVoltage(vInitTerm0 +
                                        mResistance * **mIntfCurrent);
  }

  // Setup of resistor
  if ((**mActivePower)(0, 0) != 0) {
    mSubResistor =
        std::make_shared<DP::Ph3::Resistor>(**mName + "_res", mLogLevel);
    mSubResistor->setParameters(mResistance);
    if (mReactanceInSeries) {
      mSubResistor->connect({mTerminals[0]->node(), mVirtualNodes[0]});
    } else {
      mSubResistor->connect({SimNode::GND, mTerminals[0]->node()});
    }
    mSubResistor->initialize(mFrequencies);
    mSubResistor->initializeFromNodesAndTerminals(frequency);
    //addMNASubComponent(mSubResistor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                      // MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

    if (!mReactanceInSeries) {
      **mIntfCurrent += mSubResistor->intfCurrent();
    }
  }

  // Case Inductance
  if (mReactance(0, 0) > 0) {
    mInductance = mReactance / (2 * PI * frequency);

    mSubInductor =
        std::make_shared<DP::Ph3::Inductor>(**mName + "_ind", mLogLevel);
    mSubInductor->setParameters(mInductance);
    if (mReactanceInSeries) {
      mSubInductor->connect({SimNode::GND, mVirtualNodes[0]});
    } else {
      mSubInductor->connect({SimNode::GND, mTerminals[0]->node()});
    }
    mSubInductor->initialize(mFrequencies);
    mSubInductor->initializeFromNodesAndTerminals(frequency);
    //addMNASubComponent(mSubInductor, MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     //  MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);

    if (!mReactanceInSeries) {
      **mIntfCurrent += mSubInductor->intfCurrent();
    }
  }

  // Case Capacitance
  else if (mReactance(0, 0) < 0) {
    mCapacitance = -1 / (2 * PI * frequency) * mReactance.inverse();

    mSubCapacitor =
        std::make_shared<DP::Ph3::Capacitor>(**mName + "_cap", mLogLevel);
    mSubCapacitor->setParameters(mCapacitance);
    if (mReactanceInSeries) {
      mSubCapacitor->connect({SimNode::GND, mVirtualNodes[0]});
    } else {
      mSubCapacitor->connect({SimNode::GND, mTerminals[0]->node()});
    }
    mSubCapacitor->initialize(mFrequencies);
    mSubCapacitor->initializeFromNodesAndTerminals(frequency);
    //addMNASubComponent(mSubCapacitor,
                     //  MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT,
                     //  MNA_SUBCOMP_TASK_ORDER::TASK_BEFORE_PARENT, true);
    if (!mReactanceInSeries) {
      **mIntfCurrent += mSubCapacitor->intfCurrent();
    }
  }

/*
  // Logging
   SPDLOG_LOGGER_INFO(
      mSLog,
      "\n--- Initialization from powerflow ---"
      "\nVoltage across: {:s}"
      "\nCurrent: {:s}"
      "\nTerminal 0 voltage: {:s}"
      "\nActive Power: {:s}"
      "\nReactive Power: {:s}"
      "\nResistance: {:s}"
      "\nReactance: {:s}"
      "\n--- Initialization from powerflow finished ---",
      Logger::matrixToString(**mIntfVoltage),
      Logger::matrixToString(**mIntfCurrent),
      Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
      Logger::matrixToString(**mActivePower),
      Logger::matrixToString(**mReactivePower),
      Logger::matrixToString(mResistance),
      Logger::matrixToString(mReactance));
  mSLog->flush();
  */
}

void DP::Ph3::RXLoad::mnaCompInitialize(Real omega, Real timeStep,
                                          Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  **mRightVector = Matrix::Zero(0, 0);
}

void DP::Ph3::RXLoad::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  MatrixFixedSizeComp<3, 3> conductance = Matrix::Zero(3, 3);
  conductance.real() = mResistance.inverse();
  MNAStampUtils::stampAdmittanceMatrix(
      conductance, systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1),
      terminalNotGrounded(0), terminalNotGrounded(1), mSLog);
}

void DP::Ph3::RXLoad::mnaCompAddPostStepDependencies(
    AttributeBase::List &prevStepDependencies,
    AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mIntfVoltage);
  modifiedAttributes.push_back(mIntfCurrent);
}

void DP::Ph3::RXLoad::mnaCompPostStep(Real time, Int timeStepCount,
                                        Attribute<Matrix>::Ptr &leftVector) {
  mnaCompUpdateVoltage(**leftVector);
  mnaCompUpdateCurrent(**leftVector);
}

void DP::Ph3::RXLoad::mnaCompUpdateVoltage(const Matrix &leftVector) {
  // Voltage across component is defined as V1 - V0
  **mIntfVoltage = MatrixComp::Zero(3, 1);
  if (terminalNotGrounded(1)) {
    (**mIntfVoltage)(0, 0) =
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(1, 0));
    (**mIntfVoltage)(1, 0) =
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(1, 1));
    (**mIntfVoltage)(2, 0) =
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(1, 2));
  }
  if (terminalNotGrounded(0)) {
    (**mIntfVoltage)(0, 0) =
        (**mIntfVoltage)(0, 0) -
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(0, 0));
    (**mIntfVoltage)(1, 0) =
        (**mIntfVoltage)(1, 0) -
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(0, 1));
    (**mIntfVoltage)(2, 0) =
        (**mIntfVoltage)(2, 0) -
        Math::complexFromVectorElement(leftVector, matrixNodeIndex(0, 2));
  }

  SPDLOG_LOGGER_DEBUG(mSLog, "Voltage A: {} < {}",
                      std::abs((**mIntfVoltage)(0, 0)),
                      std::arg((**mIntfVoltage)(0, 0)));
}

void DP::Ph3::RXLoad::mnaCompUpdateCurrent(const Matrix &leftVector) {
  **mIntfCurrent = mResistance.inverse() * **mIntfVoltage;

  SPDLOG_LOGGER_DEBUG(mSLog, "Current A: {} < {}",
                      std::abs((**mIntfCurrent)(0, 0)),
                      std::arg((**mIntfCurrent)(0, 0)));
}
