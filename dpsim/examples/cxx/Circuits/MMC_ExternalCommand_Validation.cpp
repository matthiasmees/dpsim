// SPDX-License-Identifier: MPL-2.0
#include "../Examples.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <dpsim-models/EMT/EMT_Ph3_SSN_MMCStation.h>

using namespace CPS;

namespace {
void require(Bool condition, const String &message) {
  if (!condition)
    throw std::runtime_error(message);
}

EMT::Ph3::SSN_MMC::Ptr makePlant(const String &name) {
  auto plant = EMT::Ph3::SSN_MMC::make(name);
  plant->setParameters(50.0, 333e3, 640e3, 0.05, 0.5, 0.01, 400, 0.02, 0.1);
  plant->setOutputCurrentController(1.0, 10.0);
  plant->setCirculatingCurrentController(1.0, 10.0);
  plant->setZeroSequenceCurrentController(1.0, 10.0);
  return plant;
}

void validatePlantInterface() {
  auto plant = makePlant("Plant");
  require(plant->controlSource() ==
              EMT::Ph3::SSN_MMC::ControlSource::InternalControllers,
          "External command became the default.");
  const Matrix defaultDerivative = plant->getStateDerivative();
  plant->setExternalDifferentialVoltageCommand(1000.0, -500.0);
  const Matrix stillDefaultDerivative = plant->getStateDerivative();
  require((defaultDerivative - stillDefaultDerivative).norm() == 0.0,
          "Inactive external command changed default MMC equations.");

  plant->setControlSource(
      EMT::Ph3::SSN_MMC::ControlSource::ExternalDifferentialVoltage);
  require(**plant->externalCommandActiveAttribute(),
          "External control source did not become active.");
  require(std::abs((**plant->externalDifferentialVoltageAttribute())(0, 0) -
                   1000.0) < 1e-12 &&
              std::abs((**plant->externalDifferentialVoltageAttribute())(1, 0) +
                       500.0) < 1e-12,
          "External command did not propagate to the plant attribute.");

  Bool rejected = false;
  try {
    plant->setExternalDifferentialVoltageCommand(
        std::numeric_limits<Real>::quiet_NaN(), 0.0);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "Non-finite MMC external command was not rejected.");

  auto initialized = makePlant("InitializedPlant");
  auto ac = SimNode<Real>::make("ac", PhaseType::ABC);
  auto dcPositive = SimNode<Real>::make("dcPositive", PhaseType::DC);
  auto dcNegative = SimNode<Real>::make("dcNegative", PhaseType::DC);
  ac->setInitialVoltage(Complex(333e3 / RMS3PH_TO_PEAK1PH, 0.0));
  dcPositive->setInitialVoltage(Complex(320e3, 0.0));
  dcNegative->setInitialVoltage(Complex(-320e3, 0.0));
  initialized->connect({ac, dcPositive, dcNegative});
  initialized->setOperatingPointInitialization(true, 25, 1e-8);
  std::static_pointer_cast<SimPowerComp<Real>>(initialized)
      ->initializeFromNodesAndTerminals(50.0);
  const Matrix initializedDefaultDerivative = initialized->getStateDerivative();
  initialized->setExternalDifferentialVoltageCommand(100e3, 10e3);
  initialized->setControlSource(
      EMT::Ph3::SSN_MMC::ControlSource::ExternalDifferentialVoltage);
  const Matrix initializedExternalDerivative =
      initialized->getStateDerivative();
  require(
      (initializedExternalDerivative - initializedDefaultDerivative).norm() >
          1e-6,
      "Initialized MMC did not consume the external command.");
  std::cout << "MMC default/external command interface passed.\n";
}

void validateTransform() {
  auto adapter = Signal::ExternallyAngledDQAdapter::make("Transform");
  const Real theta = 0.37;
  adapter->mAngle->set(theta);
  adapter->mAngularFrequency->set(2.0 * PI * 50.0);
  const Matrix voltage = adapter->dqToAbc(2.0, -0.5);
  const Matrix current = adapter->dqToAbc(0.4, -0.2);
  adapter->mVoltageAbc->set(voltage);
  adapter->mCurrentAbc->set(current);
  adapter->step();
  require(std::abs(**adapter->mVd - 2.0) < 1e-12 &&
              std::abs(**adapter->mVq + 0.5) < 1e-12 &&
              std::abs(**adapter->mId - 0.4) < 1e-12 &&
              std::abs(**adapter->mIq + 0.2) < 1e-12,
          "Externally angled dq round trip failed.");
  require(std::abs(**adapter->mActivePower - 1.35) < 1e-12,
          "Transform active-power sign/scaling failed.");
  require(std::abs(**adapter->mReactivePower - 0.3) < 1e-12,
          "Transform reactive-power sign/scaling failed.");
  std::cout << "Explicit transform boundary passed.\n";
}

void validateStationAndDependencies() {
  auto plant = makePlant("StationPlant");
  Matrix vabc(3, 1);
  const Real peak = std::sqrt(2.0 / 3.0) * 333e3;
  vabc << peak, -0.5 * peak, -0.5 * peak;
  plant->acTerminalVoltageAttribute()->set(vabc);
  plant->acTerminalCurrentAttribute()->set(Matrix::Zero(3, 1));
  plant->dcPositiveVoltageAttribute()->set(640e3);
  plant->dcNegativeVoltageAttribute()->set(0.0);
  plant->dcVoltageAttribute()->set(640e3);
  plant->dcCurrentAttribute()->set(0.0);

  auto station = EMT::Ph3::SSN_MMCStation::make("Station", plant);
  EMT::Ph3::SSN_MMCStationParameters parameters;
  parameters.nominalPower = 1e9;
  parameters.nominalAcLineLineRms = 333e3;
  parameters.nominalDcVoltage = 640e3;
  parameters.nominalFrequencyHz = 50.0;
  parameters.controllerTimeStep = 40e-6;
  parameters.armResistancePu = 0.0015;
  parameters.armInductancePu = 0.15;
  station->setParameters(parameters);
  station->setControlMode(
      EMT::Ph3::SSN_MMCStation::ControlMode::DCVoltageReactivePower);
  station->initializeFromOperatingPoint(0.0, 2.0 * PI * 50.0, 0.0, 0.0, 640e3);
  station->initialize(40e-6);
  if (!station->requestEnable())
    throw std::runtime_error(
        "Zero-error operating point did not enable bumplessly; diagnostic=" +
        std::to_string(**station->mEnableDiagnostic) +
        ", P=" + std::to_string(**station->mFilteredActivePowerPu) +
        ", Q=" + std::to_string(**station->mFilteredReactivePowerPu) +
        ", Vdc=" + std::to_string(**station->mFilteredDcVoltage));
  require(plant->controlSource() ==
              EMT::Ph3::SSN_MMC::ControlSource::ExternalDifferentialVoltage,
          "Station enable did not explicitly select external plant control.");

  const auto tasks = station->getTasks();
  require(tasks.size() == 4, "Station must expose four ordered task layers.");
  require(tasks[0]->getModifiedAttributes().size() > 0 &&
              tasks[1]->getAttributeDependencies().size() > 0 &&
              tasks[2]->getAttributeDependencies().size() > 0 &&
              tasks[3]->getAttributeDependencies().size() > 0,
          "Station tasks do not declare attribute dependencies.");
  tasks[0]->execute(0.0, 0);
  tasks[1]->execute(0.0, 0);
  tasks[2]->execute(0.0, 0);
  tasks[3]->execute(0.0, 0);
  require((**station->mPlantDifferentialVoltageCommand).allFinite(),
          "Station command is non-finite.");

  station->block();
  require(plant->controlSource() ==
              EMT::Ph3::SSN_MMC::ControlSource::InternalControllers,
          "Station block did not restore the default plant control source.");
  std::cout << "Station initialization, enable, and task chain passed.\n";
}
} // namespace

int main() {
  validatePlantInterface();
  validateTransform();
  validateStationAndDependencies();
  std::cout << "MMC external-command validation passed.\n";
  return 0;
}
