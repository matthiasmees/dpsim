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

struct Fixture {
  EMT::Ph3::SSN_MMC::Ptr plant;
  std::shared_ptr<EMT::Ph3::SSN_MMCStation> station;
  Real voltageBase;

  Fixture() {
    plant = EMT::Ph3::SSN_MMC::make("Plant");
    plant->setParameters(50.0, 333e3, 640e3, 0.05, 0.5, 0.01, 400, 0.02, 0.1);
    plant->setOutputCurrentController(1.0, 10.0);
    plant->setCirculatingCurrentController(1.0, 10.0);
    plant->setZeroSequenceCurrentController(1.0, 10.0);
    voltageBase = std::sqrt(2.0 / 3.0) * 333e3;
    Matrix voltage(3, 1);
    voltage << voltageBase, -0.5 * voltageBase, -0.5 * voltageBase;
    plant->acTerminalVoltageAttribute()->set(voltage);
    plant->acTerminalCurrentAttribute()->set(Matrix::Zero(3, 1));
    plant->dcPositiveVoltageAttribute()->set(640e3);
    plant->dcNegativeVoltageAttribute()->set(0.0);
    plant->dcVoltageAttribute()->set(640e3);
    plant->dcCurrentAttribute()->set(0.0);

    station = EMT::Ph3::SSN_MMCStation::make("Station", plant);
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
    station->initializeFromOperatingPoint(0.0, 2.0 * PI * 50.0, 0.0, 0.0,
                                          640e3);
    station->initialize(40e-6);
    require(station->requestEnable(), "Fixture station enable failed.");
  }
};

[[maybe_unused]] void openLoopValidation() {
  Fixture fixture;
  const auto tasks = fixture.station->getTasks();
  for (const auto &task : tasks)
    task->execute(0.0, 0);
  const Real expectedD =
      **fixture.station->mVdReferencePu * fixture.voltageBase;
  require(std::abs((**fixture.station->mPlantDifferentialVoltageCommand)(0, 0) -
                   expectedD) < 1e-9,
          "Per-unit to plant differential-voltage conversion failed.");
  const Real expectedInsertion =
      2.0 * expectedD / **fixture.plant->dcVoltageAttribute();
  require(std::abs(**fixture.station->mModulationMagnitude -
                   std::abs(expectedInsertion)) < 1e-12,
          "DQsym modulation and plant voltage-command scaling disagree.");
  require((**fixture.station->mConverterPhaseCommand).allFinite(),
          "Open-loop phase command is non-finite.");
  std::cout << "MMCStation open-loop signal integration passed.\n";
}

[[maybe_unused]] void currentValidation() {
  Fixture fixture;
  const auto tasks = fixture.station->getTasks();
  fixture.station->mIdReferencePu->set(0.1);
  tasks[2]->execute(0.0, 1);
  require(**fixture.station->mVdReferencePu > 1.0,
          "Positive d-current error did not increase d-voltage command.");
  fixture.station->mIdReferencePu->set(10.0);
  tasks[2]->execute(0.0, 2);
  require(**fixture.station->mCurrentDSaturated &&
              **fixture.station->mVdReferencePu == 2.0,
          "Current-loop upper saturation failed.");
  require(std::isfinite(**fixture.station->mVdReferencePu),
          "Current-loop output is non-finite.");
  std::cout << "MMCStation current-controller signal integration passed.\n";
}

[[maybe_unused]] void outerValidation() {
  Fixture fixture;
  const auto tasks = fixture.station->getTasks();
  fixture.station->block();
  fixture.station->setControlMode(
      EMT::Ph3::SSN_MMCStation::ControlMode::ActivePowerReactivePower);
  fixture.station->initializeFromOperatingPoint(0.0, 2.0 * PI * 50.0, 0.0, 0.0,
                                                640e3);
  fixture.station->initialize(40e-6);
  require(fixture.station->requestEnable(), "P-mode re-enable failed.");
  fixture.station->mActivePowerReferencePu->set(-0.01);
  for (UInt index = 0; index < 20; ++index) {
    tasks[0]->execute(0.0, index);
    tasks[1]->execute(0.0, index);
  }
  require(**fixture.station->mIdReferencePu < 0.0,
          "Active-power cascade did not produce the expected negative id_ref.");

  fixture.station->block();
  fixture.station->setControlMode(
      EMT::Ph3::SSN_MMCStation::ControlMode::DCVoltageReactivePower);
  fixture.station->initializeFromOperatingPoint(0.0, 2.0 * PI * 50.0, 0.0, 0.0,
                                                640e3);
  fixture.station->initialize(40e-6);
  require(fixture.station->requestEnable(), "Vdc-mode re-enable failed.");
  fixture.plant->dcVoltageAttribute()->set(646.4e3);
  for (UInt index = 0; index < 200; ++index) {
    tasks[0]->execute(0.0, index);
    tasks[1]->execute(0.0, index);
  }
  require(**fixture.station->mIdReferencePu < 0.0,
          "Positive DC overvoltage did not produce the SSN_MMC plant's "
          "required negative id_ref.");
  std::cout << "MMCStation outer-controller signal integration passed.\n";
}
} // namespace

int main() {
#if defined(MMCSTATION_OPEN_LOOP)
  openLoopValidation();
#elif defined(MMCSTATION_CURRENT)
  currentValidation();
#elif defined(MMCSTATION_OUTER)
  outerValidation();
#else
#error "Select one MMCStation validation scenario."
#endif
  return 0;
}
