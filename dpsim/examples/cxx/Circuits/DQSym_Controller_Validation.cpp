// SPDX-License-Identifier: MPL-2.0
#include "../Examples.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <dpsim-models/Signal/DQSymControllerBlocks.h>

using namespace CPS;

namespace {
void near(Real actual, Real expected, Real tolerance, const String &name) {
  const Real error = std::abs(actual - expected);
  std::cout << name << ": actual=" << actual << ", expected=" << expected
            << ", error=" << error << '\n';
  if (!std::isfinite(actual) || error > tolerance)
    throw std::runtime_error(name + " failed.");
}

void validatePI() {
  auto pi = Signal::DQSymPIController::make("PI");
  pi->setParameters(2.0, 4.0, -1.0, 1.0);
  pi->setInitialState(0.0);
  pi->initialize(0.01);

  pi->mError->set(0.2);
  pi->step();
  near(**pi->mUnsaturatedOutput, 0.404, 1e-12,
       "PI proportional and trapezoidal first half-step");

  pi->mError->set(10.0);
  pi->step();
  near(**pi->mOutput, 1.0, 0.0, "PI upper saturation");
  const Real heldState = **pi->mIntegratorState;
  for (UInt i = 0; i < 20; ++i)
    pi->step();
  near(**pi->mIntegratorState, heldState, 0.0,
       "PI conditional anti-windup hold");

  pi->mError->set(-0.5);
  pi->step();
  if (!(**pi->mOutput < 1.0))
    throw std::runtime_error("PI did not recover from upper saturation.");

  pi->mEnable->set(false);
  const Real disabledState = **pi->mIntegratorState;
  pi->mError->set(1.0);
  for (UInt i = 0; i < 10; ++i)
    pi->step();
  near(**pi->mIntegratorState, disabledState, 0.0, "PI disabled state hold");

  pi->mError->set(0.0);
  pi->mEnable->set(true);
  pi->step();
  near(**pi->mOutput, disabledState, 1e-12,
       "PI state-hold bumpless re-enable at zero error");
}

void validateFilter() {
  auto filter = Signal::DQSymSecondOrderFilter::make("Filter");
  const Real frequency = 10.0;
  const Real omega = 2.0 * PI * frequency;
  const Real timeStep = 1e-5;
  filter->setParameters(frequency, 1.0);
  filter->setInitialValue(0.0);
  filter->initialize(timeStep);
  filter->mInput->set(1.0);
  const Real finalTime = 0.05;
  for (UInt i = 0; i < static_cast<UInt>(finalTime / timeStep); ++i)
    filter->step();
  const Real analytical =
      1.0 - (1.0 + omega * finalTime) * std::exp(-omega * finalTime);
  near(**filter->mOutput, analytical, 3e-4,
       "critically damped second-order filter");

  filter->setInitialValue(3.0);
  filter->initialize(timeStep);
  filter->step();
  near(**filter->mOutput, 3.0, 1e-12, "filter steady-state initialization");
}

void validateCurrentController() {
  auto controller = Signal::DQSymCurrentController::make("Current");
  Signal::DQSymCurrentControllerParameters parameters;
  parameters.kp = 0.6;
  parameters.ki = 6.0;
  parameters.rFeedforwardPu = 0.0015 / 2.0;
  parameters.lFeedforwardPu = 0.15 / 2.0;
  parameters.nominalFrequencyHz = 50.0;
  parameters.lowerLimitPu = -2.0;
  parameters.upperLimitPu = 2.0;
  controller->setParameters(parameters);
  controller->setInitialIntegratorStates(0.0, 0.0);
  controller->initialize(40e-6);
  controller->mIdReference->set(0.2);
  controller->mIqReference->set(-0.1);
  controller->mId->set(0.2);
  controller->mIq->set(-0.1);
  controller->mVd->set(1.0);
  controller->mVq->set(0.0);
  controller->mFrequencyHz->set(50.0);
  controller->step();

  near(**controller->mVdReference,
       1.0 + parameters.rFeedforwardPu * 0.2 + parameters.lFeedforwardPu * 0.1,
       1e-12, "d-axis feedforward and negative q cross-coupling");
  near(**controller->mVqReference,
       parameters.rFeedforwardPu * -0.1 + parameters.lFeedforwardPu * 0.2,
       1e-12, "q-axis feedforward and positive d cross-coupling");

  controller->mIdReference->set(10.0);
  controller->step();
  near(**controller->mVdReference, 2.0, 0.0,
       "current-controller d-axis saturation");
  if (!(**controller->mDSaturated))
    throw std::runtime_error("Current-controller saturation flag not set.");
}

void validateOuterControllers() {
  auto active = Signal::DQSymOuterController::make("Active");
  active->setParameters(Signal::DQSymOuterLoopType::ActivePowerToDcVoltage,
                        0.5 / 3.0, 1.0, 0.8, 1.2, 1.0, 640e3, 1.0);
  active->mReference->set(0.4);
  active->mMeasurement->set(0.5);
  active->initialize(40e-6);
  if (!(**active->mOutput > 640e3))
    throw std::runtime_error(
        "DQsym active-power error did not increase Vdc_ref.");

  auto dcVoltage = Signal::DQSymOuterController::make("DcVoltage");
  dcVoltage->setParameters(Signal::DQSymOuterLoopType::DcVoltageToDCurrent, 4.0,
                           100.0, -2.0, 2.0, 640e3, 1.0, 0.0);
  dcVoltage->mReference->set(640e3);
  dcVoltage->mMeasurement->set(646.4e3);
  dcVoltage->initialize(40e-6);
  if (!(**dcVoltage->mOutput > 0.0))
    throw std::runtime_error(
        "DQsym positive DC-voltage error did not produce positive id_ref.");

  auto reactive = Signal::DQSymOuterController::make("Reactive");
  reactive->setParameters(Signal::DQSymOuterLoopType::ReactivePowerToQCurrent,
                          0.5 / 3.0, 1.0, -0.25, 0.25, 1.0, 1.0, 0.0);
  reactive->mReference->set(0.1);
  reactive->mMeasurement->set(0.0);
  reactive->initialize(40e-6);
  if (!(**reactive->mOutput > 0.0))
    throw std::runtime_error(
        "DQsym positive Q error did not produce positive PI output.");

  reactive->mEnable->set(false);
  const Real held = **reactive->mIntegratorState;
  reactive->mMeasurement->set(-1.0);
  reactive->step();
  near(**reactive->mIntegratorState, held, 0.0,
       "outer-loop disabled state hold");
}

void validateModulation() {
  auto modulation = Signal::DQSymModulation::make("Modulation");
  const Real vdc = 640e3;
  const Real vac = 333e3;
  modulation->setParameters(vdc, vac);
  modulation->mDcVoltage->set(vdc);
  modulation->mAngle->set(0.0);

  const Real inverseScale = vdc * 0.5 * std::sqrt(3.0 / 2.0) / vac;
  modulation->mVdCommand->set(inverseScale);
  modulation->mVqCommand->set(0.0);
  modulation->step();
  near(**modulation->mDCommand, 1.0, 1e-12, "DQsym modulation scaling");
  near((**modulation->mAbcCommand)(0, 0), 1.0, 1e-12,
       "inverse-dq phase a at zero angle");
  near((**modulation->mAbcCommand)(1, 0), -0.5, 1e-12,
       "inverse-dq phase b at zero angle");
  near((**modulation->mAbcCommand)(2, 0), -0.5, 1e-12,
       "inverse-dq phase c at zero angle");

  modulation->mVdCommand->set(100.0 * inverseScale);
  modulation->step();
  near(**modulation->mDCommand, 2.0, 0.0, "DQsym modulation axis saturation");
}
} // namespace

int main() {
  validatePI();
  validateFilter();
  validateOuterControllers();
  validateCurrentController();
  validateModulation();
  std::cout << "All independently verifiable DQsym controller validations "
               "passed.\n";
  return 0;
}
