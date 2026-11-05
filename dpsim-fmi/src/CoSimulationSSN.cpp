// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/CoSimulationSSN.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;
namespace DPsim::FMI {
namespace {
std::size_t indexOf(const std::vector<std::string> &names,
                    const std::string &name) {
  const auto found = std::find(names.begin(), names.end(), name);
  if (found == names.end())
    throw std::invalid_argument("FMI variable does not exist: " + name);
  return static_cast<std::size_t>(found - names.begin());
}
} // namespace

CoSimulationSSN::CoSimulationSSN(String name, CoSimulation::Ptr model,
                                 CoSimulationElectricalPort port)
    : MNASimPowerComp(std::move(name), true, true), mModel(std::move(model)),
      mPort(std::move(port)), mSign(mPort.currentSignToPassive) {
  if (!mModel)
    throw std::invalid_argument("CoSimulationSSN requires an imported FMU");
  if (mSign != 1 && mSign != -1)
    throw std::invalid_argument("Current sign must be +1 or -1");
  const auto channels = mPort.voltageInputs.size();
  if ((channels != 1 && channels != 3) ||
      mPort.currentOutputs.size() != channels)
    throw std::invalid_argument(
        "CoSimulationSSN requires one or three terminal channels");
  const auto &description = mModel->description();
  std::vector<bool> usedInput(description.inputNames.size(), false);
  std::vector<bool> usedOutput(description.outputNames.size(), false);
  for (const auto &name : mPort.voltageInputs) {
    const auto index = indexOf(description.inputNames, name);
    if (usedInput[index])
      throw std::invalid_argument("terminal input mappings must be distinct");
    usedInput[index] = true;
    mVoltageIndices.push_back(index);
  }
  for (const auto &name : mPort.currentOutputs) {
    const auto index = indexOf(description.outputNames, name);
    if (usedOutput[index])
      throw std::invalid_argument("terminal output mappings must be distinct");
    usedOutput[index] = true;
    mCurrentIndices.push_back(index);
  }
  mPhaseType = channels == 3 ? PhaseType::ABC : PhaseType::Single;
  setTerminalNumber(2);
  mInputs = mAttributes->create<Matrix>("fmi.inputs");
  mOutputs = mAttributes->create<Matrix>("fmi.outputs");
  **mIntfVoltage = Matrix::Zero(channels, 1);
  **mIntfCurrent = Matrix::Zero(channels, 1);
  mStepAttribute = mAttributes->create<Real>("cosim_do_steps");
  mInternalStepAttribute = mAttributes->create<Real>("fmi_internal_steps");
  mControls.resize(description.inputNames.size());
  mInputFunctions.resize(description.inputNames.size());
  mInputEvents.resize(description.inputNames.size());
  **mInputs = mModel->inputs();
  **mOutputs = mModel->outputs();
  for (std::size_t index = 0; index < description.inputNames.size(); ++index) {
    auto input = mAttributes->createDynamic<Real>(
        "fmi.input." + description.inputNames[index]);
    input->setReference(mInputs->deriveCoeff<Real>(index, 0));
    if (!usedInput[index]) {
      mControls[index] =
          mAttributes->create<Real>("control." + description.inputNames[index]);
      **mControls[index] = description.initialInputs(index, 0);
    }
  }
  for (std::size_t index = 0; index < description.outputNames.size(); ++index) {
    auto output = mAttributes->createDynamic<Real>(
        "fmi.output." + description.outputNames[index]);
    output->setReference(mOutputs->deriveCoeff<Real>(index, 0));
  }
  for (std::size_t phase = 0; phase < channels; ++phase)
    (**mIntfCurrent)(phase, 0) =
        mSign * (**mOutputs)(mCurrentIndices[phase], 0);
  **mStepAttribute = 0;
  **mInternalStepAttribute = 0;
}

std::shared_ptr<CoSimulationSSN>
CoSimulationSSN::importFMU(String name, const std::string &path,
                           CoSimulationElectricalPort port,
                           CoSimulationOptions options) {
  return CoSimulationSSN::make(std::move(name),
                               std::make_shared<CoSimulation>(path, options),
                               std::move(port));
}

SimPowerComp<Real>::Ptr CoSimulationSSN::clone(String name) {
  throw std::logic_error(
      "A running Co-Simulation FMU cannot be cloned implicitly; import a "
      "separate FMU instance for each component");
}

void CoSimulationSSN::connect(SimNode<Real>::Ptr node) {
  if (!node)
    throw std::invalid_argument("cannot connect a Co-Simulation FMU to null");
  SimPowerComp<Real>::connect({SimNode<Real>::GND, std::move(node)});
}

void CoSimulationSSN::initialize(Matrix frequencies) {
  SimPowerComp<Real>::initialize(std::move(frequencies));
}

void CoSimulationSSN::validateTerminals() {
  for (UInt terminalIndex = 0; terminalIndex < 2; ++terminalIndex) {
    const auto terminalNode = node(terminalIndex);
    if (!terminalNode ||
        (!terminalNode->isGround() && terminalNode->phaseType() != mPhaseType))
      throw std::invalid_argument(
          "Co-Simulation FMU terminal phase type does not match its port");
  }
}

void CoSimulationSSN::initializeFromNodesAndTerminals(Real) {
  validateTerminals();
  const auto channels = static_cast<Eigen::Index>(mVoltageIndices.size());
  Matrix terminalVoltage = Matrix::Zero(channels, 1);
  if (!node(1)->isGround())
    terminalVoltage += initialVoltage(1).real();
  if (!node(0)->isGround())
    terminalVoltage -= initialVoltage(0).real();
  if (channels == 3)
    terminalVoltage *= RMS3PH_TO_PEAK1PH;
  if (terminalVoltage.allFinite()) {
    **mIntfVoltage = terminalVoltage;
    for (Eigen::Index phase = 0; phase < channels; ++phase)
      (**mInputs)(mVoltageIndices[phase], 0) = terminalVoltage(phase, 0);
  }
  for (std::size_t phase = 0; phase < mCurrentIndices.size(); ++phase)
    (**mIntfCurrent)(phase, 0) =
        mSign * (**mOutputs)(mCurrentIndices[phase], 0);
}

void CoSimulationSSN::setInputFunction(const std::string &name,
                                       std::function<double(double)> function,
                                       std::vector<double> events) {
  const auto index = indexOf(mModel->description().inputNames, name);
  if (mControls[index].getPtr() == nullptr || !function)
    throw std::invalid_argument(
        "Co-Simulation input functions require an external control input");
  const double initial = function(0);
  const double declared = mModel->description().initialInputs(index, 0);
  if (!std::isfinite(initial) ||
      std::abs(initial - declared) > 1e-12 * (1 + std::abs(initial)))
    throw std::invalid_argument(
        "control function at t=0 must match the FMI input start value");
  std::sort(events.begin(), events.end());
  if (std::adjacent_find(events.begin(), events.end()) != events.end())
    throw std::invalid_argument("input discontinuity times must be unique");
  for (double event : events)
    if (!std::isfinite(event) || event <= 0)
      throw std::invalid_argument(
          "input discontinuity times must be finite and positive");
  mInputFunctions[index] = std::move(function);
  mInputEvents[index] = std::move(events);
}

void CoSimulationSSN::mnaCompInitialize(Real, Real step,
                                        Attribute<Matrix>::Ptr) {
  validateTerminals();
  if (!std::isfinite(step) || step <= 0)
    throw std::invalid_argument("Co-Simulation communication step is invalid");
  mTimeStep = step;
  for (const auto &events : mInputEvents)
    for (double event : events)
      if (std::abs(event / step - std::round(event / step)) > 1e-7)
        throw std::invalid_argument(
            "Co-Simulation input events must align with the DPsim timestep");
  updateMatrixNodeIndices();
}

void CoSimulationSSN::mnaCompApplyRightSideVectorStamp(Matrix &vector) {
  for (UInt channel = 0; channel < mVoltageIndices.size(); ++channel) {
    const UInt row0 = mPhaseType == PhaseType::ABC ? matrixNodeIndex(0, channel)
                                                   : matrixNodeIndex(0);
    const UInt row1 = mPhaseType == PhaseType::ABC ? matrixNodeIndex(1, channel)
                                                   : matrixNodeIndex(1);
    if (terminalNotGrounded(0))
      Math::setVectorElement(vector, row0, (**mIntfCurrent)(channel, 0));
    if (terminalNotGrounded(1))
      Math::setVectorElement(vector, row1, -(**mIntfCurrent)(channel, 0));
  }
}

void CoSimulationSSN::updateTerminalVoltage(const Matrix &solution) {
  const auto channels = static_cast<Eigen::Index>(mVoltageIndices.size());
  **mIntfVoltage = Matrix::Zero(channels, 1);
  for (UInt channel = 0; channel < mVoltageIndices.size(); ++channel) {
    const UInt row0 = mPhaseType == PhaseType::ABC ? matrixNodeIndex(0, channel)
                                                   : matrixNodeIndex(0);
    const UInt row1 = mPhaseType == PhaseType::ABC ? matrixNodeIndex(1, channel)
                                                   : matrixNodeIndex(1);
    if (terminalNotGrounded(1))
      (**mIntfVoltage)(channel, 0) =
          Math::realFromVectorElement(solution, row1);
    if (terminalNotGrounded(0))
      (**mIntfVoltage)(channel, 0) -=
          Math::realFromVectorElement(solution, row0);
  }
}

void CoSimulationSSN::mnaCompUpdateVoltage(const Matrix &solution) {
  updateTerminalVoltage(solution);
}

void CoSimulationSSN::mnaCompPreStep(Real time, Int) {
  if (!std::isfinite(time) || time < -1e-12)
    throw std::invalid_argument("Co-Simulation received an invalid DPsim time");
  const double communicationTime = std::max(0.0, time);
  const double inputTime =
      time > 1e-12 ? std::max(0.0, time - mTimeStep) : communicationTime;
  Matrix values = **mInputs;
  for (std::size_t index = 0; index < mInputFunctions.size(); ++index) {
    if (mControls[index].getPtr()) {
      if (mInputFunctions[index]) {
        const double value = mInputFunctions[index](inputTime);
        if (!std::isfinite(value))
          throw std::runtime_error(
              "Co-Simulation control returned non-finite input");
        values(index, 0) = value;
      } else {
        values(index, 0) = **mControls[index];
      }
    }
  }
  for (std::size_t phase = 0; phase < mVoltageIndices.size(); ++phase)
    values(mVoltageIndices[phase], 0) = (**mIntfVoltage)(phase, 0);

  // The first DPsim solve is the initialized t=0 point. Each later pre-step
  // advances the FMU over the preceding interval with the last accepted PCC
  // voltage held constant; the completed output current is stamped now.
  if (time > 1e-12) {
    const double startTime = time - mTimeStep;
    if (std::abs(startTime - mModel->time()) >
        1e-10 * std::max(1.0, std::abs(startTime)))
      throw std::runtime_error(
          "Co-Simulation time grid diverged from the DPsim communication grid");
    mModel->setInputs(values);
    mModel->doStep(startTime, mTimeStep);
    **mOutputs = mModel->outputs();
    **mInputs = values;
    for (std::size_t phase = 0; phase < mCurrentIndices.size(); ++phase)
      (**mIntfCurrent)(phase, 0) =
          mSign * (**mOutputs)(mCurrentIndices[phase], 0);
    **mStepAttribute += 1;
    **mInternalStepAttribute = 1;
  } else {
    mModel->setInputs(values);
    **mInputs = values;
    **mInternalStepAttribute = 0;
  }
  **mRightVector = Matrix::Zero((**mRightVector).rows(), 1);
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

void CoSimulationSSN::mnaCompPostStep(Real, Int,
                                      Attribute<Matrix>::Ptr &solution) {
  updateTerminalVoltage(**solution);
}

void CoSimulationSSN::mnaCompAddPreStepDependencies(
    AttributeBase::List &previous, AttributeBase::List &dependencies,
    AttributeBase::List &modified) {
  previous.push_back(mIntfVoltage);
  for (const auto &control : mControls)
    if (control.getPtr())
      dependencies.push_back(control);
  modified.push_back(mRightVector);
  modified.push_back(mInputs);
  modified.push_back(mOutputs);
  modified.push_back(mIntfCurrent);
  modified.push_back(mStepAttribute);
  modified.push_back(mInternalStepAttribute);
}

void CoSimulationSSN::mnaCompAddPostStepDependencies(
    AttributeBase::List &, AttributeBase::List &dependencies,
    AttributeBase::List &modified, Attribute<Matrix>::Ptr &solution) {
  dependencies.push_back(solution);
  modified.push_back(mIntfVoltage);
}

} // namespace DPsim::FMI
