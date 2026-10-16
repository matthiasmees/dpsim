// SPDX-License-Identifier: MPL-2.0
#include <algorithm>
#include <cmath>
#include <dpsim-fmi/NonlinearSSN.h>
#include <limits>
#include <set>
#include <stdexcept>

using namespace CPS;
namespace DPsim::FMI {
namespace {
std::size_t indexOf(const std::vector<std::string> &names,
                    const std::string &name) {
  auto found = std::find(names.begin(), names.end(), name);
  if (found == names.end())
    throw std::invalid_argument("FMI variable does not exist: " + name);
  return found - names.begin();
}
} // namespace
NonlinearSSN::NonlinearSSN(String name, ModelExchange::Ptr model,
                           const ElectricalPort3Ph &port,
                           NewtonSettings settings)
    : MNASimPowerComp(name, true, true), mModel(std::move(model)),
      mSettings(settings), mSign(port.currentSignToPassive) {
  if (!mModel)
    throw std::invalid_argument(
        "NonlinearSSN requires a Model Exchange instance");
  if (mSign != 1 && mSign != -1)
    throw std::invalid_argument("Current sign must be +1 or -1");
  for (double t :
       {settings.absoluteStateTolerance, settings.relativeStateTolerance,
        settings.absoluteCurrentTolerance, settings.relativeCurrentTolerance,
        settings.relativeJacobianPerturbation})
    if (!std::isfinite(t) || t <= 0)
      throw std::invalid_argument(
          "Newton tolerances/perturbation must be finite and positive");
  if (settings.maxIterations < 1 || settings.maxIterations > 40)
    throw std::invalid_argument(
        "Newton iteration limit must be between 1 and 40");
  mPhaseType = PhaseType::ABC;
  setTerminalNumber(2);
  const auto &d = mModel->description();
  std::set<std::size_t> uniqueVoltage, uniqueCurrent;
  mInitialVoltage.resize(3, 1);
  for (int p = 0; p < 3; ++p) {
    mVoltageIndices[p] = indexOf(d.inputNames, port.voltageInputs[p]);
    mCurrentRows[p] = indexOf(d.outputNames, port.currentOutputs[p]);
    uniqueVoltage.insert(mVoltageIndices[p]);
    uniqueCurrent.insert(mCurrentRows[p]);
    mInitialVoltage(p, 0) = d.initialInputs(mVoltageIndices[p], 0);
  }
  if (uniqueVoltage.size() != 3 || uniqueCurrent.size() != 3)
    throw std::invalid_argument(
        "Electrical abc mappings must contain distinct variable names");
  mState = mAttributes->create<Matrix>("x");
  mInputs = mAttributes->create<Matrix>("fmi.inputs");
  mOutputs = mAttributes->create<Matrix>("fmi.outputs");
  mIterationAttribute = mAttributes->create<Real>("newton_iterations");
  mResidualAttribute = mAttributes->create<Real>("newton_residual");
  mControls.resize(d.inputNames.size());
  mInputFunctions.resize(d.inputNames.size());
  mInputEvents.resize(d.inputNames.size());
  for (std::size_t k = 0; k < d.inputNames.size(); ++k) {
    auto a = mAttributes->createDynamic<Real>("fmi.input." + d.inputNames[k]);
    a->setReference(mInputs->deriveCoeff<Real>(k, 0));
    if (!uniqueVoltage.count(k)) {
      mControls[k] = mAttributes->create<Real>("control." + d.inputNames[k]);
      **mControls[k] = d.initialInputs(k, 0);
    }
  }
  for (std::size_t k = 0; k < d.outputNames.size(); ++k) {
    auto a = mAttributes->createDynamic<Real>("fmi.output." + d.outputNames[k]);
    a->setReference(mOutputs->deriveCoeff<Real>(k, 0));
  }
  for (std::size_t k = 0; k < d.stateNames.size(); ++k) {
    auto a = mAttributes->createDynamic<Real>("fmi.state." + d.stateNames[k]);
    a->setReference(mState->deriveCoeff<Real>(k, 0));
  }
  mW = Matrix::Zero(3, 3);
  mHistory = Matrix::Zero(3, 1);
  initializeValues();
}
void NonlinearSSN::initializeValues() {
  **mState = mModel->description().initialState;
  **mInputs = mModel->description().initialInputs;
  **mIntfVoltage = mInitialVoltage;
  auto e = mModel->evaluate(**mState, **mInputs, 0);
  **mOutputs = e.second;
  **mIntfCurrent = currents(e.second);
  **mIterationAttribute = 0;
  **mResidualAttribute = 0;
}
void NonlinearSSN::initialize(Matrix frequencies) {
  SimPowerComp<Real>::initialize(std::move(frequencies));
  initializeValues();
}
void NonlinearSSN::validateTerminals() {
  for (UInt p = 0; p < 2; ++p) {
    const auto n = node(p);
    if (!n || (!n->isGround() && n->phaseType() != PhaseType::ABC))
      throw std::invalid_argument(
          "NonlinearSSN requires ABC terminal nodes or ground");
  }
}
void NonlinearSSN::initializeFromNodesAndTerminals(Real) {
  validateTerminals();
  initializeValues();
}
SimPowerComp<Real>::Ptr NonlinearSSN::clone(String) {
  throw std::logic_error(
      "A persistent nonlinear FMU cannot be shared by clones; load a separate "
      "instance if its capability permits it");
}
void NonlinearSSN::setInputFunction(const std::string &name,
                                    std::function<double(double)> function,
                                    std::vector<double> events) {
  const auto k = indexOf(mModel->description().inputNames, name);
  if (!mControls[k].getPtr() || !function)
    throw std::invalid_argument(
        "Input function requires an external control input");
  const double value = function(0);
  if (!std::isfinite(value) ||
      std::abs(value - mModel->description().initialInputs(k, 0)) >
          1e-12 * (1 + std::abs(value)))
    throw std::invalid_argument(
        "Input function must match FMI initialization at zero");
  std::sort(events.begin(), events.end());
  for (double t : events)
    if (!std::isfinite(t) || t <= 0)
      throw std::invalid_argument("Input events must be finite positive times");
  if (std::adjacent_find(events.begin(), events.end()) != events.end())
    throw std::invalid_argument("Duplicate input event");
  mInputFunctions[k] = std::move(function);
  mInputEvents[k] = std::move(events);
}
Matrix NonlinearSSN::inputsAt(const Matrix &voltage, double time,
                              bool left) const {
  Matrix values = mModel->description().initialInputs;
  for (std::size_t k = 0; k < mControls.size(); ++k)
    if (mControls[k].getPtr()) {
      if (mInputFunctions[k]) {
        double sample = time;
        for (double e : mInputEvents[k])
          if (std::abs(time - e) < 1e-10) {
            sample = left ? std::nextafter(
                                e, -std::numeric_limits<double>::infinity())
                          : e;
            break;
          }
        values(k, 0) = mInputFunctions[k](sample);
      } else
        values(k, 0) = **mControls[k];
    }
  for (int p = 0; p < 3; ++p)
    values(mVoltageIndices[p], 0) = voltage(p, 0);
  return values;
}
Matrix NonlinearSSN::currents(const Matrix &outputs) const {
  Matrix current(3, 1);
  for (int p = 0; p < 3; ++p)
    current(p, 0) = mSign * outputs(mCurrentRows[p], 0);
  return current;
}
Matrix NonlinearSSN::voltageFromSolution(const Matrix &solution) {
  Matrix v = Matrix::Zero(3, 1);
  for (UInt p = 0; p < 3; ++p) {
    if (terminalNotGrounded(1))
      v(p, 0) += Math::realFromVectorElement(solution, matrixNodeIndex(1, p));
    if (terminalNotGrounded(0))
      v(p, 0) -= Math::realFromVectorElement(solution, matrixNodeIndex(0, p));
  }
  return v;
}
bool NonlinearSSN::eventAt(double time) const {
  for (const auto &events : mInputEvents)
    for (double e : events)
      if (std::abs(e - time) < 1e-10)
        return true;
  return false;
}
void NonlinearSSN::mnaCompInitialize(Real, Real step, Attribute<Matrix>::Ptr) {
  validateTerminals();
  mStep = step;
  if (!std::isfinite(step) || step <= 0)
    throw std::invalid_argument("Invalid nonlinear integration step");
  for (const auto &events : mInputEvents)
    for (double e : events)
      if (std::abs(e / step - std::round(e / step)) > 1e-7)
        throw std::invalid_argument(
            "Scheduled input events must align with the fixed timestep");
  updateMatrixNodeIndices();
  mVariableSystemMatrixEntries.clear();
  for (UInt a = 0; a < 2; ++a)
    for (UInt b = 0; b < 2; ++b)
      if (terminalNotGrounded(a) && terminalNotGrounded(b))
        for (UInt i = 0; i < 3; ++i)
          for (UInt j = 0; j < 3; ++j)
            mVariableSystemMatrixEntries.emplace_back(matrixNodeIndex(a, i),
                                                      matrixNodeIndex(b, j));
  mOldState = **mState;
  mOldDerivative = mModel->evaluate(mOldState, **mInputs, 0).first;
  mTrialState = mOldState;
  mTrialVoltage = **mIntfVoltage;
  mNewInputs = **mInputs;
  auto e = mModel->evaluate(mTrialState, mNewInputs, 0);
  mTrialDerivative = e.first;
  mTrialOutputs = e.second;
  // Initial stamp is only used to establish a complete sparse pattern.
  mTime = 0;
  mnaPrepareNonlinearIteration();
}
void NonlinearSSN::mnaCompPreStep(Real time, Int) {
  mTime = time;
  const double event = mModel->nextTimeEvent();
  if (std::isfinite(event) && event <= time + 1e-10 &&
      std::abs(event - time) > 1e-10)
    throw std::runtime_error(
        "Scheduled FMI time event must align with the fixed timestep");
  mIterations = 0;
  mOldState = **mState;
  const Matrix oldInputs = inputsAt(**mIntfVoltage, time - mStep, false);
  if ((oldInputs - **mInputs).norm() > 1e-12)
    mModel->inputEvent(mOldState, oldInputs, time - mStep);
  mOldDerivative = mModel->evaluate(mOldState, oldInputs, time - mStep).first;
  mTrialState = mOldState + mStep * mOldDerivative;
  mTrialVoltage = **mIntfVoltage;
  mNewInputs = inputsAt(mTrialVoltage, time, true);
  auto e = mModel->evaluate(mTrialState, mNewInputs, time);
  mTrialDerivative = e.first;
  mTrialOutputs = e.second;
  mnaPrepareNonlinearIteration();
}
void NonlinearSSN::mnaPrepareNonlinearIteration() {
  const Eigen::Index n = mTrialState.rows();
  Matrix a(n, n), b(n, 3), c(3, n), d(3, 3);
  for (Eigen::Index j = 0; j < n; ++j) {
    const double eps = mSettings.relativeJacobianPerturbation *
                       std::max(std::abs(mTrialState(j, 0)),
                                mModel->description().stateNominals(j, 0));
    Matrix plus = mTrialState, minus = mTrialState;
    plus(j, 0) += eps;
    minus(j, 0) -= eps;
    const auto p = mModel->evaluate(plus, mNewInputs, mTime),
               m = mModel->evaluate(minus, mNewInputs, mTime);
    a.col(j) = (p.first - m.first) / (2 * eps);
    c.col(j) = (currents(p.second) - currents(m.second)) / (2 * eps);
  }
  for (int j = 0; j < 3; ++j) {
    const double eps = mSettings.relativeJacobianPerturbation *
                       std::max(1.0, std::abs(mTrialVoltage(j, 0)));
    Matrix plus = mNewInputs, minus = mNewInputs;
    plus(mVoltageIndices[j], 0) += eps;
    minus(mVoltageIndices[j], 0) -= eps;
    const auto p = mModel->evaluate(mTrialState, plus, mTime),
               m = mModel->evaluate(mTrialState, minus, mTime);
    b.col(j) = (p.first - m.first) / (2 * eps);
    d.col(j) = (currents(p.second) - currents(m.second)) / (2 * eps);
  }
  const Matrix residual =
      mTrialState - mOldState - mStep * (mOldDerivative + mTrialDerivative) / 2;
  const Matrix left = Matrix::Identity(n, n) - mStep * a / 2;
  if (n) {
    const auto lu = left.fullPivLu();
    if (!lu.isInvertible())
      throw std::runtime_error("Singular FMI implicit state Jacobian");
    mCorrection = lu.solve(-residual);
    mSensitivity = lu.solve(mStep * b / 2);
  } else {
    mCorrection = Matrix::Zero(0, 1);
    mSensitivity = Matrix::Zero(0, 3);
  }
  mW = d + c * mSensitivity;
  mHistory = currents(mTrialOutputs) + c * mCorrection - mW * mTrialVoltage;
  if (!mW.allFinite() || !mHistory.allFinite())
    throw std::runtime_error("Nonfinite nonlinear Norton stamp");
  mLinearizationState = mTrialState;
  mLinearizationVoltage = mTrialVoltage;
  **mRightVector = Matrix::Zero((**mRightVector).rows(), 1);
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}
bool NonlinearSSN::mnaCheckNonlinearConvergence(const Matrix &solution) {
  ++mIterations;
  if (mIterations > mSettings.maxIterations)
    throw std::runtime_error("FMI Newton convergence limit at t=" +
                             std::to_string(mTime));
  mTrialVoltage = voltageFromSolution(solution);
  mTrialState = mLinearizationState + mCorrection +
                mSensitivity * (mTrialVoltage - mLinearizationVoltage);
  mNewInputs = inputsAt(mTrialVoltage, mTime, true);
  auto e = mModel->evaluate(mTrialState, mNewInputs, mTime);
  mTrialDerivative = e.first;
  mTrialOutputs = e.second;
  const Matrix r =
      mTrialState - mOldState - mStep * (mOldDerivative + mTrialDerivative) / 2;
  mStateResidual = 0;
  for (Eigen::Index j = 0; j < r.rows(); ++j) {
    const double scale =
        mSettings.absoluteStateTolerance +
        mSettings.relativeStateTolerance *
            std::max({std::abs(mOldState(j, 0)), std::abs(mTrialState(j, 0)),
                      mModel->description().stateNominals(j, 0)});
    mStateResidual = std::max(mStateResidual, std::abs(r(j, 0)) / scale);
  }
  const Matrix actual = currents(e.second),
               predicted = mW * mTrialVoltage + mHistory;
  mCurrentResidual = 0;
  for (int j = 0; j < 3; ++j) {
    const double scale =
        mSettings.absoluteCurrentTolerance +
        mSettings.relativeCurrentTolerance *
            std::max({1.0, std::abs(actual(j, 0)), std::abs(predicted(j, 0))});
    mCurrentResidual = std::max(
        mCurrentResidual, std::abs(actual(j, 0) - predicted(j, 0)) / scale);
  }
  return mStateResidual <= 1 && mCurrentResidual <= 1;
}
void NonlinearSSN::mnaCompApplySystemMatrixStamp(SparseMatrixRow &matrix) {
  MNAStampUtils::stampConductanceMatrix(
      mW, matrix, matrixNodeIndex(0), matrixNodeIndex(1),
      terminalNotGrounded(0), terminalNotGrounded(1), mSLog);
}
void NonlinearSSN::mnaCompApplyRightSideVectorStamp(Matrix &vector) {
  for (UInt j = 0; j < 3; ++j) {
    if (terminalNotGrounded(0))
      Math::setVectorElement(vector, matrixNodeIndex(0, j), mHistory(j, 0));
    if (terminalNotGrounded(1))
      Math::setVectorElement(vector, matrixNodeIndex(1, j), -mHistory(j, 0));
  }
}
void NonlinearSSN::mnaCompPostStep(Real time, Int, Attribute<Matrix>::Ptr &) {
  if (mIterations == 0 || mStateResidual > 1 || mCurrentResidual > 1)
    throw std::runtime_error("MNA solver did not converge the nonlinear FMI "
                             "component before post-step");
  // Select the converged endpoint without another provisional getter.
  mModel->prepareStep(mTrialState, mNewInputs, time);
  const Matrix rightInputs = inputsAt(mTrialVoltage, time, false);
  if (eventAt(time))
    mModel->completeStep(rightInputs);
  else
    mModel->completeStep();
  auto e = mModel->acceptedValues();
  const Matrix acceptedCurrent = currents(mTrialOutputs);
  const Matrix eventCurrent = currents(e.second);
  for (int p = 0; p < 3; ++p) {
    const double scale = mSettings.absoluteCurrentTolerance +
                         mSettings.relativeCurrentTolerance *
                             std::max({1.0, std::abs(acceptedCurrent(p, 0)),
                                       std::abs(eventCurrent(p, 0))});
    if (std::abs(acceptedCurrent(p, 0) - eventCurrent(p, 0)) > scale)
      throw std::runtime_error(
          "FMI event changed terminal current; instantaneous algebraic event "
          "iteration is unsupported");
  }
  **mState = mModel->acceptedState();
  **mInputs = rightInputs;
  **mOutputs = e.second;
  **mIntfVoltage = mTrialVoltage;
  **mIntfCurrent = currents(e.second);
  **mIterationAttribute = mIterations;
  **mResidualAttribute = std::max(mStateResidual, mCurrentResidual);
}
void NonlinearSSN::mnaCompAddPreStepDependencies(
    AttributeBase::List &previous, AttributeBase::List &dependencies,
    AttributeBase::List &modified) {
  previous.push_back(mState);
  previous.push_back(mInputs);
  previous.push_back(mIntfVoltage);
  for (auto &control : mControls)
    if (control.getPtr())
      dependencies.push_back(control);
  modified.push_back(mRightVector);
}
void NonlinearSSN::mnaCompAddPostStepDependencies(
    AttributeBase::List &, AttributeBase::List &dependencies,
    AttributeBase::List &modified, Attribute<Matrix>::Ptr &solution) {
  dependencies.push_back(solution);
  for (auto &a : {mState, mInputs, mOutputs, mIntfVoltage, mIntfCurrent})
    modified.push_back(a);
  modified.push_back(mIterationAttribute);
  modified.push_back(mResidualAttribute);
}
} // namespace DPsim::FMI
