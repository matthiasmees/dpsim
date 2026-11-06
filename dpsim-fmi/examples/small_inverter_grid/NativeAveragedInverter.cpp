// SPDX-License-Identifier: MPL-2.0
#include "NativeAveragedInverter.h"

#include <algorithm>
#include <cmath>
#include <dpsim-models/MNAStampUtils.h>
#include <stdexcept>

using namespace CPS;
namespace SmallGrid {
namespace {
constexpr double kActivePowerW = 125e6;
constexpr double kReactivePowerVar = 80.881e6;
constexpr double kNominalLineVoltageV = 345e3;
constexpr double kCurrentLimitPu = 1.20;
} // namespace

NativeAveragedInverter::NativeAveragedInverter(String name)
    : MNASimPowerComp<Real>(name, true, true) {
  mPhaseType = PhaseType::ABC;
  setTerminalNumber(2);
  mState = mAttributes->create<Matrix>("current_state");
  mIterationAttribute = mAttributes->create<Real>("newton_iterations");
  mResidualAttribute = mAttributes->create<Real>("newton_residual");
  mInitialVoltage = Matrix::Zero(3, 1);
  mOldState = Matrix::Zero(3, 1);
  mOldTarget = Matrix::Zero(3, 1);
  mTrialState = Matrix::Zero(3, 1);
  mTrialVoltage = Matrix::Zero(3, 1);
  mW = Matrix::Zero(3, 3);
  mHistory = Matrix::Zero(3, 1);
  **mState = Matrix::Zero(3, 1);
  **mIntfVoltage = Matrix::Zero(3, 1);
  **mIntfCurrent = Matrix::Zero(3, 1);
  **mIterationAttribute = 0;
  **mResidualAttribute = 0;
}

SimPowerComp<Real>::Ptr NativeAveragedInverter::clone(String name) {
  return NativeAveragedInverter::make(std::move(name));
}

void NativeAveragedInverter::connect(SimNode<Real>::Ptr node) {
  if (!node)
    throw std::invalid_argument("Cannot connect native inverter to null node");
  SimPowerComp<Real>::connect({SimNode<Real>::GND, std::move(node)});
}

void NativeAveragedInverter::validateTerminals() {
  for (UInt terminalIndex = 0; terminalIndex < 2; ++terminalIndex) {
    const auto connectedNode = node(terminalIndex);
    if (!connectedNode || (!connectedNode->isGround() &&
                           connectedNode->phaseType() != PhaseType::ABC))
      throw std::invalid_argument(
          "Native averaged inverter requires ABC terminal nodes or ground");
  }
}

Matrix NativeAveragedInverter::passiveCurrentReference(const Matrix &voltage) {
  if (voltage.rows() != 3 || voltage.cols() != 1 || !voltage.allFinite())
    throw std::invalid_argument("Native inverter voltage must be finite abc");
  constexpr double root3over2 = 0.86602540378443864676;
  const double alpha =
      (2.0 / 3.0) * (voltage(0, 0) - 0.5 * voltage(1, 0) - 0.5 * voltage(2, 0));
  const double beta =
      (2.0 / 3.0) * root3over2 * (voltage(1, 0) - voltage(2, 0));
  const double voltageSquared = alpha * alpha + beta * beta;
  Matrix passive = Matrix::Zero(3, 1);
  if (!std::isfinite(voltageSquared) || voltageSquared < 1.0)
    return passive;

  const double factor = 2.0 / (3.0 * voltageSquared);
  double injectionAlpha =
      factor * (kActivePowerW * alpha + kReactivePowerVar * beta);
  double injectionBeta =
      factor * (kActivePowerW * beta - kReactivePowerVar * alpha);
  const double ratedPeak = std::sqrt(2.0) *
                           std::hypot(kActivePowerW, kReactivePowerVar) /
                           (std::sqrt(3.0) * kNominalLineVoltageV);
  const double maximum = kCurrentLimitPu * ratedPeak;
  const double magnitude = std::hypot(injectionAlpha, injectionBeta);
  if (magnitude > maximum) {
    const double scale = maximum / magnitude;
    injectionAlpha *= scale;
    injectionBeta *= scale;
  }
  const double ia = injectionAlpha;
  const double ib = -0.5 * injectionAlpha + root3over2 * injectionBeta;
  const double ic = -0.5 * injectionAlpha - root3over2 * injectionBeta;
  // The FMU exposes passive-sign current, opposite to positive power injection.
  passive << -ia, -ib, -ic;
  return passive;
}

Matrix NativeAveragedInverter::endpointCurrent(const Matrix &voltage) const {
  const double denominator = 2.0 * kTauSeconds + mStep;
  const double oldWeight = (2.0 * kTauSeconds - mStep) / denominator;
  const double targetWeight = mStep / denominator;
  return oldWeight * mOldState +
         targetWeight * (mOldTarget + passiveCurrentReference(voltage));
}

Matrix NativeAveragedInverter::voltageFromSolution(const Matrix &solution) {
  Matrix voltage = Matrix::Zero(3, 1);
  if (terminalNotGrounded(1))
    for (UInt phase = 0; phase < 3; ++phase)
      voltage(phase, 0) +=
          Math::realFromVectorElement(solution, matrixNodeIndex(1, phase));
  if (terminalNotGrounded(0))
    for (UInt phase = 0; phase < 3; ++phase)
      voltage(phase, 0) -=
          Math::realFromVectorElement(solution, matrixNodeIndex(0, phase));
  return voltage;
}

void NativeAveragedInverter::initializeFromNodesAndTerminals(Real) {
  validateTerminals();
  MatrixComp terminalVoltage = MatrixComp::Zero(3, 1);
  if (!node(1)->isGround())
    terminalVoltage += initialVoltage(1);
  if (!node(0)->isGround())
    terminalVoltage -= initialVoltage(0);
  if (terminalVoltage.allFinite())
    mInitialVoltage = (RMS3PH_TO_PEAK1PH * terminalVoltage).real();
  const Matrix initialCurrent = passiveCurrentReference(mInitialVoltage);
  **mState = initialCurrent;
  **mIntfVoltage = mInitialVoltage;
  **mIntfCurrent = initialCurrent;
}

void NativeAveragedInverter::initialize(Matrix frequencies) {
  SimPowerComp<Real>::initialize(std::move(frequencies));
  initializeFromNodesAndTerminals(0);
}

void NativeAveragedInverter::mnaCompInitialize(Real, Real timeStep,
                                               Attribute<Matrix>::Ptr) {
  validateTerminals();
  if (!std::isfinite(timeStep) || timeStep <= 0)
    throw std::invalid_argument("Invalid native inverter timestep");
  mStep = timeStep;
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
  mOldTarget = passiveCurrentReference(**mIntfVoltage);
  mTrialVoltage = **mIntfVoltage;
  mTrialState = mOldState;
  prepareTangent();
}

void NativeAveragedInverter::mnaCompPreStep(Real, Int) {
  mIterations = 0;
  mResidual = 0;
  mOldState = **mState;
  mTrialVoltage = **mIntfVoltage;
  mOldTarget = passiveCurrentReference(**mIntfVoltage);
  mTrialState = endpointCurrent(mTrialVoltage);
  prepareTangent();
}

void NativeAveragedInverter::prepareTangent() {
  for (UInt phase = 0; phase < 3; ++phase) {
    const double perturbation =
        kJacobianPerturbation *
        std::max(1.0, std::abs(mTrialVoltage(phase, 0)));
    Matrix plus = mTrialVoltage, minus = mTrialVoltage;
    plus(phase, 0) += perturbation;
    minus(phase, 0) -= perturbation;
    mW.col(phase) =
        (endpointCurrent(plus) - endpointCurrent(minus)) / (2.0 * perturbation);
  }
  mHistory = mTrialState - mW * mTrialVoltage;
  if (!mW.allFinite() || !mHistory.allFinite())
    throw std::runtime_error("Nonfinite native inverter Norton stamp");
  **mRightVector = Matrix::Zero((**mRightVector).rows(), 1);
  mnaCompApplyRightSideVectorStamp(**mRightVector);
}

bool NativeAveragedInverter::mnaCheckNonlinearConvergence(
    const Matrix &solution) {
  if (++mIterations > kMaximumNewtonIterations)
    throw std::runtime_error("Native inverter Newton convergence limit");
  mTrialVoltage = voltageFromSolution(solution);
  mTrialState = endpointCurrent(mTrialVoltage);
  const Matrix predicted = mW * mTrialVoltage + mHistory;
  mResidual = 0;
  for (UInt phase = 0; phase < 3; ++phase) {
    const double scale = kAbsoluteCurrentTolerance +
                         kRelativeCurrentTolerance *
                             std::max({1.0, std::abs(mTrialState(phase, 0)),
                                       std::abs(predicted(phase, 0))});
    mResidual =
        std::max(mResidual,
                 std::abs(mTrialState(phase, 0) - predicted(phase, 0)) / scale);
  }
  return mResidual <= 1.0;
}

void NativeAveragedInverter::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &matrix) {
  MNAStampUtils::stampConductanceMatrix(
      mW, matrix, matrixNodeIndex(0), matrixNodeIndex(1),
      terminalNotGrounded(0), terminalNotGrounded(1), mSLog);
}

void NativeAveragedInverter::mnaCompApplyRightSideVectorStamp(Matrix &vector) {
  for (UInt phase = 0; phase < 3; ++phase) {
    if (terminalNotGrounded(0))
      Math::setVectorElement(vector, matrixNodeIndex(0, phase),
                             mHistory(phase, 0));
    if (terminalNotGrounded(1))
      Math::setVectorElement(vector, matrixNodeIndex(1, phase),
                             -mHistory(phase, 0));
  }
}

void NativeAveragedInverter::mnaCompUpdateVoltage(const Matrix &solution) {
  **mIntfVoltage = voltageFromSolution(solution);
}

void NativeAveragedInverter::mnaCompUpdateCurrent(const Matrix &) {
  **mIntfCurrent = mTrialState;
}

void NativeAveragedInverter::mnaCompPostStep(Real, Int,
                                             Attribute<Matrix>::Ptr &) {
  if (mIterations == 0 || mResidual > 1.0)
    throw std::runtime_error(
        "Native inverter did not converge before post-step");
  **mState = mTrialState;
  **mIntfVoltage = mTrialVoltage;
  **mIntfCurrent = mTrialState;
  **mIterationAttribute = mIterations;
  **mResidualAttribute = mResidual;
}

void NativeAveragedInverter::mnaCompAddPreStepDependencies(
    AttributeBase::List &previous, AttributeBase::List &dependencies,
    AttributeBase::List &modified) {
  previous.push_back(mState);
  previous.push_back(mIntfVoltage);
  modified.push_back(mRightVector);
}

void NativeAveragedInverter::mnaCompAddPostStepDependencies(
    AttributeBase::List &, AttributeBase::List &dependencies,
    AttributeBase::List &modified, Attribute<Matrix>::Ptr &solution) {
  dependencies.push_back(solution);
  for (auto &attribute : {mState, mIntfVoltage, mIntfCurrent})
    modified.push_back(attribute);
  modified.push_back(mIterationAttribute);
  modified.push_back(mResidualAttribute);
}
} // namespace SmallGrid
