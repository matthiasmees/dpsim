// SPDX-License-Identifier: MPL-2.0
#include <dpsim-fmi/LinearSSN.h>
#include <stdexcept>

using namespace CPS;
namespace DPsim::FMI {
namespace {
class OutputProjection final : public AttributeUpdateTaskBase<Real> {
  Matrix mC;
  Matrix mD;
  Attribute<Matrix>::Ptr mState, mVoltage;

public:
  OutputProjection(Matrix c, Matrix d, Attribute<Matrix>::Ptr state,
                   Attribute<Matrix>::Ptr voltage)
      : mC(std::move(c)), mD(std::move(d)), mState(state), mVoltage(voltage) {}
  void executeUpdate(std::shared_ptr<Real> &value) override {
    *value = (mC * (**mState) + mD * (**mVoltage))(0, 0);
  }
  AttributeBase::List getDependencies() override { return {mState, mVoltage}; }
};
} // namespace

LinearSSN::LinearSSN(String name, LinearModel::Ptr model, Logger::Level level)
    : VTypeSSNComp(name, name, model ? model->voltageInputs.size() : 1,
                   model ? model->voltageInputs.size() : 1, level),
      mModel(std::move(model)) {
  if (!mModel)
    throw std::invalid_argument("FMI SSN requires an imported model");
  if (!mModel->A.allFinite() || !mModel->B.allFinite() ||
      !mModel->C.allFinite() || !mModel->D.allFinite() ||
      !mModel->initialState.allFinite() || !mModel->initialVoltage.allFinite())
    throw std::invalid_argument("FMI SSN model contains non-finite values");
  const auto n = mModel->A.rows();
  const auto p = static_cast<Eigen::Index>(mModel->outputNames.size());
  const auto channels = static_cast<Eigen::Index>(mModel->voltageInputs.size());
  if ((channels != 1 && channels != 3) ||
      mModel->currentOutputs.size() != static_cast<std::size_t>(channels))
    throw std::invalid_argument(
        "FMI SSN requires one or three voltage/current channels");
  for (const auto row : mModel->currentOutputs)
    if (row >= mModel->outputNames.size())
      throw std::invalid_argument("FMI SSN current output index is invalid");
  if (mModel->A.cols() != n || mModel->B.rows() != n ||
      mModel->B.cols() != channels || mModel->C.rows() != p ||
      mModel->C.cols() != n || mModel->D.rows() != p ||
      mModel->D.cols() != channels || mModel->initialState.rows() != n ||
      mModel->initialState.cols() != 1 ||
      mModel->initialVoltage.rows() != channels ||
      mModel->initialVoltage.cols() != 1 ||
      mModel->stateNames.size() != static_cast<std::size_t>(n))
    throw std::invalid_argument("FMI SSN model dimensions are inconsistent");
  mPhaseType = channels == 3 ? PhaseType::ABC : PhaseType::Single;
  setTerminalNumber(2);
  Matrix currentC(channels, n), currentD(channels, channels);
  for (Eigen::Index channel = 0; channel < channels; ++channel) {
    currentC.row(channel) = mModel->C.row(mModel->currentOutputs[channel]);
    currentD.row(channel) = mModel->D.row(mModel->currentOutputs[channel]);
  }
  setParameters(mModel->A, mModel->B, currentC, currentD);
  **mX = mModel->initialState;
  **mIntfVoltage = mModel->initialVoltage;
  **mIntfCurrent = mC * (**mX) + mD * (**mIntfVoltage);

  // Derived attributes explicitly depend on the state and terminal-voltage
  // attributes, so loggers/consumers are scheduled after the SSN post-step.
  for (Eigen::Index row = 0; row < p; ++row) {
    auto output = mAttributes->createDynamic<Real>("fmi.output." +
                                                   mModel->outputNames[row]);
    auto dynamic =
        std::dynamic_pointer_cast<AttributeDynamic<Real>>(output.getPtr());
    dynamic->addTask(UpdateTaskKind::UPDATE_ON_GET,
                     std::make_shared<OutputProjection>(mModel->C.row(row),
                                                        mModel->D.row(row), mX,
                                                        mIntfVoltage));
  }
  for (Eigen::Index channel = 0; channel < channels; ++channel) {
    auto input = mAttributes->createDynamic<Real>(
        "fmi.input." + mModel->voltageInputs[channel]);
    input->setReference(mIntfVoltage->deriveCoeff<Real>(channel, 0));
  }
  for (Eigen::Index index = 0; index < n; ++index) {
    auto state = mAttributes->createDynamic<Real>("fmi.state." +
                                                  mModel->stateNames[index]);
    state->setReference(mX->deriveCoeff<Real>(index, 0));
  }
}

SimPowerComp<Real>::Ptr LinearSSN::clone(String name) {
  return LinearSSN::make(name, mModel, mLogLevel);
}
std::vector<String> LinearSSN::getLocalStateNames() const {
  return mModel->stateNames;
}
void LinearSSN::validateTerminals() const {
  for (UInt index = 0; index < 2; ++index) {
    const auto n = const_cast<LinearSSN *>(this)->node(index);
    if (!n)
      throw std::invalid_argument("FMI LinearSSN terminal node has expired");
    if (!n->isGround() && n->phaseType() != mPhaseType)
      throw std::invalid_argument(
          "FMI LinearSSN terminal phase type does not match its channel count");
  }
}
MatrixComp LinearSSN::buildInitialInputFromNodes(Real) {
  validateTerminals();
  // FMI initialization and the SSN history must use the same instantaneous
  // voltage. A phasor-based steady-state guess must not replace this value.
  return mModel->initialVoltage.cast<Complex>();
}
void LinearSSN::initialize(Matrix frequencies) {
  SimPowerComp<Real>::initialize(std::move(frequencies));
  **mX = mModel->initialState;
  **mIntfVoltage = mModel->initialVoltage;
  **mIntfCurrent = mC * (**mX) + mD * (**mIntfVoltage);
}
void LinearSSN::initializeFromNodesAndTerminals(Real frequency) {
  **mX = mModel->initialState;
  **mIntfVoltage = buildInitialInputFromNodes(frequency).real();
  **mIntfCurrent = mC * (**mX) + mD * (**mIntfVoltage);
}
void LinearSSN::mnaCompApplySystemMatrixStamp(SparseMatrixRow &matrix) {
  validateTerminals();
  if (mPhaseType == PhaseType::Single)
    MNAStampUtils::stampConductance(mW(0, 0), matrix, matrixNodeIndex(0),
                                    matrixNodeIndex(1), terminalNotGrounded(0),
                                    terminalNotGrounded(1), mSLog);
  else
    MNAStampUtils::stampConductanceMatrix(
        mW, matrix, matrixNodeIndex(0), matrixNodeIndex(1),
        terminalNotGrounded(0), terminalNotGrounded(1), mSLog);
}
void LinearSSN::mnaCompApplyRightSideVectorStamp(Matrix &vector) {
  for (Eigen::Index channel = 0; channel < mW.rows(); ++channel) {
    const auto phase = static_cast<UInt>(channel);
    if (terminalNotGrounded(0))
      Math::setVectorElement(vector, matrixNodeIndex(0, phase),
                             mYHist(channel, 0));
    if (terminalNotGrounded(1))
      Math::setVectorElement(vector, matrixNodeIndex(1, phase),
                             -mYHist(channel, 0));
  }
}
void LinearSSN::mnaCompUpdateVoltage(const Matrix &solution) {
  for (Eigen::Index channel = 0; channel < mW.rows(); ++channel) {
    const auto phase = static_cast<UInt>(channel);
    (**mIntfVoltage)(channel, 0) = 0;
    if (terminalNotGrounded(1))
      (**mIntfVoltage)(channel, 0) =
          Math::realFromVectorElement(solution, matrixNodeIndex(1, phase));
    if (terminalNotGrounded(0))
      (**mIntfVoltage)(channel, 0) -=
          Math::realFromVectorElement(solution, matrixNodeIndex(0, phase));
  }
}
} // namespace DPsim::FMI
