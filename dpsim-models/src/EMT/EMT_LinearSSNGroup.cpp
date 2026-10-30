// SPDX-License-Identifier: MPL-2.0
#include <dpsim-models/EMT/EMT_LinearSSNGroup.h>

#include <dpsim-models/EMT/EMT_Ph3_SSN_ImpedanceLoad.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_PiLine.h>
#include <dpsim-models/EMT/EMT_Ph3_SSN_Transformer.h>
#include <dpsim-models/MNAStampUtils.h>

#include <stdexcept>

using namespace CPS;

namespace {
using PhaseVector = Eigen::Matrix<Real, 3, 1>;

Matrix modelHistory(const EMT::LinearSSNModel &model,
                    Eigen::Ref<const Matrix> state) {
  if (model.kind == EMT::LinearSSNModel::Kind::PiLine) {
    const PhaseVector current = state.block<3, 1>(0, 0);
    const PhaseVector inductorVoltage = state.block<3, 1>(3, 0);
    const PhaseVector capVoltage0 = state.block<3, 1>(6, 0);
    const PhaseVector capCurrent0 = state.block<3, 1>(9, 0);
    const PhaseVector capVoltage1 = state.block<3, 1>(12, 0);
    const PhaseVector capCurrent1 = state.block<3, 1>(15, 0);
    const PhaseVector series =
        model.seriesSystemInverse *
        (model.inductorAdmittance * inductorVoltage + current);
    Matrix result(6, 1);
    result.block<3, 1>(0, 0) =
        series - (-capCurrent0 - model.capAdmittance * capVoltage0);
    result.block<3, 1>(3, 0) =
        -series - (-capCurrent1 - model.capAdmittance * capVoltage1);
    return result;
  }
  if (model.kind == EMT::LinearSSNModel::Kind::ImpedanceLoad) {
    Matrix result = Matrix::Zero(3, 1);
    for (UInt phase = 0; phase < 3; ++phase) {
      Real loadHistory =
          model.hasLoadInductance
              ? model.inductorAdmittanceScalar * state.coeff(phase, 0) +
                    state.coeff(3 + phase, 0)
              : 0.;
      Real shuntHistory =
          model.hasLoadCapacitance
              ? -state.coeff(9 + phase, 0) -
                    model.capacitorAdmittanceScalar * state.coeff(6 + phase, 0)
              : 0.;
      result(phase) = -(loadHistory + shuntHistory);
    }
    return result;
  }

  const PhaseVector h0 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(0, 0) +
                          state.block<3, 1>(3, 0));
  const PhaseVector h1 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(6, 0) +
                          state.block<3, 1>(9, 0));
  const PhaseVector hMag =
      model.magnetizingInductorAdmittance * state.block<3, 1>(12, 0) +
      state.block<3, 1>(15, 0);
  const PhaseVector midpointSource = model.midpointInverse * (-h0 + h1 + hMag);
  Matrix result(6, 1);
  result.block<3, 1>(0, 0) = h0 + model.seriesAdmittance * midpointSource;
  result.block<3, 1>(3, 0) =
      model.ratioFromReference * (model.seriesAdmittance * midpointSource - h1);
  if (model.referenceTerminal == 1)
    result.block<3, 1>(0, 0).swap(result.block<3, 1>(3, 0));
  return result;
}

Matrix modelNext(const EMT::LinearSSNModel &model,
                 Eigen::Ref<const Matrix> state,
                 Eigen::Ref<const Matrix> input) {
  if (model.kind == EMT::LinearSSNModel::Kind::PiLine) {
    const PhaseVector v0 = input.block<3, 1>(0, 0);
    const PhaseVector v1 = input.block<3, 1>(3, 0);
    const PhaseVector interfaceVoltage = v1 - v0;
    const PhaseVector series =
        model.seriesSystemInverse *
        (model.inductorAdmittance * state.block<3, 1>(3, 0) +
         state.block<3, 1>(0, 0));
    const PhaseVector current =
        model.seriesAdmittance * interfaceVoltage + series;
    const PhaseVector cap0 = model.capAdmittance * v0 +
                             (-state.block<3, 1>(9, 0) -
                              model.capAdmittance * state.block<3, 1>(6, 0));
    const PhaseVector cap1 = model.capAdmittance * v1 +
                             (-state.block<3, 1>(15, 0) -
                              model.capAdmittance * state.block<3, 1>(12, 0));
    Matrix result(18, 1);
    result.block<3, 1>(0, 0) = current;
    result.block<3, 1>(3, 0) =
        interfaceVoltage - model.seriesResistance * current;
    result.block<3, 1>(6, 0) = v0;
    result.block<3, 1>(9, 0) = cap0;
    result.block<3, 1>(12, 0) = v1;
    result.block<3, 1>(15, 0) = cap1;
    return result;
  }
  if (model.kind == EMT::LinearSSNModel::Kind::ImpedanceLoad) {
    Matrix result = state;
    for (UInt phase = 0; phase < 3; ++phase) {
      const Real voltage = input.coeff(phase, 0);
      if (model.hasLoadInductance) {
        result(phase) = voltage;
        result(3 + phase) =
            model.inductorAdmittanceScalar * voltage +
            (model.inductorAdmittanceScalar * state.coeff(phase, 0) +
             state.coeff(3 + phase, 0));
      }
      if (model.hasLoadCapacitance) {
        result(6 + phase) = voltage;
        result(9 + phase) =
            model.capacitorAdmittanceScalar * voltage -
            state.coeff(9 + phase, 0) -
            model.capacitorAdmittanceScalar * state.coeff(6 + phase, 0);
      }
    }
    return result;
  }

  PhaseVector v0 = input.block<3, 1>(0, 0);
  PhaseVector v1 = input.block<3, 1>(3, 0);
  const UInt ref = model.referenceTerminal;
  const PhaseVector referenceVoltage = ref == 0 ? v0 : v1;
  const PhaseVector otherVoltage = ref == 0 ? v1 : v0;
  const PhaseVector referredOther = model.ratioFromReference * otherVoltage;
  const PhaseVector h0 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(0, 0) +
                          state.block<3, 1>(3, 0));
  const PhaseVector h1 = model.seriesSystemInverse *
                         (model.inductorAdmittance * state.block<3, 1>(6, 0) +
                          state.block<3, 1>(9, 0));
  const PhaseVector hMag =
      model.magnetizingInductorAdmittance * state.block<3, 1>(12, 0) +
      state.block<3, 1>(15, 0);
  const PhaseVector midpoint =
      model.midpointInverse *
      (model.seriesAdmittance * referenceVoltage +
       model.seriesAdmittance * referredOther - h0 + h1 + hMag);
  const PhaseVector total0 = midpoint - referenceVoltage;
  const PhaseVector current0 = model.seriesAdmittance * total0 + h0;
  const PhaseVector total1 = referredOther - midpoint;
  const PhaseVector current1 = model.seriesAdmittance * total1 + h1;
  Matrix result(18, 1);
  result.block<3, 1>(0, 0) = total0 - model.seriesResistanceHalf * current0;
  result.block<3, 1>(3, 0) = current0;
  result.block<3, 1>(6, 0) = total1 - model.seriesResistanceHalf * current1;
  result.block<3, 1>(9, 0) = current1;
  result.block<3, 1>(12, 0) = -midpoint;
  result.block<3, 1>(15, 0) =
      model.magnetizingInductorAdmittance * (-midpoint) + hMag;
  return result;
}

void addBlock(Matrix &target, const std::vector<UInt> &rows,
              const std::vector<UInt> &cols, const Matrix &block) {
  for (UInt i = 0; i < rows.size(); ++i)
    for (UInt j = 0; j < cols.size(); ++j)
      target(rows[i], cols[j]) += block(i, j);
}
} // namespace

EMT::LinearSSNGroup::LinearSSNGroup(String name, ComponentList components,
                                    SimNode::List allNodes,
                                    SimNode::List retainedNodes,
                                    Logger::Level level)
    : MNASimPowerComp<Real>(name, true, true, level),
      mComponents(std::move(components)), mAllNodes(std::move(allNodes)) {
  mPhaseType = PhaseType::ABC;
  const UInt retainedCount = static_cast<UInt>(retainedNodes.size());
  setTerminalNumber(retainedCount);
  connect(std::move(retainedNodes));
  for (UInt i = 0; i < retainedCount; ++i)
    mRetainedNames.push_back(node(i)->name());
  for (UInt i = 0; i < mRetainedNames.size(); ++i)
    mRetainedNodeIndex[mRetainedNames[i]] = i;
  for (UInt i = 0; i < mAllNodes.size(); ++i)
    mLocalNodeIndex[mAllNodes[i]->name()] = i;
  mState = mAttributes->create<Matrix>("group_state");
  mState->get().resize(0, 1);
}

void EMT::LinearSSNGroup::initializeFromNodesAndTerminals(Real frequency) {
  for (const auto &component : mComponents)
    component->initializeFromNodesAndTerminals(frequency);
}

void EMT::LinearSSNGroup::buildRealization(Real omega, Real timeStep,
                                           Attribute<Matrix>::Ptr leftVector) {
  mModels.clear();
  for (const auto &component : mComponents) {
    component->mnaInitialize(omega, timeStep, leftVector);
    LinearSSNModel model;
    if (auto line = std::dynamic_pointer_cast<Ph3::SSN::PiLine>(component))
      line->exportFusedModel(model);
    else if (auto load =
                 std::dynamic_pointer_cast<Ph3::SSN::ImpedanceLoad>(component))
      load->exportFusedModel(model);
    else if (auto transformer =
                 std::dynamic_pointer_cast<Ph3::SSN::Transformer>(component))
      transformer->exportFusedModel(model);
    else
      throw std::invalid_argument("Unsupported component in LinearSSNGroup.");
    mModels.push_back(std::move(model));
  }

  const UInt localNodeCount = static_cast<UInt>(mAllNodes.size()) * 3;
  const UInt stateCount = [&] {
    UInt n = 0;
    for (const auto &model : mModels)
      n += model.state.rows();
    return n;
  }();
  Matrix G = Matrix::Zero(localNodeCount, localNodeCount);
  Matrix H = Matrix::Zero(localNodeCount, stateCount);
  Matrix F = Matrix::Zero(stateCount, stateCount);
  Matrix B = Matrix::Zero(stateCount, localNodeCount);
  mStateOffsets.clear();
  mStateOffsets.reserve(mModels.size());
  mModelInputIndices.clear();
  mModelInputIndices.reserve(mModels.size());
  Matrix initialState = Matrix::Zero(stateCount, 1);
  UInt stateOffset = 0;
  for (const auto &model : mModels) {
    mStateOffsets.push_back(stateOffset);
    std::vector<UInt> inputIndices(model.terminalCount * 3,
                                   std::numeric_limits<UInt>::max());
    for (UInt p = 0; p < inputIndices.size(); ++p) {
      const UInt terminal = p / 3;
      if (!model.nodeNames[p].empty() && !model.grounded[terminal])
        inputIndices[p] = mLocalNodeIndex.at(model.nodeNames[p]) * 3 + p % 3;
    }
    mModelInputIndices.push_back(std::move(inputIndices));
    initialState.block(stateOffset, 0, model.state.rows(), 1) = model.state;
    const UInt inputCount = model.terminalCount * 3;
    Matrix zeroState = Matrix::Zero(model.state.rows(), 1);
    Matrix zeroInput = Matrix::Zero(inputCount, 1);
    Matrix localF = Matrix::Zero(model.state.rows(), model.state.rows());
    Matrix localB = Matrix::Zero(model.state.rows(), inputCount);
    Matrix localH = Matrix::Zero(inputCount, model.state.rows());
    for (UInt i = 0; i < model.state.rows(); ++i) {
      Matrix basis = Matrix::Zero(model.state.rows(), 1);
      basis(i) = 1.;
      localF.col(i) = modelNext(model, basis, zeroInput);
      localH.col(i) = modelHistory(model, basis);
    }
    for (UInt i = 0; i < inputCount; ++i) {
      Matrix basis = Matrix::Zero(inputCount, 1);
      basis(i) = 1.;
      localB.col(i) = modelNext(model, zeroState, basis);
    }
    F.block(stateOffset, stateOffset, model.state.rows(), model.state.rows()) =
        localF;
    for (UInt p = 0; p < inputCount; ++p) {
      const String &name = model.nodeNames[p];
      const UInt terminal = p / 3;
      if (name.empty() || model.grounded[terminal])
        continue;
      const UInt base = mLocalNodeIndex.at(name) * 3 + (p % 3);
      B.block(stateOffset, base, model.state.rows(), 1) += localB.col(p);
    }
    for (UInt p = 0; p < inputCount; ++p) {
      const String &rowName = model.nodeNames[p];
      const UInt terminal = p / 3;
      if (rowName.empty() || model.grounded[terminal])
        continue;
      const UInt row = mLocalNodeIndex.at(rowName) * 3 + (p % 3);
      for (UInt q = 0; q < model.state.rows(); ++q)
        H(row, stateOffset + q) += localH(p, q);
    }
    for (UInt p = 0; p < inputCount; ++p) {
      const UInt terminal = p / 3;
      if (model.nodeNames[p].empty() || model.grounded[terminal])
        continue;
      const UInt row = mLocalNodeIndex.at(model.nodeNames[p]) * 3 + (p % 3);
      for (UInt q = 0; q < inputCount; ++q) {
        const UInt qTerminal = q / 3;
        if (model.nodeNames[q].empty() || model.grounded[qTerminal])
          continue;
        const UInt col = mLocalNodeIndex.at(model.nodeNames[q]) * 3 + (q % 3);
        G(row, col) += model.portAdmittance(p, q);
      }
    }
    stateOffset += model.state.rows();
  }
  mState->get() = initialState;

  mTerminalIndices.clear();
  for (const auto &name : mRetainedNames) {
    const UInt localBase = mLocalNodeIndex.at(name) * 3;
    for (UInt phase = 0; phase < 3; ++phase)
      mTerminalIndices.push_back(localBase + phase);
  }
  std::vector<Bool> isTerminal(localNodeCount, false);
  for (UInt index : mTerminalIndices)
    isTerminal[index] = true;
  mInternalIndices.clear();
  for (UInt i = 0; i < localNodeCount; ++i)
    if (!isTerminal[i])
      mInternalIndices.push_back(i);

  Matrix Gii = Matrix::Zero(mInternalIndices.size(), mInternalIndices.size());
  Matrix Gin = Matrix::Zero(mInternalIndices.size(), mTerminalIndices.size());
  Matrix Gni = Matrix::Zero(mTerminalIndices.size(), mInternalIndices.size());
  Matrix Gnn = Matrix::Zero(mTerminalIndices.size(), mTerminalIndices.size());
  Matrix Hi = Matrix::Zero(mInternalIndices.size(), stateCount);
  Matrix Hn = Matrix::Zero(mTerminalIndices.size(), stateCount);
  Matrix Bi = Matrix::Zero(stateCount, mInternalIndices.size());
  Matrix Bn = Matrix::Zero(stateCount, mTerminalIndices.size());
  for (UInt i = 0; i < mInternalIndices.size(); ++i) {
    Hi.row(i) = H.row(mInternalIndices[i]);
    for (UInt j = 0; j < mInternalIndices.size(); ++j)
      Gii(i, j) = G(mInternalIndices[i], mInternalIndices[j]);
    for (UInt j = 0; j < mTerminalIndices.size(); ++j)
      Gin(i, j) = G(mInternalIndices[i], mTerminalIndices[j]);
    Bi.col(i) = B.col(mInternalIndices[i]);
  }
  for (UInt i = 0; i < mTerminalIndices.size(); ++i) {
    Hn.row(i) = H.row(mTerminalIndices[i]);
    for (UInt j = 0; j < mInternalIndices.size(); ++j)
      Gni(i, j) = G(mTerminalIndices[i], mInternalIndices[j]);
    for (UInt j = 0; j < mTerminalIndices.size(); ++j)
      Gnn(i, j) = G(mTerminalIndices[i], mTerminalIndices[j]);
    Bn.col(i) = B.col(mTerminalIndices[i]);
  }
  Eigen::PartialPivLU<Matrix> factor(Gii);
  mInternalFactor = std::move(factor);
  mF = F.sparseView();
  mBi = Bi.sparseView();
  mBn = Bn.sparseView();
  mHi = Hi.sparseView();
  mHn = Hn.sparseView();
  mGin = Gin.sparseView();
  mGni = Gni.sparseView();
  mDd = mInternalIndices.empty() ? Gnn : Gnn - Gni * mInternalFactor.solve(Gin);
  mFullHistoryScratch = Matrix::Zero(localNodeCount, 1);
  mInternalHistoryScratch = Matrix::Zero(mInternalIndices.size(), 1);
  mInternalVoltageScratch = Matrix::Zero(mInternalIndices.size(), 1);
  mLocalVoltageScratch = Matrix::Zero(localNodeCount, 1);
  mTerminalInputScratch = Matrix::Zero(mTerminalIndices.size(), 1);
  mTerminalHistoryScratch = Matrix::Zero(mTerminalIndices.size(), 1);
  mStateUpdateScratch = Matrix::Zero(stateCount, 1);
  mReady = true;
}

void EMT::LinearSSNGroup::mnaCompInitialize(Real omega, Real timeStep,
                                            Attribute<Matrix>::Ptr leftVector) {
  updateMatrixNodeIndices();
  buildRealization(omega, timeStep, leftVector);
}

void EMT::LinearSSNGroup::assembleFullHistory(Matrix &fullHistory) const {
  fullHistory.setZero();
  const Matrix &state = mState->get();
  for (UInt modelIndex = 0; modelIndex < mModels.size(); ++modelIndex) {
    const auto &model = mModels[modelIndex];
    const Matrix localHistory =
        modelHistory(model, state.block(mStateOffsets[modelIndex], 0,
                                        model.state.rows(), 1));
    const auto &indices = mModelInputIndices[modelIndex];
    for (UInt p = 0; p < indices.size(); ++p)
      if (indices[p] != std::numeric_limits<UInt>::max())
        fullHistory(indices[p]) += localHistory(p);
  }
}

void EMT::LinearSSNGroup::currentHistory(Matrix &result) const {
  if (mInternalIndices.empty()) {
    result = mHn * mState->get();
    return;
  }

  // H was assembled from the component history maps during realization.
  // Evaluate the precomputed maps directly instead of rebuilding every
  // component history and scattering through the full network vector.
  mInternalHistoryScratch = mHi * mState->get();
  mInternalVoltageScratch = mInternalFactor.solve(mInternalHistoryScratch);
  result = mHn * mState->get();
  result -= mGni * mInternalVoltageScratch;
}

Matrix EMT::LinearSSNGroup::terminalVoltage(const Matrix &leftVector) {
  updateMatrixNodeIndices();
  Matrix voltage = Matrix::Zero(mTerminalIndices.size(), 1);
  for (UInt terminal = 0; terminal < mRetainedNames.size(); ++terminal)
    for (UInt phase = 0; phase < 3; ++phase)
      voltage(terminal * 3 + phase) = Math::realFromVectorElement(
          leftVector, matrixNodeIndex(terminal, phase));
  return voltage;
}

void EMT::LinearSSNGroup::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {
  updateMatrixNodeIndices();
  if (!mReady)
    throw std::logic_error("LinearSSNGroup used before initialization.");
  for (UInt i = 0; i < mRetainedNames.size() * 3; ++i)
    for (UInt j = 0; j < mRetainedNames.size() * 3; ++j)
      Math::addToMatrixElement(systemMatrix, matrixNodeIndex(i / 3, i % 3),
                               matrixNodeIndex(j / 3, j % 3), mDd(i, j));
}

void EMT::LinearSSNGroup::mnaCompApplyRightSideVectorStamp(
    Matrix &rightVector) {
  (void)rightVector;
}

void EMT::LinearSSNGroup::mnaAddCompactRightSideVectorStamp(
    Matrix &rightVector) {
  updateMatrixNodeIndices();
  currentHistory(mTerminalHistoryScratch);
  for (UInt i = 0; i < mTerminalHistoryScratch.rows(); ++i)
    Math::addToVectorElement(rightVector, matrixNodeIndex(i / 3, i % 3),
                             mTerminalHistoryScratch(i));
}

void EMT::LinearSSNGroup::mnaCompPreStep(Real, Int) {}

void EMT::LinearSSNGroup::mnaCompPostStep(Real, Int,
                                          Attribute<Matrix>::Ptr &leftVector) {
  const Matrix input = terminalVoltage(**leftVector);
  if (mInternalIndices.empty()) {
    mState->get() = mF * mState->get() + mBn * input;
    return;
  }
  const Matrix internal =
      mInternalFactor.solve(mHi * mState->get() - mGin * input);
  mState->get() = mF * mState->get() + mBi * internal + mBn * input;
}

void EMT::LinearSSNGroup::mnaCompAddPreStepDependencies(
    AttributeBase::List &prevStepDependencies, AttributeBase::List &,
    AttributeBase::List &modifiedAttributes) {
  prevStepDependencies.push_back(mState);
  modifiedAttributes.push_back(mRightVector);
}

void EMT::LinearSSNGroup::mnaCompAddPostStepDependencies(
    AttributeBase::List &, AttributeBase::List &attributeDependencies,
    AttributeBase::List &modifiedAttributes,
    Attribute<Matrix>::Ptr &leftVector) {
  attributeDependencies.push_back(leftVector);
  modifiedAttributes.push_back(mState);
}
