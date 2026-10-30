/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim/MNASolverDirect.h>

#include <dpsim-models/MathUtils.h>
#include <dpsim/SequentialScheduler.h>

#include <algorithm>
#include <fstream>
#include <numeric>

using namespace DPsim;
using namespace CPS;

namespace DPsim {

template <typename VarType>
MnaSolverDirect<VarType>::MnaSolverDirect(String name, CPS::Domain domain,
                                          CPS::Logger::Level logLevel)
    : MnaSolver<VarType>(name, domain, logLevel) {
#ifdef WITH_KLU
  mImplementationInUse = DirectLinearSolverImpl::KLU;
#elif defined(WITH_SPARSE)
  mImplementationInUse = DirectLinearSolverImpl::SparseLU;
#else
  mImplementationInUse = DirectLinearSolverImpl::DenseLU;
#endif
}

template <typename VarType>
void MnaSolverDirect<VarType>::switchedMatrixEmpty(std::size_t index) {
  mSwitchedMatrices[std::bitset<SWITCH_NUM>(index)][0].setZero();
}

template <typename VarType>
void MnaSolverDirect<VarType>::switchedMatrixEmpty(std::size_t swIdx,
                                                   Int freqIdx) {
  mSwitchedMatrices[std::bitset<SWITCH_NUM>(swIdx)][freqIdx].setZero();
}

template <typename VarType>
void MnaSolverDirect<VarType>::switchedMatrixStamp(
    std::size_t index, std::vector<std::shared_ptr<CPS::MNAInterface>> &comp) {
  auto bit = std::bitset<SWITCH_NUM>(index);
  auto &sys = mSwitchedMatrices[bit][0];
  for (auto component : comp) {
    component->mnaApplySystemMatrixStamp(sys);
  }
  for (UInt i = 0; i < mSwitches.size(); ++i)
    mSwitches[i]->mnaApplySwitchSystemMatrixStamp(bit[i], sys, 0);

  // Compute LU-factorization for system matrix
  mDirectLinearSolvers[bit][0]->preprocessing(sys,
                                              mListVariableSystemMatrixEntries);
  auto start = std::chrono::steady_clock::now();
  mDirectLinearSolvers[bit][0]->factorize(sys);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<Real> diff = end - start;
  mFactorizeTimes.push_back(diff.count());
}

template <typename VarType>
void MnaSolverDirect<VarType>::stampVariableSystemMatrix() {

  this->mDirectLinearSolverVariableSystemMatrix =
      createDirectSolverImplementation(mSLog);
  // TODO: a direct linear solver configuration is only applied if system matrix recomputation is used
  this->mDirectLinearSolverVariableSystemMatrix->setConfiguration(
      mConfigurationInUse);
  if (auto dense = std::dynamic_pointer_cast<DenseLUAdapter>(
          mDirectLinearSolverVariableSystemMatrix))
    dense->setDetailedProfiling(Solver::mLogSolveTimes);

  SPDLOG_LOGGER_INFO(mSLog,
                     "Number of variable Elements: {}"
                     "\nNumber of MNA components: {}",
                     mVariableComps.size(), mMNAComponents.size());

  // Build base matrix with only static elements
  mBaseSystemMatrix.setZero();
  for (auto statElem : mMNAComponents)
    statElem->mnaApplySystemMatrixStamp(mBaseSystemMatrix);
  SPDLOG_LOGGER_INFO(mSLog, "Base matrix with only static elements: {}",
                     Logger::matrixToString(mBaseSystemMatrix));
  mSLog->flush();

  // Continue from base matrix
  mVariableSystemMatrix = mBaseSystemMatrix;

  // Now stamp initial state of variable elements and switches into matrix
  SPDLOG_LOGGER_INFO(mSLog, "Stamping variable elements");
  for (auto varElem : mMNAIntfVariableComps)
    varElem->mnaApplySystemMatrixStamp(mVariableSystemMatrix);

  SPDLOG_LOGGER_INFO(mSLog, "Initial system matrix with variable elements {}",
                     Logger::matrixToString(mVariableSystemMatrix));
  /* TODO: find replacement for flush() */
  mSLog->flush();

  // Calculate factorization of current matrix
  mDirectLinearSolverVariableSystemMatrix->preprocessing(
      mVariableSystemMatrix, mListVariableSystemMatrixEntries);
  CPS::Math::prepareSparseMatrixValueCache(mVariableSystemMatrix);

  auto start = std::chrono::steady_clock::now();
  mDirectLinearSolverVariableSystemMatrix->factorize(mVariableSystemMatrix);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<Real> diff = end - start;
  mFactorizeTimes.push_back(diff.count());
}

template <typename VarType>
void MnaSolverDirect<VarType>::solveWithSystemMatrixRecomputation(
    Real time, Int timeStepCount) {
  auto profileStart = std::chrono::steady_clock::now();
  // Reset source vector
  mRightSideVector.setZero();
  auto profileEnd = std::chrono::steady_clock::now();
  if (Solver::mLogSolveTimes)
    recordProfileMeasurement("rhs_zero", profileStart, profileEnd);

  // Add together the right side vector (computed by the components'
  // pre-step tasks)
  profileStart = std::chrono::steady_clock::now();
  addRightVectorStamps(mRightSideVector);
  profileEnd = std::chrono::steady_clock::now();
  if (Solver::mLogSolveTimes)
    recordProfileMeasurement("rhs_assembly", profileStart, profileEnd);

  // Get switch and variable comp status and update system matrix and lu factorization accordingly
  profileStart = std::chrono::steady_clock::now();
  mVariableComponentChanged = hasVariableComponentChanged();
  profileEnd = std::chrono::steady_clock::now();
  if (Solver::mLogSolveTimes)
    recordProfileMeasurement("variable_change_detection", profileStart,
                             profileEnd);
  if (mVariableComponentChanged)
    recomputeSystemMatrix(time);

  // Calculate new solution vector
  auto start = std::chrono::steady_clock::now();
  mDirectLinearSolverVariableSystemMatrix->solveInto(mRightSideVector,
                                                     *mLeftSideVectorValue);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<Real> diff = end - start;
  mSolveTimes.push_back(diff.count());
  if (Solver::mLogSolveTimes)
    recordProfileMeasurement("numeric_solve", start, end);

  // TODO split into separate task? (dependent on x, updating all v attributes)
  profileStart = std::chrono::steady_clock::now();
  for (UInt nodeIdx = 0; nodeIdx < mNumNetNodes; ++nodeIdx)
    mNodes[nodeIdx]->mnaUpdateVoltage(*mLeftSideVectorValue);
  profileEnd = std::chrono::steady_clock::now();
  if (Solver::mLogSolveTimes)
    recordProfileMeasurement("node_voltage_propagation", profileStart,
                             profileEnd);

  // Components' states will be updated by the post-step tasks
}

template <typename VarType>
void MnaSolverDirect<VarType>::recomputeSystemMatrix(Real time) {
  const auto profileStart = std::chrono::steady_clock::now();
  // The sparse structure is fixed after initialization. Reset only the
  // coefficient values instead of copying/reallocating the full sparse
  // matrix on every timestep.
  std::copy_n(mBaseSystemMatrix.valuePtr(), mBaseSystemMatrix.nonZeros(),
              mVariableSystemMatrix.valuePtr());

  // Now stamp variable elements and switches into matrix
  for (auto comp : mMNAIntfVariableComps)
    comp->mnaApplySystemMatrixStamp(mVariableSystemMatrix);
  const auto profileEnd = std::chrono::steady_clock::now();
  if (Solver::mLogSolveTimes)
    recordProfileMeasurement("matrix_reset_and_variable_stamping", profileStart,
                             profileEnd);

  // Refactorization of matrix assuming that structure remained
  // constant by omitting analyzePattern
  auto start = std::chrono::steady_clock::now();
  mDirectLinearSolverVariableSystemMatrix->partialRefactorize(
      mVariableSystemMatrix, mListVariableSystemMatrixEntries);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<Real> diff = end - start;
  mRecomputationTimes.push_back(diff.count());
  if (Solver::mLogSolveTimes)
    recordProfileMeasurement("partial_refactorization", start, end);
  ++mNumRecomputations;
}

template <typename VarType>
void MnaSolverDirect<VarType>::extractStateSpace(Real time) {
  if (!mStateSpaceExtractor)
    return;

  if (mSystemMatrixRecomputationEnabled) {
    mStateSpaceExtractor->extract(*mDirectLinearSolverVariableSystemMatrix,
                                  mVariableComponentChanged,
                                  mVariableComponentChanged, time);
  } else {
    mStateSpaceExtractor->extract(
        *mDirectLinearSolvers[mCurrentSwitchStatus][0], false,
        mSystemMatrixChanged, time);
  }
}

template <> void MnaSolverDirect<Real>::createEmptySystemMatrix() {
  if (mSwitches.size() > SWITCH_NUM)
    throw SystemError("Too many Switches.");

  if (mSystemMatrixRecomputationEnabled) {
    mBaseSystemMatrix =
        SparseMatrix(mNumMatrixNodeIndices, mNumMatrixNodeIndices);
    mVariableSystemMatrix =
        SparseMatrix(mNumMatrixNodeIndices, mNumMatrixNodeIndices);
  } else {
    for (std::size_t i = 0; i < (1ULL << mSwitches.size()); i++) {
      auto bit = std::bitset<SWITCH_NUM>(i);
      mSwitchedMatrices[bit].push_back(
          SparseMatrix(mNumMatrixNodeIndices, mNumMatrixNodeIndices));
      mDirectLinearSolvers[bit].push_back(
          createDirectSolverImplementation(mSLog));
    }
  }
}

template <> void MnaSolverDirect<Complex>::createEmptySystemMatrix() {
  if (mSwitches.size() > SWITCH_NUM)
    throw SystemError("Too many Switches.");

  if (mFrequencyParallel) {
    for (UInt i = 0; i < std::pow(2, mSwitches.size()); ++i) {
      for (Int freq = 0; freq < mSystem.mFrequencies.size(); ++freq) {
        auto bit = std::bitset<SWITCH_NUM>(i);
        mSwitchedMatrices[bit].push_back(SparseMatrix(
            2 * (mNumMatrixNodeIndices), 2 * (mNumMatrixNodeIndices)));
        mDirectLinearSolvers[bit].push_back(
            createDirectSolverImplementation(mSLog));
      }
    }
  } else if (mSystemMatrixRecomputationEnabled) {
    mBaseSystemMatrix =
        SparseMatrix(2 * (mNumMatrixNodeIndices), 2 * (mNumMatrixNodeIndices));
    mVariableSystemMatrix =
        SparseMatrix(2 * (mNumMatrixNodeIndices), 2 * (mNumMatrixNodeIndices));
  } else {
    for (std::size_t i = 0; i < (1ULL << mSwitches.size()); i++) {
      auto bit = std::bitset<SWITCH_NUM>(i);
      mSwitchedMatrices[bit].push_back(SparseMatrix(
          2 * (mNumTotalMatrixNodeIndices), 2 * (mNumTotalMatrixNodeIndices)));
      mDirectLinearSolvers[bit].push_back(
          createDirectSolverImplementation(mSLog));
    }
  }
}

template <typename VarType>
std::shared_ptr<CPS::Task> MnaSolverDirect<VarType>::createSolveTask() {
  prepareCompactRightVectorStamps();
  return std::make_shared<MnaSolverDirect<VarType>::SolveTask>(*this);
}

template <typename VarType>
std::shared_ptr<CPS::Task> MnaSolverDirect<VarType>::createSolveTaskRecomp() {
  prepareCompactRightVectorStamps();
  return std::make_shared<MnaSolverDirect<VarType>::SolveTaskRecomp>(*this);
}

template <typename VarType>
void MnaSolverDirect<VarType>::prepareCompactRightVectorStamps() {
  std::unordered_map<const Matrix *, CPS::MNAInterface *> compactStamps;
  auto collect = [&compactStamps](const auto &components) {
    for (const auto &component : components) {
      if (!component->mnaHasCompactRightSideVectorStamp())
        continue;
      auto stamp = component->getRightVector();
      if (stamp->get().size() != 0)
        compactStamps.emplace(&stamp->get(), component.get());
    }
  };
  collect(mMNAComponents);
  collect(mMNAIntfVariableComps);

  // Resolve compact-vs-dense accumulation once. The solve loop then follows
  // this ordered list directly instead of hashing every RHS stamp each step.
  mRightVectorStampOperations.clear();
  mRightVectorStampOperations.reserve(mRightVectorStamps.size());
  for (const Matrix *stamp : mRightVectorStamps) {
    const auto compact = compactStamps.find(stamp);
    if (compact != compactStamps.end())
      mRightVectorStampOperations.push_back({nullptr, compact->second});
    else
      mRightVectorStampOperations.push_back({stamp, nullptr});
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::addRightVectorStamps(Matrix &rightVector) {
  using Clock = std::chrono::steady_clock;
  const bool profile = Solver::mLogSolveTimes;
  const auto start = profile ? Clock::now() : Clock::time_point{};
  std::chrono::nanoseconds compactTime{};
  for (const auto &operation : mRightVectorStampOperations) {
    if (operation.compact) {
      const auto compactStart = profile ? Clock::now() : Clock::time_point{};
      operation.compact->mnaAddCompactRightSideVectorStamp(rightVector);
      if (profile)
        compactTime += std::chrono::duration_cast<std::chrono::nanoseconds>(
            Clock::now() - compactStart);
      ++mCompactRightVectorStampCalls;
    } else {
      rightVector += *operation.dense;
      ++mDenseRightVectorStampCalls;
    }
  }
  if (profile) {
    const auto end = Clock::now();
    recordProfileMeasurement("rhs_compact_stamp_callbacks", start,
                             start + compactTime);
    recordProfileMeasurement("rhs_dense_stamp_accumulation",
                             start + compactTime, end);
  }
}

template <typename VarType>
std::shared_ptr<CPS::Task>
MnaSolverDirect<VarType>::createSolveTaskHarm(UInt freqIdx) {
  return std::make_shared<MnaSolverDirect<VarType>::SolveTaskHarm>(*this,
                                                                   freqIdx);
}

template <typename VarType>
std::shared_ptr<CPS::Task>
MnaSolverDirect<VarType>::createStateSpaceExtractionTask() {
  return std::make_shared<
      typename MnaSolverDirect<VarType>::StateSpaceExtractionTask>(*this);
}

template <typename VarType>
std::shared_ptr<CPS::Task> MnaSolverDirect<VarType>::createLogTask() {
  return std::make_shared<MnaSolverDirect<VarType>::LogTask>(*this);
}

template <typename VarType>
void MnaSolverDirect<VarType>::solve(Real time, Int timeStepCount) {
  const auto previousSwitchStatus = mCurrentSwitchStatus;
  // Reset source vector
  mRightSideVector.setZero();

  // Add together the right side vector (computed by the components' pre-step tasks)
  addRightVectorStamps(mRightSideVector);

  if (!mIsInInitialization)
    MnaSolver<VarType>::updateSwitchStatus();

  if (mSwitchedMatrices.size() > 0) {
    std::chrono::steady_clock::time_point start;
    if (Solver::mLogSolveTimes)
      start = std::chrono::steady_clock::now();

    mDirectLinearSolvers[mCurrentSwitchStatus][0]->solveInto(
        mRightSideVector, *mLeftSideVectorValue);

    if (Solver::mLogSolveTimes) {
      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<Real> diff = end - start;
      mSolveTimes.push_back(diff.count());
    }
  }

  // CHECK: Is this really required? Or can operations actually become part of
  // correctorStep and mnaPostStep?
  for (auto syncGen : mSyncGen)
    syncGen->updateVoltage(*mLeftSideVectorValue);

  // Reset number of iterations
  mIter = 0;

  // Additional solve steps for iterative models
  if (mSyncGen.size() > 0) {
    UInt numCompsRequireIter;
    do {
      // count synchronous generators that require iteration
      numCompsRequireIter = 0;
      for (auto syncGen : mSyncGen)
        if (syncGen->requiresIteration())
          numCompsRequireIter++;

      // recompute solve step if at least one component demands iteration
      if (numCompsRequireIter > 0) {
        mIter++;

        // Reset source vector
        mRightSideVector.setZero();

        if (!mIsInInitialization)
          MnaSolver<VarType>::updateSwitchStatus();

        for (auto syncGen : mSyncGen)
          syncGen->correctorStep();

        // Add together the right side vector (computed by the components' pre-step tasks)
        addRightVectorStamps(mRightSideVector);

        if (mSwitchedMatrices.size() > 0) {
          auto start = std::chrono::steady_clock::now();
          mDirectLinearSolvers[mCurrentSwitchStatus][0]->solveInto(
              mRightSideVector, *mLeftSideVectorValue);
          auto end = std::chrono::steady_clock::now();
          std::chrono::duration<Real> diff = end - start;
          mSolveTimes.push_back(diff.count());
        }

        // CHECK: Is this really required? Or can operations actually become part of
        // correctorStep and mnaPostStep?
        for (auto syncGen : mSyncGen)
          syncGen->updateVoltage(*mLeftSideVectorValue);
      }
    } while (numCompsRequireIter > 0);
  }

  mSystemMatrixChanged = (mCurrentSwitchStatus != previousSwitchStatus);

  // TODO split into separate task? (dependent on x, updating all v attributes)
  for (UInt nodeIdx = 0; nodeIdx < mNumNetNodes; ++nodeIdx)
    mNodes[nodeIdx]->mnaUpdateVoltage(*mLeftSideVectorValue);

  // Components' states will be updated by the post-step tasks
}

template <typename VarType>
void MnaSolverDirect<VarType>::solveWithHarmonics(Real time, Int timeStepCount,
                                                  Int freqIdx) {
  mRightSideVectorHarm[freqIdx].setZero();

  // Sum of right side vectors (computed by the components' pre-step tasks)
  for (auto stamp : mRightVectorStamps)
    mRightSideVectorHarm[freqIdx] += stamp->col(freqIdx);

  mDirectLinearSolvers[mCurrentSwitchStatus][freqIdx]->solveInto(
      mRightSideVectorHarm[freqIdx], **mLeftSideVectorHarm[freqIdx]);
}

template <typename VarType> void MnaSolverDirect<VarType>::logSystemMatrices() {
  if (mFrequencyParallel) {
    for (UInt i = 0; i < mSwitchedMatrices[std::bitset<SWITCH_NUM>(0)].size();
         ++i) {
      SPDLOG_LOGGER_INFO(mSLog, "System matrix for frequency: {:d} \n{:s}", i,
                         Logger::matrixToString(
                             mSwitchedMatrices[std::bitset<SWITCH_NUM>(0)][i]));
    }

    for (UInt i = 0; i < mRightSideVectorHarm.size(); ++i) {
      SPDLOG_LOGGER_INFO(mSLog, "Right side vector for frequency: {:d} \n{:s}",
                         i, Logger::matrixToString(mRightSideVectorHarm[i]));
    }

  } else if (mSystemMatrixRecomputationEnabled) {
    SPDLOG_LOGGER_INFO(mSLog, "Summarizing matrices: ");
    SPDLOG_LOGGER_INFO(mSLog, "Base matrix with only static elements: {}",
                       Logger::matrixToString(mBaseSystemMatrix));
    SPDLOG_LOGGER_INFO(mSLog, "Initial system matrix with variable elements {}",
                       Logger::matrixToString(mVariableSystemMatrix));
    SPDLOG_LOGGER_INFO(mSLog, "Right side vector: {}",
                       Logger::matrixToString(mRightSideVector));
  } else {
    if (mSwitches.size() < 1) {
      SPDLOG_LOGGER_INFO(mSLog, "System matrix: \n{}",
                         mSwitchedMatrices[std::bitset<SWITCH_NUM>(0)][0]);
    } else {
      SPDLOG_LOGGER_INFO(mSLog, "Initial switch status: {:s}",
                         mCurrentSwitchStatus.to_string());

      for (auto sys : mSwitchedMatrices) {
        SPDLOG_LOGGER_INFO(mSLog, "Switching System matrix {:s} \n{:s}",
                           sys.first.to_string(),
                           Logger::matrixToString(sys.second[0]));
      }
    }
    SPDLOG_LOGGER_INFO(mSLog, "Right side vector: \n{}", mRightSideVector);
  }
}

template <typename VarType> void MnaSolverDirect<VarType>::logLUTimes() {
  logFactorizationTime();
  logRecomputationTime();
  logSolveTime();
  for (const auto &entry : mProfileMeasurements) {
    const auto &measurement = entry.second;
    const double averageNs =
        measurement.count == 0
            ? 0.
            : static_cast<double>(measurement.total.count()) /
                  measurement.count;
    SPDLOG_LOGGER_INFO(mSLog,
                       "DEEP_TIMING_PHASE,{},{:.9f},{:.3f},{:.3f},{:.3f},{}",
                       entry.first, measurement.total.count() * 1e-9,
                       averageNs * 1e-3, measurement.minimum.count() * 1e-3,
                       measurement.maximum.count() * 1e-3, measurement.count);
  }
  std::ofstream phaseFile(CPS::Logger::logDir() + "/" + this->mName +
                          "_SolverPhaseTimes.csv");
  phaseFile << "phase,total_s,avg_us,min_us,max_us,calls\n";
  for (const auto &entry : mProfileMeasurements) {
    const auto &measurement = entry.second;
    const double count = static_cast<double>(measurement.count);
    phaseFile << entry.first << "," << measurement.total.count() * 1e-9 << ","
              << (count ? measurement.total.count() / count * 1e-3 : 0.) << ","
              << measurement.minimum.count() * 1e-3 << ","
              << measurement.maximum.count() * 1e-3 << "," << measurement.count
              << "\n";
  }
  phaseFile << "rhs_dense_stamp_component_calls,0,0,0,0,"
            << mDenseRightVectorStampCalls << "\n";
  phaseFile << "rhs_compact_stamp_component_calls,0,0,0,0,"
            << mCompactRightVectorStampCalls << "\n";
  for (const auto &operation : mRightVectorStampOperations)
    if (operation.compact)
      operation.compact->mnaWriteProfile(phaseFile);
  if (auto dense = std::dynamic_pointer_cast<DenseLUAdapter>(
          mDirectLinearSolverVariableSystemMatrix))
    dense->writeDetailedProfile(phaseFile);
  if (auto dense = std::dynamic_pointer_cast<DenseLUAdapter>(
          mDirectLinearSolverVariableSystemMatrix)) {
    std::ofstream diagnostics(CPS::Logger::logDir() + "/" + this->mName +
                              "_SolverDiagnostics.csv");
    dense->writeDetailedDiagnostics(diagnostics);
  }
  std::ofstream matrixFile(CPS::Logger::logDir() + "/" + this->mName +
                           "_SolverMatrixStats.csv");
  matrixFile << "rows,cols,nonzeros,structural_density,variable_entries\n";
  const double matrixSize = static_cast<double>(mVariableSystemMatrix.rows());
  const double matrixEntries = matrixSize * matrixSize;
  matrixFile << mVariableSystemMatrix.rows() << ","
             << mVariableSystemMatrix.cols() << ","
             << mVariableSystemMatrix.nonZeros() << ","
             << (matrixEntries
                     ? mVariableSystemMatrix.nonZeros() / matrixEntries
                     : 0.)
             << "," << mListVariableSystemMatrixEntries.size() << "\n";
  if (!mFactorizeTimes.empty()) {
    const Real total = std::accumulate(mFactorizeTimes.begin(),
                                       mFactorizeTimes.end(), Real(0.));
    const auto bounds =
        std::minmax_element(mFactorizeTimes.begin(), mFactorizeTimes.end());
    phaseFile << "initial_factorization," << total << ","
              << total * 1e6 / mFactorizeTimes.size() << ","
              << *bounds.first * 1e6 << "," << *bounds.second * 1e6 << ","
              << mFactorizeTimes.size() << "\n";
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::recordProfileMeasurement(
    const String &phase, std::chrono::steady_clock::time_point start,
    std::chrono::steady_clock::time_point end) {
  const auto elapsed =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  auto &measurement = mProfileMeasurements[phase];
  if (measurement.count == 0) {
    measurement.minimum = elapsed;
    measurement.maximum = elapsed;
  } else {
    measurement.minimum = std::min(measurement.minimum, elapsed);
    measurement.maximum = std::max(measurement.maximum, elapsed);
  }
  measurement.total += elapsed;
  ++measurement.count;
}

template <typename VarType> void MnaSolverDirect<VarType>::logSolveTime() {
  Real solveSum = 0.0;
  Real solveMax = 0.0;
  for (auto meas : mSolveTimes) {
    solveSum += meas;
    if (meas > solveMax)
      solveMax = meas;
  }
  SPDLOG_LOGGER_INFO(mSLog, "Cumulative solve times: {:.12f}", solveSum);
  SPDLOG_LOGGER_INFO(mSLog, "Average solve time: {:.12f}",
                     solveSum / static_cast<double>(mSolveTimes.size()));
  SPDLOG_LOGGER_INFO(mSLog, "Maximum solve time: {:.12f}", solveMax);
  SPDLOG_LOGGER_INFO(mSLog, "Number of solves: {:d}", mSolveTimes.size());
}

template <typename VarType>
void MnaSolverDirect<VarType>::logFactorizationTime() {
  for (auto meas : mFactorizeTimes) {
    SPDLOG_LOGGER_INFO(mSLog, "LU factorization time: {:.12f}", meas);
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::logRecomputationTime() {
  Real recompSum = 0.0;
  Real recompMax = 0.0;
  for (auto meas : mRecomputationTimes) {
    recompSum += meas;
    if (meas > recompMax)
      recompMax = meas;
  }
  // Sometimes, refactorization is not used
  if (mRecomputationTimes.size() != 0) {
    SPDLOG_LOGGER_INFO(mSLog, "Cumulative refactorization times: {:.12f}",
                       recompSum);
    SPDLOG_LOGGER_INFO(mSLog, "Average refactorization time: {:.12f}",
                       recompSum / ((double)mRecomputationTimes.size()));
    SPDLOG_LOGGER_INFO(mSLog, "Maximum refactorization time: {:.12f}",
                       recompMax);
    SPDLOG_LOGGER_INFO(mSLog, "Number of refactorizations: {:d}",
                       mRecomputationTimes.size());
  }
}

template <typename VarType>
std::shared_ptr<DirectLinearSolver>
MnaSolverDirect<VarType>::createDirectSolverImplementation(
    CPS::Logger::Log mSLog) {
  switch (this->mImplementationInUse) {
  case DirectLinearSolverImpl::DenseLU:
    return std::make_shared<DenseLUAdapter>(mSLog);
#ifdef WITH_SPARSE
  case DirectLinearSolverImpl::SparseLU:
    return std::make_shared<SparseLUAdapter>(mSLog);
#endif // WITH_SPARSE
#ifdef WITH_KLU
  case DirectLinearSolverImpl::KLU:
    return std::make_shared<KLUAdapter>(mSLog);
#endif
#ifdef WITH_CUDA
  case DirectLinearSolverImpl::CUDADense:
    return std::make_shared<GpuDenseAdapter>(mSLog);
#ifdef WITH_CUDA_SPARSE
  case DirectLinearSolverImpl::CUDASparse:
    return std::make_shared<GpuSparseAdapter>(mSLog);
#endif
#ifdef WITH_MAGMA
  case DirectLinearSolverImpl::CUDAMagma:
    return std::make_shared<GpuMagmaAdapter>(mSLog);
#endif
#endif
  default:
    throw CPS::SystemError("unsupported linear solver implementation.");
  }
}

template <typename VarType>
void MnaSolverDirect<VarType>::setDirectLinearSolverImplementation(
    DirectLinearSolverImpl implementation) {
  this->mImplementationInUse = implementation;
}

template <typename VarType>
void MnaSolverDirect<VarType>::setDirectLinearSolverConfiguration(
    DirectLinearSolverConfiguration &configuration) {
  this->mConfigurationInUse = configuration;
}

} // namespace DPsim

template class DPsim::MnaSolverDirect<Real>;
template class DPsim::MnaSolverDirect<Complex>;
