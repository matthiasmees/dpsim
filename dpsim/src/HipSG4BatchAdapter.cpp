/* SPDX-License-Identifier: MPL-2.0 */
#include <dpsim/HipSG4BatchAdapter.h>
#include <dpsim/HipSG4BatchKernels.h>

#include <map>
#include <stdexcept>

namespace DPsim {
struct HipSG4BatchAdapter::DeviceState {
  void *batch = nullptr;
};

HipSG4BatchAdapter::HipSG4BatchAdapter(
    std::vector<std::shared_ptr<CPS::EMT::Ph3::SynchronGenerator4OrderVBR>>
        generators)
    : mGenerators(std::move(generators)), mDevice(new DeviceState) {
  if (mGenerators.empty()) {
    delete mDevice;
    mDevice = nullptr;
    throw std::invalid_argument("HIP SG4 batch cannot be empty");
  }
  std::vector<CPS::EMT::Ph3::HipSG4State> states;
  states.reserve(mGenerators.size());
  for (const auto &generator : mGenerators) {
    if (!generator || !generator->supportsHipBatch()) {
      delete mDevice;
      mDevice = nullptr;
      throw std::invalid_argument(
          "HIP SG4 batch contains an unsupported generator");
    }
    states.push_back(generator->exportHipState());
  }
  try {
    mDevice->batch = hipSG4BatchCreate(states.data(), states.size());
  } catch (...) {
    delete mDevice;
    mDevice = nullptr;
    throw;
  }
  for (const auto &generator : mGenerators)
    generator->mDeviceManaged = true;
}

HipSG4BatchAdapter::~HipSG4BatchAdapter() {
  for (const auto &generator : mGenerators)
    if (generator) {
      generator->mDeviceManaged = false;
      generator->mDeviceMatrixManaged = false;
      generator->mDeviceRhsManaged = false;
    }
  if (mDevice) {
    hipSG4BatchDestroy(mDevice->batch);
    delete mDevice;
  }
}

void HipSG4BatchAdapter::preStep(CPS::Real time) {
  std::vector<CPS::EMT::Ph3::HipSG4State> states(mGenerators.size());
  hipSG4BatchPreStep(mDevice->batch, time,
                     mDeviceAssemblyConfigured ? nullptr : states.data());
  if (!mDeviceAssemblyConfigured)
    for (std::size_t index = 0; index < mGenerators.size(); ++index)
      mGenerators[index]->applyHipPreState(states[index], time);
}

void HipSG4BatchAdapter::postStep(const CPS::Matrix &leftVector) {
  std::vector<CPS::EMT::Ph3::HipSG4State> states(mGenerators.size());
  hipSG4BatchPostStep(mDevice->batch, leftVector.data(),
                      static_cast<std::size_t>(leftVector.size()),
                      states.data());
  for (std::size_t index = 0; index < mGenerators.size(); ++index)
    mGenerators[index]->applyHipPostState(states[index]);
}

void HipSG4BatchAdapter::configureDeviceAssembly(
    const SparseMatrix &systemMatrix, const SparseMatrix &baseSystemMatrix,
    double *deviceMatrixValues, double *deviceRightHandSide) {
  if (mDeviceAssemblyConfigured)
    throw std::logic_error("HIP SG4 device assembly already configured");
  SparseMatrix compressed = systemMatrix;
  compressed.makeCompressed();
  SparseMatrix compressedBase = baseSystemMatrix;
  compressedBase.makeCompressed();
  if (compressed.rows() != compressedBase.rows() ||
      compressed.cols() != compressedBase.cols() ||
      compressed.nonZeros() != compressedBase.nonZeros())
    throw std::runtime_error("HIP SG4 baseline matrix structure mismatch");
  std::vector<int> positions;
  std::map<int, double> restoreEntries;
  positions.reserve(mGenerators.size() * 9);
  for (const auto &generator : mGenerators) {
    const auto state = generator->exportHipState();
    for (int rowPhase = 0; rowPhase < 3; ++rowPhase) {
      const int row = state.busRows[rowPhase];
      if (state.grounded || row < 0)
        throw std::runtime_error("HIP SG4 device assembly requires an "
                                 "ungrounded generator");
      for (int columnPhase = 0; columnPhase < 3; ++columnPhase) {
        const int column = state.busRows[columnPhase];
        int position = -1;
        for (int index = compressed.outerIndexPtr()[row];
             index < compressed.outerIndexPtr()[row + 1]; ++index)
          if (compressed.innerIndexPtr()[index] == column) {
            position = index;
            break;
          }
        if (position < 0)
          throw std::runtime_error("HIP SG4 conductance entry missing from "
                                   "sparse matrix pattern");
        positions.push_back(position);
        restoreEntries[position] = compressedBase.valuePtr()[position];
      }
    }
  }
  std::vector<int> restorePositions;
  std::vector<double> restoreValues;
  restorePositions.reserve(restoreEntries.size());
  restoreValues.reserve(restoreEntries.size());
  for (const auto &[position, value] : restoreEntries) {
    restorePositions.push_back(position);
    restoreValues.push_back(value);
  }
  hipSG4BatchConfigureAssembly(
      mDevice->batch, positions.data(), deviceMatrixValues, deviceRightHandSide,
      restorePositions.data(), restoreValues.data(), restorePositions.size());
  for (const auto &generator : mGenerators) {
    generator->mDeviceMatrixManaged = true;
    generator->mDeviceRhsManaged = true;
  }
  mDeviceAssemblyConfigured = true;
}

void HipSG4BatchAdapter::accumulateDeviceRhs() {
  if (!mDeviceAssemblyConfigured)
    throw std::logic_error("HIP SG4 device assembly is not configured");
  hipSG4BatchAccumulateRhs(mDevice->batch);
}

void HipSG4BatchAdapter::scatterDeviceMatrix() {
  if (!mDeviceAssemblyConfigured)
    throw std::logic_error("HIP SG4 device assembly is not configured");
  hipSG4BatchScatterMatrix(mDevice->batch);
}

void HipSG4BatchAdapter::postStepDevice(const CPS::Matrix &leftVector,
                                        const double *deviceSolution,
                                        bool synchronizeStateToHost) {
  if (!mDeviceAssemblyConfigured || !deviceSolution)
    throw std::logic_error("HIP SG4 device solution is unavailable");
  std::vector<CPS::EMT::Ph3::HipSG4State> states;
  if (synchronizeStateToHost)
    states.resize(mGenerators.size());
  hipSG4BatchPostStepDevice(mDevice->batch, deviceSolution,
                            static_cast<std::size_t>(leftVector.size()),
                            synchronizeStateToHost ? states.data() : nullptr);
  if (synchronizeStateToHost)
    for (std::size_t index = 0; index < mGenerators.size(); ++index)
      mGenerators[index]->applyHipPostState(states[index]);
}
} // namespace DPsim
