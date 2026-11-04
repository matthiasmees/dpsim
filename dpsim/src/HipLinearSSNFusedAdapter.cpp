/* SPDX-License-Identifier: MPL-2.0 */
#include <dpsim/HipLinearSSNFusedAdapter.h>
#include <dpsim/HipLinearSSNFusedKernels.h>

#include <stdexcept>

namespace DPsim {
struct HipLinearSSNFusedAdapter::DeviceState {
  void *models = nullptr;
  std::size_t count = 0;
  double *deviceRightHandSide = nullptr;
  bool directAssembly = false;
};

HipLinearSSNFusedAdapter::HipLinearSSNFusedAdapter(
    CPS::EMT::LinearSSNFused &component)
    : mComponent(&component), mDevice(new DeviceState) {
  const auto &models = component.models();
  if (models.empty()) {
    delete mDevice;
    mDevice = nullptr;
    throw std::invalid_argument("HIP fused SSN batch cannot be empty");
  }
  try {
    mDevice->models = hipLinearSSNCreate(models.data(), models.size());
  } catch (...) {
    delete mDevice;
    mDevice = nullptr;
    throw;
  }
  mDevice->count = models.size();
  component.setDeviceManaged(true);
}

HipLinearSSNFusedAdapter::~HipLinearSSNFusedAdapter() {
  if (mComponent) {
    mComponent->setDeviceManaged(false);
    mComponent->setDeviceRhsManaged(false);
  }
  if (mDevice) {
    hipLinearSSNDestroy(mDevice->models);
    delete mDevice;
  }
}

void HipLinearSSNFusedAdapter::updateHistories() {
  if (mDevice->directAssembly) {
    hipLinearSSNUpdateHistoriesDevice(mDevice->models, mDevice->count);
    return;
  }
  std::vector<double> packed(mDevice->count * 15, 0.0);
  hipLinearSSNUpdateHistories(mDevice->models, mDevice->count, packed.data());
  const auto &models = mComponent->models();
  std::vector<CPS::Matrix> histories;
  histories.reserve(models.size());
  for (std::size_t index = 0; index < models.size(); ++index) {
    const std::size_t count =
        models[index].kind == CPS::EMT::LinearSSNModel::Kind::PiLine
            ? 9
            : (models[index].kind ==
                       CPS::EMT::LinearSSNModel::Kind::ImpedanceLoad
                   ? 6
                   : 15);
    CPS::Matrix history(static_cast<Eigen::Index>(count), 1);
    for (std::size_t value = 0; value < count; ++value)
      history(static_cast<Eigen::Index>(value)) = packed[index * 15 + value];
    histories.push_back(std::move(history));
  }
  mComponent->setDeviceHistories(histories);
}

void HipLinearSSNFusedAdapter::updateStates(const CPS::Matrix &leftVector) {
  hipLinearSSNUpdateStates(mDevice->models, mDevice->count, leftVector.data(),
                           static_cast<std::size_t>(leftVector.size()));
}

void HipLinearSSNFusedAdapter::configureDeviceRhs(double *deviceRightHandSide) {
  if (!deviceRightHandSide)
    throw std::invalid_argument("HIP fused SSN RHS buffer is required");
  mDevice->deviceRightHandSide = deviceRightHandSide;
  mDevice->directAssembly = true;
  mComponent->setDeviceRhsManaged(true);
}

void HipLinearSSNFusedAdapter::accumulateDeviceRhs() {
  if (!mDevice->directAssembly || !mDevice->deviceRightHandSide)
    throw std::logic_error("HIP fused SSN device RHS is not configured");
  hipLinearSSNAccumulateRhs(mDevice->models, mDevice->deviceRightHandSide);
}

void HipLinearSSNFusedAdapter::updateStatesDevice(
    const CPS::Matrix &leftVector, const double *deviceSolution) {
  if (!mDevice->directAssembly || !deviceSolution)
    throw std::logic_error("HIP fused SSN device solution is not configured");
  hipLinearSSNUpdateStatesDevice(mDevice->models, mDevice->count,
                                 deviceSolution,
                                 static_cast<std::size_t>(leftVector.size()));
}
} // namespace DPsim
