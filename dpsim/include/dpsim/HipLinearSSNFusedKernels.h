/* SPDX-License-Identifier: MPL-2.0 */
#pragma once

#include <cstddef>

#include <dpsim-models/EMT/EMT_LinearSSNModel.h>

namespace DPsim {
void *hipLinearSSNCreate(const CPS::EMT::LinearSSNModel *models,
                         std::size_t modelCount);
void hipLinearSSNUpdateHistories(void *deviceModels, std::size_t modelCount,
                                 double *hostHistoryRows);
void hipLinearSSNUpdateHistoriesDevice(void *deviceModels,
                                       std::size_t modelCount);
void hipLinearSSNAccumulateRhs(void *deviceModels, double *deviceRightHandSide);
void hipLinearSSNUpdateStates(void *deviceModels, std::size_t modelCount,
                              const double *hostSolution,
                              std::size_t solutionSize);
void hipLinearSSNUpdateStatesDevice(void *deviceModels, std::size_t modelCount,
                                    const double *deviceSolution,
                                    std::size_t solutionSize);
void hipLinearSSNDestroy(void *deviceModels) noexcept;
} // namespace DPsim
