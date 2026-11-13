/* SPDX-License-Identifier: MPL-2.0 */
#pragma once

#include <cstddef>

#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator4OrderVBR.h>

namespace DPsim {
void *hipSG4BatchCreate(const CPS::EMT::Ph3::HipSG4State *states,
                        std::size_t count);
void hipSG4BatchPreStep(void *batch, double time,
                        CPS::EMT::Ph3::HipSG4State *hostStates);
void hipSG4BatchScatterMatrix(void *batch);
void hipSG4BatchConfigureAssembly(void *batch, const int *hostMatrixPositions,
                                  double *deviceMatrixValues,
                                  double *deviceRightHandSide,
                                  const int *hostRestorePositions,
                                  const double *hostRestoreValues,
                                  std::size_t restoreCount);
void hipSG4BatchAccumulateRhs(void *batch);
void hipSG4BatchPostStep(void *batch, const double *hostSolution,
                         std::size_t solutionSize,
                         CPS::EMT::Ph3::HipSG4State *hostStates);
void hipSG4BatchPostStepDevice(void *batch, const double *deviceSolution,
                               std::size_t solutionSize,
                               CPS::EMT::Ph3::HipSG4State *hostStates);
void hipSG4BatchDestroy(void *batch) noexcept;
} // namespace DPsim
