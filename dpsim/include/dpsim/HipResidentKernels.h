/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <cstddef>

#include <hip/hip_runtime_api.h>

namespace DPsim {
/// Clear a device RHS vector and accumulate indexed component contributions.
/// All pointers must refer to device memory. The operation is enqueued on the
/// supplied stream and does not synchronize with the host.
void hipClearAndAccumulateIndexedRhs(double *deviceRightHandSide,
                                     std::size_t rowCount,
                                     const double *deviceContributions,
                                     const int *deviceRows,
                                     std::size_t contributionCount,
                                     hipStream_t stream = nullptr);

/// Multiply a device vector elementwise by a device scaling vector.
void hipMultiplyVectorInPlace(double *deviceVector, const double *deviceScale,
                              std::size_t elementCount,
                              hipStream_t stream = nullptr);

/// Update selected values in a device-resident CSR value array.
void hipScatterSparseValues(double *deviceValues, const int *devicePositions,
                            const double *deviceUpdates,
                            std::size_t updateCount,
                            hipStream_t stream = nullptr);

/// Execute the complete preconditioned BiCGSTAB iteration on one HIP kernel
/// launch. The status is a device integer: positive means converged, zero
/// means the iteration limit was reached, and negative means breakdown.
void hipSolveBiCGStabResident(
    const int *deviceRowPointers, const int *deviceColumnIndices,
    const double *deviceValues, const double *deviceRowScale,
    const double *deviceRightHandSide, double *deviceSolution,
    double *deviceResidual, double *deviceShadowResidual,
    double *deviceDirection, double *deviceProduct, double *deviceIntermediate,
    double *deviceIntermediateProduct, int dimension, int maximumIterations,
    double relativeTolerance, int *deviceStatus, hipStream_t stream = nullptr);
} // namespace DPsim
