/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cstdint>
#include <memory>
#include <vector>

namespace DPsim {

// Experimental direct solve for CSR MNA input. Apple MPS exposes LU only for
// dense matrices, so this adapter expands the GPU-resident CSR values into a
// dense MPS matrix before factoring it. It is not a sparse LU factorization.
class MetalSparseLUAdapter {
public:
  MetalSparseLUAdapter(void *device, uint32_t order,
                       const std::vector<uint32_t> &rowOffsets,
                       const std::vector<uint32_t> &columnIndices);
  ~MetalSparseLUAdapter();
  MetalSparseLUAdapter(const MetalSparseLUAdapter &) = delete;
  MetalSparseLUAdapter &operator=(const MetalSparseLUAdapter &) = delete;

  void encodeSolve(void *commandBuffer, void *sparseValues, void *rhs,
                   void *solution);
  void checkStatus() const;

private:
  struct Impl;
  std::unique_ptr<Impl> mImpl;
};

} // namespace DPsim
