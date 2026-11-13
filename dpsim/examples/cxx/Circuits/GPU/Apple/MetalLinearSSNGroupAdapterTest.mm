/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <dpsim/GPU/Apple/MetalLinearSSNGroupAdapter.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  try {
    using namespace CPS;
    EMT::LinearSSNGroup::PrecomputedStateMaps model;
    model.historyMap.resize(2, 3);
    model.historyMap.insert(0, 0) = 1.;
    model.historyMap.insert(0, 1) = 2.;
    model.historyMap.insert(1, 1) = -1.;
    model.historyMap.insert(1, 2) = 1.;
    model.historyMap.makeCompressed();

    model.stateMap.resize(3, 3);
    model.stateMap.insert(0, 0) = .5;
    model.stateMap.insert(1, 0) = 1.;
    model.stateMap.insert(1, 1) = .25;
    model.stateMap.insert(2, 1) = 1.;
    model.stateMap.insert(2, 2) = .5;
    model.stateMap.makeCompressed();

    model.inputMap.resize(3, 2);
    model.inputMap.insert(0, 0) = 1.;
    model.inputMap.insert(1, 1) = 2.;
    model.inputMap.insert(2, 0) = 1.;
    model.inputMap.insert(2, 1) = -1.;
    model.inputMap.makeCompressed();

    model.terminalAdmittance = Matrix::Identity(2, 2);
    model.initialState.resize(3, 1);
    model.initialState << 2., 3., 4.;
    model.terminalMatrixIndices = {1, 4};

    DPsim::MetalLinearSSNGroupAdapter adapter(model);
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    id<MTLBuffer> rhs = [device newBufferWithLength:5 * sizeof(float)
                                             options:MTLResourceStorageModeShared];
    id<MTLBuffer> rhs2 = [device newBufferWithLength:5 * sizeof(float)
                                              options:MTLResourceStorageModeShared];
    id<MTLBuffer> solution = [device newBufferWithLength:5 * sizeof(float)
                                                  options:MTLResourceStorageModeShared];
    id<MTLBuffer> solution2 = [device newBufferWithLength:5 * sizeof(float)
                                                   options:MTLResourceStorageModeShared];
    if (!rhs || !rhs2 || !solution || !solution2)
      throw std::runtime_error("test buffer allocation failed");
    std::fill_n(static_cast<float *>(rhs.contents), 5, 0.f);
    std::fill_n(static_cast<float *>(rhs2.contents), 5, 0.f);
    const float input[5] = {0.f, 5.f, 0.f, 0.f, 7.f};
    std::copy(input, input + 5, static_cast<float *>(solution.contents));
    const float input2[5] = {0.f, 2.f, 0.f, 0.f, 1.f};
    std::copy(input2, input2 + 5, static_cast<float *>(solution2.contents));

    id<MTLCommandQueue> queue = (id<MTLCommandQueue>)adapter.commandQueue();
    id<MTLCommandBuffer> command = [queue commandBuffer];
    adapter.encodeHistoryAndRhs(command, rhs);
    adapter.encodeStateUpdate(command, solution);
    adapter.encodeHistoryAndRhs(command, rhs2);
    adapter.encodeStateUpdate(command, solution2);
    [command commit];
    [command waitUntilCompleted];
    if (command.error)
      throw std::runtime_error(command.error.localizedDescription.UTF8String);

    const float expectedRhs[5] = {0.f, 8.f, 0.f, 0.f, 1.f};
    const float expectedRhs2[5] = {0.f, 39.5f, 0.f, 0.f, -13.75f};
    const float expectedState[3] = {5.f, 12.1875f, 19.25f};
    const auto *actualRhs = static_cast<const float *>(rhs.contents);
    const auto *actualRhs2 = static_cast<const float *>(rhs2.contents);
    const auto actualState = adapter.copyStateToHost();
    for (int i = 0; i < 5; ++i)
      if (std::abs(actualRhs[i] - expectedRhs[i]) > 1e-5f)
        throw std::runtime_error("GPU group history/RHS differs from reference");
    for (int i = 0; i < 5; ++i)
      if (std::abs(actualRhs2[i] - expectedRhs2[i]) > 1e-5f)
        throw std::runtime_error("GPU resident second-step history differs from reference");
    for (int i = 0; i < 3; ++i)
      if (std::abs(actualState[i] - expectedState[i]) > 1e-5f) {
        std::cerr << "state[" << i << "] actual=" << actualState[i]
                  << " expected=" << expectedState[i] << '\n';
        throw std::runtime_error("GPU group state update differs from reference");
      }

    std::cout << "Metal linear-group history/RHS and state-update smoke test passed\n";
    [rhs release];
    [rhs2 release];
    [solution release];
    [solution2 release];
    [device release];
    [pool drain];
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    [pool drain];
    return 1;
  }
}
