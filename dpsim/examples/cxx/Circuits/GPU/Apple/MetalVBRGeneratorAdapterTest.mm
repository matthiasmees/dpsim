/* SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
 * SPDX-License-Identifier: MPL-2.0
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <dpsim/GPU/Apple/MetalVBRGeneratorAdapter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using Vec3 = std::array<double, 3>;
using Mat3 = std::array<double, 9>;

Mat3 abcToDq(double theta) {
  Mat3 result{};
  const double angle[3] = {theta, theta - 2. * M_PI / 3.,
                           theta + 2. * M_PI / 3.};
  for (int phase = 0; phase < 3; ++phase) {
    result[phase] = 2. / 3. * std::cos(angle[phase]);
    result[3 + phase] = -2. / 3. * std::sin(angle[phase]);
    result[6 + phase] = 1. / 3.;
  }
  return result;
}

Mat3 dqToAbc(double theta) {
  Mat3 result{};
  const double angle[3] = {theta, theta - 2. * M_PI / 3.,
                           theta + 2. * M_PI / 3.};
  for (int phase = 0; phase < 3; ++phase) {
    result[3 * phase] = std::cos(angle[phase]);
    result[3 * phase + 1] = -std::sin(angle[phase]);
    result[3 * phase + 2] = 1.;
  }
  return result;
}

Mat3 multiply(const Mat3 &a, const Mat3 &b) {
  Mat3 result{};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      for (int k = 0; k < 3; ++k)
        result[3 * r + c] += a[3 * r + k] * b[3 * k + c];
  return result;
}

Vec3 multiply(const Mat3 &a, const Vec3 &x) {
  Vec3 result{};
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      result[r] += a[3 * r + c] * x[c];
  return result;
}

void requireClose(double actual, double expected, const char *label) {
  if (std::abs(actual - expected) > 2e-5 * std::max(1., std::abs(expected))) {
    std::cerr << label << " actual=" << actual << " expected=" << expected << '\n';
    throw std::runtime_error("Metal VBR generator kernel differs from reference");
  }
}
} // namespace

int main() {
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  try {
    using Model = DPsim::MetalVBRGeneratorAdapter::Model;
    Model::MetalDeviceData model{};
    model.voltageDq0 = {1., .2, 0.};
    model.currentDq0 = {.1, .2, 0.};
    model.transientVoltageDq0 = {.8, .7, 0.};
    model.conductanceDq0 = {1., 0., 0., 0., 2., 0., 0., 0., 3.};
    model.matrixNodeIndices = {0, 1, 2};
    model.timeStep = .001;
    model.baseVoltage = 1.;
    model.baseCurrent = 1.;
    model.baseImpedance = 1.;
    model.baseOmegaElec = 2. * M_PI * 60.;
    model.baseOmegaMech = 1.;
    model.inertia = 4.;
    model.damping = 1.;
    model.ldTransient = .2;
    model.lqTransient = .3;
    model.adTransient = .01;
    model.bdTransient = .1;
    model.aqTransient = .02;
    model.bqTransient = .2;
    model.dqTransient = .03;
    model.fieldVoltage = .4;
    model.previousFieldVoltage = .4;
    model.mechanicalTorque = .5;
    model.previousMechanicalTorque = .5;
    model.electricalTorque = .2;
    model.rotorSpeed = 1.;
    model.rotorAngle = .1;
    model.loadAngle = .2;
    model.conductanceAbc = multiply(
        multiply(dqToAbc(model.rotorAngle), model.conductanceDq0),
        abcToDq(model.rotorAngle));

    DPsim::MetalVBRGeneratorAdapter::StampPositions slots{};
    for (uint32_t i = 0; i < 9; ++i) slots.positions[i] = i;
    DPsim::MetalVBRGeneratorAdapter adapter({model}, {slots});

    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    id<MTLBuffer> matrix = [device newBufferWithLength:9 * sizeof(float)
                                               options:MTLResourceStorageModeShared];
    id<MTLBuffer> rhs = [device newBufferWithLength:3 * sizeof(float)
                                             options:MTLResourceStorageModeShared];
    id<MTLBuffer> solution = [device newBufferWithLength:3 * sizeof(float)
                                                  options:MTLResourceStorageModeShared];
    if (!device || !matrix || !rhs || !solution)
      throw std::runtime_error("test Metal allocation failed");
    std::fill_n(static_cast<float *>(matrix.contents), 9, 0.f);
    std::fill_n(static_cast<float *>(rhs.contents), 3, 0.f);
    const float vabc[3] = {1.1f, -.3f, .2f};
    std::copy(vabc, vabc + 3, static_cast<float *>(solution.contents));

    const double dt = model.timeStep;
    const double torque = model.voltageDq0[0] * model.currentDq0[0] +
                          model.voltageDq0[1] * model.currentDq0[1];
    const double speed = model.rotorSpeed + dt / (2. * model.inertia) *
        (model.previousMechanicalTorque - torque -
         model.damping * (model.rotorSpeed - 1.));
    const double angle = model.rotorAngle + dt * speed * model.baseOmegaMech;
    const double delta = model.loadAngle + dt * (speed - 1.) * model.baseOmegaMech;
    const Vec3 edq{model.voltageDq0[0] - model.currentDq0[1] * model.lqTransient,
                   model.voltageDq0[1] + model.currentDq0[0] * model.ldTransient,
                   0.};
    const Mat3 toAbc = dqToAbc(angle);
    const Mat3 toDq = abcToDq(angle);
    const Mat3 conductanceAbc =
        multiply(multiply(toAbc, model.conductanceDq0), toDq);
    const Vec3 edqHistory{
        model.adTransient * model.currentDq0[1] + model.bdTransient * edq[0],
        model.aqTransient * model.currentDq0[0] + model.bqTransient * edq[1] +
            model.dqTransient * (model.previousFieldVoltage + model.fieldVoltage),
        0.};
    Vec3 eabc = multiply(toAbc, edqHistory);
    for (double &value : eabc) value *= model.baseVoltage;
    const Vec3 source = multiply(conductanceAbc, eabc);
    const Vec3 v{vabc[0], vabc[1], vabc[2]};
    const Vec3 drop = multiply(conductanceAbc, v);
    Vec3 iabc{};
    for (int i = 0; i < 3; ++i) iabc[i] = source[i] - drop[i];
    const Vec3 expectedVoltageDq = multiply(toDq, v);
    const Vec3 expectedCurrentDq = multiply(toDq, iabc);

    id<MTLCommandQueue> queue = (id<MTLCommandQueue>)adapter.commandQueue();
    id<MTLCommandBuffer> command = [queue commandBuffer];
    adapter.encodePreStep(command, matrix, rhs, dt);
    adapter.encodePostStep(command, solution);
    [command commit];
    [command waitUntilCompleted];
    if (command.error)
      throw std::runtime_error(command.error.localizedDescription.UTF8String);

    const auto *matrixActual = static_cast<const float *>(matrix.contents);
    const auto *rhsActual = static_cast<const float *>(rhs.contents);
    for (int i = 0; i < 9; ++i)
      requireClose(matrixActual[i], conductanceAbc[i], "Gabc");
    for (int i = 0; i < 3; ++i)
      requireClose(rhsActual[i], source[i], "Norton source");
    const auto state = adapter.copyStateToHost().front();
    for (int i = 0; i < 3; ++i) {
      requireClose(state.voltageDq0[i], expectedVoltageDq[i], "Vdq0");
      requireClose(state.currentDq0[i], expectedCurrentDq[i], "Idq0");
    }
    requireClose(state.electricalTorque, torque, "electrical torque");
    requireClose(state.rotorSpeed, speed, "rotor speed");
    requireClose(state.rotorAngle, angle, "rotor angle");
    requireClose(state.loadAngle, delta, "load angle");
    for (int i = 0; i < 3; ++i)
      requireClose(state.transientVoltageDq0[i], edq[i], "transient voltage");

    std::cout << "Metal SG4 VBR pre/post-step smoke test passed\n";
    [matrix release]; [rhs release]; [solution release]; [device release];
    [pool drain];
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    [pool drain];
    return 1;
  }
}
