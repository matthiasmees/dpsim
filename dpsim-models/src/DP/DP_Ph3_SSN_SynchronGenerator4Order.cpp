// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0
#include <dpsim-models/DP/DP_Ph3_SSN_SynchronGenerator4Order.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace CPS;
namespace {
constexpr Int kEd = 0, kEq = 1, kOmega = 2, kTheta = 3;
const Complex alpha = std::polar(1., 2. * PI / 3.);
MatrixComp unpack3(const Matrix &u) {
  MatrixComp v(3, 1);
  for (Int p = 0; p < 3; ++p)
    v(p, 0) = Complex(u(2 * p, 0), u(2 * p + 1, 0));
  return v;
}
Matrix pack3(const MatrixComp &y) {
  Matrix out(6, 1);
  for (Int p = 0; p < 3; ++p) {
    out(2 * p, 0) = y(p, 0).real();
    out(2 * p + 1, 0) = y(p, 0).imag();
  }
  return out;
}
} // namespace

DP::Ph3::SSN::SynchronGenerator4Order::SynchronGenerator4Order(
    String uid, String name, Logger::Level level)
    : MixedVTypeVariableSSNComp(uid, name, 4, 0, level) {
  **mX = Matrix::Zero(4, 1);
  setParameters(Matrix::Zero(4, 4), Matrix::Zero(4, 6), Matrix::Zero(6, 4),
                Matrix::Zero(6, 6));
}

void DP::Ph3::SSN::SynchronGenerator4Order::setParametersOperationalPerUnit(
    Real nominalPower, Real nominalVoltage, Real nominalFrequency, Real inertia,
    Real resistance, Real xdTransient, Real xqTransient, Real xd, Real xq,
    Real td0Transient, Real tq0Transient, Real mechanicalDamping) {
  if (nominalPower <= 0. || nominalVoltage <= 0. || nominalFrequency <= 0. ||
      inertia <= 0. || td0Transient <= 0. || tq0Transient <= 0. ||
      xdTransient <= 0. || xqTransient <= 0.)
    throw std::invalid_argument("Invalid fourth-order SSN machine parameters.");
  mNominalPower = nominalPower;
  mNominalVoltage = nominalVoltage;
  mFrequency = nominalFrequency;
  mBaseV = nominalVoltage * 0.8164965809277260327;
  mBaseI = 2. * nominalPower / (3. * mBaseV);
  mBaseZ = mBaseV / mBaseI;
  mOmega = 2. * PI * nominalFrequency;
  mH = inertia;
  mRa = resistance;
  mXd = xd;
  mXq = xq;
  mXdp = xdTransient;
  mXqp = xqTransient;
  mTd0p = td0Transient;
  mTq0p = tq0Transient;
  mDamping = mechanicalDamping;
  mConfigured = true;
}

void DP::Ph3::SSN::SynchronGenerator4Order::setPowerflowInitialValues(
    Real activePowerPu, Real reactivePowerPu, Real voltageMagnitudePu,
    Real voltageAngleRad) {
  mPm = activePowerPu;
  mQm = reactivePowerPu;
  mVpu = voltageMagnitudePu;
  mVangle = voltageAngleRad;
  mInitialValuesConfigured = true;
}

void DP::Ph3::SSN::SynchronGenerator4Order::evaluate(const Matrix &x,
                                                     const Matrix &u, Matrix &f,
                                                     Matrix &y) const {
  const MatrixComp phases = unpack3(u);
  const Complex v =
      (phases(0, 0) + alpha * phases(1, 0) + std::conj(alpha) * phases(2, 0)) /
      (3. * mBaseV);
  const Complex vdq = v * std::polar(1., -x(kTheta, 0));
  const Real vd = vdq.real(), vq = vdq.imag();
  const Real dd = vd - x(kEd, 0), dq = vq - x(kEq, 0),
             det = mRa * mRa + mXdp * mXqp;
  const Real id = (mRa * dd + mXqp * dq) / det,
             iq = (-mXdp * dd + mRa * dq) / det;
  const Complex idq(id, iq), iph = mBaseI * idq * std::polar(1., x(kTheta, 0));
  MatrixComp currents(3, 1);
  currents(0, 0) = iph;
  currents(1, 0) = iph * std::conj(alpha);
  currents(2, 0) = iph * alpha;
  y = pack3(currents);
  f.setZero(4, 1);
  f(kEd, 0) = (-(x(kEd, 0) - mEd0) - (mXq - mXqp) * (iq - mIq0)) / mTq0p;
  f(kEq, 0) = (-(x(kEq, 0) - mEq0) + (mXd - mXdp) * (id - mId0)) / mTd0p;
  const Real pe = vd * id + vq * iq;
  f(kOmega, 0) = (mPm + pe - mDamping * (x(kOmega, 0) - 1.)) / (2. * mH);
  f(kTheta, 0) = mOmega * (x(kOmega, 0) - 1.);
}

void DP::Ph3::SSN::SynchronGenerator4Order::linearize(const Matrix &x,
                                                      const Matrix &u,
                                                      Matrix &A, Matrix &B,
                                                      Matrix &C,
                                                      Matrix &D) const {
  A.setZero(4, 4);
  B.setZero(4, 6);
  C.setZero(6, 4);
  D.setZero(6, 6);
  Matrix xp = x, xm = x, up = u, um = u, fp, fm, yp, ym;
  for (Int j = 0; j < 4; ++j) {
    const Real h = 1e-6 * std::max(1., std::abs(x(j, 0)));
    xp = x;
    xm = x;
    xp(j, 0) += h;
    xm(j, 0) -= h;
    evaluate(xp, u, fp, yp);
    evaluate(xm, u, fm, ym);
    A.col(j) = (fp - fm) / (2 * h);
    C.col(j) = (yp - ym) / (2 * h);
  }
  for (Int j = 0; j < 6; ++j) {
    const Real h = 1e-6 * std::max(1., std::abs(u(j, 0)));
    up = u;
    um = u;
    up(j, 0) += h;
    um(j, 0) -= h;
    evaluate(x, up, fp, yp);
    evaluate(x, um, fm, ym);
    B.col(j) = (fp - fm) / (2 * h);
    D.col(j) = (yp - ym) / (2 * h);
  }
}

Bool DP::Ph3::SSN::SynchronGenerator4Order::updateComponentParameters() {
  linearize(**mX, pack3(**mIntfVoltage), mA, mB, mC, mD);
  Matrix f, y;
  evaluate(**mX, pack3(**mIntfVoltage), f, y);
  setStateOffset(f - mA * (**mX) - mB * pack3(**mIntfVoltage));
  setOutputOffset(y - mC * (**mX) - mD * pack3(**mIntfVoltage));
  return true;
}

void DP::Ph3::SSN::SynchronGenerator4Order::initializeState(
    const MatrixComp &phaseVoltage) {
  if (!mConfigured || !mInitialValuesConfigured)
    throw std::logic_error(
        "Configure SSN generator parameters and PF values first.");
  const Complex v = mVpu * std::polar(1., mVangle);
  const Complex i = -std::conj(Complex(mPm, mQm) / v);
  const Complex e = v - Complex(mRa, .5 * (mXdp + mXqp)) * i;
  const Real theta = std::arg(e) - PI / 2.;
  const Complex rot = std::polar(1., -theta);
  const Complex vdq = v * rot, idq = i * rot;
  (**mX)(kTheta, 0) = theta;
  (**mX)(kOmega, 0) = 1.;
  (**mX)(kEd, 0) = vdq.real() - mRa * idq.real() + mXqp * idq.imag();
  (**mX)(kEq, 0) = vdq.imag() - mXdp * idq.real() - mRa * idq.imag();
  mEd0 = (**mX)(kEd, 0);
  mEq0 = (**mX)(kEq, 0);
  mId0 = idq.real();
  mIq0 = idq.imag();
  (void)phaseVoltage;
}

void DP::Ph3::SSN::SynchronGenerator4Order::initializeFromNodesAndTerminals(
    Real) {
  if (!mConfigured)
    throw std::logic_error("SSN generator parameters are unset.");
  MatrixComp voltage = MatrixComp::Zero(3, 1);
  voltage(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1);
  voltage(1, 0) = voltage(0, 0) * SHIFT_TO_PHASE_B;
  voltage(2, 0) = voltage(0, 0) * SHIFT_TO_PHASE_C;
  // SystemTopology transfers solved SG power onto terminal 0 for the matching
  // power-flow component name, including the slack-bus active-power result.
  const Complex solvedPower = -terminal(0)->singlePower();
  if (std::isfinite(solvedPower.real()) && std::isfinite(solvedPower.imag()) &&
      std::abs(solvedPower) > 1.) {
    mPm = solvedPower.real() / mNominalPower;
    mQm = solvedPower.imag() / mNominalPower;
  }
  mVpu = std::abs(voltage(0, 0)) / mBaseV;
  mVangle = std::arg(voltage(0, 0));
  **mIntfVoltage = voltage;
  initializeState(voltage);
  updateComponentParameters();
  Matrix f, y;
  evaluate(**mX, pack3(voltage), f, y);
  **mIntfCurrent = unpack3(y);
}

void DP::Ph3::SSN::SynchronGenerator4Order::updateLogAttributes(
    const Matrix &) const {}
Matrix DP::Ph3::SSN::SynchronGenerator4Order::getState() const { return **mX; }
