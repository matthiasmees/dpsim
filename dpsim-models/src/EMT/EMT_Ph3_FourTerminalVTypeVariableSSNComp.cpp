// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#include <dpsim-models/EMT/EMT_Ph3_FourTerminalVTypeVariableSSNComp.h>
#include <dpsim-models/MathUtils.h>

using namespace CPS;

EMT::Ph3::FourTerminalVTypeVariableSSNComp::FourTerminalVTypeVariableSSNComp(
    String uid, String name, Logger::Level logLevel)
    : VTypeVariableSSNComp(uid, name, 6, 6, logLevel),
      mTerminalIncidence(Matrix::Zero(6, 12)) {

  mPhaseType = PhaseType::ABC;
  setTerminalNumber(4);

  const Matrix I = Matrix::Identity(3, 3);

  // u = [v2-v0; v3-v1]
  mTerminalIncidence.block(0, 0, 3, 3) = -I;
  mTerminalIncidence.block(0, 6, 3, 3) = I;

  mTerminalIncidence.block(3, 3, 3, 3) = -I;
  mTerminalIncidence.block(3, 9, 3, 3) = I;
}

const Matrix &
EMT::Ph3::FourTerminalVTypeVariableSSNComp::terminalIncidenceMatrix() const {
  return mTerminalIncidence;
}

MatrixComp
EMT::Ph3::FourTerminalVTypeVariableSSNComp::buildInitialInputFromNodes(Real) {

  MatrixComp u = MatrixComp::Zero(6, 1);

  Complex v20 = RMS3PH_TO_PEAK1PH * initialSingleVoltage(2) -
                RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);

  u(0, 0) = v20;
  u(1, 0) = v20 * SHIFT_TO_PHASE_B;
  u(2, 0) = v20 * SHIFT_TO_PHASE_C;

  Complex v31 = RMS3PH_TO_PEAK1PH * initialSingleVoltage(3) -
                RMS3PH_TO_PEAK1PH * initialSingleVoltage(1);

  u(3, 0) = v31;
  u(4, 0) = v31 * SHIFT_TO_PHASE_B;
  u(5, 0) = v31 * SHIFT_TO_PHASE_C;

  return u;
}

void EMT::Ph3::FourTerminalVTypeVariableSSNComp::mnaCompUpdateVoltage(
    const Matrix &leftVector) {

  Matrix vTerm = Matrix::Zero(12, 1);

  for (UInt t = 0; t < 4; ++t) {
    if (!terminalNotGrounded(t))
      continue;
    for (UInt p = 0; p < 3; ++p)
      vTerm(3 * t + p, 0) =
          Math::realFromVectorElement(leftVector, matrixNodeIndex(t, p));
  }

  **mIntfVoltage = mTerminalIncidence * vTerm;
}

void EMT::Ph3::FourTerminalVTypeVariableSSNComp::mnaCompApplySystemMatrixStamp(
    SparseMatrixRow &systemMatrix) {

  Matrix yTerm = mTerminalIncidence.transpose() * mW * mTerminalIncidence;

  for (UInt tr = 0; tr < 4; ++tr) {
    if (!terminalNotGrounded(tr))
      continue;
    for (UInt pr = 0; pr < 3; ++pr) {
      UInt r = 3 * tr + pr;
      for (UInt tc = 0; tc < 4; ++tc) {
        if (!terminalNotGrounded(tc))
          continue;
        for (UInt pc = 0; pc < 3; ++pc) {
          UInt c = 3 * tc + pc;
          Math::addToMatrixElement(systemMatrix, matrixNodeIndex(tr, pr),
                                   matrixNodeIndex(tc, pc), yTerm(r, c));
        }
      }
    }
  }
}

void EMT::Ph3::FourTerminalVTypeVariableSSNComp::
    mnaCompApplyRightSideVectorStamp(Matrix &rightVector) {

  Matrix iHist = mTerminalIncidence.transpose() * mYHist;

  for (UInt t = 0; t < 4; ++t) {
    if (!terminalNotGrounded(t))
      continue;

    for (UInt p = 0; p < 3; ++p) {
      UInt idx = 3 * t + p;
      Math::setVectorElement(rightVector, matrixNodeIndex(t, p),
                             -iHist(idx, 0));
    }
  }
}
