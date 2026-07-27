// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <dpsim-models/EMT/EMT_VTypeVariableSSNComp.h>

namespace CPS {
namespace EMT {
namespace Ph3 {

class FourTerminalVTypeVariableSSNComp : public VTypeVariableSSNComp {
public:
  FourTerminalVTypeVariableSSNComp(String uid, String name,
                                   Logger::Level logLevel = Logger::Level::off);

  FourTerminalVTypeVariableSSNComp(String name,
                                   Logger::Level logLevel = Logger::Level::off)
      : FourTerminalVTypeVariableSSNComp(name, name, logLevel) {}

protected:
  MatrixComp buildInitialInputFromNodes(Real frequency);
  const Matrix &terminalIncidenceMatrix() const;

  void mnaCompApplySystemMatrixStamp(SparseMatrixRow &systemMatrix) override;
  void mnaCompApplyRightSideVectorStamp(Matrix &rightVector) override;
  void mnaCompUpdateVoltage(const Matrix &leftVector) override;

private:
  Matrix mTerminalIncidence;
};

} // namespace Ph3
} // namespace EMT
} // namespace CPS
