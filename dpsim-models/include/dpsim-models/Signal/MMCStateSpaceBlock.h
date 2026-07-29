// SPDX-FileCopyrightText: 2026 Institute for Automation of Complex Power Systems, EONERC, RWTH Aachen University
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <dpsim-models/Definitions.h>

namespace CPS {
namespace Signal {

/// Dense affine state-space representation evaluated at one operating point:
///
///   x_dot = A*x + B*u + E
///   y     = C*x + D*u + F
struct MMCLinearization {
  Matrix A;
  Matrix B;
  Matrix C;
  Matrix D;
  Matrix E;
  Matrix F;
};

/// Sparse counterpart used for block composition and diagnostics.
/// E and F remain dense column vectors.
struct MMCSparseLinearization {
  SparseMatrixRow A;
  SparseMatrixRow B;
  SparseMatrixRow C;
  SparseMatrixRow D;
  Matrix E;
  Matrix F;
};

/// Base interface for every MMC controller and signal-processing block.
///
/// Every derived class exposes the complete nonlinear equations and its own
/// local affine state-space model. The derived class owns no DPsim integration
/// task: its state is a slice of the enclosing MMC or HVDCLink state vector.
class MMCStateSpaceBlock {
public:
  virtual ~MMCStateSpaceBlock() = default;

  virtual UInt stateSize() const = 0;
  virtual UInt inputSize() const = 0;
  virtual UInt outputSize() const = 0;
  virtual std::vector<String> stateNames() const = 0;

  /// Evaluate the complete nonlinear block equations in one call.
  virtual void evaluate(const Matrix &x, const Matrix &u, Matrix &dx,
                        Matrix &y) const = 0;

  /// Explicit equation accessors used by analysis, testing, and composition.
  virtual void evaluateStateDerivative(const Matrix &x, const Matrix &u,
                                       Matrix &dx) const = 0;
  virtual void evaluateOutput(const Matrix &x, const Matrix &u,
                              Matrix &y) const = 0;

  /// Calculate local A, B, C, D Jacobians at (x,u).
  virtual void calculateNumericalJacobians(const Matrix &x, const Matrix &u,
                                           Matrix &A, Matrix &B, Matrix &C,
                                           Matrix &D, Real relativeStep = 1e-6,
                                           Real absoluteStep = 1e-8) const = 0;

  /// Build the complete local affine model at (x,u).
  virtual void buildStateSpaceModel(const Matrix &x, const Matrix &u, Matrix &A,
                                    Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                                    Matrix &F, Real relativeStep = 1e-6,
                                    Real absoluteStep = 1e-8) const = 0;

  /// Return all six dense matrices as one value.
  virtual MMCLinearization
  getStateSpaceModel(const Matrix &x, const Matrix &u, Real relativeStep = 1e-6,
                     Real absoluteStep = 1e-8) const = 0;

  /// Return a sparse A/B/C/D representation for block assembly.
  virtual MMCSparseLinearization
  getSparseStateSpaceModel(const Matrix &x, const Matrix &u,
                           Real relativeStep = 1e-6, Real absoluteStep = 1e-8,
                           Real sparseTolerance = 1e-14) const = 0;

protected:
  void validateDimensions(const Matrix &x, const Matrix &u) const {
    if (x.rows() != static_cast<Eigen::Index>(stateSize()) || x.cols() != 1)
      throw std::invalid_argument(
          "Signal block state vector has invalid dimensions.");
    if (u.rows() != static_cast<Eigen::Index>(inputSize()) || u.cols() != 1)
      throw std::invalid_argument(
          "Signal block input vector has invalid dimensions.");
    if (!x.allFinite() || !u.allFinite())
      throw std::invalid_argument("Signal block received NaN or Inf.");
  }

  void validateLinearizationSteps(Real relativeStep, Real absoluteStep) const {
    if (!std::isfinite(relativeStep) || relativeStep <= 0.0 ||
        !std::isfinite(absoluteStep) || absoluteStep <= 0.0)
      throw std::invalid_argument(
          "State-space linearization steps must be finite and positive.");
  }

  /// Shared central-difference implementation. Each concrete class exposes
  /// this operation through its own calculateNumericalJacobians() method.
  void calculateNumericalJacobiansGeneric(const Matrix &x, const Matrix &u,
                                          Matrix &A, Matrix &B, Matrix &C,
                                          Matrix &D, Real relativeStep,
                                          Real absoluteStep) const {
    validateDimensions(x, u);
    validateLinearizationSteps(relativeStep, absoluteStep);

    A = Matrix::Zero(stateSize(), stateSize());
    B = Matrix::Zero(stateSize(), inputSize());
    C = Matrix::Zero(outputSize(), stateSize());
    D = Matrix::Zero(outputSize(), inputSize());

    Matrix dxPlus = Matrix::Zero(stateSize(), 1);
    Matrix dxMinus = Matrix::Zero(stateSize(), 1);
    Matrix yPlus = Matrix::Zero(outputSize(), 1);
    Matrix yMinus = Matrix::Zero(outputSize(), 1);

    for (UInt column = 0; column < stateSize(); ++column) {
      const Real step =
          absoluteStep + relativeStep * std::max(1.0, std::abs(x(column, 0)));
      Matrix xPlus = x;
      Matrix xMinus = x;
      xPlus(column, 0) += step;
      xMinus(column, 0) -= step;

      evaluate(xPlus, u, dxPlus, yPlus);
      evaluate(xMinus, u, dxMinus, yMinus);
      A.col(column) = (dxPlus - dxMinus) / (2.0 * step);
      C.col(column) = (yPlus - yMinus) / (2.0 * step);
    }

    for (UInt column = 0; column < inputSize(); ++column) {
      const Real step =
          absoluteStep + relativeStep * std::max(1.0, std::abs(u(column, 0)));
      Matrix uPlus = u;
      Matrix uMinus = u;
      uPlus(column, 0) += step;
      uMinus(column, 0) -= step;

      evaluate(x, uPlus, dxPlus, yPlus);
      evaluate(x, uMinus, dxMinus, yMinus);
      B.col(column) = (dxPlus - dxMinus) / (2.0 * step);
      D.col(column) = (yPlus - yMinus) / (2.0 * step);
    }
  }

  /// Shared affine-offset implementation matching the requested DPsim style.
  void buildStateSpaceModelGeneric(const Matrix &x, const Matrix &u, Matrix &A,
                                   Matrix &B, Matrix &C, Matrix &D, Matrix &E,
                                   Matrix &F, Real relativeStep,
                                   Real absoluteStep) const {
    calculateNumericalJacobians(x, u, A, B, C, D, relativeStep, absoluteStep);

    Matrix stateDerivative = Matrix::Zero(stateSize(), 1);
    Matrix output = Matrix::Zero(outputSize(), 1);
    evaluateStateDerivative(x, u, stateDerivative);
    evaluateOutput(x, u, output);

    // Local affine offsets:
    //
    //   E = f(x0,u0) - A*x0 - B*u0
    //   F = g(x0,u0) - C*x0 - D*u0
    E = stateDerivative - A * x - B * u;
    F = output - C * x - D * u;
  }

  MMCLinearization getStateSpaceModelGeneric(const Matrix &x, const Matrix &u,
                                             Real relativeStep,
                                             Real absoluteStep) const {
    MMCLinearization model;
    buildStateSpaceModel(x, u, model.A, model.B, model.C, model.D, model.E,
                         model.F, relativeStep, absoluteStep);
    return model;
  }

  MMCSparseLinearization
  getSparseStateSpaceModelGeneric(const Matrix &x, const Matrix &u,
                                  Real relativeStep, Real absoluteStep,
                                  Real sparseTolerance) const {
    if (!std::isfinite(sparseTolerance) || sparseTolerance < 0.0)
      throw std::invalid_argument(
          "Sparse state-space tolerance must be finite and non-negative.");

    const MMCLinearization dense =
        getStateSpaceModel(x, u, relativeStep, absoluteStep);
    MMCSparseLinearization sparse;
    sparse.A = dense.A.sparseView(0.0, sparseTolerance);
    sparse.B = dense.B.sparseView(0.0, sparseTolerance);
    sparse.C = dense.C.sparseView(0.0, sparseTolerance);
    sparse.D = dense.D.sparseView(0.0, sparseTolerance);
    sparse.E = dense.E;
    sparse.F = dense.F;
    return sparse;
  }
};

/// Maps one local block into a contiguous slice of the enclosing state vector.
struct MMCStateSlice {
  UInt first = 0;
  UInt count = 0;

  Matrix get(const Matrix &global) const {
    if (count == 0)
      return Matrix::Zero(0, 1);
    return global.block(first, 0, count, 1);
  }

  void set(Matrix &global, const Matrix &local) const {
    if (local.rows() != static_cast<Eigen::Index>(count) || local.cols() != 1)
      throw std::invalid_argument(
          "Local state does not match assigned state slice.");
    if (count > 0)
      global.block(first, 0, count, 1) = local;
  }
};

} // namespace Signal
} // namespace CPS
