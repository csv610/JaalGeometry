///////////////////////////////////////////////////////////////////////////////
// This file is part of ShapeOp, a lightweight C++ library
// for static and dynamic geometry processing.
//
// Copyright (C) 2014 Sofien Bouaziz <sofien.bouaziz@gmail.com>
// Copyright (C) 2014 LGG EPFL
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
///////////////////////////////////////////////////////////////////////////////

#pragma once

///////////////////////////////////////////////////////////////////////////////
/** \file
This file contains all the linear system solvers of the ShapeOp library.*/
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <unsupported/Eigen/src/IterativeSolvers/MINRES.h>

namespace ShapeOp {
///////////////////////////////////////////////////////////////////////////////
/** \brief Base class of any sparse linear system solver. This class defines the main functionalities of the ShapeOp sparse linear system solvers (Ax = b).*/
class  LSSolver {
 public:
  virtual ~LSSolver() {};
  /** \brief Initialize the linear system solver using the sparse matrix A.*/
  virtual void initialize(const SparseMatrix &A, unsigned int iteration = 1) = 0;
  /** \brief Solve the linear system Ax = b.*/
  virtual VectorX solve(const VectorX &b, const VectorX &x0) const = 0;
  /** \brief Reports whether previous computation was successful.*/
  virtual Eigen::ComputationInfo info() const = 0;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Sparse linear system solver based on Cholesky. This class implements a sparse linear system solver based on the Cholesky LDL^T algorithm from Eigen.*/
class  SimplicialLDLTSolver : public LSSolver {
 public:
  virtual ~SimplicialLDLTSolver() {};
  /** \brief Prefactorize the sparse matrix (A = LDL^T).*/
  virtual void initialize(const SparseMatrix &A, unsigned int iteration) override final;
  /** \brief Solve the linear system by applying twice backsubstitution.*/
  virtual VectorX solve(const VectorX &b, const VectorX &x0) const override final;
  /** \brief Reports whether previous computation was successful.*/
  virtual Eigen::ComputationInfo info() const override final;
 private:
  Eigen::SimplicialLDLT<SparseMatrix> solver_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Sparse linear system solver based on CG. This class implements a sparse linear system solver based on the CG algorithm from Eigen.*/
class  CGSolver : public LSSolver {
 public:
  virtual ~CGSolver() {};
  /** \brief Initialize PCG.*/
  virtual void initialize(const SparseMatrix &A, unsigned int iteration) override final;
  /** \brief Solve the linear system by applying CG.*/
  virtual VectorX solve(const VectorX &b, const VectorX &x0) const override final;
  /** \brief Reports whether previous computation was successful.*/
  virtual Eigen::ComputationInfo info() const override final;
 private:
  Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower, Eigen::IncompleteLUT<Scalar> > solver_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Sparse linear system solver based on MINRES. This class implements a sparse linear system solver based on the MINRES algorithm from Eigen.*/
class  MINRESSolver : public LSSolver {
 public:
  virtual ~MINRESSolver() {};
  /** \brief Initialize MINRES.*/
  virtual void initialize(const SparseMatrix &A, unsigned int iteration) override final;
  /** \brief Solve the linear system by applying MINRES.*/
  virtual VectorX solve(const VectorX &b, const VectorX &x0) const override final;
  /** \brief Reports whether previous computation was successful.*/
  virtual Eigen::ComputationInfo info() const override final;
 private:
  Eigen::MINRES<SparseMatrix, Eigen::Lower, Eigen::IncompleteLUT<Scalar> > solver_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Sparse linear system solver based on successive over-relaxation (SOR).*/
class  SORSolver : public LSSolver {
 public:
  /** \brief Solver Constructor.
   \param relaxation The relaxation factor.
  */
  SORSolver(ShapeOp::Scalar relaxation = 1.6);
  virtual ~SORSolver() {};
  /** \brief Initialize SOR.*/
  virtual void initialize(const SparseMatrix &A, unsigned int iteration) override final;
  /** \brief Solve the linear system by applying SOR.*/
  virtual VectorX solve(const VectorX &b, const VectorX &x0) const override final;
  /** \brief Reports whether previous computation was successful.*/
  virtual Eigen::ComputationInfo info() const override final;
 private:
  SparseMatrixT<Eigen::RowMajor> A_;
  ShapeOp::Scalar relaxation_;
  unsigned int iteration_;
};
///////////////////////////////////////////////////////////////////////////////
