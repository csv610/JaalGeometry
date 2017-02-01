/*
 * A Simple Geometric Model for Elastic Deformations
 * Isaac Chao, Ulrich Pinkall, Patrick Sanan, Peter Schröder
 * ACM Transactions on Graphics, 29(4), 2010, 38:1-38:6.
 *
 * TODO: implemented just a gradient descent; use full Newton's method
 *
 */

#pragma once

#include "Mesh.hpp"

#include "DDG_Mesh.hpp"
#include "DDG_Complex.hpp"
#include "DDG_DenseMatrix.hpp"
#include "DDG_SparseMatrix.hpp"
#include "DDG_DiscreteExteriorCalculus.hpp"
#include "DDG_PolarDecomposition2x2.hpp"

namespace DDG
{
class JMeshElasticDeformation
{
public:
    JMeshElasticDeformation() {
        numIters = 5;
    }

    void setSource( const JMeshPtr &m);
    void setTarget( const JMeshPtr &m);

    void setNumIterations(int n) {
        numIters = n;
    }

    JMeshPtr getInterpolatedMesh( double t);

private:
    JMeshPtr  srcMesh, dstMesh, interpolatedMesh;
    DDG::Mesh source, target, mesh;
    int  numIters;

    SparseMatrix<Complex> src_L;
    SparseMatrix<Complex> tgt_L;

    void interpolate(const double t, const Mesh& source, const Mesh& target,
                     Mesh& mesh, int max_iters = 2, double tolerance = 1.0e-8);
    void init(const double step, const Mesh& source, const Mesh& target,
              Mesh& mesh);
    void computeRotation(const Mesh& meshA, const Mesh& meshB,
                         DenseMatrix<Complex>& angle) const;
    void computeLaplacian(const Mesh& mesh,
                          SparseMatrix<Complex>& L) const;
    void computeDivergence(const Mesh& mesh,
                           const DenseMatrix<Complex>& angle,
                           DenseMatrix<Complex>& rhs) const;
    void assign2DPositions(const DenseMatrix<Complex>& x, Mesh& mesh);
    void get2DPositions(const Mesh& mesh, DenseMatrix<Complex>& x) const;
};
}
