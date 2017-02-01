/*
 * HOT: Hodge-Optimized Triangulations
 * Patrick Mullen, Pooran Memari, Fernando de Goes, Mathieu Desbrun.
 * Transactions on Graphics (SIGGRAPH) 2011.
 */

#pragma once

#include "Mesh.hpp"

#include "DDG_Mesh.hpp"
#include "DDG_DenseMatrix.hpp"
#include "DDG_SparseMatrix.hpp"
#include "DDG_DiscreteExteriorCalculus.hpp"

namespace DDG
{
class JMeshHot2
{
public:
    void setMesh( const JMeshPtr &m);
    JMeshPtr getDualMesh();

protected:
    JMeshPtr   jmesh;
    DDG::Mesh  mesh;

    void optimizeWeights();
    void buildRhs(DenseMatrix<Real>& rhs);
    void assignSolution(const DenseMatrix<Real>& x);
    JNodePtr getNewNode( const Vector &p);
};
}
