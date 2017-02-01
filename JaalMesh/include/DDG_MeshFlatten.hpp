/*
 * Spectral Conformal Parameterization
 * Patrick Mullen, Yiying Tong, Pierre Alliez, Mathieu Desbrun.
 * Symposium of Geometry Processing, 2008.
 */


#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Mesh.hpp"

#include "DDG_Mesh.hpp"
#include "DDG_Complex.hpp"
#include "DDG_DenseMatrix.hpp"
#include "DDG_SparseMatrix.hpp"
#include "DDG_DiscreteExteriorCalculus.hpp"

namespace DDG
{
class JMeshFlatten
{
public:

    void setMesh( const JMeshPtr &m);
    void execute();

protected:
    JMeshPtr jmesh;
    DDG::Mesh mesh;
    void buildEnergy(SparseMatrix<Complex>& A) const;
    void assignSolution(const DenseMatrix<Complex>& x);
    void normalizeMesh(const double scale);
};
}

#endif
