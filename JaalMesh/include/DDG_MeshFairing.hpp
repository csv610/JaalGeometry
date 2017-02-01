/*
 * Implicit Fairing of Arbitrary Meshes using Diffusion and Curvature Flow
 * Mathieu Desbrun, Mark Meyer, Peter Schröder, Alan H. Barr
 * ACM Siggraph '99 Proceedings
 */

#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Mesh.hpp"

#include "DDG_Mesh.hpp"
#include "DDG_Real.hpp"
#include "DDG_DenseMatrix.hpp"
#include "DDG_SparseMatrix.hpp"
#include "DDG_DiscreteExteriorCalculus.hpp"

namespace DDG
{
class JMeshFairing
{
public:
    void setMesh( const JMeshPtr &m);

    void setTimeStep( double dt) {
        timeStep = dt;
    }
    void nextStep();

protected:
    JMeshPtr  jmesh;
    double    timeStep;

    DDG::Mesh mesh;
    /*
        DDG::SparseMatrix<Real> star0;
        DDG::SparseMatrix<Real> star1;
        DDG::SparseMatrix<Real> d0;
    */
    DDG::DenseMatrix<Real> x;

    void getPositions();
    void setPositions();
    void updateMesh();
};
}

#endif
