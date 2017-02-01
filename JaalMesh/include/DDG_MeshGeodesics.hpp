/*
 * Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow
 * Keenan Crane, Clarisse Weischedel, Max Wardetzky
 * To appear at ACM Transactions on Graphics
 *
 * TODO: pre-factorize matrices
 */

#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

#include "Mesh.hpp"

#include "DDG_Mesh.hpp"
#include "DDG_Real.hpp"
#include "DDG_Utility.hpp"
#include "DDG_DenseMatrix.hpp"
#include "DDG_SparseMatrix.hpp"
#include "DDG_DiscreteExteriorCalculus.hpp"


namespace DDG
{
class JMeshGeodesics
{
public:
    void setMesh( const JMeshPtr &m);
    void setSource( const JNodeSequence &n) {
        srcNodes = n;
    }
    void setTimeStep( double t) {
        dt = t;
    }

    int  execute();
private:
    JMeshPtr jmesh;
    DDG::Mesh mesh;
    double   dt;
    JNodeSequence srcNodes;

    SparseMatrix<Real> star0;
    SparseMatrix<Real> L;

    int  builImpulseSignal(DenseMatrix<Real>& x);
    void computeVectorField(const DenseMatrix<Real>& u);
    void computeDivergence(DenseMatrix<Real>& div);
    void assignDistance(const DenseMatrix<Real>& phi);
    void setMinToZero(DenseMatrix<Real>& phi);
};
}
#endif
