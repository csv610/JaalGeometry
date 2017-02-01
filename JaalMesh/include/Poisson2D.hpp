#pragma once

#ifndef POISSON2D_H
#define POISSON2D_H

#define EIGNN_SUPERLU_SUPPORT


#include "FEM2D.hpp"

#include <map>
#include "FEMSpace.hpp"
#include "MeshTangle.hpp"

#include "MeshExporter.hpp"
#include "SparseMatrix.hpp"

///////////////////////////////////////////////////////////////////////////////
//
// Description:  For the 2D triangle mesh, this class calculates the displacement,
//               and stresses.
// What input is needed:
//	(A)  A triangle mesh and mesh need not be simplicial in the sense that
//           triangle elements can overlap. But there are some conditiosns:
//	     (1) The mesh shhould be sliver free.
//           (2) No element should tangle with boundary edges.
//      (B)  Boundary conditions:
//           Force boundary conditions: Force applied to the geometric edges.
//           Displacemennt specfied to some to the boundary nodes. Most of the
//           the values are zero.
//
// License:  It is free and has absolutely no restrictions whatsoever.
///////////////////////////////////////////////////////////////////////////////
//
// The initial code was written in Matlab, but was converted to C++ because of
// my inertia to learn new language and (II) I want to experiment with new C++
// features.
//          Chaman Singh Verma
//          Department of computer Sciences.
//          University of Wisconsin,
//          Madison
// The original Matlab code was provided by:
//          Prof. Suresh Krishnan,
//          Department of Mechanical Engineering.
//          University of Wisconsin, Madison
//          Madison.
//
//////////////////////////////////////////////////////////////////////////////

class JPoisson2D : public JFEM2D {
public:
    JPoisson2D() ;

    void setDirichletValue( const JNodePtr &v, double u);
    void setDirichletValue( const JEdgePtr &e, double u);
    void setDirichletValue( const JEdgeSequence &ve, double u);

    void setNeumannValue( const JNodePtr &v, double u);
    void setNeumannValue( const JEdgePtr &e, double u);
    void setNeumannValue( const JEdgeSequence &ve, double u);

    void getField( vector<double> &u);

    int getValueAt( int id, double &val)
    {
        val = u[id];
        return 0;
    }

    int solve();


    void saveAs( const string &filename, const string &varname);

    void  integrate_self_tangle(const JFacePtr &face, MatrixXd &Ke, vector<double> &fe);

private:
//  JMeshTangle meshtangle;
//  vector< vector<double> >  NEdge, gradU; //  Basis function and gradient at Gauss points for an edge.

    vector<double> u;             // Displacement vector in Ku = f
    vector<double> f, freduced;  // force vector in Ku = f;
    vector<vector<double> > NFace; // Basis function for a face.
    vector<MatrixXd> gradUV;     // Gradident of Basis function for each Gauss point..

    void initParams();

    void BMatrix( const JFacePtr &face, const Point2D &uv, MatrixXd &B);
    void BMatrix( const JFacePtr &face, const Point2D &uv, MatrixXd &B, double &dJ);
    void BMatrix(const double *uv, const Matrix2d &invJ, MatrixXd &B);
    void BMatrix(const MatrixXd &gradUV, const Matrix2d &invJ, MatrixXd &Bg);

    int   buildLinearSystem();
    void  assembleK();
    void  applyNeumannBC();
    void  applyDirichletBC();
    void  applyBC();

    // Integration of a triangles within the overlap region ...
    int   integrate(const JFacePtr &face, const JFacePtr &f, MatrixXd &K);

    void  integrate( const JFacePtr &face1, const JFacePtr &face2,
                     const vector<Point2D> &triPnts, MatrixXd &Ke);

    //   Integration over the element for K = sum( B^T*D*B dx)
    void  integrate(const JFacePtr &face, MatrixXd &K, vector<double> &f);

    //  Integrate over the boundry edge for the right hand side
    //  vector f = sum ( N f dx )
    void  integrate(const JEdgePtr &e, double val, vector<double> &f);

    int  solveLinearSystem();

    void assembleTangledK( const JFacePtr &f0, const JFacePtr &f1);
    void assembleTangledK();
};

//////////////////////////////////////////////////////////////////////////////

#endif
