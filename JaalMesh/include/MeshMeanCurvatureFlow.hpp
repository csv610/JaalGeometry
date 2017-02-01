#pragma once

#include <iostream>

#ifdef USE_IGL
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/repdiag.h>
#include <igl/matlab/matlabinterface.h>
#endif

#include "Mesh.hpp"
#include "MeshAffineTransforms.hpp"
#include "basic_math.hpp"
#include "AllTriMeshGenerator.hpp"
#include "DelaunayMesh.hpp"
#include "MeshMatrix.hpp"

#include "DDG_MeshFairing.hpp"

class JMeshMeanCurvatureFlow
{
public:
    static const int IGL_METHOD    = 0;
    static const int KEENAN_METHOD = 1;
    static const int CG_SOLVER     = 0;
    static const int HSC_SOLVER    = 1;

    JMeshMeanCurvatureFlow();
    ~JMeshMeanCurvatureFlow();

    void setMesh( const JMeshPtr &m);

    // We always work with derived surface mesh of the mesh...
    JMeshPtr getSurfaceMesh() const { return triMesh; }

    // There are two method which are almost similar (1) Cotan   (2) Keenan
    void setAlgorithm( int a);

    // For fast convergence, perform "intrinsic Delaunay" tringulation on the input
    // triangle mesh...
    void setIntrinsicDelaunayMesh( bool v)
    {
        intrinsicDelaunayMesh = v;
    }

    // Bring the model to its original position ...
    void restart();

    // Calculate the position of the model after "n" iterations...
    void   nextStep(int n = 1);

    // Enable/disable time steps ...
    void   setAdaptiveTimeSteps( bool v ) {
        adaptiveTimeSteps = v;
    }

    //  Provide the time step manually ....
    void   setTimeStep( double t) {
        timeStep = t;
    }

    double getMaximumDisplacement() const {
        return maxDistance;
    }

    // If the model has euler characterisitic of 2, and the model is almost spherical,
    // directly project the vertices on the sphere...
    void   projectOnSphere();

    //  The shape of the model is continuously monitored on the witness nodes. If no
    //  Witness nodes are provide then we choose 1% nodes randomly from the input nodes...
    //  The input can also come from the "Mesh Sampling" procedures...
    void   setWitnessNode(JNodeSequence &w) {
        witnessNodes = w;
    }

    // Provides displacement of each vertex from its original position ....
    vector<double> getDisplacements();


    void setImproveInternalQuality( bool q ) { improveInternalQuality = q; }

private:
    JMeshPtr  orgMesh, triMesh;
    JMeshAffineTransform affine;
    JNodeSequence triBoundNodes;
    JNodeSequence orgBoundNodes;
    JNodeSequence witnessNodes;

    int  algorithm = KEENAN_METHOD;
    int  linearSolver = CG_SOLVER;

    // For back-up coordinates ...
    vector<size_t>  l2g;
    vector<double>  orgCoords;

    // For Cotangent method from IGL ....
    Eigen::MatrixXd V,U;
    Eigen::MatrixXi F;
    Eigen::VectorXd dblA;
    Eigen::SparseMatrix<double> L;
    Eigen::SparseMatrix<double> M;
    Eigen::MatrixXd BC;
    Eigen::RowVector3d centroid;

    // Keenan Crane's method...
    DDG::JMeshFairing meshFair;

#ifdef USE_IGL
    Engine  *engine;   // Matlab ....
#endif

    bool   adaptiveTimeSteps = 0;
    bool   intrinsicDelaunayMesh = 0;
    bool   improveInternalQuality = 0;
    int    currStep;
    int    numNodesMoved;
    double scaleFactor = 1.0;
    double timeStep;
    double maxDistance;
    double surfArea[2];
    double sphRadius;

    void  getSurfTriMesh();
    void  updateMesh( const Eigen::MatrixXd &nodeCoords);
    void  getLMatrix();
    void  isSphere();
    void  initMesh();     // Initialization of IGL method ...

    void  solveSystem( const Eigen::SparseMatrix<double> &A,
                       const Eigen::MatrixXd &b, Eigen::MatrixXd &x);
    void  getTrucatedLaplaceMatrix();
    int   stdPCG( const string &s);
};

