///////////////////////////////////////////////////////////////////////////////
// Copyright 2013 - Christian Schüller 2013, schuellc@inf.ethz.ch
// Interactive Geometry Lab - ETH Zurich
//
// Heavily modified: Chaman Singh Verma
// Univ of Wisconsin, Madison, Nov 2013
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "Mesh.hpp"
#include "MeshTopology.hpp"
#include "AllTriMeshGenerator.hpp"

#include <LIMSolverInterface.h>

using namespace Jaal;
using namespace Eigen;

class LIMSolver;
class DeformableMesh;
class TriangleMesh;
class TetrahedronMesh;

class JLocallyInjectiveMap
{
//  typedef int IndexType;
public:
    static const int DIRICHLET_ENERGY              = 0;
    static const int LAPLACIAN_ENERGY              = 1;
    static const int GREEN_STRAIN_ENERGY           = 2;
    static const int AS_RIGID_AS_POSSIBLE_ENERGY   = 3;
    static const int LEAST_SQUARE_CONFORMAL_ENERGY = 4;
    static const int POISSON_ENERGY                = 5;
    static const int IDENTITY_ENERGY               = 6;
    static const int UNIFORM_LAPLACE_ENERGY        = 7;

    static int  getEnergyType( const string &s);

    static const int ALPHA_UPDATE         = 0;
    static const int BARRIER              = 1;
    static const int BARRIER_COMPENSATION = 2;
    static const int LOG_BARRIER          = 3;
    static const int NEOHOOKEAN_BARRIER   = 4;
    static const int SUB_STEPPING         = 5;

    JLocallyInjectiveMap();

    ~JLocallyInjectiveMap();

    void setMesh( const JMeshPtr &m );

    void clear();

    // At the rest position (t=0) specify nodes that are allowed to move.  A user specifies
    // the desired locations of these source nodes. The solver will do it best to relocate
    // the source vertices to their target positions. But in many cases, it may not be possible
    // to meet the desired locatiom, as it may cause inversion of cetain elements. Either
    // we may restrict vertex movement till an optimal position is reached i.e the maxumum
    //  distant osition beyon whoch there will be inverted elements or use may be allowed to
    // proceed with creation of inverted eleme nts.

    void setEnergyType( int e );

    void setBarriers(bool v) {
        enableBarriers = v;
    }

    void setBarrierType(int t, bool v) {
        if( t == LOG_BARRIER) enableLogBarriers = v;
        if( t == NEOHOOKEAN_BARRIER) enableNeoHookeanBarriers = v;
    }

    void setBarrierCompensation(bool v, double val = 1) {
        enableBarrierCompensation = v;
    }

    void setAlphaUpdate( bool v) {
        enableAlphaUpdate = v;
    }

    void setSubstepping( bool v) {
        enableSubstepping = v;
    }

    /*
        void setAlpha( double a) {
            if( solver ) solver->Alpha = a;
            alpha = a;
        }

        double getAlpha() const {
            return alpha;
        }

        void setAlphaRatio( double a) {
            if( solver ) solver->AlphaRatio = a;
            alphaRatio = a;
        }

        double getAlphaRatio() const {
            return alphaRatio;
        }
    */

    void setBeta( double b) {
//      if( solver ) solver->Beta = b;
        beta  = b;
    }

    double getBeta() const {
        return beta;
    }

    void setGamma( double g) {
//      if( solver ) solver->Gamma = g;
        gamma = g;
    }

    double getGamma() const {
        return gamma;
    }

    void setMaxIterations( int v ) {
        maxIterations = v;
    }

    int getMaxIterations() const {
        return maxIterations;
    }

    void setSmallestArea( double a = -1) {
        smallestArea = a;
    }

    // If you wish to go for incremental solver, so that we can visualize intermediate steps,
    // initialize the data first ...
    int initSolve();

    //  Call as many time as you want after initializing the sttructure....
    int stepSolve();

    /*
        // We can set a trajectory of the displacements ...
        void setBoundaryPreservation( bool v) { enableBoundaryPreservation = v; }
        void translate( JNodeSequence &nodeSeq, const Array3D &vec);
        void rotate( JNodeSequence &nodeSeq, double angle, const Array3D *center = nullptr);
        void setDisplacements( const JNodeSequence &seq, bool tofro_ = 0);
    */
    int  solve();

    void saveAnimation( bool v) {
        saveAnim = v;
    }
    void setOriginalMesh() const;
    void setDeformedMesh() const;

    double getMaxDistance() const;

    void setConstraintsPosition(const vector<int>& constraintVertices,
                                const Eigen::Matrix<double,Dynamic,3>& positions);
    void setConstraints(const vector<int>& constraintVertices);

protected:
    JMeshPtr inMesh;   // Input mesh...
    JMeshPtr simplicialMesh;   // Input mesh...

    LIMData *limData;

    vector<int> borderVertices;
    Eigen::Matrix<double,Eigen::Dynamic,3> deformedVertices;
    Eigen::Matrix<double,Eigen::Dynamic,3> initialVertices;
    Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> elements;
    Eigen::SparseMatrix<double> constraintMatrix;
    Eigen::Matrix<double,Eigen::Dynamic,1> constraintTargets;
    Eigen::Matrix<double,Eigen::Dynamic,1> gradients;
    int  dim;
    int  elemType;
    int  energyType;

    bool enableBarriers;
    bool enableSubstepping;
    bool enableAlphaUpdate;
    bool enableOutput;
    bool enableLogBarriers;
    bool enableNeoHookeanBarriers;
    bool enableBarrierCompensation;
    bool findLocalMinima;
    bool enableBoundaryPreservation;
    bool derivedMesh = 0;

    double positionalConstraintError;
    double barrierWeight;

    double smallestArea;
    double alpha, alphaRatio, beta, gamma;
    double error, tolerance;

    bool   saveAnim;
    int    maxIterations, iterCounter;

    void   createSimplicialMesh();

    void initParams();
    void initMesh();
    void setTargets();

    bool hasConverged() const;
    void updateMesh();
    int  getEnergy() const;
    void displaceConstraints( const Point3D &tail, const Point3D &head);

    /*
    //  bool showInvertedElements;
    //  LIMSolver *solver;
    //  DeformableMesh *deformMesh;
        void initEnergy();
        void initSolver2D();
        void initSolver3D();
        void restart();
    */
};

