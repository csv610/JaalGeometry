#pragma once

#include "Mesh.hpp"
#include "MeshMatrix.hpp"
#include "MeshAffineTransforms.hpp"
#include "MeshPartitioner.hpp"
#include "AllTriMeshGenerator.hpp"

#include <Eigen/Core>

#ifdef USE_IGL
#include <igl/boundary_loop.h>
#include <igl/lscm.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/arap.h>
#include <igl/harmonic.h>
//#include <igl/svd3x3/arap.h>
#endif

// Code borrowed from : http://www.riken.jp/brict/Yoshizawa/Research/Param.html

class JSurfaceParameterization
{
public:
    // Weighttype ...
    static const int SHAPE_PRESERVING = 0;
    static const int TUTTE            = 1;
    static const int HARMONIC_MAP     = 2;
    static const int INTRINSIC_MAP    = 3;
    static const int MEAN_VALUE       = 4;
    static const int LEAST_SQUARE_CONFORMAL   = 5;
    static const int AS_RIGID_AS_POSSIBLE  = 6;

    // BoundaryType ...
    static const int SQUARE           = 0;
    static const int CIRCLE           = 1;
    static const int NATURAL_BOUNDARY = 2;

    JSurfaceParameterization();

    int  setMesh( const JMeshPtr &m);

    void setWeight( int w) {
        weighttype = w;
    }
    void setBoundary(int b) {
        boundarytype = b;
    }
    void setStartCornerID( int id) {
        vertexCornerID = id;
    }
    void setGamma(double g) {
        gammaP = g;
    }

    void setNumIterations( int n) {
        numIters = n;
    }
 
    void setRescale( bool v) { rescaleByArea = v; }

    void setSmoothBoundary( int n) { smoothIterations = n;}
    void untangle();

    void genCharts();
    vector<JMeshPtr> getSubMeshes() const { return submeshes; }
    vector<JMeshPtr> getInvertedUVMeshes() const { return invertedUVMeshes; }

    vector<pair<double,double>> getIrregularity();

private:
    JMeshPtr  mesh;
    vector<JMeshPtr> submeshes;
    vector<JMeshPtr> invertedUVMeshes;

    int    weighttype;
    int    boundarytype;
    int    numIters;
    int    smooth;
    int    vertexCornerID;
    int    paramtype;
    double gammaP;
    double intrinsiclambda;
    bool   rescaleByArea = 0;

    int    smoothIterations = 0;

    boost::scoped_ptr< JMetisPartitioner> meshPart;

    bool isDisk();
    void normalize( const JMeshPtr &m);
    void setUVMesh( const JMeshPtr &submesh);

    int  writePly2File( const JMeshPtr &mesh, const string &file);
    JMeshPtr readPly2File( const JMeshPtr &mesh, const string &file);
    void LeastSquareConformalMap( const JMeshPtr &submesh);
    void AsRigidAsPossible( const JMeshPtr &submesh);
    JMeshPtr getTriMesh( const JMeshPtr &m);
    void rescaleUVMesh( const JMeshPtr &xyzMesh, const JMeshPtr &uvMesh, int boundarytype);
    double getIrregularity( const JMeshPtr &m);
};

