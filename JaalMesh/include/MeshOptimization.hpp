#pragma once

#include "Mesh.hpp"
#include "LloydOptimizer.hpp"
#include "MeshGeometry.hpp"
#include "MeshTopology.hpp"
#include "Mesquite_all_headers.hpp"
#include "LocallyInjectiveMap.hpp"
#include "AllTetMeshGenerator.hpp"
#include "MeshQuality.hpp"

using namespace std;

struct JMeshOptimization {
    static JLogger *logger;
    JMeshOptimization()
    {
        mesh = nullptr;
        geomDim = 0;
        meshDim = 0;
    }

    void setMesh( const JMeshPtr &m);

    JMeshPtr mesh;
    int   geomDim;
    int   meshDim;
};

struct JMeshTopologyOptimization : public JMeshOptimization {
    int bandwidth_reduction(const JMeshPtr  &m);
};

///////////////////////////////////////////////////////////////////////////////////

class JMeshGeometricOptimization: public JMeshOptimization {

public:
    JMeshGeometricOptimization() {}

    void setBoundaryPreserve(bool v)
    {
        preserve_boundary = v;
    }

    void useSimplicial(bool v)
    {
        use_simplicial = v;
    }

    void setNumIterations(int n )
    {
        numIter = n;
    }
    int getNumIterations() const {
        return numIter;
    }

    void addConstraints( const JNodeSequence &v, int groupID = 0);
    void removeAllConstraints();


protected:
    int  numIter = 1;
    bool preserve_boundary = 1;
    bool use_simplicial    = 1;
    vector<Point3D> newCoords;

    void addBoundaryConstraints( const JMeshPtr &m);
};

///////////////////////////////////////////////////////////////////////////////////
class JMeshNonlinearOptimization : public JMeshGeometricOptimization {
public:
    static const int QUASI_NEWTON       = 0;
    static const int STEEPEST_DESCENT   = 1;
    static const int TRUST_REGION       = 2;
    static const int FEASIBLE_NEWTON    = 3;
    static const int CONJUGATE_GRADIENT = 4;
    static const int SMART_LAPLACIAN    = 5;

    static const int INVERSE_MEAN_RATIO = 0;
    static const int MEAN_RATIO         = 1;
    static const int CONDITION_NUMBER   = 2;
    static const int EDGE_LENGTH        = 3;

    static const int LINF_NORM          = 0;
    static const int LMEAN_NORM         = 1;
    static const int LP_NORM            = 2;


    JMeshNonlinearOptimization()
    { }

    void setQualityMetric( int q )
    {
        qualMetric = 0;
        if( q  < 4) qualMetric = q;
    }

    void setAlgorithm( int a )
    {
        algorithm = 0;
        if( a  < 6) algorithm = a;
    }
    void setNorm( int n ) {
        norm = n;
    }
    void setNormVal( int n ) {
        normVal = n;
    }

    int improveQuality(bool convert_2_simplex = 1);
    int untangle();

    void setAveragingMethod( Mesquite2::QualityMetric::AveragingMethod a) {
        average_method = a;
    }
    void setPatchType ( int a ) {
        patch_type = a;
    }

    int  shapeOptimize();

private:
    int algorithm  =  QUASI_NEWTON;
    int qualMetric =  INVERSE_MEAN_RATIO;
    int normVal    =  2;
    int norm       =  LP_NORM;

    Mesquite2::QualityMetric::AveragingMethod  average_method = Mesquite2::QualityMetric::AveragingMethod::SUM;
    int patch_type     = 0;  // local method ...

    vector<int>     vfixed;
    vector<size_t>  vNodes;
    vector<int>     etopo;
    vector<size_t>  l2g;
    vector<double>  vCoords;

    int execute(const JMeshPtr &m);
    int qualityOptimize2D();
    int qualityOptimize3D();
    int untangle( Mesquite2::ArrayMesh *m);

    Mesquite2::ArrayMesh* jaal_to_mesquite( const JMeshPtr &m);
};

//////////////////////////////////////////////////////////////////////////////////////

class JCyclicQuadMeshOptimizer : public JMeshGeometricOptimization
{
public:
    int  smoothAll();
    int  smooth( const JNodeSequence &v);
private:
    int atomicOp( const JNodePtr &vtx);
};

//////////////////////////////////////////////////////////////////////////////////////

struct JBoundedDistortionOptimizer : public JMeshGeometricOptimization
{
    int optimize( const JMeshPtr &m, double K);
};
