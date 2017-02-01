#pragma once

#include "Mesh.hpp"
#include "AllTetMeshGenerator.hpp"
#include "AllTriMeshGenerator.hpp"
#include "MeshQuality.hpp"
#include "MeshOptimization.hpp"
#include "MeshLaplacian.hpp"

#ifdef USE_IGL
#include <igl/point_mesh_squared_distance.h>
#endif

class JMeshInterpolation
{
public:
    const static int   HARMONIC_MAP = 0;
    const static int   LOCALLY_INJECTIVE_MAP = 1;
    const static int   ELASTIC_MAP = 2;

    JMeshInterpolation() {
        method = 0;
        initialized = 0;
    }

    void setSource( const JMeshPtr &s);
    void setTarget( const JMeshPtr &s);
    void setMethod( int m) {
        method = m;
    }
    JMeshPtr getInterpolatedMesh(double t);

    vector<double>  getHausdorffDistance(int dir);

private:
    int        method;
    JMeshPtr   srcMesh, dstMesh;
    JMeshPtr   simplicialMesh;
    JMeshPtr   deformedMesh;
    int        initialized;
    size_t     numSurfNodes, numSurfFaces;

    int   init();
    void  deform2D();
    void  deform3D();
    int   updatemesh();
};
