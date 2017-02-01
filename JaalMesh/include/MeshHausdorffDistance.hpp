#pragma once

#include "Mesh.hpp"
#include "MeshMatrix.hpp"
#include "AllTetMeshGenerator.hpp"
#include "AllTriMeshGenerator.hpp"

#ifdef USE_IGL
#include <igl/point_mesh_squared_distance.h>
#endif

class JMeshHausdorffDistance
{
public:
    JMeshHausdorffDistance() {}

    void setSource( const JMeshPtr &m) {
        srcMesh = m;
    }
    void setTarget( const JMeshPtr &m) {
        dstMesh = m;
    }

    vector<double>  getDistance(int dir);

private:
    JMeshPtr  srcMesh, dstMesh;
};
