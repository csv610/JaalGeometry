#pragma once

#include "Mesh.hpp"
#include "LocallyInjectiveMap.hpp"
#include "MeshLaplacian.hpp"

class JMeshOptBoundaryLayer
{
public:
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }
    void setOffset( double d)        {
        offset = d;
    }

    void setNormals();
    void update();

private:
    JMeshPtr mesh;
    JNodeSequence boundNodes, layerNodes;
    double       offset = 1.0E-06;
    bool         normals_available = 0;
    JLocallyInjectiveMap lim;

};
