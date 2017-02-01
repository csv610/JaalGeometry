#pragma once


#include "Mesh.hpp"
#include "NearestNeighbours.hpp"
#include "MeshRefine.hpp"
#include "MeshAffineTransforms.hpp"

class JMeshNormalClusters
{
public:
    void setMesh( const JMeshPtr &m) {mesh = m; }

    void setBinLevels(int n) { numBinLevels = n; sphMesh.reset(); }
    void createClusters();

    JMeshPtr getQuantizedSphere();  // Returns the spehre used for quantizing normals...

private:
    JMeshPtr mesh;
    JMeshPtr sphMesh;

    JNearestNeighbours jann;

    int  numBinLevels = 4;
    void buildSphere();
};
