#pragma once

#include "Mesh.hpp"

#include "MeshPartitioner.hpp"
#include "MeshTopology.hpp"
#include "SpectralClustering.h"

class JMeshFacesClustering
{
public:
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    void setSeeds( const JFaceSequence &f) {
        seeds = f;
    }

    JFaceSequence getSeeds() { return seeds; }

    JFaceSequence getSeeds() const { return seeds; }

    void expand();
    void removeZigZagInterfaces();

    void getSpectralClusters(int n);
private:
    JMeshPtr mesh;
    JFaceSequence seeds;
    JFacePtr   getCenter( const JFaceSequence &s);
    std::pair<int,int> max_frequency( const vector<int> &v);
    void  getPatchCenters();
};




