#pragma once

#ifndef MESHSEG_H
#define MESHSEG_H

#include "Mesh.hpp"
using namespace Jaal;

//Description: If by some means, on a model features curves are given, mesh
//             segment will cluster faces/cells in clusters akin to Mesh
// paritioners where reverse thing happens. Face/cells are paritioned, and
// we search for boundary curves. How and what makes edges features, depends
// on various user defined criteria.

class MeshSegmentation {
public:
    MeshSegmentation() {
        mesh = NULL;
    }

    MeshSegmentation( JMeshPtr m ) {
        mesh = m;
    }

    void setMesh( JMeshPtr m ) {
        mesh = m;
    }

    // In the graph traversal algorithm, the stopper is used to decide
    // when further prapogation must cease. For example, user may indicte
    // some edges (in 2D)  to be sharp or some geometric properties.
    // When an progation is reached to the edge which has "stopper" attribute
    // then face adjacent to it is not queued. If no stopper is used, then
    // the propagation stops at the boundaries if they present.

    // Note: this class is used in Motorcycle graph where separtices are
    // defined as stopper.

    void addStopper( const string &s ) {
        stoppers.push_back(s);
    }

    // The algorithm will assign an unique ID to a component ( which is a
    // cluster of faces(in 3D cells) enclosed by edges (or cells) which are
    // indicated by the stopper attributes.
    int getPartitions();
    int getNumPartitions();
    int getNumFeatureCurves();

private:
    JMeshPtr mesh;
    vector<string> stoppers;

    bool hasCollided(const JEdgePtr edge) {
        for( size_t i = 0; i < stoppers.size(); i++)
            if( edge->hasAttribute( stoppers[i] ) ) return 1;
        return 0;
    }
};

#endif
