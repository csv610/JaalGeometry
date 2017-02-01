#pragma once

#ifndef RANGE_SEARCH_H
#define RANGE_SEARCH_H

#include "Mesh.hpp"

using namespace Jaal;

class JRangeSearch
{
public:
    JRangeSearch() {
        mesh = nullptr;
    }

    void setMesh(JMeshPtr m) {
        mesh = m;
    }

    void isSimplex() const;

    void setRay( double *p0, double *p1);
    void setQueryRegion( double *p0, double *p1, double *p2);
    void setQueryRegion( double *p0, double *p1, double *p2, double *p3);

    void setQueryRegion( JFacePtr f);

    JFaceSequence getOverlappedFaces() {
        return overlapFaces;
    }
    JEdgeSequence getOverlappedEdges();
    JNodeSequence getOverlappedNodes();

private:
    JMeshPtr mesh;
    JNodeSequence overlapNodes;
    JEdgeSequence overlapEdges;
    JFaceSequence overlapFaces;

    void bfs_search( JFacePtr f);
    void bfs_search( JNodePtr v);
};

#endif
