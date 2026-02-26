#pragma once

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

#include<list>

#include "Mesh.hpp"
#include "polypartition.hpp"

using namespace std;
using namespace Jaal;

class JPolyPartitioner {
public:
    static const int EAR_CLIPPED_POLYGONS     = 0;
    static const int EDGE_LENGTH_OPT_POLYGONS = 1;
    static const int MONOTONE_POLYGONS        = 2;
    static const int CONVEX_HM_POLYGONS       = 3;
    static const int CONVEX_OPT_POLYGONS      = 4;

    JPolyPartitioner() {}

    void setMesh( const JMeshPtr &m) { mesh = m; }

    JMeshPtr getPartitions(int algo = CONVEX_OPT_POLYGONS);

protected:
    JMeshPtr  mesh, partmesh;
    JNodeSequence boundnodes;
    const JNodePtr &searchNode( const Point2D &xy);
    int createPoly( const JEdgeSequence &edges, TPPLPoly &poly, bool hole);
    int createPolyList( list<TPPLPoly> &polys);
    int createFace( TPPLPoly &p);
};
