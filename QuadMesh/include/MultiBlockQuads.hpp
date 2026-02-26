#pragma once

#ifndef MULTIBLOCK_QUAD_H
#define MULTIBLOCK_QUAD_H

#include "Mesh.hpp"

class JMultiBlockQuadMesh
{
public:
    JMultiBlockQuadMesh() {
        baseQuadmesh = nullptr;
        newQuadmesh  = nullptr;
        minDivisdions = 2;
    }

    void setBaseQuadMesh( Mesh *q)   {
        baseQuadmesh = q;
    }
    void setMinEdgeDivisions( int n) {
        minDivisions = n;
    }

    Mesh *getRefinedMesh() const;
private:
    Mesh *baseQuadmesh;
    Mesh *newQuadmesh;
    void refineEdge( Edge *e);
    void refineQuad( Face *f);
};
#endif
