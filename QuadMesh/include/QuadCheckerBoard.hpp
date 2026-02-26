#pragma once

#ifndef QUADTEXTURE
#define QUADTEXTURE

#include "Mesh.hpp"

class JQuadCheckerBoard
{
public:
    JQuadCheckerBoard() {}
    void setMesh(const JMeshPtr &m) {
        mesh = m;
    }

    void genPattern();
private:
    JMeshPtr mesh;
};

#endif
