#pragma once

#ifndef QUADPARAMLINES_H
#define QUADPARAMLINES_H

#include "Mesh.hpp"
#include "MeshTopology.hpp"

namespace Jaal {
class  QuadParametricLines {
public:
    QuadParametricLines() {
    }

    int setLines( JMeshPtr m );

private:
    void basicOp( JEdgePtr edge);

    deque<JEdgePtr> edgeQ;
    JMeshPtr mesh;
    JNodeSequence sequence;
};

}

#endif

