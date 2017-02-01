#pragma once

#include "Mesh.hpp"
#include "MSTQuadMesher.hpp"
#include "MeshRefine.hpp"
#include "Singlet.hpp"

class JMeshMinSingularity
{
public:
    void setMesh( const JMeshPtr &m);
    void setEdgeLength( double l ) {
        refEdgeLength = l;
    }

    JMeshPtr  refineAll();
private:
    JMeshPtr mesh;
    JMeshPtr newmesh;
    double   refEdgeLength;
    JMSTQuadMesher polymesher;

    void refine(const JFacePtr &f);
    void refine(const JEdgePtr &e);
    bool isBoundaryEven(const JFacePtr &f);
    void remeshBoundNode(const JNodePtr &v);
};

