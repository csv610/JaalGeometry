#pragma once

#include "Mesh.hpp"

#include "MeshSlicer.hpp"
#include "MeshSkeleton.hpp"

class JRingQuadMesh
{
public:
    void setMesh( const JMeshPtr &m);
    void setSkeleton( const vector<EdgeSequence> &ve);

    void setImproveMedialCenters( bool v);

    vector<JEdgeSequence> getRings();
private:
    JMeshPtr mesh;
    vector<JEdgeSequence> rings;
    vector<JEdgeSequence> skeletons;

    bool     improveCenters = 0;
    void     createRings();
}
