#pragma once

#include "Mesh.hpp"
#include "BoundingBox.hpp"
#include "MeshRefine.hpp"
#include <bitset>

class JMeshOctree
{
public:
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    void setNumOfPointsPerCell(int n) {
        maxPointsPerCell = n;
    }
    JMeshPtr  getVoxels();
private:
    class ONode;
    typedef boost::shared_ptr<ONode> ONodePtr;

    struct ONode
    {
        JCellPtr         hex;
        JNodeSequence     nodes;
        vector<ONodePtr> children;
    };

    void      collectLeaf( const ONodePtr &parent);
    void      split(ONodePtr &p);

    JMeshPtr mesh;
    int maxPointsPerCell;
    std::deque<ONodePtr>  nodeQ;
    std::list<ONodePtr>   leafs;
};
