
#include "MeshOctree.hpp"

///////////////////////////////////////////////////////////////////////////
void JMeshOctree :: collectLeaf( const ONodePtr &parent)
{
    if( parent->children.empty() ) {
        leafs.push_back(parent);
        return;
    }
    for( const ONodePtr &child: parent->children)
        collectLeaf(child);
}

///////////////////////////////////////////////////////////////////////////

void JMeshOctree::split(ONodePtr &parent)
{
    if( parent->nodes.size() < size_t(maxPointsPerCell)) return;

    JNodeSequence newnodes;
    JCellSequence newcells;
    JHexRefiner::refine18( parent->hex, newnodes, newcells);
    assert( newcells.size() == 8);
    Point3D pCenter = newcells[0]->getNodeAt(6)->getXYZCoords();

    vector<ONodePtr> children(8);
    for( int i = 0; i < 8; i++) {
        children[i] = boost::shared_ptr<ONode>( new ONode);
        children[i]->hex = newcells[i];
    }

    std::bitset<3> cellbits;
    for( const JNodePtr &vtx : parent->nodes) {
        cellbits.reset();
        const Point3D &xyz = vtx->getXYZCoords();
        if( xyz[0] >= pCenter[0] ) cellbits.set(0);
        if( xyz[1] >= pCenter[1] ) cellbits.set(1);
        if( xyz[2] >= pCenter[2] ) cellbits.set(2);
        int cellid = cellbits.to_ulong();
        assert( cellid >= 0 && cellid < 8);
        children[cellid]->nodes.push_back(vtx);
    }
    parent->children = children;
    parent->nodes.clear();
}

///////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshOctree :: getVoxels()
{
    JMeshPtr  voxels;
    if( mesh == nullptr) return nullptr;

    JNodeSequence nodes = mesh->getNodes();

    ONodePtr parent = boost::shared_ptr<ONode>(new ONode);
    parent->nodes = nodes;

    JBoundingBox box = JNodeGeometry::getBoundingBox(nodes);
//    box.expandBy(1.0001);

    JNodeSequence hexnodes(8);
    for( int i = 0; i < 8; i++) {
        hexnodes[i] = JNode::newObject();
        const Point3D &p = box.getCorner(i);
        hexnodes[i]->setXYZCoords(p);
    }
    nodeQ.clear();
    parent->hex = JHexahedron::newObject(hexnodes);
    nodeQ.push_back(parent);

    while(!nodeQ.empty()) {
        ONodePtr currNode = nodeQ.front();
        nodeQ.pop_front();
        split(currNode);
        for( const ONodePtr &child : currNode->children)
            nodeQ.push_back(child);
    }

    leafs.clear();
    collectLeaf( parent );

    JCellSequence hexcells;
    for( const ONodePtr &n : leafs)
        hexcells.push_back(n->hex);

    JMeshTopology::getEntitySet(hexcells, hexnodes);
    voxels = JMesh::newObject();
    voxels->addObjects(hexnodes);
    voxels->addObjects(hexcells);
    cout << "Stage 3 " << hexnodes.size() << "  " << hexcells.size() << endl;
    return voxels;
}
