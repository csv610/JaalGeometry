///////////////////////////////////////////////////////////////////////////////
//
//  Description:  Construct Dual Graph of from a Mesh. The dual of a triangle
//  ************  is its centroid.
//  Dual Graph can be represented as Node-Node Adjacency or through explcitly
//  creating and storing the edges. When Adjacency is set to (1,0)  Edges are created
//  and stored. When Adjacency is set to (0,0)  Node-Node relationships are stored.
//
//  Chaman Singh Verma
//  Department of Computer Sciences,
//  University of Wisconsin Madison.
//  Date 15th Jan 2010.
//
//  License:  Free to download and modify.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include "Mesh.hpp"
#include "MeshTopology.hpp"
#include "DDG_MeshHot2.hpp"

namespace Jaal {

class JMeshDualGraph
{
public:
    static const int DUAL_AT_CENTER         = 0;   // Cheapest one
    static const int DUAL_AT_CIRCUMCENTER   = 1;
    static const int DUAL_AT_CIRCUMCENTER_OR_CENTER  = 2;
    static const int DUAL_AT_INCENTER       = 3;
    static const int DUAL_AT_ORTHOCENTER    = 4;
    static const int DUAL_AT_HODGE_CENTER   = 5;
    static const int DUAL_AT_RANDOM         = 6;

    static JNodePtr newObject(const JEdgePtr &c);
    static JNodePtr newObject(const JFacePtr &f);
    static JNodePtr newObject(const JCellPtr &c);

    JMeshDualGraph() {
        node_on_boundary = 0;
    }

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }
    void setDualPosition(int p ) {
        dualPos = p;
    }
    void setBoundaryNodes( bool p) {
        node_on_boundary = p;
    }

    JMeshPtr getGraph();

private:
    JMeshPtr mesh, graph;
    bool node_on_boundary;

    int  dualPos;
    void build2d();
    void build3d();
    JMeshPtr getHodgeDual();
};

//////////////////////////////////////////////////////////////////////////////

inline
JNodePtr JMeshDualGraph::newObject(const JEdgePtr &edge)
{
    JNodePtr dualnode;
    int err = edge->getAttribute("DualNode", dualnode);
    if( err ) {
        dualnode = JNode::newObject();
        edge->setAttribute("DualNode", dualnode);
        dualnode->setAttribute("PrimalEdge", dualnode);
    }

    Point3D p3d;
    edge->getAvgXYZ( p3d );
    dualnode->setXYZCoords(p3d);

    return dualnode;
}

//////////////////////////////////////////////////////////////////////////////

inline
JNodePtr JMeshDualGraph::newObject(const JFacePtr &face)
{
    JNodePtr  dualnode;
    int err = face->getAttribute("DualNode", dualnode);
    if( err ) {
        dualnode = JNode::newObject();
        face->setAttribute("DualNode", dualnode);
        dualnode->setAttribute("PrimalFace", face);
    }

    Point3D p3d;
    face->getAvgXYZ( p3d );
    dualnode->setXYZCoords(p3d);

    return dualnode;
}

//////////////////////////////////////////////////////////////////////////////

inline
JNodePtr JMeshDualGraph::newObject(const JCellPtr &cell)
{
    JNodePtr  dualnode;
    int err = cell->getAttribute("DualNode", dualnode);

    if( !err) {
        dualnode = JNode::newObject();
        cell->setAttribute("DualNode", dualnode);
        dualnode->setAttribute("PrimalCell", cell);
    }

    Point3D p3d;
    cell->getAvgXYZ( p3d );
    dualnode->setXYZCoords(p3d);

    return dualnode;
}

}
