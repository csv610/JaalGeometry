#include "MeshRefine.hpp"
using namespace Jaal;

int JTriRefiner :: refineNode3(const JNodePtr &vertex, const JEdgePtr &edge1,
                               const JEdgePtr &edge2)
{
    JNodeSequence nodes(13);

    nodes[0] = vertex;
    nodes[1] = edge1->getNodeAt(0);
    nodes[2] = edge1->getNodeAt(1);

    if( edge1->getNodeAt(0) == edge2->getNodeAt(0) )
        nodes[3] = edge2->getNodeAt(1);

    if( edge1->getNodeAt(0) == edge2->getNodeAt(1) )
        nodes[3] = edge2->getNodeAt(0);

    JFacePtr f1 = JSimplex::getFaceOf(vertex, nodes[1], nodes[2] );
    JEdge::getRelations(edge1, faceneighs);
    JFacePtr of1;
    if( faceneighs[0] == f1 ) of1 = faceneighs[1];
    if( faceneighs[1] == f1 ) of1 = faceneighs[0];
    nodes[4] = JTriangle::getOppositeNode(of1, nodes[1], nodes[2] );
    if( nodes[4]->getNumRelations(0) >= 7) return 1;

    JFacePtr f2 = JSimplex::getFaceOf(vertex,nodes[1], nodes[3] );
    JEdge::getRelations(edge2, faceneighs);
    JFacePtr of2;
    if( faceneighs[0] == f2 ) of2 = faceneighs[1];
    if( faceneighs[1] == f2 ) of2 = faceneighs[0];
    nodes[5] = JTriangle::getOppositeNode(of2, nodes[1], nodes[3]);
    if( nodes[5]->getNumRelations(0) >= 7) return 1;

    nodes[6]  = JNodeGeometry::getMidNode( nodes[1], nodes[2] );
    nodes[7]  = JNodeGeometry::getMidNode( nodes[1], nodes[3] );
    nodes[8]  = JNodeGeometry::getMidNode( nodes[0], nodes[1], 0.33 );
    nodes[9]  = JNodeGeometry::getMidNode( nodes[0], nodes[1], 0.66 );
    nodes[10] = JFaceGeometry::getCentroid( vertex, nodes[1], nodes[2]);
    nodes[11] = JFaceGeometry::getCentroid( vertex, nodes[1], nodes[3]);
    nodes[12] = JFaceGeometry::getCentroid( vertex, nodes[1], nodes[6]);

    JFaceSequence newFaces(13);
    newFaces[0]  = JTriangle::newObject( nodes[1], nodes[4], nodes[6] );
    newFaces[1]  = JTriangle::newObject( nodes[2], nodes[6], nodes[4] );
    newFaces[3]  = JTriangle::newObject( nodes[1], nodes[7], nodes[5] );
    newFaces[4]  = JTriangle::newObject( nodes[3], nodes[5], nodes[7] );
    newFaces[5]  = JTriangle::newObject( nodes[0], nodes[10], nodes[2] );
    newFaces[6]  = JTriangle::newObject( nodes[2], nodes[10], nodes[6] );
    newFaces[7]  = JTriangle::newObject( nodes[0], nodes[3], nodes[11] );
    newFaces[8]  = JTriangle::newObject( nodes[3], nodes[7], nodes[11] );
    newFaces[9]  = JTriangle::newObject( nodes[0], nodes[11], nodes[8] );
    newFaces[10] = JTriangle::newObject( nodes[8], nodes[11], nodes[9] );
    newFaces[11] = JTriangle::newObject( nodes[7], nodes[9], nodes[11] );
    newFaces[12] = JTriangle::newObject( nodes[1], nodes[9], nodes[7] );
    newFaces[13] = JTriangle::newObject( nodes[0], nodes[8], nodes[10] );
    newFaces[14] = JTriangle::newObject( nodes[8], nodes[12], nodes[10] );
    newFaces[15] = JTriangle::newObject( nodes[8], nodes[9], nodes[12] );
    newFaces[16] = JTriangle::newObject( nodes[1], nodes[12], nodes[9] );
    newFaces[17] = JTriangle::newObject( nodes[1], nodes[6], nodes[12] );

    if( outmesh )  {
        outmesh->addObject(nodes[6] );
        outmesh->addObject(nodes[7] );
        outmesh->addObject(nodes[8] );
        outmesh->addObject(nodes[9] );
        outmesh->addObject(nodes[10] );
        outmesh->addObject(nodes[11] );
        outmesh->addObject(nodes[12] );
        outmesh->addObjects( newFaces);
    }

    f1->setStatus( JMeshEntity::REMOVE );
    f2->setStatus( JMeshEntity::REMOVE );
    of1->setStatus( JMeshEntity::REMOVE );
    of2->setStatus( JMeshEntity::REMOVE );
    edge1->setStatus( JMeshEntity::REMOVE);
    edge2->setStatus( JMeshEntity::REMOVE);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////

int JTriRefiner :: refineNode3(const JNodePtr &vertex)
{
    if( !vertex->isActive() )  return 1;
    JNode::getRelations(vertex, nodeneighs);
    if( nodeneighs.size() != 3 ) return 2;

    for( int i = 0; i < 3; i++)
        if( nodeneighs[i]->getNumRelations(0) == 7 ) return 3;

    for( int i = 0; i < 3; i++) {
        JEdgePtr edge1 = JSimplex::getEdgeOf( nodeneighs[(i+0)%3], nodeneighs[(i+1)%3] );
        JEdgePtr edge2 = JSimplex::getEdgeOf( nodeneighs[(i+1)%3], nodeneighs[(i+2)%3] );
        int err = refineNode3(vertex, edge1, edge2);
        if( !err ) return 0;
    }
    return 1;

}

////////////////////////////////////////////////////////////////////////////////////

int JTriRefiner :: refineNode4(const JNodePtr &vertex, const JEdgePtr &edge)
{
    JNodeSequence nodes(12);

    nodes[0] = vertex;
    nodes[1] = edge->getNodeAt(0);
    nodes[2] = edge->getNodeAt(1);

    JNode::getRelations(vertex, faceneighs);
    if( faceneighs.size() != 4 ) return 1;

    JEdgeSequence edges(4);
    for( int i = 0; i < 4; i++)
        edges[i] = JTriangle::getOppositeEdge( faceneighs[i], vertex);
    JEdgeTopology::getChain(edges, edge);
    nodes[3] = edges[2]->getNodeAt(0);
    nodes[4] = edges[3]->getNodeAt(0);

    JFacePtr tri1 = JSimplex::getFaceOf(vertex, nodes[0], nodes[1] );
    JFaceSequence edgefaces;
    JEdge::getRelations(edge, edgefaces);
    if( edgefaces.size() != 2 ) return 3;

    JFacePtr tri2;
    if( edgefaces[0] == tri1 ) tri2 = edgefaces[1];
    if( edgefaces[1] == tri1 ) tri2 = edgefaces[0];

    nodes[5] = JTriangle::getOppositeNode( tri2, nodes[1], nodes[2] );
    if( nodes[5]->getNumRelations(0) >= 7) return 4;

    nodes[6]  = JNodeGeometry::getMidNode( nodes[1], nodes[2] );
    nodes[7]  = JNodeGeometry::getMidNode( nodes[0], nodes[2] );
    nodes[8]  = JNodeGeometry::getMidNode( nodes[0], nodes[1] );
    nodes[9]  = JNodeGeometry::getMidNode( nodes[6], nodes[7] );
    nodes[10] = JNodeGeometry::getMidNode( nodes[7], nodes[8] );
    nodes[11] = JNodeGeometry::getMidNode( nodes[6], nodes[8] );

    JFaceSequence newFaces(16);
    newFaces[0] = JTriangle::newObject( nodes[3], nodes[0], nodes[7] );
    newFaces[1] = JTriangle::newObject( nodes[3], nodes[7], nodes[2]);

    newFaces[2] = JTriangle::newObject( nodes[4], nodes[1], nodes[8]);
    newFaces[3] = JTriangle::newObject( nodes[4], nodes[8], nodes[0]);

    newFaces[4] = JTriangle::newObject( nodes[5], nodes[2], nodes[6]);
    newFaces[5] = JTriangle::newObject( nodes[5], nodes[6], nodes[1]);

    newFaces[6] = JTriangle::newObject( nodes[2], nodes[9], nodes[6]);
    newFaces[7] = JTriangle::newObject( nodes[2], nodes[7], nodes[9]);

    newFaces[8] = JTriangle::newObject( nodes[0], nodes[8], nodes[10]);
    newFaces[9] = JTriangle::newObject( nodes[0], nodes[10], nodes[7]);

    newFaces[10] = JTriangle::newObject( nodes[1], nodes[6], nodes[11]);
    newFaces[11] = JTriangle::newObject( nodes[1], nodes[11], nodes[8]);

    newFaces[12] = JTriangle::newObject( nodes[6], nodes[9], nodes[11]);
    newFaces[13] = JTriangle::newObject( nodes[7], nodes[10], nodes[9]);
    newFaces[14] = JTriangle::newObject( nodes[8], nodes[11], nodes[10]);
    newFaces[15] = JTriangle::newObject( nodes[9], nodes[10], nodes[11]);

    // Some existing entities should go away ...
    for( int i = 0; i < 4; i++)
        faceneighs[i]->setStatus( JMeshEntity::REMOVE);

    tri2->setStatus(JMeshEntity::REMOVE);
    edge->setStatus(JMeshEntity::REMOVE);

    for( int i = 0; i < 4; i++) {
        JEdgePtr e = JSimplex::getEdgeOf( vertex, nodes[i+1]);
        e->setStatus(JMeshEntity::REMOVE);
    }
    vertex->setStatus(JMeshEntity::REMOVE);

    if( outmesh ) {
        outmesh->addObject( nodes[6] );
        outmesh->addObject( nodes[7] );
        outmesh->addObject( nodes[8] );
        outmesh->addObject( nodes[9] );
        outmesh->addObject( nodes[10] );
        outmesh->addObject( nodes[11] );
        outmesh->addObjects( newFaces );
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////

int JTriRefiner :: refine3( const JFacePtr &oldface)
{
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    if( oldface == nullptr) return 1;
    if( !oldface->isActive() )  return 2;
    if( oldface->getTypeID() != JFace::TRIANGLE) return 2;

    JNodePtr vcenter;
    if( oldface->hasAttribute("Steiner") ) {
        oldface->getAttribute("Steiner", vcenter);
    } else {
        vcenter = JNode::newObject();
        Point3D pc;
        oldface->getAvgXYZ( pc );
        vcenter->setXYZCoords( pc );
    }

    const JNodePtr &v0 = oldface->getNodeAt( 0 );
    const JNodePtr &v1 = oldface->getNodeAt( 1 );
    const JNodePtr &v2 = oldface->getNodeAt( 2 );

    assert(vcenter);
    if( uvCoords ) {
        Point2D uv;
        oldface->getAvgUV( uv );
        vcenter->setAttribute("UVCoords", uv);
    }

    // One new node and three new triangles are created...
    newNodes.resize(1);
    newNodes[0] = vcenter;
    newFaces.resize(3);
    newFaces[0] = JTriangle::newObject(vcenter, v0, v1);
    newFaces[1] = JTriangle::newObject(vcenter, v1, v2);
    newFaces[2] = JTriangle::newObject(vcenter, v2, v0);

    newEdges.resize(3);
    newEdges[0] = JSimplex::getEdgeOf(vcenter, v0, 1);
    newEdges[1] = JSimplex::getEdgeOf(vcenter, v1, 1);
    newEdges[2] = JSimplex::getEdgeOf(vcenter, v2, 1);

    // Old face is removed ...
    oldface->setStatus( JMeshEntity::REMOVE);

    if( outmesh ) {
        outmesh->addObjects( newNodes );
        outmesh->addObjects( newEdges );
        outmesh->addObjects( newFaces );
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
int JTriRefiner::refine4( const JFacePtr &oldface)
{
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    if( oldface == nullptr) return 1;
    if( !oldface->isActive() )  return 2;
    if( oldface->getTypeID() != JFace::TRIANGLE) return 3;

    string eAttrib = "Steiner";

    newNodes.reserve(3);
    newEdges.reserve(9);

    Point2D uv;
    // Get a new node on the edge0
    JEdgeSequence triedges(3);
    JNodeSequence edgenodes(3);
    JFaceSequence edgefaces;

    for( int i = 0; i < 3; i++) {
        triedges[i] = oldface->getEdgeAt(i);
        if( !triedges[i]->hasAttribute(eAttrib) ) {
            const JNodePtr &v0 = oldface->getNodeAt(i);
            const JNodePtr &v1 = oldface->getNodeAt(i+1);
            edgenodes[i] = JNodeGeometry::getMidNode(v0,v1);
            triedges[i]->setAttribute(eAttrib, edgenodes[i]);
            newNodes.push_back(edgenodes[i]);
            newEdges.push_back(JSimplex::getEdgeOf(v0, edgenodes[i],1));
            newEdges.push_back(JSimplex::getEdgeOf(v1, edgenodes[i],1));
            if( uvCoords) {
                triedges[i]->getAvgUV(uv);
                edgenodes[i]->setAttribute("UVCoords", uv);
            }
        }
        triedges[i]->getAttribute( eAttrib, edgenodes[i]);
        if( selective_refinement ) {
            JEdge::getRelations( triedges[i], edgefaces);
            for( size_t j = 0; j < edgefaces.size(); j++)
                inConsistentSet.insert( edgefaces[j] );
        }
    }

    const JNodePtr &v0  = oldface->getNodeAt(0);
    const JNodePtr &v1  = oldface->getNodeAt(1);
    const JNodePtr &v2  = oldface->getNodeAt(2);

    JNodePtr ev0 = edgenodes[0];
    JNodePtr ev1 = edgenodes[1];
    JNodePtr ev2 = edgenodes[2];

    // Four new triangles are generated.
    newFaces.resize(4);
    newFaces[0] = JTriangle::newObject( v0, ev0, ev2 );
    newFaces[1] = JTriangle::newObject( v1, ev1, ev0 );
    newFaces[2] = JTriangle::newObject( v2, ev2, ev1 );
    newFaces[3] = JTriangle::newObject( ev0, ev1, ev2 );

    // Three new edges are created... All are internal..
    newEdges.push_back(JSimplex::getEdgeOf(ev0,ev1,1));
    newEdges.push_back(JSimplex::getEdgeOf(ev1,ev2,1));
    newEdges.push_back(JSimplex::getEdgeOf(ev2,ev0,1));

    // Old face is removed ...
    oldface->setStatus( JMeshEntity::REMOVE);

    for( int i = 0; i < 3; i++) {
        if( triedges[i]->getNumRelations(2) == 0)
            triedges[i]->setStatus( JMeshEntity::REMOVE);
    }

    if( outmesh ) {
        outmesh->addObjects( newNodes);
        outmesh->addObjects( newEdges);
        outmesh->addObjects( newFaces);
    }

    if( selective_refinement )
        inConsistentSet.erase( oldface );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JTriRefiner::refine6( const JFacePtr &oldface)
{
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    if( oldface == nullptr) return 1;
    if( !oldface->isActive() )  return 2;
    if( oldface->getTypeID() != JFace::TRIANGLE) return 3;

    JNodePtr v0;
    Point2D uv;

    if( oldface->hasAttribute("Steiner") )
        oldface->getAttribute("Steiner", v0);
    else {
        v0 = JNode::newObject();
        Point3D pc;
        oldface->getAvgXYZ( pc );
        v0->setXYZCoords( pc );
        newNodes.push_back(v0);
        if( uvCoords) {
            oldface->getAvgUV(uv);
            v0->setAttribute("UVCoords", uv);
        }
    }

    string eAttrib = "Steiner";

    JNodeSequence edgenodes(3);
    JEdgeSequence triedges(3);
    JFaceSequence edgefaces;

    for( int j = 0; j < 3; j++) {
        triedges[j] = oldface->getEdgeAt(j);
        if( !triedges[j]->hasAttribute(eAttrib) ) {
            const JNodePtr &v0 = oldface->getNodeAt(j);
            const JNodePtr &v1 = oldface->getNodeAt(j+1);
            edgenodes[j] = JNodeGeometry::getMidNode(v0,v1);
            triedges[j]->setAttribute(eAttrib, edgenodes[j]);
            newNodes.push_back(edgenodes[j]);
            newEdges.push_back(JSimplex::getEdgeOf(v0, edgenodes[j],1));
            newEdges.push_back(JSimplex::getEdgeOf(v1, edgenodes[j],1));
            if( uvCoords) {
                triedges[j]->getAvgUV(uv);
                edgenodes[j]->setAttribute("UVCoords", uv);
            }
        }
        triedges[j]->getAttribute( eAttrib, edgenodes[j]);
        if( selective_refinement ) {
            JEdge::getRelations( triedges[j], edgefaces);
            for( size_t i = 0; i < edgefaces.size(); i++)
                 inConsistentSet.insert( edgefaces[i] );
        }
    }

    // Old face is removed ...
    oldface->setStatus( JMeshEntity::REMOVE);
    const JNodePtr &v1 = oldface->getNodeAt(0);
    const JNodePtr &v2 = oldface->getNodeAt(1);
    const JNodePtr &v3 = oldface->getNodeAt(2);

    JNodePtr ev0 = edgenodes[0];
    JNodePtr ev1 = edgenodes[1];
    JNodePtr ev2 = edgenodes[2];

    newFaces.resize(6);
    newFaces[0] = JTriangle::newObject( v0, v1,  ev0 );
    newFaces[1] = JTriangle::newObject( v0, ev0, v2 );
    newFaces[2] = JTriangle::newObject( v0, v2,  ev1 );
    newFaces[3] = JTriangle::newObject( v0, ev1, v3 );
    newFaces[4] = JTriangle::newObject( v0, v3,  ev2 );
    newFaces[5] = JTriangle::newObject( v0, ev2, v1 );

    newEdges.push_back(JSimplex::getEdgeOf(v0, ev0,1));
    newEdges.push_back(JSimplex::getEdgeOf(v0, ev1,1));
    newEdges.push_back(JSimplex::getEdgeOf(v0, ev2,1));
    newEdges.push_back(JSimplex::getEdgeOf(v0, v1, 1));
    newEdges.push_back(JSimplex::getEdgeOf(v0, v2, 1));
    newEdges.push_back(JSimplex::getEdgeOf(v0, v3, 1));

    for( int i = 0; i < 3; i++) {
        if( triedges[i]->getNumRelations(2) == 0)
            triedges[i]->setStatus( JMeshEntity::REMOVE);
    }

    if( outmesh ) {
        outmesh->addObjects( newNodes );
        outmesh->addObjects( newEdges );
        outmesh->addObjects( newFaces );
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
int JTriRefiner:: refine(const JFacePtr &face, const JEdgePtr &edge)
{
    const JNodePtr &v0 = edge->getNodeAt(0);
    const JNodePtr &v1 = edge->getNodeAt(1);

    const JNodePtr &v2 = JTriangle::getOppositeNode( face, v0, v1);
    if( v2 == nullptr) return 1;

    JNodePtr vm;
    int err = edge->getAttribute("Steiner", vm);
    if( err) {
        vm = JNodeGeometry::getMidNode( v0, v1);
        edge->setAttribute("Steiner", vm);
        if( outmesh ) outmesh->addObject(vm);
    }

    JFacePtr t1 = JTriangle::newObject( v0, vm, v2);
    JFacePtr t2 = JTriangle::newObject( vm, v1, v2);

    JEdgeSequence newedges;
    JEdgePtr   newedge;

    newedge = JSimplex::getEdgeOf(v0, vm);
    if( newedge == nullptr) {
        newedge = JEdge::newObject(v0,vm);
        newedges.push_back(newedge);
    }

    newedge = JSimplex::getEdgeOf(vm, v1);
    if( newedge == nullptr) {
        newedge = JEdge::newObject(vm,v1);
        newedges.push_back(newedge);
    }

    int bid;
    err = edge->getAttribute("Boundary", bid);
    if( !err) {
        for( const JEdgePtr &e : newedges)
            edge->setAttribute("Boundary", bid);
    }

    newedge = JEdge::newObject(vm,v2);
    newedges.push_back(newedge);

    outmesh->addObjects(newedges);
    outmesh->addObject(t1);
    outmesh->addObject(t2);

    face->setStatus( JMeshEntity::REMOVE);
    edge->setStatus( JMeshEntity::REMOVE);

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////

int JTriRefiner:: refine(const JFacePtr &face, int type)
{
    int err = 1;
    switch( type) {
    case 13:
        err = refine3( face );
        break;
    case 14:
        err = refine4( face );
        break;
    case 16:
        err = refine6( face );
        break;
    }

    if( make_consistent )
        makeConsistent();
    else
        upgradeInconsistent();

    return err;
}

////////////////////////////////////////////////////////////////////////////////
int JTriRefiner:: refine(JFaceSequence &faces2refine, int type)
{
    int nSize = faces2refine.size();

    switch( type) {
    case 13:
        for( int i = 0; i < nSize; i++)
            refine3( faces2refine[i]);
        break;
    case 14:
        for( int i = 0; i < nSize; i++)
            refine4( faces2refine[i]);
        break;
    case 16:
        for( int i = 0; i < nSize; i++)
            refine6( faces2refine[i]);
        break;
    }

    if( make_consistent )
        makeConsistent();
    else
        upgradeInconsistent();

    return 0;

}
////////////////////////////////////////////////////////////////////////////////

int JTriRefiner:: refineAll(int type)
{
    if( mesh == nullptr ) return 1;

    size_t numFaces = mesh->getSize(2);

    if(type == 13) {
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            refine3(face);
        }
    }

    if(type == 14) {
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            refine4( face);
        }
    }

    if(type == 16) {
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            refine6( face);
        }
    }

    if( make_consistent )
        makeConsistent();
    else
        upgradeInconsistent();

    mesh->deleteEdgeAttribute("Steiner");
    mesh->pruneAll();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JTriRefiner:: refineObtuse( const JFacePtr &face )
{
    if( face == nullptr )  return 1;
    if( !face->isActive()) return 2;
}
////////////////////////////////////////////////////////////////////////////////

int JTriRefiner:: refineObtuseTriangles()
{
    if( mesh == nullptr ) return 1;

    size_t numFaces = mesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        refine3( face);
    }

    mesh->deleteEdgeAttribute("Steiner");
    mesh->pruneAll();

    return 0;
}

