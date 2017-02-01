#include "Doublet.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
JLogger* JDoublet :: logger = JLogger::getInstance();

///////////////////////////////////////////////////////////////////////////////

int  JDoublet :: getSize() {
    searchDoublets();
    return doublets.size();
}

///////////////////////////////////////////////////////////////////////////////

bool
JDoublet::isDoublet (const JNodePtr &vertex)
{
    assert(vertex);

    if( !vertex->isActive()   ) return 0;
    if ( vertex->isBoundary() ) return 0;

    int numfaces = vertex->getNumRelations(2);

    if( numfaces == 0) {
        cout << "Warning: Doublet not determined, may be vertex-face relations absent" << endl;
    }

    if( numfaces == 2) return 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JDoublet:: searchDoublets()
{
    doublets.clear();

    if( mesh == nullptr) return;
    if( mesh->getTopology()->getDimension() != 2 ) {
        logger->setWarn("Prensetly doublet search is only in 2D mesh " );
        return;
    }

    //
    ///////////////////////////////////////////////////////////////////////////
    // An interior doublet is a vertex, which is shared by two face neighbours.
    // They are undesirables in the quadmesh as it would mean the angle is 180
    // between some adjacent edges...
    //
    ///////////////////////////////////////////////////////////////////////////
    size_t numnodes = mesh->getSize(0);

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() && JDoublet::isDoublet(vertex)) {
            doublets.push_back(vertex);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

JNodeSequence JDoublet::getDoublets()
{
    searchDoublets();
    return doublets;
}

///////////////////////////////////////////////////////////////////////////////

JNodeSequence JDoublet::getDoublets( const JNodeSequence &nodes)
{
    JNodeSequence doublets;

    if( mesh->getTopology()->getDimension() != 2 ) {
        logger->setWarn("Prensetly doublet search is only in 2D mesh " );
        return doublets;
    }

    if( !mesh->getTopology()->isBoundaryKnown() ) {
        logger->setWarn("Doublets require boundary identifiction ");
        mesh->getTopology()->searchBoundary();
    }

    //
    ///////////////////////////////////////////////////////////////////////////
    // An interior doublet is a vertex, which is shared by two face neighbours.
    // They are undesirables in the quadmesh as it would mean the angle is 180
    // between some adjacent edges...
    //
    ///////////////////////////////////////////////////////////////////////////
//    size_t numnodes = mesh->getSize(0);

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    for( const JNodePtr &vertex: nodes) {
        if( vertex->isActive() && JDoublet::isDoublet(vertex)) {
            doublets.push_back(vertex);
        }
    }
    return doublets;
}

///////////////////////////////////////////////////////////////////////////////

int
JDoublet::remove( const JNodePtr &v0)
{
    if( !v0->isActive() ) return 1;

    JNode::getRelations(v0, faceneighs );
    if (faceneighs.size() != 2) return 1;

    JFacePtr face0 = faceneighs[0];
    JFacePtr face1 = faceneighs[1];
    assert(face0 != face1 );

    if( !face0->isActive() ) return 1;
    if( !face1->isActive() ) return 1;

    int pos = face0->getPosOf(v0);
    assert( pos >= 0);

    const JNodePtr &v1 = face0->getNodeAt( pos + 1 );
    const JNodePtr &o1 = face0->getNodeAt( pos + 2 );
    const JNodePtr &v2 = face0->getNodeAt( pos + 3 );

    assert( face1->hasNode(v0));
    assert( face1->hasNode(v1));
    assert( face1->hasNode(v2));

    pos = face1->getPosOf(v0);
    const JNodePtr &o2 = face1->getNodeAt( pos + 2 );

    JFacePtr qface = JQuadrilateral::newObject( v1, o1, v2, o2);
    mesh->addObject(qface);

    pos = face0->getPosOf(v0);
    const JEdgePtr &edge1 = face0->getEdgeAt(pos);
    const JEdgePtr &edge2 = face0->getEdgeAt(pos+3);

    face0->setStatus(JMeshEntity::REMOVE);
    face1->setStatus(JMeshEntity::REMOVE);
    edge1->setStatus(JMeshEntity::REMOVE);
    edge2->setStatus(JMeshEntity::REMOVE);

    v0->setStatus(JMeshEntity::REMOVE);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JDoublet::removeAll( )
{
    size_t numnodes = mesh->getSize(0);
    mesh->buildRelations(0,2);
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() && JDoublet::isDoublet(vertex)) {
            remove(vertex);
        }
    }
    mesh->pruneAll();
    mesh->enumerate(0);
    mesh->enumerate(1);
    mesh->enumerate(2);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr
JDoublet::insert(const JFacePtr &face, const JNodePtr &v0, const JNodePtr &v2)
{
    //Create new vertex at the center of (v0,v2)
    Point3D p3d;
    p3d = JNodeGeometry::getMidPoint(v0, v2);

    JNodePtr doublet = JNode::newObject();
    doublet->setXYZCoords(p3d);

    int pos = face->getPosOf( v0 );
    if( pos < 0) return nullptr;

    assert( v2 == face->getNodeAt( (pos+2) ) );

    JNodePtr v1 = face->getNodeAt( pos+1 );
    JNodePtr v3 = face->getNodeAt( pos+3 );

    //  Creating a doublet in the mesh changes:
    //  (1)  insert new node
    //  (2)  one old face is removed
    //  (3)  two new faces inserted.

    mesh->addObject(doublet);

    JFacePtr newquad1 = JQuadrilateral::newObject( doublet, v0, v1, v2);
    mesh->addObject(newquad1);

    JFacePtr newquad2 = JQuadrilateral::newObject( doublet, v2, v3, v0 );
    mesh->addObject(newquad2);

    face->setStatus( JMeshEntity::REMOVE );

    return doublet;
}

///////////////////////////////////////////////////////////////////////////////

