#include "QuadChord.hpp"


///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    double elen = mesh->getGeometry()->getMeanEdgeLength();
    mesh->setAttribute("TargetEdgeLength", elen);
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: divide( const JEdgePtr &edge)
{
    if( edge->isBoundary() ) return;

    boost::array<JNodePtr,2> twonodes;
    const JNodePtr &v0 = edge->getNodeAt(0);
    const JNodePtr &v1 = edge->getNodeAt(1);
    Point3D pmid;

    if( !edge->hasAttribute("TwoNodes")) {
        twonodes[0] =  JNodeGeometry::getMidNode(v0, v1, 1.0/3.0);
        twonodes[1] =  JNodeGeometry::getMidNode(v0, v1, 2.0/3.0);
        edge->setAttribute("TwoNodes", twonodes);
        mesh->addObject( twonodes[0] );
        mesh->addObject( twonodes[1] );
    } else {
        edge->getAttribute("TwoNodes", twonodes);
        pmid = JNodeGeometry::getMidPoint(v0, v1, 1.0/3.0);
        twonodes[0]->setXYZCoords(pmid);
        pmid = JNodeGeometry::getMidPoint(v0, v1, 2.0/3.0);
        twonodes[1]->setXYZCoords(pmid);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: getEdgeNodes(const JNodePtr &v0, const JNodePtr &v1, JNodePtr &v2, JNodePtr &v3)
{
    JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
    assert(edge);

    boost::array<JNodePtr,2> twonodes;
    int err = edge->getAttribute("TwoNodes", twonodes);
    assert(!err);
    if( edge->getNodeAt(0) == v0 && edge->getNodeAt(1) == v1) {
        v2 = twonodes[0];
        v3 = twonodes[1];
        return;
    }

    if( edge->getNodeAt(0) == v1 && edge->getNodeAt(1) == v0) {
        v2 = twonodes[1];
        v3 = twonodes[0];
        return;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: getShapeFuncs(double r, double s, double phi[4] )
{
    phi[0] = 0.25*(1-r)*(1-s);
    phi[1] = 0.25*(1+r)*(1-s);
    phi[2] = 0.25*(1+r)*(1+s);
    phi[3] = 0.25*(1-r)*(1+s);
}
///////////////////////////////////////////////////////////////////////////////

JNodePtr JDiceQuadChord :: getNodes(const JNodePtr &v0, const JNodePtr &v1,
                                    const JNodePtr &v2, const JNodePtr &v3,
                                    double r, double s)
{
    vector<Point3D>  corners(4);
    corners[0] = v0->getXYZCoords();
    corners[1] = v1->getXYZCoords();
    corners[2] = v2->getXYZCoords();
    corners[3] = v3->getXYZCoords();

    double phi[4];
    getShapeFuncs(r, s, phi);

    Point3D xyz;
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    for( int i = 0; i < 3; i++) {
        for( int j = 0; j < 4; j++) {
            xyz[i] += phi[j]*corners[j][i];
        }
    }
    JNodePtr newnode = JNode::newObject();
    newnode->setXYZCoords(xyz);
    return newnode;
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: edge1( const JFacePtr &face)
{
    int side = -1;
    int ncount = 0;
    for( int i = 0; i < 4; i++) {
        const JEdgePtr &edge = face->getEdgeAt(i);
        if( edge->hasAttribute("TwoNodes")) {
            side = i;
            ncount++;
        }
    }
    assert( ncount == 1);

    JNodeSequence nodes(8);

    nodes[0] = face->getNodeAt(side+0);
    nodes[3] = face->getNodeAt(side+1);
    nodes[7] = face->getNodeAt(side+2);
    nodes[6] = face->getNodeAt(side+3);

    getEdgeNodes( nodes[0], nodes[3], nodes[1], nodes[2]);
    nodes[4] = getNodes( nodes[0], nodes[3], nodes[7], nodes[6], -0.50, 0.0);
    nodes[5] = getNodes( nodes[0], nodes[3], nodes[7], nodes[6],  0.50, 0.0);

    JFaceSequence newfaces(4);
    newfaces[0] = JQuadrilateral::newObject(nodes[0], nodes[1], nodes[4], nodes[6] );
    newfaces[1] = JQuadrilateral::newObject(nodes[1], nodes[2], nodes[5], nodes[4] );
    newfaces[2] = JQuadrilateral::newObject(nodes[2], nodes[3], nodes[7], nodes[5] );
    newfaces[3] = JQuadrilateral::newObject(nodes[4], nodes[5], nodes[7], nodes[6] );

    mesh->addObject(nodes[4]);
    mesh->addObject(nodes[5]);

    mesh->addObjects(newfaces);
    for( const JFacePtr &f : newfaces)
        newFaces.push_back(f);

    face->setStatus(JMeshEntity::REMOVE);
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: edge2( const JFacePtr &face)
{
    if( !face->isActive() ) return;

    int edgeid[] = {0, 0, 0, 0};
    int ncount = 0;
    for( int i = 0; i < 4; i++) {
        const JEdgePtr &edge = face->getEdgeAt(i);
        if(edge->hasAttribute("TwoNodes")) {
            edgeid[i] = 1;
            ncount++;
        }
    }
    // Only two edges must be divided ...
    assert( ncount == 2);

    const JNodePtr &v0 = face->getNodeAt(0);
    const JNodePtr &v1 = face->getNodeAt(1);
    const JNodePtr &v2 = face->getNodeAt(2);
    const JNodePtr &v3 = face->getNodeAt(3);

    JNodeSequence nodes(11);
    JFaceSequence  newfaces;

    int  id = 0;
    if( edgeid[0] && edgeid[2] ) id = 20;
    if( edgeid[1] && edgeid[3] ) id = 31;

    int side;
    switch(id)
    {
    case 20:
        nodes[0] = v0;
        nodes[3] = v1;
        nodes[4] = v3;
        nodes[7] = v2;
        getEdgeNodes( v0, v1, nodes[1], nodes[2]);
        getEdgeNodes( v3, v2, nodes[5], nodes[6]);
        newfaces.push_back(JQuadrilateral::newObject(nodes[0], nodes[1], nodes[5], nodes[4] ));
        newfaces.push_back(JQuadrilateral::newObject(nodes[1], nodes[2], nodes[6], nodes[5] ));
        newfaces.push_back(JQuadrilateral::newObject(nodes[2], nodes[3], nodes[7], nodes[6] ));
        break;
    case 31:
        nodes[0] = v0;
        nodes[1] = v1;
        nodes[6] = v3;
        nodes[7] = v2;
        getEdgeNodes( v0, v3, nodes[2], nodes[4]);
        getEdgeNodes( v1, v2, nodes[3], nodes[5]);
        newfaces.push_back( JQuadrilateral::newObject(nodes[0], nodes[1], nodes[3], nodes[2] ));
        newfaces.push_back( JQuadrilateral::newObject(nodes[2], nodes[3], nodes[5], nodes[4] ));
        newfaces.push_back( JQuadrilateral::newObject(nodes[4], nodes[5], nodes[7], nodes[6] ));
        break;
    default:
        if( edgeid[0] && edgeid[1] ) side = 0;
        if( edgeid[1] && edgeid[2] ) side = 1;
        if( edgeid[2] && edgeid[3] ) side = 2;
        if( edgeid[3] && edgeid[0] ) side = 3;
        nodes[0]  =  face->getNodeAt(side+0);
        nodes[3]  =  face->getNodeAt(side+1);
        nodes[10] =  face->getNodeAt(side+2);
        nodes[9]  =  face->getNodeAt(side+3);
        getEdgeNodes( nodes[0], nodes[3],  nodes[1], nodes[2]);
        getEdgeNodes( nodes[3], nodes[10], nodes[6], nodes[8]);
        nodes[4] = getNodes( nodes[0], nodes[3], nodes[10], nodes[9], -0.50, -0.50);
        nodes[5] = getNodes( nodes[0], nodes[3], nodes[10], nodes[9],  0.50, -0.50);
        nodes[7] = getNodes( nodes[0], nodes[3], nodes[10], nodes[9],  0.50,  0.50);
        newfaces.push_back( JQuadrilateral::newObject(nodes[0], nodes[1], nodes[4], nodes[9] ));
        newfaces.push_back( JQuadrilateral::newObject(nodes[1], nodes[2], nodes[5], nodes[4] ));
        newfaces.push_back( JQuadrilateral::newObject(nodes[2], nodes[3], nodes[6], nodes[5] ));
        newfaces.push_back( JQuadrilateral::newObject(nodes[5], nodes[6], nodes[8], nodes[7] ));
        newfaces.push_back( JQuadrilateral::newObject(nodes[7], nodes[8], nodes[10], nodes[9]));
        newfaces.push_back( JQuadrilateral::newObject(nodes[4], nodes[5], nodes[7], nodes[9]));
        mesh->addObject(nodes[4]);
        mesh->addObject(nodes[5]);
        mesh->addObject(nodes[7]);
        break;
    }
    mesh->addObjects(newfaces);
    for( const JFacePtr &f : newfaces)
        newFaces.push_back(f);
    face->setStatus(JMeshEntity::REMOVE);
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: edge3( const JFacePtr &face)
{
    if( !face->isActive() ) return;

    int side = -1;
    int numsides = 0;
    for( int i = 0; i < 4; i++) {
        const JEdgePtr &edge = face->getEdgeAt(i);
        if( !edge->hasAttribute("TwoNodes")) {
            side = i;
            numsides++;
        }
    }
    if( numsides != 1) return;

    JNodeSequence nodes(16);

    nodes[0]  = face->getNodeAt(side+1);
    nodes[3]  = face->getNodeAt(side+2);
    nodes[15] = face->getNodeAt(side+3);
    nodes[12] = face->getNodeAt(side+4);

    getEdgeNodes( nodes[0],  nodes[3],  nodes[1],  nodes[2]);
    getEdgeNodes( nodes[3],  nodes[15], nodes[7],  nodes[11]);
    getEdgeNodes( nodes[12], nodes[15], nodes[13], nodes[14]);

    nodes[4]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12], -0.5, -0.5);
    nodes[5]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12],  0.0, -0.5);
    nodes[6]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12],  0.5, -0.5);
    nodes[8]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12], -0.5,  0.5);
    nodes[9]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12],  0.0,  0.5);
    nodes[10] = getNodes( nodes[0], nodes[3], nodes[15], nodes[12],  0.5,  0.5);

    JFaceSequence newfaces(10);
    newfaces[0]  = JQuadrilateral::newObject( nodes[0],  nodes[1],  nodes[5],  nodes[4] );
    newfaces[1]  = JQuadrilateral::newObject( nodes[1],  nodes[2],  nodes[6],  nodes[5] );
    newfaces[2]  = JQuadrilateral::newObject( nodes[2],  nodes[3],  nodes[7],  nodes[6] );
    newfaces[3]  = JQuadrilateral::newObject( nodes[4],  nodes[5],  nodes[9],  nodes[8] );
    newfaces[4]  = JQuadrilateral::newObject( nodes[5],  nodes[6],  nodes[10], nodes[9] );
    newfaces[5]  = JQuadrilateral::newObject( nodes[6],  nodes[7],  nodes[11], nodes[10] );
    newfaces[6]  = JQuadrilateral::newObject( nodes[8],  nodes[9],  nodes[13], nodes[12] );
    newfaces[7]  = JQuadrilateral::newObject( nodes[9],  nodes[10], nodes[14], nodes[13] );
    newfaces[8]  = JQuadrilateral::newObject( nodes[10], nodes[11], nodes[15], nodes[14] );
    newfaces[9]  = JQuadrilateral::newObject( nodes[0],  nodes[4],  nodes[8],  nodes[12] );

    mesh->addObject( nodes[4] );
    mesh->addObject( nodes[5] );
    mesh->addObject( nodes[6] );

    mesh->addObject( nodes[8] );
    mesh->addObject( nodes[9] );
    mesh->addObject( nodes[10] );

    mesh->addObjects( newfaces );
    for( const JFacePtr &f : newfaces)
        newFaces.push_back(f);

    face->setStatus(JMeshEntity::REMOVE);
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: edge4( const JFacePtr &face)
{
    if( !face->isActive() ) return;

    JNodeSequence nodes(16);
    nodes[0]  = face->getNodeAt(0);
    nodes[3]  = face->getNodeAt(1);
    nodes[15] = face->getNodeAt(2);
    nodes[12] = face->getNodeAt(3);

    getEdgeNodes( nodes[0],   nodes[3],   nodes[1],   nodes[2]);
    getEdgeNodes( nodes[3],   nodes[15],  nodes[7],   nodes[11]);
    getEdgeNodes( nodes[12],  nodes[15],  nodes[13],  nodes[14]);
    getEdgeNodes( nodes[0],   nodes[12],  nodes[4],   nodes[8]);

    nodes[5]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12], -0.5, -0.5);
    nodes[6]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12],  0.5, -0.5);
    nodes[9]  = getNodes( nodes[0], nodes[3], nodes[15], nodes[12], -0.5,  0.5);
    nodes[10] = getNodes( nodes[0], nodes[3], nodes[15], nodes[12],  0.5,  0.5);

    JFaceSequence newfaces(9);
    int index = 0;
    for( int j = 0; j < 3; j++) {
        for( int i = 0; i < 3; i++) {
            JNodePtr v0 = nodes[ 4*j + i];
            JNodePtr v1 = nodes[ 4*j + i+1];
            JNodePtr v2 = nodes[ 4*j + i+5];
            JNodePtr v3 = nodes[ 4*j + i+4];
            newfaces[index++] = JQuadrilateral::newObject(v0,v1,v2,v3);
        }
    }

    mesh->addObject( nodes[5]  );
    mesh->addObject( nodes[6]  );
    mesh->addObject( nodes[9]  );
    mesh->addObject( nodes[10] );

    mesh->addObjects( newfaces);
    for( const JFacePtr &f : newfaces)
        newFaces.push_back(f);

    face->setStatus(JMeshEntity::REMOVE);
}

///////////////////////////////////////////////////////////////////////////////

void JDiceQuadChord :: setChord( const JQuadChordPtr &c)
{
    newFaces.clear();

    chord = c;
    JFaceSequence faces = chord->getFaces();
    size_t numfaces = faces.size();

    JEdgeSequence edges = chord->getEdges();
    size_t numedges = edges.size();

    double maxEdgeLength = 0.0;

    JMeshQuality quality;
    if( complete_dice ) {
        for( size_t i = 0; i < numedges; i++)
            divide( edges[i] );
    } else {
        mesh->getAttribute("TargetEdgeLength", maxEdgeLength);
        if( chord->isCyclic() || numedges < 6) {
            for( size_t i = 0; i < numedges; i++) {
                double val = JEdgeGeometry::getLength( edges[i] );
                if( val >  1.2*maxEdgeLength ) divide( edges[i] );
            }
        } else {
            for( size_t i = 3; i < numedges-3; i++) {
                double val = JEdgeGeometry::getLength( edges[i] );
                if( val >  1.1*maxEdgeLength ) divide( edges[i] );
            }
        }
    }

    for( const JFacePtr &face : faces) {
        int ncount = 0;
        for( int i = 0; i < 4; i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            if( edge->hasAttribute("TwoNodes") ) ncount++;
        }
        switch( ncount )
        {
        case 1:
            edge1(face);
            break;
        case 2:
            edge2(face);
            break;
        case 3:
            edge3(face);
            break;
        case 4:
            edge4(face);
            break;
        }
    }

    for( const JEdgePtr &edge : edges) {
        if( edge->hasAttribute("TwoNodes") )
            edge->setStatus(JMeshEntity::REMOVE);
    }

    mesh->pruneEdges();
    mesh->pruneFaces();
    mesh->enumerate(0);
    mesh->enumerate(1);
    mesh->enumerate(2);
    chord->clear();

    JLloydMeshOptimizer lloyd;
    lloyd.setMesh(mesh);
    lloyd.setNumIterations(100);
    lloyd.smoothAll();
}

