#include "PenroseTiling.hpp"

///////////////////////////////////////////////////////////////////////////////

JTrianglePtr PenroseTiling :: getCanonical( int type , double l)
{
    double angle;

    if( type < 1 || type > 2 ) return NULL;

    if( type == 1) angle = M_PI/10.0;
    if( type == 2) angle = 3.0*M_PI/10.0;

    JNodePtr v0 = JNode::newObject();
    JNodePtr v1 = JNode::newObject();
    JNodePtr v2 = JNode::newObject();
    Point3D xyz;

    xyz[0] = 0.0;
    xyz[1] = l*cos(angle);
    xyz[2] = 0.0;
    v0->setXYZCoords(xyz);

    xyz[0] = -l*sin(angle);
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    v1->setXYZCoords(xyz);

    xyz[0] = l*sin(angle);
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    v2->setXYZCoords(xyz);

    JTrianglePtr tri = JTriangle::newObject(v0,v1,v2);
    return tri;
}

///////////////////////////////////////////////////////////////////////////////

void PenroseTiling :: refine36(JFacePtr intri )
{
    if( !intri->isActive() ) return;

    JNodePtr vA = intri->getNodeAt(0);
    JNodePtr vB = intri->getNodeAt(1);
    JNodePtr vC = intri->getNodeAt(2);

    JNodePtr vP;
    JEdgePtr edge = JSimplex::getEdgeOf(vA, vB);
    if( edge->hasAttribute("Steiner") )
        edge->getAttribute("Steiner", vP);
    else {
        vP = JNodeGeometry::getMidNode(vA, vB, 1.0/goldenRatio);
        edge->setAttribute("Steiner", vP);
        mesh->addObject(vP);
    }

    JFacePtr tri1 = JTriangle::newObject(vC, vP, vB);
    tri1->setAttribute("Tile", 0);

    JFacePtr tri2 = JTriangle::newObject(vP, vC, vA);
    tri2->setAttribute("Tile", 1);

    mesh->addObject(tri1);
    mesh->addObject(tri2);
    faceQ.push_back(tri1);
    faceQ.push_back(tri2);
    intri->setStatus( JMeshEntity::REMOVE);
}

///////////////////////////////////////////////////////////////////////////////

void PenroseTiling :: refine108( JFacePtr intri)
{
    if( !intri->isActive() ) return;

    JNodePtr vA = intri->getNodeAt(0);
    JNodePtr vB = intri->getNodeAt(1);
    JNodePtr vC = intri->getNodeAt(2);

    JNodePtr vQ;
    JEdgePtr edge = JSimplex::getEdgeOf(vA, vB);
    if( edge->hasAttribute("Steiner") )
        edge->getAttribute("Steiner", vQ);
    else {
        vQ = JNodeGeometry::getMidNode(vB, vA, 1.0/goldenRatio);
        edge->setAttribute("Steiner", vQ);
        mesh->addObject(vQ);
    }

    JNodePtr vR;
    edge = JSimplex::getEdgeOf(vB, vC);
    if( edge->hasAttribute("Steiner") )
        edge->getAttribute("Steiner", vR);
    else {
        vR = JNodeGeometry::getMidNode(vB, vC, 1.0/goldenRatio);
        edge->setAttribute("Steiner", vR);
        mesh->addObject(vR);
    }


    JFacePtr tri1 = JTriangle::newObject(vR, vQ, vA);
    tri1->setAttribute("Tile", 0);

    JFacePtr tri2 = JTriangle::newObject(vQ, vR, vB);
    tri2->setAttribute("Tile", 1);

    JFacePtr tri3 = JTriangle::newObject(vR, vC, vA);
    tri3->setAttribute("Tile", 1);

    mesh->addObject(tri1);
    mesh->addObject(tri2);
    mesh->addObject(tri3);

    faceQ.push_back(tri1);
    faceQ.push_back(tri2);
    faceQ.push_back(tri3);
    intri->setStatus( JMeshEntity::REMOVE);
}

////////////////////////////////////////////////////////////////////////////////

void PenroseTiling :: refine( JFacePtr face )
{
    if( !face->isActive() ) return;

    JTrianglePtr tri  = JTriangle::down_cast(face);

    int val = 1;
    tri->getAttribute("Tile", val);

    switch(val)
    {
    case 0:
        refine36(tri);
    case 1:
        refine108(tri);
    }
}

////////////////////////////////////////////////////////////////////////////////

void PenroseTiling :: execute()
{
    assert( mesh );

    Point3D xyz;
    JNodePtr v0 = JNode::newObject();
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;
    v0->setXYZCoords(xyz);
    mesh->addObject( v0);

    double dtheta = 2.0*M_PI/10.0;
    JNodeSequence ringNodes(10);
    for( int i = 0; i < 10; i++) {
        JNodePtr vtx = JNode::newObject();
        xyz[0] = initRadius*cos(i*dtheta);
        xyz[1] = initRadius*sin(i*dtheta);
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
        mesh->addObject(vtx);
        ringNodes[i] =  vtx;
    }

    JTrianglePtr t;
    for( int i = 0; i < 10; i++) {
        JNodePtr v1 = ringNodes[i];
        JNodePtr v2 = ringNodes[(i+1)%10];
        if( i%2 == 0)
            t = JTriangle::newObject(v0,v1,v2);
        else
            t = JTriangle::newObject(v0,v2,v1);
        t->setAttribute("Tile", 0);
        mesh->addObject(t);
    }

    for( int i = 0; i < numLevels; i++) {
        size_t numfaces = mesh->getSize(2);
        for( size_t j = 0; j < numfaces; j++) {
            refine(mesh->getFaceAt(j) );
        }
    }

}

