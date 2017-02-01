#include "Mesh.hpp"
#include "AllTriMeshGenerator.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "AllTetMeshGenerator.hpp"
#include "AllHexMeshGenerator.hpp"
#include "MeshRefine.hpp"

using namespace Jaal;

JMeshPtr AllTriMeshGenerator:: getSierpinski(int nlevels )
{
    JMeshPtr mesh = JMesh::newObject();

    Point3D p3d;
    JNodePtr v0 = JNode::newObject();
    p3d[0] = -1.0;
    p3d[1] =  0.0;
    p3d[2] =  0.0;
    v0->setXYZCoords(p3d);

    JNodePtr v1 = JNode::newObject();
    p3d[0] =  1.0;
    p3d[1] =  0.0;
    p3d[2] =  0.0;
    v1->setXYZCoords(p3d);

    JNodePtr v2 = JNode::newObject();
    p3d[0] =  0.0;
    p3d[1] =  2.0;
    p3d[2] =  0.0;
    v2->setXYZCoords(p3d);

    JFacePtr baseTri = JTriangle::newObject(v0,v1,v2);
    mesh->addObject(v0);
    mesh->addObject(v1);
    mesh->addObject(v2);
    mesh->addObject(baseTri);

    JNodeSequence newnodes;
    JFaceSequence newfaces, oldfaces;

    oldfaces.push_back( baseTri );

    JTriRefiner  refiner;
    refiner.setMesh(mesh);
    for( int i = 0; i < nlevels; i++) {
        size_t numfaces = mesh->getSize(2);
        for( size_t j = 0; j < numfaces; j++) {
            JFacePtr oldface = mesh->getFaceAt(j);
            if( oldface->isActive() ) {
                refiner.refine4(oldface);
                newfaces = refiner.getNewFaces();
                newfaces[3]->setStatus(JMeshEntity::REMOVE );
            }
        }
        mesh->deleteEdgeAttribute("Steiner");
    }
    return mesh;
}

///////////////////////////////////////////////////////////////////////////////

void sierpinski_refine_quad( const JFacePtr &face, JNodeSequence &newnodes, JFaceSequence &newfaces)
{
    newnodes.clear();
    newfaces.clear();

    bool hasAttrib;
    int  ori;

    JEdgePtr edge;
    JNodeSequence twonodes(2);
    JNodeSequence qnodes(16);

    JNodePtr v0 = face->getNodeAt(0);
    JNodePtr v1 = face->getNodeAt(1);
    JNodePtr v2 = face->getNodeAt(2);
    JNodePtr v3 = face->getNodeAt(3);

    qnodes[0]  = v0;
    qnodes[3]  = v1;
    qnodes[15] = v2;
    qnodes[12] = v3;

    // Nodes at first edge
    edge = face->getEdgeAt(0);
    hasAttrib = edge->hasAttribute("Steiner");
    if( hasAttrib ) {
        edge->getAttribute("Steiner", twonodes);
        assert( twonodes.size() ==  2);
        ori = face->getOrientation( edge );
        if( ori > 0) {
            qnodes[1] =  twonodes[0];
            qnodes[2] =  twonodes[1];
        }
        if( ori < 0 ) {
            qnodes[1] =  twonodes[1];
            qnodes[2] =  twonodes[0];
        }
    } else {
        twonodes[0] = JNodeGeometry::getMidNode( v0, v1, 1.0/3.0);
        twonodes[1] = JNodeGeometry::getMidNode( v0, v1, 2.0/3.0);
        edge->setAttribute("Steiner", twonodes);
        qnodes[1] =  twonodes[0];
        qnodes[2] =  twonodes[1];
        newnodes.push_back( twonodes[0] );
        newnodes.push_back( twonodes[1] );
    }

    // Nodes at second  edge
    edge = face->getEdgeAt(1);
    hasAttrib = edge->hasAttribute("Steiner");
    if( hasAttrib ) {
        edge->getAttribute("Steiner", twonodes);
        assert( twonodes.size() ==  2);
        ori = face->getOrientation( edge );
        if( ori > 0) {
            qnodes[7]  =  twonodes[0];
            qnodes[11] =  twonodes[1];
        }
        if( ori < 0 ) {
            qnodes[7]  =  twonodes[1];
            qnodes[11] =  twonodes[0];
        }
    } else {
        twonodes[0] = JNodeGeometry::getMidNode( v1, v2, 1.0/3.0);
        twonodes[1] = JNodeGeometry::getMidNode( v1, v2, 2.0/3.0);
        edge->setAttribute("Steiner", twonodes);
        qnodes[7]  =  twonodes[0];
        qnodes[11] =  twonodes[1];
        newnodes.push_back( twonodes[0] );
        newnodes.push_back( twonodes[1] );
    }

    // Nodes at Third  edge
    edge = face->getEdgeAt(2);
    hasAttrib = edge->hasAttribute("Steiner");
    if( hasAttrib ) {
        edge->getAttribute("Steiner", twonodes);
        assert( twonodes.size() ==  2);
        ori = face->getOrientation( edge );
        if( ori > 0) {
            qnodes[14]  =  twonodes[0];
            qnodes[13]  =  twonodes[1];
        }
        if( ori < 0 ) {
            qnodes[14]  =  twonodes[1];
            qnodes[13]  =  twonodes[0];
        }
    } else {
        twonodes[0] = JNodeGeometry::getMidNode( v2, v3, 1.0/3.0);
        twonodes[1] = JNodeGeometry::getMidNode( v2, v3, 2.0/3.0);
        edge->setAttribute("Steiner", twonodes);
        qnodes[14]  =  twonodes[0];
        qnodes[13]  =  twonodes[1];
        newnodes.push_back( twonodes[0] );
        newnodes.push_back( twonodes[1] );
    }

    // Nodes at Forth  edge
    edge = face->getEdgeAt(3);
    hasAttrib = edge->hasAttribute("Steiner");
    if( hasAttrib ) {
        edge->getAttribute("Steiner", twonodes);
        assert( twonodes.size() ==  2);
        ori = face->getOrientation( edge );
        if( ori > 0) {
            qnodes[8]  =  twonodes[0];
            qnodes[4]  =  twonodes[1];
        }
        if( ori < 0 ) {
            qnodes[8]  =  twonodes[1];
            qnodes[4]  =  twonodes[0];
        }
    } else {
        twonodes[0] = JNodeGeometry::getMidNode( v3, v0, 1.0/3.0);
        twonodes[1] = JNodeGeometry::getMidNode( v3, v0, 2.0/3.0);
        edge->setAttribute("Steiner", twonodes);
        qnodes[8]  =  twonodes[0];
        qnodes[4]  =  twonodes[1];
        newnodes.push_back( twonodes[0] );
        newnodes.push_back( twonodes[1] );
    }
    qnodes[5] = JNode::newObject();
    qnodes[6] = JNode::newObject();
    qnodes[9] = JNode::newObject();
    qnodes[10] = JNode::newObject();

    newnodes.push_back( qnodes[5] );
    newnodes.push_back( qnodes[6] );
    newnodes.push_back( qnodes[9] );
    newnodes.push_back( qnodes[10] );

    set_tfi_coords(1, 1, 4, 4, qnodes);
    set_tfi_coords(2, 1, 4, 4, qnodes);
    set_tfi_coords(2, 2, 4, 4, qnodes);
    set_tfi_coords(1, 2, 4, 4, qnodes);

    newfaces.resize(8);
    newfaces[0] = JQuadrilateral::newObject( qnodes[0], qnodes[1], qnodes[5], qnodes[4]);
    newfaces[1] = JQuadrilateral::newObject( qnodes[1], qnodes[2], qnodes[6], qnodes[5]);
    newfaces[2] = JQuadrilateral::newObject( qnodes[2], qnodes[3], qnodes[7], qnodes[6]);

    newfaces[3] = JQuadrilateral::newObject( qnodes[4], qnodes[5], qnodes[9], qnodes[8]);
    newfaces[4] = JQuadrilateral::newObject( qnodes[6], qnodes[7], qnodes[11], qnodes[10]);

    newfaces[5] = JQuadrilateral::newObject( qnodes[8], qnodes[9], qnodes[13], qnodes[12]);
    newfaces[6] = JQuadrilateral::newObject( qnodes[9], qnodes[10], qnodes[14], qnodes[13]);
    newfaces[7] = JQuadrilateral::newObject( qnodes[10], qnodes[11], qnodes[15], qnodes[14]);
    face->setStatus( JMeshEntity::REMOVE);
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr AllQuadMeshGenerator:: getSierpinski( int nlevels )
{
    JMeshPtr mesh = JMesh::newObject();
    mesh->setAdjTable(0,1,1);

    Point3D p3d;
    JNodePtr v0 = JNode::newObject();
    p3d[0] = -1.0;
    p3d[1] = -1.0;
    p3d[2] =  0.0;
    v0->setXYZCoords(p3d);

    JNodePtr v1 = JNode::newObject();
    p3d[0] =  1.0;
    p3d[1] = -1.0;
    p3d[2] =  0.0;
    v1->setXYZCoords(p3d);

    JNodePtr v2 = JNode::newObject();
    p3d[0] =  1.0;
    p3d[1] =  1.0;
    p3d[2] =  0.0;
    v2->setXYZCoords(p3d);

    JNodePtr v3 = JNode::newObject();
    p3d[0] = -1.0;
    p3d[1] =  1.0;
    p3d[2] =  0.0;
    v3->setXYZCoords(p3d);

    JFacePtr baseQuad = JQuadrilateral::newObject(v0,v1,v2,v3);
    mesh->addObject(v0);
    mesh->addObject(v1);
    mesh->addObject(v2);
    mesh->addObject(v3);
    mesh->addObject(baseQuad);

    JNodeSequence newnodes;
    JFaceSequence newfaces, oldfaces;

    oldfaces.push_back( baseQuad );
    mesh->deleteEdgeAttribute("Steiner");

    for( int i = 0; i < nlevels; i++) {
        size_t numfaces = mesh->getSize(2);
        for( size_t j = 0; j < numfaces; j++) {
            JFacePtr oldface = mesh->getFaceAt(j);
            if( oldface->isActive() ) {
                sierpinski_refine_quad( oldface, newnodes, newfaces);
                mesh->addObjects( newnodes );
                mesh->addObjects( newfaces );
            }
        }
        mesh->deleteEdgeAttribute("Steiner");
    }
    return mesh;
}

//////////////////////////////////////////////////////////////////////////////////////

void sierpinski_refine_tet( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.clear();

    JNodeSequence qnodes(10);

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);

    qnodes[0] = v0;
    qnodes[1] = v1;
    qnodes[2] = v2;
    qnodes[3] = v3;

    qnodes[4] = JNodeGeometry::getMidNode(v0, v1);
    qnodes[5] = JNodeGeometry::getMidNode(v0, v2);
    qnodes[6] = JNodeGeometry::getMidNode(v0, v3);

    qnodes[7] = JNodeGeometry::getMidNode(v1, v2);
    qnodes[8] = JNodeGeometry::getMidNode(v2, v3);
    qnodes[9] = JNodeGeometry::getMidNode(v1, v3);

    newnodes.resize(6);
    newnodes[0]  = qnodes[4];
    newnodes[1]  = qnodes[5];
    newnodes[2]  = qnodes[6];
    newnodes[3]  = qnodes[7];
    newnodes[4]  = qnodes[8];
    newnodes[5]  = qnodes[9];

    newcells.resize(4);
    newcells[0] = JTetrahedron::newObject( qnodes[0], qnodes[4], qnodes[5], qnodes[6] );
    newcells[1] = JTetrahedron::newObject( qnodes[1], qnodes[9], qnodes[7], qnodes[4] );
    newcells[2] = JTetrahedron::newObject( qnodes[2], qnodes[5], qnodes[7], qnodes[8] );
    newcells[3] = JTetrahedron::newObject( qnodes[3], qnodes[6], qnodes[8], qnodes[9] );
}

//////////////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllTetMeshGenerator :: getSierpinski( int nlevels )
{
    JMeshPtr mesh = JMesh::newObject();

    /*
         Cell *baseTet = Tetrahedron::getCanonical();

         mesh->addObjects(baseTet->getNodes() );
         mesh->addObject(baseTet);

         JNodeSequence newnodes;
         JCellSequence newcells, oldcells;

         oldcells.push_back( baseTet );

         for( int i = 0; i < nlevels; i++) {
              size_t numcells = mesh->getSize(3);
              for( size_t j = 0; j < numcells; j++) {
                   Cell *oldcell = mesh->getCellAt(j);
                   if( oldcell->isActive() ) {
                        sierpinski_refine_tet( oldcell, newnodes, newcells);
                        mesh->addObjects( newnodes );
                        mesh->addObjects( newcells );
                        oldcell->setStatus(JMeshEntity::REMOVE );
                   }
              }
         }
         mesh->getTopology()->collect_faces();
    */
    return mesh;
}

///////////////////////////////////////////////////////////////////////////////


JMeshPtr AllHexMeshGenerator :: getSierpinski( int nlevels )
{
    JMeshPtr mesh = JMesh::newObject();

    /*
         Cell *baseHex = Hexahedron::getCanonical();

         mesh->addObjects(baseHex->getNodes() );
         mesh->addObject(baseHex);

         JNodeSequence newnodes;
         JCellSequence newcells, oldcells;

         oldcells.push_back( baseHex );

         for( int i = 0; i < nlevels; i++) {
              size_t numcells = mesh->getSize(3);
              for( size_t j = 0; j < numcells; j++) {
                   Cell *oldcell = mesh->getCellAt(j);
                   if( oldcell->isActive() ) {
                        HexRefiner::refine(oldcell, newnodes, newcells);
                        mesh->addObjects( newnodes );
                        mesh->addObjects( newcells );
                        cout << newcells.size() << endl;
    //             newcells[13]->setStatus( MeshEntity::REMOVE );
                        oldcell->setStatus(JMeshEntity::REMOVE );
                   }
              }
         }
    */
    return mesh;
}

///////////////////////////////////////////////////////////////////////////////
