#include "MeshRefine.hpp"
#include "BernHexOps.hpp"

using namespace Jaal;

int JHexRefiner :: refine18( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{

    newnodes.clear();
    newcells.clear();

    const JHexahedronPtr hex = JHexahedron::down_cast( cell );
    if( hex == NULL ) return 1;
    if( !cell->isActive() ) return 1;

    JNodePtr v0 = hex->getNodeAt(0);
    JNodePtr v1 = hex->getNodeAt(1);
    JNodePtr v2 = hex->getNodeAt(2);
    JNodePtr v3 = hex->getNodeAt(3);
    JNodePtr v4 = hex->getNodeAt(4);
    JNodePtr v5 = hex->getNodeAt(5);
    JNodePtr v6 = hex->getNodeAt(6);
    JNodePtr v7 = hex->getNodeAt(7);

    JNodeSequence nodes(27);

    nodes[0] = v0;
    nodes[2] = v1;
    nodes[8] = v2;
    nodes[6] = v3;

    nodes[18] = v4;
    nodes[20] = v5;
    nodes[26] = v6;
    nodes[24] = v7;

    JEdgePtr edge;
    JNodePtr node;

    edge  = hex->getEdgeOf(v0,v1);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode(v0,v1);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[1]);

    edge  = hex->getEdgeOf(v1,v2);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v1, v2);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[5]);

    edge  = hex->getEdgeOf(v2,v3);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v2, v3);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[7]);

    edge  = hex->getEdgeOf(v3,v0);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v3, v0);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[3]);

    edge  = hex->getEdgeOf(v0,v4);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v0,v4);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[9]);

    edge  = hex->getEdgeOf(v1,v5);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode(v1,v5);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[11]);

    edge  = hex->getEdgeOf(v3,v7);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v3,v7);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[15]);

    edge  = hex->getEdgeOf(v2,v6);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode(v2, v6);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[17]);

    edge  = hex->getEdgeOf(v4,v5);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v4,v5);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[19]);

    edge  = hex->getEdgeOf(v5,v6);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v5,v6);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[23]);

    edge  = hex->getEdgeOf(v6,v7);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( v6,v7);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[25]);

    edge  = hex->getEdgeOf(v4,v7);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode(v4,v7);
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", nodes[21]);

    JFacePtr face;

    face  = hex->getFaceOf(v0,v1,v2,v3);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", nodes[4] );

    face  = hex->getFaceOf(v4,v5,v6,v7);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", nodes[22] );

    face  = hex->getFaceOf(v0,v4,v7,v3);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", nodes[12] );

    face  = hex->getFaceOf(v1,v5,v6,v2);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", nodes[14] );

    face  = hex->getFaceOf(v0,v1,v5,v4);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", nodes[10] );

    face  = hex->getFaceOf(v3,v7,v6,v2);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", nodes[16] );

    nodes[13] = JCell::getCentroid( hex );
    newnodes.push_back( nodes[13] );

    JNodeSequence hexnodes(8);

    newcells.resize(8);
    int offset, index = 0;
    for( int k = 0; k < 2; k++) {
        for( int j = 0; j < 2; j++) {
            for( int i = 0; i < 2; i++) {
                offset = 9*k + 3*j + i;
                hexnodes[0] = nodes[offset];
                hexnodes[1] = nodes[offset+1];
                offset = 9*k + 3*(j+1) + i;
                hexnodes[2] = nodes[offset+1];
                hexnodes[3] = nodes[offset];
                offset = 9*(k+1) + 3*j + i;
                hexnodes[4] = nodes[offset];
                hexnodes[5] = nodes[offset+1];
                offset = 9*(k+1) + 3*(j+1) + i;
                hexnodes[6] = nodes[offset+1];
                hexnodes[7] = nodes[offset];
                newcells[index++] =  JHexahedron::newObject( hexnodes );
            }
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refine18(const JMeshPtr &hexmesh )
{
    if( hexmesh == NULL ) return 1;

    size_t numCells = hexmesh->getSize(3);

    // Delete any "Steiner" Attributes..
    hexmesh->deleteEdgeAttribute("Steiner");
    hexmesh->deleteFaceAttribute("Steiner");

    JNodeSequence newnodes;
    JCellSequence newcells;

    for( size_t i = 0; i < numCells; i++) {
        JCellPtr cell = hexmesh->getCellAt(i);
        int err = refine18( cell, newnodes, newcells );
        if( !err ) {
            hexmesh->addObjects( newnodes );
            hexmesh->addObjects( newcells );
            cell->setStatus( JMeshEntity::REMOVE );
        }
    }
    hexmesh->pruneCells();

    hexmesh->deleteEdgeAttribute("Steiner");
    hexmesh->deleteFaceAttribute("Steiner");

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int JHexRefiner :: refine17( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.clear();
    if( !cell->isActive() ) return 1;

    const JHexahedronPtr hex = JHexahedron::down_cast( cell );
    if( hex == NULL ) return 1;

    return JBernHexOps::Op1_7(hex, newnodes, newcells ) ;
}

////////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refine17( const JMeshPtr &hexmesh )
{
    size_t numCells = hexmesh->getSize(3);

    JNodeSequence newnodes;
    JCellSequence newcells;

    for( size_t i = 0; i < numCells; i++) {
        JCellPtr cell = hexmesh->getCellAt(i);
        int err = refine17( cell, newnodes, newcells );
        if( !err ) {
            hexmesh->addObjects( newnodes );
            hexmesh->addObjects( newcells );
            cell->setStatus( JMeshEntity::REMOVE );
        }
    }
    hexmesh->pruneCells();

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode0( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr  node;
    JEdgePtr  edge;
    JFacePtr  face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v0,v1);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v0,v3);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v0,v4);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v0,v3,v7,v4);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v0,v1,v5,v4);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v0,v1,v2,v3);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v0, v8, v13, v9, v10, v12, v14, v11);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v8, v1, v2, v13, v12, v5, v6, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v10, v12, v14, v11, v4, v5, v6, v7);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v9, v13, v2, v3, v11, v14, v6, v7);
    newcells.push_back( newcell );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode1( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr node;
    JEdgePtr edge;
    JFacePtr face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v0,v1);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v1,v2);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v1,v5);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v1,v5,v6,v2);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid(face);
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v0,v1,v5,v4);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v0,v1,v2,v3);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v8, v1, v9, v13, v12, v10, v11, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v13, v9, v2, v3, v14, v11, v6, v7);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v0, v8, v13, v3, v4, v12, v14, v7);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v12, v10, v11, v14, v4, v5, v6, v7);
    newcells.push_back( newcell );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode2( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr node;
    JEdgePtr edge;
    JFacePtr face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v2,v3);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v1,v2);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v2,v6);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v1,v2,v6,v5);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid(face);
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v2,v3,v7,v6);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v0,v1,v2,v3);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v13, v9, v2, v8, v14, v11, v10, v12);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v9, v13, v0, v1, v11, v14, v4, v5);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v8, v3, v0, v13, v12, v7, v4, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v10, v12, v14, v11, v6, v7, v4, v5);
    newcells.push_back( newcell );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode3( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr node;
    JEdgePtr edge;
    JFacePtr face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v3,v0);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v3,v2);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v3,v7);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v3,v7,v6,v2);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid(face);
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v3,v7,v4,v0);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v3,v0,v1,v2);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v3, v8, v13, v9, v10, v12, v14, v11);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v9, v13, v1, v2, v11, v14, v5, v6);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v8, v0, v1, v13, v12, v4, v5, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v10, v12, v14, v11, v7, v4, v5, v6);
    newcells.push_back( newcell );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode4( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr node;
    JEdgePtr edge;
    JFacePtr  face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v4,v5);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v4,v7);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v4,v0);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v4,v0,v3,v7);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid(face);
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v4,v0,v1,v5);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v4,v5,v6,v7);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14;
    v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v4, v8, v13, v9, v10, v12, v14, v11);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v9, v13, v6, v7, v11, v14, v2, v3);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v8, v5, v6, v13, v12, v1, v2, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v10, v12, v14, v11, v0, v1, v2, v3);
    newcells.push_back( newcell );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode5( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{
    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr node;
    JEdgePtr edge;
    JFacePtr face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v5,v6);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v5,v4);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v5,v1);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v5,v1,v0,v4);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid(face);
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v5,v1,v2,v6);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v5,v6,v7,v4);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14;
    v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v5, v8, v13, v9, v10, v12, v14, v11);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v9, v13, v7, v4, v11, v14, v3, v0);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v8, v6, v7, v13, v12, v2, v3, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v10, v12, v14, v11, v1, v2, v3, v0);
    newcells.push_back( newcell );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode6( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{

    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr node;
    JEdgePtr edge;
    JFacePtr face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v6,v7);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v6,v5);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v6,v2);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v6,v2,v1,v5);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid(face);
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v6,v2,v3,v7);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v6,v7,v4,v5);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14;
    v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v6, v8, v13, v9, v10, v12, v14, v11);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v9, v13, v4, v5, v11, v14, v0, v1);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v8, v7, v4, v13, v12, v3, v0, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v10, v12, v14, v11, v2, v3, v0, v1);
    newcells.push_back( newcell );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JHexRefiner :: refineNode7( const JCellPtr &cell, JNodeSequence &newnodes, JCellSequence &newcells)
{

    if( !cell->isActive() ) return 1;

    JNodePtr v0 = cell->getNodeAt(0);
    JNodePtr v1 = cell->getNodeAt(1);
    JNodePtr v2 = cell->getNodeAt(2);
    JNodePtr v3 = cell->getNodeAt(3);
    JNodePtr v4 = cell->getNodeAt(4);
    JNodePtr v5 = cell->getNodeAt(5);
    JNodePtr v6 = cell->getNodeAt(6);
    JNodePtr v7 = cell->getNodeAt(7);

    JNodePtr node;
    JEdgePtr edge;
    JFacePtr face;
    const JHexahedronPtr hex = JHexahedron::down_cast( cell );

    JNodePtr v8;
    edge  = hex->getEdgeOf(v7,v4);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v8);

    JNodePtr v9;
    edge  = hex->getEdgeOf(v7,v6);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v9);

    JNodePtr v10;
    edge  = hex->getEdgeOf(v7,v3);
    if( !edge->hasAttribute("Steiner") ) {
        node = JNodeGeometry::getMidNode( edge->getNodeAt(0), edge->getNodeAt(1) );
        edge->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    edge->getAttribute("Steiner", v10);

    JNodePtr v11;
    face  = hex->getFaceOf(v7,v3,v2,v6);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid(face);
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v11);

    JNodePtr v12;
    face  = hex->getFaceOf(v7,v3,v0,v4);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v12);

    JNodePtr v13;
    face  = hex->getFaceOf(v7,v4,v5,v6);
    if( !face->hasAttribute("Steiner") ) {
        node = JFaceGeometry::getCentroid( face );
        face->setAttribute("Steiner", node);
        newnodes.push_back( node );
    }
    face->getAttribute("Steiner", v13);

    JNodePtr v14;
    v14 = JCell::getCentroid( hex );
    newnodes.push_back( v14 );

    JCellPtr newcell;
    newcell = JHexahedron::newObject( v7, v8, v13, v9, v10, v12, v14, v11);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v9, v13, v5, v6, v11, v14, v1, v2);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v8, v4, v5, v13, v12, v0, v1, v14);
    newcells.push_back( newcell );

    newcell = JHexahedron::newObject( v10, v12, v14, v11, v3, v0, v1, v2);
    newcells.push_back( newcell );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

