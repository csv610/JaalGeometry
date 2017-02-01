#include "BernHexOps.hpp"

using namespace Jaal;
using namespace std;

///////////////////////////////////////////////////////////////////////////////////
void JBernHexOps :: getCanonical17( JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.clear();

    JHexahedronPtr hex = JHexahedron::getCanonical();
    newnodes = hex->getNodes();
    newcells.push_back(hex);
}

///////////////////////////////////////////////////////////////////////////////

void JBernHexOps :: getCanonical26( JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.resize(2);

    JQuadrilateralPtr quad1 = JQuadrilateral::getCanonical(1.0);
    JQuadrilateralPtr quad2 = JQuadrilateral::getCanonical(2.0);
    JQuadrilateralPtr quad3 = JQuadrilateral::getCanonical(1.0);

    for( int i = 0; i < 4; i++) newnodes.push_back( quad1->getNodeAt(i) );
    for( int i = 0; i < 4; i++) newnodes.push_back( quad2->getNodeAt(i) );
    for( int i = 0; i < 4; i++) newnodes.push_back( quad3->getNodeAt(i) );

    Point3D p3d;
    for( int i = 0; i < 4; i++) {
        JNodePtr v  = quad1->getNodeAt(i);
        p3d = v->getXYZCoords();
        p3d[2] = -1.0;
        v->setXYZCoords( p3d );
    }

    for( int i = 0; i < 4; i++) {
        JNodePtr v  = quad3->getNodeAt(i);
        p3d = v->getXYZCoords();
        p3d[2] = 1.0;
        v->setXYZCoords( p3d );
    }
    newcells[0] = JHexahedron::newObject( quad1, quad2 );
    newcells[1] = JHexahedron::newObject( quad2, quad3 );
}

////////////////////////////////////////////////////////////////////////////////

void JBernHexOps :: getCanonical314_5( JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.resize(3);

    double dtheta = 2.0*M_PI/6;

    Point3D xyz;
    JNodePtr vtx;

    vtx = JNode::newObject();
    xyz[0] =  0.0;
    xyz[1] =  0.0;
    xyz[2] = -0.1;
    vtx->setXYZCoords(xyz);

    newnodes.push_back(vtx);
    for( int i = 0; i < 6; i++) {
        vtx = JNode::newObject();
        xyz[0] = cos( i*dtheta );
        xyz[1] = sin( i*dtheta );
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
        newnodes.push_back(vtx);
    }

    vtx = JNode::newObject();
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 1.1;
    vtx->setXYZCoords(xyz);

    newnodes.push_back(vtx);
    for( int i = 0; i < 6; i++) {
        vtx = JNode::newObject();
        xyz[0] = cos( i*dtheta );
        xyz[1] = sin( i*dtheta );
        xyz[2] = 1.0;
        vtx->setXYZCoords(xyz);
        newnodes.push_back(vtx);
    }

    JQuadrilateralPtr quad1 = JQuadrilateral::newObject( newnodes[0], newnodes[1],
                              newnodes[2], newnodes[3] );
    JQuadrilateralPtr quad2 = JQuadrilateral::newObject( newnodes[0], newnodes[3],
                              newnodes[4], newnodes[5] );
    JQuadrilateralPtr quad3 = JQuadrilateral::newObject( newnodes[0], newnodes[5],
                              newnodes[6], newnodes[1] );

    JQuadrilateralPtr quad4 = JQuadrilateral::newObject( newnodes[7], newnodes[8],
                              newnodes[9], newnodes[10] );
    JQuadrilateralPtr quad5 = JQuadrilateral::newObject( newnodes[7], newnodes[10],
                              newnodes[11], newnodes[12] );
    JQuadrilateralPtr quad6 = JQuadrilateral::newObject( newnodes[7], newnodes[12],
                              newnodes[13], newnodes[8] );

    newcells.resize(3);
    newcells[0] = JHexahedron::newObject( quad1, quad4 );
    newcells[1] = JHexahedron::newObject( quad2, quad5 );
    newcells[2] = JHexahedron::newObject( quad3, quad6 );

}

//////////////////////////////////////////////////////////////////////////////////////

void JBernHexOps :: getCanonical316_5( JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.resize(3);

    JQuadrilateralPtr quad1 = JQuadrilateral::getCanonical(1.0);
    JQuadrilateralPtr quad2 = JQuadrilateral::getCanonical(2.0);
    JQuadrilateralPtr quad3 = JQuadrilateral::getCanonical(2.0);
    JQuadrilateralPtr quad4 = JQuadrilateral::getCanonical(1.0);

    for( int i = 0; i < 4; i++) newnodes.push_back( quad1->getNodeAt(i) );
    for( int i = 0; i < 4; i++) newnodes.push_back( quad2->getNodeAt(i) );
    for( int i = 0; i < 4; i++) newnodes.push_back( quad3->getNodeAt(i) );
    for( int i = 0; i < 4; i++) newnodes.push_back( quad4->getNodeAt(i) );

    JNodePtr vtx;
    Point3D p3d;
    for( int i = 0; i < 4; i++) {
        vtx  = quad1->getNodeAt(i);
        p3d = vtx->getXYZCoords();
        p3d[2] = -1.5;
        vtx->setXYZCoords( p3d );

        vtx  = quad2->getNodeAt(i);
        p3d = vtx->getXYZCoords();
        p3d[2] = -0.5;
        vtx->setXYZCoords( p3d );

        vtx  = quad3->getNodeAt(i);
        p3d = vtx->getXYZCoords();
        p3d[2] = 0.5;
        vtx->setXYZCoords( p3d );

        vtx  = quad4->getNodeAt(i);
        p3d = vtx->getXYZCoords();
        p3d[2] = 1.5;
        vtx->setXYZCoords( p3d );
    }
    newcells[0] = JHexahedron::newObject( quad1, quad2 );
    newcells[1] = JHexahedron::newObject( quad2, quad3 );
    newcells[2] = JHexahedron::newObject( quad3, quad4 );
}

////////////////////////////////////////////////////////////////////////////////////////
void JBernHexOps :: getCanonical416_4( JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.resize(3);

    double dtheta = 2.0*M_PI/6;

    Point3D xyz;
    JNodePtr vtx;

    vtx = JNode::newObject();
    xyz[0] =  -0.2;
    xyz[1] =   0.0;
    xyz[2] =  -0.0;
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    vtx = JNode::newObject();
    xyz[0] =   0.2;
    xyz[1] =   0.0;
    xyz[2] =  -0.0;
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    for( int i = 0; i < 6; i++) {
        vtx = JNode::newObject();
        xyz[0] = cos( i*dtheta );
        xyz[1] = sin( i*dtheta );
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
        newnodes.push_back(vtx);
    }

    vtx = JNode::newObject();
    xyz[0] = -0.2;
    xyz[1] = 0.0;
    xyz[2] = 1.5;
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    vtx = JNode::newObject();
    xyz[0] = 0.2;
    xyz[1] = 0.0;
    xyz[2] = 1.5;
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    for( int i = 0; i < 6; i++) {
        vtx = JNode::newObject();
        xyz[0] = cos( i*dtheta );
        xyz[1] = sin( i*dtheta );
        xyz[2] = 1.5;
        vtx->setXYZCoords(xyz);
        newnodes.push_back(vtx);
    }

    JQuadrilateralPtr quad1 = JQuadrilateral::newObject( newnodes[1], newnodes[7],
                              newnodes[2], newnodes[3] );
    JQuadrilateralPtr quad2 = JQuadrilateral::newObject( newnodes[0], newnodes[1],
                              newnodes[3], newnodes[4] );
    JQuadrilateralPtr quad3 = JQuadrilateral::newObject( newnodes[0], newnodes[4],
                              newnodes[5], newnodes[6] );
    JQuadrilateralPtr quad4 = JQuadrilateral::newObject( newnodes[0], newnodes[6],
                              newnodes[7], newnodes[1] );

    JQuadrilateralPtr quad5 = JQuadrilateral::newObject( newnodes[9], newnodes[15],
                              newnodes[10], newnodes[11] );
    JQuadrilateralPtr quad6 = JQuadrilateral::newObject( newnodes[8], newnodes[9],
                              newnodes[11], newnodes[12] );
    JQuadrilateralPtr quad7 = JQuadrilateral::newObject( newnodes[8], newnodes[12],
                              newnodes[13], newnodes[14] );
    JQuadrilateralPtr quad8 = JQuadrilateral::newObject( newnodes[8], newnodes[14],
                              newnodes[15], newnodes[9] );

    newcells.resize(4);
    newcells[0] = JHexahedron::newObject( quad1, quad5 );
    newcells[1] = JHexahedron::newObject( quad2, quad6 );
    newcells[2] = JHexahedron::newObject( quad3, quad7 );
    newcells[3] = JHexahedron::newObject( quad4, quad8 );
}

///////////////////////////////////////////////////////////////////////////////

void JBernHexOps :: getCanonical415_4( JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.resize(3);

    double dtheta = 2.0*M_PI/6;

    Point3D xyz;
    JNodePtr vtx;

    // On the lower plane
    vtx = JNode::newObject();
    xyz[0] =  0.0;
    xyz[1] =  0.0;
    xyz[2] = -0.0;
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);
    for( int i = 0; i < 6; i++) {
        vtx = JNode::newObject();
        xyz[0] = cos( i*dtheta );
        xyz[1] = sin( i*dtheta );
        xyz[2] = 0.0;
        vtx->setXYZCoords(xyz);
        newnodes.push_back(vtx);
    }

    // On the Upper plane
    vtx = JNode::newObject();
    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0;
    vtx->setXYZCoords(xyz);

    newnodes.push_back(vtx);
    for( int i = 0; i < 6; i++) {
        vtx = JNode::newObject();
        xyz[0] = cos( i*dtheta );
        xyz[1] = sin( i*dtheta );
        xyz[2] = 1.0;
        vtx->setXYZCoords(xyz);
        newnodes.push_back(vtx);
    }

    JQuadrilateralPtr quad1 = JQuadrilateral::newObject( newnodes[0], newnodes[1],
                              newnodes[2], newnodes[3] );
    JQuadrilateralPtr quad2 = JQuadrilateral::newObject( newnodes[0], newnodes[3],
                              newnodes[4], newnodes[5] );
    JQuadrilateralPtr quad3 = JQuadrilateral::newObject( newnodes[0], newnodes[5],
                              newnodes[6], newnodes[1] );

    JQuadrilateralPtr quad4 = JQuadrilateral::newObject( newnodes[7], newnodes[8],
                              newnodes[9], newnodes[10] );
    JQuadrilateralPtr quad5 = JQuadrilateral::newObject( newnodes[7], newnodes[10],
                              newnodes[11], newnodes[12] );
    JQuadrilateralPtr quad6 = JQuadrilateral::newObject( newnodes[7], newnodes[12],
                              newnodes[13], newnodes[8] );

    newcells.resize(4);
    newcells[0] = JHexahedron::newObject( quad1, quad4 );
    newcells[1] = JHexahedron::newObject( quad2, quad5 );
    newcells[2] = JHexahedron::newObject( quad3, quad6 );

    vtx = JNode::newObject();
    xyz[0] =  0.0;
    xyz[1] =  0.0;
    xyz[2] = -0.5;
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    JNodeSequence hexnodes(8);

    hexnodes[0] = newnodes[0];
    hexnodes[1] = newnodes[1];
    hexnodes[2] = newnodes[2];
    hexnodes[3] = newnodes[3];

    hexnodes[4] = newnodes[5];
    hexnodes[5] = newnodes[6];
    hexnodes[6] = newnodes[14];
    hexnodes[7] = newnodes[4];

    newcells[3] = JHexahedron::newObject( hexnodes);
}

///////////////////////////////////////////////////////////////////////////////
int JBernHexOps::searchPattern_7_1( vector<JHexahedronPtr> &cellseq)
{
    cellseq.clear();

    if( mesh == NULL ) return 0;

    if( mesh->getAdjTable(0,3) == 0) mesh->buildRelations(0,3);

    size_t numCells = mesh->getSize(3);

    for( size_t i = currCellID; i < numCells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        bool found = 1;
        for( int j = 0; j < JHexahedron::NumNodes; j++) {
            JNodePtr vtx = cell->getNodeAt(j);
            if( vtx->getNumRelations(3) != 4 ) {
                found = 0;
                break;
            }
        }

        if( found )  {
            cellseq.push_back(JHexahedron::down_cast(cell));
            currCellID = i+1;
        }
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////


int JBernHexOps::Op1_7( const JHexahedronPtr hex, JNodeSequence &newnodes, JCellSequence &newcells)
{
    newcells.resize(7);
    double xc[8], yc[8], zc[8];

    for( int i = 0; i < 8; i++) {
        const Point3D &p3d = hex->getNodeAt(i)->getXYZCoords();
        xc[i] = p3d[0];
        yc[i] = p3d[1];
        zc[i] = p3d[2];
    }

    newnodes.resize(8);

    JNodePtr vtx;
    Point3D p3d;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( -0.5, -0.5, -0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( -0.5, -0.5, -0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( -0.5, -0.5, -0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[0]  = vtx;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( 0.5, -0.5, -0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( 0.5, -0.5, -0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( 0.5, -0.5, -0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[1]  = vtx;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( 0.5, 0.5, -0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( 0.5, 0.5, -0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( 0.5, 0.5, -0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[2]  = vtx;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( -0.5, 0.5, -0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( -0.5, 0.5, -0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( -0.5, 0.5, -0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[3]  = vtx;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( -0.5, -0.5, 0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( -0.5, -0.5, 0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( -0.5, -0.5, 0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[4]  = vtx;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( 0.5, -0.5, 0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( 0.5, -0.5, 0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( 0.5, -0.5, 0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[5]  = vtx;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( 0.5, 0.5, 0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( 0.5, 0.5, 0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( 0.5, 0.5, 0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[6]  = vtx;

    vtx = JNode::newObject();
    p3d[0] = TFI::trilinear_interpolation( -0.5, 0.5, 0.5, xc);
    p3d[1] = TFI::trilinear_interpolation( -0.5, 0.5, 0.5, yc);
    p3d[2] = TFI::trilinear_interpolation( -0.5, 0.5, 0.5, zc);
    vtx->setXYZCoords(p3d);
    newnodes[7]  = vtx;

    JHexahedronPtr hex0 = JHexahedron::newObject( newnodes);

    JNodeSequence hexnodes(8);

    hexnodes[0] =  hex->getNodeAt(0);
    hexnodes[1] =  hex->getNodeAt(4);
    hexnodes[2] =  hex->getNodeAt(7);
    hexnodes[3] =  hex->getNodeAt(3);
    hexnodes[4] =  newnodes[0];
    hexnodes[5] =  newnodes[4];
    hexnodes[6] =  newnodes[7];
    hexnodes[7] =  newnodes[3];
    JHexahedronPtr hex1 = JHexahedron::newObject( hexnodes );

    hexnodes[0] =  hex->getNodeAt(1);
    hexnodes[1] =  hex->getNodeAt(5);
    hexnodes[2] =  hex->getNodeAt(6);
    hexnodes[3] =  hex->getNodeAt(2);
    hexnodes[4] =  newnodes[1];
    hexnodes[5] =  newnodes[5];
    hexnodes[6] =  newnodes[6];
    hexnodes[7] =  newnodes[2];
    JHexahedronPtr hex2 = JHexahedron::newObject( hexnodes);

    hexnodes[0] =  hex->getNodeAt(0);
    hexnodes[1] =  hex->getNodeAt(1);
    hexnodes[2] =  hex->getNodeAt(5);
    hexnodes[3] =  hex->getNodeAt(4);
    hexnodes[4] =  newnodes[0];
    hexnodes[5] =  newnodes[1];
    hexnodes[6] =  newnodes[5];
    hexnodes[7] =  newnodes[4];
    JHexahedronPtr hex3 = JHexahedron::newObject( hexnodes);

    hexnodes[0] =  hex->getNodeAt(3);
    hexnodes[1] =  hex->getNodeAt(2);
    hexnodes[2] =  hex->getNodeAt(6);
    hexnodes[3] =  hex->getNodeAt(7);
    hexnodes[4] =  newnodes[3];
    hexnodes[5] =  newnodes[2];
    hexnodes[6] =  newnodes[6];
    hexnodes[7] =  newnodes[7];
    JHexahedronPtr hex4 = JHexahedron::newObject( hexnodes );

    hexnodes[0] =  hex->getNodeAt(0);
    hexnodes[1] =  hex->getNodeAt(1);
    hexnodes[2] =  hex->getNodeAt(2);
    hexnodes[3] =  hex->getNodeAt(3);
    hexnodes[4] =  newnodes[0];
    hexnodes[5] =  newnodes[1];
    hexnodes[6] =  newnodes[2];
    hexnodes[7] =  newnodes[3];
    JHexahedronPtr hex5 = JHexahedron::newObject( hexnodes );

    hexnodes[0] =  hex->getNodeAt(4);
    hexnodes[1] =  hex->getNodeAt(5);
    hexnodes[2] =  hex->getNodeAt(6);
    hexnodes[3] =  hex->getNodeAt(7);
    hexnodes[4] =  newnodes[4];
    hexnodes[5] =  newnodes[5];
    hexnodes[6] =  newnodes[6];
    hexnodes[7] =  newnodes[7];
    JHexahedronPtr hex6 = JHexahedron::newObject( hexnodes );

    newcells.resize(7);
    newcells[0] = hex0;
    newcells[1] = hex1;
    newcells[2] = hex2;
    newcells[3] = hex3;
    newcells[4] = hex4;
    newcells[5] = hex5;
    newcells[6] = hex6;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
#ifdef CSV

int  JBernHexOps::Op2_6( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                         JNodeSequence &newNodes, JCellSequence &newCells)
{
    newNodes.clear();
    newCells.clear();

    assert( hex1 && hex2 );

    if( !hex1->isActive() )  return 1;
    if( !hex2->isActive() )  return 1;

    JFaceSequence comm_faces;
    Cell::get_shared_entities( hex1, hex2, comm_faces);
    if( comm_faces.size() != 1) return 1;

    // Let us assume that the common face shared by two cells is face1,
    // and face0, and face2 are the other two opposite faces in hex-I and
    // hex-II respectively.
    JFacePtr face1 =  comm_faces[0];
    /*
         Face *face0 =  hex1->getOppositeFace(face1);
         Face *face2 =  hex2->getOppositeFace(face1);
    */

    // Subdivide the face1 into 5 parts. (Quad 1->5 ). The four new
    // nodes will be created..

    /*
         JFaceSequence qfaces;
         QuadRefiner::refine15(face1, newNodes, qfaces);

         newCells.resize(6);
         JNodeSequence nodes(16), hexnodes(8);

         nodes[4]  = face1->getNodeAt(0);
         nodes[5]  = face1->getNodeAt(1);
         nodes[6]  = face1->getNodeAt(2);
         nodes[7]  = face1->getNodeAt(3);

         nodes[2]  = hex1->getDiagonalNode( nodes[4] );
         nodes[3]  = hex1->getDiagonalNode( nodes[5] );
         nodes[0]  = hex1->getDiagonalNode( nodes[6] );
         nodes[1]  = hex1->getDiagonalNode( nodes[7] );

         nodes[10]  = hex2->getDiagonalNode( nodes[4] );
         nodes[11]  = hex2->getDiagonalNode( nodes[5] );
         nodes[8]   = hex2->getDiagonalNode( nodes[6] );
         nodes[9]   = hex2->getDiagonalNode( nodes[7] );

         nodes[12] = newNodes[0];
         nodes[13] = newNodes[1];
         nodes[14] = newNodes[2];
         nodes[15] = newNodes[3];

         // First Cell...
         hexnodes[0] = nodes[0];
         hexnodes[1] = nodes[1];
         hexnodes[2] = nodes[2];
         hexnodes[3] = nodes[3];
         hexnodes[4] = nodes[12];
         hexnodes[5] = nodes[13];
         hexnodes[6] = nodes[14];
         hexnodes[7] = nodes[15];
         newCells[0] = Hexahedron::newObject( hexnodes );

         hexnodes[0] = nodes[12];
         hexnodes[1] = nodes[13];
         hexnodes[2] = nodes[14];
         hexnodes[3] = nodes[15];
         hexnodes[4] = nodes[8];
         hexnodes[5] = nodes[9];
         hexnodes[6] = nodes[10];
         hexnodes[7] = nodes[11];
         newCells[1] = Hexahedron::newObject( hexnodes );

         hexnodes[0] = nodes[1];
         hexnodes[1] = nodes[13];
         hexnodes[2] = nodes[9];
         hexnodes[3] = nodes[5];
         hexnodes[4] = nodes[2];
         hexnodes[5] = nodes[14];
         hexnodes[6] = nodes[10];
         hexnodes[7] = nodes[6];
         newCells[2] = Hexahedron::newObject( hexnodes );

         hexnodes[0] = nodes[3];
         hexnodes[1] = nodes[15];
         hexnodes[2] = nodes[11];
         hexnodes[3] = nodes[7];
         hexnodes[4] = nodes[0];
         hexnodes[5] = nodes[12];
         hexnodes[6] = nodes[8];
         hexnodes[7] = nodes[4];
         newCells[3] = Hexahedron::newObject( hexnodes );

         hexnodes[0] = nodes[0];
         hexnodes[1] = nodes[12];
         hexnodes[2] = nodes[8];
         hexnodes[3] = nodes[4];
         hexnodes[4] = nodes[1];
         hexnodes[5] = nodes[13];
         hexnodes[6] = nodes[9];
         hexnodes[7] = nodes[5];
         newCells[4] = Hexahedron::newObject( hexnodes );

         hexnodes[0] = nodes[3];
         hexnodes[1] = nodes[15];
         hexnodes[2] = nodes[11];
         hexnodes[3] = nodes[7];
         hexnodes[4] = nodes[2];
         hexnodes[5] = nodes[14];
         hexnodes[6] = nodes[10];
         hexnodes[7] = nodes[6];
         newCells[5] = Hexahedron::newObject( hexnodes );
    */

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int JBernHexOps::Op314_5( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                          const JHexahedronPtr hex3, JNodeSequence &newNodes,
                          JCellSequence &newCells)
{
    newNodes.clear();
    newCells.clear();

    if( !hex1->isActive() )  return 1;
    if( !hex2->isActive() )  return 1;
    if( !hex3->isActive() )  return 1;

    JFaceSequence comm_faces;

    // The face shared by Hex1 and Hex2 is F(7,0,3,10)
    Cell::get_shared_entities( hex1, hex2, comm_faces);
    if( comm_faces.size() != 1 ) return 1;
    Quadrilateral *face70310 = Quadrilateral::down_cast( comm_faces[0] );

    // The face shared by Hex2 and Hex3 is F(7,0,5,12)
    Cell::get_shared_entities( hex2, hex3, comm_faces);
    if( comm_faces.size() != 1 ) return 1;
    Quadrilateral *face70512 = Quadrilateral::down_cast( comm_faces[0] );

    // The face shared by Hex3 and Hex1 is F(7,0,1,8)
    Cell::get_shared_entities( hex3, hex1, comm_faces);
    if( comm_faces.size() != 1 ) return 1;
    Quadrilateral *face7018 = Quadrilateral::down_cast( comm_faces[0] );

    // All the cells must be the common Edge(0,7)
    JEdgeSequence comm_edge1, comm_edge2;
    Face::get_shared_entities( face70310, face70512, comm_edge1);
    if( comm_edge1.size() != 1 ) return 1;

    Face::get_shared_entities( face70512, face7018, comm_edge2);
    if( comm_edge2.size() != 1 ) return 1;

    if( comm_edge1[0] != comm_edge2[0] ) return 1;

    JEdgePtr edge07 = comm_edge1[0];

//   Now identify all the nodes
    JNodeSequence nodes(16);

    nodes[0] = edge07->getNodeAt(0);
    nodes[7] = edge07->getNodeAt(1);

    nodes[1] = face7018->getDiagonalNode(nodes[7]);
    nodes[8] = face7018->getDiagonalNode(nodes[0]);
    nodes[2] = hex1->getDiagonalNode(nodes[7]);
    nodes[9] = hex1->getDiagonalNode(nodes[0]);

    nodes[3]  = face70310->getDiagonalNode(nodes[7]);
    nodes[10] = face70310->getDiagonalNode(nodes[0]);
    nodes[4]  = hex2->getDiagonalNode(nodes[7]);
    nodes[11] = hex2->getDiagonalNode(nodes[0]);

    nodes[5]  = face70512->getDiagonalNode(nodes[7]);
    nodes[12] = face70512->getDiagonalNode(nodes[0]);
    nodes[6]  = hex3->getDiagonalNode(nodes[7]);
    nodes[13] = hex3->getDiagonalNode(nodes[0]);

    nodes[14] = JNode::getMidNode( nodes[0], nodes[7], 0.25);
    nodes[15] = JNode::getMidNode( nodes[0], nodes[7], 0.75);

    newNodes.resize(2);
    newCells.resize(5);

    newNodes[0] = nodes[14];
    newNodes[1] = nodes[15];

    JNodeSequence hexnodes(8);

    hexnodes[0] = nodes[14];
    hexnodes[1] = nodes[6];
    hexnodes[2] = nodes[1];
    hexnodes[3] = nodes[2];
    hexnodes[4] = nodes[15];
    hexnodes[5] = nodes[13];
    hexnodes[6] = nodes[8];
    hexnodes[7] = nodes[9];
    newCells[0] = Hexahedron::newObject( hexnodes );

    hexnodes[0] = nodes[14];
    hexnodes[1] = nodes[2];
    hexnodes[2] = nodes[3];
    hexnodes[3] = nodes[4];
    hexnodes[4] = nodes[15];
    hexnodes[5] = nodes[9];
    hexnodes[6] = nodes[10];
    hexnodes[7] = nodes[11];
    newCells[1] = Hexahedron::newObject( hexnodes );

    hexnodes[0] = nodes[14];
    hexnodes[1] = nodes[4];
    hexnodes[2] = nodes[5];
    hexnodes[3] = nodes[6];
    hexnodes[4] = nodes[15];
    hexnodes[5] = nodes[11];
    hexnodes[6] = nodes[12];
    hexnodes[7] = nodes[13];
    newCells[2] = Hexahedron::newObject( hexnodes );

    hexnodes[0] = nodes[0];
    hexnodes[1] = nodes[1];
    hexnodes[2] = nodes[2];
    hexnodes[3] = nodes[3];
    hexnodes[4] = nodes[5];
    hexnodes[5] = nodes[6];
    hexnodes[6] = nodes[14];
    hexnodes[7] = nodes[4];
    newCells[3] = Hexahedron::newObject( hexnodes );

    hexnodes[0] = nodes[7];
    hexnodes[1] = nodes[8];
    hexnodes[2] = nodes[9];
    hexnodes[3] = nodes[10];
    hexnodes[4] = nodes[12];
    hexnodes[5] = nodes[13];
    hexnodes[6] = nodes[15];
    hexnodes[7] = nodes[11];
    newCells[4] = Hexahedron::newObject( hexnodes );

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JBernHexOps:: Op316_5( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                           const JHexahedronPtr hex3, JCellSequence &newCells)
{
    JHexahedron bHex;   // Bottom Hex
    JHexahedron mHex;   // Middile hex
    JHexahedron tHex;   // Top Hex

    //
    // It is unlikely or incorrect that the three input hex elements are correctly given, so check and
    // find out the order.
    //

    JFaceSequence botfaces, topfaces;

    if( mHex == NULL ) {
        Cell::get_shared_entities( hex1, hex2, topfaces);
        Cell::get_shared_entities( hex1, hex3, botfaces);
        if( (topfaces.size() == 1) && (botfaces.size() == 1 )) {
            tHex = hex2;
            mHex = hex1;
            bHex = hex3;
        }
    }

    if( mHex == NULL ) {
        Cell::get_shared_entities( hex2, hex1, topfaces);
        Cell::get_shared_entities( hex2, hex3, botfaces);
        if( (topfaces.size() == 1) && (botfaces.size() == 1 )) {
            tHex = hex1;
            mHex = hex2;
            bHex = hex3;
        }
    }

    if( mHex == NULL ) {
        Cell::get_shared_entities( hex3, hex1, topfaces);
        Cell::get_shared_entities( hex3, hex2, botfaces);
        if( (topfaces.size() == 1) && (botfaces.size() == 1 )) {
            tHex = hex1;
            mHex = hex3;
            bHex = hex2;
        }
    }

    if( mHex == NULL ) return 1;

    JNodeSequence nodes(16);

    nodes[4] = topfaces[0]->getNodeAt(0);
    nodes[5] = topfaces[0]->getNodeAt(1);
    nodes[6] = topfaces[0]->getNodeAt(2);
    nodes[7] = topfaces[0]->getNodeAt(3);

    nodes[0] = tHex->getDiagonalNode( nodes[6] );
    nodes[1] = tHex->getDiagonalNode( nodes[7] );
    nodes[2] = tHex->getDiagonalNode( nodes[4] );
    nodes[3] = tHex->getDiagonalNode( nodes[5] );

    nodes[8]  = mHex->getDiagonalNode( nodes[6] );
    nodes[9]  = mHex->getDiagonalNode( nodes[7] );
    nodes[10] = mHex->getDiagonalNode( nodes[4] );
    nodes[11] = mHex->getDiagonalNode( nodes[5] );

    nodes[12] = bHex->getDiagonalNode( nodes[10] );
    nodes[13] = bHex->getDiagonalNode( nodes[11] );
    nodes[14] = bHex->getDiagonalNode( nodes[8] );
    nodes[15] = bHex->getDiagonalNode( nodes[9] );

    newCells.resize(5);
    JNodeSequence hexnodes(8);

    hexnodes[0] = nodes[0];
    hexnodes[1] = nodes[1];
    hexnodes[2] = nodes[2];
    hexnodes[3] = nodes[3];
    hexnodes[4] = nodes[12];
    hexnodes[5] = nodes[13];
    hexnodes[6] = nodes[14];
    hexnodes[7] = nodes[15];
    newCells[0] = Hexahedron::newObject(hexnodes);

    hexnodes[0] = nodes[4];
    hexnodes[1] = nodes[8];
    hexnodes[2] = nodes[11];
    hexnodes[3] = nodes[7];
    hexnodes[4] = nodes[0];
    hexnodes[5] = nodes[12];
    hexnodes[6] = nodes[15];
    hexnodes[7] = nodes[3];
    newCells[1] = Hexahedron::newObject(hexnodes);

    hexnodes[0] = nodes[9];
    hexnodes[1] = nodes[5];
    hexnodes[2] = nodes[6];
    hexnodes[3] = nodes[10];
    hexnodes[4] = nodes[13];
    hexnodes[5] = nodes[1];
    hexnodes[6] = nodes[2];
    hexnodes[7] = nodes[14];
    newCells[2] = Hexahedron::newObject(hexnodes);

    hexnodes[0] = nodes[11];
    hexnodes[1] = nodes[10];
    hexnodes[2] = nodes[6];
    hexnodes[3] = nodes[7];
    hexnodes[4] = nodes[15];
    hexnodes[5] = nodes[14];
    hexnodes[6] = nodes[2];
    hexnodes[7] = nodes[3];
    newCells[3] = Hexahedron::newObject(hexnodes);

    hexnodes[0] = nodes[9];
    hexnodes[1] = nodes[8];
    hexnodes[2] = nodes[4];
    hexnodes[3] = nodes[5];
    hexnodes[4] = nodes[13];
    hexnodes[5] = nodes[12];
    hexnodes[6] = nodes[0];
    hexnodes[7] = nodes[1];
    newCells[4] = Hexahedron::newObject(hexnodes);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////

int JBernHexOps:: Op416_4( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                           const JHexahedronPtr hex3, const JHexahedronPtr hex4,
                           JCellSequence &newcells)
{
    vector<JHexahedronPtr> groupCells(4);

    groupCells[0] = hex1;
    groupCells[1] = hex2;
    groupCells[2] = hex3;
    groupCells[3] = hex4;

    map<JEdgePtr,vector<JHexahedronPtr> > edgecells;

    for( int j = 0; j < 4; j++) {
        JHexahedronPtr currCell = groupCells[j];
        for( int i = 0; i < Hexahedron::NumEdges; i++) {
            JEdgePtr edge = currCell->getEdgeAt(i);
            edgecells[edge].push_back(currCell);
        }
    }

    map<JEdgePtr, vector<const Hexahedron*> >::const_iterator it;

    Edge *edge08 = NULL, *edge19 = NULL;

    int ncount = 0;
    for( it = edgecells.begin(); it != edgecells.end(); ++it) {
        size_t  ncells = it->second.size();
        if( ncells == 3 ) {
            if( ncount == 0) edge19 = it->first;
            if( ncount == 1) edge08 = it->first;
            ncount++;
        }
    }
    if( ncount != 2 || edge08 == NULL || edge19 == NULL ) return 1;

    Face *face0198 = Simplex::getFaceOf( edge08->getNodeAt(0), edge08->getNodeAt(1),
                                         edge19->getNodeAt(0), edge19->getNodeAt(1), 1);

    if( face0198 == NULL ) return 2;

    vector<const Hexahedron*> neighCells;
    for( int j = 0; j < 4; j++) {
        const Hexahedron *currCell = groupCells[j];
        for( int i = 0; i < Hexahedron::NumFaces; i++) {
            Face *f = currCell->getFaceAt(i);
            if( f == face0198) neighCells.push_back(currCell);
        }
    }

    if( neighCells.size() != 2 ) return 3;

    groupCells[0] = neighCells[0];
    groupCells[1] = neighCells[1];

    for( int i = 0; i < 3; i++) {
        const Hexahedron *c = edgecells[edge08][i];
        if( c != groupCells[0]  && c != groupCells[1] ) {
            groupCells[2] = c;
            break;
        }
    }

    for( int i = 0; i < 3; i++) {
        const Hexahedron *c = edgecells[edge19][i];
        if( c != groupCells[0]  && c != groupCells[1] ) {
            groupCells[3] = c;
            break;
        }
    }

    JNodeSequence nodes(16);

    nodes[0] = edge08->getNodeAt(0);
    nodes[8] = edge08->getNodeAt(1);
    nodes[1] = edge19->getNodeAt(0);
    nodes[9] = edge19->getNodeAt(1);

    nodes[3]  =  groupCells[0]->getDiagonalNode( nodes[8] );
    nodes[4]  =  groupCells[0]->getDiagonalNode( nodes[9] );
    nodes[11] =  groupCells[0]->getDiagonalNode( nodes[0] );
    nodes[12] =  groupCells[0]->getDiagonalNode( nodes[1] );

    nodes[6]  =  groupCells[1]->getDiagonalNode( nodes[9] );
    nodes[7]  =  groupCells[1]->getDiagonalNode( nodes[8] );
    nodes[14] =  groupCells[1]->getDiagonalNode( nodes[1] );
    nodes[15] =  groupCells[1]->getDiagonalNode( nodes[0] );

    nodes[5]  =  groupCells[2]->getDiagonalNode( nodes[8] );
    nodes[13] =  groupCells[2]->getDiagonalNode( nodes[0] );

    nodes[2]  =  groupCells[3]->getDiagonalNode( nodes[9] );
    nodes[10] =  groupCells[3]->getDiagonalNode( nodes[1] );

    JNodeSequence  hexnodes(8);
    newcells.resize(4);
    hexnodes[0] = nodes[5];
    hexnodes[1] = nodes[2];
    hexnodes[2] = nodes[3];
    hexnodes[3] = nodes[4];
    hexnodes[4] = nodes[13];
    hexnodes[5] = nodes[10];
    hexnodes[6] = nodes[11];
    hexnodes[7] = nodes[12];
    newcells[0] = Hexahedron::newObject( hexnodes);

    hexnodes[0] = nodes[6];
    hexnodes[1] = nodes[7];
    hexnodes[2] = nodes[2];
    hexnodes[3] = nodes[5];
    hexnodes[4] = nodes[14];
    hexnodes[5] = nodes[15];
    hexnodes[6] = nodes[10];
    hexnodes[7] = nodes[13];
    newcells[1] = Hexahedron::newObject( hexnodes );

    hexnodes[0] = nodes[13];
    hexnodes[1] = nodes[12];
    hexnodes[2] = nodes[8];
    hexnodes[3] = nodes[14];
    hexnodes[4] = nodes[10];
    hexnodes[5] = nodes[11];
    hexnodes[6] = nodes[9];
    hexnodes[7] = nodes[15];
    newcells[2] = Hexahedron::newObject( hexnodes );

    hexnodes[0] = nodes[5];
    hexnodes[1] = nodes[4];
    hexnodes[2] = nodes[0];
    hexnodes[3] = nodes[6];
    hexnodes[4] = nodes[2];
    hexnodes[5] = nodes[3];
    hexnodes[6] = nodes[1];
    hexnodes[7] = nodes[7];
    newcells[3] = Hexahedron::newObject( hexnodes );

    return 0;
}
#endif
