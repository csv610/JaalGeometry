#include "OBBTree.hpp"

#ifdef CSV

////////////////////////////////////////////////////////////////

int  JOBBTree :: getMaxLengthSide( const JHexahedronPtr &hex)
{
    double xlen = JNodeGeometry::getLength(hex->getNodeAt(0), hex->getNodeAt(1));
    double ylen = JNodeGeometry::getLength(hex->getNodeAt(0), hex->getNodeAt(3));
    double zlen = JNodeGeometry::getLength(hex->getNodeAt(0), hex->getNodeAt(4));

    int side = 0;
    double maxlen = xlen;

    if( ylen > maxlen ) {
        maxlen = ylen;
        side   = 1;
    }

    if( zlen > maxlen ) {
        maxlen = zlen;
        side   = 2;
    }
    return side;
}

/////////////////////////////////////////////////////////////////
int  JOBBTree  :: getPlane( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2,
                            vector<double> &plane)
{
    plane.resize(4);
    const Point3D &p0 = v0->getXYZCoords();
    const Point3D &p1 = v1->getXYZCoords();
    const Point3D &p2 = v2->getXYZCoords();
    Vec3D normal  = JTriGeometry::getNormal(p0, p1, p2);
    double d = normal[0]*p0[0] + normal[1]*p0[1] + normal[2]*p0[2];
    plane[0] = normal[0];
    plane[1] = normal[1];
    plane[2] = normal[2];
    plane[3] = -d;
}
/////////////////////////////////////////////////////////////////

int  JOBBTree :: xsplitHex(const JHexahedronPtr &hex,
                           JHexahedronPtr &lhex, JHexahedronPtr &rhex,
                           vector<double> &plane)
{
    JNodeSequence vold = hex->getNodes();

    JNodePtr v01 = JNodeGeometry::getMidNode( vold[0], vold[1] );
    JNodePtr v23 = JNodeGeometry::getMidNode( vold[2], vold[3] );
    JNodePtr v45 = JNodeGeometry::getMidNode( vold[4], vold[5] );
    JNodePtr v67 = JNodeGeometry::getMidNode( vold[6], vold[7] );

    lhex = JHexahedron::newObject( vold[0], v01, v23, vold[3],
                                  vold[4], v45, v67, vold[7] );

    boxes.push_back(lhex);

    v01 = JNodeGeometry::getMidNode( vold[0], vold[1] );
    v23 = JNodeGeometry::getMidNode( vold[2], vold[3] );
    v45 = JNodeGeometry::getMidNode( vold[4], vold[5] );
    v67 = JNodeGeometry::getMidNode( vold[6], vold[7] );

    rhex = JHexahedron::newObject( v01, vold[1], vold[2], v23,
                                  v45, vold[5], vold[6], v67 );
    boxes.push_back(rhex);

    getPlane( v01, v23, v67, plane);

    return 0;
}
/////////////////////////////////////////////////////////////////

int  JOBBTree :: ysplitHex(const JHexahedronPtr &hex,
                           JHexahedronPtr &lhex, JHexahedronPtr &rhex,
                           vector<double> &plane)
{
    JNodeSequence vold = hex->getNodes();

    JNodePtr v03 = JNodeGeometry::getMidNode( vold[0], vold[3] );
    JNodePtr v12 = JNodeGeometry::getMidNode( vold[1], vold[2] );
    JNodePtr v47 = JNodeGeometry::getMidNode( vold[4], vold[7] );
    JNodePtr v56 = JNodeGeometry::getMidNode( vold[5], vold[6] );

    lhex = JHexahedron::newObject( vold[0], vold[1], v12, v03,
                                  vold[4], vold[5], v56, v47 );
    boxes.push_back(lhex);

    v03 = JNodeGeometry::getMidNode( vold[0], vold[3] );
    v12 = JNodeGeometry::getMidNode( vold[1], vold[2] );
    v47 = JNodeGeometry::getMidNode( vold[4], vold[7] );
    v56 = JNodeGeometry::getMidNode( vold[5], vold[6] );

    rhex = JHexahedron::newObject( v03, v12, vold[2], vold[3],
                                  v47, v56, vold[6], vold[7] );
    boxes.push_back(rhex);
    getPlane( v47, v56, v12, plane);

    return 0;
}
/////////////////////////////////////////////////////////////////

int  JOBBTree :: zsplitHex(const JHexahedronPtr &hex,
                           JHexahedronPtr &lhex, JHexahedronPtr &rhex,
                           vector<double> &plane)
{
    JNodeSequence vold = hex->getNodes();

    JNodePtr v04 = JNodeGeometry::getMidNode( vold[0], vold[4] );
    JNodePtr v37 = JNodeGeometry::getMidNode( vold[3], vold[7] );
    JNodePtr v15 = JNodeGeometry::getMidNode( vold[1], vold[5] );
    JNodePtr v26 = JNodeGeometry::getMidNode( vold[2], vold[6] );

    lhex = JHexahedron::newObject( vold[0], vold[1], vold[2], vold[3],
                                  v04, v15, v26, v37 );
    boxes.push_back(lhex);

    v04 = JNodeGeometry::getMidNode( vold[0], vold[4] );
    v37 = JNodeGeometry::getMidNode( vold[3], vold[7] );
    v15 = JNodeGeometry::getMidNode( vold[1], vold[5] );
    v26 = JNodeGeometry::getMidNode( vold[2], vold[6] );

    rhex = JHexahedron::newObject( v04, v15, v26, v37,
                                  vold[4], vold[5], vold[6], vold[7] );
    boxes.push_back(rhex);
    getPlane( v04, v15, v26, plane);
    return 0;
}
/////////////////////////////////////////////////////////////////

int JOBBTree :: splitHex(const JHexahedronPtr &hex, int side,
                         JHexahedronPtr &lhex, JHexahedronPtr &rhex,
                         vector<double> &plane)
{
    switch( side)
    {
    case 0:
        xsplitHex( hex, lhex, rhex, plane);
        break;
    case 1:
        ysplitHex( hex, lhex, rhex, plane);
        break;
    case 2:
        zsplitHex( hex, lhex, rhex, plane);
        break;
    }
}

/////////////////////////////////////////////////////////////////
bool JOBBTree :: onLeftSide( const vector<double> &plane, const Point3D &query)
{
    double a = plane[0];
    double b = plane[1];
    double c = plane[2];
    double d = plane[3];

    double x = query[0];
    double y = query[1];
    double z = query[2];

    double val =  a*x + b*y + c*z + d;

    if( val > 0.0) return  0;
    return 1;
}
/////////////////////////////////////////////////////////////////
void JOBBTree :: splitHex( const NodePtr &parent)
{
    if( parent->depth == maxDepth) return;
    if( parent->hex   == nullptr)  return;
    size_t numnodes = parent->nodes.size();

    if( numnodes < 2) return;

    int maxSide = getMaxLengthSide( parent->hex);

    JHexahedronPtr  lhex, rhex;
    JHexahedronPtr  hex = parent->hex;

    vector<double>  plane;
    splitHex(hex, maxSide, lhex, rhex, plane);

    NodePtr lChild = NodePtr( new OBBNode);
    NodePtr rChild = NodePtr( new OBBNode);

    lChild->id = nCount++;
    rChild->id = nCount++;

    parent->lChild = lChild;
    parent->rChild = rChild;
    lChild->depth = parent->depth + 1;
    rChild->depth = parent->depth + 1;

    for( const JNodePtr &vtx : parent->nodes) {
        const Point3D &query = vtx->getXYZCoords();
        if( onLeftSide(plane, query ))
            lChild->nodes.push_back(vtx);
        else
            rChild->nodes.push_back(vtx);
    }

    if( lChild->nodes.empty() ) lhex->setStatus( JMeshEntity::INACTIVE);
    if( rChild->nodes.empty() ) rhex->setStatus( JMeshEntity::INACTIVE);

    parent->nodes.clear();
    parent->hex->setStatus( JMeshEntity::INACTIVE);

    JHexahedronPtr hex1 = JMeshGeometry::getMinBox(lChild->nodes);
    for( int i = 0; i < 8; i++) {
        const Point3D &xyz = hex1->getNodeAt(i)->getXYZCoords();
        lhex->getNodeAt(i)->setXYZCoords(xyz);
    }

    JHexahedronPtr hex2 = JMeshGeometry::getMinBox(rChild->nodes);
    for( int i = 0; i < 8; i++) {
        const Point3D &xyz = hex2->getNodeAt(i)->getXYZCoords();
        rhex->getNodeAt(i)->setXYZCoords(xyz);
    }

    lChild->hex = lhex;
    rChild->hex = rhex;

    splitHex( lChild );
    splitHex( rChild );
}

/////////////////////////////////////////////////////////////////

JCellSequence JOBBTree :: getBoxes(const JNodeSequence &nodes)
{
    boxes.clear();

    NodePtr root = NodePtr(new OBBNode);

    root->id    = 0;
    root->depth = 1;
    root->nodes = nodes;
    root->hex   = JMeshGeometry::getMinBox(nodes);
    boxes.push_back(root->hex);

    nCount = 1;
    splitHex(root);

    JCellSequence leafs;
    for( int i = 0; i < boxes.size(); i++) {
        if( boxes[i]->isActive() ) leafs.push_back( boxes[i] );
    }
    return leafs;
}
/////////////////////////////////////////////////////////////////
#endif
