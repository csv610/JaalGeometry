#include "Mesh.hpp"
#include "tfiblend.hpp"
#include "MeshRefine.hpp"
#include "GeomPredicates.hpp"

using namespace Jaal;

size_t JCell::NumObjectsCreated = 0;
std::map<string,string> JCell::attribInfo;

///////////////////////////////////////////////////////////////////////////////
int JCell :: registerAttribute( const string &name, const string &type)
{
    int  found = 0;

    if( type =="int"    ) found = 1;
    if( type =="char"   ) found = 1;
    if( type =="float"  ) found = 1;
    if( type =="double" ) found = 1;
    if( type =="uchar"  ) found = 1;

    if( !found) {
        cout << "Warning: invalid attribute type " << type << endl;
        return 2;
    }

    attribInfo[name] = type;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

string JCell :: getAttributeTypeName( const string &name)
{
    string str;
    if( attribInfo.find( name ) == attribInfo.end() ) return str;
    return  attribInfo[name];
}

///////////////////////////////////////////////////////////////////////////////

JCellSequence JTetrahedron::newObjects(size_t n)
{
    JCellSequence newtets;
    if( n < 1) return newtets;
    newtets.resize(n);
    for( size_t i = 0; i < n; i++)
        newtets[i] = JTetrahedron::newObject();
    return newtets;
}

///////////////////////////////////////////////////////////////////////////////
JCellSequence JHexahedron::newObjects(size_t n)
{
    JCellSequence newhexs;
    if( n < 1) return newhexs;
    newhexs.resize(n);
    for( size_t i = 0; i < n; i++)
        newhexs[i] = JHexahedron::newObject();
    return newhexs;
}
///////////////////////////////////////////////////////////////////////////////

JCellPtr JCell :: newObject( const JNodeSequence &nodeseq )
{
    JCellPtr newentity;
    int nsize = nodeseq.size();
    switch( nsize ) {
    case 4:
        newentity = JTetrahedron::newObject(nodeseq);
        break;
    case 8:
        newentity = JHexahedron::newObject(nodeseq);
        break;
    }
    return newentity;
}
///////////////////////////////////////////////////////////////////////////////

JNodePtr
JCell::getCentroid(const JCellPtr cell)
{
    Point3D pC;
    pC[0] = 0.0;
    pC[1] = 0.0;
    pC[2] = 0.0;
    int nSize = cell->getSize(0);
    for( int i = 0; i < nSize; i++) {
        const Point3D &xyz = cell->getNodeAt(i)->getXYZCoords();
        pC[0] += xyz[0];
        pC[1] += xyz[1];
        pC[2] += xyz[2];
    }
    pC[0] /= ( double)nSize;
    pC[1] /= ( double)nSize;
    pC[2] /= ( double)nSize;

    JNodePtr vc = JNode::newObject();
    vc->setXYZCoords(pC);
    return vc;
}
///////////////////////////////////////////////////////////////////////////////

JCellPtr JCell :: explode( double alpha) const
{
    JCellPtr fe = getClone();
    JNodeSequence newnodes;

    Point3D p0, pmid;
    this->getAvgXYZ(p0);
    int nn = getSize(0);
    newnodes.resize(nn);
    for( int i = 0; i < nn; i++) {
        const Point3D &p1 = getNodeAt(i)->getXYZCoords();
        pmid[0] = (1 - alpha) * p0[0] + alpha * p1[0];
        pmid[1] = (1 - alpha) * p0[1] + alpha * p1[1];
        pmid[2] = (1 - alpha) * p0[2] + alpha * p1[2];
        JNodePtr vtx = JNode::newObject();
        vtx->setXYZCoords(pmid);
        newnodes[i] = vtx;
    }
    fe->setNodes(newnodes);

    return fe;
}
///////////////////////////////////////////////////////////////////////////////

void JCell::getRelations23( const JCellPtr &cell, JCellSequence &seq)
{
    /*
     seq.clear();

       static JFaceSequence cellfaces;
       static JCellSequence neighscells;
       static JCellSet cset;

       cset.clear();

       this->getFaces( cellfaces );
       for( size_t iface = 0; iface < cellfaces.size(); iface++) {
            cellfaces[iface]->getRelations(neighscells);
            for( size_t jcell = 0; jcell < neighscells.size(); jcell++)
                 cset.insert( neighscells[jcell] );
       }

       cset.erase(this);
       if( !cset.empty() ) {
            seq.resize( cset.size() );
            std::copy( cset.begin(), cset.end(), seq.begin() );
        }
    */
    cout << "Debug exit" << endl;
    exit(0);
}
////////////////////////////////////////////////////////////////////////

JTetrahedronPtr JTetrahedron::getCanonical()
{
    Point3D p3d;
    JNodeSequence nodes(4);

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    nodes[0] = JNode::newObject();
    nodes[0]->setXYZCoords(p3d);

    p3d[0] = 1.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    nodes[1] = JNode::newObject();
    nodes[1]->setXYZCoords(p3d);

    p3d[0] = 0.0;
    p3d[1] = 1.0;
    p3d[2] = 0.0;
    nodes[2] = JNode::newObject();
    nodes[2]->setXYZCoords(p3d);

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 1.0;
    nodes[3] = JNode::newObject();
    nodes[3]->setXYZCoords(p3d);

    JTetrahedronPtr tet = JTetrahedron::newObject();
    tet->setNodes(nodes);
    return tet;
}

JMeshPtr JTetrahedron::getSchonhardt()
{
    // Object copied from: software3d.com

    JNodeSequence nodes(6);
    Point3D xyz;

    nodes[0] = JNode::newObject();
    xyz[0] = 0.0 ;
    xyz[1] = 1.0 ;
    xyz[2] = 0.82644582514053474004;
    nodes[0]->setXYZCoords(xyz);

    nodes[1] = JNode::newObject();
    xyz[0] = 0.86602540378443864676;
    xyz[1] = -0.5 ;
    xyz[2] = 0.82644582514053474004;
    nodes[1]->setXYZCoords(xyz);

    nodes[2] = JNode::newObject();
    xyz[0] = -0.86602540378443864676;
    xyz[1] = -0.5 ;
    xyz[2] =  0.82644582514053474004;
    nodes[2]->setXYZCoords(xyz);

    nodes[3] = JNode::newObject();
    xyz[0] = -1.0 ;
    xyz[1] =  0.0 ;
    xyz[2] = -0.82644582514053474004;
    nodes[3]->setXYZCoords(xyz);

    nodes[4] = JNode::newObject();
    xyz[0] = 0.5 ;
    xyz[1] = 0.86602540378443864676;
    xyz[2] = -0.82644582514053474004;
    nodes[4]->setXYZCoords(xyz);

    nodes[5] = JNode::newObject();
    xyz[0] = 0.5 ;
    xyz[1] = -0.86602540378443864676;
    xyz[2] = -0.82644582514053474004;
    nodes[5]->setXYZCoords(xyz);

    JFaceSequence faces(8);
    faces[0] =  JTriangle::newObject( nodes[0], nodes[1], nodes[5] );
    faces[1] =  JTriangle::newObject( nodes[0], nodes[5], nodes[4] );
    faces[2] =  JTriangle::newObject( nodes[0], nodes[4], nodes[2] );
    faces[3] =  JTriangle::newObject( nodes[0], nodes[2], nodes[1] );
    faces[4] =  JTriangle::newObject( nodes[3], nodes[1], nodes[2] );
    faces[5] =  JTriangle::newObject( nodes[3], nodes[2], nodes[4] );
    faces[6] =  JTriangle::newObject( nodes[3], nodes[4], nodes[5] );
    faces[7] =  JTriangle::newObject( nodes[3], nodes[5], nodes[1] );

    JMeshPtr newmesh = JMesh::newObject();
    newmesh->addObjects( nodes );
    newmesh->addObjects( faces );
    return newmesh;
}

///////////////////////////////////////////////////////////////////////////////
int JTetrahedron :: getEdgeTopology( int id, int &v0, int &v1)
{
    switch( id) {
    case 0:
        v0 = 0;
        v1 = 1;
        break;
    case 1:
        v0 = 0;
        v1 = 2;
        break;
    case 2:
        v0 = 0;
        v1 = 3;
        break;
    case 3:
        v0 = 1;
        v1 = 2;
        break;
    case 4:
        v0 = 1;
        v1 = 3;
        break;
    case 5:
        v0 = 2;
        v1 = 3;
        break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JTetrahedron :: getPosOf(const JNodePtr &v0, const JNodePtr &v1) const
{
    int pos[2];

    pos[0] =  JCell::getPosOf(v0);
    pos[1] =  JCell::getPosOf(v1);
    sort( pos, pos + 2);

    if( pos[0] == 0 &&  pos[1] == 1) return 0;
    if( pos[0] == 0 &&  pos[1] == 2) return 1;
    if( pos[0] == 0 &&  pos[1] == 3) return 2;
    if( pos[0] == 1 &&  pos[1] == 2) return 3;
    if( pos[0] == 1 &&  pos[1] == 3) return 4;
    if( pos[0] == 2 &&  pos[1] == 3) return 5;

    cout << "Warning: Invalid edge requested from tetrahedron " << endl;
    return -1;
}
///////////////////////////////////////////////////////////////////////////////

int JTetrahedron :: getPosOf(const JEdgePtr &edge) const
{
    assert( edge != NULL );
    for( int i= 0; i < NumEdges; i++) {
        JEdgePtr e = getEdgeAt(i);
        if( e == edge ) return i;
    }
    return -1;
}
///////////////////////////////////////////////////////////////////////////////

int JTetrahedron :: getPosOf(const JFacePtr &face) const
{
    assert( face != NULL );
    for( int i= 0; i < NumFaces; i++) {
        JFacePtr f = getFaceAt(i);
        if( f == face ) return i;
    }
    return -1;
}


///////////////////////////////////////////////////////////////////////////////

int JTetrahedron :: getFaceTopology(int id, int &v0, int &v1, int &v2)
{
    assert( id >= 0 && id < 4);

    switch(id) {
    case 0:
        v0 = 1;
        v1 = 2;
        v2 = 3;
        break;
    case 1:
        v0 = 0;
        v1 = 3;
        v2 = 2;
        break;
    case 2:
        v0 = 0;
        v1 = 1;
        v2 = 3;
        break;
    case 3:
        v0 = 0;
        v1 = 2;
        v2 = 1;
        break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JEdgePtr JTetrahedron :: getEdgeAt( int i ) const
{
    int i0, i1;
    getEdgeTopology(i, i0, i1);
    JNodePtr v0 =  getNodeAt(i0);
    JNodePtr v1 =  getNodeAt(i1);
    return JSimplex::getEdgeOf(v0, v1, 1);
}

///////////////////////////////////////////////////////////////////////////////

JEdgePtr JTetrahedron :: getEdgeAt( int pos, int &ori ) const
{
    int i0, i1;
    getEdgeTopology(pos, i0, i1);
    JNodePtr v0 = getNodeAt(i0);
    JNodePtr v1 = getNodeAt(i1);
    JEdgePtr edge = JSimplex::getEdgeOf(v0, v1, 1);
    ori = edge->getOrientation(v0,v1);
    return edge;
}

///////////////////////////////////////////////////////////////////////////////
bool JTetrahedron :: hasFaceAt( int pos ) const
{
    bool build = 0;
    int i0, i1, i2;
    getFaceTopology(pos, i0, i1, i2);
    JNodePtr v0 =  getNodeAt(i0);
    JNodePtr v1 =  getNodeAt(i1);
    JNodePtr v2 =  getNodeAt(i2);
    JFacePtr f =  JSimplex::getFaceOf(v0, v1, v2, build);
    if( f ) return 1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
JFacePtr JTetrahedron :: getFaceAt( int pos ) const
{
    bool build = 1;
    int i0, i1, i2;
    getFaceTopology(pos, i0, i1, i2);
    JNodePtr v0 =  getNodeAt(i0);
    JNodePtr v1 =  getNodeAt(i1);
    JNodePtr v2 =  getNodeAt(i2);
    return JSimplex::getFaceOf(v0, v1, v2, build);
}

///////////////////////////////////////////////////////////////////////////////

JFacePtr JTetrahedron :: getFaceAt( int i, int &ori ) const
{
    bool build = 1;
    int i0, i1, i2;
    getFaceTopology(i, i0, i1, i2);
    JNodePtr v0 =  getNodeAt(i0);
    JNodePtr v1 =  getNodeAt(i1);
    JNodePtr v2 =  getNodeAt(i2);
    JFacePtr fs =    JSimplex::getFaceOf(v0, v1, v2, build);
    ori  =  fs->getOrientation(v0,v1,v2);
    return fs;
}

///////////////////////////////////////////////////////////////////////////////
int JTetrahedron :: build_lower_entities( int dim )
{
    switch(dim) {
    case 1:
        for( int i = 0; i < JTetrahedron::NumEdges; i++)
            this->getEdgeAt(i);
        break;
    case 2:
        for( int i = 0; i < JTetrahedron::NumFaces; i++)
            this->getFaceAt(i);
        break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JTetrahedron ::  getOrientation( const JEdgePtr &edge) const
{
    assert( edge  );

    int i0, i1;
    for( int i = 0; i < NumEdges; i++) {
        if( getEdgeAt(i) == edge ) {
            getEdgeTopology(i, i0, i1);
            JNodePtr ve0 = getNodeAt(i0);
            JNodePtr ve1 = getNodeAt(i1);
            JNodePtr va0 = edge->getNodeAt(0);
            JNodePtr va1 = edge->getNodeAt(1);
            if( ve0 == va0 &&  ve1 == va1 ) return  1;
            if( ve0 == va1 &&  ve1 == va0 ) return -1;
            return 0;
        }
    }
    cout << "Warning: Invalid edge search  in a tetrahedra " << endl;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JTetrahedron ::  getOrientation( const JFacePtr &face) const
{
    assert( face );
    int i0, i1, i2;
    for( int i = 0; i < NumFaces; i++) {
        if( getFaceAt(i) == face ) {
            getFaceTopology(i, i0, i1, i2);
            JNodePtr v0 = getNodeAt(i0);
            JNodePtr v1 = getNodeAt(i1);
            JNodePtr v2 = getNodeAt(i2);
            return face->getOrientation(v0,v1,v2);
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JTetrahedron :: getPosOf( const JNodePtr &v0, const JNodePtr &v1 , const JNodePtr &v2) const
{
    int pos[3];

    pos[0] =  JCell::getPosOf(v0);
    pos[1] =  JCell::getPosOf(v1);
    pos[2] =  JCell::getPosOf(v2);
    sort( pos, pos + 3);

    if( pos[0] == 1 &&  pos[1] == 2  && pos[2] == 3 ) return 0;
    if( pos[0] == 0 &&  pos[1] == 2  && pos[2] == 3 ) return 1;
    if( pos[0] == 0 &&  pos[1] == 1  && pos[2] == 3 ) return 2;
    if( pos[0] == 0 &&  pos[1] == 1  && pos[2] == 2 ) return 3;

    cout << "Warning: Invalid face requested from tetrahedron " << endl;
    return -1;
}

///////////////////////////////////////////////////////////////////////////////
JEdgePtr JTetrahedron :: getEdgeOf( const JNodePtr &v0, const JNodePtr &v1) const
{
    JEdgePtr nullPtr;
    if( !this->hasNode(v0) ) return nullPtr;
    if( !this->hasNode(v1) ) return nullPtr;
    return JSimplex::getEdgeOf(v0,v1, 1);
}

///////////////////////////////////////////////////////////////////////////////
JFacePtr JTetrahedron :: getFaceOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2) const
{
    JFacePtr nullPtr;
    if( !this->hasNode(v0) ) return nullPtr;
    if( !this->hasNode(v1) ) return nullPtr;
    if( !this->hasNode(v2) ) return nullPtr;
    return JSimplex::getFaceOf(v0,v1,v2,1);
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JTetrahedron :: getEdges() const
{
    JEdgeSequence edges(NumEdges);
    int p0, p1;
    for( int i = 0; i < NumEdges; i++) {
        getEdgeTopology(i, p0, p1);
        JNodePtr v0 = nodes[p0];
        JNodePtr v1 = nodes[p1];
        edges[i] = JSimplex::getEdgeOf(v0,v1,1);
        assert( edges[i] );
    }

    return edges;
}

///////////////////////////////////////////////////////////////////////////////

JFaceSequence JTetrahedron :: getFaces()  const
{
    JFaceSequence faces(NumFaces);

    int p0, p1, p2;
    for( int i = 0; i < 4; i++) {
        getFaceTopology( i, p0, p1, p2);
        JNodePtr v0 = getNodeAt(p0);
        JNodePtr v1 = getNodeAt(p1);
        JNodePtr v2 = getNodeAt(p2);
        faces[i]  = JSimplex::getFaceOf(v0,v1,v2,1);
    }
    return faces;
}

///////////////////////////////////////////////////////////////////////////////
int JTetrahedron :: getHexCells( vector<JHexahedron> &) const
{
    /*
         JNodeSequence v(15), vconn;

         v[0] = this->getNodeAt(0);
         v[1] = this->getNodeAt(1);
         v[2] = this->getNodeAt(2);
         v[3] = this->getNodeAt(3);

         JNodePtr vhash;
         JEdgePtr   existedge;

         vhash = min( v[0], v[1] );
         existedge =  vhash->getEdgeOf(v[0], v[1] );
         if( existedge == NULL ) return 1;
         existedge->getAttribute("Steiner", v[4] );
         if( v[4] == NULL ) return 1;

         vhash = min( v[0], v[2] );
         existedge =  vhash->getEdgeOf(v[0], v[2]);
         if( existedge == NULL ) return 1;
         existedge->getAttribute("Steiner", v[5] );
         if( v[5] == NULL ) return 1;

         vhash = min( v[0], v[3] );
         existedge =  vhash->getEdgeOf(v[0], v[3]);
         if( existedge == NULL ) return 1;
         existedge->getAttribute("Steiner", v[6] );
         if( v[6] == NULL ) return 1;

         vhash = min( v[1], v[2] );
         existedge =  vhash->getEdgeOf(v[1], v[2]);
         if( existedge == NULL ) return 1;
         existedge->getAttribute("Steiner", v[7] );
         if( v[7] == NULL ) return 1;

         vhash = min( v[2], v[3] );
         existedge =  vhash->getEdgeOf(v[2], v[3]);
         if( existedge == NULL ) return 1;
         existedge->getAttribute("Steiner", v[8] );
         if( v[8] == NULL ) return 1;

         vhash = min( v[1], v[3] );
         existedge =  vhash->getEdgeOf(v[1], v[3]);
         if( existedge == NULL ) return 1;
         existedge->getAttribute("Steiner", v[9] );
         if( v[9] == NULL ) return 1;

         JFacePtr existface;

         vconn.resize(3);
         vconn[0] = v[0];
         vconn[1] = v[1];
         vconn[2] = v[2];
         vhash = *min_element( vconn );
         existface = vhash->getFaceOf( vconn );
         if( existface == NULL ) return 1;
         existface->getAttribute("Steiner", v[10] );
         if( v[10] == NULL ) return 1;

         vconn[0] = v[0];
         vconn[1] = v[2];
         vconn[2] = v[3];
         vhash = *min_element( vconn );
         existface = vhash->getFaceOf( vconn );
         if( existface == NULL ) return 1;
         existface->getAttribute("Steiner", v[11] );
         if( v[11] == NULL ) return 1;

         vconn[0] = v[0];
         vconn[1] = v[1];
         vconn[2] = v[3];
         vhash = *min_element( vconn );
         existface = vhash->getFaceOf( vconn );
         if( existface == NULL ) return 1;
         existface->getAttribute("Steiner", v[12] );
         if( v[12] == NULL ) return 1;

         vconn[0] = v[1];
         vconn[1] = v[2];
         vconn[2] = v[3];
         vhash = *min_element( vconn) );
         existface = vhash->getFaceOf( vconn );
         if( existface == NULL ) return 1;
         existface->getAttribute("Steiner", v[13] );
         if( v[13] == NULL ) return 1;

         this->getAttribute("Steiner", v[14] );

         vconn.resize(8);
         hex.resize(4);
         vconn[0] = v[6];
         vconn[1] = v[12];
         vconn[2] = v[4];
         vconn[3] = v[0];
         vconn[4] = v[11];
         vconn[5] = v[14];
         vconn[6] = v[10];
         vconn[7] = v[5];
         hex[0].setNodes(vconn);

         vconn[0] = v[9];
         vconn[1] = v[1];
         vconn[2] = v[4];
         vconn[3] = v[12];
         vconn[4] = v[13];
         vconn[5] = v[7];
         vconn[6] = v[10];
         vconn[7] = v[14];
         hex[1].setNodes(vconn);

         vconn[0] = v[11];
         vconn[1] = v[14];
         vconn[2] = v[10];
         vconn[3] = v[5];
         vconn[4] = v[8];
         vconn[5] = v[13];
         vconn[6] = v[7];
         vconn[7] = v[2];
         hex[2].setNodes(vconn);

         vconn[0] = v[3];
         vconn[1] = v[9];
         vconn[2] = v[12];
         vconn[3] = v[6];
         vconn[4] = v[8];
         vconn[5] = v[13];
         vconn[6] = v[14];
         vconn[7] = v[11];
         hex[3].setNodes(vconn);
    */

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
double TetGeometry :: getRegularTetrahedronArea( double a )
{
    return sqrt(3.0)*a*a;
}

double TetGeometry :: getRegularTetrahedronHeight( double a )
{
    return sqrt(6.0)*a/3.0;
}

double TetGeometry :: getRegularTetrahedronVolume( double a )
{
    return sqrt(2.0)*a*a*a/12.0;
}

////////////////////////////////////////////////////////////////////////////////
void TetGeometry :: getDihedralAngles( const JCellPtr &tet, vector<double> &angle)
{

}
////////////////////////////////////////////////////////////////////////////////
double TetGeometry :: getMinDihedralAngle( const JCellPtr &tet)
{
    double angle = 0.0;
    return angle;
}
////////////////////////////////////////////////////////////////////////////////
double TetGeometry :: getMaxDihedralAngle( const JCellPtr &tet)
{
    double angle = 360.0;
    return angle;
}

Point3D TetGeometry :: getRandomPoint( const JCellPtr &cell)
{
    Point4D bary;

    double highval = 1.0;
    bary[0] = JMath::random_value(0.0, highval);
    highval -= bary[0];

    bary[1] = JMath::random_value(0.0, highval);
    highval -= bary[1];

    bary[2] = JMath::random_value(0.0, highval);
    highval -= bary[2];

    bary[3] = JMath::random_value(0.0, highval);
    highval -= bary[3];

    double sum = bary[0] + bary[1] + bary[2] + bary[3];

    bary[0] /= (double) sum;
    bary[1] /= (double) sum;
    bary[2] /= (double) sum;
    bary[3] /= (double) sum;

    assert( bary[0] >= 0.0 && bary[0] <= 1.0);
    assert( bary[1] >= 0.0 && bary[1] <= 1.0);
    assert( bary[2] >= 0.0 && bary[2] <= 1.0);
    assert( bary[3] >= 0.0 && bary[3] <= 1.0);

    Point3D xyz;
    TetGeometry::getXYZCoordinates(cell, bary, xyz);
    return xyz;
}

////////////////////////////////////////////////////////////////////////////////

bool TetGeometry :: isDegenerate( const JCellPtr &cell)
{
    assert(cell->getSize(0) == 4);
    const Point3D &pa = cell->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = cell->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = cell->getNodeAt(2)->getXYZCoords();
    const Point3D &pd = cell->getNodeAt(3)->getXYZCoords();

    double vol = orient3dexact( &pa[0], &pb[0], &pc[0], &pd[0]);
    if( vol == 0.0) return 1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int TetGeometry :: getOrientation( const JCellPtr &cell)
{
    assert(cell->getSize(0) == 4);
    const Point3D &pa = cell->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = cell->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = cell->getNodeAt(2)->getXYZCoords();
    const Point3D &pd = cell->getNodeAt(3)->getXYZCoords();

    double vol = orient3dexact( &pb[0], &pc[0], &pd[0], &pa[0]);
    if( vol < 0 ) return -1;
    if( vol > 0 ) return +1;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

bool TetGeometry :: isInside( const JCellPtr &cell, const Point3D &queryPoint, bool include_boundary)
{
    assert( cell ) ;

    assert(cell->getSize(0) == 4);
    const Point3D &pa = cell->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = cell->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = cell->getNodeAt(2)->getXYZCoords();
    const Point3D &pd = cell->getNodeAt(3)->getXYZCoords();

    double vol;
    vol = orient3dexact( &pa[0], &pc[0], &pb[0], &queryPoint[0]);
    if( vol < 0 ) return 0;

    vol = orient3dexact( &pa[0], &pd[0], &pc[0], &queryPoint[0]);
    if( vol < 0 ) return 0;

    vol = orient3dexact( &pa[0], &pb[0], &pd[0], &queryPoint[0]);
    if( vol < 0 ) return 0;

    vol = orient3dexact( &pb[0], &pc[0], &pd[0], &queryPoint[0]);
    if( vol < 0 ) return 0;

    return 1;
}
///////////////////////////////////////////////////////////////////////////////
bool TetGeometry :: isOutside( const JCellPtr &cell, const Point3D &queryPoint, bool include_boundary)
{
    assert( cell ) ;

    assert(cell->getSize(0) == 4);
    const Point3D &pa = cell->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = cell->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = cell->getNodeAt(2)->getXYZCoords();
    const Point3D &pd = cell->getNodeAt(3)->getXYZCoords();

    double vol;
    vol = orient3dexact( &pa[0], &pc[0], &pb[0], &queryPoint[0]);
    if( vol < 0 ) return 1;

    vol = orient3dexact( &pa[0], &pd[0], &pc[0], &queryPoint[0]);
    if( vol < 0 ) return 1;

    vol = orient3dexact( &pa[0], &pb[0], &pd[0], &queryPoint[0]);
    if( vol < 0 ) return 1;

    vol = orient3dexact( &pb[0], &pc[0], &pd[0], &queryPoint[0]);
    if( vol < 0 ) return 1;

    return 0;

}
///////////////////////////////////////////////////////////////////////////////

bool TetGeometry :: isOnBoundary( const JCellPtr &cell, const Point3D &queryPoint)
{
    assert(cell->getSize(0) == 4);
    const Point3D &pa = cell->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = cell->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = cell->getNodeAt(2)->getXYZCoords();
    const Point3D &pd = cell->getNodeAt(3)->getXYZCoords();

    double vol;
    vol = orient3dexact( &pa[0], &pc[0], &pb[0], &queryPoint[0]);
    if( vol == 0 ) return 1;

    vol = orient3dexact( &pa[0], &pd[0], &pc[0], &queryPoint[0]);
    if( vol == 0 ) return 1;

    vol = orient3dexact( &pa[0], &pb[0], &pd[0], &queryPoint[0]);
    if( vol == 0 ) return 1;

    vol = orient3dexact( &pb[0], &pc[0], &pd[0], &queryPoint[0]);
    if( vol == 0 ) return 1;

    return 0;
}

int TetGeometry::getCircumSphere( const JCellPtr &c, Point3D &center, double &radius)
{
    return 1;
}
///////////////////////////////////////////////////////////////////////////////////

int TetGeometry :: getInSphere( const JCellPtr &c, Point3D &center, double &radius)
{
    return 1;
}

///////////////////////////////////////////////////////////////////////////////////

int TetGeometry :: getBaryCoordinates( const JCellPtr &cell, const Point3D &xyz, Point4D &baryCoords)
{
    assert(cell->getSize(0) == 4);

    const Point3D &pa = cell->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = cell->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = cell->getNodeAt(2)->getXYZCoords();
    const Point3D &pd = cell->getNodeAt(3)->getXYZCoords();

    double vol = orient3d( &pa[0], &pc[0], &pb[0], &pd[0]);
    if( fabs(vol) < 1.0E-10) return 1;

    double vola = orient3d( &pb[0], &pc[0], &pd[0], &xyz[0]);
//   if( vola < 0.0) cout << "Warning: Barycentric coordinate is negative " << vola << endl;
    double volb = orient3d( &pa[0], &pd[0], &pc[0], &xyz[0]);
//   if( volb < 0.0) cout << "Warning: Barycentric coordinate is negative " << volb << endl;
    double volc = orient3d( &pa[0], &pb[0], &pd[0], &xyz[0]);
//   if( volc < 0.0) cout << "Warning: Barycentric coordinate is negative " << volc << endl;
    double vold = orient3d( &pa[0], &pc[0], &pb[0], &xyz[0]);
//   if( vold < 0.0) cout << "Warning: Barycentric coordinate is negative " << vold << endl;

    baryCoords[0] = vola/vol;
    baryCoords[1] = volb/vol;
    baryCoords[2] = volc/vol;
    baryCoords[3] = vold/vol;

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int TetGeometry :: getXYZCoordinates( const JCellPtr &cell, const Point4D &baryCoords, Point3D &xyz)
{
    assert(cell->getSize(0) == 4);
    Point3D points[4];

    points[0] = cell->getNodeAt(0)->getXYZCoords();
    points[1] = cell->getNodeAt(1)->getXYZCoords();
    points[2] = cell->getNodeAt(2)->getXYZCoords();
    points[3] = cell->getNodeAt(3)->getXYZCoords();

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    for( int i = 0; i < 4; i++) {
        xyz[0] += baryCoords[i]*points[i][0];
        xyz[1] += baryCoords[i]*points[i][1];
        xyz[2] += baryCoords[i]*points[i][2];
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

JHexahedronPtr JHexahedron::getCanonical( double len )
{
    Point3D p3d;
    JNodeSequence  nodes(8);

    p3d[0] = -0.5*len;
    p3d[1] = -0.5*len;
    p3d[2] = -0.5*len;
    nodes[0] = JNode::newObject();
    nodes[0]->setXYZCoords(p3d);

    p3d[0] =  0.5*len;
    p3d[1] = -0.5*len;
    p3d[2] = -0.5*len;
    nodes[1] = JNode::newObject();
    nodes[1]->setXYZCoords(p3d);

    p3d[0] =  0.5*len;
    p3d[1] =  0.5*len;
    p3d[2] = -0.5*len;
    nodes[2] = JNode::newObject();
    nodes[2]->setXYZCoords(p3d);

    p3d[0] = -0.5*len;
    p3d[1] =  0.5*len;
    p3d[2] = -0.5*len;
    nodes[3] = JNode::newObject();
    nodes[3]->setXYZCoords(p3d);

    p3d[0] = -0.5*len;
    p3d[1] = -0.5*len;
    p3d[2] =  0.5*len;
    nodes[4] = JNode::newObject();
    nodes[4]->setXYZCoords(p3d);

    p3d[0] =  0.5*len;
    p3d[1] = -0.5*len;
    p3d[2] =  0.5*len;
    nodes[5] = JNode::newObject();
    nodes[5]->setXYZCoords(p3d);

    p3d[0] = 0.5*len;
    p3d[1] = 0.5*len;
    p3d[2] = 0.5*len;
    nodes[6] = JNode::newObject();
    nodes[6]->setXYZCoords(p3d);

    p3d[0] = -0.5*len;
    p3d[1] =  0.5*len;
    p3d[2] =  0.5*len;
    nodes[7] = JNode::newObject();
    nodes[7]->setXYZCoords(p3d);

    JHexahedronPtr hex = JHexahedron::newObject();
    hex->setNodes(nodes);
    return hex;
}
////////////////////////////////////////////////////////////////////////////

JHexahedronPtr JHexahedron:: newObject( const JFacePtr &face1, const JFacePtr &face2,
                                       const JFacePtr &face3, const JFacePtr &face4,
                                       const JFacePtr &face5, const JFacePtr &face6)
{
    JNodeSequence hexnodes(8);

    hexnodes[0] =  face1->getNodeAt(0);
    hexnodes[1] =  face1->getNodeAt(1);
    hexnodes[2] =  face1->getNodeAt(2);
    hexnodes[3] =  face1->getNodeAt(3);

    JEdgeSequence myedges, comm_edges;
    myedges = face1->getEdges();

    JFaceSequence  myfaces(6), comm_faces;
    myfaces[0] =  face1;
    myfaces[1] =  face2;
    myfaces[2] =  face3;
    myfaces[3] =  face4;
    myfaces[4] =  face5;
    myfaces[5] =  face6;

    comm_faces.clear();
    for( int i = 0; i < 6; i++) {
        if( !myfaces[i]->hasNode( face1->getNodeAt(0) ) ) continue;
        if( !myfaces[i]->hasNode( face1->getNodeAt(1) ) ) continue;
        comm_faces.push_back( myfaces[i] );
    }
    assert( comm_faces.size() == 2 );

    JQuadrilateralPtr face0145;
    if( comm_faces[0] == myfaces[0] )  face0145 = JQuadrilateral::down_cast( comm_faces[1]);
    if( comm_faces[1] == myfaces[0] )  face0145 = JQuadrilateral::down_cast( comm_faces[0]);

    if( face0145 == NULL ) return NULL;

    JEdgePtr edge45 = JQuadrilateral::getOppositeEdge( face0145, face1->getEdgeAt(0) );
    assert( edge45 );

    comm_faces.clear();
    for( int i = 0; i < 6; i++) {
        if( !myfaces[i]->hasNode( edge45->getNodeAt(0) ) ) continue;
        if( !myfaces[i]->hasNode( edge45->getNodeAt(1) ) ) continue;
        comm_faces.push_back( myfaces[i] );
    }

    JQuadrilateralPtr face4567;
    if( comm_faces[0] == face0145 )  face4567 = JQuadrilateral::down_cast( comm_faces[1]);
    if( comm_faces[1] == face0145 )  face4567 = JQuadrilateral::down_cast( comm_faces[0]);
    if( face4567 == NULL ) return NULL;

    int pos;
    JNodePtr vtx;

    vtx = face1->getNodeAt(0);
    pos  = face0145->getPosOf(vtx);
    hexnodes[5] =  face0145->getNodeAt(pos + 2);

    vtx = face1->getNodeAt(1);
    pos  = face0145->getPosOf(vtx);
    hexnodes[4] =  face0145->getNodeAt(pos + 2);

    pos = face4567->getPosOf( hexnodes[4] );
    hexnodes[6] =  face4567->getNodeAt(pos + 2);

    pos = face4567->getPosOf( hexnodes[5] );
    hexnodes[7] =  face4567->getNodeAt(pos + 2);

    JNoImpl();
    return NULL;
}

///////////////////////////////////////////////////////////////////////////////

JHexahedronPtr JHexahedron:: newObject( const JNodePtr &v0, const  JNodePtr &v1,
                                       const JNodePtr &v2, const  JNodePtr &v3,
                                       const JNodePtr &v4, const  JNodePtr &v5,
                                       const JNodePtr &v6, const  JNodePtr &v7 )
{
    JHexahedronPtr hex = JHexahedron::newObject();
    hex->setNodes(v0,v1,v2,v3,v4,v5,v6,v7);
    return hex;
}

///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getEdgeTopology( int id, int &v0, int &v1)
{
    switch( id ) {
    case 0:
        v0 = 0;
        v1 = 1;
        break;
    case 1:
        v0 = 3;
        v1 = 2;
        break;
    case 2:
        v0 = 7;
        v1 = 6;
        break;
    case 3:
        v0 = 4;
        v1 = 5;
        break;
    // Around the Y -Axis in cw direction.
    case 4:
        v0 = 0;
        v1 = 3;
        break;
    case 5:
        v0 = 4;
        v1 = 7;
        break;
    case 6:
        v0 = 5;
        v1 = 6;
        break;
    case 7:
        v0 = 1;
        v1 = 2;
        break;
    // Around the z -Axis in cw direction.
    case 8:
        v0 = 0;
        v1 = 4;
        break;
    case 9:
        v0 = 1;
        v1 = 5;
        break;
    case 10:
        v0 = 2;
        v1 = 6;
        break;
    case 11:
        v0 = 3;
        v1 = 7;
        break;
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////
int JHexahedron :: getCommonEdgeID( int face1, int face2 )
{
    // Edge 0-1
    if( face1 == BOTTOM_SIDE && face2 == BACK_SIDE ) return 0;
    if( face2 == BOTTOM_SIDE && face1 == BACK_SIDE ) return 0;

    // Edge 3 2
    if( face1 == BACK_SIDE && face2 == TOP_SIDE ) return 1;
    if( face2 == BACK_SIDE && face1 == TOP_SIDE ) return 1;

    // Edge 7 6
    if( face1 == TOP_SIDE && face2 == FRONT_SIDE ) return 2;
    if( face2 == TOP_SIDE && face1 == FRONT_SIDE ) return 2;

    // Edge 4 5
    if( face1 == FRONT_SIDE && face2 == BOTTOM_SIDE ) return 3;
    if( face2 == FRONT_SIDE && face1 == BOTTOM_SIDE ) return 3;

    // Edge 0 3
    if( face1 == BACK_SIDE && face2 == LEFT_SIDE ) return 4;
    if( face2 == BACK_SIDE && face1 == LEFT_SIDE ) return 4;

    // Edge 4 7
    if( face1 == LEFT_SIDE && face2 == FRONT_SIDE ) return 5;
    if( face2 == LEFT_SIDE && face1 == FRONT_SIDE ) return 5;

    // Edge 5 6
    if( face1 == FRONT_SIDE && face2 == RIGHT_SIDE ) return 6;
    if( face2 == FRONT_SIDE && face1 == RIGHT_SIDE ) return 6;

    // Edge 1 2
    if( face1 == RIGHT_SIDE && face2 == BACK_SIDE ) return 7;
    if( face2 == RIGHT_SIDE && face1 == BACK_SIDE ) return 7;

    // Edge 0 4
    if( face1 == LEFT_SIDE && face2 == BOTTOM_SIDE ) return 8;
    if( face2 == LEFT_SIDE && face1 == BOTTOM_SIDE ) return 8;

    // Edge 1 5
    if( face1 == BOTTOM_SIDE && face2 == RIGHT_SIDE ) return 9;
    if( face2 == BOTTOM_SIDE && face1 == RIGHT_SIDE ) return 9;

    // Edge 2 6
    if( face1 == RIGHT_SIDE && face2 == TOP_SIDE ) return 10;
    if( face2 == RIGHT_SIDE && face1 == TOP_SIDE ) return 10;

    // Edge 3 7
    if( face1 == TOP_SIDE && face2 == LEFT_SIDE ) return 11;
    if( face2 == TOP_SIDE && face1 == LEFT_SIDE ) return 11;

    return -1;
}

int JHexahedron :: getPosOf( const JNodePtr &v0, const JNodePtr &v1) const
{
    int p[2];
    p[0] = JCell::getPosOf( v0 );
    p[1] = JCell::getPosOf( v1 );
    sort( p, p + 2 );

    // Around the X -Axis in cw direction.
    if( p[0] == 0 && p[1] == 1) return 0;
    if( p[0] == 2 && p[1] == 3) return 1;
    if( p[0] == 6 && p[1] == 7) return 2;
    if( p[0] == 4 && p[1] == 5) return 3;


    // Around the Y -Axis in cw direction.
    if( p[0] == 0 && p[1] == 3) return 4;
    if( p[0] == 4 && p[1] == 7) return 5;
    if( p[0] == 5 && p[1] == 6) return 6;
    if( p[0] == 1 && p[1] == 2) return 7;


    // Around the z -Axis in cw direction.
    if( p[0] == 0 && p[1] == 4) return 8;
    if( p[0] == 1 && p[1] == 5) return 9;
    if( p[0] == 2 && p[1] == 6) return 10;
    if( p[0] == 3 && p[1] == 7) return 11;

    cout << "Warning: No hex-edge exist " << endl;

    return -1;
}

//////////////////////////////////////////////////////////////////////////////////


int JHexahedron :: getPosOf(const JEdgePtr &edge) const
{
    assert( edge != NULL );
    for( int i= 0; i < NumEdges; i++) {
        JEdgePtr e = getEdgeAt(i);
        if( e == edge ) return i;
    }
    return -1;
}
///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getPosOf(const JFacePtr &face) const
{
    assert( face != NULL );
    for( int i= 0; i < NumFaces; i++) {
        JFacePtr f = getFaceAt(i);
        if( f == face ) return i;
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////
int JHexahedron :: getFaceTopology( int id, int &v0, int &v1, int &v2, int &v3)
{
    assert( id >= 0 && id < 6);

    switch( id ) {
    case 0:
        v0 = 0;
        v1 = 4 ;
        v2 = 7;
        v3 = 3 ;
        break;
    case 1:
        v0 = 1;
        v1 = 2;
        v2 = 6;
        v3 = 5;
        break;
    case 2:
        v0 = 0;
        v1 = 1;
        v2 = 5;
        v3 = 4;
        break;
    case 3:
        v0 = 2;
        v1 = 3;
        v2 = 7;
        v3 = 6;
        break;
    case 4:
        v0 = 0;
        v1 = 3;
        v2 = 2;
        v3 = 1;
        break;
    case 5:
        v0 = 4;
        v1 = 5;
        v2 = 6;
        v3 = 7;
        break;
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////
int JHexahedron :: getFaceID( const string &str)
{
    if( str == "Left")    return 0;
    if( str == "Right")   return 1;
    if( str == "Bottom")  return 2;
    if( str == "Top")     return 3;
    if( str == "Back")    return 4;
    if( str == "Front")   return 5;
    return -1;
}
//////////////////////////////////////////////////////////////////////////

int JHexahedron :: getPosOf( const JNodePtr &v0, const JNodePtr &v1,
                            const JNodePtr &v2, const JNodePtr &v3) const
{
    int v[4];
    v[0] = JCell::getPosOf( v0 );
    v[1] = JCell::getPosOf( v1 );
    v[2] = JCell::getPosOf( v2 );
    v[3] = JCell::getPosOf( v3 );
    sort(v, v+4);

    // Convention: Plane along the X-Y-Z directions...
    if(v[0] == 0 && v[1] ==  3 && v[2] == 4 && v[3] == 7 ) return 0;
    if(v[0] == 1 && v[1] ==  2 && v[2] == 5 && v[3] ==  6) return 1;

    if(v[0] == 0 && v[1] ==  1 && v[2] == 4 && v[3] ==  5) return 2;
    if(v[0] == 2 && v[1] ==  3 && v[2] == 6 && v[3] ==  7) return 3;

    if(v[0] == 0 && v[1] ==  1 && v[2] == 2 && v[3] ==  3) return 4;
    if(v[0] == 4 && v[1] ==  5 && v[2] == 6 && v[3] ==  7) return 5;

    return -1;
}

////////////////////////////////////////////////////////////////////////////
JEdgePtr JHexahedron :: getEdgeOf( const JNodePtr &v0, const JNodePtr &v1) const
{
    JEdgePtr nullPtr;
    if( !this->hasNode(v0) ) return nullPtr;
    if( !this->hasNode(v1) ) return nullPtr;
    return JSimplex::getEdgeOf(v0,v1,1);
}

///////////////////////////////////////////////////////////////////////////////
JFacePtr JHexahedron :: getFaceOf( const JNodePtr &v0, const JNodePtr &v1,
                                  const JNodePtr &v2, const JNodePtr &v3) const
{
    if( !this->hasNode(v0) ) return NULL;
    if( !this->hasNode(v1) ) return NULL;
    if( !this->hasNode(v2) ) return NULL;
    if( !this->hasNode(v3) ) return NULL;

    return JSimplex::getFaceOf(v0,v1,v2,v3,1);
}
///////////////////////////////////////////////////////////////////////////////
int JHexahedron :: getEdgesAt( const JNodePtr &vertex, JEdgeSequence &vertexEdges)
{
    vertexEdges.clear();

    const JNodePtr &v0 = this->getNodeAt(0);
    const JNodePtr &v1 = this->getNodeAt(1);
    const JNodePtr &v2 = this->getNodeAt(2);
    const JNodePtr &v3 = this->getNodeAt(3);
    const JNodePtr &v4 = this->getNodeAt(4);
    const JNodePtr &v5 = this->getNodeAt(5);
    const JNodePtr &v6 = this->getNodeAt(6);
    const JNodePtr &v7 = this->getNodeAt(7);

    if( vertex == v0 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v0,v1);
        vertexEdges[1] = this->getEdgeOf(v0,v3);
        vertexEdges[2] = this->getEdgeOf(v0,v4);
        return 0;
    }

    if( vertex == v1 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v1,v0);
        vertexEdges[1] = this->getEdgeOf(v1,v2);
        vertexEdges[2] = this->getEdgeOf(v1,v5);
        return 0;
    }

    if( vertex == v3 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v2,v3);
        vertexEdges[1] = this->getEdgeOf(v2,v1);
        vertexEdges[2] = this->getEdgeOf(v2,v6);
        return 0;
    }

    if( vertex == v3 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v3,v2);
        vertexEdges[1] = this->getEdgeOf(v3,v0);
        vertexEdges[2] = this->getEdgeOf(v3,v7);
        return 0;
    }

    if( vertex == v4 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v4,v5);
        vertexEdges[1] = this->getEdgeOf(v4,v7);
        vertexEdges[2] = this->getEdgeOf(v4,v0);
        return 0;
    }

    if( vertex == v5 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v5,v4);
        vertexEdges[1] = this->getEdgeOf(v5,v6);
        vertexEdges[2] = this->getEdgeOf(v5,v1);
        return 0;
    }

    if( vertex == v6 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v6,v7);
        vertexEdges[1] = this->getEdgeOf(v6,v5);
        vertexEdges[2] = this->getEdgeOf(v6,v2);
        return 0;
    }

    if( vertex == v7 ) {
        vertexEdges.resize(3);
        vertexEdges[0] = this->getEdgeOf(v7,v6);
        vertexEdges[1] = this->getEdgeOf(v7,v4);
        vertexEdges[2] = this->getEdgeOf(v7,v3);
        return 0;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getFacesAt( const JNodePtr &vertex, JFaceSequence &vertexFaces)
{
    vertexFaces.clear();

    const JNodePtr &v0 = this->getNodeAt(0);
    const JNodePtr &v1 = this->getNodeAt(1);
    const JNodePtr &v2 = this->getNodeAt(2);
    const JNodePtr &v3 = this->getNodeAt(3);
    const JNodePtr &v4 = this->getNodeAt(4);
    const JNodePtr &v5 = this->getNodeAt(5);
    const JNodePtr &v6 = this->getNodeAt(6);
    const JNodePtr &v7 = this->getNodeAt(7);

    if( vertex == v0 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v0, v4, v7, v3);
        vertexFaces[1] = this->getFaceOf(v0, v1, v5, v4);
        vertexFaces[2] = this->getFaceOf(v0, v1, v2, v3);
        return 0;
    }

    if( vertex == v1 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v1, v2, v6, v5);
        vertexFaces[1] = this->getFaceOf(v1, v5, v4, v0);
        vertexFaces[2] = this->getFaceOf(v1, v2, v3, v0);
        return 0;
    }

    if( vertex == v2 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v2, v6, v5, v1);
        vertexFaces[1] = this->getFaceOf(v2, v3, v7, v6);
        vertexFaces[2] = this->getFaceOf(v2, v1, v0, v3);
        return 0;
    }

    if( vertex == v3 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v3, v0, v4, v7);
        vertexFaces[1] = this->getFaceOf(v3, v7, v6, v2);
        vertexFaces[2] = this->getFaceOf(v3, v2, v1, v0);
        return 0;
    }

    if( vertex == v4 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v4, v7, v3, v0);
        vertexFaces[1] = this->getFaceOf(v4, v0, v1, v5);
        vertexFaces[2] = this->getFaceOf(v4, v5, v6, v7);
        return 0;
    }

    if( vertex == v5 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v5, v1, v2, v6);
        vertexFaces[1] = this->getFaceOf(v5, v4, v0, v1);
        vertexFaces[2] = this->getFaceOf(v5, v6, v7, v4);
        return 0;
    }

    if( vertex == v6 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v6, v5, v1, v2);
        vertexFaces[1] = this->getFaceOf(v6, v2, v3, v7);
        vertexFaces[2] = this->getFaceOf(v6, v7, v4, v5);
        return 0;
    }
    if( vertex == v7 ) {
        vertexFaces.resize(3);
        vertexFaces[0] = this->getFaceOf(v7, v3, v0, v4);
        vertexFaces[1] = this->getFaceOf(v7, v6, v2, v3);
        vertexFaces[2] = this->getFaceOf(v7, v4, v5, v6);
        return 0;
    }
    return 1;
}

////////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getFacesAt( const JEdgePtr &edge, JFaceSequence &edgeFaces)
{
    edgeFaces.clear();

    int pos = this->getPosOf( edge );

    JNodePtr v0 = this->getNodeAt(0);
    JNodePtr v1 = this->getNodeAt(1);
    JNodePtr v2 = this->getNodeAt(2);
    JNodePtr v3 = this->getNodeAt(3);
    JNodePtr v4 = this->getNodeAt(4);
    JNodePtr v5 = this->getNodeAt(5);
    JNodePtr v6 = this->getNodeAt(6);
    JNodePtr v7 = this->getNodeAt(7);

    switch( pos ) {
    case 0:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v0, v1, v5, v4);
        edgeFaces[1] = this->getFaceOf(v0, v1, v2, v3);
        break;
    case 1:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v2, v3, v7, v6);
        edgeFaces[1] = this->getFaceOf(v2, v3, v0, v1);
        break;
    case 2:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v6, v7, v3, v2);
        edgeFaces[1] = this->getFaceOf(v6, v7, v4, v5);
        break;
    case 3:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v4, v5, v1, v0);
        edgeFaces[1] = this->getFaceOf(v4, v5, v6, v7);
        break;
    case 4:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v0, v3, v7, v4);
        edgeFaces[1] = this->getFaceOf(v0, v3, v2, v1);
        break;
    case 5:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v4, v7, v3, v0);
        edgeFaces[1] = this->getFaceOf(v4, v7, v6, v5);
        break;
    case 6:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v5, v6, v2, v1);
        edgeFaces[1] = this->getFaceOf(v5, v6, v7, v4);
        break;
    case 7:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v1, v2, v6, v5);
        edgeFaces[1] = this->getFaceOf(v1, v2, v3, v0);
        break;
    case 8:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v0, v4, v7, v3);
        edgeFaces[1] = this->getFaceOf(v0, v4, v5, v1);
        break;
    case 9:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v1, v5, v6, v2);
        edgeFaces[1] = this->getFaceOf(v1, v5, v4, v0);
        break;
    case 10:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v2, v6, v5, v1);
        edgeFaces[1] = this->getFaceOf(v2, v6, v7, v3);
        break;
    case 11:
        edgeFaces.resize(2);
        edgeFaces[0] = this->getFaceOf(v3, v7, v4, v0);
        edgeFaces[1] = this->getFaceOf(v3, v7, v6, v2);
        break;
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: build_lower_entities( int dim )
{
    switch(dim) {
    case 1:
        for( int i = 0; i < JTetrahedron::NumEdges; i++)
            this->getEdgeAt(i);
        break;
    case 2:
        for( int i = 0; i < JTetrahedron::NumFaces; i++)
            this->getFaceAt(i);
        break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JHexahedron :: getOrientation(const JEdgePtr &edge) const
{
    int i0, i1;
    for( int i = 0; i < NumEdges; i++) {
        if( getEdgeAt(i) == edge ) {
            getEdgeTopology(i, i0, i1);
            JNodePtr v0 = getNodeAt(i0);
            JNodePtr v1 = getNodeAt(i1);
            return edge->getOrientation(v0,v1);
        }
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
int JHexahedron :: getOrientation( const JFacePtr &face) const
{
    if( face->getSize(0) != 4) return 0;

    int i0, i1, i2, i3;
    for( int i = 0; i < NumFaces; i++) {
        if( getFaceAt(i) == face ) {
            getFaceTopology(i, i0, i1, i2, i3);
            JNodePtr v0 = getNodeAt(i0);
            JNodePtr v1 = getNodeAt(i1);
            JNodePtr v2 = getNodeAt(i2);
            JNodePtr v3 = getNodeAt(i3);
            return face->getOrientation(v0,v1,v2, v3);
        }
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getOppositeFaceNodes( const JNodeSequence &this_facenodes,
                                        JNodeSequence &opp_facenodes) const
{
    opp_facenodes.resize(4);

    int v[4];
    v[0] = JCell::getPosOf( this_facenodes[0] );
    v[1] = JCell::getPosOf( this_facenodes[1] );
    v[2] = JCell::getPosOf( this_facenodes[2] );
    v[3] = JCell::getPosOf( this_facenodes[3] );

    std::sort(  v, v + 4 );

    // X -Axis faces ...
    if(v[0] == 0 && v[1] ==  3 && v[2] == 4 && v[3] == 7 )  {
        opp_facenodes[0] = getNodeAt(1);
        opp_facenodes[1] = getNodeAt(2);
        opp_facenodes[2] = getNodeAt(6);
        opp_facenodes[3] = getNodeAt(5);
        return 0;
    }

    if(v[0] == 1 && v[1] ==  2 && v[2] == 5 && v[3] ==  6)  {
        opp_facenodes[0] = getNodeAt(0);
        opp_facenodes[1] = getNodeAt(4);
        opp_facenodes[2] = getNodeAt(7);
        opp_facenodes[3] = getNodeAt(3);
        return 0;
    }

    // Z Axis faces ...
    if(v[0] == 0 && v[1] ==  1 && v[2] == 2 && v[3] ==  3)  {
        opp_facenodes[0] = getNodeAt(4);
        opp_facenodes[1] = getNodeAt(5);
        opp_facenodes[2] = getNodeAt(6);
        opp_facenodes[3] = getNodeAt(7);
        return 0;
    }

    if(v[0] == 4 && v[1] ==  5 && v[2] == 6 && v[3] ==  7)  {
        opp_facenodes[0] = getNodeAt(0);
        opp_facenodes[1] = getNodeAt(3);
        opp_facenodes[2] = getNodeAt(2);
        opp_facenodes[3] = getNodeAt(1);
        return 0;
    }

    // Y Axis faces ..

    if(v[0] == 0 && v[1] ==  1 && v[2] == 4 && v[3] ==  5)  {
        opp_facenodes[0] = getNodeAt(2);
        opp_facenodes[1] = getNodeAt(3);
        opp_facenodes[2] = getNodeAt(7);
        opp_facenodes[3] = getNodeAt(6);
        return 0;
    }


    if(v[0] == 2 && v[1] ==  3 && v[2] == 6 && v[3] ==  7)  {
        opp_facenodes[0] = getNodeAt(0);
        opp_facenodes[1] = getNodeAt(1);
        opp_facenodes[2] = getNodeAt(5);
        opp_facenodes[3] = getNodeAt(4);
        return 0;
    }

    cout << "Error: Invalid face " << endl;
    exit(0);

    return 1;
}
///////////////////////////////////////////////////////////////////

JFacePtr JHexahedron :: getOppositeFace( const JFacePtr &f) const
{
    JFacePtr parface;

    int v[4];

    v[0] = JCell::getPosOf( f->getNodeAt(0) );
    v[1] = JCell::getPosOf( f->getNodeAt(1) );
    v[2] = JCell::getPosOf( f->getNodeAt(2) );
    v[3] = JCell::getPosOf( f->getNodeAt(3) );
    std::sort(  v, v + 4);

    JNodeSequence facenodes(4);
    JNodePtr vhash;

    // X Axis...
    if(v[0] == 0 && v[1] ==  3 && v[2] == 4 && v[3] == 7 )  {
        facenodes[0] = getNodeAt(1);
        facenodes[1] = getNodeAt(2);
        facenodes[2] = getNodeAt(6);
        facenodes[3] = getNodeAt(5);
        vhash = *min_element(facenodes );
        parface = vhash->getFaceOf( facenodes );

#ifdef DEBUG
        v[0] = Cell::getPosOf( parface->getNodeAt(0) );
        v[1] = Cell::getPosOf( parface->getNodeAt(1) );
        v[2] = Cell::getPosOf( parface->getNodeAt(2) );
        v[3] = Cell::getPosOf( parface->getNodeAt(3) );
        std::sort(  v, v + 4);
        assert( v[0] == 1 && v[1] == 2 && v[2] == 5 && v[3] == 6);
#endif
        assert(parface);
        return parface;
    }

    if(v[0] == 1 && v[1] ==  2 && v[2] == 5 && v[3] ==  6)  {
        facenodes[0] = getNodeAt(0);
        facenodes[1] = getNodeAt(4);
        facenodes[2] = getNodeAt(7);
        facenodes[3] = getNodeAt(3);
        vhash = *min_element(facenodes );
        parface = vhash->getFaceOf( facenodes );

#ifdef DEBUG
        v[0] = Cell::getPosOf( parface->getNodeAt(0) );
        v[1] = Cell::getPosOf( parface->getNodeAt(1) );
        v[2] = Cell::getPosOf( parface->getNodeAt(2) );
        v[3] = Cell::getPosOf( parface->getNodeAt(3) );
        std::sort(  v, v + 4);
        assert( v[0] == 0 && v[1] == 3 && v[2] == 4 && v[3] == 7);
#endif

        assert(parface);
        return parface;
    }

    // Z Axis ...

    if(v[0] == 0 && v[1] ==  1 && v[2] == 2 && v[3] ==  3)  {
        facenodes[0] = getNodeAt(4);
        facenodes[1] = getNodeAt(5);
        facenodes[2] = getNodeAt(6);
        facenodes[3] = getNodeAt(7);
        vhash = *min_element(facenodes );
        parface = vhash->getFaceOf( facenodes );
#ifdef DEBUG
        v[0] = Cell::getPosOf( parface->getNodeAt(0) );
        v[1] = Cell::getPosOf( parface->getNodeAt(1) );
        v[2] = Cell::getPosOf( parface->getNodeAt(2) );
        v[3] = Cell::getPosOf( parface->getNodeAt(3) );
        std::sort(  v, v + 4);
        assert( v[0] == 4 && v[1] == 5 && v[2] == 6 && v[3] == 7);
#endif
        assert(parface);
        return parface;
    }

    if(v[0] == 4 && v[1] ==  5 && v[2] == 6 && v[3] ==  7)  {
        facenodes[0] = getNodeAt(0);
        facenodes[1] = getNodeAt(3);
        facenodes[2] = getNodeAt(2);
        facenodes[3] = getNodeAt(1);
        vhash = *min_element( facenodes );
        parface = vhash->getFaceOf( facenodes );
#ifdef DEBUG
        v[0] = Cell::getPosOf( parface->getNodeAt(0) );
        v[1] = Cell::getPosOf( parface->getNodeAt(1) );
        v[2] = Cell::getPosOf( parface->getNodeAt(2) );
        v[3] = Cell::getPosOf( parface->getNodeAt(3) );
        std::sort(  v, v + 4);
        assert( v[0] == 0 && v[1] == 1 && v[2] == 2 && v[3] == 3);
#endif

        assert(parface);
        return parface;
    }

    // Y Axis ...

    if(v[0] == 0 && v[1] ==  1 && v[2] == 4 && v[3] ==  5)  {
        facenodes[0] = getNodeAt(2);
        facenodes[1] = getNodeAt(3);
        facenodes[2] = getNodeAt(7);
        facenodes[3] = getNodeAt(6);
        vhash = *min_element( facenodes );
        parface = vhash->getFaceOf( facenodes );
#ifdef DEBUG
        v[0] = Cell::getPosOf( parface->getNodeAt(0) );
        v[1] = Cell::getPosOf( parface->getNodeAt(1) );
        v[2] = Cell::getPosOf( parface->getNodeAt(2) );
        v[3] = Cell::getPosOf( parface->getNodeAt(3) );
        std::sort(  v, v + 4);
        assert( v[0] == 2 && v[1] == 3 && v[2] == 6 && v[3] == 7);
#endif
        assert(parface);
        return parface;
    }

    if(v[0] == 2 && v[1] ==  3 && v[2] == 6 && v[3] ==  7)  {
        facenodes[0] = getNodeAt(0);
        facenodes[1] = getNodeAt(1);
        facenodes[2] = getNodeAt(5);
        facenodes[3] = getNodeAt(4);
        vhash = *min_element(facenodes );
        parface = vhash->getFaceOf( facenodes );
#ifdef DEBUG
        v[0] = Cell::getPosOf( parface->getNodeAt(0) );
        v[1] = Cell::getPosOf( parface->getNodeAt(1) );
        v[2] = Cell::getPosOf( parface->getNodeAt(2) );
        v[3] = Cell::getPosOf( parface->getNodeAt(3) );
        std::sort(  v, v + 4);
        assert( v[0] == 0 && v[1] == 1 && v[2] == 4 && v[3] == 5);
#endif
        assert(parface);
        return parface;
    }

    cout << "Error: Invalid face " << parface << endl;

    return parface;
}
///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getNeighbors( const JEdgePtr &edge, JFaceSequence &efaces)
{
    efaces.resize(2);

    JNodePtr v0 = edge->getNodeAt(0);
    JNodePtr v1 = edge->getNodeAt(1);
    JNodePtr v2, v3, v4, v5;

    int pos0   = JCell::getPosOf( v0);
    int pos1   = JCell::getPosOf( v1);
    int p0     = std::min( pos0, pos1);
    int p1     = std::max( pos0, pos1);

    if( p0 == 0 && p1 == 1 ) {
        v2 = getNodeAt(2);
        v3 = getNodeAt(3);
        v4 = getNodeAt(4);
        v5 = getNodeAt(5);
    }

    if( p0 == 0 && p1 == 3 ) {
        v2 = getNodeAt(1);
        v3 = getNodeAt(2);
        v4 = getNodeAt(4);
        v5 = getNodeAt(7);
    }

    if( p0 == 0 && p1 == 4 ) {
        v2 = getNodeAt(3);
        v3 = getNodeAt(7);
        v4 = getNodeAt(1);
        v5 = getNodeAt(5);
    }

    if( p0 == 1 && p1 == 2 ) {
        v2 = getNodeAt(0);
        v3 = getNodeAt(3);
        v4 = getNodeAt(5);
        v5 = getNodeAt(6);
    }

    if( p0 == 1 && p1 == 5 ) {
        v2 = getNodeAt(0);
        v3 = getNodeAt(4);
        v4 = getNodeAt(2);
        v5 = getNodeAt(6);
    }

    if( p0 == 2 && p1 == 3 ) {
        v2 = getNodeAt(0);
        v3 = getNodeAt(1);
        v4 = getNodeAt(6);
        v5 = getNodeAt(7);
    }

    if( p0 == 2 && p1 == 6 ) {
        v2 = getNodeAt(1);
        v3 = getNodeAt(5);
        v4 = getNodeAt(3);
        v5 = getNodeAt(7);
    }

    if( p0 == 3 && p1 == 7 ) {
        v2 = getNodeAt(0);
        v3 = getNodeAt(4);
        v4 = getNodeAt(2);
        v5 = getNodeAt(6);
    }

    if( p0 == 4 && p1 == 5 ) {
        v2 = getNodeAt(0);
        v3 = getNodeAt(1);
        v4 = getNodeAt(6);
        v5 = getNodeAt(7);
    }

    if( p0 == 4 && p1 == 7 ) {
        v2 = getNodeAt(0);
        v3 = getNodeAt(3);
        v4 = getNodeAt(5);
        v5 = getNodeAt(6);
    }

    if( p0 == 5 && p1 == 6 ) {
        v2 = getNodeAt(1);
        v3 = getNodeAt(2);
        v4 = getNodeAt(4);
        v5 = getNodeAt(7);
    }

    if( p0 == 6 && p1 == 7 ) {
        v2 = getNodeAt(2);
        v3 = getNodeAt(3);
        v4 = getNodeAt(4);
        v5 = getNodeAt(5);
    }

    efaces[0] = this->getFaceOf(v0,v1,v2,v3);
    efaces[1] = this->getFaceOf(v0,v1,v4,v5);
    assert( efaces[0] );
    assert( efaces[1] );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JCell::remove_unattached_lower_entities()
{
    if( this->getStatus() != JMeshEntity::REMOVE) return 1;

    int nfaces = this->getSize(2);
    for( int i = 0; i < nfaces; i++) {
        JFacePtr face = this->getFaceAt(i);
        if( face->getNumHigherRelations() == 0)  {
            face->setStatus( JMeshEntity::REMOVE );
            face->remove_unattached_lower_entities();
        }
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JCell :: get_shared_entities( const JCellPtr &cell1, const JCellPtr &cell2, JNodeSequence &nodes)
{
    nodes.clear();
    int nsize = cell1->getSize(0);
    for( int i = 0; i < nsize; i++) {
        JNodePtr vtx = cell1->getNodeAt(i);
        if( cell2->hasNode( vtx ) ) nodes.push_back(vtx);
    }
    return nodes.size();
}

///////////////////////////////////////////////////////////////////////////////

int JCell :: get_shared_entities( const JCellPtr &cell1, const JCellPtr &cell2, JEdgeSequence &comm_edges)
{
    comm_edges.clear();

    JEdgeSequence cell1edges = cell1->getEdges();
    JEdgeSequence cell2edges = cell2->getEdges();

    boost::sort( cell1edges );
    boost::sort( cell2edges );

    boost::set_intersection( cell1edges, cell2edges, back_inserter(comm_edges) );

    return comm_edges.size();
}

///////////////////////////////////////////////////////////////////////////////

int JCell :: get_shared_entities( const JCellPtr &cell1, const JCellPtr &cell2, JFaceSequence &comm_faces)
{
    comm_faces.clear();

    JFaceSequence cell1faces = cell1->getFaces();
    JFaceSequence cell2faces = cell2->getFaces();

    boost::sort( cell1faces );
    boost::sort( cell2faces );

    boost::set_intersection( cell1faces, cell2faces, back_inserter(comm_faces) );

    return comm_faces.size();
}

///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: rotate_common_face( const JMeshPtr &mesh, const JHexahedronPtr &hex1, const JHexahedronPtr &hex2, int dir)
{
    assert( mesh->getAdjTable(2,3) );

    JFaceSequence comm_faces;
    get_shared_entities( hex1, hex2, comm_faces);
    if( comm_faces.size() != 1) return 1;

    const JNodePtr &v0 = comm_faces[0]->getNodeAt(0);
    const JNodePtr &v1 = comm_faces[0]->getNodeAt(1);
    const JNodePtr &v2 = comm_faces[0]->getNodeAt(2);
    const JNodePtr &v3 = comm_faces[0]->getNodeAt(3);

    const JNodePtr &v4 = hex1->getDiagonalNode(v2);
    const JNodePtr &v5 = hex1->getDiagonalNode(v3);
    const JNodePtr &v6 = hex1->getDiagonalNode(v0);
    const JNodePtr &v7 = hex1->getDiagonalNode(v1);

    const JNodePtr &v8  = hex2->getDiagonalNode(v2);
    const JNodePtr &v9  = hex2->getDiagonalNode(v3);
    const JNodePtr &v10 = hex2->getDiagonalNode(v0);
    const JNodePtr &v11 = hex2->getDiagonalNode(v1);

    /*
         if( dir == 1) {
              hex1->setNodes( v6, v9, v10, v2, v7, v8, v11, v3 );
              hex2->setNodes( v6, v5, v1, v9, v7, v4, v0, v8 );
         }

         if( dir == -1) {
              hex1->setNodes( v6, v2, v10, v5, v7, v3, v11, v4 );
              hex2->setNodes( v5, v10, v9, v1, v4, v11, v8, v0 );
         }
    */

    if( dir == 1) {
        hex1->setNodes( v4, v5, v1, v9, v7,  v6, v2,  v10 );
        hex2->setNodes( v4, v9, v8, v0, v7, v10, v11, v3 );
    }

    if( dir == -1) {
        hex1->setNodes( v4, v5, v8, v0, v7, v6, v11, v3  );
        hex2->setNodes( v5, v1, v9, v8, v6, v2, v10, v11 );
    }

    return 0;

}

////////////////////////////////////////////////////////////////////////////////

JEdgeSequence JHexahedron :: getEdges()  const
{
    JEdgeSequence edges(NumEdges);

    int p0, p1;
    for( int i = 0; i < NumEdges; i++) {
        getEdgeTopology(i, p0, p1);
        JNodePtr v0 = nodes[p0];
        JNodePtr v1 = nodes[p1];
        edges[i] = JSimplex::getEdgeOf(v0,v1,1);
        assert( edges[i] );
    }
    return edges;
}

////////////////////////////////////////////////////////////////////////////////////
JFaceSequence JHexahedron :: getFaces() const
{
    JFaceSequence faces(6);

    int p0, p1, p2, p3;
    for( int i = 0; i < 6; i++) {
        getFaceTopology( i, p0, p1, p2, p3);
        JNodePtr v0 = getNodeAt(p0);
        JNodePtr v1 = getNodeAt(p1);
        JNodePtr v2 = getNodeAt(p2);
        JNodePtr v3 = getNodeAt(p3);
        faces[i]  = JSimplex::getFaceOf(v0,v1,v2,v3,1);
        assert( faces[i] );
    }
    return faces;
}
///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getTetrahedra( vector<JCellPtr> &tets) const
{
    tets.resize(6);
    JNodeSequence tnodes(4);

    tnodes[0] = getNodeAt(0);
    tnodes[1] = getNodeAt(4);
    tnodes[2] = getNodeAt(5);
    tnodes[3] = getNodeAt(7);
    tets[0] = JTetrahedron::newObject(tnodes);

    tnodes[0] = getNodeAt(0);
    tnodes[1] = getNodeAt(5);
    tnodes[2] = getNodeAt(6);
    tnodes[3] = getNodeAt(7);
    tets[1] = JTetrahedron::newObject(tnodes);

    tnodes[0] = getNodeAt(0);
    tnodes[1] = getNodeAt(6);
    tnodes[2] = getNodeAt(3);
    tnodes[3] = getNodeAt(7);
    tets[2] = JTetrahedron::newObject(tnodes);

    tnodes[0] = getNodeAt(0);
    tnodes[1] = getNodeAt(1);
    tnodes[2] = getNodeAt(6);
    tnodes[3] = getNodeAt(5);
    tets[3] = JTetrahedron::newObject(tnodes);

    tnodes[0] = getNodeAt(0);
    tnodes[1] = getNodeAt(6);
    tnodes[2] = getNodeAt(1);
    tnodes[3] = getNodeAt(3);
    tets[4] = JTetrahedron::newObject(tnodes);

    tnodes[0] = getNodeAt(2);
    tnodes[1] = getNodeAt(1);
    tnodes[2] = getNodeAt(6);
    tnodes[3] = getNodeAt(3);
    tets[5] = JTetrahedron::newObject(tnodes);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: get_topological_parallel_edges(const JEdgePtr &edge, JEdgeSequence &paredges) const
{
    paredges.clear();

    JNodePtr v0 = edge->getNodeAt(0);
    JNodePtr v1 = edge->getNodeAt(1);
    assert(v0 != nullptr  && v1 != nullptr);

    int p0  = min( JCell::getPosOf(v0), JCell::getPosOf(v1) );
    int p1  = max( JCell::getPosOf(v0), JCell::getPosOf(v1) );


    // Give the edges in clockwise direction ....
    assert( p0 >= 0 && p0 < 8);
    assert( p1 >= 0 && p1 < 8);

    paredges.resize(3);

    if( p0 == 0 && p1 == 1 ) {
        v0 = getNodeAt(4);
        v1 = getNodeAt(5);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(6);
        v1 = getNodeAt(7);
        paredges[1] = getEdgeOf(v0,v1);

        v0 = getNodeAt(2);
        v1 = getNodeAt(3);
        paredges[2] = getEdgeOf(v0,v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );

        return 0;
    }

    if( p0 == 4 && p1 == 5 ) {
        v0 = getNodeAt(6);
        v1 = getNodeAt(7);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(2);
        v1 = getNodeAt(3);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(0);
        v1 = getNodeAt(1);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if(  p0 == 6 && p1 == 7 ) {
        v0 = getNodeAt(2);
        v1 = getNodeAt(3);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(0);
        v1 = getNodeAt(1);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(4);
        v1 = getNodeAt(5);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if(  p0 == 2 && p1 == 3 ) {
        v0 = getNodeAt(0);
        v1 = getNodeAt(1);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(4);
        v1 = getNodeAt(5);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(6);
        v1 = getNodeAt(7);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if( p0 == 0 && p1 == 3 ) {
        v0 = getNodeAt(1);
        v1 = getNodeAt(2);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(5);
        v1 = getNodeAt(6);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(4);
        v1 = getNodeAt(7);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if( p0 == 1 && p1 == 2) {
        v0 = getNodeAt(5);
        v1 = getNodeAt(6);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(4);
        v1 = getNodeAt(7);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(0);
        v1 = getNodeAt(3);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if(  p0 == 5 && p1 == 6 ) {
        v0 = getNodeAt(4);
        v1 = getNodeAt(7);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(0);
        v1 = getNodeAt(3);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(1);
        v1 = getNodeAt(2);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if( p0 == 4 && p1 == 7 ) {
        v0 = getNodeAt(0);
        v1 = getNodeAt(3);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(1);
        v1 = getNodeAt(2);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(5);
        v1 = getNodeAt(6);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if( p0 == 0 && p1 == 4 ) {
        v0 = getNodeAt(3);
        v1 = getNodeAt(7);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(2);
        v1 = getNodeAt(6);
        paredges[1] = getEdgeOf(v0, v1);

        v0 = getNodeAt(1);
        v1 = getNodeAt(5);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if( p0 == 3 && p1 == 7 ) {
        v0 = getNodeAt(2);
        v1 = getNodeAt(6);
        paredges[0] = getEdgeOf(v0, v1);

        v0 = getNodeAt(1);
        v1 = getNodeAt(5);
        paredges[1] = getEdgeOf(v0,v1);

        v0 = getNodeAt(0);
        v1 = getNodeAt(4);
        paredges[2] = getEdgeOf(v0, v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if( p0 == 2 && p1 == 6 ) {
        v0 = getNodeAt(1);
        v1 = getNodeAt(5);
        paredges[0] = getEdgeOf(v0,v1);

        v0 = getNodeAt(0);
        v1 = getNodeAt(4);
        paredges[1] = getEdgeOf(v0,v1);

        v0 = getNodeAt(3);
        v1 = getNodeAt(7);
        paredges[2] = getEdgeOf(v0,v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }

    if( p0 == 1 && p1 == 5 ) {
        v0 = getNodeAt(0);
        v1 = getNodeAt(4);
        paredges[0] = getEdgeOf(v0,v1);

        v0 = getNodeAt(3);
        v1 = getNodeAt(7);
        paredges[1]= getEdgeOf(v0,v1);

        v0 = getNodeAt(2);
        v1 = getNodeAt(6);
        paredges[2] = getEdgeOf(v0,v1);

        assert( paredges[0] != NULL );
        assert( paredges[1] != NULL );
        assert( paredges[2] != NULL );
        return 0;
    }
    cout << p0 << " " << p1 << endl;

    cout << "Error: No hexahedron edge " << endl;
    exit(0);
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getSignature( const JEdgePtr &edge) const
{
    int p[2];
    p[0] = JCell::getPosOf( edge->getNodeAt(0));
    p[1] = JCell::getPosOf( edge->getNodeAt(1));
    sort( p, p+2);

    if( p[0] == 0 && p[1] == 1 ) return 10;
    if( p[0] == 2 && p[1] == 3 ) return 23;
    if( p[0] == 6 && p[1] == 7 ) return 67;
    if( p[0] == 4 && p[1] == 5 ) return 45;

    if( p[0] == 0 && p[1] == 3 ) return 30;
    if( p[0] == 1 && p[1] == 2 ) return 12;
    if( p[0] == 5 && p[1] == 6 ) return 56;
    if( p[0] == 4 && p[1] == 7 ) return 47;

    if( p[0] == 0 && p[1] == 4 ) return 40;
    if( p[0] == 3 && p[1] == 7 ) return 37;
    if( p[0] == 2 && p[1] == 6 ) return 26;
    if( p[0] == 1 && p[1] == 5 ) return 15;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JHexahedron :: getSignature( const JFacePtr &face) const
{
    int p[4];
    p[0] = JCell::getPosOf( face->getNodeAt(0));
    p[1] = JCell::getPosOf( face->getNodeAt(1));
    p[2] = JCell::getPosOf( face->getNodeAt(2));
    p[3] = JCell::getPosOf( face->getNodeAt(3));
    sort( p, p+4);

    if( p[0] == 0 && p[1] == 1  && p[2] == 2 && p[3] == 3) return 1230;
    if( p[0] == 1 && p[1] == 2  && p[2] == 5 && p[3] == 6) return 1256;
    if( p[0] == 0 && p[1] == 1  && p[2] == 4 && p[3] == 5) return 1450;
    if( p[0] == 2 && p[1] == 3  && p[2] == 6 && p[3] == 7) return 2367;
    if( p[0] == 0 && p[1] == 3  && p[2] == 4 && p[3] == 7) return 3470;
    if( p[0] == 4 && p[1] == 5  && p[2] == 6 && p[3] == 6) return 4567;
    return 0;
}

/////////////////////////////////////////////////////////////////////////////

int JHexahedron :: getCyclicFaces( const JEdgePtr &edge, JFaceSequence &cfaces) const
{
    cfaces.clear();
    JNodePtr v0 = getNodeAt(0);
    JNodePtr v1 = getNodeAt(1);
    JNodePtr v2 = getNodeAt(2);
    JNodePtr v3 = getNodeAt(3);
    JNodePtr v4 = getNodeAt(4);
    JNodePtr v5 = getNodeAt(5);
    JNodePtr v6 = getNodeAt(6);
    JNodePtr v7 = getNodeAt(7);

    JFacePtr f1230 = getFaceOf(v0, v1, v2, v3);
    assert(f1230);
    JFacePtr f4567 = getFaceOf(v4, v5, v6, v7);
    assert(f4567);

    JFacePtr f3470 = getFaceOf(v0, v3, v4, v7);
    assert(f3470);
    JFacePtr f1256 = getFaceOf(v1, v2, v5, v6);
    assert(f1256);

    JFacePtr f1450 = getFaceOf(v0, v1, v4, v5);
    assert(f1450);
    JFacePtr f2367 = getFaceOf(v2, v3, v6, v7);
    assert(f2367);

    int sig = getSignature(edge);

    // Revolve around X-Axis and get all the faces ...
    if( sig == 10 || sig == 23 || sig == 45 || sig == 67 ) {
        cfaces.resize(4);
        cfaces[0] =  f1450;
        cfaces[1] =  f4567;
        cfaces[2] =  f2367;
        cfaces[3] =  f1230;
        return 0;
    }

    // Revolve around Y Axis and get all the faces ...
    if( sig == 12 || sig == 30 || sig == 47 || sig == 56 ) {
        cfaces.resize(4);
        cfaces[0] =  f1230;
        cfaces[1] =  f1256;
        cfaces[2] =  f4567;
        cfaces[3] =  f3470;
        return 0;
    }

    if( sig == 15 || sig == 26 || sig == 37 || sig == 40 ) {
        cfaces.resize(4);
        cfaces[0] =  f3470;
        cfaces[1] =  f2367;
        cfaces[2] =  f1256;
        cfaces[3] =  f1450;
        return 0;
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

inline
JEdgePtr JHexahedron :: getEdgeAt(int id) const
{
    assert( id >=0 && id < 12);
    int p0, p1;
    getEdgeTopology(id, p0, p1);
    const JNodePtr &v0 = getNodeAt(p0);
    const JNodePtr &v1 = getNodeAt(p1);
    return JSimplex::getEdgeOf( v0, v1, 1);
}
////////////////////////////////////////////////////////////////////////////////

inline
JEdgePtr JHexahedron :: getEdgeAt(int id, int &ori) const
{
    assert( id >=0 && id < 12);
    int p0, p1;
    getEdgeTopology(id, p0, p1);
    const JNodePtr &v0 = getNodeAt(p0);
    const JNodePtr &v1 = getNodeAt(p1);
    JEdgePtr edge = JSimplex::getEdgeOf( v0, v1, 1);
    ori = edge->getOrientation( v0, v1);
    return edge;
}

////////////////////////////////////////////////////////////////////////////////

inline
bool JHexahedron :: hasFaceAt(int id) const
{
    assert( id >=0 && id < 6);

    int p0, p1, p2, p3;
    getFaceTopology(id, p0, p1, p2, p3);

    const JNodePtr &v0 = getNodeAt(p0);
    const JNodePtr &v1 = getNodeAt(p1);
    const JNodePtr &v2 = getNodeAt(p2);
    const JNodePtr &v3 = getNodeAt(p3);

    JFacePtr f = JSimplex::getFaceOf( v0, v1, v2, v3, 0);
    if( f ) return 1;
    return 0;
}

inline
JFacePtr JHexahedron :: getFaceAt(int id) const
{
    assert( id >=0 && id < 6);

    int p0, p1, p2, p3;
    getFaceTopology(id, p0, p1, p2, p3);

    const JNodePtr &v0 = getNodeAt(p0);
    const JNodePtr &v1 = getNodeAt(p1);
    const JNodePtr &v2 = getNodeAt(p2);
    const JNodePtr &v3 = getNodeAt(p3);

    return JSimplex::getFaceOf( v0, v1, v2, v3, 1);
}
////////////////////////////////////////////////////////////////////////////////

inline
JFacePtr JHexahedron :: getFaceAt(int id, int &ori) const
{
    assert( id >=0 && id < 6);

    int p0, p1, p2, p3;
    getFaceTopology(id, p0, p1, p2, p3);

    const JNodePtr &v0 = getNodeAt(p0);
    const JNodePtr &v1 = getNodeAt(p1);
    const JNodePtr &v2 = getNodeAt(p2);
    const JNodePtr &v3 = getNodeAt(p3);

    JFacePtr fs = JSimplex::getFaceOf( v0, v1, v2, v3, 1);
    ori  = fs->getOrientation(v0,v1,v2,v3);
    return fs;
}

////////////////////////////////////////////////////////////////////////////////

void JHexahedron :: tesseract(  JNodeSequence &newnodes, JCellSequence &newcells)
{
    JNodePtr vtx;
    newnodes.reserve(8);
    newcells.reserve(5);

    double xc[8], yc[8], zc[8];

    Point3D xyz;
    for( int i = 0; i < 8; i++) {
        xyz   =  nodes[i]->getXYZCoords();
        xc[i] =  xyz[0];
        yc[i] =  xyz[1];
        zc[i] =  xyz[2];
    }

    double r,s, t;

    r = -0.5;
    s = -0.5;
    t = -0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    r =  0.5;
    s = -0.5;
    t = -0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    r =  0.5;
    s =  0.5;
    t = -0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    r = -0.5;
    s =  0.5;
    t = -0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);


    r = -0.5;
    s = -0.5;
    t =  0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    r =  0.5;
    s = -0.5;
    t =  0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    r =  0.5;
    s =  0.5;
    t =  0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    r = -0.5;
    s =  0.5;
    t =  0.5;
    xyz[0] =  TFI::trilinear_interpolation(r, s, t, xc);
    xyz[1] =  TFI::trilinear_interpolation(r, s, t, yc);
    xyz[2] =  TFI::trilinear_interpolation(r, s, t, zc);
    vtx = JNode::newObject();
    vtx->setXYZCoords(xyz);
    newnodes.push_back(vtx);

    JHexahedronPtr hex;
    JNodeSequence hexnodes(8);

    hexnodes[0] = nodes[0];
    hexnodes[1] = nodes[3];
    hexnodes[2] = nodes[7];
    hexnodes[3] = nodes[4];
    hexnodes[4] = newnodes[0];
    hexnodes[5] = newnodes[3];
    hexnodes[6] = newnodes[7];
    hexnodes[7] = newnodes[4];
    hex = JHexahedron::newObject();
    hex->setNodes(  hexnodes );
    newcells.push_back( hex);


    hexnodes[0] = nodes[1];
    hexnodes[1] = nodes[2];
    hexnodes[2] = nodes[6];
    hexnodes[3] = nodes[5];
    hexnodes[4] = newnodes[1];
    hexnodes[5] = newnodes[2];
    hexnodes[6] = newnodes[6];
    hexnodes[7] = newnodes[5];
    hex = JHexahedron::newObject();
    hex->setNodes(  hexnodes );
    newcells.push_back( hex);

    hexnodes[0] = nodes[0];
    hexnodes[1] = nodes[1];
    hexnodes[2] = nodes[5];
    hexnodes[3] = nodes[4];
    hexnodes[4] = newnodes[0];
    hexnodes[5] = newnodes[1];
    hexnodes[6] = newnodes[5];
    hexnodes[7] = newnodes[4];
    hex = JHexahedron::newObject();
    hex->setNodes(  hexnodes );
    newcells.push_back( hex);


    hexnodes[0] = nodes[2];
    hexnodes[1] = nodes[3];
    hexnodes[2] = nodes[7];
    hexnodes[3] = nodes[6];
    hexnodes[4] = newnodes[2];
    hexnodes[5] = newnodes[3];
    hexnodes[6] = newnodes[7];
    hexnodes[7] = newnodes[6];
    hex = JHexahedron::newObject();
    hex->setNodes(  hexnodes );
    newcells.push_back( hex);

    hexnodes[0] = nodes[0];
    hexnodes[1] = nodes[1];
    hexnodes[2] = nodes[2];
    hexnodes[3] = nodes[3];
    hexnodes[4] = newnodes[0];
    hexnodes[5] = newnodes[1];
    hexnodes[6] = newnodes[2];
    hexnodes[7] = newnodes[3];
    hex = JHexahedron::newObject();
    hex->setNodes(  hexnodes );
    newcells.push_back( hex);

    hexnodes[0] = newnodes[4];
    hexnodes[1] = newnodes[5];
    hexnodes[2] = newnodes[6];
    hexnodes[3] = newnodes[7];
    hexnodes[4] = nodes[4];
    hexnodes[5] = nodes[5];
    hexnodes[6] = nodes[6];
    hexnodes[7] = nodes[7];
    hex = JHexahedron::newObject();
    hex->setNodes(  hexnodes );
    newcells.push_back( hex);
}

///////////////////////////////////////////////////////////////////////////////
void JHexahedron :: get_parallel_face_diagonal( const JNodePtr &v0, const JNodePtr &v1,
        JNodePtr &v2, JNodePtr &v3)
{
    v2.reset();
    v3.reset();
    int pos0 = JCell::getPosOf(v0);
    int pos1 = JCell::getPosOf(v1);

    // On face 1;
    if( pos0 == 0  && pos1 == 7 ) {
        v2 = getNodeAt( 1 );
        v3 = getNodeAt( 6 );
        return;
    }
    if( pos0 == 7  && pos1 == 0 ) {
        v2 = getNodeAt( 6 );
        v3 = getNodeAt( 1 );
        return;
    }
    if( pos0 == 3  && pos1 == 4 ) {
        v2 = getNodeAt( 2 );
        v3 = getNodeAt( 5 );
        return;
    }
    if( pos0 == 4  && pos1 == 3 ) {
        v2 = getNodeAt( 5 );
        v3 = getNodeAt( 2 );
        return;
    }

    // On face 2;
    if( pos0 == 1  && pos1 == 6 ) {
        v2 = getNodeAt( 0 );
        v3 = getNodeAt( 7 );
        return;
    }
    if( pos0 == 6  && pos1 == 1 ) {
        v2 = getNodeAt( 7 );
        v3 = getNodeAt( 0 );
        return;
    }
    if( pos0 == 2  && pos1 == 5 ) {
        v2 = getNodeAt( 3 );
        v3 = getNodeAt( 4 );
        return;
    }
    if( pos0 == 5  && pos1 == 2 ) {
        v2 = getNodeAt( 4 );
        v3 = getNodeAt( 3 );
        return;
    }

    // On face 3;
    if( pos0 == 0  && pos1 == 5 ) {
        v2 = getNodeAt( 3 );
        v3 = getNodeAt( 6 );
        return;
    }
    if( pos0 == 5  && pos1 == 0 ) {
        v2 = getNodeAt( 6 );
        v3 = getNodeAt( 3 );
        return;
    }
    if( pos0 == 1  && pos1 == 4 ) {
        v2 = getNodeAt( 2 );
        v3 = getNodeAt( 7 );
        return;
    }
    if( pos0 == 4  && pos1 == 1 ) {
        v2 = getNodeAt( 7 );
        v3 = getNodeAt( 2 );
        return;
    }

    // Nodes on face 4
    if( pos0 == 3  && pos1 == 6 ) {
        v2 = getNodeAt( 0 );
        v3 = getNodeAt( 5 );
        return;
    }
    if( pos0 == 6  && pos1 == 3 ) {
        v2 = getNodeAt( 5 );
        v3 = getNodeAt( 0 );
        return;
    }

    if( pos0 == 2  && pos1 == 7 ) {
        v2 = getNodeAt( 1 );
        v3 = getNodeAt( 4 );
        return;
    }
    if( pos0 == 7  && pos1 == 2 ) {
        v2 = getNodeAt( 4 );
        v3 = getNodeAt( 1 );
        return;
    }

    // Nodes on face 5
    if( pos0 == 0  && pos1 == 2 ) {
        v2 = getNodeAt( 4 );
        v3 = getNodeAt( 6 );
        return;
    }
    if( pos0 == 2  && pos1 == 0 ) {
        v2 = getNodeAt( 6 );
        v3 = getNodeAt( 4 );
        return;
    }

    if( pos0 == 1  && pos1 == 3 ) {
        v2 = getNodeAt( 5 );
        v3 = getNodeAt( 7 );
        return;
    }
    if( pos0 == 3  && pos1 == 1 ) {
        v2 = getNodeAt( 7 );
        v3 = getNodeAt( 5 );
        return;
    }

    // Nodes on face 6
    if( pos0 == 4  && pos1 == 6 ) {
        v2 = getNodeAt( 0 );
        v3 = getNodeAt( 2 );
        return;
    }
    if( pos0 == 6  && pos1 == 4 ) {
        v2 = getNodeAt( 2 );
        v3 = getNodeAt( 0 );
        return;
    }
    if( pos0 == 5  && pos1 == 7 ) {
        v2 = getNodeAt( 1 );
        v3 = getNodeAt( 3 );
        return;
    }
    if( pos0 == 7  && pos1 == 5 ) {
        v2 = getNodeAt( 3 );
        v3 = getNodeAt( 1 );
        return;
    }
}
///////////////////////////////////////////////////////////////////////////////
JNodePtr JHexahedron :: getDiagonalNode( const JNodePtr &vtx) const
{

    int pos = JCell::getPosOf(vtx);
    if( pos < 0 || pos > 7) return NULL;

    switch( pos) {
    case 0:
        return getNodeAt(6);
    case 1:
        return getNodeAt(7);
    case 2:
        return getNodeAt(4);
    case 3:
        return getNodeAt(5);
    case 4:
        return getNodeAt(2);
    case 5:
        return getNodeAt(3);
    case 6:
        return getNodeAt(0);
    case 7:
        return getNodeAt(1);
    }

    return NULL;
}

////////////////////////////////////////////////////////////////////////////////

JEdgePtr JHexahedron :: getDiagonalEdge( const JEdgePtr &edge) const
{
    JEdgeSequence paredges;
    get_topological_parallel_edges(edge, paredges);
    if( paredges.size() != 3 ) return NULL;
    return paredges[1];
}


////////////////////////////////////////////////////////////////////////////////
int JHexahedron :: dice( JEdgePtr &edge, int npieces,
                        JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.clear();

    JNoImpl();
// if( !hasEdge(edge)) return 1;

    if( npieces < 2) return 1;

    JEdgeSequence paredges;
    get_topological_parallel_edges( edge, paredges);
    assert( paredges.size() == 3);

    paredges.insert(paredges.begin(), edge);

    int nNodes = npieces-1;

    JNodeSequence edgenodes[4];

    for( int i = 0; i < 4; i++) {
        if( !paredges[i]->hasAttribute("Steiner") ) {
            paredges[i]->getMidNodes( nNodes, edgenodes[i]);
            paredges[i]->setAttribute("Steiner", edgenodes[i] );
            for( size_t j = 0; j < edgenodes[i].size(); j++)
                newnodes.push_back( edgenodes[i][j] );
        }
        paredges[i]->getAttribute("Steiner", edgenodes[i]);
    }

    JNodeSequence currnodes(4);

    int refdir = getOrientation(edge);

    for( int i = 0; i < 4; i++) {
        int dir =  getOrientation( paredges[i] );
        if( dir*refdir == -1)  {
            currnodes[i] = paredges[i]->getNodeAt(1);
            edgenodes[i].push_back( paredges[i]->getNodeAt(0));
            std::reverse( edgenodes[i].begin(), edgenodes[i].end() );
        } else {
            currnodes[i] = paredges[i]->getNodeAt(0);
            edgenodes[i].push_back( paredges[i]->getNodeAt(1));
        }
    }

    newcells.resize( npieces);
    JNodeSequence hnodes(8);
    for( int i = 0; i < npieces; i++) {
        hnodes[0]  = currnodes[0];
        hnodes[1]  = currnodes[1];
        hnodes[2]  = currnodes[2];
        hnodes[3]  = currnodes[3];
        hnodes[4]  = edgenodes[0][i];
        hnodes[5]  = edgenodes[1][i];
        hnodes[6]  = edgenodes[2][i];
        hnodes[7]  = edgenodes[3][i];
        JHexahedronPtr h = JHexahedron::newObject();
        h->setNodes( hnodes );
        newcells[i] = h;
        currnodes[0] = hnodes[4];
        currnodes[1] = hnodes[5];
        currnodes[2] = hnodes[6];
        currnodes[3] = hnodes[7];
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////////////
int JHexRefiner::refine(const JCellPtr &hex, JNodeSequence &newnodes, JCellSequence &newcells)
{
    newnodes.clear();
    newcells.clear();

    if( !hex->isActive() ) return 1;

    JNodeSequence hexnodes(27);

    JNodePtr v0 = hex->getNodeAt(0);
    JNodePtr v1 = hex->getNodeAt(1);
    JNodePtr v2 = hex->getNodeAt(2);
    JNodePtr v3 = hex->getNodeAt(3);
    JNodePtr v4 = hex->getNodeAt(4);
    JNodePtr v5 = hex->getNodeAt(5);
    JNodePtr v6 = hex->getNodeAt(6);
    JNodePtr v7 = hex->getNodeAt(7);

    hexnodes[0] = hex->getNodeAt(0);
    hexnodes[2] = hex->getNodeAt(1);
    hexnodes[8] = hex->getNodeAt(2);
    hexnodes[6] = hex->getNodeAt(3);

    hexnodes[18] = hex->getNodeAt(4);
    hexnodes[20] = hex->getNodeAt(5);
    hexnodes[26] = hex->getNodeAt(6);
    hexnodes[24] = hex->getNodeAt(7);

    JNodePtr vertex;
    JEdgePtr edge;
    JFacePtr face;

    edge = JSimplex::getEdgeOf(v0,v1);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v0,v1);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[1] );

    edge = JSimplex::getEdgeOf(v0,v3);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v0,v3);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[3] );

    edge = JSimplex::getEdgeOf(v1,v2);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v1,v2);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[5] );

    edge = JSimplex::getEdgeOf(v2,v3);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v2,v3);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[7] );

    edge = JSimplex::getEdgeOf(v0,v4);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v0,v4);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[9] );

    edge = JSimplex::getEdgeOf(v1,v5);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v1,v5);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[11] );

    edge = JSimplex::getEdgeOf(v3,v7);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v3,v7);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[15] );

    edge = JSimplex::getEdgeOf(v2,v6);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v2,v6);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[17] );

    edge = JSimplex::getEdgeOf(v4,v5);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v4,v5);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[19] );

    edge = JSimplex::getEdgeOf(v4,v7);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v4,v7);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[21] );

    edge = JSimplex::getEdgeOf(v5,v6);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v5,v6);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[23] );

    edge = JSimplex::getEdgeOf(v6,v7);
    assert(edge);
    if( !edge->hasAttribute("Steiner") ) {
        vertex = JNodeGeometry::getMidNode(v6,v7);
        newnodes.push_back(vertex);
        edge->setAttribute("Steiner", vertex);
    }
    edge->getAttribute("Steiner", hexnodes[25] );

    Point3D xyz;
    face = JSimplex::getFaceOf(v0,v1,v2,v3);
    assert(face);
    if( !face->hasAttribute("Steiner") ) {
        face->getAvgXYZ(xyz);
        vertex = JNode::newObject();
        vertex->setXYZCoords(xyz);
        newnodes.push_back(vertex);
        face->setAttribute("Steiner", vertex);
    }
    face->getAttribute("Steiner", hexnodes[4] );

    face = JSimplex::getFaceOf(v4,v5,v6,v7);
    assert(face);
    if( !face->hasAttribute("Steiner") ) {
        face->getAvgXYZ(xyz);
        vertex = JNode::newObject();
        vertex->setXYZCoords(xyz);
        newnodes.push_back(vertex);
        face->setAttribute("Steiner", vertex);
    }
    face->getAttribute("Steiner", hexnodes[22] );

    face = JSimplex::getFaceOf(v0,v4,v7,v3);
    assert(face);
    if( !face->hasAttribute("Steiner") ) {
        face->getAvgXYZ(xyz);
        vertex = JNode::newObject();
        vertex->setXYZCoords(xyz);
        newnodes.push_back(vertex);
        face->setAttribute("Steiner", vertex);
    }
    face->getAttribute("Steiner", hexnodes[12] );

    face = JSimplex::getFaceOf(v1,v2,v6,v5);
    assert(face);
    if( !face->hasAttribute("Steiner") ) {
        face->getAvgXYZ(xyz);
        vertex = JNode::newObject();
        vertex->setXYZCoords(xyz);
        newnodes.push_back(vertex);
        face->setAttribute("Steiner", vertex);
    }
    face->getAttribute("Steiner", hexnodes[14] );

    face = JSimplex::getFaceOf(v0,v1,v5,v4);
    assert(face);
    if( !face->hasAttribute("Steiner") ) {
        face->getAvgXYZ(xyz);
        vertex = JNode::newObject();
        vertex->setXYZCoords(xyz);
        newnodes.push_back(vertex);
        face->setAttribute("Steiner", vertex);
    }
    face->getAttribute("Steiner", hexnodes[10] );

    face = JSimplex::getFaceOf(v2,v3,v7,v6);
    assert(face);
    if( !face->hasAttribute("Steiner") ) {
        face->getAvgXYZ(xyz);
        vertex = JNode::newObject();
        vertex->setXYZCoords(xyz);
        newnodes.push_back(vertex);
        face->setAttribute("Steiner", vertex);
    }
    face->getAttribute("Steiner", hexnodes[16] );

    hex->getAvgXYZ(xyz);
    vertex = JNode::newObject();
    vertex->setXYZCoords(xyz);
    newnodes.push_back(vertex);
    hexnodes[13] = vertex;

    int nx = 3;
    int ny = 3;
    int index, offset;
    int hexid = 0;

    JNodeSequence newhexnodes(8);
    newcells.resize(8);

    for( int k = 0; k < 2; k++) {
        for( int j = 0; j < 2; j++) {
            for( int i = 0; i < 2; i++) {
                offset = k*nx*ny + j*nx + i;
                index  = offset;
                newhexnodes[0] = hexnodes[index];

                index = offset + 1;
                newhexnodes[1] = hexnodes[index];

                index = offset + 1 + nx;
                newhexnodes[2] = hexnodes[index];

                index = offset + nx;
                newhexnodes[3] = hexnodes[index];

                offset = (k+1)*nx*ny + j*nx + i;
                index  = offset;
                newhexnodes[4] = hexnodes[index];

                index = offset + 1;
                newhexnodes[5] = hexnodes[index];

                index = offset + 1 + nx;
                newhexnodes[6] = hexnodes[index];

                index = offset + nx;
                newhexnodes[7] = hexnodes[index];
                JHexahedronPtr h = JHexahedron::newObject();
                h->setNodes( newhexnodes );
                newcells[hexid++] = h;
            }
        }
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

JTriangularPrismPtr JTriangularPrism::getCanonical()
{
    JNodeSequence  vnodes(6);

    for( int i = 0; i < 6; i++) vnodes[i] = JNode::newObject();

    Point3D p3d;

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    vnodes[0]->setXYZCoords( p3d );

    p3d[0] = 1.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    vnodes[1]->setXYZCoords( p3d );

    p3d[0] = 0.0;
    p3d[1] = 1.0;
    p3d[2] = 0.0;
    vnodes[2]->setXYZCoords( p3d );

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 1.0;
    vnodes[3]->setXYZCoords( p3d );

    p3d[0] = 1.0;
    p3d[1] = 0.0;
    p3d[2] = 1.0;
    vnodes[4]->setXYZCoords( p3d );

    p3d[0] = 0.0;
    p3d[1] = 1.0;
    p3d[2] = 1.0;
    vnodes[5]->setXYZCoords( p3d );

    JTriangularPrismPtr tp( new JTriangularPrism());
    tp->setNodes(vnodes);
    return tp;
}

///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JTriangularPrism :: getEdges() const
{
    JEdgeSequence edges(NumEdges);
    JNodePtr v0, v1;

    // First Edge ...
    v0 = nodes[0];
    v1 = nodes[1];
    edges[0] = JSimplex::getEdgeOf(v0, v1,1);

    // 2nd  Edge ...
    v0 = nodes[1];
    v1 = nodes[2];
    edges[1] = JSimplex::getEdgeOf(v0, v1,1);

    // 3rd Edge ...
    v0 = nodes[0];
    v1 = nodes[2];
    edges[2]  = JSimplex::getEdgeOf(v0, v1,1);

    // 4th Edge ...
    v0 = nodes[3];
    v1 = nodes[4];
    edges[3] = JSimplex::getEdgeOf(v0, v1,1);

    // 5th Edge ..
    v0 = nodes[4];
    v1 = nodes[5];
    edges[4] = JSimplex::getEdgeOf(v0, v1, 1);

    // 6th Edge ...
    v0 = nodes[3];
    v1 = nodes[5];
    edges[5] = JSimplex::getEdgeOf(v0, v1, 1);

    // 7th Edge ...
    v0 = nodes[0];
    v1 = nodes[3];
    edges[6] = JSimplex::getEdgeOf(v0, v1, 1);

    // 8th Edge ...
    v0 = nodes[1];
    v1 = nodes[4];
    edges[7]  = JSimplex::getEdgeOf(v0, v1,1);

    // 9th Edge ...
    v0 = nodes[2];
    v1 = nodes[5];
    edges[8] = JSimplex::getEdgeOf(v0, v1, 1);

    return edges;
}

///////////////////////////////////////////////////////////////////////////////

JFaceSequence JTriangularPrism :: getFaces()  const
{
    JFaceSequence faces(NumFaces);

    JNodePtr v0, v1, v2, v3;

    v0 = nodes[0];
    v1 = nodes[1];
    v2 = nodes[2];
    faces[0] = JSimplex::getFaceOf(v0,v1,v2, 1);

    // Triangle opposite vertex "1"
    v0 = nodes[3];
    v1 = nodes[4];
    v2 = nodes[5];
    faces[1] = JSimplex::getFaceOf(v0,v1,v2, 1);

    // Quad  opposite vertex "2"
    v0 = nodes[0];
    v1 = nodes[1];
    v2 = nodes[4];
    v3 = nodes[3];
    faces[2] = JSimplex::getFaceOf(v0,v1,v2, v3, 1);

    v0 = nodes[1];
    v1 = nodes[2];
    v2 = nodes[5];
    v3 = nodes[4];
    faces[3] = JSimplex::getFaceOf(v0,v1,v2, v3, 1);

    v0 = nodes[0];
    v1 = nodes[2];
    v2 = nodes[5];
    v3 = nodes[3];
    faces[4] = JSimplex::getFaceOf(v0,v1,v2, v3, 1);
    return faces;
}

//////////////////////////////////////////////////////////////////////////////////

int JTriangularPrism :: build_lower_entities( int )
{
    JExit();
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int JPolyhedron :: build_lower_entities( int  )
{
    JExit();
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////
bool JCellGeometry ::  isInverted( const JCellPtr &c)
{
    JMeshQuality q;
    double  J = q.getJacobian(c);
    if( J < 0.0) return 1;
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////
