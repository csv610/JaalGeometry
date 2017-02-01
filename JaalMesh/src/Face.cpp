#include <iomanip>

#include "Mesh.hpp"
#include "basic_math.hpp"
#include "MeshAffineTransforms.hpp"

using namespace std;
using namespace Jaal;

size_t JFace::NumObjectsCreated = 0;
std::map<string,string> JFace::attribInfo;

///////////////////////////////////////////////////////////////////////////////

int JFace :: registerAttribute( const string &name, const string &type)
{
    int  found = 0;

    if( type =="int"    ) found = 1;
    if( type =="char"   ) found = 1;
    if( type =="float"  ) found = 1;
    if( type =="double" ) found = 1;
    if( type =="uchar"  ) found = 1;

    if( !found) {
        cout << "Warning: invalid attribute type " << endl;
        return 2;
    }

    attribInfo[name] = type;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

string JFace :: getAttributeTypeName( const string &name)
{
    string str;
    if( attribInfo.find( name ) == attribInfo.end() ) return str;
    return  attribInfo[name];
}
///////////////////////////////////////////////////////////////////////////////

void
JFace::getRelations02( const JFacePtr &iface, JFaceSequence &faceneighs)
{
    faceneighs.clear();

    if( iface == nullptr) return;
    JFaceSequence vneighs;
    int nSize = iface->getSize(0);

    for (int i = 0; i < nSize; i++) {
        const JNodePtr &vi = iface->getNodeAt(i);
        JNode::getRelations(vi, vneighs);
        assert( !vneighs.empty() ) ;
        for( const JFacePtr  &jface : vneighs) {
            if (iface != jface ) {
                if (find(faceneighs.begin(), faceneighs.end(), jface) == faceneighs.end())
                    faceneighs.push_back(jface);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void
JFace::getRelations12( const JFacePtr &iface, JFaceSequence &faceneighs)
{
    faceneighs.clear();
    if( iface == nullptr) return;

    JFaceSequence fneighs;
    int nSize = iface->getSize(0);
    for (int i = 0; i < nSize; i++) {
        const JEdgePtr &edge = iface->getEdgeAt(i);
        JEdge::getRelations(edge, fneighs);
        for( const JFacePtr &jface: fneighs) {
            if ( jface != iface ) {
                if (find(faceneighs.begin(), faceneighs.end(), jface) == faceneighs.end())
                    faceneighs.push_back(jface);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void
JFace::getRelations( const JFacePtr &iface, JCellSequence &cellneighs)
{
    cellneighs.clear();
    if( iface == nullptr) return;
    iface->getRelations_(cellneighs);
}
///////////////////////////////////////////////////////////////////////////////
int JFaceGeometry :: getOrientation2D( const JFacePtr &face)
{
    int np = face->getSize(0);
    vector<double> x(np), y(np);
    for( int i = 0; i < np; i++) {
        const Point3D &p = face->getNodeAt(i)->getXYZCoords();
        x[i] = p[0];
        y[i] = p[1];
    }
    double area = JGeometry::getSignedArea( &x[0], &y[0], np);

    if( area > 0.0) return  1;
    if( area < 0.0) return -1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool JFaceGeometry :: isInverted(const JFacePtr &face)
{
    JMeshQuality quality;
    if( JFaceGeometry :: getSignedArea(face) < 0.0) return 1;
    /*
        double val = quality.getScaledJacobian(face);
        if( val < 0.0) return 1;
    */

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
Point3D JTriGeometry :: getRandomPoint( const JFacePtr &face)
{
    assert( face->getSize(0) == 3);

    Point3D bary;

    double highval = 1.0;
    bary[0] = JMath::random_value(0.0, highval);
    highval -= bary[0];

    bary[1] = JMath::random_value(0.0, highval);
    highval -= bary[1];

    bary[2] = JMath::random_value(0.0, highval);
    highval -= bary[2];

    double sum = bary[0] + bary[1] + bary[2];

    bary[0] /= (double) sum;
    bary[1] /= (double) sum;
    bary[2] /= (double) sum;

    assert( bary[0] >= 0.0 && bary[0] <= 1.0);
    assert( bary[1] >= 0.0 && bary[1] <= 1.0);
    assert( bary[2] >= 0.0 && bary[2] <= 1.0);

    Point3D xyz;
    JTriGeometry::getXYZCoordinates(face, bary, xyz);
    return xyz;
}

double JTriGeometry :: getSignedArea(const double *p0, const double *p1, const double *p2)
{
    double x[3], y[3];

    x[0] = p0[0];
    x[1] = p1[0];
    x[2] = p2[0];

    y[0] = p0[1];
    y[1] = p1[1];
    y[2] = p2[1];
    return JGeometry::getSignedArea(x,y,3);
}
////////////////////////////////////////////////////////////////////////////////
int JTriGeometry :: getXYZCoordinates( const JFacePtr &face, const Point3D &baryCoords, Point3D &xyz)
{
    assert(face->getSize(0) == 3);
    Point3D points[3];

    points[0] = face->getNodeAt(0)->getXYZCoords();
    points[1] = face->getNodeAt(1)->getXYZCoords();
    points[2] = face->getNodeAt(2)->getXYZCoords();

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    for( int i = 0; i < 3; i++) {
        xyz[0] += baryCoords[i]*points[i][0];
        xyz[1] += baryCoords[i]*points[i][1];
        xyz[2] += baryCoords[i]*points[i][2];
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JTriGeometry :: getXYZCoordinates( const double *p0, const double *p1, const double *p2,
                                       const double *uvw, double *xyz)
{
    xyz[0] = uvw[0]*p0[0] + uvw[1]*p1[0] + uvw[2]*p2[0];
    xyz[1] = uvw[0]*p0[1] + uvw[1]*p1[1] + uvw[2]*p2[1];
    xyz[2] = uvw[0]*p0[2] + uvw[1]*p1[2] + uvw[2]*p2[2];

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JTriGeometry :: getBaryCoords( const double *p0, const double *p1, const double *p2,
                                   const double *xy, double *uv)
{
    double totalArea = JTriGeometry::getSignedArea( p0, p1, p2);
    assert( fabs(totalArea) > 1.0E-10);

    uv[1] = JTriGeometry::getSignedArea( xy, p2, p0)/totalArea;
    uv[2] = JTriGeometry::getSignedArea( xy, p0, p1)/totalArea;
    uv[0] = 1.0- uv[1] - uv[2];

    if( uv[1] < 0.0 || uv[1] > 1.0) return 1;
    if( uv[2] < 0.0 || uv[2] > 1.0) return 1;

    return 0;
}
////////////////////////////////////////////////////////////////////////////////
int JTriGeometry :: getUVCoords( const double *p0, const double *p1, const double *p2,
                                 const double *xy, double *uv)
{
    double eps = 1.0E-10;
    double totalArea = JTriGeometry::getSignedArea( p0, p1, p2);
    assert( fabs(totalArea) > 1.0E-10);

    uv[0] = JTriGeometry::getSignedArea( xy, p2, p0)/totalArea;
    uv[1] = JTriGeometry::getSignedArea( xy, p0, p1)/totalArea;

    if( uv[0] < -eps || uv[0] > 1.0+eps) return 1;
    if( uv[1] < -eps || uv[1] > 1.0+eps) return 1;

    double u3 = 1 - uv[0] - uv[1];
    if( u3 < -eps ||   u3 > 1.0 + eps ) return 1;

    uv[0] = fabs(uv[0]);
    uv[1] = fabs(uv[1]);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JTriGeometry :: getBaryCoordinates( const JFacePtr &face, const Point3D &xy, Point3D &uvw)
{
    assert(face->getSize(0) == 3);

    const Point3D &pa = face->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = face->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = face->getNodeAt(2)->getXYZCoords();

    JTriGeometry :: getBaryCoords( &pa[0], &pb[0],&pc[0], &xy[0],  &uvw[0]);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

bool JTriGeometry :: isInside( const JFacePtr &face, const Point3D &queryPoint, bool include_boundary)
{
    assert( face ) ;
    assert(face->getSize(0) == 3);

    const Point3D &pa = face->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = face->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = face->getNodeAt(2)->getXYZCoords();

    double val1 = orient2dexact( &pa[0], &pb[0], &queryPoint[0]);
    if( val1 < 0 ) return 0;

    double val2 = orient2dexact( &pb[0], &pc[0], &queryPoint[0]);
    if( val2 < 0 ) return 0;

    double val3 = orient2dexact( &pc[0], &pa[0], &queryPoint[0]);
    if( val3 < 0 ) return 0;

    if( include_boundary ) {
        if( val1 == 0.0 || val2 == 0.0  || val3 == 0.0) return 1;
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

bool JTriGeometry :: isOutside( const JFacePtr &face, const Point3D &queryPoint, bool include_boundary)
{
    assert( face ) ;

    assert(face->getSize(0) == 3);

    const Point3D &pa = face->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = face->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = face->getNodeAt(2)->getXYZCoords();

    double val1 = orient2dexact( &pa[0], &pb[0], &queryPoint[0]);
    if( val1 < 0 ) return 1;

    double val2  = orient2dexact( &pb[0], &pc[0], &queryPoint[0]);
    if( val2 < 0 ) return 1;

    double val3 = orient2dexact( &pc[0], &pa[0], &queryPoint[0]);
    if( val3 < 0 ) return 1;

    if( include_boundary ) {
        if( val1 == 0.0 || val2 == 0.0  || val3 == 0.0) return 1;
    }

    return 0;

}
///////////////////////////////////////////////////////////////////////////////

bool JTriGeometry :: isOnBoundary( const JFacePtr &face, const Point3D &queryPoint)
{
    assert(face->getSize(0) == 3);

    const Point3D &pa = face->getNodeAt(0)->getXYZCoords();
    const Point3D &pb = face->getNodeAt(1)->getXYZCoords();
    const Point3D &pc = face->getNodeAt(2)->getXYZCoords();

    double vol;
    vol = orient2dexact( &pa[0], &pb[0], &queryPoint[0]);
    if( vol == 0 ) return 1;

    vol = orient2dexact( &pb[0], &pc[0], &queryPoint[0]);
    if( vol == 0 ) return 1;

    vol = orient2dexact( &pc[0], &pa[0], &queryPoint[0]);
    if( vol == 0 ) return 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JTriGeometry :: getMaxAngle( const JFacePtr &face, double &angle, int &pos)
{
    assert( face->getSize(0) == 3);

    const Point3D &p0 = face->getNodeAt(0)->getXYZCoords();
    const Point3D &p1 = face->getNodeAt(1)->getXYZCoords();
    const Point3D &p2 = face->getNodeAt(2)->getXYZCoords();

    double a   = JMath::length(p1,p2);
    double b   = JMath::length(p0,p2);
    double c   = JMath::length(p0,p1);

    double maxlen = a;
    pos = 0;
    if( b > maxlen) {
        pos = 1;
        maxlen = b;
    }

    if( c > maxlen) {
        pos = 2;
        maxlen = c;
    }

    double t;
    switch(pos)
    {
    case 0:
        t = (b*b + c*c - a*a)/(2.0*b*c);
        break;
    case 1:
        t = (a*a + c*c - b*b)/(2.0*a*c);
        break;
    case 2:
        t = (a*a + b*b - c*c)/(2.0*a*b);
        break;
    }
    if( t >  1.0) t =  1.0;
    if( t < -1.0) t = -1.0;
    angle  = 180.0*acos(t)/M_PI;
}
///////////////////////////////////////////////////////////////////////////////
bool JFace::lexiCompare(const JFacePtr face1, const JFacePtr face2)
{
    int num1 = face1->getSize(0);
    int num2 = face2->getSize(0);

    if( num1 < num2) return 1;

    vector<size_t> id1(num1);
    for( int i = 0; i < num1; i++)
        id1[i] = face1->getNodeAt(i)->getID();
    boost::sort( id1 );

    vector<size_t> id2(num2);
    for( int i = 0; i < num2; i++)
        id2[i] = face2->getNodeAt(i)->getID();
    boost::sort(id2);

    for( int i = 0; i < num1; i++)
        if( id2[i] > id1[i] ) return 1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JFaceGeometry::getBoundedSide( const JFacePtr &face, const Point3D &ptest)
{
    if( face->getSize(0) == 3 ) {
        const Point3D &p0 = face->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = face->getNodeAt(1)->getXYZCoords();
        const Point3D &p2 = face->getNodeAt(2)->getXYZCoords();
        return JGeometry::getBoundedSide(&p0[0], &p1[0], &p2[0], &ptest[0] );
    }

    cout << "Fatal error: not yet implemented " << __LINE__ << endl;
    exit(0);

}

////////////////////////////////////////////////////////////////////////////////
JFaceSequence JTriangle::newObjects(size_t n)
{
    JFaceSequence faces;
    if( n < 1) return faces;
    faces.resize(n);
    for( size_t i = 0; i < n; i++)
        faces[i] = JTriangle::newObject();
    return faces;
}
////////////////////////////////////////////////////////////////////////////////

JFaceSequence JQuadrilateral::newObjects(size_t n)
{
    JFaceSequence faces;
    if( n < 1) return faces;
    faces.resize(n);
    for( size_t i = 0; i < n; i++)
        faces[i] = JQuadrilateral::newObject();
    return faces;
}

///////////////////////////////////////////////////////////////////////////////

JFaceSequence JPolygon::newObjects(size_t n)
{
    JFaceSequence faces;
    if( n < 1) return faces;
    faces.resize(n);
    for( size_t i = 0; i < n; i++)
        faces[i] = JPolygon::newObject();
    return faces;
}
///////////////////////////////////////////////////////////////////////////////

void
set_tfi_coords(int i, int j, int nx, int ny, JNodeSequence &qnodes)
{
    assert( qnodes.size() == (size_t)nx*ny );
    int offset;

    offset = 0;
    const Point3D &v00 = qnodes[offset]->getXYZCoords();

    offset = i;
    const Point3D &vr0 = qnodes[offset]->getXYZCoords();

    offset = (nx - 1);
    const Point3D &v10 = qnodes[offset]->getXYZCoords();

    offset = j*nx;
    const Point3D &v0s = qnodes[offset]->getXYZCoords();

    offset = j * nx + (nx - 1);
    const Point3D &v1s = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx;
    const Point3D &v01 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + i;
    const Point3D &vr1 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + (nx - 1);
    const Point3D &v11 = qnodes[offset]->getXYZCoords();

    Point3D vrs;

    double dr = 2.0 / (double) (nx - 1);
    double ds = 2.0 / (double) (ny - 1);

    double r = -1.0 + i*dr;
    double s = -1.0 + j*ds;
    for (int k = 0; k < 3; k++) {
        vrs[k] = TFI::transfinite_blend(r, s,
                                        v00[k], v10[k], v11[k], v01[k],
                                        vr0[k], v1s[k], vr1[k], v0s[k]);
    }
    offset = j * nx + i;
    qnodes[offset]->setXYZCoords(vrs);
}

////////////////////////////////////////////////////////////////////////////////

int JQuadGeometry :: getBilinearCoords( const JFacePtr face, const Point2D &paramCoord, Point3D &xyz)
{
    assert( face->getSize(0) == 4);
    double x[4], y[4], z[4];

    for( int i = 0; i < 4; i++) {
        const Point3D &p = face->getNodeAt(i)->getXYZCoords();
        x[i] = p[0];
        y[i] = p[1];
        z[i] = p[2];
    }
    xyz[0] = TFI::bilinear_interpolation( paramCoord[0], paramCoord[1], x );
    xyz[1] = TFI::bilinear_interpolation( paramCoord[0], paramCoord[1], y );
    xyz[2] = TFI::bilinear_interpolation( paramCoord[0], paramCoord[1], z );
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int JQuadGeometry :: getUVCoords( const JFacePtr face, const Point3D &xyz, Point2D &paramCoord)
{
    double xmin = -1.0;
    double ymin = -1.0;
    double xmax =  1.0;
    double ymax =  1.0;
    double xmid, ymid;

    double dist[4], minval;
    Point3D p;

    Point2D uv;
    for( int i = 0; i < 100; i++) {
        xmid = 0.5*(xmax+xmin);
        ymid = 0.5*(ymax+ymin);
        uv[0] = xmid;
        uv[1] = ymid;
        getBilinearCoords( face, uv, p);
        double currdist = JMath::length(p, xyz);
        if( currdist  < 1E-10) {
            paramCoord = uv;
            return 0;
        }


        uv[0] = 0.5*(xmin+xmid);
        uv[1] = 0.5*(ymin+ymid);
        getBilinearCoords( face, uv, p);
        dist[0] = JMath::length(p, xyz);
        if( dist[0]  < 1E-10) {
            paramCoord = uv;
            return 0;
        }

        uv[0] = 0.5*(xmax + xmid);
        uv[1] = 0.5*(ymin + ymid);
        getBilinearCoords( face, uv, p);
        dist[1] = JMath::length(p, xyz);
        if( dist[1]  < 1E-10) {
            paramCoord = uv;
            return 0;
        }

        uv[0] = 0.5*(xmax + xmid);
        uv[1] = 0.5*(ymax + ymid);
        getBilinearCoords( face, uv, p);
        dist[2] = JMath::length(p, xyz);
        if( dist[2]  < 1E-10) {
            paramCoord = uv;
            return 0;
        }

        uv[0] = 0.5*(xmin + xmid);
        uv[1] = 0.5*(ymax + ymid);
        getBilinearCoords( face, uv, p);
        dist[3] = JMath::length(p, xyz);
        if( dist[3]  < 1E-10) {
            paramCoord = uv;
            return 0;
        }

        minval = dist[0];
        int minpos = 0;
        for( int j = 1; j < 4; j++) {
            if( dist[j] < minval) {
                minval = dist[j];
                minpos = j;
            }
        }

        switch( minpos) {
        case 0:
            xmax = xmid;
            ymax = ymid;
            break;
        case 1:
            xmin = xmid;
            ymax = ymid;
            break;
        case 2:
            xmin = xmid;
            ymin = ymid;
            break;
        case 3:
            ymin = ymid;
            xmax = xmid;
            break;
        }
        cout << "MinPos " << minpos << endl;
        cout << xmin << " " << ymin << "  " << xmax << "  " << ymax << endl;
        cout << i << "  " << xmid << "  " << ymid << "  " << currdist << endl;
        getchar();
    }
    cout << "Warning: Parametric coordinate not found correctly " << minval << endl;
    return 1;
}

////////////////////////////////////////////////////////////////////////////////
int JFace :: setStartNode( const JNodePtr &p)
{
    int pos = getPosOf(p);

    if( pos >= 0) {

        int nsize = this->getSize(0);
        JNodeSequence rotnodes(nsize);
        for( int i = 0; i < nsize; i++)
            rotnodes[i] = this->getNodeAt(pos+i);

        nodes = rotnodes;
        return 0;
    }
    return 1;
}
////////////////////////////////////////////////////////////////////////////////

int JFace :: remove_unattached_lower_entities()
{
    if( this->getStatus() != JMeshEntity::REMOVE ) return 1;
    int nedges = this->getSize(1);
    for( int i = 0; i < nedges; i++) {
        const JEdgePtr &edge = this->getEdgeAt(i);
        if( edge->getNumHigherRelations() == 0)  {
            edge->setStatus( JMeshEntity::REMOVE );
            edge->remove_unattached_lower_entities();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

JEdgePtr JFace:: getEdgeOf( const JNodePtr &v0, const JNodePtr &v1 )
{
    if( !this->hasNode(v0)  ) return nullptr;
    if( !this->hasNode(v1)  ) return nullptr;

    int  nn = getSize(0);
    JEdgePtr e;
    for( int i = 0; i < nn; i++) {
        if( getNodeAt(i) == v0 && getNodeAt(i+1) == v1 ) {
            e =  JSimplex::getEdgeOf(v0,v1,1);
            return e;
        }
        if( getNodeAt(i) == v1 && getNodeAt(i+1) == v0 ) {
            e =  JSimplex::getEdgeOf(v0,v1,1);
            return e;
        }
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
int JFace :: getOrientation(const JEdgePtr &rhs) const
{
    const JNodePtr &vr0 = rhs->getNodeAt(0);
    const JNodePtr &vr1 = rhs->getNodeAt(1);
    if( !hasNode(vr0) ) return 0;
    if( !hasNode(vr1) ) return 0;
    int nSize = nodes.size();
    for( int i = 0; i < nSize; i++) {
        JEdgePtr lhs  = getEdgeAt(i);
        if( lhs == rhs) {
            JNodePtr v0 = getNodeAt(i);
            JNodePtr v1 = getNodeAt(i+1);
            if( v0 == vr0 && v1 == vr1 ) return  1;
            if( v0 == vr1 && v1 == vr0 ) return -1;
            return 0;
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JFace :: getOrientation( const JNodePtr &n0, const JNodePtr &n1, const JNodePtr &n2) const
{
    if( getSize(0) != 3 ) return 0;

    if( !hasNode(n0)  )   return 0;
    if( !hasNode(n1)  )   return 0;
    if( !hasNode(n2)  )   return 0;

    int pos = JSimplex::getPosOf(n0);
    const JNodePtr &v0 = getNodeAt( pos );
    const JNodePtr &v1 = getNodeAt( pos+1 );
    const JNodePtr &v2 = getNodeAt( pos+2 );

    if( v0 == n0  && v1 == n1  && v2 == n2 ) return  1;
    if( v0 == n0  && v1 == n2  && v2 == n1 ) return -1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JFace :: getOrientation( const JNodePtr &n0, const JNodePtr &n1,
                             const JNodePtr &n2, const JNodePtr &n3) const
{
    if( getSize(0) != 4 ) return 0;

    if( !hasNode(n0)  )   return 0;
    if( !hasNode(n1)  )   return 0;
    if( !hasNode(n2)  )   return 0;
    if( !hasNode(n3)  )   return 0;

    int pos = JSimplex::getPosOf(n0);

    const JNodePtr &v0 = getNodeAt( pos);
    const JNodePtr &v1 = getNodeAt( pos+1);
    const JNodePtr &v2 = getNodeAt( pos+2);
    const JNodePtr &v3 = getNodeAt( pos+3);

    if( v0 == n0 && v1 == n1 && v2 == n2 && v3 == n3 ) return   1;
    if( v0 == n0 && v1 == n3 && v2 == n2 && v3 == n1 ) return  -1;

    return 0;
}


///////////////////////////////////////////////////////////////////////////////
JFacePtr JFace :: getProduct( int type )
{
    JFacePtr newentity;
    switch( type ) {
    case JFace::TRIANGLE:
        newentity = JTriangle::newObject();
        break;
    case JFace::QUADRILATERAL:
        newentity = JQuadrilateral::newObject();
        break;
    case JFace::POLYGON:
        newentity = JPolygon::newObject();
        break;
    }
    return newentity;
}

///////////////////////////////////////////////////////////////////////////////
JFacePtr JFace :: newObject( const JNodeSequence &nodeseq )
{
    JFacePtr newentity;
    int nsize = nodeseq.size();
    switch( nsize ) {
    case 3:
        newentity = JTriangle::newObject(nodeseq);
        break;
    case 4:
        newentity = JQuadrilateral::newObject(nodeseq);
        break;
    default:
        newentity = JPolygon::newObject(nodeseq);
        break;
    }
    return newentity;
}
///////////////////////////////////////////////////////////////////////////////

void JFace::getSharedEntities( const JFacePtr &face1,
                               const JFacePtr &face2, JNodeSequence &nodes)
{
    nodes.clear();
    int nsize = face1->getSize(0);
    for( int i = 0; i < nsize; i++) {
        const JNodePtr &vtx = face1->getNodeAt(i);
        if( face2->hasNode( vtx ) ) nodes.push_back(vtx);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JFace:: getSharedEntities( const JFacePtr &face1,
                                const JFacePtr &face2, JEdgeSequence &comm_edges)
{
    JEdgeSequence face1edges, face2edges;

    face1edges = face1->getEdges();
    face2edges = face2->getEdges();

    boost::sort( face1edges );
    boost::sort( face2edges );

    comm_edges.clear();
    boost::set_intersection( face1edges, face2edges, back_inserter(comm_edges) );
}

///////////////////////////////////////////////////////////////////////////////

int JFace :: build_lower_entities(int )
{
    int numedges = this->getSize(1);
    JNodePtr vhash = nullptr;
    for( int i = 0; i < numedges; i++) {
        const JEdgePtr &fedge = this->getEdgeAt(i);
        vhash = fedge->getHashNode();
        if( vhash ) vhash->attach(fedge);
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JFacePtr JFace :: explode( double alpha) const
{
    JFacePtr fe = getClone();
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

double JFaceGeometry:: getAngleAt( const JFacePtr &face, const JNodePtr &v, int measure)
{
    int pos = face->getPosOf(v);
    assert( pos >= 0);
    int nnodes = face->getSize(0);
    double sum = 0.0;
    for( int i = 0; i < nnodes-2; i++) {
        const Point3D &p0 = face->getNodeAt(pos+i+ 0)->getXYZCoords();
        const Point3D &p1 = face->getNodeAt(pos+i+ 1)->getXYZCoords();
        const Point3D &p2 = face->getNodeAt(pos+i+nnodes-1)->getXYZCoords();
        sum += JMath::getTriAngle(p0, p1, p2, measure);
    }
    return sum;
}

///////////////////////////////////////////////////////////////////////////////

bool
JFaceGeometry::isConvex( const JFacePtr &face)
{
    int nnodes = face->getSize(0);

    if( nnodes < 4) return 1;

    int pos = JFaceGeometry::reflexAngleAt(face);
    if( pos < 0) return 1;
    return 0;

//  if( nnodes == 4 ) return QuadGeometry::isConvex(face);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool
JQuadGeometry::isConvex(const Point3D &p0, const Point3D &p1,
                        const Point3D &p2, const Point3D &p3)
{
    double qarea = JQuadGeometry::getArea(p0, p1, p2, p3);

    double tarea1, tarea2;
    tarea1 = JTriGeometry::getArea(p0, p1, p2);
    tarea2 = JTriGeometry::getArea(p0, p2, p3);
    if (fabs(tarea1 + tarea2 - qarea) > 1.0E-10) return 0;

    tarea1 = JTriGeometry::getArea(p0, p1, p3);
    tarea2 = JTriGeometry::getArea(p1, p2, p3);
    if (fabs(tarea1 + tarea2 - qarea) > 1.0E-10) return 0;

    return 1;
}

////////////////////////////////////////////////////////////////////////////////

bool JQuadGeometry :: isConvex( const JFacePtr face)
{
    assert( face->getSize(0) == 4);
    return JQuadGeometry::isConvex( face->getNodeAt(0)->getXYZCoords(),
                                    face->getNodeAt(1)->getXYZCoords(),
                                    face->getNodeAt(2)->getXYZCoords(),
                                    face->getNodeAt(3)->getXYZCoords());
}


////////////////////////////////////////////////////////////////////////////////

double
JFaceGeometry::getAspectRatio( const JFacePtr &face)
{
    int nSize = face->getSize(0);

    double minlen = MAXDOUBLE;
    double maxlen = 0.0;

    for (int i = 0; i < nSize; i++) {
        JNodePtr v0 = face->getNodeAt(i);
        JNodePtr v1 = face->getNodeAt(i+1);
        double len2 = JNodeGeometry::getLength2(v0, v1);
        if (len2 > maxlen) maxlen = len2;
        if (len2 < minlen) minlen = len2;
    }

    return sqrt(minlen / maxlen);
}

////////////////////////////////////////////////////////////////////////////////

bool
JFace::has_boundary_edge() const
{
    int nSize = nodes.size();
    JFaceSequence neighs;
    for (int i = 0; i < nSize; i++) {
        const JNodePtr &v0 = getNodeAt(i);
        const JNodePtr &v1 = getNodeAt(i + 1);
        if (v0->isBoundary() && v1->isBoundary()) {
            JEdgePtr edge = getEdgeAt(i);
//             Mesh::getRelations112(v0, v1, neighs);
            JEdge::getRelations(edge, neighs);
            if (neighs.size() == 1) return 1;
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////

int JFaceGeometry :: getDualPosition( const JFacePtr face, Point3D &xyz)
{
    int numnodes = face->getSize(0);
    if( numnodes < 3) return 1;

    // For general Polygon, we return the average node.
    if( numnodes > 3 ) {
        face->getAvgXYZ(xyz);
        return 0;
    }

    const JNodePtr &v0 = face->getNodeAt(0);
    const JNodePtr &v1 = face->getNodeAt(1);
    const JNodePtr &v2 = face->getNodeAt(2);

    Point3D pa = v0->getXYZCoords();
    Point3D pb = v1->getXYZCoords();
    Point3D pc = v2->getXYZCoords();

    // For any obtuse triangle, the circumcenter lies outside the triangle.
    // but return the mid node of the longest edge...

    double angle = 0.0;
    int pos = JMath::getMaxTriAngle(pa, pb, pc, angle);
    if( angle > 90.0) {
        switch(pos) {
        case 0:
            xyz = JNodeGeometry::getMidPoint(v1,v2);
            break;
        case 1:
            xyz = JNodeGeometry::getMidPoint(v2,v0);
            break;
        case 2:
            xyz = JNodeGeometry::getMidPoint(v0,v1);
            break;
        }
        return 0;
    }
    //
    // If the triangle is acute the circumcenter will lies within the
    // triangle: Bring any 3D triangle to a canonical oientation (z = 0)
    //
    Point3D center, param;
    TriCircumCenter3D( &pa[0], &pb[0], &pc[0], &center[0], &param[0] );

    unused_parameter(param);

    xyz[0] = center[0];
    xyz[1] = center[1];
    xyz[2] = center[2];

#ifdef DEBUG
    double r0 = JMath::length(pa, center);
    double r1 = JMath::length(pb, center);
    double r2 = JMath::length(pc, center);

    assert( fabs(r0 -r1 ) < 1.0E-06);
    assert( fabs(r1 -r2 ) < 1.0E-06);
    assert( fabs(r2 -r0 ) < 1.0E-06);
#endif

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JTrianglePtr JTriangle::getCanonical( double len)
{
    Point3D p3d;

    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    JNodePtr v0 =  JNode::newObject();
    v0->setXYZCoords(p3d);
    v0->setID(0);

    p3d[0] = len;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    JNodePtr v1 =  JNode::newObject();
    v1->setXYZCoords(p3d);
    v1->setID(1);

    p3d[0] = 0.0;
    p3d[1] = len;
    p3d[2] = 0.0;
    JNodePtr v2 =  JNode::newObject();
    v2->setXYZCoords(p3d);
    v2->setID(2);

    JTrianglePtr t = JTriangle::newObject(v0,v1,v2);
    return t;
}

////////////////////////////////////////////////////////////////////////////////

Point3D JTriGeometry::getConformal( const Point3D &pa, const Point3D &pb, const Point3D &pc, const Point3D &angles)
{
    Point3D pC;
    /*
        // //////////////////////////////////////////////////////////////////////////
        // Given three points and three corresponding angles, calculate the ideal position of the
        // third point so that the triangles have given angles. The length of the AB is fixed.
        // //////////////////////////////////////////////////////////////////////////

            assert( fabs(angles[0] + angles[1] + angles[2] - M_PI) < 1.0E-06);

            Point3D p0, p1;
            if( orient2d(  &pa[0], &pb[0], &pc[0] ) > 0) {
                p0 = pa;
                p1 = pb;
            } else  {
                p0 = pb;
                p1 = pa;
            }
            p1[0] -= p0[0];
            p1[1] -= p0[1];
            p1[2] -= p0[2];

            double rotangle = atan2( p1[1], p1[0] );

            AffineTransform::zRotatePoint(p1, -rotangle);

            assert( p1[0] > 0.0 && fabs(p1[1]) < 1.0E-06);

            double ab = JMath::length(pa, pb);
            double b  = ab*sin( angles[1] )/sin( angles[2] );

            Point3D pC;
            pC[0] = b*cos( angles[0] );
            pC[1] = b*sin( angles[0] );
            pC[2] = 0.0;

            AffineTransform::zRotatePoint(pC, rotangle);

            pC[0] += p0[0];
            pC[1] += p0[1];
            pC[2] = 0.5*(pa[2]  + pb[2]);
            return pC;
    */
    return pC;

}

///////////////////////////////////////////////////////////////////////////////

Point3D JTriGeometry :: getIdealPosition( const JFacePtr &face, const JNodePtr &vtx, double tipangle)
{
    //////////////////////////////////////////////////////////////////////////
    // Given a triangular face, locate the new position of the vertex so that it has the
    // given angle at that vertex...
    //////////////////////////////////////////////////////////////////////////

    assert( face->getSize(0) == 3 );
    tipangle *= M_PI/180.0;

    double alpha = 0.5*(M_PI - tipangle);
    int pos = face->getPosOf(vtx);
    Array3D angles;

    angles[0] =  alpha;
    angles[1] =  alpha;
    angles[2] =  tipangle;
    return JTriGeometry::getConformal( face->getNodeAt(pos+1)->getXYZCoords(),
                                       face->getNodeAt(pos+2)->getXYZCoords(),
                                       face->getNodeAt(pos+0)->getXYZCoords(),
                                       angles);
}

///////////////////////////////////////////////////////////////////////////////

JEdgePtr JTriangle :: getOppositeEdge( const JFacePtr &face, const JNodePtr &vtx)
{
    assert( face->getSize(0) == 3 );
    int pos = face->getPosOf(vtx);
    const JNodePtr &v0 = face->getNodeAt(pos + 1);
    const JNodePtr &v1 = face->getNodeAt(pos + 2);
    return JSimplex::getEdgeOf(v0,v1);
}
///////////////////////////////////////////////////////////////////////////////
Vec3D
JFaceGeometry::getNormal(const JFacePtr &face)
{
    Vec3D n;

    cout << "Exit " << endl;
    exit(0);
    return n;
}
///////////////////////////////////////////////////////////////////////////////
Vec3D
JTriGeometry::getNormal(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    Vec3D normal, p1p0, p2p0;
    JMath::make_vector(p2, p0, p2p0);
    JMath::make_vector(p1, p0, p1p0);
    JMath::cross_product(p2p0, p1p0, normal);

    double mag = JMath::magnitude(normal);
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;

    return normal;
}

////////////////////////////////////////////////////////////////////////////////

Vec3D
JTriGeometry::getNormal(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2)
{
    Vec3D n;
    JMath::normal(v0->getXYZCoords(), v1->getXYZCoords(), v2->getXYZCoords(), n);
    return n;
}

////////////////////////////////////////////////////////////////////////////////

double
JTriGeometry::getArea(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    ////////////////////////////////////////////////////////////////////////////////
    // Ref: http://mathworld.wolfram.com/HeronsFormula.html
    //
    // Heron's formular  A = sqrt(s(s-a)(s-b)(s-c)) is expensive because if require
    // three square roots for each a,b, and c.
    // Instead we will use alternate formula of "Heron" which avoids three
    // expensive square roots.
    //
    // Sorting is done to reduce the truncation errors. Very similar to Kahan's
    // Original idea.
    ////////////////////////////////////////////////////////////////////////////////

    double a  = JMath::length(p1, p2);
    double b  = JMath::length(p2, p0);
    double c  = JMath::length(p0, p1);
    double s  = 0.5*(a+b+c);
    double area = sqrt(s*(s-a)*(s-b)*(s-c));
    return area;
/*
    double d[3];
    d[0] = JMath::length2(p1, p2);
    d[1] = JMath::length2(p2, p0);
    d[2] = JMath::length2(p0, p1);

    std::sort(d, d + 3); // May be we should have optimized version than STL one

    double a2 = d[0];
    double b2 = d[1];
    double c2 = d[2];

    double area = 0.25 * sqrt(4 * a2 * b2 - (a2 + b2 - c2)*(a2 + b2 - c2));
    return area;
*/
}

///////////////////////////////////////////////////////////////////////////////////
JQuadrilateralPtr JQuadrilateral::getCanonical(const Point3D &center, double len)
{
    Point3D p3d;

    p3d[0] = center[0] -0.5*len;
    p3d[1] = center[1] -0.5*len;
    p3d[2] = center[2];
    JNodePtr v0 =  JNode::newObject();
    v0->setXYZCoords(p3d);
    v0->setID(0);

    p3d[0] =  center[0] + 0.5*len;
    p3d[1] =  center[1] - 0.5*len;
    p3d[2] =  center[2];
    JNodePtr v1 =  JNode::newObject();
    v1->setXYZCoords(p3d);
    v1->setID(1);

    p3d[0] = center[0] + 0.5*len;
    p3d[1] = center[1] + 0.5*len;
    p3d[2] = center[2];
    JNodePtr v2 =  JNode::newObject();
    v2->setXYZCoords(p3d);
    v2->setID(2);

    p3d[0] = center[0] - 0.5*len;
    p3d[1] = center[1] + 0.5*len;
    p3d[2] = center[2];
    JNodePtr v3 =  JNode::newObject();
    v3->setXYZCoords(p3d);
    v3->setID(3);

    JQuadrilateralPtr q = JQuadrilateral::newObject(v0,v1,v2, v3);
    return q;
}


JQuadrilateralPtr JQuadrilateral::getCanonical(double len)
{
    Point3D  center;
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;
    return getCanonical(center, len);
}

///////////////////////////////////////////////////////////////////////////////

JQuadrilateralPtr JQuadrilateral::getCanonical(double lx, double ly)
{
    Point3D p3d;

    p3d[0] = -0.5*lx;
    p3d[1] = -0.5*ly;
    p3d[2] = 0.0;
    JNodePtr v0 =  JNode::newObject();
    v0->setXYZCoords(p3d);
    v0->setID(0);

    p3d[0] =  0.5*lx;
    p3d[1] = -0.5*ly;
    p3d[2] =  0.0;
    JNodePtr v1 =  JNode::newObject();
    v1->setXYZCoords(p3d);
    v1->setID(1);

    p3d[0] = 0.5*lx;
    p3d[1] = 0.5*ly;
    p3d[2] = 0.0;
    JNodePtr v2 =  JNode::newObject();
    v2->setXYZCoords(p3d);
    v2->setID(2);

    p3d[0] = -0.5*lx;
    p3d[1] =  0.5*ly;
    p3d[2] =  0.0;
    JNodePtr v3 =  JNode::newObject();
    v3->setXYZCoords(p3d);
    v3->setID(3);

    JQuadrilateralPtr q = JQuadrilateral::newObject(v0,v1,v2, v3);
    return q;
}

///////////////////////////////////////////////////////////////////////////////

JPolygonPtr JPolygon::getCanonical(int n, double len)
{
    assert( n > 2);

    JNodeSequence vnodes(n);
    double dtheta =  2*M_PI/(double)n;
    double radius = 0.5*len/sin( 0.5*dtheta);

    Point3D p3d;
    for( int i = 0; i < n; i++) {
        p3d[0] = radius*cos( i*dtheta);
        p3d[1] = radius*sin( i*dtheta);
        p3d[2] = 0.0;
        JNodePtr v0 = JNode::newObject();
        v0->setXYZCoords(p3d);
        v0->setID(i);
        vnodes[i] = v0;
    }
    JPolygonPtr p = JPolygon::newObject(vnodes);
    return p;
}

///////////////////////////////////////////////////////////////////////////////


int
JQuadrilateral::quad_tessalate(const JNodeSequence &orgNodes, JNodeSequence &rotatedNodes)
{
    double A013 = JTriGeometry::getArea(orgNodes[0]->getXYZCoords(),
                                        orgNodes[1]->getXYZCoords(),
                                        orgNodes[3]->getXYZCoords());

    double A123 = JTriGeometry::getArea(orgNodes[1]->getXYZCoords(),
                                        orgNodes[2]->getXYZCoords(),
                                        orgNodes[3]->getXYZCoords());

    double A012 = JTriGeometry::getArea(orgNodes[0]->getXYZCoords(),
                                        orgNodes[1]->getXYZCoords(),
                                        orgNodes[2]->getXYZCoords());

    double A023 = JTriGeometry::getArea(orgNodes[0]->getXYZCoords(),
                                        orgNodes[2]->getXYZCoords(),
                                        orgNodes[3]->getXYZCoords());

    rotatedNodes.resize(4);
    if (fabs(A013) + fabs(A123) < fabs(A012) + fabs(A023)) {
        rotatedNodes[0] = orgNodes[1];
        rotatedNodes[1] = orgNodes[2];
        rotatedNodes[2] = orgNodes[3];
        rotatedNodes[3] = orgNodes[0];
        return 1;
    }
    rotatedNodes[0] = orgNodes[0];
    rotatedNodes[1] = orgNodes[1];
    rotatedNodes[2] = orgNodes[2];
    rotatedNodes[3] = orgNodes[3];
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JQuadrilateral::triangulate( JFaceSequence &trifaces, JNodeSequence &newnodes, int n) const
{
    trifaces.clear();
    newnodes.clear();

    if( n == 2 ) {
        JNodeSequence rotatedNodes(4);
        quad_tessalate(nodes, rotatedNodes);
        trifaces.resize(2);
        trifaces[0] = JTriangle::newObject(rotatedNodes[0], rotatedNodes[1],
                                           rotatedNodes[2] );
        trifaces[1] = JTriangle::newObject(rotatedNodes[0], rotatedNodes[2],
                                           rotatedNodes[3] );
        return 0;
    }

    Point3D pc;
    if( n == 4 ) {
        JNodePtr vc = JNode::newObject();
        getAvgXYZ(pc);
        vc->setXYZCoords(pc);
        newnodes.resize(1);
        newnodes[0] = vc;
        trifaces.resize(4);
        trifaces[0] = JTriangle::newObject(nodes[0], nodes[1], vc);
        trifaces[1] = JTriangle::newObject(nodes[1], nodes[2], vc);
        trifaces[2] = JTriangle::newObject(nodes[2], nodes[3], vc);
        trifaces[3] = JTriangle::newObject(nodes[3], nodes[0], vc);
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
bool JQuadrilateral :: isPlanar() const
{
    cout << "Not done " << endl;
    /*
       const Point3D &p0 = getNodeAt(0)->getXYZCoords();
       const Point3D &p1 = getNodeAt(1)->getXYZCoords();
       const Point3D &p2 = getNodeAt(2)->getXYZCoords();
       const Point3D &p3 = getNodeAt(3)->getXYZCoords();

       double angle1, angle2;
       angle1 = JMath::getAngle( p0, p1, p2, ANGLE_IN_DEGREES);
       angle2 = JMath::getAngle( p0, p2, p3, ANGLE_IN_DEGREES);
       if( fabs(angle1-angle2) > 1.0) return 0;

       angle1 = JMath::getAngle( p1, p2, p3, ANGLE_IN_DEGREES);
       angle2 = JMath::getAngle( p1, p3, p0, ANGLE_IN_DEGREES);
       if( fabs(angle1-angle2) > 1.0) return 0;
    */
    exit(0);

    return 1;

}
///////////////////////////////////////////////////////////////////////////////
JNodePtr JQuadrilateral :: getDiagonalNode( const JFacePtr &face, const JNodePtr &vtx)
{
    assert( face->getSize(0) == 4);
    int pos = face->getPosOf(vtx);
    assert( pos >= 0);
    return face->getNodeAt( pos + 2 );
}
///////////////////////////////////////////////////////////////////////////////

JNodePtr
JTriangle::getOppositeNode(const JFacePtr &tri, const JNodePtr &n1, const JNodePtr &n2)
{
    assert( tri->getSize(0) == 3 );

    JNodePtr tn0 = tri->getNodeAt(0);
    JNodePtr tn1 = tri->getNodeAt(1);
    JNodePtr tn2 = tri->getNodeAt(2);

    if (tn0 == n1 && tn1 == n2)
        return tn2;
    if (tn0 == n2 && tn1 == n1)
        return tn2;

    if (tn1 == n1 && tn2 == n2)
        return tn0;
    if (tn1 == n2 && tn2 == n1)
        return tn0;

    if (tn2 == n1 && tn0 == n2)
        return tn1;
    if (tn2 == n2 && tn0 == n1)
        return tn1;

    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr
JFace::getDiagonalNode(const JFacePtr &face, const JNodePtr &n1)
{
    int pos = face->getPosOf(n1);
    if (pos >= 0) {
        int size = face->getSize(0);
        if (size == 4)
            return face->getNodeAt(pos + 2);
    }
    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////

int JQuadGeometry::isCyclic(const Point3D &p0, const Point3D &p1, const Point3D &p2,
                            const Point3D &p3)
{
    //
    // http://en.wikipedia.org/wiki/Ptolemy%27s_theorem
    //

    double d02 = JMath::length(p0, p2);
    double d13 = JMath::length(p1, p3);
    double d01 = JMath::length(p0, p1);
    double d12 = JMath::length(p1, p2);
    double d23 = JMath::length(p2, p3);
    double d30 = JMath::length(p3, p0);

    double lval = d02*d13;
    double rval = d01 * d23 + d12*d30;

    if (fabs(lval - rval) < 1.0E-15) return 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void
JFace::getOppositeNodes(const JFacePtr &quad, JNodePtr &n1, JNodePtr &n2,
                        JNodePtr &n3, JNodePtr &n4)
{
    const JNodePtr  &qn0 = quad->getNodeAt(0);
    const JNodePtr  &qn1 = quad->getNodeAt(1);
    const JNodePtr  &qn2 = quad->getNodeAt(2);
    const JNodePtr  &qn3 = quad->getNodeAt(3);

    if ((qn0 == n1 && qn1 == n2) || (qn0 == n2 && qn1 == n1)) {
        n3 = qn2;
        n4 = qn3;
        return;
    }

    if ((qn1 == n1 && qn2 == n2) || (qn1 == n2 && qn2 == n1)) {
        n3 = qn0;
        n4 = qn3;
        return;
    }

    if ((qn2 == n1 && qn3 == n2) || (qn2 == n2 && qn3 == n1)) {
        n3 = qn0;
        n4 = qn1;
        return;
    }

    if ((qn3 == n1 && qn0 == n2) || (qn3 == n2 && qn0 == n1)) {
        n3 = qn1;
        n4 = qn2;
        return;
    }

    cout << " Warning: You should not come here " << endl;
    cout << " search for " << n1 << "  " << n2 << endl;
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

JQuadrilateralPtr
JQuadrilateral::newObject(const JFacePtr &tri1, const JFacePtr &tri2)
{
    if( tri1->getTypeID() != JFace::TRIANGLE || tri2->getTypeID() != JFace::TRIANGLE ) {
        cout << "Warning: Need a triangle to create a quad " << endl;
        return nullptr;
    }

    JEdgeSequence edges;
    JFace::getSharedEntities(tri1, tri2, edges);
    if( edges.size() != 1 ) {
        cout << "Fatal Error:  Faces don't share a common edge " << endl;
        return nullptr;
    }

    const JNodePtr &v0   = edges[0]->getNodeAt(0);
    const JNodePtr &v1   = edges[0]->getNodeAt(1);
    const JNodePtr &ot1  = JTriangle::getOppositeNode(tri1, v0, v1);
    const JNodePtr &ot2  = JTriangle::getOppositeNode(tri2, v0, v1);

    JQuadrilateralPtr qface = JQuadrilateral::newObject(ot1, v0, ot2, v1 );
    return qface;
}

///////////////////////////////////////////////////////////////////////////////

int
JQuadrilateral::hexagon_2_quads(const JNodeSequence &hexnodes, JFaceSequence &quads, int offset)
{
#ifdef DEBUG
    if( hexnodes.size() != 6) {
        cout << "Fatal Error: Given polygon is not hex " << endl;
        exit(0);
    }
#endif

    JFacePtr face1 = JQuadrilateral::newObject();
    JFacePtr face2 = JQuadrilateral::newObject();

    JNodeSequence nodes(4);

    for (int i = 0; i < 3; i++) {
        nodes[0] = hexnodes[ (i + offset + 0) % 6];
        nodes[1] = hexnodes[ (i + offset + 1) % 6];
        nodes[2] = hexnodes[ (i + offset + 2) % 6];
        nodes[3] = hexnodes[ (i + offset + 3) % 6];
        face1->setNodes(nodes);

        nodes[0] = hexnodes[ (i + offset + 3) % 6];
        nodes[1] = hexnodes[ (i + offset + 4) % 6];
        nodes[2] = hexnodes[ (i + offset + 5) % 6];
        nodes[3] = hexnodes[ (i + offset + 6) % 6];
        face2->setNodes(nodes);
        bool convex1 = JQuadGeometry::isConvex(face1);
        bool convex2 = JQuadGeometry::isConvex(face2);
        if (convex1 && convex2) {
            quads.resize(2);
            quads[0] = face1;
            quads[1] = face2;
            return 0;
        }
    }
    quads.clear();
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
double
JQuadGeometry::getArea(const Point3D &p0, const Point3D &p1,
                       const Point3D &p2, const Point3D &p3)
{
    /////////////////////////////////////////////////////////////////////////////
    // For explanation of some amazing proofs and theorems, please refer to the
    // following site. This implementation is based on this article.
    //
    // http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#Quadrilaterals
    //
    // For Bretschneider's Formula: refer to
    // http://mathworld.wolfram.com/BretschneidersFormula.html
    // Given a general quadrilateral with sides of lengths a, b, c, and d, the area
    // is given by
    //              K = 1/4sqrt(4p^2q^2-(b^2+d^2-a^2-c^2)^2)
    //
    // It seems that the first method is better because it doesn't require
    // expensive 6 lenghts. ( 4 sides + 2 diagonal ). But a proper analysis needs
    // be done. But probably, the most amazing thing about the formula is that
    // it is valid for 3D quadrilateral and for both convex and concave.
    // But why it handles concavity, I am not quite sure.
    //
    // Chaman Singh Verma
    // 16th Feb 2010.
    //////////////////////////////////////////////////////////////////////////////

    Vec3D d0d1, v2v0, v3v1;
    JMath::make_vector(p2, p0, v2v0);
    JMath::make_vector(p3, p1, v3v1);
    JMath::cross_product(v2v0, v3v1, d0d1);

    double area = 0.5 * JMath::magnitude(d0d1);
    return area;
}
////////////////////////////////////////////////////////////////////////////////
void
JQuadGeometry::bilinear_weights(double xi, double eta, vector<double> &weight)
{
    weight.resize(4);

#ifdef DEBUG
    assert(xi >= -1.0 && xi <= 1.0);
    assert(eta >= -1.0 && eta <= 1.0);
#endif

    double coeff = 1.0 / 4.0;

    weight[0] = coeff * (1.0 - xi)*(1.0 - eta);
    weight[1] = coeff * (1.0 + xi)*(1.0 - eta);
    weight[2] = coeff * (1.0 + xi)*(1.0 + eta);
    weight[3] = coeff * (1.0 - xi)*(1.0 + eta);
}

/////////////////////////////////////////////////////////////////////////////////////

double
JSimplex::linear_interpolation(const vector<double> &x, const vector<double> &w)
{
    size_t numNodes = x.size();
    assert(w.size() == numNodes);

    double sum = 0.0;
    for (size_t i = 0; i < numNodes; i++) sum += x[i] * w[i];

    return sum;
}

/////////////////////////////////////////////////////////////////////////////////////

int JFaceGeometry::reflexAngleAt( const JFacePtr &face)
{
    int nsize = face->getSize(0);

    if( nsize == 3) return -1;

    double x[5], y[5], z[5], triarea;
    JNodePtr vertex;

    Point3D xyz;

    for (int i = 0; i < nsize; i++) {
        vertex = face->getNodeAt(i);
        xyz = vertex->getXYZCoords();
        x[0] = xyz[0];
        y[0] = xyz[1];
        z[0] = xyz[2];

        vertex = face->getNodeAt(i + 1 );
        xyz = vertex->getXYZCoords();
        x[1] = xyz[0];
        y[1] = xyz[1];
        z[1] = xyz[2];

        vertex = face->getNodeAt(i + nsize - 1 );
        xyz = vertex->getXYZCoords();
        x[2] = xyz[0];
        y[2] = xyz[1];
        z[2] = xyz[2];

        triarea = PolygonArea3D(3, x, y, z);

        if (triarea < 0.0) return i;
    }

    return -1;
}

/////////////////////////////////////////////////////////////////////////////////////

bool
JFaceGeometry::isSimple( const JFacePtr &face)
{
    int nnodes =  face->getSize(0);
    if( nnodes == 3 ) return 1;

    double d0, d1, d2, d3;
    if( nnodes == 4 ) {
        const Point3D &v0 = face->getNodeAt(0)->getXYZCoords();
        const Point3D &v1 = face->getNodeAt(1)->getXYZCoords();
        const Point3D &v2 = face->getNodeAt(2)->getXYZCoords();
        const Point3D &v3 = face->getNodeAt(3)->getXYZCoords();

        d0  = JNodeGeometry::getOrientation( v0, v1, v2);
        d1  = JNodeGeometry::getOrientation( v0, v1, v3);

        d2  = JNodeGeometry::getOrientation( v2, v3, v0);
        d3  = JNodeGeometry::getOrientation( v2, v3, v1);

        if( (d0*d1 < 0.0) && (d2*d3 < 0.0 ) ) return 0;

        d0  = JNodeGeometry::getOrientation( v1, v2, v0);
        d1  = JNodeGeometry::getOrientation( v1, v2, v3);

        d2  = JNodeGeometry::getOrientation( v0, v3, v1);
        d3  = JNodeGeometry::getOrientation( v0, v3, v2);

        if( (d0*d1 < 0.0) && (d2*d3 < 0.0 ) ) return 0;

        return 1;
    }

    cout << "Warning: General Polygons not supported yet " << endl;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

double JFaceGeometry :: getSignedArea( const JFacePtr &face)
{
    int np = face->getSize(0);
    vector<double> x(np), y(np);
    for( int i = 0; i < np; i++) {
        const Point3D &p = face->getNodeAt(i)->getXYZCoords();
        x[i] = p[0];
        y[i] = p[1];
    }
    double area = JGeometry::getSignedArea( &x[0], &y[0], np);
    return area;
}
///////////////////////////////////////////////////////////////////////////////

double JFaceGeometry :: getArea( const JFacePtr &face)
{
    int nnodes = face->getSize(0);
    if (nnodes == 3) {
        return JTriGeometry::getArea(face->getNodeAt(0)->getXYZCoords(),
                                     face->getNodeAt(1)->getXYZCoords(),
                                     face->getNodeAt(2)->getXYZCoords());
    }

    if (nnodes == 4) {
        return JQuadGeometry::getArea(face->getNodeAt(0)->getXYZCoords(),
                                      face->getNodeAt(1)->getXYZCoords(),
                                      face->getNodeAt(2)->getXYZCoords(),
                                      face->getNodeAt(3)->getXYZCoords());
    }
    cout << "Redo it " << __LINE__ << endl;
    exit(0);
}

////////////////////////////////////////////////////////////////////////////////
double JFaceGeometry :: getArea( const JFaceSequence &seq)
{
    double sum = 0.0;
    for( const JFacePtr &f : seq) sum += getArea(f);
    return sum;
}
////////////////////////////////////////////////////////////////////////////////
void JFaceGeometry :: getAngles( const JFacePtr &face, vector<double> &angles, int measure)
{
    angles.clear();
    int nnodes = face->getSize(0);
    if( nnodes < 3) return;
    angles.resize(nnodes);
    for( int i = 0; i < nnodes; i++)
        angles[i] = JFaceGeometry::getAngleAt(face, face->getNodeAt(i), measure);
}

////////////////////////////////////////////////////////////////////////////////
double JFaceGeometry :: getMaxAngle( const JFacePtr &face, int measure)
{
    int nnodes = face->getSize(0);
    assert( nnodes > 2 );
    double maxAngle = 0.0;
    for( int i = 0; i < nnodes; i++)  {
        double angle = JFaceGeometry::getAngleAt(face, face->getNodeAt(i), measure);
        maxAngle  = max(angle, maxAngle);
    }
    return maxAngle;
}

////////////////////////////////////////////////////////////////////////////////

double JFaceGeometry :: getMinAngle( const JFacePtr &face, int measure)
{
    int nnodes = face->getSize(0);
    assert( nnodes > 2 );
    double minAngle = std::numeric_limits<double>::max();
    for( int i = 0; i < nnodes; i++)  {
        double angle = JFaceGeometry::getAngleAt(face, face->getNodeAt(i), measure);
        minAngle  = min(angle, minAngle);
    }
    return minAngle;
}

////////////////////////////////////////////////////////////////////////////////

void
JFaceGeometry::getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, Point3D &pc)
{
    Point3D xyz;
    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;

    xyz = v0->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v1->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v2->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    pc[0] /= 3.0;
    pc[1] /= 3.0;
    pc[2] /= 3.0;
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr
JFaceGeometry::getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2)
{
    Point3D xyz;
    JFaceGeometry::getCentroid(v0, v1, v2, xyz);

    JNodePtr vc = JNode::newObject();
    vc->setXYZCoords(xyz);
    return vc;
}
///////////////////////////////////////////////////////////////////////////////

void
JFaceGeometry::getCentroid(const JFacePtr &face, Point3D &pC)
{
    pC[0] = 0.0;
    pC[1] = 0.0;
    pC[2] = 0.0;
    int nSize = face->getSize(0);
    for( int i = 0; i < nSize; i++) {
        const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
        pC[0] += xyz[0];
        pC[1] += xyz[1];
        pC[2] += xyz[2];
    }
    pC[0] /= ( double)nSize;
    pC[1] /= ( double)nSize;
    pC[2] /= ( double)nSize;
}
///////////////////////////////////////////////////////////////////////////////

JNodePtr
JFaceGeometry::getCentroid(const JFacePtr &face)
{
    Point3D pC;
    getCentroid(face, pC);
    JNodePtr vc = JNode::newObject();
    vc->setXYZCoords(pC);
    return vc;
}
///////////////////////////////////////////////////////////////////////////////

void
JFaceGeometry::getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3, Point3D &pc)
{
    Point3D  xyz;

    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;

    xyz = v0->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v1->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v2->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v3->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    pc[0] /= 4.0;
    pc[1] /= 4.0;
    pc[2] /= 4.0;
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr
JFaceGeometry::getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3)
{
    Point3D xyz;
    getCentroid( v0, v1, v2, v3, xyz);
    JNodePtr vc = JNode::newObject();
    vc->setXYZCoords(xyz);
    return vc;
}

///////////////////////////////////////////////////////////////////////////////

void
JFaceGeometry::getCentroid(const JNodePtr &v0, const JNodePtr &v1,
                           const JNodePtr &v2, const JNodePtr &v3,
                           const JNodePtr &v4, Point3D &pc)
{
    Point3D  xyz;

    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;

    xyz = v0->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v1->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v2->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v3->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    xyz = v4->getXYZCoords();
    pc[0] += xyz[0];
    pc[1] += xyz[1];
    pc[2] += xyz[2];

    pc[0] /= 5.0;
    pc[1] /= 5.0;
    pc[2] /= 5.0;
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr
JFaceGeometry::getCentroid(const JNodePtr &v0, const JNodePtr &v1,
                           const JNodePtr &v2, const JNodePtr &v3,
                           const JNodePtr &v4)
{
    Point3D xyz;
    getCentroid( v0, v1, v2, v3, v4, xyz);
    JNodePtr vc = JNode::newObject();
    vc->setXYZCoords(xyz);
    return vc;
}

///////////////////////////////////////////////////////////////////////////////


int
JFaceGeometry::is_3_sided_convex_loop_quad_meshable(const int *segments, int *partsegments)
{
    double M[6][6], rhs[6];

    // octave:1> A = [ 0  0 -1  1  0  0;
    //                -1  0  0  0  1  0;
    //                 0 -1  0  0  0  1;
    //		       1  0  0  1  0  0;
    //                 0  1  0  0  1  0;
    //                 0  0  1  0  0  1]
    //A =
    //
    //   0   0  -1   1   0   0
    //  -1   0   0   0   1   0
    //   0  -1   0   0   0   1
    //   1   0   0   1   0   0
    //   0   1   0   0   1   0
    //   0   0   1   0   0   1
    //  octave:2> inv(A)
    //  ans =

    //  -0.5  -0.5   0.5   0.5   0.5  -0.5
    //   0.5  -0.5  -0.5  -0.5   0.5   0.5
    //  -0.5   0.5  -0.5   0.5  -0.5   0.5
    //   0.5   0.5  -0.5   0.5  -0.5   0.5
    //  -0.5   0.5   0.5   0.5   0.5  -0.5
    //   0.5  -0.5   0.5  -0.5   0.5   0.5

    //   First Row
    //  -0.5  -0.5   0.5   0.5   0.5  -0.5
    M[0][0] = -0.5;
    M[0][1] = -0.5;
    M[0][2] = 0.5;
    M[0][3] = 0.5;
    M[0][4] = 0.5;
    M[0][5] = -0.5;

    //   Second Row
    //   0.5  -0.5  -0.5  -0.5   0.5   0.5
    M[1][0] = 0.5;
    M[1][1] = -0.5;
    M[1][2] = -0.5;
    M[1][3] = -0.5;
    M[1][4] = 0.5;
    M[1][5] = 0.5;

    //  Third Row
    //  -0.5   0.5  -0.5   0.5  -0.5   0.5
    M[2][0] = -0.5;
    M[2][1] = 0.5;
    M[2][2] = -0.5;
    M[2][3] = 0.5;
    M[2][4] = -0.5;
    M[2][5] = 0.5;

    // Forth Row
    //   0.5   0.5  -0.5   0.5  -0.5   0.5
    M[3][0] = 0.5;
    M[3][1] = 0.5;
    M[3][2] = -0.5;
    M[3][3] = 0.5;
    M[3][4] = -0.5;
    M[3][5] = 0.5;

    //   Fifth Row
    //  -0.5   0.5   0.5   0.5   0.5  -0.5
    M[4][0] = -0.5;
    M[4][1] = 0.5;
    M[4][2] = 0.5;
    M[4][3] = 0.5;
    M[4][4] = 0.5;
    M[4][5] = -0.5;

    //  Sixth Row
    //   0.5  -0.5   0.5  -0.5   0.5   0.5
    M[5][0] = 0.5;
    M[5][1] = -0.5;
    M[5][2] = 0.5;
    M[5][3] = -0.5;
    M[5][4] = 0.5;
    M[5][5] = 0.5;

    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = 0.0;

    rhs[3] = segments[0];
    rhs[4] = segments[1];
    rhs[5] = segments[2];

    vector<int> x(6);
    for (int i = 0; i < 6; i++) {
        double sum = 0.0;
        for (int j = 0; j < 6; j++)
            sum += M[i][j] * rhs[j];
        x[i] = (int) sum;
    }

    for (int i = 0; i < 6; i++)
        if (x[i] < 1) return 0;

    if (x[0] + x[3] != rhs[3]) return 0;
    partsegments[0] = x[0];
    partsegments[1] = x[3];

    if (x[1] + x[4] != rhs[4]) return 0;
    partsegments[2] = x[1];
    partsegments[3] = x[4];

    if (x[2] + x[5] != rhs[5]) return 0;
    partsegments[4] = x[2];
    partsegments[5] = x[5];

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int
JFaceGeometry::is_5_sided_convex_loop_quad_meshable(const int *segments, int *partSegments)
{
    double M[10][10], rhs[10];

    //  Equations:
    //      b0 -a2   = 0
    //      b1 -a3   = 0
    //      b2 -a4   = 0
    //      b3 -a0   = 0
    //      b4 -a1   = 0
    //      a0 + b0  = s0
    //      a1 + b1  = s1
    //      a2 + b2  = s2
    //      a3 + b3  = s3
    //      a4 + b4  = s4

    // For more details; See Guy Bunin's paper.
    // M =
    //   a0  a1 a2  a3  a4  b0   b1 b2   b3  b4
    //   0   0  -1   0   0   1   0   0   0   0
    //   0   0   0  -1   0   0   1   0   0   0
    //   0   0   0   0  -1   0   0   1   0   0
    //  -1   0   0   0   0   0   0   0   1   0
    //   0  -1   0   0   0   0   0   0   0   1
    //   1   0   0   0   0   1   0   0   0   0
    //   0   1   0   0   0   0   1   0   0   0
    //   0   0   1   0   0   0   0   1   0   0
    //   0   0   0   1   0   0   0   0   1   0
    //   0   0   0   0   1   0   0   0   0   1

    // octave:38> inv(M)
    // ans =
    //  -0.5   0.5   0.5  -0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    //  -0.5  -0.5   0.5   0.5  -0.5   0.5   0.5  -0.5  -0.5   0.5
    //  -0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    //   0.5  -0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    //   0.5   0.5  -0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    //   0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    //   0.5   0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    //   0.5   0.5   0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    //  -0.5   0.5   0.5   0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    //  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5   0.5


    //   First Row
    //  -0.5   0.5   0.5  -0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    M[0][0] = -0.5;
    M[0][1] = 0.5;
    M[0][2] = 0.5;
    M[0][3] = -0.5;
    M[0][4] = -0.5;
    M[0][5] = 0.5;
    M[0][6] = -0.5;
    M[0][7] = -0.5;
    M[0][8] = 0.5;
    M[0][9] = 0.5;

    //  -0.5  -0.5   0.5   0.5  -0.5   0.5   0.5  -0.5  -0.5   0.5
    M[1][0] = -0.5;
    M[1][1] = -0.5;
    M[1][2] = 0.5;
    M[1][3] = 0.5;
    M[1][4] = -0.5;
    M[1][5] = 0.5;
    M[1][6] = 0.5;
    M[1][7] = -0.5;
    M[1][8] = -0.5;
    M[1][9] = 0.5;


    //  -0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    M[2][0] = -0.5;
    M[2][1] = -0.5;
    M[2][2] = -0.5;
    M[2][3] = 0.5;
    M[2][4] = 0.5;
    M[2][5] = 0.5;
    M[2][6] = 0.5;
    M[2][7] = 0.5;
    M[2][8] = -0.5;
    M[2][9] = -0.5;

    //   0.5  -0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    M[3][0] = 0.5;
    M[3][1] = -0.5;
    M[3][2] = -0.5;
    M[3][3] = -0.5;
    M[3][4] = 0.5;
    M[3][5] = -0.5;
    M[3][6] = 0.5;
    M[3][7] = 0.5;
    M[3][8] = 0.5;
    M[3][9] = -0.5;

    //   0.5   0.5  -0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    M[4][0] = 0.5;
    M[4][1] = 0.5;
    M[4][2] = -0.5;
    M[4][3] = -0.5;
    M[4][4] = -0.5;
    M[4][5] = -0.5;
    M[4][6] = -0.5;
    M[4][7] = 0.5;
    M[4][8] = 0.5;
    M[4][9] = 0.5;

    //   0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    M[5][0] = 0.5;
    M[5][1] = -0.5;
    M[5][2] = -0.5;
    M[5][3] = 0.5;
    M[5][4] = 0.5;
    M[5][5] = 0.5;
    M[5][6] = 0.5;
    M[5][7] = 0.5;
    M[5][8] = -0.5;
    M[5][9] = -0.5;

    //   0.5   0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    M[6][0] = 0.5;
    M[6][1] = 0.5;
    M[6][2] = -0.5;
    M[6][3] = -0.5;
    M[6][4] = 0.5;
    M[6][5] = -0.5;
    M[6][6] = 0.5;
    M[6][7] = 0.5;
    M[6][8] = 0.5;
    M[6][9] = -0.5;

    //   0.5   0.5   0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    M[7][0] = 0.5;
    M[7][1] = 0.5;
    M[7][2] = 0.5;
    M[7][3] = -0.5;
    M[7][4] = -0.5;
    M[7][5] = -0.5;
    M[7][6] = -0.5;
    M[7][7] = 0.5;
    M[7][8] = 0.5;
    M[7][9] = 0.5;

    //  -0.5   0.5   0.5   0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    M[8][0] = -0.5;
    M[8][1] = 0.5;
    M[8][2] = 0.5;
    M[8][3] = 0.5;
    M[8][4] = -0.5;
    M[8][5] = 0.5;
    M[8][6] = -0.5;
    M[8][7] = -0.5;
    M[8][8] = 0.5;
    M[8][9] = 0.5;

    //  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5   0.5
    M[9][0] = -0.5;
    M[9][1] = -0.5;
    M[9][2] = 0.5;
    M[9][3] = 0.5;
    M[9][4] = 0.5;
    M[9][5] = 0.5;
    M[9][6] = 0.5;
    M[9][7] = -0.5;
    M[9][8] = -0.5;
    M[9][9] = 0.5;

    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = 0.0;
    rhs[3] = 0.0;
    rhs[4] = 0.0;

    rhs[5] = segments[0];
    rhs[6] = segments[1];
    rhs[7] = segments[2];
    rhs[8] = segments[3];
    rhs[9] = segments[4];

    vector<int> x(10);
    for (int i = 0; i < 10; i++) {
        double sum = 0.0;
        for (int j = 0; j < 10; j++)
            sum += M[i][j] * rhs[j];
        x[i] = (int) sum;
    }

    for (int i = 0; i < 10; i++)
        if (x[i] < 1) return 0;

    if (x[0] != x[8]) return 0;
    if (x[0] + x[5] != rhs[5]) return 0;
    partSegments[0] = x[0];
    partSegments[1] = x[5];

    if (x[1] != x[9]) return 0;
    if (x[1] + x[6] != rhs[6]) return 0;
    partSegments[2] = x[1];
    partSegments[3] = x[6];

    if (x[2] != x[5]) return 0;
    if (x[2] + x[7] != rhs[7]) return 0;
    partSegments[4] = x[2];
    partSegments[5] = x[7];

    if (x[3] != x[6]) return 0;
    if (x[3] + x[8] != rhs[8]) return 0;
    partSegments[6] = x[3];
    partSegments[7] = x[8];

    if (x[4] != x[7]) return 0;
    if (x[4] + x[9] != rhs[9]) return 0;
    partSegments[8] = x[4];
    partSegments[9] = x[9];

    return 1;
}

JEdgePtr JQuadrilateral :: getOppositeEdge( const JFacePtr &face, const JEdgePtr &e)
{
    assert( face->getSize(0) == 4);

    JNodePtr v0 = e->getNodeAt(0);
    JNodePtr v1 = e->getNodeAt(1);

    int pos0   = face->getPosOf( v0 );
    int pos1   = face->getPosOf( v1 );

    JEdgePtr oedge;

    if( (pos0 == 0 && pos1 == 1 )  || ( pos0 == 1 && pos1 == 0) ) {
        v0 = face->getNodeAt(2);
        v1 = face->getNodeAt(3);
        oedge = face->getEdgeOf(v0,v1);
        return oedge;
    }

    if( (pos0 == 1 && pos1 == 2 )  || ( pos0 == 2 && pos1 == 1) ) {
        v0 = face->getNodeAt(3);
        v1 = face->getNodeAt(0);
        oedge = face->getEdgeOf(v0,v1);
        return oedge;
    }

    if( (pos0 == 2 && pos1 == 3 )  || ( pos0 == 3 && pos1 == 2) ) {
        v0 = face->getNodeAt(0);
        v1 = face->getNodeAt(1);
        oedge = face->getEdgeOf(v0,v1);
        return oedge;
    }

    if( (pos0 == 0 && pos1 == 3 )  || ( pos0 == 3 && pos1 == 0) ) {
        v0 = face->getNodeAt(1);
        v1 = face->getNodeAt(2);
        oedge = face->getEdgeOf(v0,v1);
        return oedge;
    }
    cout << "Warning: No edge present" << endl;
    return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////

JEdgeSequence JPolygon :: getCircle( const Point2D &center,  double radius, int np)
{
    double dtheta = 2.0*M_PI/(double)np;

    JNodeSequence nodes(np);
    Point3D xyz;
    for( int i = 0; i < np; i++) {
        xyz[0] = center[0] + radius*cos(i*dtheta);
        xyz[1] = center[0] + radius*sin(i*dtheta);
        xyz[2] = 0.0;
        JNodePtr vtx = JNode::newObject();
        vtx->setXYZCoords(xyz);
        nodes[i] = vtx;
    }

    JEdgeSequence edges(np);
    for( int i = 0; i < np; i++) {
        JNodePtr v0 = nodes[i];
        JNodePtr v1 = nodes[(i+1)%np];
        edges[i]   = JEdge::newObject(v0,v1);
    }
    return edges;
}


/*
double JFaceGeometry :: getArea2D(const double *x, const double *y, int n)
{
    // Given "n" points ( 2D ) on a polygon, (P(n+1) = p(0)) this function
    // return the signed area..
    assert( n >= 3);

    double sum = 0.0;
    for( int i  = 0; i < n; i++) {
        double x0 = x[i];
        double y0 = y[i];
        double x1 = x[(i+1)%n];
        double y1 = y[(i+1)%n];
        sum += x0*y1 - x1*y0;
    }
    return 0.5*sum;
}
*/

/*
double JFaceGeometry ::  getArea2D( const vector<Point2D> &poly)
{
    int np = poly.size();
    vector<double> x(np), y(np);
    for( int i = 0; i < np; i++) {
        x[i] = poly[i][0];
        y[i] = poly[i][1];
    }

    return getArea2D( &x[0], &y[0], np);
}
*/

/*
void JFaceGeometry :: getCentroid( const double *x, const double *y, int n, double *center )
{
    assert( n >= 3);

    // Formula from wikipedia. Note that centroid does not depend on the
    // orientation of the polygon. The division by Area will take care
    // of correct value.

    // For convex bodies, centroid is always inside the region.

    double cx = 0.0;
    double cy = 0.0;
    double cf;

    for( int i  = 0; i < n; i++) {
        cf  = x[i]*y[(i+1)%n] - x[(i+1)%n]*y[i];
        cx +=  (x[i]+x[(i+1)%n])*cf;
        cy +=  (y[i]+y[(i+1)%n])*cf;
    }

    double A = getArea2D(x, y, n);

    center[0] = cx/(6.0*A);
    center[1] = cy/(6.0*A);
}
*/

/*
int
JFace::get_topological_neighbors( JNodeSequence &vneighs)
{
    vneighs.clear();
    JNodeSet vset;

    int nnodes = nodes.size();

    for (int i = 0; i < nnodes; i++) {
        JNodePtr vertex = nodes[i];
        vertex->getRelations( vneighs );
        int numneighs = vneighs.size();
        for (int j = 0; j < numneighs; j++)
            vset.insert(vneighs[j]);
    }

    for (int i = 0; i < nnodes; i++)
        vset.erase(nodes[i]);

    JNodeSequence vresult;
    if (!vset.empty()) {
        JNodeSet::const_iterator it;
        size_t index = 0;
        for (it = vset.begin(); it != vset.end(); ++it)
            vresult[index++] = *it;
    }
    return 0;
}
*/


/*
int
Face::get_interior_angles( vector<double> &vals ) const
{
     vals.clear();

     JNodeSequence rotatedNodes;
     map<Vertex*, double> mapangles;

     int nsize = nodes.size();
     for (int i = 0; i < nsize; i++)
          mapangles[ nodes[i] ] = 0.0;
     Quadrilateral::quad_tessalate(nodes, rotatedNodes);

     Array3D tangles;

     const Point3D &p0 = rotatedNodes[0]->getXYZCoords();
     const Point3D &p1 = rotatedNodes[1]->getXYZCoords();
     const Point3D &p2 = rotatedNodes[2]->getXYZCoords();
     const Point3D &p3 = rotatedNodes[3]->getXYZCoords();

     JMath::getTriAngles(p0, p1, p2, tangles);
     mapangles[ rotatedNodes[0]] += tangles[0];
     mapangles[ rotatedNodes[1]] += tangles[1];
     mapangles[ rotatedNodes[2]] += tangles[2];

     JMath::getTriAngles(p0, p2, p3, tangles);
     mapangles[ rotatedNodes[0]] += tangles[0];
     mapangles[ rotatedNodes[2]] += tangles[1];
     mapangles[ rotatedNodes[3]] += tangles[2];

     vals.resize(nsize);
     for (int i = 0; i < nsize; i++)
          vals[i] = mapangles[ nodes[i] ];
     return 0;
}
*/
bool JFaceGeometry::intersectPredicate2d(const JFacePtr &face1, const JFacePtr &face2)
{
    /*
        if( face1 == face2) return 0;
        // Cheapest check: Bounding Boxes must intersect...
        const BoundingBox &box1 = face1->getBoundingBox();
        const BoundingBox &box2 = face2->getBoundingBox();
        if( BoundingBox::intersect(box1,box2) == 0) return 0;

        // Check if the edges are intersecting ....
        int numEdges1 = face1->getSize(1);
        int numEdges2 = face2->getSize(1);

        for(int i = 0; i < numEdges1; i++) {
            const JEdgePtr &ei = face1->getEdgeAt(i);
            for( int j = 0; j < numEdges2; j++) {
                const JEdgePtr &ej = face2->getEdgeAt(j);
                if( JEdgeGeometry::intersectPredicate2d(ei, ej) ) return 1;
            }
        }

        // Check if one face is totally inside the other ...
        int nCount, numnodes;

        numnodes = face2->getSize(0);
        nCount = 0;
        for( int i = 0; i < numnodes; i++) {
            int ori = JFaceGeometry::getBoundedSide( face1, face2->getNodeAt(i)->getXYZCoords() );
            if( ori == GeomOrient::INSIDE) nCount++;
        }
        if( nCount == numnodes) return 1;

        numnodes = face1->getSize(0);
        nCount = 0;
        for( int i = 0; i < numnodes; i++) {
            int ori = JFaceGeometry::getBoundedSide( face2, face1->getNodeAt(i)->getXYZCoords() );
            if( ori == GeomOrient::INSIDE) nCount++;
        }
        if( nCount == numnodes) return 1;
    */
    cout << "REIMPLEMENT " << __LINE__ << endl;
    exit(0);
    return 0;
}

