#include <iomanip>

#include "Mesh.hpp"
#include "basic_math.hpp"

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

size_t JNode::GlobalID = 0;
size_t JNode::NumObjectsCreated = 0;
std::map<string,string> JNode::attribInfo;

void JMeshEntity::init_random_number()
{
    cout << "Rnadom number" << endl;

    srand48( time(0));
}

//EntityFeatureSet* Node :: featureSet = new EntityFeatureSet;

///////////////////////////////////////////////////////////////////////////////
int JNode :: registerAttribute( const string &name, const string &type)
{
    int  found = 0;

    if( type == "int"    ) found = 1;
    if( type == "char"   ) found = 1;
    if( type == "float"  ) found = 1;
    if( type == "double" ) found = 1;
    if( type == "uchar"  ) found = 1;
    if( !found) {
        cout << "Warning: invalid attribute type " << endl;
        return 2;
    }
    attribInfo[name] = type;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

string JNode :: getAttributeTypeName( const string &name)
{
    string str;
    if( attribInfo.find( name ) == attribInfo.end() ) return str;
    return  attribInfo[name];
}
///////////////////////////////////////////////////////////////////////////////

JNodeSequence JNode :: newObjects( size_t n)
{
    assert(n > 0);
    JNodeSequence newvec(n);
    for( size_t i = 0; i < n; i++)
        newvec[i] = JNode::newObject();
    return newvec;
}

///////////////////////////////////////////////////////////////////////////////

void JNode :: getRelations( const JNodePtr &vtx, JNodeSequence &seq)
{
    seq.clear();
    if( vtx == nullptr) return;
    vtx->getRelations_(seq);
}

void JNode :: getRelations( const JNodePtr &vtx, JEdgeSequence &seq)
{
    seq.clear();
    if( vtx == nullptr) return;
    vtx->getRelations_(seq);
}

void JNode :: getRelations( const JNodePtr &vtx, JFaceSequence &seq)
{
    seq.clear();
    if( vtx == nullptr) return;
    vtx->getRelations_(seq);
}

void JNode :: getRelations( const JNodePtr &vtx, JCellSequence &seq)
{
    seq.clear();
    if( vtx == nullptr) return;
    vtx->getRelations_(seq);
}

int JNode::get_ideal_face_degree(int n) const
{
    /*
        if( n == 3 ) {
            if (!isBoundary()) return 6;

            if (getSpanAngle() <= 90.0 ) return 1;
            if (getSpanAngle() <= 220.0) return 3;
            if (getSpanAngle() <= 300.0) return 4;
            return 5;
        }

        if( n == 4 ) {
            if (!isBoundary()) return 4;

            if (getSpanAngle() <= 100 )  return 1;
            if (getSpanAngle() <= 220.0) return 2;
            if (getSpanAngle() <= 300.0) return 3;
            return 4;
        }
    */

    cout << "Error: Ideal vertex degree only for Quad right now " << endl;
    exit(0);
}
///////////////////////////////////////////////////////////////////////////////

bool JNode :: isInternal() const
{
    if( this->isBoundary() ) return 0;
    if( this->hasAttribute("Interface"))  return 0;
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
Point3D
JNodeGeometry ::getMidPoint(const JNodePtr &v0, const JNodePtr &v1, double alpha)
{

    const Point3D &p0 = v0->getXYZCoords();
    const Point3D &p1 = v1->getXYZCoords();

    Point3D pmid;
    pmid[0] = (1 - alpha) * p0[0] + alpha * p1[0];
    pmid[1] = (1 - alpha) * p0[1] + alpha * p1[1];
    pmid[2] = (1 - alpha) * p0[2] + alpha * p1[2];
    return pmid;
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr
JNodeGeometry::getMidNode(const JNodePtr &v0, const JNodePtr &v1, double alpha)
{
    Point3D xyz = JNodeGeometry::getMidPoint(v0, v1, alpha);

    JNodePtr vmid = JNode::newObject();
    vmid->setXYZCoords(xyz);
    return vmid;
}

///////////////////////////////////////////////////////////////////////////////

double
JNodeGeometry::getLength(const JNodePtr &v0, const JNodePtr &v1)
{
    const Point3D &p0 = v0->getXYZCoords();
    const Point3D &p1 = v1->getXYZCoords();

    double dx = p0[0] - p1[0];
    double dy = p0[1] - p1[1];
    double dz = p0[2] - p1[2];

    return sqrt(dx * dx + dy * dy + dz * dz);
}

///////////////////////////////////////////////////////////////////////////////

double
JNodeGeometry::getLength2(const JNodePtr &v0, const JNodePtr &v1)
{
    const Point3D &p0 = v0->getXYZCoords();
    const Point3D &p1 = v1->getXYZCoords();

    double dx = p0[0] - p1[0];
    double dy = p0[1] - p1[1];
    double dz = p0[2] - p1[2];

    return dx * dx + dy * dy + dz * dz;
}

///////////////////////////////////////////////////////////////////////////////

double JNodeGeometry:: getOrientation( const Point3D &p0, const Point3D &p1, const Point3D &qpoint)
{
    double x[5], y[5], z[5];

    x[0] = p0[0];
    x[1] = p1[0];
    x[2] = qpoint[0];

    y[0] = p0[1];
    y[1] = p1[1];
    y[2] = qpoint[1];

    z[0] = p0[2];
    z[1] = p1[2];
    z[2] = qpoint[2];

    return PolygonArea3D(3, x, y, z);
}
/////////////////////////////////////////////////////////////////////////////////////
JBoundingBox JNodeGeometry :: getBoundingBox( const JNodeSequence &nodes)
{
    JBoundingBox box;
    if( nodes.empty() ) return box;
    Point3D xyz = nodes[0]->getXYZCoords();
    double xmin = xyz[0];
    double ymin = xyz[1];
    double zmin = xyz[2];
    double xmax = xmin;
    double ymax = ymin;
    double zmax = ymin;

    for( const JNodePtr &vtx : nodes) {
        const Point3D &p3d = vtx->getXYZCoords();
        xmin  = min( xmin, p3d[0] );
        ymin  = min( ymin, p3d[1] );
        zmin  = min( zmin, p3d[2] );

        xmax  = max( xmax, p3d[0] );
        ymax  = max( ymax, p3d[1] );
        zmax  = max( zmax, p3d[2] );
    }
    Point3D pmin, pmax;
    pmin[0] = xmin;
    pmin[1] = ymin;
    pmin[2] = zmin;

    pmax[0] = xmax;
    pmax[1] = ymax;
    pmax[2] = zmax;

    box.setPoints( pmin, pmax);
    return box;
}

/////////////////////////////////////////////////////////////////////////////////////
double JNodeGeometry :: getSpanAngleAt( const JNodePtr &vtx, int measure)
{
    JFaceSequence vfaces;
    JNode::getRelations( vtx, vfaces );

    double sum = 0.0;
    for( const JFacePtr &face: vfaces)  {
        sum += JFaceGeometry::getAngleAt(face, vtx, measure);
    }
    return sum;
}
/////////////////////////////////////////////////////////////////////////////////////

/*
bool Node::isFeature() const
{
    set<string>::const_iterator it;
    for( it = featureSet->attributes.begin(); it != featureSet->attributes.end(); ++it)
	if( this->hasAttribute(*it) ) return 1;
    return 0;
}
*/


