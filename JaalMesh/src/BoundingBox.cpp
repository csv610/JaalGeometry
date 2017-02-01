#include "BoundingBox.hpp"
#include "GeomPredicates.hpp"

using namespace Jaal;

Point3D JBoundingBox :: getCorner( int ic ) const
{
    double xmin = lowerCorner[0];
    double ymin = lowerCorner[1];
    double zmin = lowerCorner[2];

    double xmax = upperCorner[0];
    double ymax = upperCorner[1];
    double zmax = upperCorner[2];

    Point3D p3d;

    switch(ic) {
    case 0:
        p3d[0] = xmin;
        p3d[1] = ymin;
        p3d[2] = zmin;
        break;
    case 1:
        p3d[0] = xmax;
        p3d[1] = ymin;
        p3d[2] = zmin;
        break;
    case 2:
        p3d[0] = xmax;
        p3d[1] = ymax;
        p3d[2] = zmin;
        break;
    case 3:
        p3d[0] = xmin;
        p3d[1] = ymax;
        p3d[2] = zmin;
        break;
    case 4:
        p3d[0] = xmin;
        p3d[1] = ymin;
        p3d[2] = zmax;
        break;
    case 5:
        p3d[0] = xmax;
        p3d[1] = ymin;
        p3d[2] = zmax;
        break;
    case 6:
        p3d[0] = xmax;
        p3d[1] = ymax;
        p3d[2] = zmax;
        break;
    case 7:
        p3d[0] = xmin;
        p3d[1] = ymax;
        p3d[2] = zmax;
    }
    return p3d;
}
///////////////////////////////////////////////////////////////////////////////

void JBoundingBox :: expandBy(double l)
{
    double dl = l*getLength(0);
    lowerCorner[0] -= 0.5*dl;
    upperCorner[0] += 0.5*dl;

    dl = l*getLength(1);
    lowerCorner[1] -= 0.5*dl;
    upperCorner[1] += 0.5*dl;

    dl = l*getLength(2);
    lowerCorner[2] -= 0.5*dl;
    upperCorner[2] += 0.5*dl;
}
///////////////////////////////////////////////////////////////////////////////


void JBoundingBox :: getEdges(vector<Point3D> &edges) const
{
    double xmin = lowerCorner[0];
    double ymin = lowerCorner[1];
    double zmin = lowerCorner[2];

    double xmax = upperCorner[0];
    double ymax = upperCorner[1];
    double zmax = upperCorner[2];

    vector<Point3D> corners(8);

    Point3D p3d;

    p3d[0] = xmin;
    p3d[1] = ymin;
    p3d[2] = zmin;
    corners[0] = p3d;

    p3d[0] = xmax;
    p3d[1] = ymin;
    p3d[2] = zmin;
    corners[1] = p3d;

    p3d[0] = xmax;
    p3d[1] = ymax;
    p3d[2] = zmin;
    corners[2] = p3d;

    p3d[0] = xmin;
    p3d[1] = ymax;
    p3d[2] = zmin;
    corners[3] = p3d;

    p3d[0] = xmin;
    p3d[1] = ymin;
    p3d[2] = zmax;
    corners[4] = p3d;

    p3d[0] = xmax;
    p3d[1] = ymin;
    p3d[2] = zmax;
    corners[5] = p3d;

    p3d[0] = xmax;
    p3d[1] = ymax;
    p3d[2] = zmax;
    corners[6] = p3d;

    p3d[0] = xmin;
    p3d[1] = ymax;
    p3d[2] = zmax;
    corners[7] = p3d;

    edges.resize(24);

    edges[0] =  corners[0];
    edges[1] =  corners[1];

    edges[2] =  corners[1];
    edges[3] =  corners[2];

    edges[4] =  corners[2];
    edges[5] =  corners[3];

    edges[6] =  corners[3];
    edges[7] =  corners[0];

    edges[8] =  corners[0];
    edges[9] =  corners[4];

    edges[10] =  corners[1];
    edges[11] =  corners[5];

    edges[12] =  corners[3];
    edges[13] =  corners[7];

    edges[14] =  corners[2];
    edges[15] =  corners[6];

    edges[16] =  corners[4];
    edges[17] =  corners[5];

    edges[18] =  corners[5];
    edges[19] =  corners[6];

    edges[20] =  corners[6];
    edges[21] =  corners[7];

    edges[22] =  corners[7];
    edges[23] =  corners[4];
}

///////////////////////////////////////////////////////////////////////////////

int JBoundingBox :: getOrientation( const JBoundingBox &box, const Point3D &p)
{
    int   dim = box.getDimension();
    const Point3D &lowerCorner = box.getLower();
    const Point3D &upperCorner = box.getUpper();

    if( p[0] < lowerCorner[0] ) return JGeomPredicates::OUTSIDE;
    if( p[0] > upperCorner[0] ) return JGeomPredicates::OUTSIDE;

    if( p[1] < lowerCorner[1] ) return JGeomPredicates::OUTSIDE;
    if( p[1] > upperCorner[1] ) return JGeomPredicates::OUTSIDE;

    if( p[2] < lowerCorner[2] ) return JGeomPredicates::OUTSIDE;
    if( p[2] > upperCorner[2] ) return JGeomPredicates::OUTSIDE;

    if( dim == 2) {
        if( p[0] > lowerCorner[0] && p[0] < upperCorner[0] &&
                p[1] > lowerCorner[1] && p[0] < upperCorner[1] ) return JGeomPredicates::INSIDE;
    }

    if( dim == 3) {
        if( p[0] > lowerCorner[0] && p[0] < upperCorner[0] &&
                p[1] > lowerCorner[1] && p[0] < upperCorner[1] &&
                p[2] > lowerCorner[2] && p[0] < upperCorner[2] ) return JGeomPredicates::INSIDE;
    }

    return JGeomPredicates::BOUNDARY;
}

///////////////////////////////////////////////////////////////////////////////

bool JBoundingBox :: intersect( const JBoundingBox &box1, const JBoundingBox &box2)
{
    const Point3D &p0 = box1.getLower();
    const Point3D &p1 = box1.getUpper();
    const Point3D &p2 = box2.getLower();
    const Point3D &p3 = box2.getUpper();

    for( int i = 0; i < 3; i++) {
        if( p2[i] > p1[i] ) return 0;
        if( p3[i] < p0[i] ) return 0;
    }

    return 1;
}
///////////////////////////////////////////////////////////////////////////////
void JBoundingBox :: setUnion(const JBoundingBox &box)
{
    Point3D xyz;

    xyz = box.getCorner(0);
    lowerCorner[0] = min( lowerCorner[0], xyz[0] );
    lowerCorner[1] = min( lowerCorner[1], xyz[1] );
    lowerCorner[2] = min( lowerCorner[2], xyz[2] );

    xyz = box.getCorner(6);
    upperCorner[0] = max( upperCorner[0], xyz[0] );
    upperCorner[1] = max( upperCorner[1], xyz[1] );
    upperCorner[2] = max( upperCorner[2], xyz[2] );
}

///////////////////////////////////////////////////////////////////////////////

/*
ostream& operator << (ostream &os, const Jaal::BoundingBox &box)
{
     Point3D p3d;

     p3d = box.getCorner(0);
     os << "{" << p3d[0] << " " << p3d[1] << " " << p3d[2] << "},";
     p3d = box.getCorner(6);
     os << "{" << p3d[0] << " " << p3d[1] << " " << p3d[2] << "}";
     os << endl;

     return os;
}
*/
///////////////////////////////////////////////////////////////////////////////
