#include <iomanip>

#include "Mesh.hpp"
#include "basic_math.hpp"

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int
JSimplex ::getAvgXYZ( Point3D &pc) const
{
    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;

    int nsize = nodes.size();
    for (int inode = 0; inode < nsize; inode++) {
        JNodePtr v = nodes[inode];
        const Point3D &p3d = v->getXYZCoords();
        pc[0] += p3d[0];
        pc[1] += p3d[1];
        pc[2] += p3d[2];
    }

    double multby = 1.0/(double)nsize;
    pc[0] *=  multby;
    pc[1] *=  multby;
    pc[2] *=  multby;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int
JSimplex ::getAvgUV( Point2D &uv) const
{
    Point2D  p2d;

    uv[0]  = 0.0;
    uv[1]  = 0.0;
    p2d[0] = 0.0;
    p2d[1] = 0.0;
    int nsize = nodes.size();
    for (int inode = 0; inode < nsize; inode++) {
        JNodePtr v = nodes[inode];
        if( !v->hasAttribute("UVCoords")) return 1;
        v->getAttribute("UVCoords", p2d);
        uv[0] += p2d[0];
        uv[1] += p2d[1];
    }

    double multby = 1.0/(double)nsize;
    uv[0] *=  multby;
    uv[1] *=  multby;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

JBoundingBox
JSimplex::getBoundingBox() const
{
    JBoundingBox box;

    size_t numnodes = nodes.size();
    if( numnodes < 2) return box;

    double xmin, xmax, ymin, ymax, zmin, zmax;
    Point3D xyz;
    xyz = nodes[0]->getXYZCoords();

    xmin = xyz[0];
    xmax = xyz[0];
    ymin = xyz[1];
    ymax = xyz[1];
    zmin = xyz[2];
    zmax = xyz[2];

    for (size_t i = 0; i < numnodes; i++) {
        JNodePtr v = nodes[i];
        if( v->isActive() ) {
            xyz = nodes[i]->getXYZCoords();
            xmin = min(xmin, xyz[0]);
            xmax = max(xmax, xyz[0]);
            ymin = min(ymin, xyz[1]);
            ymax = max(ymax, xyz[1]);
            zmin = min(zmin, xyz[2]);
            zmax = max(zmax, xyz[2]);
        }
    }

    xyz[0] = xmin;
    xyz[1] = ymin;
    xyz[2] = zmin;
    box.setLower(xyz);

    xyz[0] = xmax;
    xyz[1] = ymax;
    xyz[2] = zmax;
    box.setUpper(xyz);

    return box;
}

/////////////////////////////////////////////////////////////////////////////////
bool JSimplex :: hasUniqueNodes() const
{
    int nn = nodes.size();
    for( int i = 0; i < nn; i++)
        for( int j = i+1; j < nn; j++)
            if( nodes[i] == nodes[j] ) return 0;
    return 1;
}

