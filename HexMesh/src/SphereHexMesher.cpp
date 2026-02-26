#include "SphereHexMesher.hpp"

///////////////////////////////////////////////////////////////////////

void JSphereHexMesher :: setMesh(const JMeshPtr &m)
{
    inMesh = m;
}

///////////////////////////////////////////////////////////////////////

void JSphereHexMesher :: setSphere(const JMeshPtr &m)
{
    sphMesh = m;
}

///////////////////////////////////////////////////////////////////////

void JSphereHexMesher :: setCubeCells( int *cellDim)
{
    int nodeDim[3];
    nodeDim[0] = max(2, cellDim[0] + 1);
    nodeDim[1] = max(2, cellDim[1] + 1);
    nodeDim[2] = max(2, cellDim[2] + 1);

    double len[3];
    len[0] = 2.0;
    len[1] = 2.0;
    len[2] = 2.0;

    double org[3];
    org[0] = -1.0;
    org[1] = -1.0;
    org[2] = -1.0;
    double radius = sqrt(2.0);

    if( sphMesh ) {
        JBoundingBox box = sphMesh->getGeometry()->getBoundingBox();
        len[0] = box.getLength(0);
        len[1] = box.getLength(1);
        len[2] = box.getLength(2);
        Point3D pC = box.getCorner(0);
        org[0] = pC[0];
        org[1] = pC[1];
        org[2] = pC[2];
        radius   =  0.5*len[0];

    }
    refHexMesh = AllHexMeshGenerator::getStructuredMesh(nodeDim, len, org);

    size_t numnodes = refHexMesh->getSize(0);

    Point3D xyz;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = refHexMesh->getNodeAt(i);
        if( vtx->isBoundary() ) {
            xyz = vtx->getXYZCoords();
            double len = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
            double nx   = xyz[0]/len;
            double ny   = xyz[1]/len;
            double nz   = xyz[2]/len;
            xyz[0]      = radius*nx;
            xyz[1]      = radius*ny;
            xyz[2]      = radius*nz;
            vtx->setXYZCoords(xyz);
        }
    }
}
