#include "MeshNormalClusters.hpp"

///////////////////////////////////////////////////////////////////////////////////

void JMeshNormalClusters :: buildSphere()
{
    JHexahedronPtr  hex = JHexahedron::getCanonical();

    sphMesh = JMesh::newObject();
    sphMesh->addObjects( hex->getNodes() );
    sphMesh->addObjects( hex->getFaces() );

    JQuadRefiner refiner;
    refiner.setMesh(sphMesh);
    for( int i = 0; i < numBinLevels-1; i++) {
        refiner.refineAll(14);
    }

    JMeshAffineTransform  affine;
    affine.setMesh(sphMesh);
    affine.toCenter();

    size_t numNodes = sphMesh->getSize(0);
    Point3D xyz;
    for( int i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = sphMesh->getNodeAt(i);
        Point3D xyz = vtx->getXYZCoords();
        double  x = xyz[0];
        double  y = xyz[1];
        double  z = xyz[2];
        double l = sqrt(x*x + y*y + z*z);
        xyz[0] = x/l;
        xyz[1] = y/l;
        xyz[2] = z/l;
        vtx->setXYZCoords(xyz);
    }

    JMeshIO mio;
    mio.saveAs( sphMesh, "sphtemplate.off");

    int  numFaces = sphMesh->getSize(2);
    JNodeSequence faceCenters;
    faceCenters.resize(numFaces);

    for( int  i = 0; i < numFaces; i++) {
        const JFacePtr &f = sphMesh->getFaceAt(i);
        f->setAttribute("Partition", i);
        JFaceGeometry::getCentroid(f, xyz);
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];
        double l = sqrt(x*x + y*y + z*z);
        xyz[0] = x/l;
        xyz[1] = y/l;
        xyz[2] = z/l;
        faceCenters[i] = JNode::newObject();
        faceCenters[i]->setID(i);
        faceCenters[i]->setXYZCoords(xyz);
    }

    jann.setCloud(faceCenters);
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshNormalClusters :: getQuantizedSphere()
{
    if( sphMesh == nullptr) buildSphere();
    return sphMesh;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalClusters :: createClusters()
{
    if( sphMesh == nullptr) buildSphere();

    Vec3F normal;
    Point3D query;
    size_t numFaces = mesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
         const JFacePtr &face = mesh->getFaceAt(i);
         face->getAttribute("Normal", normal);
         query[0] = normal[0];
         query[1] = normal[1];
         query[2] = normal[2];
         JNodePtr sphNode = jann.getNearest(query);
         int partID = sphNode->getID();
         face->setAttribute("Partition", partID);
     }
}

