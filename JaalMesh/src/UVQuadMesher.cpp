#include "UVQuadMesher.hpp"

////////////////////////////////////////////////////////////////////////////////

void JUVQuadMesher :: remesh()
{
    if( mesh == nullptr) return;
    int err = mesh->getAttribute("UVMesh", uvMesh);

    if( err ) {
        cout << "Warning: No UV mesh associated with the mesh " << endl;
        return;
    }

    size_t numfaces = uvMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = uvMesh->getFaceAt(i);
        double area = JFaceGeometry::getSignedArea(f);
        if( area < 0.0) {
            cout << "Fatal Error: The UV mesh has an inverted element" << endl;
            return;
        }
    }
     
    newUVMesh = uvMesh->deepCopy();
    JAlphaMSTQuadMesh amst;
    amst.setMesh(newUVMesh);
    amst.remeshAll();

    JMeshPtr qmesh = uvMesh->deepCopy();
    AllTriMeshGenerator alltri;
    JMeshPtr uvTriMesh = alltri.getFromQuadMesh(qmesh);

    JPointLocation pLocator;
    pLocator.setMesh(uvTriMesh);

    Point3D xyz;
    Point3D bary;
    Point3D pQuery;

    newMesh = JMesh::newObject();

    // Update the interior coordinates;
    size_t numNodes = newUVMesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++)  {
        const JNodePtr &vi = newUVMesh->getNodeAt(i);
        if( !vi->isBoundary() ) {
            pQuery = vi->getXYZCoords();
            JFacePtr uvface = pLocator.searchFace(pQuery, 1);
            const Point3D &p0 = uvface->getNodeAt(0)->getXYZCoords();
            const Point3D &p1 = uvface->getNodeAt(1)->getXYZCoords();
            const Point3D &p2 = uvface->getNodeAt(2)->getXYZCoords();
            JTriGeometry::getBaryCoords( &p0[0], &p1[0], &p2[0], &pQuery[0],  &bary[3]);
            int fid = uvface->getID();
            JTriGeometry::getXYZCoordinates(uvMesh->getFaceAt(fid), bary, xyz);
            JNodePtr newnode = JNode::newObject();
            newnode->setXYZCoords(xyz);
            newMesh->addObject(newnode);
        } else 
            newMesh->addObject(mesh->getNodeAt(i));
      }
}

////////////////////////////////////////////////////////////////////////////////
