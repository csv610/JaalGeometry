#include "MeshExtrude.hpp"

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshExtrude :: getMesh( double distance, int ntimes)
{
    if( mesh == nullptr) return nullptr;

    JMeshPtr newMesh = JMesh::newObject();
    newMesh->addObject(mesh);

    JMeshPtr lastlayer = mesh;
    JEdgeSequence boundEdges[2];
    size_t numNodes = mesh->getSize(0);
    size_t numFaces = mesh->getSize(2);
    JQuadrilateralPtr newQuad;
    JNodeSequence qConnect(4);

    for( int i = 0; i < ntimes; i++) {
        lastlayer->getTopology()->getBoundary( boundEdges[0] );
        JMeshPtr newlayer = lastlayer->deepCopy();
        newlayer->getTopology()->getBoundary( boundEdges[1] );
        for( size_t j = 0; j < numNodes; j++) {
            const JNodePtr &node0 = lastlayer->getNodeAt(j);
            const JNodePtr &node1 = newlayer->getNodeAt(j);
            Point3D p0 = node0->getXYZCoords();
            p0[2] += distance;
            node1->setXYZCoords(p0);
        }

        newMesh->addObject( newlayer );
        size_t numEdges = boundEdges[0].size();
        for( size_t j = 0; j < numEdges; j++) {
            const JEdgePtr &edge0 = boundEdges[0][j];
            const JEdgePtr &edge1 = boundEdges[1][j];
            qConnect[0] = edge0->getNodeAt(0);
            qConnect[1] = edge0->getNodeAt(1);
            qConnect[2] = edge1->getNodeAt(1);
            qConnect[3] = edge1->getNodeAt(0);
            newQuad =  JQuadrilateral::newObject(qConnect);
            newMesh->addObject(newQuad);
        }
        for( size_t j = 0; j < numFaces; j++) {
            const JFacePtr &face = newlayer->getFaceAt(j);
            face->reverse();
         }
         lastlayer = newlayer;
    }
    return newMesh;
}

////////////////////////////////////////////////////////////////////////////////
