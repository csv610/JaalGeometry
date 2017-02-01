#include "SimpleLaplaceSkeleton.hpp"

///////////////////////////////////////////////////////////////////////////////
void JSimpleLaplacianSkeleton :: setMesh( const JMeshPtr &m)
{
    jmesh = m->deepCopy();
    if( jmesh == nullptr)  return;

    jmesh->buildRelations(0,0);
    jmesh->buildRelations(0,2);
    jmesh->enumerate(0);
    jmesh->enumerate(2);

    size_t numNodes = jmesh->getSize(0);
    size_t numFaces = jmesh->getSize(2);

    mesh.nodes.resize(numNodes);

    JNodeSequence vneighs;
    JFaceSequence fneighs;
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = jmesh->getNodeAt(i);
        assert( vtx->getID() == i);
        JNode::getRelations(vtx, vneighs);
        for( const JNodePtr &node : vneighs) {
            size_t j = node->getID();
            mesh.nodes[i].adjNodes.push_back(j);
        }
        JNode::getRelations(vtx, fneighs);
        for( const JFacePtr &face : fneighs) {
            size_t j = face->getID();
            mesh.nodes[i].adjFaces.push_back(j);
        }
        mesh.nodes[i].xyz = vtx->getXYZCoords();
        mesh.nodes[i].lambda = 1.0;
    }

    mesh.faces.resize(numFaces);
    for( size_t i = 0; i < numFaces; i++) {
         const JFacePtr &face = jmesh->getFaceAt(i);
         for( int j = 0; j < 3; j++) {
             const JNodePtr &node = face->getNodeAt(j);
             mesh.faces[i].nodes[j] = node->getID();
         }
    }

}
///////////////////////////////////////////////////////////////////////////////

void JSimpleLaplacianSkeleton :: atomicOp(size_t vi)
{
    newCoords[vi][0] = 0.0;
    newCoords[vi][1] = 0.0;
    newCoords[vi][2] = 0.0;

    double xi = mesh.nodes[vi].xyz[0];
    double yi = mesh.nodes[vi].xyz[1];
    double zi = mesh.nodes[vi].xyz[2];

    int nSize = mesh.nodes[vi].adjNodes.size();
    for( int j = 0; j < nSize; j++) {
        int vj = mesh.nodes[vi].adjNodes[j];
        double dx = mesh.nodes[vj].xyz[0] - xi;
        double dy = mesh.nodes[vj].xyz[1] - yi;
        double dz = mesh.nodes[vj].xyz[2] - zi;
        double lambda = mesh.nodes[vi].lambda;
        newCoords[vi][0] +=  xi + lambda*dx;
        newCoords[vi][1] +=  yi + lambda*dy;
        newCoords[vi][2] +=  zi + lambda*dz;
    }
    newCoords[vi][0] /= (double)nSize;
    newCoords[vi][1] /= (double)nSize;
    newCoords[vi][2] /= (double)nSize;

    double dx = newCoords[vi][0] - xi;
    double dy = newCoords[vi][1] - yi;
    double dz = newCoords[vi][2] - zi;
 
    // Nodes must move inside..
    double nx = -mesh.nodes[vi].normal[0];
    double ny = -mesh.nodes[vi].normal[1];
    double nz = -mesh.nodes[vi].normal[2];
    double dl  = dx*nx + dy*ny + dz*nz;

    newCoords[vi][0] = xi + dl*nx;
    newCoords[vi][1] = yi + dl*ny;
    newCoords[vi][2] = zi + dl*nz;
//    mesh.nodes[vi].distance = dl;
}

///////////////////////////////////////////////////////////////////////////////

void JSimpleLaplacianSkeleton :: applyOneStep()
{
    static size_t ncount = 0;
    cout << "Laplacian Contraction " << ncount++ << endl;

    int numNodes = mesh.nodes.size();
    newCoords.resize(numNodes);

    setNodesNormal();

    for( size_t i = 0; i < numNodes; i++)
         mesh.nodes[i].lambda = 1.0;

    for( size_t i = 0; i < numNodes; i++) 
	atomicOp(i);

/*
    setNodesArea();
    double maxVal = 0.0;
    double minVal = std::numeric_limits<double>::max();

    for( size_t i = 0; i < numNodes; i++) {
         minVal = std::min( minVal, mesh.nodes[i].area);
         maxVal = std::max( minVal, mesh.nodes[i].area);
    }

    normArea.resize(numNodes);
    for( size_t i = 0; i < numNodes; i++)
         normArea[i] = mesh.nodes[i].area/maxVal;

    for( size_t i = 0; i < numNodes; i++)
         mesh.nodes[i].lambda = 1.0/(normArea[i] + 1.0E-06);

    setNodesAngle();
    for( size_t i = 0; i < numNodes; i++) {
         if( mesh.nodes[i].maxAngle > 170) 
             mesh.nodes[i].lambda = 0.001;
         else 
             mesh.nodes[i].lambda = 1.000;
    }

    for( size_t i = 0; i < numNodes; i++) atomicOp(i);
*/

    for( size_t i = 0; i < numNodes; i++)  {
         mesh.nodes[i].xyz =  newCoords[i];
         JNodePtr vtx   = jmesh->getNodeAt(i);
         vtx->setXYZCoords(newCoords[i] );
     }
}

///////////////////////////////////////////////////////////////////////////////
void JSimpleLaplacianSkeleton :: setNodesNormal()
{
    int numFaces = mesh.faces.size();
    for( size_t i = 0; i < numFaces; i++)
         setFaceNormal(i);

    int numNodes = mesh.nodes.size();
    for( size_t i = 0; i < numNodes; i++)
         setNodeNormal(i);
   
}
///////////////////////////////////////////////////////////////////////////////
void JSimpleLaplacianSkeleton :: setNodesArea()
{
    int numFaces = mesh.faces.size();
    for( size_t i = 0; i < numFaces; i++)
         setFaceArea(i);

    int numNodes = mesh.nodes.size();
    for( size_t i = 0; i < numNodes; i++)
         setNodeArea(i);
   
}
///////////////////////////////////////////////////////////////////////////////
