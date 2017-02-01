#include "NearestNeighbours.hpp"

JNearestNeighbours :: JNearestNeighbours()
{
}
//////////////////////////////////////////////////////////////////////////

JNearestNeighbours :: ~JNearestNeighbours()
{
    if( dataPts )    annDeallocPts(dataPts);
    if( queryPoint ) annDeallocPt(queryPoint);
}


//////////////////////////////////////////////////////////////////////////
void JNearestNeighbours :: setMesh( const JMeshPtr &mesh)
{
    inCloud.clear();
    kdTree.reset();
    if( mesh) inCloud = mesh->getNodes();
}
//////////////////////////////////////////////////////////////////////////
void JNearestNeighbours :: setCloud( const JNodeSequence &seq)
{
    kdTree.reset();
    inCloud = seq;
}

void JNearestNeighbours :: buildTree()
{
    if( inCloud.empty()) return;

    int index = 0;
    Point3D xyz;

    size_t numNodes = inCloud.size();
    dataPts = annAllocPts( numNodes, 3);
    queryPoint = annAllocPt(3);

    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = inCloud[i];
        vtx->setID(index);
        xyz = vtx->getXYZCoords();
        dataPts[index][0] = xyz[0];
        dataPts[index][1] = xyz[1];
        dataPts[index][2] = xyz[2];
        index++;
    }
    kdTree.reset(new ANNkd_tree(dataPts, index, 3));
}

//////////////////////////////////////////////////////////////////////////

JNodeSequence  JNearestNeighbours :: getNearest(const Point3D &qPoint, int k)
{
    JNodeSequence nodes;
    if( k  < 1) return nodes;
    if( kdTree == nullptr) buildTree();

    queryPoint[0] = qPoint[0];
    queryPoint[1] = qPoint[1];
    queryPoint[2] = qPoint[2];

    nnIdx = new ANNidx[k];
    dists = new ANNdist[k];

    double eps = 0.0;
    kdTree->annkSearch( queryPoint, k, nnIdx, dists, eps);

    nodes.resize(k);
    for( int i = 0; i  < k; i++)  {
        nodes[i] = inCloud[nnIdx[i]];
    }

    delete [] nnIdx;
    delete [] dists;

    return nodes;
}

//////////////////////////////////////////////////////////////////////////
JNodePtr  JNearestNeighbours :: getNearest(const Point3D &qPoint)
{
    JNodeSequence  nodes = getNearest(qPoint, 1);
    return nodes[0];
}

//////////////////////////////////////////////////////////////////////////
