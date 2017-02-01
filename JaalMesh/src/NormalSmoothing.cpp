#include "NormalSmoothing.hpp"

////////////////////////////////////////////////////////////////////////////////////////
void JMeshNormalSmoothing :: atomicOp( const JNodePtr &apex)
{
    JNodeSequence vneighs;
    JNode::getRelations( apex, vneighs);

    double xsum = 0.0;
    double ysum = 0.0;
    double zsum = 0.0;
    int nSize = vneighs.size();

    if( nSize == 0) cout << "Warning: Node is isolated " << endl;

    Vec3F normal;
    for( const JNodePtr &vtx: vneighs) {
        vtx->getAttribute("Normal", normal);
        xsum += normal[0];
        ysum += normal[1];
        zsum += normal[2];
    }
    xsum /= ( double) nSize;
    ysum /= ( double) nSize;
    zsum /= ( double) nSize;
    apex->setAttribute("Normal", normal);
}
////////////////////////////////////////////////////////////////////////////////////////

void JMeshNormalSmoothing :: atomicOp( const JFacePtr &face)
{
    double xsum = 0.0;
    double ysum = 0.0;
    double zsum = 0.0;
    int nSize = face->getSize(0);

    Vec3F normal;
    for( int i = 0; i < nSize; i++) {
        const JNodePtr &vtx = face->getNodeAt(i);
        vtx->getAttribute("Normal", normal);
        xsum += normal[0];
        ysum += normal[1];
        zsum += normal[2];
    }
    xsum /= ( double) nSize;
    ysum /= ( double) nSize;
    zsum /= ( double) nSize;
    face->setAttribute("Normal", normal);
}

////////////////////////////////////////////////////////////////////////////////////////

void JMeshNormalSmoothing :: execute()
{
    if( mesh == nullptr) return;

    mesh->buildRelations(0,0);

    size_t numnodes = mesh->getSize(0);
    for( int i = 0; i < numIterations; i++) {
        for( size_t j = 0; j < numnodes; j++) atomicOp( mesh->getNodeAt(j) );
    }

    size_t numfaces = mesh->getSize(2);
    for( size_t j = 0; j < numfaces; j++) atomicOp( mesh->getFaceAt(j) );
}

////////////////////////////////////////////////////////////////////////////////////////
