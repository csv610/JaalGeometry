#include "RingQuadMesh.hpp"

void JRingQuadMesh :: setMesh(const JMeshPtr &m)
{
   mesh = m;
}

//////////////////////////////////////////////////////////////////////
void JRingQuadMesh :: setSkeleton(const vector<JEdgeSequence> &e)
{
   skeleton = e;
}

void JRingQuadMesh :: createRings( const JEdgeSequence &branch)
{
}

//////////////////////////////////////////////////////////////////////
void JRingQuadMesh :: createRings()
{
    if( skeleton.empty() ) {

    }
}
//////////////////////////////////////////////////////////////////////
