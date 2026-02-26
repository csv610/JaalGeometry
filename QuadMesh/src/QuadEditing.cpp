#include "QuadEditing.hpp"

JEdgePtr QuadEdit :: getAnyEdge(int nv0, int nv1)
{
    size_t minindex = 0;
    size_t numedges = mesh->getSize(1);
    size_t index = JMath::random_value(minindex,numedges-1);
    JEdgePtr edge;

    for( size_t i = 0; i < numedges; i++) {
        edge = mesh->getEdgeAt( (index+i)%numedges);
        if( edge->isActive()  && !edge->isBoundary() ) {
            JNodePtr v0 = edge->getNodeAt(0);
            JNodePtr v1 = edge->getNodeAt(1);
            int n0 = v0->getNumRelations(2);
            int n1 = v1->getNumRelations(2);
            if( n0 == nv0 && n1 == nv1 ) return edge;
            if( n1 == nv0 && n0 == nv1 ) return edge;
        }
    }
    return NULL;
}

/*
int QuadEdit :: getAnyPair( int nv0, int nv1, NodePair &vpair)
{
    return 0;
}
*/
