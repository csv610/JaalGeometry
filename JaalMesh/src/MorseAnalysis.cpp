#include "MorseAnalysis.hpp"


void JMorseAnalysis :: setHeightDirection( int dir)
{
    if( mesh == NULL ) return;

    size_t nSize =  mesh->getSize(0);
    if( nSize == 0) return;

    for( size_t i = 0; i < nSize; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            const Point3D  &xyz  = vtx->getXYZCoords();
            vtx->setAttribute("MorseScalar", xyz[dir] );
        }
    }
}

