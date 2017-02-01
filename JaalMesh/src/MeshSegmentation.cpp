#include "MeshSegmentation.hpp"

int MeshSegmentation ::getPartitions()
{
    if( mesh == NULL ) return 1;

    int relexist2 = mesh->buildRelations(1, 2);

    size_t nSize = mesh->getSize(2);
    for( size_t i = 0; i < nSize; i++) {
        JFacePtr f = mesh->getFaceAt(i);
        f->setAttribute( "Partition", 0);
        f->setID( i );
    }

    deque<JFacePtr> faceQ;
    JFaceSequence fneighs;
    int partid = 1;
    int pid = 0;

    while(1) {
        JFacePtr seedface = NULL;
        for( size_t i = 0; i < nSize; i++) {
            JFacePtr f = mesh->getFaceAt(i);
            if( f->isActive()  ) {
                f->getAttribute("Partition", pid);
                if( pid == 0) {
                    seedface = f;
                    break;
                }
            }
        }
        if( seedface == NULL ) break;

        faceQ.clear();
        faceQ.push_back(seedface);

        while(!faceQ.empty() ) {
            JFacePtr currface = faceQ.front();
            faceQ.pop_front();
            if( currface->isActive() ) {
                currface->getAttribute( "Partition", pid) ;
                if( pid == 0) {
                    currface->setAttribute( "Partition", partid) ;
                    int nnodes = currface->getSize(0);
                    for( int i = 0; i < nnodes; i++) {
                        JEdgePtr edge = currface->getEdgeAt(i);
                        if( !hasCollided(edge) ) {
                            JEdge::getRelations(edge, fneighs);
                            for( size_t j = 0; j < fneighs.size(); j++) {
                                fneighs[j]->getAttribute("Partition", pid);
                                if( pid == 0) faceQ.push_back( fneighs[j] );
                            }
                        }
                    }
                }
            }
        }
        partid++;
    }

    if (!relexist2) mesh->clearRelations(1, 2);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////
