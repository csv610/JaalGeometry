#include "RangeSearch.hpp"

void JRangeSearch :: setQueryRegion( JFacePtr f)
{
    if( mesh == nullptr ) return;
    if( mesh->getAdjTable(0,2) == 0) mesh->buildRelations(0,2);

    bfs_search(f);
}

/////////////////////////////////////////////////////////////////////////////////

void JRangeSearch :: bfs_search(JFacePtr seedface)
{
    overlapFaces.clear();

    deque<JNodePtr> nodeQ;
    JNodeSet nodesVisited;

    JFaceSet overlapSet;
    overlapSet.insert(seedface);

    int numnodes = seedface->getSize(0);
    for( int i = 0; i < numnodes; i++)
        nodeQ.push_back( seedface->getNodeAt(i) );

    JFaceSequence faceneighs;
    while( !nodeQ.empty() ) {
        JNodePtr currnode = nodeQ.front();
        nodeQ.pop_front();
        if( nodesVisited.find(currnode) == nodesVisited.end() ) {
            JNode::getRelations(currnode, faceneighs);
            for( JFacePtr nextface: faceneighs) {
                if( overlapSet.find(nextface) == overlapSet.end() ) {
                    if( JFaceGeometry::intersectPredicate2d( seedface, nextface ) ) {
                        overlapSet.insert( nextface);
                        numnodes = nextface->getSize(0);
                        for( int i = 0; i < numnodes; i++) {
                            JNodePtr vtx = nextface->getNodeAt(i);
                            if( nodesVisited.find(vtx) == nodesVisited.end() )
                                nodeQ.push_back(vtx);
                        }
                    }
                }
            }
            nodesVisited.insert( currnode);
        }
    }

    if(overlapSet.size() > 1)
    {
        overlapFaces.resize( overlapSet.size() );
        copy( overlapSet.begin(), overlapSet.end(), overlapFaces.begin() );
    }

}
////////////////////////////////////////////////////////////////////////////////
void JRangeSearch :: bfs_search(JNodePtr moveVertex)
{
    overlapFaces.clear();
    JFaceSequence faceneighs;
    JNode::getRelations(moveVertex, faceneighs);

    JFaceSet overlapSet;
    for( size_t i = 0;  i < faceneighs.size(); i++) {
        bfs_search(faceneighs[i] );
        copy( overlapFaces.begin(), overlapFaces.end(), inserter( overlapSet, overlapSet.end() ) );
    }

    if( !overlapSet.empty() ) {
        overlapFaces.resize( overlapSet.size() );
        copy( overlapSet.begin(), overlapSet.end(), overlapFaces.begin() );
    }
}
////////////////////////////////////////////////////////////////////////////////
