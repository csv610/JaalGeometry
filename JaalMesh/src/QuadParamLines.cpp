#include "QuadParamLines.hpp"

using namespace Jaal;

void
QuadParametricLines :: basicOp(JEdgePtr edge)
{
    //////////////////////////////////////////////////////////////////////////
    //                   **********************
    //                   *         *          *
    //                   *         * Next     *
    //                   *         *          *
    //           Avoid   **********************  Avoid
    //                   *         *          *
    //                   *         * Current  *
    //                   *         *          *
    //                   *         *          *
    //                   **********************
    //                            Source
    // A Source vertex and Current edge is chosen.
    // We want to avoid two edges and want to select "Next" edge.
    //////////////////////////////////////////////////////////////////////////
    if( edge->getVisitBit() ) return;

//     Vertex *v0, *v1, *v2, *v3, *v4;

    bool param;
    int err = edge->getAttribute("ParametricLine", param);
    assert( err == 0);

    edge->setVisitBit(1);

    JFaceSequence adjFaces;
    JEdge::getRelations(edge, adjFaces);

    int numNeighs = adjFaces.size();
    JEdgePtr nextEdge;
    for( int i = 0; i < numNeighs; i++) {
        int pos = adjFaces[i]->getPosOf( edge );
        assert( pos >= 0);
        nextEdge = adjFaces[i]->getEdgeAt(pos+1);
        if( !nextEdge->getVisitBit() ) {
            nextEdge->setAttribute("ParametricLine", !param);
            edgeQ.push_back(nextEdge);
        }

        nextEdge = adjFaces[i]->getEdgeAt(pos+2);
        if( !nextEdge->getVisitBit() ) {
            nextEdge->setAttribute("ParametricLine", param);
            edgeQ.push_back(nextEdge);
        }

        nextEdge = adjFaces[i]->getEdgeAt(pos+3);
        if( !nextEdge->getVisitBit() ) {
            nextEdge->setAttribute("ParametricLine", !param);
            edgeQ.push_back(nextEdge);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

int QuadParametricLines :: setLines(JMeshPtr m)
{
    mesh = m;

    int nTopo = mesh->getTopology()->getElementsType(2);
    if (nTopo != JFace::QUADRILATERAL) {
        mesh->getLogger()->setWarn("For parametric lines, the mesh must be all-quads ");
        return 1;
    }

//   mesh->getTopology()->collect_edges();

    mesh->buildRelations(1, 2);
    mesh->buildRelations(0, 2);

    size_t numedges = mesh->getSize(1);

    JEdgePtr seedEdge = NULL;
    for( size_t i = 0; i < numedges; i++) {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            JNodePtr v0 = edge->getNodeAt(0);
            if( v0->getNumRelations(2) == 4 ) {
                seedEdge = edge;
                break;
            }
            JNodePtr v1 = edge->getNodeAt(1);
            if( v1->getNumRelations(2) == 4 ) {
                seedEdge = edge;
                break;
            }
        }
    }
    edgeQ.clear();
    mesh->deleteEdgeAttribute("ParametricLine" );

    mesh->setVisitBits(1,0);
    bool param = 0;
    seedEdge->setAttribute("ParametricLine", param);

    if( seedEdge )  {
        edgeQ.push_back( seedEdge);
        while( !edgeQ.empty() ) {
            JEdgePtr curredge = edgeQ.front();
            edgeQ.pop_front();
            basicOp( curredge );
        }
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

