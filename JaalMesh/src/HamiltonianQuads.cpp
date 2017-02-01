#include "AllQuadMeshGenerator.hpp"

JMeshPtr AllQuadMeshGenerator :: getHamiltonianQuads()
{
    triMesh->pruneAll();
    JEdgeSequence orgEdges = triMesh->getEdges();

    JTriRefiner trefine;
    trefine.setMesh(triMesh);
    trefine.refineAll(13);
    quadMesh = triMesh;

    JFaceSequence efaces;
    JNodePtr v0, v1, ot1, ot2;

    size_t numedges = orgEdges.size();
    for( size_t iedge = 0; iedge < numedges; iedge++) {
        const JEdgePtr &edge = orgEdges[iedge];
        JEdge::getRelations( edge, efaces );
        if( efaces.size() == 2 ) {
            if( efaces[0]->isActive() && efaces[1]->isActive() ) {
                v0  = edge->getNodeAt(0);
                v1  = edge->getNodeAt(1);
                ot1 = JTriangle::getOppositeNode( efaces[0], v0, v1);
                ot2 = JTriangle::getOppositeNode( efaces[1], v0, v1);
                JFacePtr newquad =  JQuadrilateral::newObject(ot1, v0, ot2, v1);
                quadMesh->addObject(newquad);
                edge->setStatus( JMeshEntity::REMOVE );
                efaces[0]->setStatus( JMeshEntity::REMOVE );
                efaces[1]->setStatus( JMeshEntity::REMOVE );
            }
        }
    }

    numedges = quadMesh->getSize(1);
    for( size_t iedge = 0; iedge < numedges; iedge++) {
        const JEdgePtr &edge = quadMesh->getEdgeAt(iedge);
        if( edge->isActive() ) {
            JEdge::getRelations( edge, efaces );
            if( efaces.size() == 2 ) {
                if( efaces[0]->isActive() && efaces[1]->isActive() ) {
                    if( efaces[0]->getSize(0) == 3 && efaces[1]->getSize(0) == 3 ) {
                        v0  = edge->getNodeAt(0);
                        v1  = edge->getNodeAt(1);
                        ot1 = JTriangle::getOppositeNode( efaces[0], v0, v1);
                        ot2 = JTriangle::getOppositeNode( efaces[1], v0, v1);
                        JFacePtr newquad =  JQuadrilateral::newObject(ot1, v0, ot2, v1);
                        quadMesh->addObject(newquad);
                        edge->setStatus( JMeshEntity::REMOVE );
                        efaces[0]->setStatus( JMeshEntity::REMOVE );
                        efaces[1]->setStatus( JMeshEntity::REMOVE );
                    }
                }
            }
        }
    }
    quadMesh->pruneEdges();
    quadMesh->pruneFaces();

    JDoublet doublet;
    doublet.setMesh(quadMesh);
    doublet.removeAll();


    return quadMesh;
}

