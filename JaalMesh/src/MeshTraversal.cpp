#include "MeshTraversal.hpp"

/*
void JMeshTraversal :: getTopologicalSequence(const JFacePtr &seedface, JFaceSequence &sequence, size_t maxSize, int method)
{
    if( method == BREADTH_FIRST_ORDER ) {
        mesh->buildRelations(1, 2);

        bfs.clear();

        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            JFacePtr f = mesh->getFaceAt(i);
            f->setVisitBit(0);
        }

        list<JFacePtr> faceQ;
        faceQ.push_back(seedface);

        JFaceSequence neighs;

        while(1) {
            while(!faceQ.empty() ) {
                JFacePtr currface = faceQ.front();
                faceQ.pop_front();
                if( currface->isActive() ) {
                    currface->getRelations12(neighs);
                    for( size_t i = 0; i < neighs.size(); i++)
                        if( !neighs[i]->getVisitBit() ) faceQ.push_back(neighs[i] );
                    if( !currface->getVisitBit() ) {
                        bfs.push_back(currface);
                        if( bfs.size() == maxSize ) return;
                        currface->setVisitBit(1);
                    }
                }
            }

            for( size_t i = 0; i < numfaces; i++) {
                JFacePtr f = mesh->getFaceAt(i);
                if( !f->getVisitBit() && f->isActive() ) {
                    faceQ.push_back( f );
                    break;
                }
            }
            if( faceQ.empty() ) break;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void JMeshTraversal :: getTopologicalSequence(const JNodePtr &seednode, JFaceSequence &sequence, size_t maxSize, int method)
{

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
void JMeshTraversal :: init02()
{
    double infty =  0.999*std::numeric_limits<double>::max();

    mesh->buildRelations(0, 2);

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        v->setVisitBit(0);
        v->setAttribute("Distance", infty);
    }

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr f = mesh->getFaceAt(i);
        f->setVisitBit(0);
    }

    initialized = 1;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void JMeshTraversal :: getGeometricSequence(const JNodePtr &seednode, JFaceSequence &sequence, double maxDistance, int method)
{
    double infty =  0.999*std::numeric_limits<double>::max();

    sequence.clear();

    if( !initialized) init02();

    if( method == BREADTH_FIRST_ORDER ) {
        JFaceSequence faces, facesChanged;
        JNodeSequence nodesChanged;

        list<JNodePtr> nodeQ;
        nodeQ.push_back(seednode);
        double idist, jdist;
        idist = 0.0;
        seednode->setAttribute("Distance", idist);

        while(!nodeQ.empty() ) {
            JNodePtr inode = nodeQ.front();
            nodeQ.pop_front();
            if( inode->getVisitBit() == 1) continue;
            inode->getAttribute("Distance", idist);
            if( inode->isActive() ) {
                JNode::getRelations(inode, faces);
                for( JFacePtr face : faces)  {
                    int  inrange = 1;
                    int  numnodes = face->getSize(0);
                    for( int i = 0; i < numnodes; i++) {
                        JNodePtr jnode = face->getNodeAt(i);
                        if( inode != jnode) {
                            jnode->getAttribute("Distance", jdist);
                            double dij  = JNodeGeometry::getLength(inode, jnode);
                            double dmin =  min( idist + dij, jdist);
                            if(  dmin < maxDistance)  {
                                jnode->setAttribute("Distance", dmin);
                                nodeQ.push_back( jnode);
                            }  else
                                inrange = 0;
                        }
                    }
                    if( inrange) sequence.push_back( face );
                    face->setVisitBit(1);
                    facesChanged.push_back(face);
                }
                inode->setVisitBit(1);
                nodesChanged.push_back(inode);
            }

            // Bring back the changed so that you can call this function for other seeds. Initialization
            // is expensive ....
            for( JNodePtr vtx: nodesChanged )  {
                vtx->setAttribute("Distance", infty);
                vtx->setVisitBit(0);
            }
            for( JFacePtr face : facesChanged)
                face->setVisitBit(0);
        }
    }
}
