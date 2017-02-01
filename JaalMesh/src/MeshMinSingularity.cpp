#include "MeshMinSingularity.hpp"

void JMeshMinSingularity :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    refEdgeLength = mesh->getGeometry()->getMeanEdgeLength();
}

////////////////////////////////////////////////////////////////////////////////
void JMeshMinSingularity :: remeshBoundNode( const JNodePtr &apex)
{
    if( apex->getNumRelations(2) < 6) return;
    JEdgeSequence edges;
    mesh->getTopology()->getRim(apex, edges);

    int numedges = edges.size();
    if( numedges < 6) return;

    if( numedges%2 ) {
        cout << "Warning: number of boundary edges is not even" << endl;
        return;
    }

    JEdgeTopology::getChain(edges, apex);
    if( edges[0]->getNodeAt(0) != apex ) {
        cout << "Warning: Wrong statring node of the boundary loop " << endl;
        return;
    }

    if( edges[numedges]->getNodeAt(1) != apex ) {
        cout << "Warning: Wrong ending node of the boundary loop " << endl;
        return;
    }
}
////////////////////////////////////////////////////////////////////////////////

void JMeshMinSingularity :: refine( const JEdgePtr &edge)
{
    if( edge->isBoundary() ) return;
    if( newmesh == nullptr)  return;

    double elen = JEdgeGeometry::getLength(edge);
    int    numSegments = max(1.0, elen/refEdgeLength);

    if( numSegments%2 == 0) numSegments++;

    JNodeSequence nodes;
    JEdgeGeometry::generateLinearNodes( edge, numSegments+1, nodes);

    if( nodes.empty() ) return;

    newmesh->addObjects(nodes);
}

////////////////////////////////////////////////////////////////////////////////

void JMeshMinSingularity :: refine( const JFacePtr &face)
{
    int numedges = face->getSize(0);
    if( numedges > 6) return;

    if( newmesh == nullptr) return;

    int nCount  = 0;
    vector<int>   segments(numedges);
    JNodeSequence  nodes, srcnodes;

    int numedgenodes;
    for( int i = 0; i < numedges; i++) {
        const JEdgePtr &edge = face->getEdgeAt(i);
        const JNodePtr &v0   = edge->getNodeAt(0);
        const JNodePtr &v1   = edge->getNodeAt(1);
        srcnodes.push_back(face->getNodeAt(i));
        if( edge->isBoundary() )  {
            segments[i]  = 1;
            nCount++;
        } else {
            int err = edge->getAttribute("Steiner", nodes);
            if( !err ) {
                numedgenodes = nodes.size();
                segments[i] = numedgenodes+1;
                nCount +=  segments[i];
                int ori = face->getOrientation(edge);
                if( ori == 1) {
                    for( size_t j = 0; j < numedgenodes; j++)
                        srcnodes.push_back(nodes[j] );
                } else {
                    for( size_t j = 0; j < numedgenodes; j++)
                        srcnodes.push_back(nodes[numedgenodes-1-j] );
                }
            } else {
                segments[i] = 1;
                nCount++;
            }
        }
    }

    if( nCount%2 ) {
        cout << "Fatal error: Number of segments on a face not even" << endl;
        newmesh.reset();
        return;
    }

    JMeshPtr patchmesh =  polymesher.getPatch(srcnodes, segments);
    if( patchmesh == nullptr) return;

    JFaceSequence newfaces = patchmesh->getFaces();

    JNodeSequence newnodes;
    patchmesh->getTopology()->getSubmeshInternal(newfaces, newnodes );

    newmesh->addObjects( newnodes );
    newmesh->addObjects( newfaces );

    int pid = face->getID();
    for( const JFacePtr &face : newfaces)
        face->setAttribute("Partition", pid);

    for( const JFacePtr &face : newfaces) {
        double area = JFaceGeometry::getSignedArea(face);
        if( area < 0.0)  face->reverse();
    }

    face->setStatus(JMeshEntity::REMOVE);
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshMinSingularity :: refineAll()
{
    if( mesh == nullptr) return nullptr;

    newmesh = JMesh::newObject();

    size_t numnodes = mesh->getSize(0);
    for(size_t i = 0; i < numnodes; i++)
        newmesh->addObject( mesh->getNodeAt(i) );

    size_t numedges = mesh->getSize(1);
    for(size_t i = 0; i < numedges; i++)
        refine(mesh->getEdgeAt(i));

    size_t numfaces  = mesh->getSize(2);
    int    numColors = numfaces-1;
    for(size_t i = 0; i < numfaces; i++)
        refine(mesh->getFaceAt(i));

    for(size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->hasAttribute("Steiner") )
            edge->setStatus(JMeshEntity::REMOVE);
    }

    newmesh->getTopology()->searchBoundary();

    JSinglet singlet;
    singlet.setMesh(newmesh);
//  singlet.removeAll();

    // We want to separate faces shared by high degree vertex..
    JFaceSequence faceneighs;
    newmesh->buildRelations(0,2);
    numnodes = newmesh->getSize(0);
    JNodeSequence boundnodes;
    newmesh->getTopology()->getBoundary(boundnodes);

    JNodeSequence highdegree;
    for( const JNodePtr &vtx : boundnodes) {
        if( vtx->getNumRelations(2) > 5)
            highdegree.push_back(vtx);
    }

    newmesh->getTopology()->collectEdges();

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(newmesh);
    mopt.setBoundaryPreserve(1);
    mopt.setNumIterations(1000);
    mopt.improveQuality(1);

    return newmesh;
}
////////////////////////////////////////////////////////////////////////////////
