#include "TriDecimator.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////
//
bool JTriDecimator :: isCollapsable(const JEdgePtr &edge, bool check567 )
{
    if( mesh == nullptr || edge == nullptr ) return 1;
    if( !edge->isActive()  ) return 0;

    JEdge::getRelations(edge, faceneighs);
    if( faceneighs.size() != 2 ) return 0;

    JNodePtr v0 = edge->getNodeAt(0);
    JNodePtr v1 = edge->getNodeAt(1);
    if( v0->isBoundary() || v1->isBoundary() ) return 0;

    JNode::getRelations(v0, f0neighs);
    remove(f0neighs.begin(), f0neighs.end(), faceneighs[0] );
    remove(f0neighs.begin(), f0neighs.end(), faceneighs[1] );

    int f0size = f0neighs.size()-2;
    if( f0size < 0) return 0;
    for( int i = 0; i < f0size; i++)
        if( f0neighs[i]->hasNode(v1) ) return 0;

    JNode::getRelations(v1, f1neighs);
    remove(f1neighs.begin(), f1neighs.end(), faceneighs[0]);
    remove(f1neighs.begin(), f1neighs.end(), faceneighs[1]);

    int f1size = f1neighs.size()-2;
    if( f1size < 0) return 0;
    for( int i = 0; i < f1size; i++)
        if( f1neighs[i]->hasNode(v0) ) return 0;

    if( check567 ) {
        JNodePtr vo1 =  JTriangle::getOppositeNode(faceneighs[0], v0, v1);
        if( vo1->getNumRelations(2) <= 6) return 0;

        JNodePtr vo2 =  JTriangle::getOppositeNode(faceneighs[1], v0, v1);
        if( vo2->getNumRelations(2) <= 6) return 0;

        if( f0size + f1size > 7) return 0;
    }

    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
int JTriDecimator :: collapse(const JEdgePtr &edge)
{
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    if(!isCollapsable( edge) ) return 1;

    JNodePtr v0 = edge->getNodeAt(0);
    JNodePtr v1 = edge->getNodeAt(1);

    JNodePtr vmid = JNodeGeometry::getMidNode(v0, v1);
    mesh->addObject(vmid);
    newNodes.push_back(vmid);

    int f0size = f0neighs.size()-2;

    for( int i = 0; i < f0size; i++) {
        JFacePtr face = f0neighs[i]->getClone();
        face->replace(v0, vmid);
        newFaces.push_back(face);
    }

    int f1size = f1neighs.size()-2;
    for( int i = 0; i < f1size; i++) {
        JFacePtr face = f1neighs[i]->getClone();
        face->replace(v1, vmid);
        newFaces.push_back(face);
    }

    mesh->addObjects(newFaces);

    JEdgePtr fedge;
    int pos;
    for( size_t i = 0; i < f0neighs.size(); i++) {
        pos = f0neighs[i]->getPosOf(v0);
        fedge = f0neighs[i]->getEdgeAt(pos);
        fedge->setStatus( JMeshEntity::REMOVE);
        fedge = f0neighs[i]->getEdgeAt(pos+2);
        fedge->setStatus( JMeshEntity::REMOVE);
        f0neighs[i]->setStatus(JMeshEntity::REMOVE);
    }

    for( size_t i = 0; i < f1neighs.size(); i++) {
        pos = f1neighs[i]->getPosOf(v1);
        fedge = f1neighs[i]->getEdgeAt(pos);
        fedge->setStatus( JMeshEntity::REMOVE);
        fedge = f1neighs[i]->getEdgeAt(pos+2);
        fedge->setStatus( JMeshEntity::REMOVE);
        f1neighs[i]->setStatus(JMeshEntity::REMOVE);
    }

    for( size_t i = 0; i < newFaces.size(); i++) {
        pos = newFaces[i]->getPosOf(vmid);
        fedge = newFaces[i]->getEdgeAt(pos);
        newEdges.push_back(fedge);
        fedge = newFaces[i]->getEdgeAt(pos+2);
        newEdges.push_back(fedge);
    }

    JEdgeSequence::iterator it;
    it  = std::unique( newEdges.begin(), newEdges.end() );
    newEdges.erase( it, newEdges.end() );

    mesh->addObjects(newEdges);

    faceneighs[0]->setStatus( JMeshEntity::REMOVE);
    faceneighs[1]->setStatus( JMeshEntity::REMOVE);
    edge->setStatus( JMeshEntity::REMOVE);
    v0->setStatus( JMeshEntity::REMOVE);
    v1->setStatus( JMeshEntity::REMOVE);

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void JTriDecimator :: updateRelations( const JFacePtr &face)
{
    int nd = face->getSize(0);
    for( int i = 0; i < nd; i++) {
        const JEdgePtr &edge =  face->getEdgeAt(i);
        const JNodePtr &v0   =  edge->getNodeAt(0);
        const JNodePtr &v1   =  edge->getNodeAt(1);
        v0->addRelation(v1);
        v1->addRelation(v0);
        v0->addRelation(face);
        v1->addRelation(face);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////
int JTriDecimator :: remove_degree3_node( const JNodePtr &vertex)
{
    JNodeSequence nodeneighs;
    JFaceSequence faceneighs;

    // Node(0,2), vertex(0,0) relations must be present...
    newFaces.clear();

    if( !vertex->isActive()  ) return 1;
    if( vertex->isBoundary() ) return 2;
    if( vertex->getNumRelations(2) != 3 ) return 3;

    nodeneighs.clear();
    JNode::getRelations(vertex, nodeneighs);
    if( nodeneighs.size() != 3 ) return 4;

    // Since every triangle vertex degree is going to be reduced by one, if you delete this
    // vertex, then you will create another vertex with degree three and this thing may not
    // terminate, therefore, use the refinement method for such nodes...
    for( int i = 0; i < 3; i++)
        if( nodeneighs[i]->getNumRelations(2) < 5) return 5;

    JNode::getRelations(vertex, faceneighs);
    for( int j = 0; j < 3; j++) {
        JFacePtr face = faceneighs[j];
        face->setStatus( JMeshEntity::REMOVE);
        JEdgePtr edge = JSimplex::getEdgeOf(vertex, nodeneighs[j] );
        if(edge) edge->setStatus( JMeshEntity::REMOVE);
    }
    vertex->setStatus( JMeshEntity::REMOVE);
    newFaces.resize(1);
    newFaces[0] = JTriangle::newObject( nodeneighs );
    if( mesh) mesh->addObject( newFaces[0] );
    updateRelations( newFaces[0] );

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

int JTriDecimator :: remove_degree4_node( const JNodePtr &vertex)
{
    // Node(0,2), vertex(0,0) relations must be present...
    newFaces.clear();
    newEdges.clear();

    if( !vertex->isActive()  ) return 1;
    if( vertex->isBoundary() ) return 2;
    if( vertex->getNumRelations(2) != 4 ) return 3;

    nodeneighs.clear();
    JNode::getRelations(vertex, nodeneighs);
    if( nodeneighs.size() != 4 ) return 4;

    JNode::getRelations(vertex, faceneighs);

    assert(faceneighs.size() == 4);
    coveredges.resize(4);
    for( int i = 0; i < 4; i++) {
        coveredges[i] =  JTriangle::getOppositeEdge(faceneighs[i], vertex);
    }

    // Assuming that the all the input triangles are consistently oriented, we shall attempt to make
    // output triangles consistent wi the input triangles also.
    int pos = faceneighs[0]->getPosOf(vertex);
    JNodePtr start_vertex = faceneighs[0]->getNodeAt(pos+1);

    JEdgeTopology::getChain( coveredges, start_vertex);

    JNodePtr v0 = coveredges[0]->getNodeAt(0);
    JNodePtr v1 = coveredges[1]->getNodeAt(0);
    JNodePtr v2 = coveredges[2]->getNodeAt(0);
    JNodePtr v3 = coveredges[3]->getNodeAt(0);
    JNodePtr retain[2];
    newFaces.clear();
    if( v1->getNumRelations(2) >= 6 && v3->getNumRelations(2) >= 6)  {
        newFaces.resize(2);
        newFaces[0] = JTriangle::newObject( v0, v1, v2);
        newFaces[1] = JTriangle::newObject( v0, v2, v3);
        retain[0]   = v0;
        retain[1]   = v2;
    }

    if( newFaces.empty() ) {
        if( v0->getNumRelations(2) >= 6 && v2->getNumRelations(2) >= 6)  {
            newFaces.resize(2);
            newFaces[0] = JTriangle::newObject( v0, v1, v3);
            newFaces[1] = JTriangle::newObject( v1, v2, v3);
            retain[0]   = v1;
            retain[1]   = v3;
        }
    }

    if( newFaces.empty() ) return 1;

    JEdgePtr oldedges[4];
    oldedges[0] = JSimplex::getEdgeOf(vertex, v0);
    oldedges[1] = JSimplex::getEdgeOf(vertex, v1);
    oldedges[2] = JSimplex::getEdgeOf(vertex, v2);
    oldedges[3] = JSimplex::getEdgeOf(vertex, v3);
    for( int j = 0; j < 4; j++) {
        JFacePtr face = faceneighs[j];
        face->setStatus( JMeshEntity::REMOVE);
        oldedges[j]->setStatus( JMeshEntity::REMOVE);
    }
    vertex->setStatus( JMeshEntity::REMOVE);
    newEdges.resize(1);
    newEdges[0] = JSimplex::getEdgeOf(retain[0], retain[1], 1);

    assert( newEdges[0] );

    if( mesh ) {
        mesh->addObjects( newFaces );
        mesh->addObjects( newEdges );
    }
    updateRelations( newFaces[0] );
    updateRelations( newFaces[1] );
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void JTriDecimator :: remove_degree3_nodes()
{
    JFaceSequence newfaces;  // Collect all new faces from each operation..
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            if( !vertex->isBoundary() ) {
                int nsize = vertex->getNumRelations(2);
                if( nsize == 3 )  {
                    int err = remove_degree3_node( vertex);
                    if( !err)
                        boost::copy( newFaces, std::back_inserter(newfaces));
                }
            }
        }
    }
    newFaces = newfaces;
}

//////////////////////////////////////////////////////////////////////////////////////////////

void JTriDecimator :: remove_degree4_nodes()
{
    JFaceSequence newfaces;  // Collect all new faces from each operation..
    JEdgeSequence newedges;  // Collect all new faces from each operation..

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            if( !vertex->isBoundary() ) {
                int nsize = vertex->getNumRelations(2);
                if( nsize == 4 )  {
                    int err = remove_degree4_node( vertex);
                    if( !err)  {
                        boost::copy( newEdges, std::back_inserter(newedges));
                        boost::copy( newFaces, std::back_inserter(newfaces));
                    }
                }
            }
        }
    }
    newFaces = newfaces;
    newEdges = newedges;
}

////////////////////////////////////////////////////////////////////////////////
void JTriDecimator :: removeLowDegreeNodes()
{
    newEdges.clear();
    newFaces.clear();

    if( mesh == nullptr) return;
    if( mesh->getAdjTable(0,2) == 0) mesh->buildRelations(0,2);
    if( mesh->getAdjTable(0,0) == 0) mesh->buildRelations(0,0);

    remove_degree3_nodes();
    remove_degree4_nodes();
}
//////////////////////////////////////////////////////////////////////////////////////////////////

void JTriDecimator :: smooth_local_patch(const JNodePtr &vapex, const JNodeSequence &neighs)
{
    double xsum = 0;
    double ysum = 0;
    double zsum = 0;
    int nSize = neighs.size();
    for( int i = 0; i < nSize; i++) {
        const Point3D &xyz = neighs[i]->getXYZCoords();
        xsum += xyz[0];
        ysum += xyz[1];
        zsum += xyz[2];
    }
    xsum /= (double)nSize;
    ysum /= (double)nSize;
    zsum /= (double)nSize;

    vapex->setXYZCoords(xsum, ysum, zsum);
}
////////////////////////////////////////////////////////////////////////////////

int JTriDecimator :: remove_above_degree8_node( const JNodePtr &vertex)
{
    // Node(0,2), vertex(0,0) relations must be present...
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    if( !vertex->isActive()  ) return 1;

//   We do not want to mess with the boundary nodes and don't want to change the connectivity too.
    if( vertex->isBoundary() ) return 2;

//   Node degree eight is handled separately..
    if( vertex->getNumRelations(2) < 9) return 3;

    JNode::getRelations(vertex, faceneighs);

    int numedges = faceneighs.size();
    coveredges.resize(numedges);
    for( int i = 0; i < numedges; i++) {
        coveredges[i] =  JTriangle::getOppositeEdge(faceneighs[i], vertex);
    }

//  Assuming that the input mesh is consistently oriented, we want that new faces are
//  also consistently oriented with the input mesh, so that we do not have to check
//  and validate the mesh every time. To do this, follow these steps.
//  Step-I:  the first face is consistently oriented, make the first edge correctly
//           oriented (i.e. counter-clockwise, if the input triangle is counter-clockwise.
//
    int edgesign = faceneighs[0]->getOrientation( coveredges[0] );
    if( edgesign == -1) coveredges[0]->reverse();

    JEdgeTopology::getChain( coveredges);
    JEdgeTopology::getChainNodes( coveredges, nodeneighs);

    JNodePtr start_vertex;
    for( int i = 0; i < numedges; i++) {
        int d1 = nodeneighs[i]->getNumRelations(2);
        int d2 = nodeneighs[(i+5)%numedges]->getNumRelations(2);
        if( d1 < 8 && d2 < 8 ) {
            start_vertex = nodeneighs[i];
            break;
        }
    }
    if( start_vertex == nullptr ) return 1;
    JEdgeTopology::getChain( coveredges, start_vertex);

    // Generate a new node ...
    JNodePtr vnew = JNode::newObject();
    vnew->setXYZCoords( vertex->getXYZCoords());
    vnew->setID( mesh->getSize(0) + 1);

    newNodes.resize(1);
    newNodes[0] = vnew;
    mesh->addObject(vnew);

    // The new node will have vertex degree of 7. First five segments
    // are glued to the new vertex.
    JNodeSequence newneighs(7);
    newneighs[0] = coveredges[0]->getNodeAt(0);
    newneighs[1] = coveredges[1]->getNodeAt(0);
    newneighs[2] = coveredges[2]->getNodeAt(0);
    newneighs[3] = coveredges[3]->getNodeAt(0);
    newneighs[4] = coveredges[4]->getNodeAt(0);
    newneighs[5] = coveredges[5]->getNodeAt(0);
    newneighs[6] = vnew;

    // Except the first five segments, all will remain with the old
    // vertex.
    JNodeSequence oldneighs( numedges-4);
    oldneighs[0] = coveredges[0]->getNodeAt(0);
    int index = 1;
    for( int i = 5; i < numedges; i++)
        oldneighs[index++] = coveredges[i]->getNodeAt(0);

    // Create new edges for the new vertex ...
    newEdges.resize(7);
    for( int i = 0; i < 6; i++)
        newEdges[i]  = JEdge::newObject( vnew, newneighs[i]);
    newEdges[6]  = JEdge::newObject( vnew, vertex);
    mesh->addObjects(newEdges);

    // New faces are inserted at the new vertex. The new vertex will have
    // seven neighbours.
    newFaces.resize(7);
    newFaces[0]  = JTriangle::newObject( vnew, vertex, newneighs[0]);
    for( int i = 0; i < 5; i++)
        newFaces[i+1]  = JTriangle::newObject( vnew,  newneighs[i], newneighs[i+1]);
    newFaces[6]  = JTriangle::newObject( vertex, vnew, newneighs[5]);
    mesh->addObjects( newFaces);
    assert( vnew->getNumRelations(2) == 7 );

    // Five faces at the old faces will be removed ....
    for( int i = 0; i < 5; i++) {
        JFacePtr face = JSimplex::getFaceOf( vertex, newneighs[i], newneighs[i+1] );
        assert(face);
        face->setStatus( JMeshEntity::REMOVE);
    }

    // Four edges are the old node will be removed ...
    for( int i = 0; i < 4; i++) {
        JEdgePtr edge = JSimplex::getEdgeOf( vertex, newneighs[i+1]);
        assert(edge);
        edge->setStatus( JMeshEntity::REMOVE);
    }

    // Keeping the local boundary fixed, we can relocate the two internal nodes
    // to improve the elements quality..
    for( int i = 0; i < 3; i++) {
        smooth_local_patch( vnew,   newneighs);
        smooth_local_patch( vertex, oldneighs);
    }

    if(!mesh->getTopology()->isConsistent() ) {
        cout << "Warning: Mesh is not consistent " << endl;
    }

    for( int i = 0; i < 7; i++) updateRelations( newFaces[i] );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JTriDecimator :: remove_degree8_node( const JNodePtr &vertex)
{
// Relations(0,2), Relations(0,0) relations must be present...
    newNodes.clear();
    newEdges.clear();
    newFaces.clear();

    if( !vertex->isActive()  ) return 1;
//  We do not to mess with the boundary nodes ...
    if( vertex->isBoundary() ) return 2;
    if( vertex->getNumRelations(2) != 8) return 3;

//  Collect all the edges covering the "vertex" at center...
    JNode::getRelations(vertex, faceneighs);
    int numedges = faceneighs.size();

    coveredges.resize(numedges);
    for( int i = 0; i < numedges; i++) {
        coveredges[i] =  JTriangle::getOppositeEdge(faceneighs[i], vertex);
    }

//  Assuming that the input mesh is consistently oriented, we want that new faces are
//  also consistently oriented with the input mesh, so that we do not have to check
//  and validate the mesh every time. To do this, follow these steps.
//  Step-I:  the first face is consistently oriented, make the first edge correctly
//           oriented (i.e. counter-clockwise, if the input triangle is counter-clockwise.
//
    int edgesign = faceneighs[0]->getOrientation( coveredges[0] );
    if( edgesign == -1) coveredges[0]->reverse();

//  Step 2: Keeping the first edge orientation fixed, get the entire edge consistently
//  oriented...

    JEdgeTopology::getChain( coveredges);
    JEdgeTopology::getChainNodes( coveredges, nodeneighs);

//  Now we should select two important nodes i.e. nodes[0] and nodes[5]. The degree of
//  these two nodes increases by one, therefore, we search for those pairs which have
//  vertex degrees less than 8, After operation, these two vertices will have degree = 7.

    JNodePtr start_vertex;

    for( int i = 0; i < 8; i++) {
        int d0  = nodeneighs[i]->getNumRelations(2);
        int d5  = nodeneighs[(i+5)%8]->getNumRelations(2);
        if( d0 < 8 && d5 < 8) {
            start_vertex = nodeneighs[i];
            break;
        }
    }
    if( start_vertex == nullptr) return 1;

//  Orient the coveredges so that the first vertex is the "start_vertex"...
    JEdgeTopology::getChain( coveredges, start_vertex);

//  One new vertex is created which will have 7 neigbours and the original
//  center vertex will have degree equal to 5.
//  We do not know the good position of the new vertex, so initialize it with
//  the position of the old center vertex, we will locally optimize the position
//  of the both nodes in the end...
    JNodePtr vnew = JNode::newObject();
    vnew->setXYZCoords( vertex->getXYZCoords());
    vnew->setID( mesh->getSize(0) + 1);

    newNodes.resize(1);
    newNodes[0] = vnew;
    mesh->addObject(vnew);

    JNodeSequence newneighs(7);
    newneighs[0] = coveredges[0]->getNodeAt(0);
    newneighs[1] = coveredges[1]->getNodeAt(0);
    newneighs[2] = coveredges[2]->getNodeAt(0);
    newneighs[3] = coveredges[3]->getNodeAt(0);
    newneighs[4] = coveredges[4]->getNodeAt(0);
    newneighs[5] = coveredges[5]->getNodeAt(0);
    newneighs[6] = vnew;

    JNodeSequence oldneighs(5);
    oldneighs[0] = coveredges[0]->getNodeAt(0);
    oldneighs[1] = coveredges[5]->getNodeAt(0);
    oldneighs[2] = coveredges[6]->getNodeAt(0);
    oldneighs[3] = coveredges[7]->getNodeAt(0);
    oldneighs[4] = vnew;

//  Creating new edges ...
    newEdges.resize(7);
    for( int i = 0; i < 6; i++)
        newEdges[i]  = JEdge::newObject( vnew, newneighs[i]);
    newEdges[6]  = JEdge::newObject( vnew, vertex);
    mesh->addObjects(newEdges);

//  Creating new faces ....
    newFaces.resize(7);
    newFaces[0]  = JTriangle::newObject( vnew, vertex, newneighs[0]);
    for( int i = 0; i < 5; i++)
        newFaces[i+1]  = JTriangle::newObject( vnew,  newneighs[i], newneighs[i+1]);
    newFaces[6]  = JTriangle::newObject( vertex, vnew, newneighs[5]);
    mesh->addObjects( newFaces);

//  Removing some of old faces ...
    for( int i = 0; i < 5; i++) {
        JFacePtr face = JSimplex::getFaceOf( vertex, newneighs[i], newneighs[i+1] );
        face->setStatus( JMeshEntity::REMOVE);
    }

//  Reoving some of old edges ....
    for( int i = 0; i < 4; i++) {
        JEdgePtr edge = JSimplex::getEdgeOf( vertex, newneighs[i+1]);
        edge->setStatus( JMeshEntity::REMOVE);
    }

    for( int i = 0; i < 3; i++) {
        smooth_local_patch( vnew,   newneighs);
        smooth_local_patch( vertex, oldneighs);
    }

    assert( vnew->getNumRelations(2)   == 7);
    assert( vertex->getNumRelations(2) == 5);

    for( int i = 0; i < 7; i++) updateRelations(newFaces[i]);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////
int JTriDecimator :: removeHighDegreeNodes()
{
    if( mesh == nullptr) return 1;
    if( mesh->getAdjTable(0,2) == 0) mesh->buildRelations(0,2);
    if( mesh->getAdjTable(0,0) == 0) mesh->buildRelations(0,0);

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            int nd = vtx->getNumRelations(2);
            if( nd >  8) remove_above_degree8_node( vtx );
            if( nd == 8) remove_degree8_node( vtx );
        }
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////
