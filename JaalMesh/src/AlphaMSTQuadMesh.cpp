#include "AlphaMSTQuadMesh.hpp"

void JAlphaMSTQuadMesh :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    nearSearch.reset( new JNearestNeighbours);
    JNodeSequence bnodes;
    mesh->getTopology()->getBoundary(bnodes);
    nearSearch->setCloud(bnodes);
    double elen  =  mesh->getGeometry()->getMeanEdgeLength();
    meanArea = elen*elen;
    mesh->buildRelations(0,2);
}

//////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMesh :: clear()
{
    patchNodes.clear();
    newNodes.clear();
    patchEdges.clear();
    newEdges.clear();
    patchFaces.clear();
    newFaces.clear();
}

//////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMesh :: setCircle(const JCircle &c)
{
    patchNodes.clear();
    patchFaces.clear();
    patchEdges.clear();
    circle = c;
}
//////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMesh :: setCenter(const Point3D &center)
{
    patchNodes.clear();
    patchFaces.clear();
    patchEdges.clear();

    if( nearSearch == nullptr) return;
    JNodePtr  vtx = nearSearch->getNearest(center);
    double    r   = JMath::length( vtx->getXYZCoords(), center);
    circle.setCenter(center);
    circle.setRadius(r);
}

//////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMesh :: buildPatch()
{
    if( mesh == nullptr) return;

    Point3D  center = circle.getCenter();
    double   radius = circle.getRadius();

    patchNodes.clear();

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        const Point3D  &p0  = vtx->getXYZCoords();
        double dist         = JMath::length(p0, center);
        if( dist < radius) patchNodes.push_back(vtx);
    }
    boost::sort( patchNodes);

    patchFaces.clear();

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        bool allinside = 1;
        for( int j = 0; j < nnodes; j++) {
            if( boost::binary_search(patchNodes, face->getNodeAt(j) ) == 0) {
                allinside = 0;
                break;
            }
        }
        if( allinside ) patchFaces.push_back(face);
    }

    mesh->getTopology()->getSubmeshBoundary( patchFaces, patchEdges);

    JNodeSequence nodes;
    JMeshTopology::getEntitySet( patchFaces, nodes);
    numSingularities = 0;
    for( const JNodePtr  &vtx : nodes) {
        if( !vtx->isBoundary() && vtx->getNumRelations(2) != 4)
            numSingularities++;
    }

    JEdgeTopology::getChain(patchEdges);
    JEdgeTopology::getChainNodes(patchEdges, patchNodes);
}

//////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMesh :: remeshPatch()
{
    if( patchEdges.empty() ) return;

    int numBoundEdges = patchEdges.size();
    if( numBoundEdges%2 ) {
        cout << "Warning: number of boundary edges not even " << endl;
        return;
    }
    vector<int> segments(4);
    segments[0] = numBoundEdges/4;
    segments[1] = numBoundEdges/4;
    segments[2] = numBoundEdges/4;
    segments[3] = numBoundEdges/4;

    int nleft   = numBoundEdges%4;
    assert( nleft <= 2);
    if( nleft ) {
        segments[0] += 1;
        nleft--;
    }
    if( nleft ) {
        segments[2] += 1;
        nleft--;
    }
    assert( nleft == 0);

    adaptFactor = 2.0*M_PI*circle.getRadius()/(numBoundEdges*expectedEdgeLength);

    JMSTQuadMesher polymesher;
    JMeshPtr meshtemplate = polymesher.getPatch(patchNodes, segments);
//  JMeshPtr meshtemplate = polymesher.getAdaptivePatch(patchNodes, segments, adaptFactor);

    // Entities to be removed from the mesh ..
    JEdgeSequence oldEdges;
    JNodeSequence oldNodes;
    mesh->getTopology()->getSubmeshInternal(patchFaces, oldEdges, oldNodes);

    for( const JFacePtr &face: patchFaces)
        face->setStatus(JMeshEntity::REMOVE);

    for( const JEdgePtr &edge: oldEdges)
        edge->setStatus(JMeshEntity::REMOVE);

    for( const JNodePtr &vtx: oldNodes)
        vtx->setStatus(JMeshEntity::REMOVE);

    newFaces = meshtemplate->getFaces();
    newNodes.clear();
    newEdges.clear();

    // Entiies to be added to the mesh ...
    meshtemplate->getTopology()->getSubmeshInternal(newFaces, newEdges, newNodes);

    mesh->addObjects(newNodes);
    mesh->addObjects(newEdges);
    mesh->addObjects(newFaces);

    // Check if any doublet is formed at the patch boundary. Remove it ...
    mesh->buildRelations(0,2);
    JFaceSequence faceneighs;

    JDoublet doublet;
    doublet.setMesh(mesh);
    for( const JNodePtr &vtx : patchNodes) {
        JNode::getRelations(vtx, faceneighs);
        if( faceneighs.empty() ) vtx->setStatus( JMeshEntity::REMOVE);
        if( faceneighs.size() == 2 ) doublet.remove(vtx);
    }
    mesh->pruneAll();

    // Everything is done. Clear the patch for the new patch.

    patchNodes.clear();
    patchEdges.clear();
    patchFaces.clear();

    // Smooth the patch with Laplacian. Domain must be convex ....
    JLaplaceMeshSmoother smooth;
    smooth.setMesh(mesh);
    smooth.setNumIterations(100);
    smooth.smooth(newNodes);

    for( const JFacePtr &face : newFaces ) {
        double area = JFaceGeometry::getSignedArea(face);
        if( area < 0.0) face->reverse();
    }

    // Report how many singularities are left in the patch..
    JNodeSequence nodes;
    JMeshTopology::getEntitySet( newFaces, nodes);
    numSingularities = 0;
    for( const JNodePtr  &vtx : nodes) {
        if( !vtx->isBoundary() && vtx->getNumRelations(2) != 4)
            numSingularities++;
    }
}
//////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMesh :: remeshAll()
{
/*
    if( mesh == nullptr) return;
  
    vector<JEdgeSequence> allboundedges;
    mesh->getTopology()->getBoundary(allboundedges);
    if( allbounedges.size() != 1) {
        cout << "Warning: There should be only one boundary loop " << endl;
        return;
    }
    patchEdges = allboundedges[0];

    if( patchEdges.empty() ) return;

    int numBoundEdges = patchEdges.size();
    if( numBoundEdges%2 ) {
        cout << "Warning: number of boundary edges not even " << endl;
        return;
    }
    vector<int> segments(4);
    segments[0] = numBoundEdges/4;
    segments[1] = numBoundEdges/4;
    segments[2] = numBoundEdges/4;
    segments[3] = numBoundEdges/4;

    int nleft   = numBoundEdges%4;
    assert( nleft <= 2);
    if( nleft ) {
        segments[0] += 1;
        nleft--;
    }
    if( nleft ) {
        segments[2] += 1;
        nleft--;
    }
    assert( nleft == 0);

    adaptFactor = 1.0;

    JMSTQuadMesher polymesher;
    JMeshPtr meshtemplate = polymesher.getPatch(patchNodes, segments);
//  JMeshPtr meshtemplate = polymesher.getAdaptivePatch(patchNodes, segments, adaptFactor);

    // Entities to be removed from the mesh ..
    JEdgeSequence oldEdges;
    JNodeSequence oldNodes;
    mesh->getTopology()->getSubmeshInternal(patchFaces, oldEdges, oldNodes);

    for( const JFacePtr &face: patchFaces)
        face->setStatus(JMeshEntity::REMOVE);

    for( const JEdgePtr &edge: oldEdges)
        edge->setStatus(JMeshEntity::REMOVE);

    for( const JNodePtr &vtx: oldNodes)
        vtx->setStatus(JMeshEntity::REMOVE);

    newFaces = meshtemplate->getFaces();
    newNodes.clear();
    newEdges.clear();

    // Entiies to be added to the mesh ...
    meshtemplate->getTopology()->getSubmeshInternal(newFaces, newEdges, newNodes);

    mesh->addObjects(newNodes);
    mesh->addObjects(newEdges);
    mesh->addObjects(newFaces);

    // Check if any doublet is formed at the patch boundary. Remove it ...
    mesh->buildRelations(0,2);
    JFaceSequence faceneighs;

    JDoublet doublet;
    doublet.setMesh(mesh);
    for( const JNodePtr &vtx : patchNodes) {
        JNode::getRelations(vtx, faceneighs);
        if( faceneighs.empty() ) vtx->setStatus( JMeshEntity::REMOVE);
        if( faceneighs.size() == 2 ) doublet.remove(vtx);
    }
    mesh->pruneAll();

    // Everything is done. Clear the patch for the new patch.

    patchNodes.clear();
    patchEdges.clear();
    patchFaces.clear();

    // Smooth the patch with Laplacian. Domain must be convex ....
    JLaplaceMeshSmoother smooth;
    smooth.setMesh(mesh);
    smooth.setNumIterations(100);
    smooth.smooth(newNodes);

    for( const JFacePtr &face : newFaces ) {
        double area = JFaceGeometry::getSignedArea(face);
        if( area < 0.0) face->reverse();
    }

    // Report how many singularities are left in the patch..
    JNodeSequence nodes;
    JMeshTopology::getEntitySet( newFaces, nodes);
    numSingularities = 0;
    for( const JNodePtr  &vtx : nodes) {
        if( !vtx->isBoundary() && vtx->getNumRelations(2) != 4)
            numSingularities++;
    }
*/
}
//////////////////////////////////////////////////////////////////////
