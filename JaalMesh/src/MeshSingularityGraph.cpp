#include "MeshSingularityGraph.hpp"

using namespace Jaal;

JLogger* JMeshSingularityGraph :: logger = JLogger::getInstance();

/////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraph :: clear()
{
    if( mesh == nullptr) return;
    mesh->deleteNodeAttribute("Interface");
    mesh->deleteEdgeAttribute("Interface");
    mesh->deleteFaceAttribute("Partition");
    mesh->deleteCellAttribute("Partition");
}

/////////////////////////////////////////////////////////////////////////////

int
JMeshSingularityGraph ::advance_single_quad_step()
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
    JNodePtr v0, v1, v2, v3, v4;

    JNodeSet vset;

    size_t index = nodeSeq.size();
    v0 = nodeSeq[index - 2];
    v1 = nodeSeq[index - 1];
    JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
    assert(edge);
    edge->setAttribute("Interface", 1);

    JNodeSequence vneighs;
    JNode::getRelations( v1, vneighs );
    if (vneighs.size() != 4) return 1; // Terminate at singular node

    JFaceSequence adjFaces;
    JEdge::getRelations(edge, adjFaces);
    assert(edge);
    if( adjFaces.size() != 2 ) return 1;
    v2 = JFace::getDiagonalNode(adjFaces[0], v0);
    v3 = JFace::getDiagonalNode(adjFaces[1], v0);

    vset.clear();
    vset.insert(vneighs[0]);
    vset.insert(vneighs[1]);
    vset.insert(vneighs[2]);
    vset.insert(vneighs[3]);
    vset.erase(v0);
    vset.erase(v2);
    vset.erase(v3);
    assert(vset.size() == 1);
    v4 = *vset.begin();
    edge = JSimplex::getEdgeOf(v1,v4);
    assert(edge);
    edge->setAttribute("Interface", 1);
    if( v4 == nodeSeq.front() ) return 1; // Terminate at self loop...

    if( stopAt == STOP_AT_CROSSINGS)  {
        if( v4->hasAttribute("Interface") ) return 1;
    }

    nodeSeq.push_back(v4);
    v4->setAttribute("Interface", 1);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshSingularityGraph ::create_new_branch( const JNodePtr &node1, const JNodePtr &node2)
{
    int pid = 1;
    nodeSeq.resize(2);
    nodeSeq[0] = node1;
    nodeSeq[1] = node2;

    node1->setAttribute("Interface", pid);
    node2->setAttribute("Interface", pid);

    while (1) {
        int terminate = advance_single_quad_step();
        if ( terminate ) break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void  JMeshSingularityGraph :: setAllPaths()
{
    if( mesh == nullptr) return;

    clear();

    diskTopology = !mesh->getTopology()->isClosed();

    int ntopo = mesh->getTopology()->getElementsType(2);
    if( ntopo == JFace::TRIANGLE) setTriPaths();
    if( ntopo == JFace::QUADRILATERAL) setQuadPaths();
}

///////////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraph :: setPaths( const JNodeSequence &nodes)
{
    if( mesh == nullptr) return;

    int ntopo = mesh->getTopology()->getElementsType(2);
    if( ntopo == 3) setTriPaths( nodes);
    if( ntopo == 4) setQuadPaths(nodes);
}

///////////////////////////////////////////////////////////////////////////////////

int JMeshSingularityGraph :: makeCyclicNodes(const JNodePtr &vertex)
{
    if( !vertex->isActive() ) return 1;
    if( diskTopology && vertex->isBoundary() ) return 2;
    if( vertex->getNumRelations(2) != 6 ) return 3;

    JNode::getRelations(vertex, faceSeq);

    edgeSeq.resize(6);
    for( int i = 0; i < 6; i++)  {
        edgeSeq[i] = JTriangle::getOppositeEdge(faceSeq[i], vertex);
    }

    JEdgeTopology::getChain(edgeSeq);

    nodeSeq.resize(6);
    for( int i = 0; i < 6; i++) {
        const JNodePtr &v0 = edgeSeq[i]->getNodeAt(1);
        const JNodePtr &v1 = edgeSeq[(i+1)%6]->getNodeAt(0);
        if( v0 != v1 ) return 2;
        nodeSeq[i] = v0;
    }

    // Now the vertices are clockwise or anticlockwise order.
    vertex->setRelations(nodeSeq);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraph :: advance_tri_path( const JNodePtr &v0, const JNodePtr &v1)
{
    bool val = 1;
    JEdgePtr edge = JSimplex::getEdgeOf(v0, v1);
    edge->setAttribute("Interface", val);

    if( diskTopology && v1->isBoundary() ) return;

    if( v1->getNumRelations(2) != 6 ) return;
    if( v1->hasAttribute("Interface") ) return;

    v0->setAttribute("Interface", val);
    v1->setAttribute("Interface", val);

    JNode::getRelations(v1, nodeSeq);
    assert( nodeSeq.size() == 6 );

    int pos = -1;
    for( int i = 0; i < 6; i++) {
        if( nodeSeq[i] == v0) {
            pos = i;
            break;
        }
    }
    assert(pos >= 0);

    JNodePtr nextVertex = nodeSeq[(pos+3)%6];
    advance_tri_path(v1, nextVertex);
}

///////////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraph :: setTriPath( const JNodePtr &vtx)
{
    JNodeSequence vneighs;
    JNode::getRelations(vtx, vneighs);
        for( size_t j = 0; j< vneighs.size(); j++)
            advance_tri_path( vtx, vneighs[j] );
}
///////////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraph :: setTriPaths( const JNodeSequence &singularNodes)
{
    for( const JNodePtr &vtx: singularNodes)  setTriPath(vtx);
}
///////////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraph :: setTriPaths()
{
    size_t numnodes = mesh->getSize(0);
    JNodeSequence irregularNodes;

    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vertex = mesh->getNodeAt(i);
        if( vertex->getNumRelations(2) != 6 ) 
            irregularNodes.push_back(vertex);
    }
    setTriPaths(irregularNodes);
}

///////////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraph :: setQuadPath( const JNodePtr &vtx)
{
    JNodeSequence vnodes;
    JNode::getRelations( vtx, vnodes );
    if( vnodes.empty() ) {
        mesh->buildRelations(0,0);
        JNode::getRelations( vtx, vnodes );
    }

    if( vnodes.empty() ) {
        cout << "Warning: Vertex had no neighbours " << endl;
        return;
    }

    for( const JNodePtr &vneigh: vnodes) 
        create_new_branch(vtx, vneigh );
}
///////////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraph :: setQuadPaths( const JNodeSequence &singularNodes)
{
    for( const JNodePtr &vtx: singularNodes) setQuadPath(vtx);
}
///////////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraph :: setQuadPaths()
{
    if( mesh == nullptr ) return;
    singularNodes.clear();

    size_t numnodes = mesh->getSize(0);

    JFaceSequence vfaces;
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        JNode::getRelations( vtx, vfaces );
        if( diskTopology &&  !vtx->isBoundary()  ) {
            if( vfaces.size() != 4 ) singularNodes.push_back(vtx);
        } else {
            if( vfaces.size() != 4 ) singularNodes.push_back(vtx);
        }

    }
    setQuadPaths(singularNodes);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraph :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    mesh->buildRelations(0,2);
    mesh->buildRelations(0,0);
    mesh->getTopology()->searchBoundary();
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JMeshSingularityGraph :: getEdges()
{
    JEdgeSequence eseq;

    if( mesh == nullptr ) return eseq;

    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
         const JEdgePtr &e = mesh->getEdgeAt(i);
         if( e->hasAttribute("Interface") )
             eseq.push_back(e);
    }
    return eseq;
}
///////////////////////////////////////////////////////////////////////////////

JNodeSequence JMeshSingularityGraph :: getNodes()
{
   return singularNodes;
}

///////////////////////////////////////////////////////////////////////////////

/*
JNodeSequence JMotorcycleGraph :: getJunctions()
{
    return singularNodes;
}
///////////////////////////////////////////////////////////////////////////////
int JMotorcycleGraph :: getPartitions( const JNodeSequence &nodes)
{
    if( mesh == nullptr ) {
        logger->setWarn("A null mesh object passed to Motorcycle graph");
        return 1;
    }

    logger->setInfo("Calculating Motorcycle graph paths");
    getPaths(nodes);

    JMeshPartitioner  mpart;
    mpart.setMesh(mesh);
    mpart.searchRegions();
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JMotorcycleGraph :: getPartitions()
{
    if( mesh == nullptr ) return 1;

    setAllPaths();

    JMeshPartitioner  mpart;
    mpart.setMesh(mesh);
    mpart.searchRegions();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMotorcycleGraph :: optimize()
{
    cout << "Not yet implemented " << endl;
    return 0;
}
*/
