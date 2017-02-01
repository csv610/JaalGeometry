#include "MeshSkeleton.hpp"

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeleton :: clear()
{
    faceCentroid.clear();
    nodeSDF.clear();
    leafNodes.clear();
    junctionNodes.clear();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeleton :: storeTriMesh()
{
    if( inMesh == nullptr) return;

    int dim = inMesh->getTopology()->getDimension();
    JMeshPtr  surfMesh;

    // If the mesh have 3 cells, extract the boundary mesh ...
    if( dim == 3)
        surfMesh = inMesh->getTopology()->getSurfaceMesh();
    else
        surfMesh = inMesh;

    int nc =  surfMesh->getTopology()->getNumComponents();
    if( nc  > 1)  {
        cout << "Erro: Mesh skeleton requires single component mesh" << endl;
        return;
    }

    AllTriMeshGenerator alltri;
    int elemType = surfMesh->getTopology()->getElementsType(2);
    if( elemType == JFace::QUADRILATERAL) {
        JMeshPtr tmpmesh = surfMesh->deepCopy();
        triMesh =  alltri.getFromQuadMesh(tmpmesh,4);
    } else if( elemType == JFace::TRIANGLE)
        triMesh = surfMesh->deepCopy();

    string name = inMesh->getName() + std::to_string(numRefine);
    triMesh->setName(name);

    if( numRefine) {
        // Refine each triangle into 4 subtriangles ...
        JTriRefiner trefine;
        trefine.setMesh(triMesh);
        for( int i = 0; i < numRefine; i++)
            trefine.refineAll(14);
    }

    // Calculate the centroid of each triangle ....
    size_t numTris = triMesh->getSize(2);
    faceCentroid.resize( numTris);
    for( size_t i = 0; i < numTris; i++) {
        const JFacePtr &tri = triMesh->getFaceAt(i);
        JFaceGeometry::getCentroid( tri, faceCentroid[i]);
    }

    // Enumerate the nodes and store the file into the "off" format.
    triMesh->enumerate(0);
    JMeshIO meshio;
    meshio.saveAs(triMesh, "tmp.off");
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeleton :: initMesh()
{
    if( inMesh == nullptr) return;

    storeTriMesh();

#ifdef USE_CGAL
    // Read the mesh from the file ...
    std::ifstream input( "tmp.off");
    input >> tmesh;
    mcs.reset( new Skeletonization(tmesh));
#endif
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshSkeleton :: getWorkingMesh()
{

#ifdef USE_CGAL
    if( mcs == nullptr) initMesh();
    return triMesh;
#endif
    return nullptr;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeleton :: setNodeRadius(const JNodePtr &vtx)
{
    vector<size_t> surfnodes;
    vtx->getAttribute("SurfaceNodes", surfnodes);

    double minDist = std::numeric_limits<double>::max();
    double dist;
    for(size_t id : surfnodes) {
        const JNodePtr &vs = triMesh->getNodeAt(id);
        dist =  JNodeGeometry::getLength( vtx, vs);
        minDist =  min(dist, minDist);
    }
    vtx->setAttribute("Radius", minDist);

}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeleton :: genSkeleton()
{
    if( skelGraph ) return;
    
#ifdef USE_CGAL
    if( mcs == nullptr) initMesh();

    mcs->set_quality_speed_tradeoff(speedQuality);
    mcs->set_medially_centered_speed_tradeoff(centerQuality);

    mcs->contract_until_convergence();

    Skeleton skeleton;
    mcs->convert_to_skeleton(skeleton);

    // Output all the edges of the skeleton.
    std::ofstream output("skel.cgal");
    SkelPolylines display(skeleton,output);
    CGAL::split_graph_into_polylines(skeleton, display);
    output.close();

    skelGraph = JMesh::newObject();
    Point3D xyz;
    BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton)) {
        xyz[0] = skeleton[v].point.x();
        xyz[1] = skeleton[v].point.y();
        xyz[2] = skeleton[v].point.z();
        JNodePtr v0 = JNode::newObject();
        v0->setXYZCoords(xyz);
        v0->setID(v);
        skelGraph->addObject(v0);
    }

    BOOST_FOREACH(Skeleton_edge e, edges(skeleton))
    {
        auto vsrc = source(e,skeleton);
        const Point& s = skeleton[vsrc].point;
        xyz[0] = s.x();
        xyz[1] = s.y();
        xyz[2] = s.z();
        JNodePtr v0  = skelGraph->getNodeAt(vsrc);

        auto vdst = target(e,skeleton);
        const Point& t = skeleton[vdst].point;
        xyz[0] = t.x();
        xyz[1] = t.y();
        xyz[2] = t.z();
        JNodePtr v1  = skelGraph->getNodeAt(vdst);

        JEdgePtr edge = JEdge::newObject(v0,v1);
        skelGraph->addObject(edge);
    }

    vector<size_t> surfnodes;
    BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton)) {
        JNodePtr svtx = skelGraph->getNodeAt(v);
        surfnodes.clear();
        BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
        surfnodes.push_back((size_t)vd);
        svtx->setAttribute("SurfaceNodes", surfnodes);
        setNodeRadius(svtx);
    }
    processGraph();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void JMeshSkeleton :: contractGeometry()
{
   static size_t ncount = 0;
#ifdef USE_CGAL
   if( mcs == nullptr) initMesh();

   mcs->set_quality_speed_tradeoff(speedQuality);
   mcs->set_medially_centered_speed_tradeoff(centerQuality);

   mcs->contract();
   auto P = mcs->meso_skeleton();
   std::ofstream out("meso.off");
   CGAL::set_ascii_mode(out);
   out << P << endl;

   JMeshIO mio;
   triMesh = mio.readFile("meso.off");
#endif
}
////////////////////////////////////////////////////////////////////////////////

void JMeshSkeleton :: processGraph()
{
    if( skelGraph == nullptr) return;

    map<JNodePtr, JNodeSequence> vRelations;

    size_t numEdges = skelGraph->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = skelGraph->getEdgeAt(i);
        const JNodePtr &v0 =  edge->getNodeAt(0);
        const JNodePtr &v1 =  edge->getNodeAt(1);
        vRelations[v0].push_back(v1);
        vRelations[v1].push_back(v0);
    }

    size_t nSize = skelGraph->getSize(0);
    leafNodes.clear();
    junctionNodes.clear();
    for( size_t i = 0; i < nSize; i++) {
        JNodePtr vtx  = skelGraph->getNodeAt(i);
  
        if( vRelations[vtx].size() == 1) leafNodes.push_back(vtx);
        if( vRelations[vtx].size() >  2) junctionNodes.push_back(vtx);
    }
    // Keep them sorted, used while extracting branches.
    boost::sort( leafNodes);
    boost::sort( junctionNodes);

    segmentBranches();

    segmentSurface();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeleton :: segmentBranches()
{
    if( skelGraph == nullptr) return;

    // Mark all points on the graph to 0.
    size_t numNodes = skelGraph->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &node = skelGraph->getNodeAt(i);
        node->setVisitBit(0);
    }
    // Mark junction point to 1 which will not allow the next process pass through
    // the junction nodes ...
    for( const JNodePtr &vtx : junctionNodes)
        vtx->setVisitBit(1);

    skelGraph->buildRelations(0,1);

    // Now perform BFS on the graph to identify the strongly connected branches ...
    JNodePtr seedNode, currNode, nextNode;
    deque<JNodePtr> nodeQ;
    JEdgeSequence adjEdges;
    int partID = 0;

    while(1) {
        // Identify the seed node ...
        seedNode.reset();
        for( size_t i = 0; i < numNodes; i++) {
            const JNodePtr &node = skelGraph->getNodeAt(i);
            if( !node->getVisitBit() ) {
                seedNode = node;
                break;
            }
        }
        if( seedNode == nullptr) break;

        nodeQ.clear();
        nodeQ.push_back(seedNode);
        while(!nodeQ.empty() )  {
            currNode = nodeQ.front();
            nodeQ.pop_front();
            if( !currNode->getVisitBit() ) {
                currNode->setVisitBit(1);
                currNode->setAttribute("Partition", partID);
                JNode::getRelations(currNode, adjEdges);
                for( const JEdgePtr &e : adjEdges) {
                    e->setAttribute("Partition", partID);
                    nextNode = e->getNodeAt(0);
                    if( !nextNode->getVisitBit()) nodeQ.push_back(nextNode);
                    nextNode = e->getNodeAt(1);
                    if( !nextNode->getVisitBit()) nodeQ.push_back(nextNode);
                }
            }
        }
        partID++;
    }

    int numBranches = partID;

    branches.clear();
    for( int i = 0; i < numBranches; i++)  {
        auto newbranch = getNewBranch(i);
        branches.push_back(newbranch);
    }

    boost::sort(branches, []( const JMeshSkeletonBranchPtr &a, const JMeshSkeletonBranchPtr &b)
    {
        return a->getLength() < b->getLength();
    });
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSkeleton :: segmentSurface()
{
    int branchID = -1;
    size_t numnodes = inMesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = inMesh->getNodeAt(i);
        vtx->setAttribute("Partition", branchID);
    }

    int numBranches = branches.size();
    JNodeSequence skelNodes;

    vector<size_t>  surfnodes;
    for( int i = 0; i < numBranches; i++) {
        auto branch = branches[i];
        int bid = branch->id;
        JMeshTopology::getEntitySet( branch->skelEdges, skelNodes);
        for( const JNodePtr &sv : skelNodes) {
            sv->getAttribute("SurfaceNodes", surfnodes);
            for( auto srfnode : surfnodes) {
                const JNodePtr &snode = inMesh->getNodeAt(srfnode);
                snode->setAttribute("Partition", bid);
            }
        }
    }

    vector<int> nodePart;
    size_t numfaces = inMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = inMesh->getFaceAt(i);
        branchID = -1;
        face->setAttribute("Partition", branchID);
        int nn = face->getSize(0);
        nodePart.resize(nn);
        for( int j = 0; j < nn; j++) {
            const JNodePtr &vtx = face->getNodeAt(j);
            vtx->getAttribute("Partition", nodePart[j] );
        }
        int minID = *std::min_element(nodePart.begin(), nodePart.end() );
        if( minID >= 0) {
            branchID = nodePart[0];
            bool all_equal = 1;
            for( int id : nodePart) if( id != branchID ) all_equal = 0 ;
            if( all_equal ) face->setAttribute("Partition", branchID);
        }
    }

    //
    // Faces which have not been sigend the ParittionID so far are ambigious and
    // they do not belong to any speciiic branch, assign them the highest ID i.e.
    // "numBranches"...
    //
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &node = inMesh->getNodeAt(i);
        node->getAttribute("Partition", branchID);
        if( branchID < 0) node->setAttribute("Partition", numBranches);
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = inMesh->getFaceAt(i);
        face->getAttribute("Partition", branchID);
        if( branchID < 0) face->setAttribute("Partition", numBranches);
    }
}

///////////////////////////////////////////////////////////////////////////////

JNodePtr JMeshSkeleton :: getStartNode( JEdgeSequence &edges)
{
    JNodeSequence edgeNodes;
    JMeshTopology::getEntitySet(edges, edgeNodes);
    boost::sort( edgeNodes);

    // First check if there is a leaf node on the branch.
    JNodeSequence commNodes;
    boost::set_intersection( edgeNodes, leafNodes, back_inserter(commNodes) );

    if( !commNodes.empty() ) return commNodes[0];

    commNodes.clear();
    // First check if there is a junction node on the branch.
    boost::set_intersection( edgeNodes, junctionNodes, back_inserter(commNodes) );
    if( !commNodes.empty() ) return commNodes[0];

    // In rare case, the edge could be perfectly cyclic with no leaf or junction
    // node. In such case, the starting node is arbitrary...

    return edgeNodes[0];
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeleton :: classify( JMeshSkeletonBranchPtr &branch)
{
   char  str[2];
   str[0] = '0';
   str[1] = '0';

   JNodePtr node;

/*
   node = skelEdges->front()->getNodeAt(0);
   if( boost::find(leafNodes, node)) 
       str[0] = '1';
   else if( boost::find(junctionNodes, node)) 
       str[0] = '3';

    node = skelEdges->back()->getNodeAt(1);
    if( boost::find(leafNodes, node)) 
        str[1] = '1';
    else if( boost::find(junctionNodes, node2)) 
         str[1] = '3';
*/

    branch->type = atoi(str);
}
///////////////////////////////////////////////////////////////////////////////

JMeshSkeletonBranchPtr JMeshSkeleton :: getNewBranch(int partID)
{
    if( skelGraph == nullptr) return nullptr;

    JMeshSkeletonBranchPtr branch(new JMeshSkeletonBranch);
    branch->id = partID;

    size_t numEdges = skelGraph->getSize(1);
    int id = 0;

    JEdgeSequence skelEdges;
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = skelGraph->getEdgeAt(i);
        edge->getAttribute("Partition", id);
        if( id == partID ) skelEdges.push_back(edge);
    }

    if( skelEdges.empty() ) return nullptr;

    JNodePtr startNode = getStartNode( skelEdges );
    JEdgeTopology::getChain( skelEdges, startNode);

    branch->skelEdges = skelEdges;

    classify( branch );

    return branch;
}

///////////////////////////////////////////////////////////////////////////////
bool JMeshSkeletonBranch :: isTopologicalCylinder()
{

}
///////////////////////////////////////////////////////////////////////////////

double JMeshSkeletonBranch :: getMeanRadius()
{
    /*
        JFaceSequence branchfaces = getBranchSurface(id);

        JNodeSequence surfnodes;
        JMeshTopology::getEntitySet( branchfaces, surfnodes);
        int numnodes = surfnodes.size();
        double sum = 0.0;
        for (int i = 0; i < numnodes; i++)
             sum += nodeSDF[surfnodes[i]->getID() ];

        return sum/(double)numnodes;
    */
}
///////////////////////////////////////////////////////////////////////////////

/*
JFaceSequence JMeshSkeleton :: getBranchSurface(int partID)
{
    return branch;
}
*/

void JMeshSkeletonBranch :: fitCylinder()
{
    /*
        double radius = getMeanRadius(id);

        JFaceSequence branchfaces = getBranchSurface(id);
        JNodeSequence surfnodes;
        JMeshTopology::getEntitySet( branchfaces, surfnodes);

        Point3D p0, p1;
        p0[0] = 0.0;
        p0[1] = 0.0;
        p0[2] = 0.0;
        for( const JNodePtr &vtx : surfnodes) {
              p1  = vtx->getXYZCoords();
              p0[1] = p1[1];
              double dx = p1[0] - p0[0];
              double dy = p1[1] - p0[1];
              double dz = p1[2] - p0[2];
              double dl = sqrt(dx*dx + dy*dy + dz*dz);
              p1[0] = p0[0] + radius*dx/dl;
              p1[1] = p0[1] + radius*dy/dl;
              p1[2] = p0[2] + radius*dz/dl;
             vtx->setXYZCoords(p1);
         }
    */
}
///////////////////////////////////////////////////////////////////////////////
double JMeshSkeletonBranch :: getLength()
{
    double sum = 0.0;
    double elen;
    for( const JEdgePtr &edge : skelEdges)
        sum += JEdgeGeometry::getLength(edge);
    return sum;
}
///////////////////////////////////////////////////////////////////////////////
