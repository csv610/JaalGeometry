#include "BinaryTreeMatch.hpp"
#include "StopWatch.hpp"
#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int JBinaryTreeMatch::match_boundary_triangle(JFacePtr btri)
{
    if (btri->getSize(0) != 3) return 1;

    JNodePtr bv0 = nullptr;
    JNodePtr bv1 = nullptr;
    JNodePtr bv2 = nullptr;

    for (int i = 0; i < 3; i++) {
        JNodePtr ev1 = btri->getNodeAt(i + 1);
        JNodePtr ev2 = btri->getNodeAt(i + 2);
        if (ev1->isBoundary() && ev2->isBoundary()) {
            bv0 = ev1;
            bv1 = ev2;
            bv2 = btri->getNodeAt(i);
            break;
        }
    }

    if (bv0 == nullptr || bv1 == nullptr) return 2;

    btri->setStatus( JMeshEntity::REMOVE );

    Point3D p3d = JNodeGeometry::getMidPoint(bv0, bv1);

    JNodePtr bound = JNode::newObject();
    bound->setXYZCoords(p3d);
    trimesh->addObject(bound);

    JFacePtr tri0 = JTriangle::newObject( bv0, bound, bv2 );
    trimesh->addObject(tri0);

    JFacePtr tri1 = JTriangle::newObject( bound, bv1, bv2 );
    trimesh->addObject(tri1);

    steinerNodes.push_back(bound);

    JFacePair facepair;
    facepair.first  = tri0;
    facepair.second = tri1;
    facematching.push_back(facepair);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////

void JBinaryTreeMatch::splitParent(TreeNode *parent, TreeNode *child1,
                                   TreeNode *child2)
{
    JNodePtr dnode = nullptr;
    dnode = parent->getDualNode();

    JFacePtr parentface = nullptr;
    dnode->getAttribute("PrimalFace", parentface);

    dnode = child1->getDualNode();
    JFacePtr face1 = nullptr;
    dnode->getAttribute("PrimalFace", face1);

    dnode = child2->getDualNode();
    JFacePtr face2 = nullptr;
    dnode->getAttribute("PrimalFace", face2);

    JNodeSequence connect(3);

    // Remove all existing vertex-face relations;
    JNodePtr vertex;
    for (int i = 0; i < 3; i++) {
        vertex = parentface->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face1->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face2->getNodeAt(i);
        vertex->clearRelations(2);
    }

    // Rebuild vertex-face relations...
    for (int i = 0; i < 3; i++) {
        vertex = parentface->getNodeAt(i);
        vertex->addRelation(parentface);

        vertex = face1->getNodeAt(i);
        vertex->addRelation(face1);

        vertex = face2->getNodeAt(i);
        vertex->addRelation(face2);
    }

    dnode = nullptr;
    parentface->getAttribute("DualNode", dnode);
    JNodePtr steiner = dnode->getClone();
    steiner->setID(parentface->getID());
    trimesh->addObject(steiner);
    steinerNodes.push_back(steiner);

    int edge1, edge2, edge3;

    edge1 = edge2 = edge3 = -1;
    JFaceSequence neighs;
    for (int i = 0; i < 3; i++) {
        JNodePtr v0 = parentface->getNodeAt(i + 1);
        JNodePtr v1 = parentface->getNodeAt(i + 2);
        JEdgePtr edge = JSimplex::getEdgeOf(v0,v1);
        JEdge::getRelations(edge, neighs);

        if (neighs.size() == 1)
            edge3 = i;

        if (neighs.size() == 2) {
            if (find(neighs.begin(), neighs.end(), face1) != neighs.end())
                edge1 = i;
            if (find(neighs.begin(), neighs.end(), face2) != neighs.end())
                edge2 = i;
        }
    }

    JFacePtr qface;
    JNodePtr ev0, ev1;

    // Match Child1 and One of the Split Triangle ...
    maxfaceID++;
    ev0 = parentface->getNodeAt(edge1 + 1);
    ev1 = parentface->getNodeAt(edge1 + 2);
    connect[0] = steiner;
    connect[1] = ev0;
    connect[2] = ev1;
    qface = JTriangle::newObject( connect );
    qface->setID(maxfaceID);
    JNodePtr dc1 = JMeshDualGraph::newObject(qface );
    dc1->setID(maxfaceID);
    trimesh->addObject(qface);
    steinerFaces.push_back(qface);

    dnode = nullptr;
    face1->getAttribute("DualNode", dnode);
    matchnodes( dnode, dc1);
    TreeNode *bnode1 = new TreeNode(dc1);
    bnode1->setMatchMark(1);
    bnode1->setParent(parent);
    bnode1->addChild(child1);
    btree->addNode(bnode1);

    // Match Child2 and One of the Split Triangle ...
    maxfaceID++;
    ev0 = parentface->getNodeAt(edge2 + 1);
    ev1 = parentface->getNodeAt(edge2 + 2);
    connect[0] = steiner;
    connect[1] = ev0;
    connect[2] = ev1;
    qface = JTriangle::newObject( connect );
    qface->setID(maxfaceID);
    JNodePtr dc2 = JMeshDualGraph::newObject(qface);
    dc2->setID(maxfaceID);
    trimesh->addObject(qface);
    steinerFaces.push_back(qface);
    dnode = nullptr;
    face2->getAttribute( "DualNode", dnode);
    matchnodes( dnode, dc2);

    TreeNode *bnode2 = new TreeNode(dc2);
    bnode2->setMatchMark(1);
    bnode2->setParent(parent);
    bnode2->addChild(child2);
    btree->addNode(bnode2);

    // Now Parent have different connectivity ...
    ev0 = parentface->getNodeAt(edge3 + 1);
    ev1 = parentface->getNodeAt(edge3 + 2);
    connect[0] = steiner;
    connect[1] = ev0;
    connect[2] = ev1;
    parentface->setNodes(connect);
    Point3D p3d;
    parentface->getAvgXYZ( p3d );

    JNodePtr dc3 = nullptr;
    parentface->getAttribute("DualNode", dc3);
    dc3->setXYZCoords(p3d);
    parent->addChild(bnode1);
    parent->addChild(bnode2);
    modifiedFaces.push_back(parentface);

    for (int i = 0; i < 3; i++) {
        vertex = parentface->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face1->getNodeAt(i);
        vertex->clearRelations(2);

        vertex = face2->getNodeAt(i);
        vertex->clearRelations(2);
    }

    child1->setMatchMark(1);
    child2->setMatchMark(1);

    btree->removeNode(child1);
    btree->removeNode(child2);
}

////////////////////////////////////////////////////////////////////////////////

void JBinaryTreeMatch::matchnode(TreeNode* v)
{
    TreeNode *parv = v->getParent();

    if (parv == nullptr)
        return;

    int degree = parv->getDegree();

    if (parv->isRoot() && degree == 2) {
        TreeNode *vsib = v->getSibling();
        splitParent(parv, v, vsib);
        return;
    }

    if (degree == 1 || degree == 2) {
        matchnodes(v, parv);
        return;
    }

    if ((degree == 3)) {

        if (required_topology == ALL_QUADS) {
            TreeNode *vsib = v->getSibling();
            splitParent(parv, v, vsib);
            return;
        }

        TreeNode *vsib = v->getSibling();
        JNodePtr d0 = v->getDualNode();
        JNodePtr d1 = vsib->getDualNode();
        if (d0->getNumRelations(0) < d1->getNumRelations(0)) {
            matchnodes(v, parv);
            btree->unlinkNode(vsib);
        } else {
            matchnodes(vsib, parv);
            btree->unlinkNode(v);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

TreeNode* JBinaryTreeMatch::getChildofDegreeNParent(TNodeList &levelnodes,
        int nd)
{
    TreeNode *currnode, *parent, *child;

    int ncount;
    TNodeList::const_iterator it;

    for (it = levelnodes.begin(); it != levelnodes.end(); ++it) {
        currnode = *it;
        parent = currnode->getParent();
        if (parent) {
            if (!parent->isMatched()) {
                ncount = 0;
                if (parent->getParent())
                    ncount = 1;
                for (int i = 0; i < parent->getNumChildren(); i++) {
                    child = parent->getChild(i);
                    if (!child->isMatched())
                        ncount++;
                }
                if (ncount == nd)
                    return currnode;
            }
        }
    }

    return nullptr;
}

///////////////////////////////////////////////////////////////////////////////

TreeNode *JBinaryTreeMatch::getNextNode(TNodeList &levelnodes)
{
    TreeNode *currnode = nullptr;

    if (levelnodes.empty())
        return currnode;

    TNodeList::iterator it;
    for (it = levelnodes.begin(); it != levelnodes.end(); ++it) {
        currnode = *it;
        if (currnode->isMatched())
            btree->unlinkNode(currnode);
    }

    it = remove_if(levelnodes.begin(), levelnodes.end(), already_matched);
    levelnodes.erase(it, levelnodes.end());

    TreeNode *child = nullptr;

    // High Priority: parent having degree = 1;
    child = getChildofDegreeNParent(levelnodes, 1);

    if (!child)
        child = getChildofDegreeNParent(levelnodes, 2);

    // Low Priority: parent having degree = 3;
    if (!child)
        child = getChildofDegreeNParent(levelnodes, 3);

    return child;
}

////////////////////////////////////////////////////////////////////////////////

void JBinaryTreeMatch::prunelevel(TNodeList &levelnodes)
{
    while (1) {
        TreeNode *currnode = getNextNode(levelnodes);
        if (currnode == nullptr) break;
        matchnode(currnode);
    }
}

////////////////////////////////////////////////////////////////////////////////

void JBinaryTreeMatch::percolateup()
{
    /*
        steinerNodes.clear();
        steinerFaces.clear();

        int height = btree->getHeight();
        TNodeList levelnodes;
        TNodeList::const_iterator it;

        //Reset all the Matching marks to 0;
        for (int i = 0; i < height; i++) {
            levelnodes = btree->getLevelNodes(height - i - 1);
            TreeNode *currnode;
            for (it = levelnodes.begin(); it != levelnodes.end(); ++it) {
                currnode = *it;
                currnode->setMatchMark(0);
            }
        }

        // Start Prunning the level. At most the root will be unmatched.
        for (int i = 0; i < height; i++) {
            levelnodes = btree->getLevelNodes(height - i - 1);
            prunelevel(levelnodes);
        }

        size_t numfaces = trimesh->getSize(2);
        facematching.reserve(numfaces);

        JFacePtr mateface;
        JNodePtr u, v;
        for (size_t i = 0; i < numfaces; i++) {
            JFacePtr face = trimesh->getFaceAt(i);
            face->getAttribute("DualNode", u);
            u->getAttribute("DualMate", v);
            if (v) {
                if(v > u) {
                    v->getAttribute("PrimalFace", mateface);
                    JFacePair facepair;
                    facepair.first  = face;
                    facepair.second = mateface;
                    facematching.push_back(facepair);
                }
            }
        }

        // If the root is unmatched, bring it down to a leaf and then split the
        // leaf. Do this step after the triangles have been matched.
        TreeNode *root = btree->getRoot();
        if (!root->isMatched()) {
            cout << "Warning: Boundary Triangle modified " << endl;
    #ifdef VERBOSE
            cout << "Warning: Boundary Triangle modified " << endl;
    #endif
            JNodePtr dnode = root->getDualNode();
            JFacePtr rootface = nullptr;
            dnode->getAttribute("PrimalFace", rootface);
            match_boundary_triangle(rootface);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JBinaryTreeMatch::maximum_tree_matching()
{
    // In order to insert any steiner point on the boundary triangle (at the root)
    // We should know which triangles and nodes are on the boundary. Therefore,
    // call this function to set the boundary flags. Building the relationship
    // at this stage is good as even the DualGraph construction require it.

    trimesh->getTopology()->searchBoundary();

#ifdef VERBOSE
    cout << " Creating Dual Graph ... " << endl;
#endif

    JMeshDualGraph dgrapher;
    dgrapher.setMesh(trimesh);
    JMeshPtr dgraph = dgrapher.getGraph();

    if (verbose)
        cout << " Building Binary Tree of Dual Graph ... " << endl;

    btree = new BinaryTree(dgraph);
    btree->build();

#ifdef VERBOSE
    btree->saveAs("btree");
    cout << " Tree Matching ... " << endl;
#endif

    percolateup();

    btree->deleteAll();
    dgraph->deleteAll();

    delete btree;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JBinaryTreeMatch::getQuadMesh(JMeshPtr inmesh, int topo)
{
    JMeshTopologyPtr topology = inmesh->getTopology();
    if (topology->getElementsType(2) != JFace::TRIANGLE) {
        cout << "Warning: Input mesh is not triangular " << endl;
        return nullptr;
    }

    if (!topology->isSimple()) {
        cout << "Warning: Input mesh is not simple, use edmonds' algorithm " << endl;
        return nullptr;
    }

#ifdef DEBUG
    if (inmesh->getNumOfComponents() > 1) {
        cout << "Warning: There are multiple components in the mesh" << endl;
        cout << "         Algorithm works for single component " << endl;
        return nullptr;
    }
    int euler0 = trimesh->getEulerCharacteristic();
    cout << " Input Euler # : " << euler0 << endl;
    cout << inmesh->saveAs( "Check.dat");
#endif

    trimesh = inmesh;

    required_topology = topo;

    trimesh->enumerate(2);
    maxfaceID = trimesh->getSize(2) - 1;

    ///////////////////////////////////////////////////////////////////////////
    // Generate Maximum Matching on a binary tree using Suneeta's Algorithm.
    // If the required topology is set to ALL_QUADS, steiner points( and new
    // faces) will be inserted in the input triangle mesh.  Please note that
    // this implementation doesn't produces "The" optimal soluation as
    // described in the original papers, and doesn't even guarantee that
    // the resulting quadrilaterals will be convex. This along with other
    // topological and geometric optimization are anyhow essential and
    // are carried out during the "Post Processing" step. Therefore, we
    // have sacrifised performance over quality in this implementation.
    // Roughly we can expect that about 4-5% steiner points are inserted in
    // most of the general cases.
    ///////////////////////////////////////////////////////////////////////////
    JStopWatch swatch;
    swatch.start();

    maximum_tree_matching();

    JMeshPtr quadmesh;
//    quadmesh = AllQuadMeshGenerator::collapse_matched_triangles(trimesh, facematching);

    swatch.stop();
    cout << "Info: Tri->Quad Elapsed Time " << swatch.getSeconds() << endl;

    if( quadmesh ) {
        topology = quadmesh->getTopology();
        if (!topology->isSimple())
            cout << "Warning: Quadrilateral Mesh is not simple " << endl;

        topology->getConsistent();
        if (!topology->isConsistent()) {
            cout << "Alas ! Quadrilateral Mesh is still inconsistently oriented: Check manually " << endl;
        }

        quadmesh->enumerate(0);
        quadmesh->enumerate(2);
    } else {
        cout << "Saving the modified triangle mesh: modmesh.xml" << endl;
//      trimesh->saveAs("modmesh.xml");

        cout << "Saving the binary tree " << endl;


    }

    //////////////////////////////////////////////////////////////////////////
    // Since Steiner points may be inserted in the mesh ( and new triangles).
    // Renumber all the nodes and faces for future processing.
    //////////////////////////////////////////////////////////////////////////

    trimesh->enumerate(0);
    trimesh->enumerate(2);

    return quadmesh;
}

///////////////////////////////////////////////////////////////////////////////

void JBinaryTreeMatch::match_tree_walk(BinaryTree *btree, TreeNode *parent)
{
    //
    // Brings all the internal unmatched nodes at the leaf.
    //
    if (parent == nullptr)
        return;
    if (parent->isLeaf())
        return;

    int numChildren = parent->getNumChildren();

    for (int i = 0; i < numChildren; i++) {
        TreeNode *child1 = parent->getChild(i);
        if (!btree->isMatched(parent, child1)) {
            int numGrandChildren = child1->getNumChildren();
            for (int j = 0; j < numGrandChildren; j++) {
                TreeNode *child2 = child1->getChild(j);
                if (btree->isMatched(child1, child2)) {
                    JNodePtr np = parent->getDualNode();
                    assert(np);
                    JNodePtr c1 = child1->getDualNode();
                    assert(c1);
                    JNodePtr c2 = child2->getDualNode();
                    assert(c2);
                    matchnodes(np, c1);
                    c2->setAttribute("DualMate", 0);
                    c2->setStatus(JMeshEntity::ACTIVE);
                    match_tree_walk(btree, child2);
                    return;
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
