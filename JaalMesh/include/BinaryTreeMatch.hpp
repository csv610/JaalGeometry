////////////////////////////////////////////////////////////////////////////////
//            Jaal:: Triangle-to-Quad Transformation Code
//            ********************************************
//  Description:
//  Given a triangulated mesh, convert into All-Quads or Quad-Dominated Mesh.
//  using Maximum Tree Matching Algorithm.
//
//  Maximum Cardinality using Edmond's Algorithm may give perfect matching
//  but the algorithm is too expensive for large dataset and often perfect
//  matching is not required.
//
//  If the user requires All-Quad Meshing then all the unmatched triangles
//  using Tree Matching Algorithms are split and steiner points are created.
//  In most of the cases, number of Steiner points are less than 5% of the
//  total triangles present in the original mesh.
//
//  Chaman Singh Verma
//  University of Wisconsin Madison, USA,
//  Date: 15th Jan 2010.
//
//  License: Free to download and modify.

//  For more information, please follow the link
//
//  http://pages.cs.wisc.edu/~csverma/CS899_09/JQuad.html
//
////////////////////////////////////////////////////////////////////////////////

#ifndef Tri2Quad_H
#define Tri2Quad_H

#include "Mesh.hpp"
#include "MeshDualGraph.hpp"
#include "BinaryTree.hpp"
#include "cycle.hpp"         // performance counter.

namespace Jaal {

class JBinaryTreeMatch {
public:
    const static int ALL_QUADS = 0;
    const static int QUAD_DOMINATED = 1;

    JBinaryTreeMatch()
    {
        trimesh = nullptr;
        btree = nullptr;
        verbose = 1;
        required_topology = ALL_QUADS;
    }

    const vector<JFacePair> &getMaximumDualMatching();

    JMeshPtr getQuadMesh(JMeshPtr tmesh, int topo = ALL_QUADS);

    int   getQuadMesh( vector<double> &nodes, vector<size_t> &tmesh, vector<size_t> &qmesh);

    JNodeSequence getSteinerNodes()  const;
    JFaceSequence getInsertedFaces() const;
    JFaceSequence getModifiedFaces() const;

private:
    JMeshPtr trimesh; // Input mesh.

    struct LVertex : public JNode {
        LVertex( JNodePtr v )
        {
            vertex = v;
        }
        JNodePtr vertex;
        JNodePtr mate;
        JFacePtr dual;
    };

    JFaceSequence steinerFaces, modifiedFaces;
    JNodeSequence steinerNodes;

    BinaryTree *btree;

    int required_topology;
    bool verbose;
    size_t maxfaceID;

    vector<JFacePair> facematching;

    void splitParent(JFacePtr parent, JFacePtr child1, JFacePtr child2);
    void splitParent(TreeNode *p, TreeNode *c1, TreeNode *c2);

    int match_boundary_triangle(JFacePtr face);

    void percolateup();

    void matchnode(TreeNode *v);
    void matchnodes(TreeNode *child, TreeNode *parent);
    void matchnodes(JNodePtr child, JNodePtr parent);

    TreeNode* getChildofDegreeNParent(TNodeList &levelnodes, int nd);

    TreeNode *getNextNode(TNodeList &levelnodes);
    void prunelevel(TNodeList &levelnodes);
    void maximum_tree_matching();

    void match_tree_walk(BinaryTree *tree, TreeNode *u);
};

////////////////////////////////////////////////////////////////////////////////

bool has_same_dual(const TreeNode *nd1, const TreeNode *nd2);

inline
bool already_matched(const TreeNode *node)
{
    return node->isMatched();
}

////////////////////////////////////////////////////////////////////////////////

inline
void JBinaryTreeMatch::matchnodes(JNodePtr child, JNodePtr parent)
{
    child->setAttribute("DualMate", parent);
    parent->setAttribute("DualMate", child);
    child->setStatus(JMeshEntity::REMOVE);
    parent->setStatus(JMeshEntity::REMOVE);
}

///////////////////////////////////////////////////////////////////////////////////

inline
void JBinaryTreeMatch::matchnodes(TreeNode *child, TreeNode *parent)
{
    if (parent->isMatched() && !child->isMatched()) {
        cout << "Warning: parent is already matched " << endl;
    }

    if (!child->isMatched() && !parent->isMatched())
        matchnodes(child->getDualNode(), parent->getDualNode());

    btree->unlinkNode(child);
    btree->unlinkNode(parent);
}

} // namespace Jaal


#endif
