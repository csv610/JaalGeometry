///////////////////////////////////////////////////////////////////////////////
//  Description: Builds a binary tree from a dual graph of triangualated mesh.
//  Chaman Singh Verma
//  University of Wisconsin Madison, USA
//  Date 15th Jan 2010
//
//  License: Free to distribute and modify.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include "Mesh.hpp"
#include "MeshDualGraph.hpp"

namespace Jaal {

class TreeNode;

typedef std::list<TreeNode*> TNodeList;

////////////////////////////////////////////////////////////////////////////////

class TreeNode {
public:
    TreeNode()
    {}

    TreeNode( JNodePtr n)
    {
        dualNode = n;
        levelID = -1;
        parent = nullptr;
    }

    void setLevelID(int l)
    {
        levelID = l;
    }

    int getLevelID() const
    {
        return levelID;
    }

    bool isMatched() const
    {
        return dualNode->isRemoved();
    }

    void setMatchMark(char r)
    {
        dualNode->setStatus(r);
    }

    bool isRoot() const
    {
        if (parent == nullptr)
            return 1;
        return 0;
    }

    bool isLeaf() const
    {
        if (getNumChildren() == 0)
            return 1;
        return 0;
    }

    int getDegree() const
    {
        int ncount = 0;
        if (parent)
            ncount = 1;
        ncount += getNumChildren();
        return ncount;
    }

    TreeNode * getSibling() const
    {
        if (parent) {
            for (int i = 0; i < parent->getNumChildren(); i++) {
                TreeNode *child = parent->getChild(i);
                if (child != this)
                    return child;
            }
        }
        return nullptr;
    }

    void setParent(TreeNode * p)
    {
        parent = p;
    }

    TreeNode * getParent() const
    {
        return parent;
    }

    int getNumChildren() const
    {
        int ncount = 0;
        for (size_t i = 0; i < children.size(); i++)
            if (children[i].active)
                ncount++;
        return ncount;
    }

    void addChild(TreeNode * c)
    {
        ActiveNode activenode(c);
        if (find(children.begin(), children.end(), activenode) == children.end())
            children.push_back(activenode);
    }

    void removeChild(TreeNode * child)
    {
        ActiveNode activenode(child);
        vector<ActiveNode>::iterator it;
        it = remove(children.begin(), children.end(), activenode);
        children.erase(it, children.end());
    }

    void unlinkChild(TreeNode * child)
    {
        for (size_t i = 0; i < children.size(); i++) {
            if (children[i].node == child) {
                children[i].active = 0;
                return;
            }
        }
    }

    void relinkChild(TreeNode * child)
    {
        for (size_t i = 0; i < children.size(); i++) {
            if (children[i].node == child) {
                children[i].active = 1;
                return;
            }
        }
    }

    void relinkAll()
    {
        for (size_t i = 0; i < children.size(); i++)
            children[i].active = 1;
    }

    TreeNode * getChild(int cid) const
    {
        int index = -1;
        for (size_t i = 0; i < children.size(); i++) {
            if (children[i].active) {
                index++;
                if (index == cid)
                    return children[i].node;
            }
        }
        return nullptr;
    }

    int getID() const
    {
        return dualNode->getID();
    }

    JNodePtr getDualNode() const
    {
        return dualNode;
    }

private:
    struct ActiveNode {
        ActiveNode(TreeNode *n)
        {
            active = 1;
            node = n;
        }
        bool active;
        bool operator ==(const ActiveNode &rhs) const
        {
            return node == rhs.node;
        }
        TreeNode *node;
    };
    int levelID;
    JNodePtr dualNode;
    TreeNode* parent;

    vector<ActiveNode> children;
};

////////////////////////////////////////////////////////////////////////////////

class BinaryTree {
public:

    const static int BREADTH_FIRST_TREE = 0;
    const static int DEPTH_FIRST_TREE   = 1;

    BinaryTree(JMeshPtr g)
    {
        dgraph = g;
        treetype = BREADTH_FIRST_TREE;
    }

    ~BinaryTree()
    {
        clear();
    }

    void setTreeType(int t)
    {
        treetype = t;
    }

    void build(TreeNode *r = nullptr);

    TreeNode* getRoot() const
    {
        return root;
    }

    // Remove a given node from the tree. Unlink child and parents.
    //
    void removeNode(TreeNode *node)
    {
        if (node) {
            TreeNode *parv = node->getParent();
            if (parv)
                parv->removeChild(node);
        }
    }

    void addNode(TreeNode *tnode)
    {
        tnodemap[tnode->getDualNode()] = tnode;
    }

    void unlinkNode(TreeNode *node)
    {
        if (node) {
            TreeNode *parv = node->getParent();
            if (parv)
                parv->unlinkChild(node);
        }
    }

    bool isMatched(const TreeNode *u, const TreeNode *v) const
    {
        JNodePtr du = u->getDualNode();
        JNodePtr dv = v->getDualNode();
        JNodePtr umate, vmate;
        du->getAttribute("DualMate", umate);
        dv->getAttribute("DualMate", vmate);
        if (umate == dv && vmate == du) return 1;
        return 0;
    }

    // Give nodes at a given level. Root is at level = 0,
    const TNodeList &getLevelNodes(int l) const;

    // Total Height of the tree...
    int getHeight() const;

    // How many nodes in the tree.
    size_t getSize() const;

    // Clear all the nodes created in the tree.
    void clear();

    void deleteAll();

    void relinkAll();

    // Save the tree in GraphViz data format...
    void saveAs(const string &s);

private:
    JMeshPtr dgraph;
    TreeNode *root;
    TNodeList emptylist;
    int treetype;

    typedef std::map<JNodePtr, TreeNode*> TNodeMap;
    TNodeMap  tnodemap;

    std::map<int, TNodeList> levelnodes;

    void bfs_traverse(TreeNode *parent, TNodeList &nextnodes);
    void bfs_traverse(TreeNode *parent);

    void dfs_traverse(TreeNode *parent);
    void dfs_traverse(TreeNode *parent, TNodeList &nextnodes);
};

} // namespace Jaal

