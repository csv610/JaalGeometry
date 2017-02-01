#include "BinaryTree.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
void BinaryTree::deleteAll()
{
    TNodeMap ::iterator it;
    for (it = tnodemap.begin(); it != tnodemap.end(); ++it)
        delete it->second;
    tnodemap.clear();
}

///////////////////////////////////////////////////////////////////////////////

void BinaryTree::clear()
{
    tnodemap.clear();
    levelnodes.clear();
}

///////////////////////////////////////////////////////////////////////////////
void BinaryTree::relinkAll()
{
    TNodeMap::iterator it;
    for (it = tnodemap.begin(); it != tnodemap.end(); ++it) {
        TreeNode *tnode = it->second;
        tnode->relinkAll();
    }
}

///////////////////////////////////////////////////////////////////////////////

int BinaryTree::getHeight() const
{
    map<int, TNodeList>::const_iterator it;

    int h = 0;
    for (it = levelnodes.begin(); it != levelnodes.end(); ++it)
        h = max(h, it->first);

    return h + 1;
}

///////////////////////////////////////////////////////////////////////////////

size_t BinaryTree::getSize() const
{
    return tnodemap.size();
}
///////////////////////////////////////////////////////////////////////////////

void BinaryTree::bfs_traverse(TreeNode *parent, TNodeList &nextnodes)
{
    JNodePtr dualnode = parent->getDualNode();

    if (dualnode->getVisitBit())
        return;

    JNodeSequence neighs;
    JNode::getRelations( dualnode, neighs );

    int nextlevel = parent->getLevelID() + 1;

    for (size_t i = 0; i < neighs.size(); i++) {
        dualnode = neighs[i];
        if (!dualnode->getVisitBit()) {
            TreeNode *tnode = tnodemap[dualnode];
            if (tnode->getParent() == nullptr) {
                int lid = max(nextlevel, tnode->getLevelID());
                tnode->setLevelID(lid);
                tnode->setParent(parent);
                parent->addChild(tnode);
            }
            nextnodes.push_back(tnode);

        }
    }

    dualnode = parent->getDualNode();
    dualnode->setVisitBit(1);
}

///////////////////////////////////////////////////////////////////////////////

void BinaryTree::dfs_traverse(TreeNode *parent, TNodeList &nextnodes)
{
    JNodePtr dualnode = parent->getDualNode();

    if (dualnode->getVisitBit())
        return;

    vector<TreeNode*> neighQ;
    JNodeSequence neighs;
    JNode::getRelations( dualnode, neighs );

    int nextlevel = parent->getLevelID() + 1;

    int nSize = neighs.size();
    for (int i = 0; i < nSize; i++) {
        dualnode = neighs[i];
        if (!dualnode->getVisitBit()) {
            TreeNode *tnode = tnodemap[dualnode];
            if (tnode->getParent() == nullptr) {
                int lid = max(nextlevel, tnode->getLevelID());
                tnode->setLevelID(lid);
                tnode->setParent(parent);
                parent->addChild(tnode);
            }
            neighQ.push_back(tnode);
        }
    }

    int nQ = neighQ.size();
    for (int i = 0; i < nQ; i++)
        nextnodes.push_front(neighQ[nQ - i - 1]);

    dualnode = parent->getDualNode();
    dualnode->setVisitBit(1);
}

///////////////////////////////////////////////////////////////////////////////

void BinaryTree::dfs_traverse(TreeNode *parent)
{
    TNodeList listnodes;
    listnodes.push_back(parent);

    while (!listnodes.empty()) {
        TreeNode *thisnode = listnodes.front();
        listnodes.pop_front();
        dfs_traverse(thisnode, listnodes);
    }
}

///////////////////////////////////////////////////////////////////////////////

void BinaryTree::bfs_traverse(TreeNode *parent)
{
    TNodeList listnodes;
    listnodes.push_back(parent);

    while (!listnodes.empty()) {
        TreeNode *thisnode = listnodes.front();
        listnodes.pop_front();
        bfs_traverse(thisnode, listnodes);
    }
}

///////////////////////////////////////////////////////////////////////////////

const TNodeList &BinaryTree::getLevelNodes(int level) const
{
    map<int, TNodeList>::const_iterator it;

    it = levelnodes.find(level);
    if (it == levelnodes.end())
        return emptylist;

    return it->second;
}
///////////////////////////////////////////////////////////////////////////////

void BinaryTree::build(TreeNode *r)
{
    assert( dgraph);

    dgraph->buildRelations(0,0);

    int numnodes = dgraph->getSize(0);

    for (int i = 0; i < numnodes; i++) {
        JNodePtr dualnode = dgraph->getNodeAt(i);
        dualnode->setVisitBit(0);
        TreeNode *tnode = new TreeNode(dualnode);
        assert(tnode);
        tnodemap[dualnode] = tnode;
    }

    srand( time(nullptr)); // Initialize random number

    int mindegree = numnodes;
    if (r == nullptr) {
        for (int i = 0; i < numnodes; i++) {
            JNodePtr dualnode = dgraph->getNodeAt(i);
            int degree = dualnode->getNumRelations(0);
            if (degree < mindegree) mindegree = degree;
        }

        int offset = rand()%numnodes;
        assert( offset >=0 && offset < numnodes );
        JNodePtr rootnode;
        for (int i = 0; i < numnodes; i++) {
            int j =  (offset + i)%numnodes;
            assert( j >=0 && j < numnodes );
            JNodePtr dualnode = dgraph->getNodeAt(j);
            int degree = dualnode->getNumRelations(0);
            if( degree == mindegree) {
                rootnode = dualnode;
                break;
            }
        }

        assert(rootnode);
        root = tnodemap[rootnode];
    } else
        root = r;

    assert(root != nullptr);

    root->setLevelID(0);

    if (treetype == BREADTH_FIRST_TREE)
        bfs_traverse( root);

    for (int i = 0; i < numnodes; i++) {
        JNodePtr dualnode = dgraph->getNodeAt(i);
        assert(dualnode->getVisitBit());
        TreeNode *tnode = tnodemap[dualnode];
        int lid = tnode->getLevelID();
        levelnodes[lid].push_back(tnode);
    }

    cout << " Binary Tree Build Completed : #Nodes : " << getSize() << endl;
}

///////////////////////////////////////////////////////////////////////////////

void BinaryTree::saveAs(const string &fname)
{
    string filename = fname + ".dot";
    ofstream ofile(filename.c_str(), ios::out);

    ofile << "digraph G { " << endl;
    ofile << "size = \"8,8\" " << endl;

    TNodeMap::const_iterator it;

    for (it = tnodemap.begin(); it != tnodemap.end(); ++it) {
        TreeNode *tnode = it->second;
        if (tnode->isMatched())
            ofile << tnode->getID() << "[style=filled, fillcolor=green]" << endl;
        else
            ofile << tnode->getID() << "[style=filled, fillcolor=red]" << endl;
    }

    for (it = tnodemap.begin(); it != tnodemap.end(); ++it) {
        TreeNode *tnode = it->second;
        int numChildren = tnode->getNumChildren();
        for (int j = 0; j < numChildren; j++) {
            ofile << tnode->getID() << "->" << tnode->getChild(j)->getID();
            if (isMatched(tnode, tnode->getChild(j)))
                ofile << "[color=blue]";
            else
                ofile << "[color=red]";
            ofile << ";" << endl;
        }

    }

    ofile << " } " << endl;
}

///////////////////////////////////////////////////////////////////////////////

