#pragma once

#ifndef OBBTREE_H
#define OBBTREE_H

#include "Mesh.hpp"

class JOBBTree
{
    struct OBBNode
    {
        int   id;
        int   depth;
        JNodeSequence   nodes;
        JHexahedronPtr hex;
        boost::shared_ptr<OBBNode> lChild, rChild;
    };

    typedef boost::shared_ptr<OBBNode> NodePtr;
public:
    void setMaxDepth( int d) {
        maxDepth = d;
    }

    JCellSequence getBoxes( const JNodeSequence &nodes);

private:
    int maxDepth;
    JCellSequence  boxes;
    boost::shared_ptr<OBBNode> root;
    int   nCount;

    bool onLeftSide( const vector<double> &plane, const Point3D &query);

    int  getMaxLengthSide( const JHexahedronPtr &h);

    int  getPlane( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, vector<double> &plane);
    void splitHex( const NodePtr &parent);
    int  splitHex( const JHexahedronPtr &parent, int side, JHexahedronPtr &lChild, JHexahedronPtr &rChild, vector<double> &plane);
    int  xsplitHex( const JHexahedronPtr &parent, JHexahedronPtr &lChild, JHexahedronPtr &rChild, vector<double> &plane);
    int  ysplitHex( const JHexahedronPtr &parent, JHexahedronPtr &lChild, JHexahedronPtr &rChild, vector<double> &plane);
    int  zsplitHex( const JHexahedronPtr &parent, JHexahedronPtr &lChild, JHexahedronPtr &rChild, vector<double> &plane);
};

#endif
