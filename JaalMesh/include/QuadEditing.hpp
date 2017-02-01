#ifndef QUADEDIT_H
#define QUADEDIT_H

#include "Mesh.hpp"

using namespace Jaal;

class QuadEdit
{
public:
    typedef std::pair<JNodePtr, JNodePtr> NodePair;

    QuadEdit( JMeshPtr m ) {
        mesh = m;
    }

    JEdgePtr getAnyEdge( int n0, int n1);
    int moveEdge( JEdgePtr edge);

    int getAnyPair( int n0, int n1, NodePair &p);
    int movePair( NodePair &p);

private:
    JMeshPtr mesh;
    void initialize();
};
#endif
