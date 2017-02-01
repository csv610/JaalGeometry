#pragma once

#include "Mesh.hpp"
#include "MeshOptimization.hpp"
#include "MeshRefine.hpp"
#include "SwapEdges.hpp"

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

class JSinglet : public JMeshTopologyOptimization {
public:
    static bool   isSinglet(const JNodePtr &v);
    static void   addAttribute(const JMeshPtr &m, const string &s = "Singlet" );

    JNodeSequence getSinglets();

    int getSize();
    int remove(const JNodePtr &v);
    int removeAll();

    void mergeAll();

private:
    JQuadRefiner quadRefiner;
    JNodeSequence  singlets;
    void  searchSinglets();
    JSwapQuadEdge  edgeSwap;

    int mergeOpposite( const JNodePtr &v);
};

///////////////////////////////////////////////////////////////////////////////
