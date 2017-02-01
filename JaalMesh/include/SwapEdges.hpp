#pragma once

#ifndef JEDGEFLIP_H
#define JEDGEFLIP_H

#include "Mesh.hpp"
#include "MeshOptimization.hpp"

namespace Jaal {

struct JEdgeSwap : public JMeshTopologyOptimization {
    static const int  NO_RULE_FLIP  = 0;
    static const int  DELAUNAY_FLIP = 1;
    static const int  DEGREE_REDUCING_FLIP = 2;
    static const int  ADVANCE_FRONT_FLIP = 3;

    JEdgeSwap() {}

    ~JEdgeSwap() {}

    void setMesh( const JMeshPtr &m, int r = DELAUNAY_FLIP)
    {
        mesh = m;
        fliprule = r;
    }

    void setCreaseAngle(double a)
    {
        creaseAngle = a;
    }

    void setConstraintEdges(JEdgeSequence &)
    {
//          constraint_edges.add(emesh);
    }

    size_t get_number_of_edges_flipped() const
    {
        return num_edges_flipped;
    }

protected:
    int   fliprule;
    JMeshPtr  mesh;
    JEdgePtr  swapedge;
    double creaseAngle;
    size_t num_edges_flipped;
    JFaceSequence faceNeighs;
};

//////////////////////////////////////////////////////////////////////////////

class JSwapTriEdge : public JEdgeSwap {
public:

    JSwapTriEdge()
    {
        mesh = nullptr;
        creaseAngle = 30.0;
    }

    ~JSwapTriEdge() { }


    bool isSwappable( const JEdgePtr &e, bool check567 = 0);

    // When the swap is successful: the input edge is marked for deletion and
    // a new edge is created. At present, I am creating a new edge, not updating
    // the old edge, which is inefficient, but could ease the burdens of relation
    // manangment.

    int applyAt(const JEdgePtr &e);

    int execute();

private:
    JNodePtr nodes[4];
    bool  isDart() const;
    bool  isSharp() const;
    int   commit();
};

//////////////////////////////////////////////////////////////////////////////////////////////////

class JSwapQuadEdge : public JEdgeSwap
{
    static bool is_topologically_valid_swap(int d1, int d2, int d3, int d4);

public:
    JSwapQuadEdge()
    {
        mesh =   nullptr;
        swapedge = nullptr;
        firstFace = nullptr;
    }

    bool isSwappable( const JEdgePtr &e);

    int applyAt( const JEdgePtr &e, JFacePtr firstface = nullptr);
    int applyAt( const JEdgePtr &e, const JNodePtr &v);
    int applyAt( const JNodePtr &v);

private:
    // Which one of the two faces is the first one. It is needed
    JFacePtr firstFace;
    JNodeSequence boundNodes;  // It is always going to be six nodes...
    JFaceSequence newfaces, edgefaces;

    int applyDegreeRule();
    int applyConcaveRule();   // Swap some concave edge
    int applyBoundRule();     // Swap some boundary edge
    int applySingletRule(const JNodePtr &singlet); // Force creating diagonal at singlet..
    int applyDeficientRule(const JNodePtr &v); // Force creating diagonal at deficient vertex..

    // Get Position of the vertex in the closed chain.
    int getPosOf(const JNodePtr &v) const;
    int buildBoundary();
    int makeDiagonalAt(int pos, bool bound_check = 1);
};
}

///////////////////////////////////////////////////////////////////////////////////

#endif
