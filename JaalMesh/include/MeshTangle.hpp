#pragma once

#ifndef MESHTANGLE_H
#define MESHTANGLE_H

#include "Mesh.hpp"
#include "RangeSearch.hpp"
#include "PolyBoolean.hpp"

class JMeshTangle
{
public:
    static const int  BRUTE_SEARCH   = 0;
    static const int  BREADTH_SEARCH = 1;
    static int  random_tangle( JMeshPtr mesh, size_t ntangle, vector<size_t> &permute);

    void setMesh( JMeshPtr m) {
        mesh = m;
        rangeSearch.setMesh(m);
        method = BRUTE_SEARCH;
    }

    void  clear() {
        facePairs.clear();
        edgePairs.clear();
        negativeFaces.clear();
        tangledFaces.clear();
    }

    // Search all overlaps. This function must be called when the mesh is modified.
    // All subsequent queries must be following by this call...
    void  searchOverlap();

    // How many negatively oriented faces in the mesh ?
    size_t getNumInvertedElements();

    // Given two faces, give the intesection region
    int  getIntersectionOf( const JFacePtr A, const JFacePtr B, vector<Point2D> &p);

    // How many faces are tangling with the given face ...
    JFaceSequence    getIntersectionOf( const JFacePtr f);

    // Return all the edges in the mesh which are intesected because of tangling ...
    JEdgeSequence    getIntersectEdges()  {
        return intersectEdges;
    }

    // Give all the faces which are negatively oriented in the mesh. All tangling
    // search start from these faces ...
    JFaceSequence    getNegativeFaces()   {
        return negativeFaces;
    }

    // Give all the faces which are overlapped by atleast one face ...
    JFaceSequence    getOverlapFaces()    {
        return tangledFaces;
    }

    // Return all the intersection points of edges ...
    vector<Point2D> getIntersectPoints() {
        return intersectPoints;
    }

    // Give the list of pair of faces which are entangled ....
    vector< pair<JFacePtr,JFacePtr> > getTangledFaces();

private:
    int  method;
    JMeshPtr mesh;
    JRangeSearch rangeSearch;

    struct Intersect
    {
        JEdgePtr edge;
        double  uCoord;
    };

    JPolyBoolean  polyBool;

    set< pair<JFacePtr, JFacePtr> > facePairs;  // Which pair of faces overlap..
    set< pair<JEdgePtr, JEdgePtr> > edgePairs;  // which pair of edges overlap..
    vector<Point2D> intersectPoints;      // Store interecting points of edges...
    JFaceSequence negativeFaces, tangledFaces;
    JEdgeSequence intersectEdges;

    // A brute force is expensive but reliable, It is used to verify the correctness
    // of other optimal algorithms ...
    void brute_overlap_search();
    void breadth_overlap_search();

    void searchNegativeFaces();
    void searchOverlapFaces();
    void searchOverlapEdges();
    void genIntersectPoints();
};

#endif
