#pragma once

#include "Mesh.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "DijkstraShortestPath.hpp"
#include "boost/range/algorithm/find.hpp"

class JQuadDominant2PureQuadsMesher
{
public:
    JQuadDominant2PureQuadsMesher() { }

    void setMesh(const JMeshPtr &m);
    JMeshPtr getCatmullClarkMesh();

    JFaceSequence getNewStrip();   // Find some strip.
    JFaceSequence getNewStrip(const JFacePtr &f0); // Find a strip passing thru the face.
    JFaceSequence getNewStrip(const JFacePtr &f0, const JFacePtr &f1);
    void remeshStrip();

    vector<JFaceSequence> getAllStrips();
    int remeshAllStrips();

private:
    JMeshPtr mesh;
    JMeshPtr dualGraph;
    boost::scoped_ptr<JDijkstraShortestPath> djk;

    JFaceSet nonQuadsSet;
    JFaceSequence newStrip;

    void searchNonQuads();

    void searchNewStrip();
    void searchNewStrip( const JFacePtr &f0);
    void searchNewStrip( const JFacePtr &f0, const JFacePtr &f1);
    void genDualGraph();
    bool checkStrip();

    void  setOnlyTrianglesAsNonQuads();
    void  mergeTriQuad( const JFacePtr &tri, const JFacePtr &quad);
};
