#pragma once
#include "Mesh.hpp"
#include "MeshQuality.hpp"
#include "MeshOptimization.hpp"

class JQuadChord
{
    static const int CLOSED_CHORD   = 0;
    static const int BOUNDARY_CHORD = 1;

public:
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    void clear();
    void setSeed( const JEdgePtr &edge);

    JEdgeSequence getEdges() const {
        return chordEdges;
    }
    JFaceSequence getFaces() const {
        return chordFaces;
    }

    JFaceSequence getSelfIntersectingFaces() const
    {
        return selfIntersectingFaces;
    }

    bool isSelfIntersecting() const
    {
        if( selfIntersectingFaces.empty() ) return 0;
        return 0;
    }

    bool isCyclic() const ;

    void setID( int id );


private:
    JMeshPtr   mesh;

    JEdgePtr     seedEdge;
    JFaceSequence selfIntersectingFaces;
    JEdgeSequence chordEdges;
    JFaceSequence chordFaces;
    int  chordType;

    void expand( const JEdgePtr &edge, const JFacePtr &face);
};
typedef boost::shared_ptr<JQuadChord> JQuadChordPtr;

///////////////////////////////////////////////////////////////////////////////

class JDiceQuadChord
{
public:
    void setMesh( const JMeshPtr &m);
    void setCompleteDice( bool p) {
        complete_dice = p;
    }

    void setChord( const JQuadChordPtr &c);

    JFaceSequence getNewFaces() const {
        return newFaces;
    }
private:
    JMeshPtr  mesh;
    JQuadChordPtr   chord;
    bool      complete_dice = 1;
    double    aspectRatio;
    JFaceSequence newFaces;

    void getShapeFuncs( double r, double s, double phi[4] );
    void getEdgeNodes(const JNodePtr &v0, const JNodePtr &v1, JNodePtr &v2, JNodePtr &v3);
    JNodePtr getNodes(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2,
                      const JNodePtr &v3, double r, double s);

    int  getStartIndex();
    void verifyChord();
    void divide( const JEdgePtr &e);
    void edge1( const JFacePtr &f);
    void edge2( const JFacePtr &f);
    void edge3( const JFacePtr &f);
    void edge4(const JFacePtr &f);
};

///////////////////////////////////////////////////////////////////////////////
