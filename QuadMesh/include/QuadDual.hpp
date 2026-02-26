#pragma once

#include "Mesh.hpp"

#include "QuadChord.hpp"

class  JQuadDual
{
    static JLogger* logger;
public:
    static int  NoIntersection_Refinement( const JFacePtr &face, JNodeSequence &s1, JFaceSequence &s2);

    JQuadDual( const JMeshPtr &m )
    {
        if( m->getTopology()->getDimension() == 3)
            mesh = m->getTopology()->getSurfaceMesh();
        else
            mesh = m;

        initialize();
    }

    // Start from the given edge and search for its chord.
    JQuadChordPtr getChord(const JEdgePtr &edge);

    // Start from the given faces and return two chord which passes through it.
    vector<JQuadChordPtr> getChords(JEdgeSequence &seq);
    // Start from the given faces and return two chord which passes through it.
    vector<JQuadChordPtr> getColumn(const JFacePtr &quad);

    // Start from the given faces and return two chord which passes through it.
    vector<JQuadChordPtr> getColumns(JFaceSequence  &f);

    // Get all the chord in the mesh.
    vector<JQuadChordPtr> getAllChords();

    vector<JQuadChordPtr> getAllSimpleChords() ;
    vector<JQuadChordPtr> getAllCyclicChords() ;

    int collapse_simple_cycle( const JQuadChordPtr &d);

private:
    size_t nCounter;
    void   initialize();
    JMeshPtr  mesh;
};

