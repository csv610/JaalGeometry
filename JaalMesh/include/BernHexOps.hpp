#pragma once

#ifndef BERNHEXOPS_H
#define BERNHEXOPS_H

#include "Mesh.hpp"
#include "tfiblend.hpp"
#include "MeshRefine.hpp"
#include <map>

using namespace Jaal;

class JBernHexOps {
public:
    static void getCanonical17( JNodeSequence &newnodes, JCellSequence &newcells);
    static void getCanonical26( JNodeSequence &newnodes, JCellSequence &newcells);
    static void getCanonical314_5( JNodeSequence &newnodes, JCellSequence &newcells);
    static void getCanonical316_5( JNodeSequence &newnodes, JCellSequence &newcells);
    static void getCanonical416_4( JNodeSequence &newnodes, JCellSequence &newcells);
    static void getCanonical415_4( JNodeSequence &newnodes, JCellSequence &newcells);

    JBernHexOps() {
        mesh = nullptr;
        currCellID = 0;
    }

    void setMesh( JMeshPtr m) {
        mesh = m;
        currCellID = 0;
    }

    int searchPattern_1_7( vector<JHexahedron*> &);
    static int Op1_7( const JHexahedronPtr hex,
                      JNodeSequence &newnodes, JCellSequence &newcells);

    int searchPattern_7_1( vector<JHexahedronPtr> &);
    static int Op7_1( const JHexahedron *hex, JHexahedron *newhex);

    int searchPattern_2_6( vector<JHexahedronPtr> &);
    static int Op2_6( const JHexahedronPtr hex1,  const JHexahedronPtr hex2,
                      JNodeSequence &newnodes, JCellSequence &result);

    int search_pattern_6_2( vector<JHexahedronPtr> &);
    static int Opr6_2( JCellSequence &oldcells,
                       JNodeSequence &nodesremoved,
                       JCellSequence &cellsremoved,
                       JNodeSequence &newnodes,
                       JCellSequence &newcells);

    int searchPattern_314_5( vector<JHexahedronPtr> &);
    static int Op314_5( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                        const JHexahedronPtr hex3, JNodeSequence &newnodes,
                        JCellSequence &newcells);

    int searchPattern_316_5( vector<JHexahedronPtr> &);
    static int Op316_5( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                        const JHexahedronPtr hex3, JCellSequence &newcells);

    int searchPattern_415_4( vector<JHexahedronPtr> &);
    static int Op415_4( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                        const JHexahedronPtr hex3, const JHexahedronPtr hex4,
                        JCellSequence &newcells);

    int searchPattern_416_4( vector<JHexahedronPtr> &);
    static int Op416_4( const JHexahedronPtr hex1, const JHexahedronPtr hex2,
                        const JHexahedronPtr hex3, const JHexahedronPtr hex4,
                        JCellSequence &newcells);
private:
    JMeshPtr mesh;
    size_t currCellID;

};

#endif

