#pragma once

#include "Mesh.hpp"

using namespace Jaal;

struct GeodeTemplate
{
    int unitTest();
    int refineHex( const JCellPtr hex, const JFacePtr bot,
                   JNodeSequence &newnodes,
                   JFaceSequence &newfaces,
                   JCellSequence &newcells);
private:
    void refineBottomFace();
    void refineTopFace();
    void refineSideFace();
};
