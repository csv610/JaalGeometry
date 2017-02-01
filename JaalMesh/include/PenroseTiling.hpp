#ifndef PENROSE_H
#define PENROSE_H

#include "Mesh.hpp"

using namespace Jaal;

class PenroseTiling
{
public:
    static JTrianglePtr  getCanonical( int type, double l = 1.0);

    PenroseTiling() {
        goldenRatio = 0.5*(1 + sqrt(5.0));
    }

    JMeshPtr getMesh(int n, double r = 1.0 ) {
        mesh = JMesh::newObject();
        numLevels = n;
        initRadius  = r;
        execute();
        return mesh;
    }


    JMeshPtr mesh;

private:
    int numLevels;
    double initRadius;
    double goldenRatio;

    deque<JFacePtr> faceQ;

    void execute();
    void refine( JFacePtr f);
    void refine36(JFacePtr tri);
    void refine108(JFacePtr tri);

};

#endif
