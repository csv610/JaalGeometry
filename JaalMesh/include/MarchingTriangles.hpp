#pragma once

#include "Mesh.hpp"

using namespace Jaal;

class JMarchingTriangles
{
public:
    JMarchingTriangles() {
        mesh = nullptr;
    }
    void setMesh( JMeshPtr m, const string &a)
    {
        mesh = m;
        attribname = a;
    }

    void setInsideOutside(double val);
    void getContour( double val, JEdgeSequence &eseq);
    void getContours( int nCount, vector<JEdgeSequence> &seq);
    void getContours( const vector<double> &val, vector<JEdgeSequence> &seq);
private:
    JMeshPtr mesh;
    string attribname;
    void checkTriangle(  const JFacePtr face, double scalar, JEdgeSequence &eseq);
};
