#pragma once

#include "Mesh.hpp"

class JCurveShorteningFlow
{
public:
    JCurveShorteningFlow() {};
    ~JCurveShorteningFlow() {};

    void setMesh(const JMeshPtr &m);
    void setPreserveArea( bool a) {
        preserveArea = a;
    }
    void performOneStep();
    double getCurveLength() const;
private:
    JMeshPtr mesh;
    bool  preserveArea = 0;
    vector<JEdgeSequence> boundLoops;
    vector<JNodeSequence> boundNodes;
};
