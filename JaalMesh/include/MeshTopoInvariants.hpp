#pragma once

#include "Mesh.hpp"
#include "MeshIO.hpp"

class JMeshTopoInvariants
{
public:
    void setMesh( const JMeshPtr &m);

    int getEulerCharacteristic();
    int getNumHoles();
    int getNumHandles();
    int getNumTunnels();

    vector<JMeshPtr> getHoles();
    vector<JEdgeSequence> getHandleLoops();
    vector<JEdgeSequence> getTunnelLoops();
private:
    JMeshPtr mesh;
    void checkInput();
    void getInvariant();
    vector<JEdgeSequence> handleLoops;
    vector<JEdgeSequence> tunnelLoops;
};

