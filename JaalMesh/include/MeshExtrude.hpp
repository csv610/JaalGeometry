#pragma once

#include "Mesh.hpp"

class JMeshExtrude
{
public:
    void setMesh( const JMeshPtr &m) {
         mesh = m; 
    }

    JMeshPtr getMesh( double distance, int numTimes = 1);
private:
    JMeshPtr mesh;
};
