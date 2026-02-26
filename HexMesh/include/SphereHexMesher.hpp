#pragma once

#include "Mesh.hpp"
#include "AllHexMeshGenerator.hpp"

class JSphereHexMesher
{
public:
    void setMesh( const JMeshPtr &m);
    void setSphere(const JMeshPtr &m);
    void setCubeCells(int *d);

    JMeshPtr  getHexMesh()  {
        return modelHexMesh;
    }
    JMeshPtr  getSphere()   {
        return sphMesh;
    }
    JMeshPtr  getCubeMesh() {
        return refHexMesh;
    }
private:
    JMeshPtr  inMesh;
    JMeshPtr  sphMesh, refHexMesh, modelHexMesh;
    void     projectRefMesh();
    void     genRefMesh();
};
