#pragma once

#ifdef USE_IGL
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/n_polyvector.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#endif

#include <cstdlib>
#include <iostream>
#include <vector>

#include "MeshSurfaceVectorField.hpp"

class JPolyVectorsField : public JMeshSurfaceVectorField
{
public:
    void setMesh( const JMeshPtr &m);
    void genRandomConstraints(int nRandom);
    int  genField();
};
