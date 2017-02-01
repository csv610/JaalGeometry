#pragma once

#ifdef USE_IGL
#include <igl/local_basis.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <igl/readOBJ.h>
#include <igl/copyleft/comiso/nrosy.h>
#endif

#include "MeshSurfaceVectorField.hpp"

using namespace Eigen;
using namespace Jaal;

class JNRoSyField : public JMeshSurfaceVectorField
{
public:
    JNRoSyField() {
        numVecPerFace = 4;
    }

    void genRandomConstraints( int n);

    void setMesh( const JMeshPtr &m);
    int  genField();
private:
    double edgelen;
    void representative_to_nrosy(const Eigen::MatrixXd& R, Eigen::MatrixXd& Y);
};


