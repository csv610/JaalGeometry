#pragma once

#include "Mesh.hpp"
#include "MeshMatrix.hpp"
#include <Eigen/Core>
#include <armadillo>

#ifdef USE_IGL
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/grad.h>
#include <igl/doublearea.h>
#include <igl/repdiag.h>
#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
#include <iostream>
#endif

class JMeshSpectrum
{
    typedef boost::tuple<int,int, double> MatrixData;
public:
    void setMesh(const JMeshPtr &m) {
        mesh = m;
    }

    void genEigenVectors(int n);

    vector<double> getEigenVector(int i) const;
private:
    JMeshPtr mesh;
    Eigen::MatrixXd EVec;
};
