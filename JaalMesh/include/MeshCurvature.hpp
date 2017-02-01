#pragma once

#include "Mesh.hpp"
#include "MeshMatrix.hpp"
#include "AllTriMeshGenerator.hpp"

#include <Eigen/Core>

#ifdef USE_IGL
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/principal_curvature.h>
#endif

class JMeshCurvature
{
public:
    void setMesh( const JMeshPtr &m);

    void setEvalPos( int p) {
        evalPos = p;
    }

    int  setGaussianCurvature();
    int  setMeanCurvature(int method = 0);
    void setScale( double s) {
        scale = s;
    }

    std::pair<JMeshPtr, JMeshPtr> getCurvatureDirections();

private:
    JMeshPtr mesh;     // Input can have any 2D complex..
    JMeshPtr trimesh;  // We need only the triangle mesh. If the trimesh is
                       // nullptr, then none of the methods will work.
    bool     derived;  // Have we derived a triangle mesh from the input mesh ?
                       // If yes, we will delete it in the destructor.
    double   scale;    // Scaling factor. 
    int      evalPos;  // Where the values are evaluated ? Nodes or faces ?
    JMeshPtr getVectors( const Eigen::MatrixXd &headPoints,
                         const Eigen::MatrixXd &tailPoints);
    void   assignAverage( const JFacePtr &f,  const string &s);
};
