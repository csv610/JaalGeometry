#pragma once

#include <sstream>

#ifdef USE_IGL
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/comb_cross_field.h>
#include <igl/comb_frame_field.h>
#include <igl/compute_frame_field_bisectors.h>
#include <igl/cross_field_missmatch.h>
#include <igl/cut_mesh_from_singularities.h>
#include <igl/find_cross_field_singularities.h>
#include <igl/local_basis.h>
#include <igl/rotate_vectors.h>
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/copyleft/comiso/miq.h>
#include <igl/readOBJ.h>
#endif

#include "MeshSurfaceVectorField.hpp"
#include "MeshMatrix.hpp"
//#include <qex.h>

class JMixedIntegerQuads : public JMeshSurfaceVectorField
{
public:
    void setMesh( const JMeshPtr &m);

    void genRandomConstraints( int n);

    JMeshPtr getUVMesh();
    JMeshPtr getSeams();

    std::pair<JMeshPtr,JMeshPtr> getCrossField();
    std::pair<JMeshPtr,JMeshPtr> getBisectorField();
    std::pair<JMeshPtr,JMeshPtr> getBisectorCombinedField();

private:
    // Scale for visualizing the fields
    bool extend_arrows = false;

    // Cross field
    Eigen::MatrixXd X1,X2;

// Bisector field
    Eigen::MatrixXd BIS1, BIS2;

// Combed bisector
    Eigen::MatrixXd BIS1_combed, BIS2_combed;

// Per-corner, integer mismatches
    Eigen::MatrixXi MMatch;

// Field singularities
    Eigen::VectorXi isSingularity, singularityIndex;

// Per corner seams
    Eigen::MatrixXi Seams;

// Combed field
    Eigen::MatrixXd X1_combed, X2_combed;

// Global parametrization (with seams)
    Eigen::MatrixXd UV_seams;
    Eigen::MatrixXi FUV_seams;

// Global parametrization
    Eigen::MatrixXd UV;
    Eigen::MatrixXi FUV;

// Create a texture that hides the integer translation in the parametrization
    void line_texture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_B);

//  void addEdges( Eigen::MatrixXd &v);

    void getVecField(int k);
    JMeshPtr extractQuadMesh();
};
