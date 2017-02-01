#pragma once

#ifdef USE_IGL
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/n_polyvector.h>
#include <igl/integrable_polyvector_fields.h>
#include <igl/local_basis.h>
#include <igl/avg_edge_length.h>
#include <igl/is_border_vertex.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/polyvector_field_matchings.h>
#include <igl/polyvector_field_singularities_from_matchings.h>
#include <igl/polyvector_field_cut_mesh_with_singularities.h>
#include <igl/polyvector_field_comb_from_matchings_and_cuts.h>
#include <igl/polyvector_field_poisson_reconstruction.h>
#include <igl/cut_mesh.h>
#include <igl/slice.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/matlab_format.h>
#endif

#include <iostream>
#include <fstream>

#include "MeshSurfaceVectorField.hpp"

class JIntegrableField : public JMeshSurfaceVectorField
{
public:
    void setMesh( const JMeshPtr &m);

    void genRandomConstraints( int n);
    int  smoothField();
    int  genField();

private:
// Input mesh
    std::vector<bool> V_border;
    std::vector<std::vector<int> > VF, VFi;
    std::vector<std::vector<int> > VV;
    Eigen::MatrixXi TT, TTi;
    Eigen::MatrixXi E, E2F, F2E;

// "Subdivided" mesh obtained by splitting each triangle into 3
// (only needed for display)
    Eigen::MatrixXd Vbs;
    Eigen::MatrixXi Fbs;

// Scale for visualizing the fields
    double global_scale;

// Scale for visualizing textures
    double uv_scale;

// Data for original PolyVector field
    Eigen::MatrixXd two_pv_ori; // field
    Eigen::VectorXi singularities_ori; // singularities
    Eigen::VectorXd curl_ori; // curl per edge
    Eigen::MatrixXi cuts_ori; // cut edges
    Eigen::MatrixXd two_pv_poisson_ori; // field after poisson integration
    Eigen::VectorXf poisson_error_ori; // poisson integration error
    Eigen::MatrixXd scalars_ori;
    Eigen::MatrixXd Vcut_ori;
    Eigen::MatrixXi Fcut_ori;

// Data for curl-free PolyVector field
    Eigen::MatrixXd two_pv; // field
    Eigen::VectorXi singularities; // singularities
    Eigen::VectorXd curl; // curl per edge
    Eigen::MatrixXi cuts; // cut edges
    Eigen::MatrixXd two_pv_poisson; // field after poisson integration
    Eigen::VectorXf poisson_error; // poisson integration error
    Eigen::MatrixXd scalars;
    Eigen::MatrixXd Vcut;
    Eigen::MatrixXi Fcut;

// "constraint level" flag (level=2 indicates that both directions are constrained,
// level = 1 indicates a partially constrained face, i.e. only the first vector will
// be constrained)
    Eigen::VectorXi blevel;

// percentage of constrained faces
    double constraint_percentage = 0.002;

// Random length factor
    double rand_factor = 5;

#ifdef USE_IGL
// The set of parameters for calculating the curl-free fields
    igl::integrable_polyvector_fields_parameters params;

// Solver data (needed for precomputation)
    igl::IntegrableFieldSolverData<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd> ipfdata;
#endif

//texture image
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;

    int display_mode = 1;

    int iter = 0;

    void line_texture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_B);

    void generate_constraints();
    void getCuts(const Eigen::MatrixXi &cuts);
    void getField( const Eigen::MatrixXd &field, const Eigen::RowVector3d &color);
    void colorEdgeMeshFaces(const Eigen::VectorXd &values,
                            const double &minimum, const double &maximum, Eigen::MatrixXd &C);

};
