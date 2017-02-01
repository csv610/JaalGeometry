#pragma once

#include "Mesh.hpp"
#include "basic_math.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Jaal;
using namespace Eigen;

class HeatConduction3D
{
    typedef Eigen::Matrix<double,4,4>  Mat4D;

public:
    HeatConduction3D();

    void setMesh( const JMeshPtr &m) {
        mesh = m;
        matrix_ready = 0;
    }

    void setBoundaryCondition( const JNodePtr &vtx, int , double val) {
        if( vtx == nullptr ) return;
        vtx->setAttribute("Dirichlet", val);
    }

    void setAttributeName( const string &n);

    void setMaxIterations( int n ) {
        numIterations = n;
    }

    int solve();

private:
    JMeshPtr mesh;
    int    numBoundConditions;
    int    numIterations;
    bool   matrix_ready;
    string attribname;
    Eigen::SparseMatrix<double> M, A;
    Eigen::VectorXd  rhs, temperature;

    void genVolMesh();
    void build_global_matrix();
    void get_local_tet_mat(Vec4D &x, Vec4D &y, Vec4D &z, Mat4D &lmat);
    void apply_boundary_conditions();
    void solve_linear_system();
};

//////////////////////////////////////////////////////////////////////////

struct PipeHeatConduction : public HeatConduction3D
{
 private:
      void prepareModel();
};
        


