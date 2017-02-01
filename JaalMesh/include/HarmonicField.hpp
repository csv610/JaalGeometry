#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#ifdef USE_IGL
#include <igl/matlab/matlabinterface.h>
#endif

#include "Mesh.hpp"
#include "SparseMatrix.hpp"
#include "MeshMatrix.hpp"
#include "StopWatch.hpp"

class JHarmonicField
{
public:
    static const int   CG  = 0;   // Standard Conjugate Gradient method ...
    static const int   CMG = 1;   // Combinatorial Multigrid method:  Yiannis Koutis
    static const int   HSC = 2;   // Dilip Krishnan's HSC method

    JHarmonicField();
    ~JHarmonicField();

    void setMesh( const JMeshPtr &m) {
        mesh = m;
        currStep = 0;
    }

    void setFieldName(const string &s);

    void clearConstraints() { fixedNodes.clear(); }
    void addConstraint( const JNodePtr &a, double v);
    void addConstraints( const JNodeSequence &a, double v);
    JNodeSequence getConstraints() const;

    void setSolver(int a, int niters = 1000, double t = 1.0E-06)
    {   solver = a;
        numIters = niters;
        tol = t;
        currStep = 0;
    }

    int    solveSystem();
    double getNorm()const;

private:
    JMeshPtr mesh;

#ifdef USE_IGL
    Engine  *engine;   // Matlab ....
#endif
    std::string fieldName;

    JNodeSequence fixedNodes;
    JNodeSequence freeNodes; // Are those node where we need to calculate the field.

    Eigen::SparseMatrix<double> LA, LB;

    std::map<JNodePtr, int>  nodeLocalID;

    Eigen::MatrixXd xi, yi, zi, si;
    Eigen::MatrixXd xb, yb, zb, sb;

    int currStep = 0;
    int solver   = CMG;
    int numIters = 1000;
    double tol   = 1.0E-06;
    double maxnorm;

    void initSystem();

    void initCG();
    void initHSC();
    void initCMG();

    int  stdPCG( const string &cmd);
    int  stdPCG( const string &cmd, Eigen::MatrixXd &xb, Eigen::MatrixXd &x);

    int  solveCG();
    int  solveHSC();
    int  solveCMG();
};
