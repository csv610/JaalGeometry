#pragma once

#include "Mesh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "SparseMatrix.hpp"

class JMeshEigenMatrix
{
    typedef Eigen::MatrixXd  NodeMatrix;
    typedef Eigen::MatrixXi  EdgeMatrix;
    typedef Eigen::MatrixXi  FaceMatrix;
    typedef Eigen::MatrixXi  CellMatrix;
    typedef Eigen::MatrixXi  ElemMatrix;

public:
    void setMesh(const JMeshPtr &m);

    NodeMatrix  getNodeMatrix();
    EdgeMatrix  getEdgeMatrix();
    FaceMatrix  getFaceMatrix();
    CellMatrix  getCellMatrix();
    FaceMatrix  getTriMatrix(); // Internall converts elements into triangles.
    FaceMatrix  getTetMatrix(); // Internallt convets the elements into tets.

    ElemMatrix  getElementMatrix();

    JMeshPtr     getMesh(const NodeMatrix &n, const FaceMatrix &f);
private:
    JMeshPtr mesh;
/*
    std::map<JNodePtr, size_t> nodeIDMap;
    void createNodeIDMap();
*/
};

class JMeshAdjacencyMatrix
{
public:
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    JGeneralSparseMatrix<int> getMatrix() const;
private:
    JMeshPtr mesh;
};

class JMeshLaplaceMatrix
{
public:
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }
    JGeneralSparseMatrix<double> getMatrix()const;

    // A complete matrix can be decomposed into two part as
    //      (  A   |   B  )
    //      (  0   |   I  )
    //
    //  Where A consists of unknown values matrix and B contains
    //  boundary condition mtrix therefore
    //        [A] Xu =  Bu - B*Xk
    //  Xu: Unknown values; Xk: Known values ...
    //
    JGeneralSparseMatrix<double> getInnerMatrix()const;           // Returns A Matrix
    JGeneralSparseMatrix<double> getBoundMatrix()const;  // Return B Matrix ..
private:
    JMeshPtr mesh;
};

struct JEigenMatrixAdaptor
{
    Eigen::SparseMatrix<double>  getMatrix(JGeneralSparseMatrix<double> &m, bool clearAfterCopy = 0);
    JGeneralSparseMatrix<double> getMatrix( Eigen::SparseMatrix<double> &m, bool clearAfterCopy = 0);
};
