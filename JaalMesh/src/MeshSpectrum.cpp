#include "MeshSpectrum.hpp"

#ifdef USE_IGL
#include <igl/matlab/matlabinterface.h>
#endif

////////////////////////////////////////////////////////////////////////////////

void JMeshSpectrum :: genEigenVectors(int numVectors)
{
    if( mesh == nullptr) return;

#ifdef USE_IGL

    JMeshLaplaceMatrix lpmat;
    lpmat.setMesh(mesh);
    JGeneralSparseMatrix<double> Ltmp = lpmat.getMatrix();

    JEigenMatrixAdaptor adpt;
    Eigen::SparseMatrix<double> L  = adpt.getMatrix(Ltmp, 1);

    Engine *engine;
    igl::matlab::mlinit(&engine);

    // Send Laplacian matrix to matlab
    igl::matlab::mlsetmatrix(&engine,"L",L);

    ostringstream oss;
    oss << "[EV,~] = eigs(-L," << numVectors << ",'sm')";
    string cmd = oss.str();

    // Extract the first 10 eigenvectors
    igl::matlab::mleval(&engine, cmd);
    igl::matlab::mleval(&engine, oss.str());

    // Retrieve the result
    igl::matlab::mlgetmatrix(&engine,"EV",EVec);

    igl::matlab::mlclose(&engine);
#endif
}

////////////////////////////////////////////////////////////////////////////////

vector<double> JMeshSpectrum :: getEigenVector(int id) const
{
    vector<double> vec;
    Eigen::VectorXd v = EVec.col(id);

    v = v.array() - v.minCoeff();
    v = v.array() / v.maxCoeff();

    size_t nsize = v.size();
    vec.resize(nsize);
    for( size_t i = 0; i < nsize; i++)
        vec[i] = v(i);
    return vec;
}
////////////////////////////////////////////////////////////////////////////////
