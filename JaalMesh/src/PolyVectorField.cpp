#include "PolyVectorsField.hpp"

using namespace  Eigen;

////////////////////////////////////////////////////////////////////////////////////

// Create a random set of tangent vectors
Eigen::VectorXd random_constraints(const Eigen::VectorXd& b1,
                                   const Eigen::VectorXd& b2, int n)
{
    double rand_factor = 5.0;
    Eigen::VectorXd r(n*3);
    for (unsigned i=0; i<n; ++i)
    {
        double a = (double(rand())/RAND_MAX)*2*M_PI;
        double s = 1 + ((double(rand())/RAND_MAX)) * rand_factor;
        Eigen::Vector3d t = s * (cos(a) * b1 + sin(a) * b2);
        r.block(i*3,0,3,1) = t;
    }
    return r;
}


void JPolyVectorsField :: genRandomConstraints(int nRandom)
{
    b.resize(nRandom);

    int numFaces = mesh->getSize(2);

    for( int i = 0; i < nRandom; i++)
        b[i]    =   JMath::random_value(0, numFaces-1);

    bc.resize(nRandom,3*numVecPerFace);
    for(unsigned i=0; i<b.size(); ++i)
    {
        VectorXd t = random_constraints(B1.row(b(i)),B2.row(b(i)),numVecPerFace);
        bc.row(i) = t;
    }
}

////////////////////////////////////////////////////////////////////////////////////
void JPolyVectorsField :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

#ifdef USE_IGL
    JMeshEigenMatrix mat;
    mat.setMesh(mesh);

    V = mat.getNodeMatrix();
    F = mat.getFaceMatrix();

    // Compute local basis for faces
    igl::local_basis(V,F,B1,B2,B3);

    // Compute face barycenters
    igl::barycenter(V, F, B);

    // Compute scale for visualizing fields
    global_scale =  0.5*igl::avg_edge_length(V, F);

    // Make the example deterministic
    srand(0);
#endif
}

////////////////////////////////////////////////////////////////////////////////////

int JPolyVectorsField :: genField()
{
#ifdef USE_IGL
    Eigen::MatrixXd pvf;
    igl::n_polyvector(V, F, b, bc, pvf);

    vecField = JMesh::newObject();

    for (int n=0; n< numVecPerFace; ++n)
    {
        MatrixXd VF = pvf.block(0,3*n, F.rows(), 3);
        addEdges(vecField, B - global_scale*VF, B + global_scale*VF);
    }
    return 0;
#endif
}
