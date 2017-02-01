#include "MeshHausdorffDistance.hpp"

vector<double> JMeshHausdorffDistance:: getDistance(int dir)
{
    vector<double> dist;
    if( srcMesh == nullptr || dstMesh == nullptr) return dist;

    JMeshEigenMatrix mat;

    int elemType;
    JMeshPtr asurf = srcMesh->getTopology()->getSurfaceMesh();

    elemType = asurf->getTopology()->getElementsType(2);
    JMeshPtr  trimesh1 = asurf;
    if( elemType == JFace::QUADRILATERAL) {
        AllTriMeshGenerator alltri;
        trimesh1 = alltri.getFromQuadMesh(asurf, 2);
    }

    mat.setMesh(trimesh1);
    Eigen::MatrixXd VA = mat.getNodeMatrix();
    Eigen::MatrixXi FA = mat.getFaceMatrix();

    JMeshPtr bsurf = dstMesh->getTopology()->getSurfaceMesh();
    elemType = bsurf->getTopology()->getElementsType(2);
    JMeshPtr  trimesh2 = bsurf;
    if( elemType == JFace::QUADRILATERAL) {
        AllTriMeshGenerator alltri;
        trimesh2 = alltri.getFromQuadMesh(bsurf, 2);
    }

    mat.setMesh(trimesh2);
    Eigen::MatrixXd VB = mat.getNodeMatrix();
    Eigen::MatrixXi FB = mat.getFaceMatrix();

    Eigen::Matrix<double, Eigen::Dynamic,1> sqrD;
    Eigen::Matrix<size_t, Eigen::Dynamic,1> I;
    Eigen::Matrix<double, Eigen::Dynamic,3> C;

    /*
    if( dir == 1)
        igl::point_mesh_squared_distance(VA,VB,FB,sqrD,I,C);
    else
        igl::point_mesh_squared_distance(VB,VA,FA,sqrD,I,C);

        int nsize = sqrD.size();
        vector<double> dist(nsize);
        for( size_t i = 0; i < nsize; i++)
            dist[i] = sqrt(sqrD[i] );

        return dist;
    */
}
