#include "NRoSyField.hpp"

using namespace Jaal;

#ifdef USE_IGL
using namespace igl::copyleft::comiso;
#endif


/////////////////////////////////////////////////////////////////////

void JNRoSyField::setMesh( const JMeshPtr &m)
{
#ifdef USE_IGL
    mesh = m;
    if( mesh == nullptr) return;

    JMeshEigenMatrix mat;
    mat.setMesh(mesh);

    V = mat.getNodeMatrix();
    F = mat.getFaceMatrix();

    edgelen = 0.5*igl::avg_edge_length(V, F);
    igl::barycenter(V,F,B);

    JEdgeSequence bedges;
    mesh->getTopology()->getBoundary(bedges);

    if( bedges.empty() ) return;

    int nbound = bedges.size();
    b.resize(nbound);
    bc.resize(nbound,3);

    JFaceSequence bfaces;
    Vec3D vec;
    int index = 0;
    for( const JEdgePtr &edge : bedges) {
        JEdge::getRelations(edge, bfaces);
        assert( bfaces.size() == 1);
        const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
        JMath::unit_vector(p1,p0, vec);
        b[index] =  bfaces[0]->getID();
        bc.coeffRef(index,0) = vec[0];
        bc.coeffRef(index,1) = vec[1];
        bc.coeffRef(index,2) = vec[2];
        index++;
    }
#endif
}
/////////////////////////////////////////////////////////////////////

void JNRoSyField::genRandomConstraints( int nrand)
{
    if( mesh == nullptr) return;
    assert( nrand >= 1);

    b.resize(nrand);
    bc.resize(nrand,3);

    int numFaces = mesh->getSize(2);

    Point3D p0, p1;
    for( int i = 0; i < nrand; i++) {
        b[i]    =   JMath::random_value(0, numFaces-1);
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAvgXYZ(p0);
        p1 = face->getNodeAt(0)->getXYZCoords();
        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dz = p1[2] - p0[2];
        double dl = sqrt(dx*dx + dy*dy + dz*dz);
        bc.coeffRef(i,0) =  dx/dl;
        bc.coeffRef(i,1) =  dy/dl;
        bc.coeffRef(i,2) =  dz/dl;
    }
}

/////////////////////////////////////////////////////////////////////
// Converts a representative vector per face in the full set of vectors that describe
// an N-RoSy field
/////////////////////////////////////////////////////////////////////

void JNRoSyField :: representative_to_nrosy( const Eigen::MatrixXd& R, Eigen::MatrixXd& Y)
{
#ifdef USE_IGL
    MatrixXd B1, B2, B3;
    igl::local_basis(V,F,B1,B2,B3);

    Y.resize(F.rows()*numVecPerFace,3);
    for (unsigned i=0; i<F.rows(); ++i)
    {
        double x = R.row(i) * B1.row(i).transpose();
        double y = R.row(i) * B2.row(i).transpose();
        double angle = atan2(y,x);

        for (unsigned j=0; j< numVecPerFace; ++j)
        {
            double anglej = angle + 2*M_PI*double(j)/double(numVecPerFace);
            double xj = cos(anglej);
            double yj = sin(anglej);
            Y.row(i*numVecPerFace+j) = xj * B1.row(i) + yj * B2.row(i);
        }
    }
#endif
}

/////////////////////////////////////////////////////////////////////

int JNRoSyField :: genField()
{
#ifdef USE_IGL
    MatrixXd R;
    VectorXd S;

    if( b.rows() == 0) genRandomConstraints(1);
        
    nrosy(V,F,b,bc,VectorXi(),VectorXd(),MatrixXd(), numVecPerFace, 0.5, R, S);

    MatrixXd Y;
    representative_to_nrosy(R, Y);

    MatrixXd Be(B.rows()*numVecPerFace,3);
    for(unsigned i=0; i<B.rows(); ++i)
        for(unsigned j=0; j< numVecPerFace; ++j)
            Be.row(i*numVecPerFace+j) = B.row(i);

    addEdges(Be, Be+ Y*edgelen);

    positiveSingularNodes.clear();
    negativeSingularNodes.clear();
    for (unsigned i=0; i<S.size(); ++i)
    {
        if(S(i) < -0.001) positiveSingularNodes.push_back( mesh->getNodeAt(i) );
        if(S(i) > +0.001) negativeSingularNodes.push_back( mesh->getNodeAt(i) );
    }

    constrainedFaces.clear();
    for (size_t i=0; i< b.size(); ++i)
        constrainedFaces.push_back( mesh->getFaceAt(b(i)));

    return 0;
#endif
}
