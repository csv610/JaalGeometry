#include "MeshSurfaceVectorField.hpp"

#ifdef USE_IGL
#include <igl/readDMAT.h>
#endif

////////////////////////////////////////////////////////////////////////////////////

bool JMeshSurfaceVectorField :: isPlanar( const JFacePtr &face, const Vec3D &query)
{
    Vec3D  normal  = JFaceGeometry::getNormal(face);
    double dotprod = JMath::dot_product(normal, query);
    if( fabs(dotprod) < 1.0E-06) return 1;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////////

Vec3F JMeshSurfaceVectorField :: getRandomVector( const JFacePtr &face)
{
    Point3D  pCenter;
    face->getAvgXYZ( pCenter);

    double vmin = 0.0;
    double vmax=  1.0;
    double t = JMath::random_value(vmin, vmax);

    int  imin  = 0;
    int  imax  = face->getSize(0)- 1;
    double edgenum = JMath::random_value(imin, imax);

    const JEdgePtr &edge = face->getEdgeAt(edgenum);

    const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
    const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();

    Point3D pnew;
    pnew[0] = (1-t)*p0[0] + t*p1[0];
    pnew[1] = (1-t)*p0[1] + t*p1[1];
    pnew[2] = (1-t)*p0[2] + t*p1[2];

    double dx = pnew[0] - pCenter[0];
    double dy = pnew[1] - pCenter[1];
    double dz = pnew[2] - pCenter[2];
    double dl = sqrt(dx*dx + dy*dy + dz*dz);

    Vec3F av;
    av[0] = dx/dl;
    av[1] = dy/dl;
    av[2] = dz/dl;
    return av;
}

////////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd JMeshSurfaceVectorField :: genRandomVector(const Eigen::VectorXd& b1,
        const Eigen::VectorXd& b2, int n)
{
// Create a random set of tangent vectors
    Eigen::VectorXd r(n*3);
    for (unsigned i=0; i< n; ++i)
    {
        double a = (double(rand())/RAND_MAX)*2*M_PI;
        double s = 1 + ((double(rand())/RAND_MAX)) * rand_factor;
        Eigen::Vector3d t = s * (cos(a) * b1 + sin(a) * b2);
        r.block(i*3,0,3,1) = t;
    }
    return r;
}

////////////////////////////////////////////////////////////////////////////////////
int JMeshSurfaceVectorField :: readConstraints( const string &file1, const string &file2)
{
#ifdef USE_IGL
    igl::readDMAT(file1, b);
    igl::readDMAT(file2, bc);
    return 0;
#endif
}
////////////////////////////////////////////////////////////////////////////////////
JFaceSequence JMeshSurfaceVectorField :: getConstrainedFaces()
{
    constrainedFaces.clear();
    if( mesh == nullptr) return constrainedFaces;

    for( size_t i = 0; i < b.size(); i++)
        constrainedFaces.push_back( mesh->getFaceAt( b(i) ));
    return constrainedFaces;
}
////////////////////////////////////////////////////////////////////////////////////

int JMeshSurfaceVectorField :: readVecField()
{
    ifstream infile("vecfield.dat", ios::in);
    if( infile.fail() ) return 1;

    size_t numfaces, N;
    infile >> numfaces >> N;

    /*
    vecField = JMesh::newObject();

    Point3D p0, p1;
    for( int i = 0; i < numfaces; i++) {
        infile >> p0[0] >> p0[1] >> p0[2];
        JNodePtr v0 =  JNode::newObject();
        v0->setXYZCoords(p0);
        vecField->addObject(v0);
        for( int j = 0; j < N; j++) {
            infile >> p1[0] >> p1[1] >> p1[2];
            JNodePtr v1 =  JNode::newObject();
            v1->setXYZCoords(p1);
            JEdgePtr edge = JEdge::newObject(v0,v1);
            vecField->addObject(v1);
            vecField->addObject(edge);
        }
    }

        int numSingular, id;
        double index;
        infile >> numSingular;
        if( numSingular ) {
            singularNodes.resize(numSingular);
            for( int i = 0; i < numSingular; i++) {
                infile >> id >> index;
                singularNodes[i] = mesh->getNodeAt(id);
            }
        }
    */

    return 0;
}

/////////////////////////////////////////////////////////////////////
int JMeshSurfaceVectorField :: saveConstraints()
{
    size_t numfaces = mesh->getSize(2);
    size_t nfixed = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->isActive() && f->hasAttribute("ConstraintVector") )
            nfixed++;
    }

    if( nfixed  == 0 ) {
        cout << "Warning: At least one face with vector direction must be specified " << endl;
        return 1;
    }

    ofstream ofile("fconstraints.dat", ios::out);
    if( ofile.fail() ) return 1;

    ofile << nfixed << endl;

    JEdgePtr edge;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        if( f->isActive() ) {
            int err = f->getAttribute("ConstraintVector", edge);
            if( !err) {
                const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
                const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
                double dx = p1[0] - p0[0];
                double dy = p1[1] - p0[1];
                double dz = p1[2] - p0[2];
                double dl = sqrt(dx*dx + dy*dy + dz*dz);
                ofile << f->getID() <<  " " << dx/dl << " " << dy/dl << " " << dz/dl << endl;

            }
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////

void JMeshSurfaceVectorField :: addEdges(const MatrixXd &tail, const MatrixXd &head)
{
    int numfaces = tail.rows()/numVecPerFace;

    vecFields.resize(numVecPerFace);

    for( int j = 0; j < numVecPerFace; j++)
        vecFields[j] = JMesh::newObject();

    Point3D  p0, p1;
    size_t index = 0;
    for( size_t i = 0; i < numfaces; i++) {
        for( int j = 0; j < numVecPerFace; j++) {
            JNodePtr v0 = JNode::newObject();
            p0[0]       = tail.coeff(index,0);
            p0[1]       = tail.coeff(index,1);
            p0[2]       = tail.coeff(index,2);
            v0->setXYZCoords(p0);
            vecFields[j]->addObject(v0);

            JNodePtr v1 = JNode::newObject();
            p1[0]       = head.coeff(index,0);
            p1[1]       = head.coeff(index,1);
            p1[2]       = head.coeff(index,2);
            v1->setXYZCoords(p1);
            vecFields[j]->addObject(v1);
            index++;
            JEdgePtr newedge = JEdge::newObject(v0,v1);
            vecFields[j]->addObject(newedge);
        }
    }
}

void JMeshSurfaceVectorField :: addEdges(const JMeshPtr &vmesh, const MatrixXd &tail, const MatrixXd &head)
{
    int nSize = tail.rows();

    Point3D  p0, p1;
    for( size_t i = 0; i < nSize; i++) {
        JNodePtr v0 = JNode::newObject();
        p0[0]       = tail.coeff(i,0);
        p0[1]       = tail.coeff(i,1);
        p0[2]       = tail.coeff(i,2);
        v0->setXYZCoords(p0);
        vmesh->addObject(v0);

        JNodePtr v1 = JNode::newObject();
        p1[0]       = head.coeff(i,0);
        p1[1]       = head.coeff(i,1);
        p1[2]       = head.coeff(i,2);
        v1->setXYZCoords(p1);
        vmesh->addObject(v1);
        JEdgePtr newedge = JEdge::newObject(v0,v1);
        vmesh->addObject(newedge);
    }
}
////////////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshSurfaceVectorField :: getConstrainedField()
{
}
////////////////////////////////////////////////////////////////////////////////////
