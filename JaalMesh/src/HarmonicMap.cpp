#include "HarmonicMap.hpp"

void JHarmonicMap :: init2D(const JMeshPtr &srcMesh)
{
    if( srcMesh->getTopology()->getElementsType(2) == JFace::TRIANGLE) {
        simplicialMesh = srcMesh->deepCopy();
        return;
    }

    JMeshPtr tmpmesh = srcMesh->deepCopy();
    simplicialMesh = AllTriMeshGenerator::getFromQuadMesh( tmpmesh, 4);
}

///////////////////////////////////////////////////////////////////////////////

void JHarmonicMap :: init3D(const JMeshPtr &srcMesh)
{
    if( srcMesh->getTopology()->getElementsType(3) == 8) {
        JMeshPtr tmpmesh =  srcMesh->deepCopy();
        AllTetMeshGenerator alltets;
        simplicialMesh = alltets.fromHexMesh(tmpmesh);
    } else
        simplicialMesh =  srcMesh->deepCopy();
    cout << "Converted to tet mesh " << endl;
}

///////////////////////////////////////////////////////////////////////////////

void JHarmonicMap :: setSource( const JMeshPtr &mesh)
{
    if(mesh == nullptr) return;

    int topdim = mesh->getTopology()->getDimension();
    if( topdim == 2) init2D(mesh);
    if( topdim == 3) init3D(mesh);

    JMeshEigenMatrix mat;
    mat.setMesh(simplicialMesh);
    V = mat.getNodeMatrix();

    if( topdim == 2) F = mat.getFaceMatrix();
    if( topdim == 3) F = mat.getCellMatrix();

    deformedMesh = mesh->deepCopy();
}

///////////////////////////////////////////////////////////////////////////////

void JHarmonicMap :: setTarget( const JMeshPtr &dstMesh)
{
    simplicialMesh->getTopology()->searchBoundary();
    size_t numnodes = dstMesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vsrc = simplicialMesh->getNodeAt(i);
        const JNodePtr &vdst = dstMesh->getNodeAt(i);
        if( vsrc->isBoundary() ) {
            const Point3D &p = vdst->getXYZCoords();
            vsrc->setAttribute("TargetPos", p);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void JHarmonicMap:: setBoundConditions( double t )
{
    cout << "HELLO " << endl;
    size_t numnodes = simplicialMesh->getSize(0);

    vector<int> fixID;
    vector<Point3D> fixPos;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vsrc = simplicialMesh->getNodeAt(i);
        if( vsrc->hasAttribute("TargetPos") ) {
            fixID.push_back(i);
            const Point3D &xyz = vsrc->getXYZCoords();
            fixPos.push_back(xyz);
        }
    }
    int nCount = fixID.size();
    if( nCount < 1) return;

    b.resize(nCount);
    bc.resize(nCount,3);
    for( int i = 0; i < nCount; i++) {
        b[i] = fixID[i];
        bc.coeffRef(i,0) = fixPos[i][0];
        bc.coeffRef(i,1) = fixPos[i][1];
        bc.coeffRef(i,2) = fixPos[i][2];
    }
}
///////////////////////////////////////////////////////////////////////////////
void JHarmonicMap:: solveSystem()
{
#ifdef USE_IGL
    Eigen::MatrixXd D;
    igl::harmonic(V, F, b, bc, order, D);

    size_t numnodes = deformedMesh->getSize(0);
    Point3D xyz;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = deformedMesh->getNodeAt(i);
        xyz[0] = D.coeff(i,0);
        xyz[1] = D.coeff(i,1);
        xyz[2] = D.coeff(i,2);
        vtx->setXYZCoords(xyz);
    }
#endif
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JHarmonicMap:: getDeformedMesh( double t)
{
    setBoundConditions(t);
    solveSystem();
    return deformedMesh;
}
///////////////////////////////////////////////////////////////////////////////
