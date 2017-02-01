#include "MeshInterpolation.hpp"

using namespace Eigen;
///////////////////////////////////////////////////////////////////////////////
void JMeshInterpolation :: setSource( const JMeshPtr &m)
{
    srcMesh = m;
    initialized = 0;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshInterpolation :: setTarget( const JMeshPtr &m)
{
    dstMesh = m;
    initialized = 0;
}
///////////////////////////////////////////////////////////////////////////////
int JMeshInterpolation :: init()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolation :: deform2D()
{
    if( srcMesh->getSize(2) != dstMesh->getSize(2) ) {
        cout << "Warning: number of faces on source and target meshes is not same" << endl;
        return;
    }

    if( srcMesh->getSize(0) == 0) {
        cout << "Warning: source mesh is empty" << endl;
        return;
    }

    /*
        if( !initialized) {
            simplicialMesh.reset();
            if( srcMesh->getTopology()->getElementsType(2) == 4)
                simplicialMesh =  srcMesh->deepCopy();
            AllTriMeshGenerator alltri;
        }
        initialized = 1;
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshInterpolation :: deform3D()
{
    cout << "Deformed 3D " << endl;

    if( srcMesh->getSize(3) != dstMesh->getSize(3) ) {
        cout << "Warning: number of triangles are not same on source and target mesh " << endl;
        return;
    }

    if( srcMesh->getSize(0) == 0) {
        cout << "Warning: source mesh is empty" << endl;
        return;
    }

    if( !initialized) {
        simplicialMesh.reset();
        if( srcMesh->getTopology()->getElementsType(3) == 8) {
            JMeshPtr tmpmesh =  srcMesh->deepCopy();
            AllTetMeshGenerator alltets;
            simplicialMesh = alltets.fromHexMesh(tmpmesh);
        } else
            simplicialMesh =  srcMesh->deepCopy();
        simplicialMesh->getTopology()->searchBoundary();

        size_t numnodes = srcMesh->getSize(0);
        #pragma omp parallel for
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vsrc = simplicialMesh->getNodeAt(i);
            const JNodePtr &vdst = dstMesh->getNodeAt(i);
            if( vsrc->isBoundary() ) {
                const Point3D &p = vdst->getXYZCoords();
                vsrc->setAttribute("TargetPos", p);
            }
        }
        deformedMesh = srcMesh->deepCopy();
        initialized = 1;
    }

//    if( method == HARMONIC_MAP)  getHarmonicMap();
}
///////////////////////////////////////////////////////////////////////////////
/*
void JMeshInterpolation :: getHarmonicMap()
{
    Point3D xyz;
    size_t numnodes = srcMesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vsrc = simplicialMesh->getNodeAt(i);
        int err = vsrc->getAttribute("TargetPos", xyz);
        if( !err) vsrc->setXYZCoords(xyz);
    }

    JLaplaceMeshSmoother mlap;
    mlap.setMesh(simplicialMesh);
    mlap.setNumIterations(100);
    mlap.smoothAll();
}
*/
///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshInterpolation :: getInterpolatedMesh( double t)
{
    if( srcMesh == nullptr ) return nullptr;
    if( dstMesh == nullptr ) return nullptr;

    if( srcMesh->getSize(0) != dstMesh->getSize(0) ) {
        cout << "Warning: number of nodes on source and target meshes is not same" << endl;
        return nullptr;
    }

    if( srcMesh->getTopology()->getDimension() == 2) deform2D();
    if( srcMesh->getTopology()->getDimension() == 3) deform3D();

    if( simplicialMesh ) {
        size_t numnodes = srcMesh->getSize(0);
        #pragma omp parallel for
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vsrc = simplicialMesh->getNodeAt(i);
            const JNodePtr &vdst = deformedMesh->getNodeAt(i);
            const Point3D  &xyz = vsrc->getXYZCoords();
            vdst->setXYZCoords( xyz ) ;
        }
    }

    return deformedMesh;
}
///////////////////////////////////////////////////////////////////////////////

int JMeshInterpolation :: updatemesh()
{
    size_t numnodes = srcMesh->getSize(0);

    Point3D xyz;
    int grp = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vsrc = srcMesh->getNodeAt(i);
        const JNodePtr &vdst = dstMesh->getNodeAt(i);
        Point3D p0 = vsrc->getXYZCoords();
        Point3D p1 = vdst->getXYZCoords();
        vsrc->setXYZCoords( p1 ) ;
        vsrc->setAttribute("Constraint", grp);
    }

    /*
        size_t nCount = tetmesh->getGeometry()->getNumInvertedElements();

        cout << "#Inverted Elements: " << nCount << endl;

        if( nCount) {
            cout << "Untangling ..." << endl;
            JMeshNonlinearOptimization mopt;
            mopt.setMesh(tetmesh);
            mopt.untangle();
            nCount = tetmesh->getGeometry()->getNumInvertedElements();
            cout << "#Inverted Elements: " << nCount << endl;
        }
        return nCount;
    */
}

///////////////////////////////////////////////////////////////////////////////
/*
int JMeshInterpolation :: refinemesh()
{
    addPoints.clear();

    double val;
    Point3D xyz;
    JMeshQuality mq;

    size_t numNodes = tetmesh->getSize(0);

    for( int i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = tetmesh->getNodeAt(i);
        if( !vtx->isBoundary() ) {
            xyz = vtx->getXYZCoords();
            addPoints.push_back( xyz );
        }
    }

    size_t numCells = tetmesh->getSize(3);
    for (size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = tetmesh->getCellAt(i);
        if( cell->isActive() ) {
            mq.getJacobian(cell,val);
            if( val < 0.0) {
                cell->getAvgXYZ( xyz);
                addPoints.push_back(xyz);
            }
        }
    }
    JMeshPtr newtetmesh = AllTetMeshGenerator::getQualityMesh(srcMesh, addPoints);
    JMeshIO::saveAs(newtetmesh, "tet.xml");

    tetmesh->deleteCells();
    tetmesh->deleteFaces(JMeshEntity::INTERNAL_ENTITY);
    tetmesh->deleteEdges(JMeshEntity::INTERNAL_ENTITY);
    tetmesh->deleteNodes(JMeshEntity::INTERNAL_ENTITY);

    tetmesh = newtetmesh;
}
*/

///////////////////////////////////////////////////////////////////////////////
