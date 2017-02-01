#include "MeshUntangle.hpp"
#include "StopWatch.hpp"

////////////////////////////////////////////////////////////////////////

JMeshUntangle :: JMeshUntangle()
{
    energyType = JLocallyInjectiveMap::DIRICHLET_ENERGY;
}

JMeshUntangle :: ~JMeshUntangle()
{
    inflatedMesh.reset();
}

///////////////////////////////////////////////////////////////////////////////////

void JMeshUntangle :: setMesh( const JMeshPtr &m)
{
    inflatedMesh.reset();

    srcMesh = m;
    if( srcMesh == nullptr) return;

    srcMesh->pruneAll();
    srcMesh->enumerate(0);

    srcMesh->getTopology()->getConsistent();

    entityDim = srcMesh->getTopology()->getDimension();

    // If all the elements are positive, we do not have to do anything.
    size_t nCount = countInverted( srcMesh );

    if( nCount == 0) return;
    org_simplicial = 0;

    // If more than half are inverted, then flip the mesh ...
    if( entityDim == 2) {
        if( nCount > 0.5*srcMesh->getSize(2) ) {
            size_t numfaces = srcMesh->getSize(2);
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &f = srcMesh->getFaceAt(i);
                f->reverse();
            }
        }

        // For untangling purposes, we will be working with a triangle
        // mesh so make a trimesh from a quad mesh...
        if( srcMesh->getTopology()->getElementsType(2) == JFace::TRIANGLE) {
            triMesh = srcMesh->deepCopy();
            org_simplicial = 1;
        } else {
            JMeshPtr tmpmesh = srcMesh->deepCopy();
            triMesh = AllTriMeshGenerator::getFromQuadMesh( tmpmesh, 4);
        }
        inflatedMesh = triMesh;
        JMeshGeometry jm(inflatedMesh);
        area[0] = jm.getSurfaceArea();
    }

    if( entityDim == 3) {
        if( nCount > 0.5*srcMesh->getSize(3) ) {
            size_t numcells = srcMesh->getSize(3);
            for( size_t i = 0; i < numcells; i++) {
                const JCellPtr &c = srcMesh->getCellAt(i);
                c->reverse();
            }
        }

        // For untangling purposes, we will be working with a tet
        // mesh so make a tetmesh from a hex  mesh...
        if( srcMesh->getTopology()->getElementsType(3) == JCell::HEXAHEDRON) {
            JMeshPtr tmpmesh = srcMesh->deepCopy();
            AllTetMeshGenerator  alltet;
            tetMesh = alltet.fromHexMesh(tmpmesh);
        } else {
            tetMesh = srcMesh->deepCopy();
            org_simplicial = 1;
        }

        inflatedMesh = tetMesh;
        JMeshGeometry jm(inflatedMesh);
//        volume[0] = jm.getVolume();
    }
    string name = "inflated_" + srcMesh->getName();
    inflatedMesh->setName(name);

    inflatedMesh->getTopology()->searchBoundary();
    size_t numnodes = inflatedMesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = inflatedMesh->getNodeAt(i);
        if( v->isBoundary() ) {
            const Point3D  &p  = v->getXYZCoords();
            v->setAttribute("TargetPos", p);
        }
    }

    inflatedMesh->buildRelations(0,0);
    lap.setMesh(inflatedMesh);
    lap.setNumIterations(100);
    lap.setBoundaryPreserve(1); // We will explicitly change the boundary and
    // Improve the internal nodes.
}

////////////////////////////////////////////////////////////////////////
int JMeshUntangle :: smoothCorner(const JNodePtr &vtx)
{
    JNodeSequence vneighs;
    vtx->getRelations(vtx, vneighs);

    Point3D pmid;
    int nCount = 0;

    pmid[0] = 0.0;
    pmid[1] = 0.0;
    pmid[2] = 0.0;
    for( const JNodePtr &v : vneighs) {
        if( v->isBoundary() ) {
            const Point3D &p = v->getXYZCoords();
            pmid[0] += p[0];
            pmid[1] += p[1];
            pmid[2] += p[2];
            nCount++;
        }
    }
    if( nCount == 0) return 1;
    pmid[0] /= ( double)nCount;
    pmid[1] /= ( double)nCount;
    pmid[2] /= ( double)nCount;
    vtx->setXYZCoords(pmid);
    return 0;
}

////////////////////////////////////////////////////////////////////////
int JMeshUntangle :: smoothBoundary()
{
    if( inflateAll ) {
        size_t numnodes = triMesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = triMesh->getNodeAt(i);
            if( vtx->isBoundary() ) smoothCorner( vtx);
        }
        return 0;
    }

    for( const JNodePtr &vtx : concaveCorners)
        smoothCorner( vtx);

    return 0;
}

////////////////////////////////////////////////////////////////////////
JFaceSequence JMeshUntangle :: getInvertedFaces()
{
    JFaceSequence inverted;
    size_t numfaces = srcMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = srcMesh->getFaceAt(i);
        if( JFaceGeometry::isInverted(face) ) inverted.push_back(face);
    }
    return inverted;
}
////////////////////////////////////////////////////////////////////////
JCellSequence JMeshUntangle :: getInvertedCells()
{
    JCellSequence inverted;

    size_t numcells = srcMesh->getSize(3);
    for( size_t i = 0; i < numcells; i++) {
        const JCellPtr &cell = srcMesh->getCellAt(i);
        if( JCellGeometry::isInverted(cell) ) inverted.push_back(cell);
    }
    return inverted;
}
////////////////////////////////////////////////////////////////////////

size_t JMeshUntangle :: countInverted( const JMeshPtr &m)
{
    if( m == nullptr) return 0;
    concaveCorners.clear();

    size_t nCount = 0;

    int edim = m->getTopology()->getDimension();
    if( edim == 2 ) {
        size_t numfaces = m->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = m->getFaceAt(i);
            if( JFaceGeometry::isInverted(face) ) {
                for( int j = 0; j < face->getSize(0); j++) {
                    const JNodePtr &v = face->getNodeAt(j);
                    if( v->isBoundary() ) {
                        concaveCorners.insert(v);
                    }
                }
                nCount++;
            }
        }
        return nCount;
    }

    if( edim == 3 ) {
        size_t numcells = m->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = m->getCellAt(i);
            if( JCellGeometry::isInverted(cell) ) {
                for( int j = 0; j < cell->getSize(0); j++) {
                    const JNodePtr &v = cell->getNodeAt(j);
                    if( v->isBoundary() ) concaveCorners.insert(v);
                }
                nCount++;
            }
        }
        return nCount;
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////
int JMeshUntangle :: updateSource()
{
    size_t numnodes = srcMesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const Point3D &p = inflatedMesh->getNodeAt(i)->getXYZCoords();
        const JNodePtr &vtx = srcMesh->getNodeAt(i);
        vtx->setXYZCoords(p);
    }
}
////////////////////////////////////////////////////////////////////////

int JMeshUntangle :: startBackProjection()
{
    size_t nCount = countInverted( inflatedMesh );
    if( nCount) {
        cout << "Warning: Inflated mesh contains " << nCount <<  " inverted elements; Backprojection does not start " << endl;
        return 1;
    }

    JMeshGeometry jm(inflatedMesh);
    if( entityDim == 2) area[1] = jm.getSurfaceArea();
//    if( entityDim == 3) volume[1] = jm.getVolume();

    JStopWatch swatch;
    inflatedMesh->deleteNodeAttribute("Constraint");
    JLocallyInjectiveMap lim;

    lim.setMesh(inflatedMesh);
    lim.setEnergyType( energyType );
    lim.setMaxIterations( maxBackProjectionSteps);
    lim.solve();

    double maxDist = getMaxDistance();
    if( maxDist > 1.0E-03) return 1;

    updateSource();

    return 0;
}

////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshUntangle :: getSimplicialMesh()
{
    if( srcMesh == nullptr) return nullptr;

    if( entityDim == 2) return triMesh;
    if( entityDim == 3) return tetMesh;

    return inflatedMesh;
}

////////////////////////////////////////////////////////////////////////
double JMeshUntangle :: getMaxDistance() const
{
    Point3D p0, p1;

    double maxDistance = 0.0;
    size_t numnodes = triMesh->getSize(0);
//#pragma omp parallel fohhhhhhhhhhhh
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = triMesh->getNodeAt(i);
        if( vtx->isActive() ) {
            int err = vtx->getAttribute("TargetPos", p0);
            if( !err ) {
                const Point3D  &p1  = vtx->getXYZCoords();
                double d = JMath::length(p0, p1);
                maxDistance = max(maxDistance, d);
            }
        }
    }
    return maxDistance;
}
////////////////////////////////////////////////////////////////////////

int JMeshUntangle :: startInflation()
{
    if( inflatedMesh == nullptr) return 0;
    cout << "Total number of elements in inflated mesh " << inflatedMesh->getSize(3) << endl;

    int miter = 0;
    int nCount = 0;
    for( int iter = 0; iter < maxInflationSteps; iter++) {
        lap.smoothAll();
        nCount = countInverted( inflatedMesh );
        cout << "#Inverted elements "  << nCount << endl;
        if( nCount == 0) break;
        smoothBoundary();
    }

    if( nCount ) return 1;

    offset.clear();
    size_t numnodes = srcMesh->getSize(0);
    offset.resize(numnodes);
//#pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = inflatedMesh->getNodeAt(i);
        const JNodePtr &v1 = inflatedMesh->getNodeAt(i);
        const Point3D  &p1 = v1->getXYZCoords();
        v0->setXYZCoords(p1);
        offset[i] = JNodeGeometry::getLength(srcMesh->getNodeAt(i), triMesh->getNodeAt(i));
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////
int JMeshUntangle :: optimize()
{
    JMeshNonlinearOptimization mopt;
    mopt.setNumIterations(100);
    mopt.setBoundaryPreserve(1);

    // inflatedMesh is always simplicial
    mopt.setMesh(inflatedMesh);
    mopt.improveQuality();

    updateSource();

    // Original mesh could be any mesh...
    if( !org_simplicial) {
        mopt.setMesh(srcMesh);
        mopt.useSimplicial(0);
        mopt.improveQuality();
    }
}
////////////////////////////////////////////////////////////////////////
