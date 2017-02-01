#include "PointLocation.hpp"

JPointLocation :: JPointLocation()
{
    kdTree = nullptr;
    nnIdx  = new ANNidx[1];
    dists  = new ANNdist[1];
}

//////////////////////////////////////////////////////////////////////////
void JPointLocation :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    preprocess();
}
//////////////////////////////////////////////////////////////////////////

void JPointLocation :: preprocess()
{
    if( mesh == nullptr) return;
    if( kdTree) delete kdTree;

    int dim = mesh->getTopology()->getDimension();

    int index = 0;
    Point3D xyz;

    if( dim == 2) {
        mesh->pruneFaces();
        mesh->enumerate(2);
        size_t numFaces = mesh->getSize(2);
        dataPts = annAllocPts( numFaces, 2);
        queryPoint = annAllocPt(2);

        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAvgXYZ(xyz);
                dataPts[index][0] = xyz[0];
                dataPts[index][1] = xyz[1];
                index++;
            }
        }
        kdTree = new ANNkd_tree( dataPts, index, 2);
        return;
    }

    mesh->pruneCells();
    mesh->enumerate(3);

    size_t numCells = mesh->getSize(3);
    dataPts = annAllocPts( numCells, 3);
    queryPoint = annAllocPt(3);

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            cell->getAvgXYZ(xyz);
            dataPts[index][0] = xyz[0];
            dataPts[index][1] = xyz[1];
            dataPts[index][2] = xyz[2];
            index++;
        }
    }
    kdTree = new ANNkd_tree( dataPts, index, 3);
}

//////////////////////////////////////////////////////////////////////////
JCellPtr JPointLocation :: getNextCell( const JCellPtr &cell, const JNodePtr &vtx)
{
    int pos = cell->getPosOf(vtx);
    const JFacePtr &face = cell->getFaceAt(pos);
    JCellSequence neighs;
    JFace::getRelations( face, neighs);
    assert( neighs.size() < 3);
    if( neighs.size() == 1) return nullptr;
    if( neighs[0] == cell ) return neighs[1];
    if( neighs[1] == cell ) return neighs[0];
    return nullptr;
}

//////////////////////////////////////////////////////////////////////////
JCellPtr JPointLocation :: searchCell( const Point3D &qPoint, bool include_boundary )
{
    if( mesh == nullptr ) return NULL;

    queryPoint[0] = qPoint[0];
    queryPoint[1] = qPoint[1];
    queryPoint[2] = qPoint[2];

    double eps = 0.0;
    kdTree->annkSearch( queryPoint, 1, nnIdx, dists, eps);

    int id = nnIdx[0];
    const JCellPtr &seedCell = mesh->getCellAt( id );

    deque<JCellPtr> cellQ;
    JCellSequence    neighs;
    JCellSet         cset;

    cellQ.push_back(seedCell);

    JCellPtr currCell, nextCell;
    Point4D uvw;

    int nCount = 0;
    while(!cellQ.empty() ) {
        currCell = cellQ.front();
        cellQ.pop_front();
        nCount++;
        if( TetGeometry::isInside( currCell, qPoint), include_boundary) {
            return currCell;
        }
        TetGeometry::getBaryCoordinates(currCell, qPoint, uvw);
        cset.insert( currCell );

        for( int i = 0; i < 4; i++) {
            nextCell = getNextCell( currCell, currCell->getNodeAt(i) );
            if( nextCell && (cset.find(nextCell) == cset.end()) )
                if( uvw[i] < 0.0)
                    cellQ.push_front(nextCell);
                else
                    cellQ.push_back(nextCell);
        }
    }
    return nullptr;
}

//////////////////////////////////////////////////////////////////////////
JFacePtr JPointLocation :: getNextFace( const JFacePtr &face, const JNodePtr &vtx)
{
    int pos = face->getPosOf(vtx);
    const JEdgePtr &edge = face->getEdgeAt(pos);
    JFaceSequence neighs;
    JEdge::getRelations( edge, neighs);
    assert( neighs.size() < 3);
    if( neighs.size() == 1) return nullptr;
    if( neighs[0] == face ) return neighs[1];
    if( neighs[1] == face ) return neighs[0];
    return nullptr;
}
//////////////////////////////////////////////////////////////////////////
JFacePtr JPointLocation :: searchFace( const Point3D &qPoint, bool include_boundary )
{
    if( mesh == nullptr ) return NULL;

    queryPoint[0] = qPoint[0];
    queryPoint[1] = qPoint[1];

    double eps = 0.0;
    kdTree->annkSearch( queryPoint, 1, nnIdx, dists, eps);

    int id = nnIdx[0];
    const JFacePtr &seedFace = mesh->getFaceAt( id );

    deque<JFacePtr> faceQ;
    JFaceSequence   neighs;
    JFaceSet        fset;

    faceQ.push_back(seedFace);

    JFacePtr currFace, nextFace;
    Point3D uv;

    int nCount = 0;
    while(!faceQ.empty() ) {
        currFace = faceQ.front();
        faceQ.pop_front();
        nCount++;
        if( JTriGeometry::isInside( currFace, qPoint, include_boundary)) {
            return currFace;
        }
        JTriGeometry::getBaryCoordinates(currFace, qPoint, uv);
        fset.insert( currFace );

        for( int i = 0; i < 3; i++) {
            nextFace = getNextFace( currFace, currFace->getNodeAt(i) );
            if( nextFace && (fset.find(nextFace) == fset.end()) )
                if( uv[i] < 0.0)
                    faceQ.push_front(nextFace);
                else
                    faceQ.push_back(nextFace);
        }
    }

    // Worst case scenario: scan all the faces ...
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &currFace = mesh->getFaceAt(i);
        if( JTriGeometry::isInside( currFace, qPoint, include_boundary)) {
            cout << "FOUND " << endl;
            return currFace;
        }
    }

    return nullptr;
}
