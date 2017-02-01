#include "MeshMatrix.hpp"

/////////////////////////////////////////////////////////////////////////

void JMeshEigenMatrix :: setMesh( const JMeshPtr &m)
{
    mesh = m;
}

/////////////////////////////////////////////////////////////////////////
JMeshEigenMatrix::NodeMatrix JMeshEigenMatrix :: getNodeMatrix()
{
    JMeshEigenMatrix::NodeMatrix nodes;

    size_t index = 0;
    size_t numnodes = 0;
    if( mesh ) {
        numnodes = mesh->getActiveSize(0);
        nodes.resize(numnodes,3);
        index = 0;
        size_t nsize = mesh->getSize(0);
        for( size_t i = 0; i < nsize; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isActive() ) {
                const Point3D &xyz = vtx->getXYZCoords();
                nodes(index,0) = xyz[0];
                nodes(index,1) = xyz[1];
                nodes(index,2) = xyz[2];
                index++;
            }
        }
    }

    return nodes;
}

/////////////////////////////////////////////////////////////////////////

JMeshEigenMatrix::EdgeMatrix JMeshEigenMatrix :: getEdgeMatrix()
{
    JMeshEigenMatrix::EdgeMatrix edges;

    size_t numedges = 0;
    size_t index = 0;

    if( mesh ) {
        numedges  = mesh->getActiveSize(1);
        edges.resize(numedges,2);

        size_t index = 0;
        size_t nsize = mesh->getSize(1);
        for( size_t i = 0; i < nsize; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                for( int j = 0; j < 2; j++)
                    edges(index,j) = edge->getNodeAt(j)->getID();
                index++;

            }
        }
    }

    return edges;
}

/////////////////////////////////////////////////////////////////////////

JMeshEigenMatrix::FaceMatrix JMeshEigenMatrix :: getFaceMatrix()
{
    JMeshEigenMatrix::FaceMatrix  faces;

    size_t numfaces;
    size_t index = 0;

    if( mesh ) {
        numfaces = mesh->getActiveSize(2);
        int elemType = mesh->getTopology()->getElementsType(2);
        if( elemType == JFace::TRIANGLE) 
            faces.resize(numfaces,3);

        else if( elemType == JFace::QUADRILATERAL) 
            faces.resize(numfaces,4);

        index = 0;
        size_t nsize = mesh->getSize(2);
        for( size_t i = 0; i < nsize; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                int nn = face->getSize(0);
                for( int j = 0; j < nn; j++)
                    faces(index,j) = face->getNodeAt(j)->getID();
                index++;
            }
        }
    }
    return faces;
}

/////////////////////////////////////////////////////////////////////////

JMeshEigenMatrix::CellMatrix JMeshEigenMatrix :: getCellMatrix()
{
    JMeshEigenMatrix::CellMatrix  cells;
    if( mesh == nullptr ) return cells;

    size_t numcells = mesh->getActiveSize(3);
    cells.resize(numcells,4);

    size_t index = 0;
    size_t nsize = mesh->getSize(3);
    for( size_t i = 0; i < nsize; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < cell->getSize(0); j++)
                cells(index,j) = cell->getNodeAt(j)->getID();
            index++;
        }
    }
    return cells;
}

/////////////////////////////////////////////////////////////////////////
JMeshEigenMatrix::ElemMatrix JMeshEigenMatrix :: getElementMatrix()
{
    JMeshEigenMatrix::ElemMatrix  elems;
    if( mesh == nullptr ) return elems;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim == 1) return getEdgeMatrix();
    if( topDim == 2) return getFaceMatrix();
    if( topDim == 3) return getCellMatrix();

    return elems;
}

/////////////////////////////////////////////////////////////////////////

JMeshPtr  JMeshEigenMatrix :: getMesh( const NodeMatrix &nodes, const FaceMatrix &faces)
{
    JMeshPtr  newmesh = JMesh::newObject();

    size_t numnodes = nodes.rows();
    size_t dim      = nodes.cols();

    JNodeSequence vnodes;
    if( numnodes ) {
        vnodes.resize(numnodes);
        Point3D xyz;
        xyz[0] = 0.0;
        xyz[1] = 0.0;
        xyz[2] = 0.0;
        for( size_t i = 0; i < numnodes; i++) {
            for( int j = 0; j < dim; j++)
                xyz[j] = nodes.coeff(i,j);
            vnodes[i] = JNode::newObject();
            vnodes[i]->setXYZCoords(xyz);
            vnodes[i]->setID(i);
        }
        newmesh->addObjects(vnodes);
    }

    size_t numfaces = faces.rows();
    if( numfaces == 0 ) return newmesh;

    int     nnodes  = faces.cols();

    JFacePtr dummyface;
    if( nnodes == 3) dummyface = JTriangle::newObject();
    if( nnodes == 4) dummyface = JQuadrilateral::newObject();
    JNodeSequence connect(nnodes);

    JFaceSequence newfaces(numfaces);
    for( size_t i = 0; i < numfaces; i++) {
        for( int j = 0; j < nnodes; j++)
            connect[j] =  vnodes[faces.coeff(i,j)];
        newfaces[i] = dummyface->getClone();
        newfaces[i]->setNodes( connect );
    }
    newmesh->addObjects(newfaces);
    return newmesh;
}
/////////////////////////////////////////////////////////////////////////

JGeneralSparseMatrix<int> JMeshAdjacencyMatrix::getMatrix() const
{
    JGeneralSparseMatrix<int> A;
    if( mesh == nullptr ) return A;

    size_t numnodes  = mesh->getActiveSize(0);
    size_t numedges  = mesh->getSize(1);

    A.setSize(numnodes,numnodes);
    size_t index = 0;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            int v0 = edge->getNodeAt(0)->getID();
            int v1 = edge->getNodeAt(1)->getID();
            A.setValue(v0,v1,1);
            A.setValue(v1,v0,1);
        }
    }
    return A;
}

/////////////////////////////////////////////////////////////////////////

JGeneralSparseMatrix<double> JMeshLaplaceMatrix::getMatrix() const
{
    JGeneralSparseMatrix<double> L;
    if( mesh == nullptr ) return L;

    JMeshAdjacencyMatrix meshadj;
    meshadj.setMesh(mesh);
    JGeneralSparseMatrix<int> b = meshadj.getMatrix();

    size_t numRows = b.numRows();
    L.setSize(numRows,numRows);

    JGeneralSparseMatrix<int>::SparseRow row;
    for( size_t i = 0; i < numRows; i++) {
        row = b.getRow(i);
        int sum = 0;
        for( auto &keyVal : row) {
            size_t j = keyVal.first;
            double v = keyVal.second;
            L.setValue(i,j,-v);
            sum += v;
        }
        L.setValue(i,i,sum);
    }
    return L;
}

/////////////////////////////////////////////////////////////////////////

Eigen::SparseMatrix<double> JEigenMatrixAdaptor::getMatrix( JGeneralSparseMatrix<double> &mat, bool clearAfterCopy)
{
    size_t numNonZeros = mat.getNumNonZeros();
    std::vector<Eigen::Triplet<double>>  triplets(numNonZeros);

    size_t numRows = mat.numRows();
    size_t numCols = mat.numCols();
    JGeneralSparseMatrix<double>::SparseRow row;

    size_t index = 0;
    for( size_t i = 0; i < numRows; i++) {
        row = mat.getRow(i);
        for( auto &keyVal : row) {
            size_t j   = keyVal.first;
            double val = keyVal.second;
            triplets[index++] = Eigen::Triplet<double>(i,j,val);
        }
        if( clearAfterCopy) mat.clearRow(i);
    }
    Eigen::SparseMatrix<double> egMat(numRows, numCols);
    egMat.setFromTriplets( triplets.begin(), triplets.end() );
    return egMat;
}

JGeneralSparseMatrix<double> JEigenMatrixAdaptor::getMatrix( Eigen::SparseMatrix<double> &mat, bool clearAfterCopy)
{
    JGeneralSparseMatrix<double>  gm;
    int numRows = mat.rows();
    int numCols = mat.cols();
    gm.setSize(numRows, numCols);
    for (int k=0; k<mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
            int i = it.row();   // row index
            int j = it.col();   // col index (here it is equal to k)
            double val = it.value();
            gm.setValue(i,j,val);
        }
    return gm;

}
