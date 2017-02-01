#include "AllHexMeshGenerator.hpp"

using namespace Jaal;

////////////////////////////////////////////////////////////////////////////////

JStructuredMesh3D :: JStructuredMesh3D()
{
    cellDim[0] = 1;
    cellDim[1] = 1;
    cellDim[2] = 1;

    nodeDim[0] = 2;
    nodeDim[1] = 2;
    nodeDim[2] = 2;
}
////////////////////////////////////////////////////////////////////////////////

size_t JStructuredMesh3D::getSize(int e) const
{
    if( mesh == nullptr) return 0;

    return mesh->getSize(e);
}

////////////////////////////////////////////////////////////////////////////////
void JStructuredMesh3D :: setCellDimensions( const vector<int> &dim)
{
    cellDim[0] = dim[0];
    cellDim[1] = dim[1];
    cellDim[2] = dim[2];

    nodeDim[0] = cellDim[0] + 1;
    nodeDim[1] = cellDim[1] + 1;
    nodeDim[2] = cellDim[2] + 1;

    numCells   = cellDim[0]*cellDim[1]*cellDim[2];
    numNodes   = nodeDim[0]*nodeDim[1]*nodeDim[2];
}

////////////////////////////////////////////////////////////////////////////////

JCellPtr JStructuredMesh3D :: getCellAt( int i, int j, int k) const
{
    if( mesh == nullptr) return nullptr;

    if( i < 0 || i >= cellDim[0] ) return nullptr;
    if( j < 0 || j >= cellDim[1] ) return nullptr;
    if( k < 0 || k >= cellDim[2] ) return nullptr;

    size_t pos = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
    return mesh->getCellAt(pos);
}

////////////////////////////////////////////////////////////////////////////////
JNodePtr JStructuredMesh3D :: getNodeAt( int i, int j, int k) const
{
    if( mesh == nullptr) return nullptr;

    if( i < 0 || i >= nodeDim[0] ) return nullptr;
    if( j < 0 || j >= nodeDim[1] ) return nullptr;
    if( k < 0 || k >= nodeDim[2] ) return nullptr;

    size_t pos =  k*nodeDim[0]*nodeDim[1] + j*nodeDim[0] + i;
    return mesh->getNodeAt(pos);
}
////////////////////////////////////////////////////////////////////////////////

int JStructuredMesh3D :: getIJK( const JCellPtr &cell, int &i, int &j, int &k) const
{
    size_t cellID;
    int err = cell->getAttribute("VoxelID", cellID);
    if( err ) {
        cout << "Warning: Cell Voxel ID not specified " << endl;
        return 1;
    }

    if( cellID < 0 ||  cellID >= numCells ) return 1;
    k = cellID/(cellDim[0]*cellDim[1]);
    j = (cellID - k*cellDim[0]*cellDim[1])/cellDim[0];
    i = cellID - k*cellDim[0]*cellDim[1] - j*cellDim[0];

//  assert( cellID == k*cellDim[0]*cellDim[1] + j*cellDim[0] + i);
    assert( i >= 0 || i < cellDim[0] );
    assert( j >= 0 || j < cellDim[1] );
    assert( k >= 0 || k < cellDim[2] );

    return 0;
}
////////////////////////////////////////////////////////////////////////////////
int JStructuredMesh3D :: getIJK( const JNodePtr &vtx, int &i, int &j, int &k) const
{
    size_t nodeID;
    int err = vtx->getAttribute("VoxelID", nodeID);
    if( err ) {
        cout << "Warning: Node Voxel ID not specified " << endl;
        return 1;
    }
    if( nodeID < 0 ||  nodeID >= numNodes ) return 1;

    k = nodeID%(nodeDim[0]*nodeDim[1]);
    j = (nodeID - k*nodeDim[0]*nodeDim[1])%nodeDim[0];
    i = nodeID - k*nodeDim[0]*nodeDim[1] -j*nodeDim[0];
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JStructuredMesh3D :: getRelations( const JNodePtr &vtx, JNodeSequence &neighs) const
{
    neighs.clear();

    int i, j, k;
    int err = getIJK(vtx, i, j, k);
    if( err ) return err;

    JNodePtr node;
    node = getNodeAt( i+1,j,k);
    if( node ) neighs.push_back(node);

    node = getNodeAt( i-1,j,k );
    if( node ) neighs.push_back(node);

    node = getNodeAt( i,j+1,k);
    if( node ) neighs.push_back(node);

    node = getNodeAt( i,j-1,k);
    if( node ) neighs.push_back(node);

    node = getNodeAt( i,j, k-1);
    if( node) neighs.push_back(node);

    node = getNodeAt( i,j, k+1 );
    if( node ) neighs.push_back(node);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
int JStructuredMesh3D :: getRelations( const JNodePtr &vtx, JCellSequence &neighs) const
{
    neighs.clear();

    int i, j, k;
    int err = getIJK(vtx, i, j, k);
    if( err ) return err;

    JCellPtr cell;
    cell = getCellAt( i-1,j-1,k-1);
    if( cell ) neighs.push_back(cell);

    cell = getCellAt( i,j-1,k-1 );
    if( cell ) neighs.push_back(cell);

    cell = getCellAt( i-1,j,k-1);
    if( cell) neighs.push_back(cell);

    cell = getCellAt( i,j,k-1);
    if( cell ) neighs.push_back(cell);

    cell = getCellAt( i-1,j-1,k);
    if( cell ) neighs.push_back(cell);

    cell = getCellAt( i,j-1,k );
    if( cell ) neighs.push_back(cell);

    cell = getCellAt( i-1,j,k);
    if( cell) neighs.push_back(cell);

    cell = getCellAt( i,j,k);
    if( cell) neighs.push_back(cell);

    cout << "NEIGH " << neighs.size() << endl;

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JStructuredMesh3D :: getRelations6( const JCellPtr &cell, JCellSequence &neighs) const
{
    neighs.clear();

    int i, j, k;
    int err = getIJK(cell, i, j, k);
    if( err ) return 1;

    JCellPtr nextCell;
    nextCell = getCellAt( i+1,j,k);
    if( nextCell ) neighs.push_back(nextCell);

    nextCell = getCellAt( i-1,j,k );
    if( nextCell ) neighs.push_back(nextCell);

    nextCell = getCellAt( i,j+1,k );
    if( nextCell ) neighs.push_back(nextCell);

    nextCell = getCellAt( i,j-1,k );
    if( nextCell ) neighs.push_back(nextCell);

    nextCell = getCellAt( i,j, k-1 );
    if( nextCell ) neighs.push_back(nextCell);

    nextCell = getCellAt( i,j, k+1);
    if( nextCell ) neighs.push_back(nextCell);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllHexMeshGenerator::getStructuredMesh(int *node_dim, double *length, double *origin)
{
    JMeshPtr mesh = JMesh::newObject();

    int nx = 1, ny = 1, nz = 1;

    if( node_dim ) {
        nx = node_dim[0];
        ny = node_dim[1];
        nz = node_dim[2];
    }

    double xorg = 0.0, yorg = 0.0, zorg = 0.0;
    if( origin  ) {
        xorg = origin[0];
        yorg = origin[1];
        zorg = origin[2];
    }

    double xlen = 1.0, ylen = 1.0, zlen = 1.0;
    if( length ) {
        xlen = length[0];
        ylen = length[1];
        zlen = length[2];
    }

    double  dx = 0.0, dy = 0.0, dz = 0.0;
    if( nx > 1) dx = xlen/ (double) (nx - 1);
    if( ny > 1) dy = ylen/ (double) (ny - 1);
    if (nz > 1) dz = zlen/ (double) (nz - 1);

    JNodeSequence nodes = JNode::newObjects(nx * ny * nz);

    int offset;
    Point3D xyz;
    int index = 0;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                offset = k * nx * ny + j * nx + i;
                xyz[0] = xorg + i*dx;
                xyz[1] = yorg + j*dy;
                xyz[2] = zorg + k*dz;
                nodes[offset]->setXYZCoords(xyz);
                nodes[offset]->setID(index++);
            }
        }
    }

    mesh->addObjects(nodes);

    JCellSequence cells = JHexahedron::newObjects((nx-1) * (ny-1) * (nz-1));

    JNodeSequence qnodes(8);
    index = 0;
    for (int k = 0; k < nz - 1; k++) {
        for (int j = 0; j < ny - 1; j++) {
            for (int i = 0; i < nx - 1; i++) {
                offset = k*nx*ny + j * nx + i;
                qnodes[0] = nodes[offset];
                qnodes[1] = nodes[offset + 1];
                qnodes[2] = nodes[offset + 1 + nx];
                qnodes[3] = nodes[offset + nx];
                offset = (k+1)*nx*ny + j * nx + i;
                qnodes[4] = nodes[offset];
                qnodes[5] = nodes[offset + 1];
                qnodes[6] = nodes[offset + 1 + nx];
                qnodes[7] = nodes[offset + nx];
                cells[index]->setID(index);
                cells[index]->setNodes(qnodes);
                index++;
            }
        }
    }
    mesh->addObjects(cells);

    return mesh;
}

////////////////////////////////////////////////////////////////////////////////
int AllHexMeshGenerator :: stackQuadMesh( const JMeshPtr &mesh1, const JMeshPtr &mesh2, JCellSequence &newcells)
{
    newcells.clear();

    size_t numFaces = mesh1->getSize(2);
    for( size_t j = 0; j < numFaces; j++) {
        JFacePtr f0 =  mesh1->getFaceAt(j);
        JFacePtr f1 =  mesh2->getFaceAt(j);
        JCellPtr hex = JHexahedron::newObject(f0, f1);
        newcells.push_back(hex);
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////

JMeshPtr AllHexMeshGenerator :: stackQuadMesh( vector<JMeshPtr> &meshlayers, bool closed )
{
    JMeshPtr hexmesh = JMesh::newObject();

    int nstacks = meshlayers.size();

    JNodeSequence nodes;
    for( int i = 0; i < nstacks; i++) {
        nodes = meshlayers[i]->getNodes();
        hexmesh->addObjects( nodes );
    }

    JCellSequence newcells;
    for( int i = 0; i < nstacks-1; i++) {
        stackQuadMesh( meshlayers[i], meshlayers[i+1], newcells);
        hexmesh->addObjects( newcells );
    }

    if( closed ) {
        stackQuadMesh( meshlayers[nstacks-1], meshlayers[0], newcells);
        hexmesh->addObjects( newcells );
    }

    return hexmesh;
}
////////////////////////////////////////////////////////////////////////////////

bool AllHexMeshGenerator :: isDoublet(const JNodePtr &v)
{
    if( v->isBoundary() ) return 0;

    JCellSequence cellneighs;
    JNode::getRelations(v, cellneighs);

    if( cellneighs.size() == 2 ) return 1;

    if( cellneighs.empty() )
        cout << "Warning: Cell Neighbours empty: Doublet unknown" << endl;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int AllHexMeshGenerator :: removeDoublet( const JMeshPtr &mesh, const JNodePtr &vdoublet)
{
    assert( mesh->getAdjTable( 0,3) );

    JCellSequence cellneighs;
    JNode::getRelations(vdoublet, cellneighs);

    if( cellneighs.size() != 2 ) return 1;

    JHexahedronPtr hex1 = JHexahedron::down_cast( cellneighs[0] );
    JHexahedronPtr hex2 = JHexahedron::down_cast( cellneighs[1] );

    JFaceSequence faceneighs;
    JHexahedron::get_shared_entities( hex1, hex2, faceneighs);
    boost::sort( faceneighs );

    if( faceneighs.size() != 3 ) return 1;

    JFaceSequence hexfaces, newfaces;

    hexfaces = hex1->getFaces();
    assert( hexfaces.size() == 6 );
    boost::sort( hexfaces );
    boost::set_difference( hexfaces, faceneighs, back_inserter( newfaces ) );

    hexfaces = hex2->getFaces();
    assert( hexfaces.size() == 6 );
    boost::sort( hexfaces );
    boost::set_difference( hexfaces, faceneighs, back_inserter( newfaces ));

    assert( newfaces.size() == 6 );

    if( newfaces.size() != 6 ) return 2;

    JCellPtr mergedHex = JHexahedron::newObject( newfaces[0], newfaces[1], newfaces[2],
                         newfaces[3], newfaces[4], newfaces[5] );
    if( mergedHex == NULL ) return 3;

    mesh->remove( cellneighs[0] );
    mesh->remove( cellneighs[1] );
    mesh->remove( vdoublet );
    mesh->addObject( mergedHex );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
void AllHexMeshGenerator :: searchDoublets( const JMeshPtr &mesh, JNodeSequence &doublets)
{
    doublets.clear();
    size_t nSize = mesh->getSize(0);

    if( mesh->getAdjTable(2,3) == 0)  mesh->buildRelations(2,3);

    for( size_t i = 0; i < nSize; i++) {
        JNodePtr v = mesh->getNodeAt(i);
        if( isDoublet(v) )  doublets.push_back(v);
    }
}
////////////////////////////////////////////////////////////////////////////////
JMeshPtr AllHexMeshGenerator::Tetrahedral2Hexahedral(const JMeshPtr &tetmesh)
{
    if( tetmesh->getTopology()->getElementsType(3) != JCell::TETRAHEDRON) {
        tetmesh->getLogger()->setWarn("Tet->Hex require all-tet mesh as an input ");
        return NULL;
    }

    JMeshPtr hexmesh = JMesh::newObject();
    tetmesh->getTopology()->searchBoundary();

    size_t numnodes = tetmesh->getSize(0);
    size_t numedges = tetmesh->getSize(1);
    size_t numfaces = tetmesh->getSize(2);
    size_t numcells = tetmesh->getSize(3);

    for( size_t i = 0; i < numnodes; i++)
        hexmesh->addObject( tetmesh->getNodeAt(i));

    int val = 1;
    Point3D p3d;
    JNodePtr steiner, vh;
    for( size_t i = 0; i < numedges; i++) {
        JEdgePtr edge = tetmesh->getEdgeAt(i);
        if( edge->isActive() )  {
            steiner = JNode::newObject();
            val = 1;
            steiner->setAttribute("Steiner", val);
            if( edge->isBoundary() ) {
                steiner->setAttribute("GeoDimension", 1);
            }
            edge->getAvgXYZ(p3d);
            steiner->setXYZCoords(p3d);
            edge->setAttribute("Steiner", steiner);
            hexmesh->addObject(  steiner );
            if( edge->isBoundary() ) {
//                    steiner->setBoundaryMark( max(1,edge->getBoundaryMark() ));
                cout << "Incomplete " << endl;
                exit(0);
            }
            if( edge->hasAttribute("Rigid") ) {
                JEdgePtr edge1 = JEdge::newObject( edge->getNodeAt(0), steiner);
                JEdgePtr edge2 = JEdge::newObject( edge->getNodeAt(1), steiner);
                edge1->setAttribute("Rigid", 1);
                edge2->setAttribute("Rigid", 1);
                vh = edge1->getHashNode();
                vh->attach(edge1);
                vh = edge2->getHashNode();
                vh->attach(edge2);
            }
        }
    }

    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = tetmesh->getFaceAt(i);
        if( face->isActive() ) {
            assert( face->getSize(0) == 3 ) ;  // Triangle face ...
            steiner = JNode::newObject();
            val = 2;
            steiner->setAttribute("Steiner", val);
            if( face->isBoundary() )  {
                steiner->setAttribute("GeoDimension", 2);
            }
            face->getAvgXYZ(p3d);
            steiner->setXYZCoords(p3d);
            face->setAttribute("Steiner", steiner);
            hexmesh->addObject(  steiner );
            if( face->isBoundary() ) {
//                  steiner->setBoundaryMark( max(1,face->getBoundaryMark() ));
                cout << "Incomplete " << endl;
                exit(0);
            }
        }
    }

    for( size_t i = 0; i < numcells; i++) {
        JCellPtr cell = tetmesh->getCellAt(i);
        if( cell->isActive() ) {
            steiner = JNode::newObject();
            val = 3;
            steiner->setAttribute("Steiner", val);
            cell->getAvgXYZ(p3d);
            steiner->setXYZCoords(p3d);
            cell->setAttribute("Steiner", steiner);
            hexmesh->addObject(  steiner );
        }
    }

    vector<JHexahedron> hex;
    for( size_t i = 0; i < numcells; i++) {
        JCellPtr cell = tetmesh->getCellAt(i);
        if( cell->isActive() ) {
            JTetrahedronPtr tet = JTetrahedron::down_cast(cell);
            tet->getHexCells(hex);
            for( size_t j = 0; j < hex.size(); j++) {
                JCellPtr h = hex[j].getClone();
                hexmesh->addObject(h);
            }
        }
    }

    tetmesh->deleteEdgeAttribute("Steiner");
    tetmesh->deleteFaceAttribute("Steiner");
    tetmesh->deleteCellAttribute("Steiner");

    return hexmesh;
}

////////////////////////////////////////////////////////////////////////////////

int AllHexMeshGenerator::setSurfaceQuadMesh(const JMeshPtr &m)
{
    int homog = m->getTopology()->getElementsType(2);

    if( homog != JCell::TETRAHEDRON) {
        cout << "Warning: Hexmeshing requires All-Quads surface mesh " << endl;
        return 1;
    }

    if( !m->getTopology()->isClosed() ) {
        cout << "Warning: The input mesh for hexahedral must be closed " << endl;
        return 2;
    }

#ifdef CSV
    if( !m->getTopology()->isOrientable() ) {
        cout << "Warning: The input mesh for hexahedral must be orientable " << endl;
        return 2;
    }

    if( !m->getTopology()->isSimple() ) {
        cout << "Warning: The input mesh for hexahedral must be simple" << endl;
        return 2;
    }

    surfmesh = m;

    // Step 1: Convert the Quadrilateral meshing into Triangle mesh;
    Mesh2D *trimesh;


    // Step 2: Use triangle mesh to generate tetrahedral mesh ....

    // Step 3: Use Tetrahedral mesh to get All-Hex mesh;

    // Step 4; Improve the tiopological hex mesh.
#endif
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int AllHexMeshGenerator :: rotate_boundary_edge( const JMeshPtr &mesh, const JEdgePtr &boundedge,
        JEdgeSequence &affectedEdges)
{

#ifdef LATER
    affectedEdges.clear();
    if( boundedge == NULL || !boundedge->isActive() ) return 1;
    affectedEdges.push_back( boundedge );

    if( mesh->getAdjTable(1,2) == 0) mesh->buildRelations(1,2);
    if( mesh->getAdjTable(1,3) == 0) mesh->buildRelations(1,3);
    if( mesh->getAdjTable(2,3) == 0) mesh->buildRelations(2,3);

    JNodePtr v0 = boundedge->getNodeAt(0);
    JNodePtr v1 = boundedge->getNodeAt(1);

    if( !v0->isBoundary() || !v1->isBoundary() ) {
        cout << "Info: Currently only boundary edges are rotable " << endl;
        return 1;
    }

    JCellSequence cellneighs;
    JFaceSequence faceneighs;

    JEdge::getRelations( boundedge, faceneighs);
    if( faceneighs.size() != 3 ) {
        cout << "Warning: Boundary edge should have three faces: One internal and two boundary" << endl;
        return 1;
    }

    JFacePtr partFace;
    JFaceSequence currSideFaces, nextSideFaces(2);
    JCellSequence currSideCells(2), prevSideCells(2);

    for( int i = 0; i < 3; i++) {
        if(  faceneighs[i]->isBoundary() ) currSideFaces.push_back( faceneighs[i] );
        if( !faceneighs[i]->isBoundary() ) partFace = faceneighs[i];
    }

    if( partFace == nullptr) return 1;

    if( currSideFaces.size() != 2 ) {
        cout << "Warning: Boundary edge doesn't have two boundary faces as neighs" << endl;
        return 1;
    }

    JHexahedronPtr hex;
    JEdgePtr currEdge = boundedge;
    prevSideCells[0] = nullptr;
    prevSideCells[1] = nullptr;

    while( 1  ) {
        JFace::getRelations(currSideFace[0], cellneighs);
        if( cellneighs.size() == 1) {
            currSideCells[0] = cellneighs[0];
        } else {
            if( cellneighs[0] == prevSideCells[0] ) currSideCells[0] = cellneighs[1];
            if( cellneighs[1] == prevSideCells[0] ) currSideCells[0] = cellneighs[0];
        }
        prevSideCells[0] = currSideCells[0];

        JFace::getRelations(currSideFace[1], cellneighs);
        if( cellneighs.size() == 1) {
            currSideCells[1] = cellneighs[0];
        } else {
            if( cellneighs[0] == prevSideCells[1] ) currSideCells[1] = cellneighs[1];
            if( cellneighs[1] == prevSideCells[1] ) currSideCells[1] = cellneighs[0];
        }
        prevSideCells[1] = currSideCells[1];

        hex = Hexahedron::down_cast( currSideCells[0] );
        nextSideFaces[0] = hex->getOppositeFace( currSideFaces[0] );

        hex = Hexahedron::down_cast( currSideCells[1] ) ;
        nextSideFaces[1] = hex->getOppositeFace( currSideFaces[1] );

        JQuadrilateralPtr quad = Quadrilateral::down_cast( partFace );
        currEdge = quad->getOppositeEdge(currEdge);
        affectedEdges.push_back( currEdge );

        if( currEdge->isBoundary() ) break;

        JEdge::getRelations(currEdge, faceneighs);
        if( faceneighs.size() != 4 ) {
            cout << "Info: Internal edge must have four adjacent faces " << endl;
            return 1;
        }

        JFacePtr nextFace = nullptr;
        for( int i = 0; i < 4; i++) {
            if( faceneighs[i] == partFace ||
                    faceneighs[i] == nextSideFaces[0] ||
                    faceneighs[i] == nextSideFaces[1] ) continue;
            nextFace = faceneighs[i];
            break;
        }

        if( nextFace == nullptr) {
            affectedEdges.clear();
            affectedEdges.push_back( boundedge );
            return 1;
        }
        partFace = nextFace;
        currSideFaces[0] = nextSideFaces[0];
        currSideFaces[1] = nextSideFaces[1];
    }
#endif
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
int AllHexMeshGenerator :: rotate3hex(const JEdgePtr &edge, const JMeshPtr &)
{
    // Condition : The given edge must share exactly three hexahedral cell
    JCellSequence cellneighs;
    JEdge::getRelations(edge, cellneighs);
    if( cellneighs.size() != 3 ) return 1;

    JExit();
    return 1;
}
///////////////////////////////////////////////////////////////////////////////
