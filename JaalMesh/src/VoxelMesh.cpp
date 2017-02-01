#include "AllHexMeshGenerator.hpp"
#include <boost/lexical_cast.hpp>

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JVoxelMesh ::getBackgroundMesh(int *gdim, double *glength, double *gorigin)
{
    cellDim[0] = gdim[0];
    cellDim[1] = gdim[1];
    cellDim[2] = gdim[2];
    numCells   = cellDim[0]*cellDim[1]*cellDim[2];

    nodeDim[0] = gdim[0] + 1;
    nodeDim[1] = gdim[1] + 1;
    nodeDim[2] = gdim[2] + 1;
    numNodes   = nodeDim[0]*nodeDim[1]*nodeDim[2];

    if( glength ) {
        length[0] = glength[0];
        length[1] = glength[1];
        length[2] = glength[2];
    }

    if( origin) {
        origin[0] = gorigin[0];
        origin[1] = gorigin[1];
        origin[2] = gorigin[2];
    }

    mesh  = AllHexMeshGenerator::getStructuredMesh(nodeDim, length, origin);

    if( mesh == nullptr) {
        cout << "Warning: Background mesh not created " << endl;
        return nullptr;
    }

    if( mesh->getSize(0) != numNodes) {
        cout << "Warning: Background mesh has different nodes " << endl;
        return nullptr;
    }

    if( mesh->getSize(3) != numCells) {
        cout << "Warning: Background mesh has different nodes " << endl;
        return nullptr;
    }

    // All cells are inactive in the beginning ...
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell =  mesh->getCellAt(i);
        cell->setAttribute("VoxelID", i);
        cell->setStatus(JMeshEntity::INACTIVE);
    }

    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx =  mesh->getNodeAt(i);
        vtx->setAttribute("VoxelID", i);
    }



    return mesh;
}

////////////////////////////////////////////////////////////////////////////////

boost::tuple<int,int> JVoxelMesh :: getRange( int dir) const
{
    if( dir == 0) return getXRange();
    if( dir == 1) return getYRange();
    if( dir == 2) return getZRange();
    boost::tuple<int,int>   empty;
    return empty;
}

/////////////////////////////////////////////////////////////////////////////////////////////

boost::tuple<int,int> JVoxelMesh :: getXRange() const
{
    int minVal = 0, maxVal = 0;
    for( int i = 0; i < cellDim[0]; i++) {
        for( int j = 0; j < cellDim[1]; j++) {
            for( int k = 0; k < cellDim[2]; k++) {
                size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
                const JCellPtr &cell = mesh->getCellAt(vid);
                if( cell->isActive() ) {
                    minVal = i;
                    break;
                }
            }
        }
    }

    for( int i = cellDim[0]-1; i >=0 ; i--) {
        for( int j = 0; j < cellDim[1]; j++) {
            for( int k = 0; k < cellDim[2]; k++) {
                size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
                const JCellPtr &cell = mesh->getCellAt(vid);
                if( cell->isActive() ) {
                    maxVal = i;
                    break;
                }
            }
        }
    }
    return boost::make_tuple(minVal,  maxVal);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

boost::tuple<int,int> JVoxelMesh :: getYRange() const
{

    int minVal = 0, maxVal = 0;
    for( int j = 0; j < cellDim[1]; j++) {
        for( int i = 0; i < cellDim[0]; i++) {
            for( int k = 0; k < cellDim[2]; k++) {
                size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
                const JCellPtr &cell = mesh->getCellAt(vid);
                if( cell->isActive() ) {
                    minVal = i;
                    break;
                }
            }
        }
    }

    for( int j = 0; j < cellDim[1]; j++) {
        for( int i = cellDim[0]-1; i >=0 ; i--) {
            for( int k = 0; k < cellDim[2]; k++) {
                size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
                const JCellPtr &cell = mesh->getCellAt(vid);
                if( cell->isActive() ) {
                    maxVal = i;
                    break;
                }
            }
        }
    }
    return boost::make_tuple(minVal,  maxVal);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

boost::tuple<int,int> JVoxelMesh :: getZRange() const
{
    int minVal = 0, maxVal = 0;
    for( int k = 0; k < cellDim[2]; k++) {
        for( int j = 0; j < cellDim[1]; j++) {
            for( int i = 0; i < cellDim[0]; i++) {
                size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
                const JCellPtr &cell = mesh->getCellAt(vid);
                if( cell->isActive() ) {
                    minVal = i;
                    break;
                }
            }
        }
    }

    for( int k = 0; k < cellDim[2]; k++) {
        for( int j = 0; j < cellDim[1]; j++) {
            for( int i = cellDim[0]-1; i >=0 ; i--) {
                size_t vid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
                const JCellPtr &cell = mesh->getCellAt(vid);
                if( cell->isActive() ) {
                    maxVal = i;
                    break;
                }
            }
        }
    }
    return boost::make_tuple(minVal,  maxVal);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
int JVoxelMesh ::  getNumSingularPoints() const
{
    size_t numNodes = mesh->getSize(0);
    mesh->buildRelations(0,3);

    size_t nCount = 0;
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        int numCells = vtx->getNumRelations(3);
        if( numCells == 3 || numCells == 7 ) nCount++;
    }
    return nCount;
}
/////////////////////////////////////////////////////////////////////////////////////////////

string JVoxelMesh ::  getBitString() const
{
    string str;
    if( mesh == nullptr) return str;

    ostringstream oss;

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        oss << cell->isActive();
    }
    return oss.str();
}

////////////////////////////////////////////////////////////////////////////////////////////

JMeshPtr JVoxelMesh :: getModelMesh()
{
    if( mesh == nullptr) return nullptr;

    if( modelmesh == nullptr) modelmesh  = JMesh::newObject();

    modelmesh->clearAll();
    assert( modelmesh->getSize(0) == 0);
    assert( modelmesh->getSize(1) == 0);
    assert( modelmesh->getSize(2) == 0);
    assert( modelmesh->getSize(3) == 0);

    size_t numCells = mesh->getSize(3);
    JNodeSet nodeSet;
    JEdgeSet edgeSet;
    JFaceSet faceSet;

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < 8; j++)
                nodeSet.insert( cell->getNodeAt(j) );
            for( int j = 0; j < 12; j++)
                edgeSet.insert( cell->getEdgeAt(j) );
            for( int j = 0; j < 6; j++)
                faceSet.insert( cell->getFaceAt(j) );
        }
    }

    for( const JNodePtr &vtx : nodeSet)
        modelmesh->addObject(vtx);

    for( const JEdgePtr &edge : edgeSet)
        modelmesh->addObject(edge);

    for( const JFacePtr &face : faceSet)
        modelmesh->addObject(face);

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) modelmesh->addObject(cell);
    }
    modelmesh->enumerate(0);
    modelmesh->enumerate(1);
    modelmesh->enumerate(2);
    modelmesh->enumerate(3);

    return modelmesh;
}

/////////////////////////////////////////////////////////////////////////////////
int JVoxelMesh :: saveAs(const string &fname) const
{
    ofstream ofile( fname.c_str(), ios::out);
    if( ofile.fail() ) return 1;

    ofile << "<CellDimensions>" << "  "
          << cellDim[0] << " " << cellDim[1]  << " " << cellDim[2] << " </CellDimensions>" << endl;
    ofile << "<Boxlength>" << "  "
          << length[0] << " " << length[1]  << " " << length[2] << " </Boxlength>" << endl;
    ofile << "<Origin>" << "  "
          << origin[0] << " " << origin[1]  << " " << origin[2] << " </Origin>" << endl;
    ofile << "<StorageOrder>" << "  ZYX " << "</StorageOrder>" << endl;
    ofile << "<BitString>" << endl;
    string str = getBitString();
    ofile << str << endl;
    ofile << "</BitString>" << endl;

    ofile << endl;
    return 0;
}
//////////////////////////////////////////////////////////////////////////////////

int JVoxelMesh :: readFrom( const string &fname)
{
    ifstream infile( fname.c_str(), ios::in);
    if( infile.fail() ) return 1;

    string str, voxstr;
    infile >> str;
    if( str != "<CellDimensions>" ) {
        cout <<"Warning: Cell Dimenstion  tag missing" << endl;
        return 1;
    }
    infile >> cellDim[0] >> cellDim[1] >> cellDim[2];
    infile >> str;
    if( str != "</CellDimensions>" ) {
        cout <<"Warning: Cell Dimenstion  tag missing" << endl;
        return 1;
    }
    infile >> str;
    if( str != "<Boxlength>") {
        cout <<"Warning: Box length is missing" << endl;
        return 1;
    }
    infile >> length[0] >> length[1] >> length[2];
    infile >> str;
    if( str != "</Boxlength>") {
        cout <<"Warning: Box length is missing" << endl;
        return 1;
    }

    infile >> str;
    if( str != "<Origin>") {
        cout <<"Warning: Origin is missing" << endl;
        return 1;
    }
    infile >> origin[0] >> origin[1] >> origin[2];
    infile >> str;
    if( str != "</Origin>") {
        cout <<"Warning: Box length is missing" << endl;
        return 1;
    }

    infile >> str;
    if( str != "<StorageOrder>") {
        cout <<"Warning: Storage Order is missing" << endl;
        return 1;
    }
    infile >> str;
    infile >> str;
    if( str != "</StorageOrder>") {
        cout <<"Warning: Storage Order is missing" << endl;
        return 1;
    }
    infile >> str;
    if( str != "<BitString>") {
        cout <<"Warning: BitString is  missing" << endl;
        return 1;
    }
    infile >> voxstr;

    numCells   = cellDim[0]*cellDim[1]*cellDim[2];
    if( voxstr.size() != numCells) {
        cout <<"Warning: Bitstring size not matching " << endl;
        return 1;
    }

    infile >> str;
    if( str != "</BitString>") {
        cout <<"Warning: Box length is missing" << endl;
        return 1;
    }

    nodeDim[0] = cellDim[0] + 1;
    nodeDim[1] = cellDim[1] + 1;
    nodeDim[2] = cellDim[2] + 1;
    numNodes   = nodeDim[0]*nodeDim[1]*nodeDim[2];

    mesh  = AllHexMeshGenerator::getStructuredMesh(nodeDim, length, origin);

    if( mesh == nullptr) {
        cout << "Warning: Background mesh not created " << endl;
        return 1;
    }


    if( mesh->getSize(0) != numNodes) {
        cout << "Warning: Background mesh has different nodes " << endl;
        return 1;
    }

    if( mesh->getSize(3) != numCells) {
        cout << "Warning: Background mesh has different nodes " << endl;
        return 1;
    }

    // All cells are inactive in the beginning ...
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell =  mesh->getCellAt(i);
        cell->setAttribute("VoxelID", i);
        int val = boost::lexical_cast<int> (voxstr[i]);
        if( val )
            cell->setStatus(JMeshEntity::ACTIVE);
        else
            cell->setStatus(JMeshEntity::INACTIVE);
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
bool JVoxelMesh :: isBoundary( int i, int j, int k) const
{
    if( i == 0 || i == cellDim[0] ) return 1;
    if( j == 0 || j == cellDim[1] ) return 1;
    if( k == 0 || k == cellDim[2] ) return 1;
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

int JVoxelMesh :: makeSolid()
{
    size_t numVoxels = cellDim[0]*cellDim[1]*cellDim[2];
    if( numVoxels == 0) return 0;

    for( size_t i = 0; i < numVoxels; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        cell->setVisitBit(0);
    }

    JCellSequence neighs, boundCells;

    for( int k = 0; k < cellDim[2]; k++) {
        for( int j = 0; j < cellDim[1]; j++) {
            for( int i = 0; i < cellDim[0]; i++) {
                if( isBoundary(i,j,k) ) {
                    size_t cid = k*cellDim[0]*cellDim[1] + j*cellDim[0] + i;
                    const JCellPtr &cell = mesh->getCellAt(cid);
                    boundCells.push_back(cell);
                }
            }
        }
    }

    std::deque<JCellPtr>  cellQ;
    for( size_t i = 0; i < boundCells.size(); i++) {
        if( boundCells[i]->getVisitBit() == 0) {
            cellQ.clear();
            cellQ.push_back( boundCells[i] );
            while(!cellQ.empty() ) {
                const JCellPtr &currCell = cellQ.front();
                cellQ.pop_front();
                if( currCell->getVisitBit() == 0) {
                    currCell->setVisitBit(1);
                    getRelations6(currCell, neighs);
                    for( const JCellPtr &nextCell : neighs) {
                        if( (nextCell->getVisitBit() == 0) && (!nextCell->isActive()))
                            cellQ.push_back(nextCell);
                    }
                }
            }
        }
    }

    for( size_t i = 0; i < numVoxels; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->getVisitBit() == 0)
            cell->setStatus(JMeshEntity::ACTIVE);
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
int JMonotoneVoxelizer :: getNumSingularPoints( const JCellPtr &cell) const
{
    int nCount = 0;
    for( int i = 0; i < 8; i++) {
        const JNodePtr &vtx = cell->getNodeAt(i);
        int numneighs = vtx->getNumRelations(3);
        if( numneighs == 1 || numneighs == 7) nCount++;
    }

    return nCount;
}

/////////////////////////////////////////////////////////////////////////////////////////////

int JMonotoneVoxelizer :: getChangeInSingularPointsAfterAddition( const JCellPtr &cell)
{
    if( cell->isActive() ) return 0;

    int n0 = getNumSingularPoints(cell);
    cell->setStatus( JMeshEntity::ACTIVE);
    for( int j = 0; j < 8; j++)
        cell->getNodeAt(j)->addRelation(cell);
    int n1 = getNumSingularPoints(cell);
    cell->setStatus( JMeshEntity::INACTIVE);
    return n1-n0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
int  JMonotoneVoxelizer :: getChangeInSingularPointsAfterDeletion( const JCellPtr &cell)
{
    if( !cell->isActive() ) return 0;

    int n0 = getNumSingularPoints(cell);
    cell->setStatus( JMeshEntity::INACTIVE);
    int n1 = getNumSingularPoints(cell);
    cell->setStatus( JMeshEntity::ACTIVE);
    for( int j = 0; j < 8; j++)
        cell->getNodeAt(j)->addRelation(cell);

    return n1-n0;
}
/////////////////////////////////////////////////////////////////////////////////////////////

void JMonotoneVoxelizer :: makeGroups()
{
    addGroup.clear();
    delGroup.clear();

    JMeshPtr bgmesh = voxmesh->getBackgroundMesh();
    vector<int> cellDim = voxmesh->getCellDimensions();
    size_t numCells = cellDim[0]*cellDim[1]*cellDim[2];

    JCellSequence neighs;
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = bgmesh->getCellAt(i);
        voxmesh->getRelations6(cell, neighs);
        int numactive = 0;
        for( const JCellPtr &neigh : neighs)
            if( neigh->isActive() ) numactive++;

        if(  cell->isActive() && numactive < 6) delGroup.push_back(cell);
        if( !cell->isActive() && numactive > 0) addGroup.push_back(cell);
        if( cell->isActive() ) {
            for( int j = 0; j < 8; j++)
                cell->getNodeAt(j)->addRelation(cell);
        }
    }

    for( const JCellPtr &cell: addGroup ) {
        int x = getChangeInSingularPointsAfterAddition(cell);
        cell->setAttribute("SingularIndex", x);
    }
    sort( addGroup.begin(), addGroup.end(), []( const JCellPtr &a, const JCellPtr &b)
    {
        int x, y;
        a->getAttribute("SingularIndex", x);
        b->getAttribute("SingularIndex", y);
        return x < y;
    });

    for( const JCellPtr &cell: delGroup ) {
        int x = getChangeInSingularPointsAfterDeletion(cell);
        cell->setAttribute("SingularIndex", x);
    }

    sort( delGroup.begin(), delGroup.end(), []( const JCellPtr &a, const JCellPtr &b)
    {
        int x, y;
        a->getAttribute("SingularIndex", x);
        b->getAttribute("SingularIndex", y);
        return x < y;
    });

    bgmesh->deleteCellAttribute("SingularIndex");
}

/////////////////////////////////////////////////////////////////////////////////////////////

void JMonotoneVoxelizer :: optimize()
{
    makeGroups();
    JCellPtr addcell, delcell;
    size_t nSize1, nSize2;
    nSize1 = addGroup.size();

    // First round: Which empty voxels can be filled ...

    for( size_t i = 0; i < nSize1; i++) {
        const JCellPtr &addcell = addGroup[i];
        int x = getChangeInSingularPointsAfterAddition(addcell);
        nSize2 = delGroup.size();
        for( size_t j = 0; j < nSize2; j++) {
            const JCellPtr &delcell = delGroup.front();
            delGroup.pop_front();
            int y = getChangeInSingularPointsAfterDeletion(delcell);
            if( x + y < 0) {
                addcell->setStatus(JMeshEntity::ACTIVE);
                delcell->setStatus(JMeshEntity::INACTIVE);
                for( int k = 0; k < 8; k++)
                    addcell->getNodeAt(k)->addRelation(addcell);
                break;
            } else
                delGroup.push_back(delcell);
        }
    }

    // Second round: Which filled voxels can be made empty  ....

    makeGroups();
    nSize1 = delGroup.size();

    for( size_t i = 0; i < nSize1; i++) {
        const JCellPtr &delcell = delGroup[i];
        int x = getChangeInSingularPointsAfterDeletion(delcell);
        size_t nSize2 = addGroup.size();
        for( size_t j = 0; j < nSize2; j++) {
            const JCellPtr &addcell = addGroup.front();
            addGroup.pop_front();
            int y = getChangeInSingularPointsAfterAddition(addcell);
            if( x + y < 0) {
                addcell->setStatus(JMeshEntity::ACTIVE);
                delcell->setStatus(JMeshEntity::INACTIVE);
                for( int k = 0; k < 8; k++)
                    addcell->getNodeAt(k)->addRelation(addcell);
                break;
            } else
                addGroup.push_back(addcell);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////

void JMonotoneVoxelizer :: setDistance()
{
    /*
       int numVoxels = cellDim[0]*cellDim[1]*cellDim[2];
       if( numVoxels == 0) return 0;

       for( size_t i = 0; i < numVoxels; i++)
            distance[i] =std::numeric_limits<float>::max();

       for( size_t i = 0; i < numVoxels; i++) {
           const JCellPtr &cell = mesh->getCellAt(i);
           if( cell->isActive() ) {
               distance[i] = 0.0;
               cellQ.push_back(cell);
           }
       }

       double currdist, nextdist;
       while(!cellQ.empty() ) {
              currCell = cellQ.front();
              cellQ.pop_front();
              currdist = distance[currCell->getID()];
              getRelations26(currCell, neighs);
              for( const JCellPtr &nextCell : neighs) {
                  nextdist = distance[nextCell->getID()];
                  nextdist = min( nextdist, currdist + 1.0);
                  distance[nextCell->getID()] = nextdist;
                  cellQ.push_back(nextCell);
              }
        }
    */
}
/////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////
/*
void JMonotoneVoxelizer :: setAddRank( const JCellPtr &cell)
{
    if( cell->isActive()  ) {
        addRank[0].push_back(cell);
        return;
    }

    // Temporary set the cell active and calculate number of singular points ...
    cell->setStatus( JMeshEntity::ACTIVE);
    int nsum = 0;
    for( int i = 0; i < 8; i++) {
         if( isSingular( cell->getNodeAt(i))) nsum++;
    }
    addRank[nsum].push_back(cell);

    // Reset the value ...
    cell->setStatus( JMeshEntity::INACTIVE);
}
/////////////////////////////////////////////////////////////////////////////////////////////

void JMonotoneVoxelizer :: setRemoveRank( const JCellPtr &cell)
{
    if( !cell->isActive()  ) {
        removeRank[0].push_back(cell);
        return;
    }

    // Temporary set the cell active and calculate number of singular points ...
    cell->setStatus( JMeshEntity::INACTIVE);
    int nsum = 0;
    for( int i = 0; i < 8; i++)
         if( isSingular( cell->getNodeAt(i))) nsum++;
    removeRank[nsum].push_back(cell);

    // Reset the value ...
    cell->setStatus( JMeshEntity::ACTIVE);
}
/////////////////////////////////////////////////////////////////////////////////////////////


void JMonotoneVoxelizer :: generateRanks()
{
   addRank.clear();
   removeRank.clear();
   if( voxmesh == nullptr) return;

   JMeshPtr bgmesh = voxmesh->getBackgroundMesh();
   vector<int> cellDim = voxmesh->getCellDimensions();
   size_t numCells = cellDim[0]*cellDim[1]*cellDim[2];
   for( int i = 0; i < numCells; i++) {
         const JCellPtr &cell = bgmesh->getCellAt(i);
         setAddRank(cell);
         setRemoveRank(cell);
   }
}
*/

