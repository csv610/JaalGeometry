#include "MeshPartitioner.hpp"
#include "MeshTopology.hpp"

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////
JLogger* JMeshPartitioner :: logger = JLogger::getInstance();

void JMeshPartitioner :: clear()
{
    elemPart.clear();
    nodePart.clear();
    edgeGroup.clear();

    mesh->deleteNodeAttribute("Partition");
    mesh->deleteNodeAttribute("Interface");
    mesh->deleteNodeAttribute("PartitionCorner");
    mesh->deleteEdgeAttribute("Interface");
    mesh->deleteFaceAttribute("Partition");
    mesh->deleteFaceAttribute("Interface");
    mesh->deleteCellAttribute("Partition");
}
///////////////////////////////////////////////////////////////////////////////
int JMeshPartitioner :: setElementsPartitionFromNodes()
{
    size_t numfaces = mesh->getSize(2);
    vector<int> nodeParts;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int nnodes = face->getSize(0);
        nodeParts.resize(nnodes);
        for( int j = 0; j < nnodes; j++) {
            const JNodePtr &vtx = face->getNodeAt(j);
            vtx->getAttribute("Partition", nodeParts[j] );
        }
        int faceID = JMath::high_frequency_item( nodeParts);
        face->setAttribute("Partition", faceID);
    }
}
///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner :: searchFaceComponents()
{
    std::map<int, JFaceSequence>  faceGroups;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        assert( face->hasAttribute("Partition") );
        face->setVisitBit(0);
        face->setID(i);
    }

    JFacePtr seedface, currface;
    deque<JFacePtr> faceQ;
    JFaceSequence fneighs;

    int compID = 0;
    int id1, id2;

//    size_t nCount = 0;
    while(1) {
        seedface.reset();
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->getVisitBit() == 0 && face->isActive() ) {
                seedface = face;
                break;
            }
        }
        if( seedface == nullptr) break;
        faceQ.push_back( seedface);
        int err = seedface->getAttribute("Partition", id1);
        assert( !err);
        while( !faceQ.empty() ) {
            currface = faceQ.front();
            faceQ.pop_front();
            if( currface->getVisitBit() == 0) {
                JFace::getRelations12(currface, fneighs);
                faceGroups[compID].push_back(currface);
                currface->setVisitBit(1);
                for( const JFacePtr &nextface : fneighs) {
                    nextface->getAttribute("Partition", id2);
                    if( (id1 == id2) && nextface->getVisitBit() == 0)
                        faceQ.push_back(nextface);
                }
            }
        }
        compID++;
    }

    for( auto keyVal : faceGroups) {
        compID = keyVal.first;
        for( const JFacePtr &face: faceGroups[compID])
            face->setAttribute("Partition", compID);
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////
int JMeshPartitioner :: searchComponents()
{
    return searchFaceComponents();
}
///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner ::searchRegion(const JFacePtr  &seedface)
{
    if( seedface == nullptr ) return 1;

    deque<JFacePtr> faceQ;
    JFaceSequence fneighs;

    int pid, seedid;
    int err = seedface->getAttribute( "Partition", seedid) ;
    if( err ) {
        cout << "Fatal: Seed surface partition id is not known " << endl;
        return 1;
    }

    faceQ.clear();
    faceQ.push_back(seedface);

    JFacePtr currface;
    while(!faceQ.empty() ) {
        currface = faceQ.front();
        faceQ.pop_front();
        if( !currface->isActive() ) continue;
        int numedges = currface->getSize(1);
        for( int i = 0; i < numedges; i++) {
            const JEdgePtr &edge = currface->getEdgeAt(i);
            if( !edge->hasAttribute("Interface") ) {
                JEdge::getRelations(edge, fneighs);
                for( const JFacePtr &neigh : fneighs) {
                    neigh->getAttribute("Partition", pid);
                    if( pid < 0) {
                        neigh->setAttribute( "Partition", seedid) ;
                        faceQ.push_back( neigh );
                    }
                }
            }
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner ::searchEdgesInterface()
{
    if( mesh == NULL ) return 1;

    mesh->deleteEdgeAttribute("Interface");

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() )
            if( !face->hasAttribute("Partition") ) {
                cout << "Warning: One of the active face does not have parition attribute" << endl;
                cout << "         Interface not identified" << endl;
                return 1;
            }
    }

    size_t numedges = mesh->getSize(1);
    std::map< pair<int,int>, int>  edgemap;

    JFaceSequence faceneighs;

    std::pair<int,int>  idpair;
    int part1, part2;
    int partid = 0;

    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        JEdge::getRelations(edge, faceneighs);
        if( faceneighs.size() == 2 ) {
            faceneighs[0]->getAttribute("Partition", part1);
            faceneighs[1]->getAttribute("Partition", part2);
            if( part1 != part2) {
                idpair.first  = min(part1, part2);
                idpair.second = max(part1, part2);
                if( edgemap.find(idpair) == edgemap.end() ) {
                    edgemap[idpair] = partid;
                    partid++;
                }
                edge->setAttribute("Interface", edgemap[idpair] );
            }
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner ::searchFacesInterface()
{
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner ::searchInterfaces()
{
    int topDim = mesh->getTopology()->getDimension();

    if( topDim == 2 ) searchEdgesInterface();
    if( topDim == 3 ) searchFacesInterface();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner ::searchRegions()
{
    //
    // If you are given interfaces, identify the faces/cells partitions
    //
    if( mesh == NULL ) return 1;
    logger->setInfo("Calculating regions in parition mesh");


    mesh->buildRelations(1, 2);
    mesh->deleteFaceAttribute("Interface");

    size_t nSize = mesh->getSize(2);
    int index = 0;

    int pid, partid = -1;
    for( size_t i = 0; i < nSize; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->setAttribute( "Partition", partid);
            face->setID( index++ );
        }
    }

    partid = 0;
    while(1) {
        JFacePtr seedface;
        for( size_t i = 0; i < nSize; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive()  ) {
                f->getAttribute("Partition", pid);
                if( pid < 0) {
                    seedface = f;
                    break;
                }
            }
        }
        if( seedface == nullptr ) break;
        seedface->setAttribute( "Partition", partid) ;
        searchRegion( seedface);
        partid++;
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
int JMeshPartitioner ::searchCorners()
{
    if( mesh == NULL ) return 1;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2) {
        cout << "Warning: Presently the corners are identified for surface mesh only " << endl;
        return 1;
    }

    mesh->deleteNodeAttribute("PartitionCorner");

    size_t nSize = mesh->getSize(2);
    for( size_t i = 0; i < nSize; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            if( !face->hasAttribute( "Partition") ) {
                cout << "Warning : One of the active face no partition attribute" << endl;
                cout << "          Corners not identified" << endl;
            }
        }
    }

    mesh->buildRelations(0, 2);

    nSize = mesh->getSize(0);
    JFaceSequence faces;
    set<int> aset;
    int pid;
    boost::any anyval;
    for( size_t i = 0; i < nSize; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        JNode::getRelations( vtx, faces);
        aset.clear();
        for( const JFacePtr &face : faces) {
            face->getAttribute("Partition", pid);
            aset.insert(pid);
        }
        if( aset.size() > 2 ) {
            vtx->setAttribute("PartitionCorner", anyval);
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////

void JMeshPartitioner :: getCorners( JNodeSequence  &nodes)
{
    nodes.clear();
    if( mesh == nullptr) return;

    size_t numnodes = mesh->getSize(0);

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx  = mesh->getNodeAt(i);
        if( vtx->isActive()  ) {
            if( vtx->hasAttribute("PartitionCorner") ) {
                nodes.push_back(vtx);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner :: getNumPartitions() const
{
    if( mesh == nullptr) return 0;
    set<int> pset;

    int topDim = mesh->getTopology()->getDimension();

    int pid;

    if( topDim == 2 ) {
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f  =  mesh->getFaceAt(i);
            if( f->isActive() ) {
                int err = f->getAttribute("Partition", pid );
                if( !err ) pset.insert(pid);
            }
        }
    }

    if( topDim == 3 ) {
        size_t numcells = mesh->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &c  =  mesh->getCellAt(i);
            if( c->isActive() ) {
                int err = c->getAttribute("Partition", pid );
                if( !err) pset.insert(pid);
            }
        }
    }

    return pset.size();
}

///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner :: getNumInterfaces() const
{
    if( mesh == nullptr) return 0;

    set<int> pset;

    int topDim = mesh->getTopology()->getDimension();

    int pid;

    if( topDim == 2 ) {
        size_t numedges = mesh->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge  =  mesh->getEdgeAt(i);
            int err = edge->getAttribute("Interface", pid );
            if( !err) pset.insert(pid);
        }
    }

    if( topDim == 3 ) {
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face  =  mesh->getFaceAt(i);
            int err = face->getAttribute("Interface", pid );
            if( !err) pset.insert(pid);
        }
    }
    return pset.size();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshPartitioner :: getInterface(JEdgeSequence &iedges)
{
    iedges.clear();
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive()  && edge->hasAttribute("Interface") )
            iedges.push_back(edge);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPartitioner :: getInterface(int id, JEdgeSequence &iedges)
{
    if( edgeGroup.empty() ) {
        size_t numedges = mesh->getSize(1);
        int val;
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                int err = edge->getAttribute("Interface", val);
                if( !err) edgeGroup[val].push_back(edge);
            }
        }
    }

    iedges = edgeGroup[id];
}
///////////////////////////////////////////////////////////////////////////////
int JMeshPartitioner :: getNumCorners() const
{
    if( mesh == nullptr) return 0;

    int nCount = 0;
    nCount = mesh->getNumAttributes("PartitionCorner", 0);
    return nCount;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshPartitioner :: getGraph()
{
    JMeshPtr partGraph = JMesh::newObject();

    map<int,JNodePtr> vmap;

    if( mesh == nullptr ) return partGraph;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim == 2 ) {
        size_t numedges = mesh->getSize(1);

        std::pair<int,int>  ipair;
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                int err = edge->getAttribute("Interface", ipair);
                if( !err) {
                    int p1 = ipair.first;
                    int p2 = ipair.second;
                    if( vmap.find(p1) == vmap.end() ) {
                        JNodePtr vnew = JNode::newObject();
                        vnew->setAttribute("Partition", p1);
                        vmap[p1] = vnew;
                        partGraph->addObject( vnew );
                    }

                    if( vmap.find(p2) == vmap.end() ) {
                        JNodePtr vnew = JNode::newObject();
                        vnew->setAttribute("Partition", p2);
                        vmap[p2] = vnew;
                        partGraph->addObject( vnew );
                    }
                    if( JSimplex::getEdgeOf( vmap[p1], vmap[p2], 0) == nullptr ) {
                        JEdgePtr edge = JEdge::newObject( vmap[p1], vmap[p2] );
                        partGraph->addObject(edge);
                    }
                }
            }
        }
    }
    return partGraph;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshPartitioner :: getSubMesh( int partid)
{
    if( mesh == nullptr) return nullptr;

    size_t numfaces = mesh->getSize(2);

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            if( !face->hasAttribute("Partition")  )  {
                cout << "Warning: Face does not have Partition attribute" << endl;
                return nullptr;
            }
        }
    }

    JNodeSet vset;
    JFaceSequence  pfaces;

    int pid;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int err = face->getAttribute("Partition", pid );
            if( !err) {
                if( pid == partid) {
                    pfaces.push_back(face);
                    int nnodes = face->getSize(0);
                    for( int j = 0; j < nnodes; j++)
                        vset.insert( face->getNodeAt(j));
                }
            }
        }
    }

    JMeshPtr sbmesh = JMesh::newObject();
    for(JNodePtr vertex : vset) sbmesh->addObject(vertex);

    numfaces = pfaces.size();
    for( size_t i = 0; i < numfaces; i++)
        sbmesh->addObject( pfaces[i] );
    return sbmesh;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshPartitioner :: getPartitions()
{
    if( mesh == nullptr ) return 1;

    searchInterfaces();
    searchCorners();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JMeshPartitioner :: getRegionBoundary( int partid, JEdgeSequence &edges)
{
    edges.clear();
    if( mesh == nullptr ) return 1;

    JMeshPtr submesh = this->getSubMesh(partid);
    size_t numfaces = submesh->getSize(2);

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = submesh->getFaceAt(i);
        if( !face->hasAttribute("Partition" ) ) {
            cout << "Warning: Face does not have Partition attribute" << endl;
            return 1;
        }
        int   numedges = face->getSize(1);
        for( int j = 0; j < numedges; j++) {
            const JEdgePtr &edge = face->getEdgeAt(j);
            if( edge->hasAttribute("Interface") ) {
                edges.push_back ( edge );
            }
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshPartitioner :: getPartition( int partid, JFaceSequence &faces)
{
    faces.clear();
    if( mesh ==  nullptr) return;

    int val;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int err =  face->getAttribute("Partition", val);
            if( err == 0 && partid == val )
                faces.push_back(face);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPartitioner :: removeInterface( int breakid)
{
    if( mesh == nullptr) return;

    size_t numedges = mesh->getSize(1);
    size_t numfaces = mesh->getSize(2);

    JEdgeSequence edges;
    JFaceSequence faces;
    getInterface(breakid, edges);

    int err;
    if( !edges.empty() )  {
        int id;
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            int err = edge->getAttribute("Interface", id);
            if( err == 0 && id > breakid)
                edge->setAttribute("Interface", id-1);
        }
        for( const JEdgePtr &edge : edges)
            edge->deleteAttribute("Interface");

        JFaceSequence faces;
        JEdge::getRelations(edges[0], faces);

        if( faces.size() == 1) return;

        if( faces.size() > 2) {
            cout << "Fatal error: Non manifold edge detected" << endl;
            return;
        }

        int fid1;
        err = faces[0]->getAttribute("Partition", fid1);
        if( err ) {
            cout << "Error: A face has no Partition attribute" << endl;
            return;
        }
        int fid2;

        err = faces[1]->getAttribute("Partition", fid2);
        if( err ) {
            cout << "Error: A face has no Partition attribute" << endl;
            return;
        }
        int maxid = max(fid1, fid2);
        int minid = min(fid1, fid2);

        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            int err = face->getAttribute("Partition", id);
            if( err == 0 && id > minid)
                face->setAttribute("Partition", id-1);
        }
        getPartition(maxid, faces);
        for( const JFacePtr &face : faces)
            face->setAttribute("Partition", minid);
    }
}

///////////////////////////////////////////////////////////////////////////////////
void JMeshPartitioner :: savePartitions( int indextype )
{
    int nparts = getNumPartitions();

    if( indextype == 0) {

        JMeshIO mio;
        for( int i = 0; i < nparts; i++) {
            JMeshPtr meshpart = getSubMesh(i);
            meshpart->enumerate(0);
            string  filename = "meshpart" + to_string(i) + ".off";
            mio.saveAs( meshpart, filename);
        }
        mesh->enumerate(0);  // restore back
    }
}
