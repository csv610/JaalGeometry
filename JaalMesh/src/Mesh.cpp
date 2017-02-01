#include <iomanip>

#include "Mesh.hpp"
#include "MeshImporter.hpp"
#include "MeshExporter.hpp"

#include "MeshGeometry.hpp"
#include "MeshTopology.hpp"

#include "basic_math.hpp"
#include <complex>

using namespace std;
using namespace Jaal;

size_t JMesh::nCounter = 0;

JLogger* JMesh::logger =  JLogger::getInstance();

JLogger* JMesh::getLogger()
{
    if( logger ) logger = JLogger::getInstance();
    return logger;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JMesh::newObject()
{
    JMeshPtr msh(new JMesh);

    JMeshGeometryPtr geom( new JMeshGeometry(msh));
    JMeshTopologyPtr topo( new JMeshTopology(msh));
    msh->setGeometry(geom);
    msh->setTopology(topo);

    ostringstream oss;
    oss << "Mesh" << nCounter++;
    msh->setName( oss.str() );
    return msh;
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMesh::newObject( JEdgeSequence &edges)
{
    JMeshPtr newmesh = newObject();
    JNodeSequence nodes;
    JMeshTopology::getEntitySet( edges, nodes);
    newmesh->addObjects(nodes);
    newmesh->addObjects(edges);
    return newmesh;
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMesh::newObject( JFaceSequence &faces)
{
    JMeshPtr newmesh = newObject();
    JNodeSequence nodes;
    JMeshTopology::getEntitySet( faces, nodes);
    newmesh->addObjects(nodes);
    newmesh->addObjects(faces);
    return newmesh;
}
///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMesh::newObject( JCellSequence &cells)
{
    JMeshPtr newmesh = newObject();
    JNodeSequence nodes;
    JMeshTopology::getEntitySet( cells, nodes);
    newmesh->addObjects(nodes);
    newmesh->addObjects(cells);
    return newmesh;
}

///////////////////////////////////////////////////////////////////////////////

void JMesh :: deleteNodeAttribute(const string &s)
{
    size_t nSize = nodes.size();
    for( size_t i = 0; i < nSize; i++)
        nodes[i]->deleteAttribute(s);
}
///////////////////////////////////////////////////////////////////////////////
JNodeSequence JMesh :: getNodes() const
{
    JNodeSequence seq;
    logger->setInfo("Collecting active nodes in the mesh");

    size_t nSize = nodes.size();
    if( nSize) {
        seq.reserve( nSize );

        for( size_t i = 0; i <  nSize; i++) {
            JNodePtr v = getNodeAt(i);
            if( v->isActive() ) seq.push_back(v);
        }
    }
    return seq;
}
///////////////////////////////////////////////////////////////////////////////
JNodeSequence JMesh :: getElementsNodes() const
{
    JNodeSequence seq;

    if( !cells.empty() ) {
        JMeshTopology::getEntitySet( cells, seq);
        return seq;
    }

    if( !faces.empty() ) {
        JMeshTopology::getEntitySet( faces, seq);
        return seq;
    }

    if( !edges.empty() ) {
        JMeshTopology::getEntitySet( edges, seq);
        return seq;
    }

    return nodes;
}
///////////////////////////////////////////////////////////////////////////////



void JMesh :: setVisitBits( int e, bool val)
{
    if( e == 0) {
        size_t numnodes =  nodes.size();
        for( size_t i = 0; i < numnodes; i++)
            nodes[i]->setVisitBit(val);
    }

    if( e == 1) {
        size_t numedges =  edges.size();
        for( size_t i = 0; i < numedges; i++)
            edges[i]->setVisitBit(val);
    }

    if( e == 2) {
        size_t numfaces =  faces.size();
        for( size_t i = 0; i < numfaces; i++)
            faces[i]->setVisitBit(val);
    }

    if( e == 3) {
        size_t numcells =  cells.size();
        for( size_t i = 0; i < numcells; i++)
            cells[i]->setVisitBit(val);
    }
}
//////////////////////////////////////////////////////////////////////////////////

JNodeSequence JMesh::getNodes(const string &attribname) const
{
    JNodeSequence seq;

    size_t nSize = nodes.size();

    if( nSize ) {
        seq.reserve( nSize );
        for( size_t i = 0; i <  nSize; i++) {
            JNodePtr v = getNodeAt(i);
            if( v->isActive() && v->hasAttribute(attribname) )
                seq.push_back(v);
        }
    }
    return seq;
}

///////////////////////////////////////////////////////////////////////////////

void JMesh :: initialize()
{
    logger = JLogger::getInstance();
    logger->setInfo("Initialize new mesh");

    geomModel = 0;
    name = "unknown";

    activeBit  = 1;
    status     = ACTIVE;
    topoChanged[0] = 0;
    topoChanged[1] = 0;
    topoChanged[2] = 0;
    topoChanged[3] = 0;

    geomChanged = 0;

    keep_mesh_entity[0] = 0;
    keep_mesh_entity[1] = 0;
    keep_mesh_entity[2] = 0;
    keep_mesh_entity[3] = 0;
    relManager.reset( new JMeshRelationManager(this));
    attribManager.reset( new JAttributeManager);

    srand48( time(0));   // For random entity selection ...

}
///////////////////////////////////////////////////////////////////////////////

void JMesh :: finalize()
{
    logger->setInfo("Cleating mesh data structures ");
    clearAll();
}

///////////////////////////////////////////////////////////////////////////////

size_t JMesh :: getSize(int d)
{
    switch(d) {
    case 0:
        return nodes.size();
    case 1:
        return edges.size();
    case 2:
        return faces.size();
    case 3:
        return cells.size();
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JMesh :: getEdges()
{
    JEdgeSequence seq;

    size_t nSize = edges.size();
    if( nSize ) {
        seq.reserve( nSize );

        for( size_t i = 0; i <  nSize; i++) {
            const JEdgePtr &e = getEdgeAt(i);
            if( e->isActive() ) seq.push_back(e);
        }
    }
    return seq;
}
///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JMesh :: getElementsEdges() const
{
    JEdgeSequence seq;

    if( !cells.empty() ) {
        JMeshTopology::getEntitySet( cells, seq);
        return seq;
    }

    if( !faces.empty() ) {
        JMeshTopology::getEntitySet( faces, seq);
        return seq;
    }

    return edges;
}

///////////////////////////////////////////////////////////////////////////////
void JMesh:: deleteEdgeAttribute(const string &s)
{
    size_t nSize = edges.size();
    for( size_t i = 0; i < nSize; i++)
        edges[i]->deleteAttribute(s);
}

void JMesh::generate_objects(size_t n, JNodeSequence &seq)
{
    seq.clear();
    if(n < 1) return;

    seq.resize(n);
    for (size_t i = 0; i < n; i++)
        seq[i] = JNode::newObject();
}

///////////////////////////////////////////////////////////////////////////////

int JMesh :: getPosOf( const JNodePtr &vsearch, size_t &pos) const
{
    // If the mesh is pruned and enumerated then finding position is easy
    size_t id = vsearch->getID();
    if( nodes[id] == vsearch) {
        pos = id;
        return 1;
    }

    pos = std::numeric_limits<size_t>::max() - 1;
    size_t nSize = nodes.size();
    for( size_t i = 0; i < nSize; i++) {
        if( nodes[i] == vsearch ) {
            pos =  i;
            return 1;
        }
    }

    return -1;
}

///////////////////////////////////////////////////////////////////////////////

int JMesh :: swap( const JNodePtr &v1, const JNodePtr &v2)
{
    size_t pos1, pos2;
    int    found;

    found = getPosOf(v1, pos1);
    if( found < 0) return 1;

    found = getPosOf(v2, pos2);
    if( found < 0) return 1;

    swap( nodes[pos1], nodes[pos2] );

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMesh :: getPosOf( const JEdgePtr &esearch, size_t &pos) const
{
    // If the mesh is pruned and enumerated then finding position is easy
    size_t id = esearch->getID();
    if( edges[id] == esearch) {
        pos = id;
        return 1;
    }

    pos = std::numeric_limits<size_t>::max() - 1;
    size_t nSize = edges.size();
    for( size_t i = 0; i < nSize; i++) {
        if( edges[i] == esearch ) {
            pos =  i;
            return 1;
        }
    }
    return -1;
}

///////////////////////////////////////////////////////////////////////////////

int JMesh :: swap( const JEdgePtr &e1, const JEdgePtr &e2)
{
    size_t pos1, pos2;
    int  found;

    found = getPosOf(e1, pos1);
    if( found < 0) return 1;

    found = getPosOf(e2, pos2);
    if( found < 0) return 1;

    swap( edges[pos1], edges[pos2] );

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void JMesh::generate_objects(size_t n, JEdgeSequence &seq)
{
    seq.clear();
    if( n < 1) return;

    seq.resize(n);
    for (size_t i = 0; i < n; i++)
        seq[i] = JEdge::newObject();
}

///////////////////////////////////////////////////////////////////////////////

void JMesh::generate_objects(size_t n, JFaceSequence &seq, int type)
{
    seq.clear();

    if( n < 1 ) return;

    seq.resize(n);

    switch( type ) {
    case JFace::TRIANGLE:
        for (size_t i = 0; i < n; i++)
            seq[i] = JTriangle::newObject();
        break;
    case JFace::QUADRILATERAL:
        for (size_t i = 0; i < n; i++)
            seq[i] = JQuadrilateral::newObject();
        break;
    case JFace::POLYGON:
        for (size_t i = 0; i < n; i++)
            seq[i] = JPolygon::newObject();
        break;
    default:
        cout << "No such face object available " << endl;
        exit(0);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMesh::generate_objects(size_t n, JCellSequence &seq, int type)
{
    seq.clear();

    if( n < 1 ) return;

    seq.resize(n);

    switch( type ) {
    case JCell::TETRAHEDRON:
        for (size_t i = 0; i < n; i++)
            seq[i] = JTetrahedron::newObject();
        break;
    case JCell::HEXAHEDRON:
        for (size_t i = 0; i < n; i++)
            seq[i] = JHexahedron::newObject();
        break;
    case JCell::POLYHEDRON:
        for (size_t i = 0; i < n; i++)
            seq[i] = JPolyhedron::newObject();
        break;
    default:
        cout << "No such cell object available " << endl;
        exit(0);
    }
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr
JMesh::deepCopy()
{
    logger->setInfo("Making deep copy of the mesh");
    JMeshPtr newmesh = JMesh::newObject();

    std::map<JNodePtr, JNodePtr> vmap;

    size_t numnodes = nodes.size();
    if( numnodes ) {
        newmesh->reserve(numnodes, 0);
        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vold = getNodeAt(i);
            if( vold->isActive() ) {
                JNodePtr vnew = JNode::newObject();
                vnew->setXYZCoords(vold->getXYZCoords());
                vmap[vold] = vnew;
                newmesh->addObject(vnew);
            }
        }
    }

    size_t numfaces = faces.size();
    if( numfaces ) {
        newmesh->reserve(numfaces, 2);
        JNodeSequence nodes;
        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &fold = getFaceAt(i);
            if( fold->isActive() ) {
                nodes  = fold->getNodes();
                int nv   = nodes.size();
                for (int j = 0; j < nv; j++)
                    nodes[j] = vmap[nodes[j]];
                JFacePtr fnew = fold->getClone();
                fnew->setNodes(nodes);
                newmesh->addObject(fnew);
            }
        }
    }

    size_t numcells = cells.size();
    if( numcells ) {
        newmesh->reserve(numcells, 3);
        JNodeSequence nodes;
        for (size_t i = 0; i < numcells; i++) {
            const JCellPtr &cold = getCellAt(i);
            if( cold->isActive() ) {
                nodes  = cold->getNodes();
                int nv   = nodes.size();
                for (int j = 0; j < nv; j++)
                    nodes[j] = vmap[nodes[j]];
                JCellPtr cnew = cold->getClone();
                cnew->setNodes(nodes);
                newmesh->addObject(cnew);
            }
        }
    }

    return newmesh;
}

///////////////////////////////////////////////////////////////////////////////
void JMesh::addObjects(const JNodeSequence  &seq, bool embedded )
{
    for(const JNodePtr &v : seq) addObject( v, embedded);
}

///////////////////////////////////////////////////////////////////////////////

int JMesh :: addObject(const JEdgePtr &edge, bool embedded)
{
    if( edge == nullptr ) return 1;

    if( !edge->isActive() ) {
        logger->setWarn("Only active edges can be added to the container");
        return 1;
    }

    // Attach the new edge at the hash vertex...
    const JNodePtr &vhash = edge->getHashNode();
    if( !edge->hasActiveNodes() ) {
        logger->setWarn("An active edge should have all active nodes");
        return 1;
    }

    if( vhash ) vhash->attach(edge);

    topoChanged[1] = 1;

    edge->setStatus(JMeshEntity::ACTIVE);

    if( embedded ) {
        embedded_edges.push_back(edge);
        return 0;
    }

    edges.push_back(edge);
    relManager->buildRelations(edge);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
void JMesh::addObjects(const JEdgeSequence  &seq, bool embedded)
{
    for(const JEdgePtr &e : seq) addObject( e, embedded);
}
////////////////////////////////////////////////////////////////////////////////

int JMesh::addObject(const JFacePtr &newface, bool embedded)
{
    if( newface == nullptr ) return 1;

    if( !newface->isActive() ) {
        logger->setWarn("Only an active face can be added to the container");
        return 1;
    }

    // Basic stuff ....
    if( !newface->hasActiveNodes() ) {
        logger->setWarn("An active face must have all active nodes ");
        return 1;
    }

    if( !newface->hasUniqueNodes() ) {
        cout << "Warning: Face has non-unique nodes: not added to mesh" << endl;
        int nn = newface->getSize(0);
        for( int i = 0; i < nn; i++)
            cout << newface->getNodeAt(i)->getID() << " ";
        cout << endl;
        return 1;
    }

    newface->setStatus(JMeshEntity::ACTIVE);

    if( embedded ) {
        embedded_faces.push_back(newface);
        return 0;
    }

    topoChanged[2] = 1;

    const JNodePtr &vhash = newface->getHashNode();
    if( vhash ) vhash->attach(newface);

    faces.push_back(newface);

    int numedges = newface->getSize(0);
    JEdgePtr edge;
    for( int i = 0; i < numedges; i++) {
        if( !newface->hasEdgeAt(i) ) {
            edge = newface->getEdgeAt(i);
            addObject(edge);
        } else
            edge = newface->getEdgeAt(i);

        edge->addRelation(newface);
    }

    relManager->buildRelations(newface);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
void JMesh :: addObjects( const JFaceSequence &seq)
{
    for(const JFacePtr &f : seq) addObject(f);
}
///////////////////////////////////////////////////////////////////////////////
JFaceSequence JMesh :: getFaces() const
{
    JFaceSequence seq;

    size_t nSize = faces.size();
    if( nSize ) {
        seq.reserve( nSize );
        for( size_t i = 0; i <  nSize; i++) {
            const JFacePtr &f = getFaceAt(i);
            if( f->isActive() ) seq.push_back(f);
        }
    }
    return seq;
}
////////////////////////////////////////////////////////////////////////////////

JFaceSequence JMesh :: getElementsFaces() const
{
    JFaceSequence seq;
    if( !cells.empty() ) {
        JMeshTopology::getEntitySet( cells, seq);
        return seq;
    }
    return faces;
}
///////////////////////////////////////////////////////////////////////////////

JFaceSequence JMesh :: getFaces(const string &attribname) const
{
    JFaceSequence seq;

    size_t nSize = faces.size();
    if( nSize ) {
        seq.reserve( nSize );
        for( size_t i = 0; i <  nSize; i++) {
            const JFacePtr &f = getFaceAt(i);
            if( f->isActive() && f->hasAttribute(attribname) )
                seq.push_back(f);
        }
    }
    return seq;
}

////////////////////////////////////////////////////////////////////////////////

int JMesh ::addObject(const JCellPtr &newcell)
{
    if( newcell == nullptr ) return 1;

    if( !newcell->isActive() ) {
        logger->setWarn("Only an active cell can be added to the container");
        return 1;
    }
    // Basic stuff ....
    if( !newcell->hasActiveNodes() ) {
        logger->setWarn("An active face must have all active nodes");
        return 1;
    }

    topoChanged[3] = 1;

    const JNodePtr &vhash = newcell->getHashNode();
    if( vhash ) vhash->attach(newcell);

    cells.push_back(newcell);

    int numfaces = newcell->getSize(2);

    JFacePtr face = nullptr;
    for( int i = 0; i < numfaces; i++) {
        if( !newcell->hasFaceAt(i) ) {
            face = newcell->getFaceAt(i);
            addObject(face);
        } else {
            face = newcell->getFaceAt(i);
        }
        face->addRelation(newcell);
    }
    relManager->buildRelations(newcell);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JMesh :: addObjects( const JCellSequence &seq)
{
    for(const JCellPtr &e : seq) addObject(e);
}

///////////////////////////////////////////////////////////////////////////////
int JMesh::remove(const JNodePtr &v)
{
    if (v == nullptr) return 1;
    if (v->isRemoved()) return 0;  // Already removed
    v->setStatus(JMeshEntity::REMOVE);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMesh::remove(const JEdgePtr &edge)
{
    if (edge == nullptr) return 1;
    if (edge->isRemoved()) return 0;    // Already removed ...
    topoChanged[1] = 1;

    edge->setStatus(JMeshEntity::REMOVE);
    relManager->removeRelations(edge);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMesh::remove(const JFacePtr &face)
{
    if( face == nullptr ) return 0;
    if( face->isRemoved() ) return 0;   // Already removed ...
    topoChanged[2] = 1;
    face->setStatus(JMeshEntity::REMOVE);
    relManager->removeRelations(face);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMesh :: remove(const JCellPtr &cell)
{
    if( cell == nullptr ) return 0;
    if( !cell->isActive() ) return 0;   // Already removed ...

    topoChanged[3] = 1;
    cell->setStatus(JMeshEntity::REMOVE);
    relManager->removeRelations(cell);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int  JMesh :: getAttributeType(const string &name, int entity) const
{
    assert( entity >=0 && entity < 4);

    size_t nCount = 0;

    switch( entity ) {
    case 0:
        for( const JNodePtr &vertex : nodes) {
            if( vertex->isActive() ) {
                if( !vertex->hasAttribute(name) && nCount > 1)
                    return JAttribute::SPARSE_ATTRIBUTE;
                nCount++;
            }
        }
        if( nCount ) return JAttribute::DENSE_ATTRIBUTE;
        break;
    case 1:
        for( const JEdgePtr &edge : edges) {
            if( edge->isActive() ) {
                if( !edge->hasAttribute(name) && nCount > 1)
                    return JAttribute::SPARSE_ATTRIBUTE;
                nCount++;
            }
        }
        if( nCount ) return JAttribute::DENSE_ATTRIBUTE;
        break;
    case 2:
        for( const JFacePtr &face: faces) {
            if( face->isActive() ) {
                if( !face->hasAttribute(name) && nCount > 1)
                    return JAttribute::SPARSE_ATTRIBUTE;
                nCount++;
            }
        }
        if( nCount ) return JAttribute::DENSE_ATTRIBUTE;
        break;
    case 3:
        for( const JCellPtr &cell: cells) {
            if( cell->isActive() ) {
                if( !cell->hasAttribute(name) && nCount > 1)
                    return JAttribute::SPARSE_ATTRIBUTE;
                nCount++;
            }
        }
        if( nCount ) return JAttribute::DENSE_ATTRIBUTE;
        break;
    }

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

int JMesh :: getNumAttributes(int entity ) const
{
    assert( entity >=0 && entity < 4);

    int nCount = 0;

    switch( entity ) {
    case 0:
        for( const JNodePtr &vertex: nodes) {
            if( vertex->isActive() )
                nCount = max( nCount, vertex->getNumAttributes());
        }
        break;
    case 1:
        for(const JEdgePtr &edge: edges) {
            if( edge->isActive() )
                nCount = max( nCount, edge->getNumAttributes());
        }
        break;
    case 2:
        for(const JFacePtr &face: faces) {
            if( face->isActive() )
                nCount = max( nCount, face->getNumAttributes());
        }
        break;
    case 3:
        for(const JCellPtr &cell: cells) {
            if( cell->isActive() )
                nCount = max( nCount, cell->getNumAttributes());
        }
        break;
    }
    return nCount;
}

////////////////////////////////////////////////////////////////////////////////
size_t JMesh :: getNumAttributes(const string &name, int entity ) const
{
    assert( entity >=0 && entity < 4);

    int nCount = 0;

    switch( entity ) {
    case 0:
        for(const JNodePtr &vertex: nodes) {
            if( vertex->isActive() && vertex->hasAttribute(name)) nCount++;
        }
        break;
    case 1:
        for(const JEdgePtr &edge: edges) {
            if( edge->isActive() && edge->hasAttribute(name)) nCount++;
        }
        break;
    case 2:
        for( const JFacePtr &face: faces) {
            if( face->isActive() && face->hasAttribute(name)) nCount++;
        }
        break;
    case 3:
        for( const JCellPtr &cell: cells) {
            if( cell->isActive() && cell->hasAttribute(name)) nCount++;
        }
        break;
    }
    return nCount;
}

///////////////////////////////////////////////////////////////////////////////

int JMesh :: getAttributeNames( vector<string> &vnames, int entity) const
{
    assert( entity >=0 && entity < 4);

    vnames.clear();
    set<string>  nameset;

    switch( entity ) {
    case 0:
        for( const JNodePtr &vertex: nodes) {
            if( vertex->isActive() ) {
                vertex->getAttributeNames(vnames);
                nameset.insert( vnames.begin(), vnames.end() );
            }
        }
        break;
    case 1:
        for(const JEdgePtr &edge: edges) {
            if( edge->isActive() ) {
                edge->getAttributeNames(vnames);
                nameset.insert( vnames.begin(), vnames.end() );
            }
        }
        break;
    case 2:
        for(const JFacePtr &face: faces) {
            if( face->isActive() ) {
                face->getAttributeNames(vnames);
                nameset.insert( vnames.begin(), vnames.end() );
            }
        }
        break;
    case 3:
        for(const JCellPtr &cell: cells) {
            if( cell->isActive() ) {
                cell->getAttributeNames(vnames);
                nameset.insert( vnames.begin(), vnames.end() );
            }
        }
        break;
    }

    vnames.clear();

    set<string>::const_iterator it;
    for( it = nameset.begin(); it != nameset.end(); ++it)
        vnames.push_back(*it);
    return 0;
}

///////////////////////////////////////////////////////////////////////////

bool JMesh :: hasAttribute(const string &name, int entity) const
{
    assert( entity >=0 && entity < 4);

    switch( entity ) {
    case 0:
        for(const JNodePtr &vertex: nodes) {
            if( vertex->isActive() && vertex->hasAttribute(name) ) return 1;
        }
        break;
    case 1:
        for(const JEdgePtr &edge: edges) {
            if( edge->isActive() && edge->hasAttribute(name) ) return 1;
        }
        break;
    case 2:
        for(const JFacePtr &face: faces) {
            if( face->isActive() && face->hasAttribute(name) ) return 1;
        }
        break;
    case 3:
        for( const JCellPtr &cell : cells) {
            if( cell->isActive() && cell->hasAttribute(name) ) return 1;
        }
        break;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////

size_t JMesh :: getActiveSize( int d )  const
{
    size_t nCount = 0;

    if( d == 0) {
        for(const JNodePtr &vtx: nodes)
            if( vtx->isActive() ) nCount++;
        return nCount;
    }

    if( d == 1) {
        for(const JEdgePtr &edge: edges)
            if( edge->isActive()  ) nCount++;
        return nCount;
    }

    if( d == 2) {
        for(const JFacePtr &face: faces)
            if(face->isActive() ) nCount++;
        return nCount;
    }

    if( d == 3) {
        for(const JCellPtr &cell: cells)
            if( cell->isActive() ) nCount++;
        return nCount;
    }

    return nCount;
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::enumerate(int etype)
{
    size_t index = 0;

    switch( etype ) {
    case 0:
        for( const JNodePtr &vertex: nodes)
            if( vertex->isActive() ) vertex->setID(index++);
        break;
    case 1:
        for( const JEdgePtr &edge: edges)
            if( edge->isActive() ) edge->setID(index++);
        break;
    case 2:
        for( const JFacePtr &face: faces)
            if( face->isActive() ) face->setID(index++);
        break;
    case 3:
        for( const JCellPtr &cell: cells)
            if( cell->isActive() ) cell->setID(index++);
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::progressiveSort(int etype)
{
    size_t index = 0;
    getTopology()->searchBoundary();
    if( etype == 0) {
        for( const JNodePtr &vertex: nodes)
            if( vertex->isActive() && vertex->isBoundary() ) vertex->setID(index++);
        for( const JNodePtr &vertex: nodes)
            if( vertex->isActive() && !vertex->isBoundary() ) vertex->setID(index++);
        std::sort(nodes.begin(), nodes.end(),
                  []( const JNodePtr &v1,  const JNodePtr &v2)
        {
            return v1->getID() < v2->getID();
        });
    }

    if( etype == 1) {
        for( const JEdgePtr &edge : edges)
            if( edge->isActive() && edge->isBoundary() ) edge->setID(index++);
        for( const JEdgePtr &edge : edges)
            if( edge->isActive() && !edge->isBoundary() ) edge->setID(index++);
        std::sort(edges.begin(), edges.end(),
                  []( const JEdgePtr &e1,  const JEdgePtr &e2)
        {
            return e1->getID() < e2->getID();
        });
    }

    if( etype == 3) {
        for( const JFacePtr &face : faces)
            if( face->isActive() && face->isBoundary() ) face->setID(index++);
        for( const JFacePtr &face : faces)
            if( face->isActive() && !face->isBoundary() ) face->setID(index++);
        std::sort(faces.begin(), faces.end(),
                  []( const JFacePtr &f1,  const JFacePtr &f2)
        {
            return f1->getID() < f2->getID();
        });
    }
}
///////////////////////////////////////////////////////////////////////////////
int
JMesh::sort( const JEntityCompare *cmp, int entity)
{
    /*
        if( cmp  == nullptr)  {
            cout << "Warning: No proxy class defined for sorting the nodes " << endl;
            return 1;
        }
    */

    if( entity == 0) {
        std::sort(nodes.begin(), nodes.end(),
                  [cmp]( const JNodePtr &v1,  const JNodePtr &v2)
        {
            return cmp->isLess(v1,v2);
        });
        /*
                std::sort(nodes.begin(), nodes.end(),
                          []( const JNodePtr &v1,  const JNodePtr &v2)
                                 { return v1->getID() < v2->getID();});
                std::sort(nodes.begin(), nodes.end(), JEntityIDCompare() );
        */
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

bool
JMesh::isPruned() const
{
    for (size_t i = 0; i < nodes.size(); i++)
        if (!nodes[i]->isActive()) return 0;

    for (size_t i = 0; i < edges.size(); i++)
        if (!edges[i]->isActive()) return 0;

    for (size_t i = 0; i < faces.size(); i++)
        if (!faces[i]->isActive()) return 0;

    for (size_t i = 0; i < cells.size(); i++)
        if (!cells[i]->isActive()) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::pruneNodes()
{
    logger->setInfo("Prunning node container");

    JNodeSequence::iterator vend;
    vend = remove_if(nodes.begin(), nodes.end(), EntityRemovedPred());
    nodes.erase(vend, nodes.end());

    enumerate(0);
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::pruneEdges()
{
    logger->setInfo("Prunning edge container ");

    JEdgeSequence::iterator eend;
    eend = remove_if(edges.begin(), edges.end(), EntityRemovedPred());
    edges.erase(eend, edges.end());

    enumerate(1);
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::pruneFaces()
{
    logger->setInfo("Prunning face container ");

    JFaceSequence::iterator fend;
    fend = remove_if(faces.begin(), faces.end(), EntityRemovedPred());
    faces.erase(fend, faces.end());

    enumerate(2);
}
///////////////////////////////////////////////////////////////////////////////

void
JMesh::pruneCells()
{
    logger->setInfo("Prunning cell container");
    if( cells.empty() ) return;

    JCellSequence::iterator cend;
    cend = remove_if(cells.begin(), cells.end(), EntityRemovedPred());
    cells.erase(cend, cells.end());

    enumerate(3);
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::pruneAll()
{
    logger->setInfo("Prunning entire mesh");
    pruneCells();
    pruneFaces();
    pruneEdges();
    pruneNodes();
}

///////////////////////////////////////////////////////////////////////////////
JCellSequence JMesh :: getCells() const
{
    JCellSequence seq;

    size_t nSize = cells.size();
    if( nSize) {
        seq.reserve( nSize );
        for( size_t i = 0; i <  nSize; i++) {
            const JCellPtr c = getCellAt(i);
            if( c->isActive() ) seq.push_back(c);
        }
    }
    return seq;
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::deleteSurfaceMesh()
{
    size_t numCells = cells.size();
    if( numCells) {
        cout << "Warning: Volume mesh present: surface mesh can not be removed" << endl;
        return;
    }

    size_t numFaces = faces.size();
    if( numFaces == 0) return;

    logger->setInfo("deleting volumetric mesh");

    size_t numNodes = nodes.size();
    size_t numEdges = edges.size();

    for( size_t i = 0; i < numNodes; i++) nodes[i]->setVisitBit(0);
    for( size_t i = 0; i < numEdges; i++) edges[i]->setVisitBit(0);

    JFaceSequence neighs;
    for( size_t i = 0; i < numEdges; i++) {
        JEdge::getRelations(edges[i], neighs);
        if( neighs.size() == 1) {
            edges[i]->setVisitBit(1);
            for( int j = 0; j < 2; j++) {
                const JNodePtr &vtx = edges[i]->getNodeAt(j);
                vtx->setVisitBit(1);
            }
        }
    }

    for( size_t i = 0; i < numFaces; i++) {
        faces[i]->setStatus( JMeshEntity::REMOVE);
    }

    for( size_t i = 0; i < numEdges; i++) {
        if( edges[i]->getVisitBit() == 0)
            edges[i]->setStatus( JMeshEntity::REMOVE);
    }

    for( size_t i = 0; i < numNodes; i++) {
        if( nodes[i]->getVisitBit() == 0)
            nodes[i]->setStatus( JMeshEntity::REMOVE);
    }

    pruneAll();
}
///////////////////////////////////////////////////////////////////////////////

void
JMesh::deleteVolumeMesh()
{
    size_t numCells = cells.size();
    if( numCells == 0) return;

    logger->setInfo("deleting volumetric mesh");

    size_t numNodes = nodes.size();
    size_t numEdges = edges.size();
    size_t numFaces = faces.size();

    for( size_t i = 0; i < numNodes; i++) nodes[i]->setVisitBit(0);
    for( size_t i = 0; i < numEdges; i++) edges[i]->setVisitBit(0);
    for( size_t i = 0; i < numFaces; i++) faces[i]->setVisitBit(0);

    JCellSequence neighs;
    for( size_t i = 0; i < numFaces; i++) {
        JFace::getRelations(faces[i], neighs);
        if( neighs.size() == 1) {
            faces[i]->setVisitBit(1);
            for( int j = 0; j < faces[i]->getSize(1); j++) {
                const JEdgePtr &edge = faces[i]->getEdgeAt(j);
                edge->setVisitBit(1);
            }
            for( int j = 0; j < faces[i]->getSize(0); j++) {
                const JNodePtr &vtx = faces[i]->getNodeAt(j);
                vtx->setVisitBit(1);
            }
        }
    }

    for( size_t i = 0; i < numCells; i++)
        cells[i]->setStatus( JMeshEntity::REMOVE);

    for( size_t i = 0; i < numFaces; i++) {
        if( faces[i]->getVisitBit() == 0)
            faces[i]->setStatus( JMeshEntity::REMOVE);
    }

    for( size_t i = 0; i < numEdges; i++) {
        if( edges[i]->getVisitBit() == 0)
            edges[i]->setStatus( JMeshEntity::REMOVE);
    }

    for( size_t i = 0; i < numNodes; i++) {
        if( nodes[i]->getVisitBit() == 0)
            nodes[i]->setStatus( JMeshEntity::REMOVE);
    }
    pruneAll();
}

///////////////////////////////////////////////////////////////////////////////
void JMesh::deleteNodes( int type )
{
    if( !cells.empty() ) {
        logger->setWarn("Cell container not empty: Nodes can not be deleted ");
        return;
    }

    if( !faces.empty() ) {
        logger->setWarn("Face container not empty: Nodes can not be deleted ");
        return;
    }

    if( !edges.empty() ) {
        logger->setWarn("Edge container not empty: Nodes can not be deleted " );
        return;
    }

    switch( type ) {
    case JMeshEntity::ANY_ENTITY:
        for( JNodePtr vtx : nodes)
            vtx->setStatus( JMeshEntity::REMOVE);
        break;
    case JMeshEntity::BOUNDARY_ENTITY:
        logger->setInfo("Freeing mesh boundary nodes");
        for( JNodePtr vtx : nodes)
            if( vtx->isBoundary() ) vtx->setStatus( JMeshEntity::REMOVE);
        break;
    case JMeshEntity::INTERNAL_ENTITY:
        logger->setInfo("Freeing mesh internal nodes");
        for( JNodePtr vtx : nodes)
            if( !vtx->isBoundary() ) vtx->setStatus( JMeshEntity::REMOVE);
        break;
    }
    pruneNodes();
}

///////////////////////////////////////////////////////////////////////////////

void JMesh:: deleteEdges(int type)
{
    if( !cells.empty() ) {
        logger->setWarn("Cell container not empty: edges will not be freed");
        return;
    }

    if( !faces.empty() ) {
        logger->setWarn("Face container not empty: faces can not be freed");
        return;
    }

    switch( type ) {
    case JMeshEntity::ANY_ENTITY:
        logger->setInfo("Deleting all mesh edges");
        for( JEdgePtr edge: edges)
            edge->setStatus( JMeshEntity::REMOVE);
        break;
    case JMeshEntity::BOUNDARY_ENTITY:
        logger->setInfo("Deleting mesh boundary edges");
        for( JEdgePtr edge: edges) {
            if( edge->isBoundary() ) edge->setStatus( JMeshEntity::REMOVE);
        }
        break;
    case JMeshEntity::INTERNAL_ENTITY:
        logger->setInfo("Deleting mesh internal edges");
        for( JEdgePtr edge: edges)
            if( !edge->isBoundary() ) edge->setStatus( JMeshEntity::REMOVE);
        break;
    }
    pruneEdges();
}

///////////////////////////////////////////////////////////////////////////////

void JMesh :: deleteFaceAttribute(const string &s)
{
    size_t nSize = faces.size();
    for( size_t i = 0; i < nSize; i++)
        faces[i]->deleteAttribute(s);
}
///////////////////////////////////////////////////////////////////////////////

void JMesh:: deleteFaces( int type )
{
    if( !cells.empty() ) {
        logger->setWarn("Cell container not empty: faces can not be freed");
        return;
    }

    switch( type ) {
    case JMeshEntity::ANY_ENTITY:
        logger->setInfo("Freeing all the faces in the mesh ");
        for( const JFacePtr &face: faces)
            face->setStatus( JMeshEntity::REMOVE);
        break;
    case JMeshEntity::BOUNDARY_ENTITY:
        logger->setInfo("Freeeing  all the boundary faces in the mesh ");
        for( const JFacePtr &face: faces)
            if(face->isBoundary() ) face->setStatus( JMeshEntity::REMOVE);
        break;
    case JMeshEntity::INTERNAL_ENTITY:
        logger->setInfo("Freeing all the internal faces in the mesh ");
        for( const JFacePtr &face: faces)
            if(!face->isBoundary() ) face->setStatus( JMeshEntity::REMOVE);
        break;
    }
    pruneFaces();
}

///////////////////////////////////////////////////////////////////////////////
void JMesh :: deleteCellAttribute(const string &s)
{
    size_t nSize = cells.size();
    for( size_t i = 0; i < nSize; i++)
        cells[i]->deleteAttribute(s);
}

///////////////////////////////////////////////////////////////////////////////

void JMesh:: deleteCells()
{
    logger->setInfo("Deleting all cells in the container");

    for( const JCellPtr &cell : cells)
        cell->setStatus( JMeshEntity::REMOVE);
    cells.clear();
}

///////////////////////////////////////////////////////////////////////////////

void
JMesh::deleteAll()
{
    logger->setInfo("Deleting all mesh entitites" );
    deleteCells();
    deleteFaces();
    deleteEdges();
    deleteNodes();
}

///////////////////////////////////////////////////////////////////////////////

void
Jaal::set_tfi_coords(int i, int j, int nx, int ny, JNodeSequence &qnodes)
{
    int offset;

    offset = 0;
    const Point3D &v00 = qnodes[offset]->getXYZCoords();

    offset = i;
    const Point3D &vr0 = qnodes[offset]->getXYZCoords();

    offset = (nx - 1);
    const Point3D &v10 = qnodes[offset]->getXYZCoords();

    offset = j*nx;
    const Point3D &v0s = qnodes[offset]->getXYZCoords();

    offset = j * nx + (nx - 1);
    const Point3D &v1s = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx;
    const Point3D &v01 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + i;
    const Point3D &vr1 = qnodes[offset]->getXYZCoords();

    offset = (ny - 1) * nx + (nx - 1);
    const Point3D &v11 = qnodes[offset]->getXYZCoords();

    Point3D vrs;

    double dr = 2.0 / (double) (nx - 1);
    double ds = 2.0 / (double) (ny - 1);

    double r = -1.0 + i*dr;
    double s = -1.0 + j*ds;
    for (int k = 0; k < 3; k++) {
        vrs[k] = TFI::transfinite_blend(r, s,
                                        v00[k], v10[k], v11[k], v01[k],
                                        vr0[k], v1s[k], vr1[k], v0s[k]);
    }
    offset = j * nx + i;
    qnodes[offset]->setXYZCoords(vrs);
}

///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
JFaceSequence
JMesh::filteredFaces(int facetype) const
{
    JFaceSequence::const_iterator it;
    size_t ncount = 0;
    for (it = faces.begin(); it != faces.end(); ++it) {
        Face *face = *it;
        if (face->getType() == facetype)
            ncount++;
    }

    JFaceSequence tmpfaces;
    if (ncount) {
        tmpfaces.resize(ncount);
        size_t index = 0;
        for (it = faces.begin(); it != faces.end(); ++it) {
            Face *face = *it;
            if (face->getType() == facetype)
                tmpfaces[index++] = face;
        }
    }

    return tmpfaces;
}

///////////////////////////////////////////////////////////////////////////////

Mesh *
Jaal::struct_quad_grid(int nx, int ny)
{
    Mesh *quadmesh = new Mesh;

    double dx = 2.0 / (nx - 1);
    double dy = 2.0 / (ny - 1);

    Point3D xyz;

    int index = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            xyz[0] = -1.0 + i * dx;
            xyz[1] = -1.0 + j * dy;
            xyz[2] = 0.0;
            Vertex *vnew = Vertex::newObject();
            vnew->setID(index++);
            vnew->setXYZCoords(xyz);
            quadmesh->addObject(vnew);
        }
    }

    JNodeSequence nodes(4);
    index = 0;
    Face *newquad;
    for (int j = 0; j < ny - 1; j++) {
        for (int i = 0; i < nx - 1; i++) {
            int n0 = j * nx + i;
            int n1 = n0 + 1;
            int n2 = n1 + nx;
            int n3 = n0 + nx;
            nodes[0] = quadmesh->getNodeAt(n0);
            nodes[1] = quadmesh->getNodeAt(n1);
            nodes[2] = quadmesh->getNodeAt(n2);
            nodes[3] = quadmesh->getNodeAt(n3);
            newquad = Face::newObject();
            newquad->setNodes(nodes);
            quadmesh->addObject(newquad);
        }
    }
    return quadmesh;
}

////////////////////////////////////////////////////////////////////

bool
Face::has_all_bound_nodes() const
{
    for (int i = 0; i < getSize(0); i++) {
        Vertex *v = getNodeAt(i);
        if (!v->isBoundary()) return 0;
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

int Mesh::get_breadth_first_ordered_nodes(JNodeSequence &seq, Vertex *vstart, MeshFilter *filter)
{
    assert(vstart != nullptr);

    seq.clear();

    int relexist0 = buildRelations(0, 0);

    size_t numnodes = getSize(0);

    if (numnodes == 0) return 1;

    for (size_t i = 0; i < numnodes; i++) {
        Vertex *v = getNodeAt(i);
        v->setVisitBit(0);
        v->setAttribute("Layer", 0);
    }

    if (vstart == 0) vstart = getNodeAt(0);

    seq.reserve(numnodes);

    list<Vertex*> vertexQ;
    vertexQ.push_back(vstart);
    JNodeSequence vneighs;

    int proceed = 1;
    int currlevel;
    while (!vertexQ.empty()) {
        Vertex *curr_vertex = vertexQ.front();
        vertexQ.pop_front();
        curr_vertex->getAttribute("Layer", currlevel);
        if (filter) {
            if (curr_vertex != vstart) proceed = filter->pass(curr_vertex);
        }
        if (!curr_vertex->isVisited()) {
            seq.push_back(curr_vertex);
            if (!proceed) break;
            curr_vertex->setVisitBit(1);
            curr_vertex->getRelations( vneighs );
            for (size_t i = 0; i < vneighs.size(); i++) {
                if (!vneighs[i]->isVisited()) {
                    vertexQ.push_back(vneighs[i]);
                    vneighs[i]->setAttribute("Layer", currlevel + 1);
                }
            }
        }
    }

    if (!relexist0)
        clearRelations(0, 0);

    // Free unused memory in sequence...
    if (!seq.empty()) JNodeSequence(seq).swap(seq);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int Mesh::get_depth_first_ordered_nodes(JNodeSequence &seq, Vertex *vstart, MeshFilter *filter)
{
    int relexist0 = buildRelations(0, 0);

    size_t numnodes = getSize(0);
    list<Vertex*> vertexQ;
    for (size_t i = 0; i < numnodes; i++) {
        Vertex *v = getNodeAt(i);
        v->setVisitBit(0);
    }

    seq.clear();

    if (vstart == 0) vstart = getNodeAt(0);
    vertexQ.push_back(vstart);
    JNodeSequence vneighs;

    while (!vertexQ.empty()) {
        Vertex *curr_vertex = vertexQ.front();
        vertexQ.pop_front();
        if (!curr_vertex->isVisited()) {
            seq.push_back(curr_vertex);
            curr_vertex->setVisitBit(1);
            curr_vertex->getRelations( vneighs );
            for (size_t i = 0; i < vneighs.size(); i++) {
                if (!vneighs[i]->isVisited())
                    vertexQ.push_front(vneighs[i]);
            }
        }
    }

    if (!relexist0)
        clearRelations(0, 0);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int Jaal::SurfPatch::getPosOf(const Vertex *v)
{
    for (size_t i = 0; i < bound_nodes.size(); i++)
        if (bound_nodes[i] == v) return i;

    cout << "Error: Vertex not found on the boundary " << endl;
    exit(0);

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

JNodeSequence SurfPatch::get_bound_nodes(const Vertex *src, const Vertex *dst)
{
    int start_pos = getPosOf(src);
    int end_pos = getPosOf(dst);
    int nsize = bound_nodes.size();

    if (end_pos == 0) end_pos = nsize;
    assert(end_pos > start_pos);

    JNodeSequence seq(end_pos - start_pos + 1);
    int index = 0;
    for (int i = start_pos; i <= end_pos; i++)
        seq[index++] = bound_nodes[i % nsize];

    return seq;
}

////////////////////////////////////////////////////////////////////////////////

int Jaal::SurfPatch::search_boundary()
{
    corners.clear();
    boundary.clear();
    Vertex *vertex;

    assert(!faces.empty());

    // We need to rebuild relations locally to identfy corners and boundary.
    JFaceSet::const_iterator fiter;
    std::map<Vertex*, JFaceSet> relations02;

    for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
        Face *face = *fiter;
        for (int j = 0; j < face->getSize(0); j++) {
            vertex = face->getNodeAt(j);
            relations02[vertex].insert(face);
        }
    }

    // A boundary edge must have exactly one face neighbor...
    JFaceSequence faceneighs;
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
        Face *face = *fiter;
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++) {
            Vertex *v0 = face->getNodeAt( j + 0 );
            Vertex *v1 = face->getNodeAt( j + 1 );
            faceneighs.clear();
            assert(relations02[v0].size() > 0);
            assert(relations02[v1].size() > 0);
            boost::set_intersection(relations02[v0], relations02[v1], back_inserter(faceneighs));
            if (faceneighs.size() == 1) {
                Edge newedge(v0, v1);
                boundary.push_back(newedge);
            }
        }
    }

    // Sequence the chain and start from one of the corner...
    int err = Mesh::make_chain(boundary);
    if (err) return 2;

    //
    // Identify corners in the mesh.
    // Should we check only one vertex per edge ?
    //
    JFaceSet neighs;
    int boundSize = boundary.size();
    for (int k = 0; k < boundSize; k++) {
        vertex = boundary[k].getNodeAt(0);
        neighs = relations02[vertex];
        if (neighs.size() == 1) corners.insert(vertex);

        vertex = boundary[k].getNodeAt(1);
        neighs = relations02[vertex];
        if (neighs.size() == 1) corners.insert(vertex);
    }

    // Start the chain from one of the corners.
    err = Mesh::rotate_chain(boundary, *corners.begin());
    if (err) return 3;

    // Collect the sequence of boundary nodes...,
    bound_nodes.resize( boundSize );
    for (int k = 0; k < boundSize; k++) {
        bound_nodes[k] = boundary[k].getNodeAt(0); // Only the first node.
    }

    //
    // Collect the inner nodes of the patch. These nodes will be deleted, if
    // the remesh operation is successful...
    //
    inner_nodes.clear();
    for (fiter = faces.begin(); fiter != faces.end(); ++fiter) {
        Face *face = *fiter;
        int nnodes = face->getSize(0);
        for (int j = 0; j < nnodes; j++) {
            Vertex *v = face->getNodeAt(j);
            if (find(bound_nodes.begin(), bound_nodes.end(), v) == bound_nodes.end())
                inner_nodes.insert(v);
        }
    }

    // Split the boundary loop into segments.
    // (i.e. End of the segments are the corners identified earlier )
    set_boundary_segments();

    return 0;
}

////////////////////////////////////////////////////////////////////

void Jaal::SurfPatch::set_boundary_segments()
{
    // Although this stage will not come in this algorithm...
    if (corners.size() == 0) return;

    cornerPos.resize(corners.size() + 1);

    JNodeSet::const_iterator it;
    int index = 0;
    for (it = corners.begin(); it != corners.end(); ++it) {
        cornerPos[index++] = getPosOf(*it);
    }

    cornerPos[corners.size()] = bound_nodes.size();

    boost::sort(cornerPos);

    segSize.resize(corners.size());

    for (size_t i = 0; i < corners.size(); i++)
        segSize[i] = cornerPos[(i + 1)] - cornerPos[i] + 1;
}

////////////////////////////////////////////////////////////////////

int Jaal::SurfPatch::reorient_4_sided_loop()
{
    //Always remeshable, nothing has to be done...
    if ((segSize[0] == segSize[2]) && (segSize[1] == segSize[3])) return 0;

    //////////////////////////////////////////////////////////////////////////
    // Defination:  A four sided convex loop has four segments.
    // Objectives:  A four sided convex loop must be orietned such that
    //   1.  First segment must be smaller than the third one, because
    //       we need to create a triangle patch based at segment#3.
    //
    //   2.  If there are two choices, then the side having irregular
    //       node must be given higher priority. ( Here irregulaty means
    //       that vertex valency < 4 ).
    //
    //   3.  A side having less number of nodes on the first segment than
    //       the third is given preference.
    //
    // Pre-Conditions  :  A loop must be oriented ( CW or CCW ).
    //
    // Date: 17th Nov. 2010.
    //////////////////////////////////////////////////////////////////////////

    Vertex *start_corner = nullptr;

    if (segSize[0] == segSize[2]) {
        if (min(segSize[1], segSize[3]) == 2) return 1;
        //  Either Segment 2 or 3 must be starting node
        if (segSize[1] < segSize[3])
            start_corner = bound_nodes[ cornerPos[1] ];
        else
            start_corner = bound_nodes[ cornerPos[3] ];
        start_boundary_loop_from(start_corner);
    }

    if (min(segSize[0], segSize[2]) == 2) return 1;

    if (segSize[2] < segSize[0]) {
        start_corner = bound_nodes[ cornerPos[2] ];
        start_boundary_loop_from(start_corner);
    }

    // By this stage, the loop must be reoriented correctly.
    assert(segSize[0] < segSize[2]);

    cout << " Careful to change this code " << endl;
    // Great, we found one irregular node on the first boundary...
    //   if( has_irregular_node_on_first_segment() ) return 1;

    // If the segment 2 and 3 have same size, Alas, nothing can be done.
    if (segSize[1] == segSize[3]) return 1;

    if (min(segSize[1], segSize[3]) == 2) return 1;

    if (segSize[3] < segSize[1]) {
        start_corner = bound_nodes[ cornerPos[3] ];
        start_boundary_loop_from(start_corner);
    } else {
        start_corner = bound_nodes[ cornerPos[1] ];
        start_boundary_loop_from(start_corner);
    }

    //
    // Note that we didn't check for irregular node here. So if this segment
    // has at least one irregular node, then we are lucky. Otherwise decision
    // to remesh it done based wthere remeshing will result in the reduction
    // of irregular nodes in patch.
    //
    assert(segSize[0] < segSize[2]);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void Jaal::SurfPatch::start_boundary_loop_from(Vertex *vmid)
{
    assert(corners.find(vmid) != corners.end());

    JNodeSequence::iterator middle;
    middle = find(bound_nodes.begin(), bound_nodes.end(), vmid);
    assert(middle != bound_nodes.end());

    std::rotate(bound_nodes.begin(), middle, bound_nodes.end());
    assert(bound_nodes[0] == vmid);

    set_boundary_segments();
}
///////////////////////////////////////////////////////////////////////////////

void Jaal::project_boundary_on_sphere( Mesh *mesh, double radius)
{
    size_t numNodes = mesh->getSize(0);
    Point3D p3d, pC;

    pC[0] = 0.0;
    pC[1] = 0.0;
    pC[2] = 0.0;
    size_t nCount = 0;
    for( size_t i = 0; i < numNodes; i++) {
        Vertex *vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            p3d = vtx->getXYZCoords();
            pC[0] +=  p3d[0];
            pC[1] +=  p3d[1];
            pC[2] +=  p3d[2];
            nCount++;
        }
    }
    pC[0] /= ( double) nCount;
    pC[1] /= ( double) nCount;
    pC[2] /= ( double) nCount;

    double dl;
    for( size_t i = 0; i < numNodes; i++) {
        Vertex *vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && vtx->isBoundary() ) {
            p3d = vtx->getXYZCoords();
            p3d[0] =  p3d[0] - pC[0];
            p3d[1] =  p3d[1] - pC[1];
            p3d[2] =  p3d[2] - pC[2];
            dl = sqrt(p3d[0]*p3d[0] + p3d[1]*p3d[1] + p3d[2]*p3d[2] );
            p3d[0] = pC[0] + radius*p3d[0]/dl;
            p3d[1] = pC[1] + radius*p3d[1]/dl;
            p3d[2] = pC[2] + radius*p3d[2]/dl;
            vtx->setXYZCoords(p3d);
        }
    }
}
void Mesh :: objects_from_pool( size_t n, JNodeSequence  &objects)
{
    objects.clear();

    if( n == 0) return;

    objects.reserve( n );

    size_t ncount = 0;

    /*
         while( !garbageNodes.empty() ) {
              Vertex *v = garbageNodes.front();
              garbageNodes.pop_front();
              if( v->isRemoved() ) {
                   v->reset();
                   v->setStatus( MeshEntity::INACTIVE);
                   objects.push_back(v);
                   ncount++;
                   if( ncount == n) break;
              }
         }
    */

    for( size_t i = ncount; i < n; i++) {
        Vertex *v = Vertex::newObject();
        objects.push_back(v);
        addObject( v );
    }

    assert( objects.size() == n );
}

///////////////////////////////////////////////////////////////////////////////

void Mesh :: objects_from_pool( size_t n, JFaceSequence &objects, int type)
{
    objects.clear();

    if( n == 0) return;

    objects.reserve( n );

    size_t ncount = 0;
    /*
         while( !garbageFaces.empty() ) {
              Face *f = garbageFaces.front();
              garbageFaces.pop_front();
              if( f->isRemoved() ) {
                   f->reset();
                   f->setStatus( MeshEntity::INACTIVE);
                   objects.push_back(f);
                   ncount++;
                   if( ncount == n) break;
              }
         }
    */

    for( size_t i = ncount; i < n; i++) {
        Face *f = Face::getProduct(type);
        objects.push_back(f);
        addObject(f);
    }

    assert( objects.size() == n );
}

////////////////////////////////////////////////////////////////////////////////
bool JMesh :: hasTopologyChanged() const
{
    if( topoChanged[1] || topoChanged[2] || topoChanged[3] ) {
        topology->setBoundaryKnown(0);
        return 1;
    }
    return 0;
}


#endif

