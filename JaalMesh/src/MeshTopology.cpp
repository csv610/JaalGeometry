#include <iomanip>

#include "Mesh.hpp"
#include "MeshTopology.hpp"
#include "basic_math.hpp"
#include "circumcenter.hpp"
#include "MeshMatrix.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>

using namespace std;
using namespace Jaal;
using namespace boost;

///////////////////////////////////////////////////////////////////////////////

JLogger* JMeshTopology::logger = JLogger::getInstance();

///////////////////////////////////////////////////////////////////////////////
bool
JMeshTopology ::isSimple()
{
    return 1;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshTopology :: setRCMOrdering()
{
    if( mesh == nullptr) return;

    typedef adjacency_list<vecS, vecS, undirectedS,
            property<vertex_color_t, default_color_type,
            property<vertex_degree_t,int> > > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::vertices_size_type size_type;

    int numNodes = mesh->getSize(0);
    int numEdges = mesh->getSize(1);

    Graph G(numNodes);
    for (size_t i = 0; i < numEdges; ++i) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        size_t v0 = edge->getNodeAt(0)->getID();
        size_t v1 = edge->getNodeAt(1)->getID();
        add_edge( v0, v1, G);
    }

    graph_traits<Graph>::vertex_iterator ui, ui_end;

    property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
        deg[*ui] = degree(*ui, G);

    property_map<Graph, vertex_index_t>::type
    index_map = get(vertex_index, G);

    std::vector<Vertex> inv_perm(num_vertices(G));
    std::vector<size_type> perm(num_vertices(G));

    cuthill_mckee_ordering(G, inv_perm.rbegin(), get(vertex_color, G),
                           make_degree_map(G));

    for (size_type c = 0; c != inv_perm.size(); ++c)
        perm[index_map[inv_perm[c]]] = c;

    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        size_t id = vtx->getID();
        vtx->setID( perm[id] );
    }

    JEntityCompare  *jcomp = new JEntityIDCompare;
    mesh->sort(jcomp);
    delete jcomp;
}

///////////////////////////////////////////////////////////////////////////////

vector<int> JMeshTopology::getMinNodesColor() const
{
    typedef adjacency_list<listS, vecS, undirectedS> Graph;
    typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef property_map<Graph, vertex_index_t>::const_type vertex_index_map;

    mesh->pruneNodes();
    mesh->pruneEdges();

    size_t numnodes = mesh->getSize(0);
    size_t numedges = mesh->getSize(1);

    Graph graph(numnodes);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            int v0 = edge->getNodeAt(0)->getID();
            int v1 = edge->getNodeAt(1)->getID();
            add_edge(v0,v1,graph);
        }
    }

    std::vector<vertices_size_type> color_vec(num_vertices(graph));
    iterator_property_map<vertices_size_type*, vertex_index_map>
    color(&color_vec.front(), get(vertex_index, graph));
    vertices_size_type num_colors = sequential_vertex_coloring(graph, color);

    vector<int> nodeColor( numnodes);
    graph_traits<Graph>::vertex_iterator vi, vi_end;
    for(tie(vi,vi_end) = vertices(graph); vi != vi_end; ++vi) {
        nodeColor[*vi] = color[*vi];
    }
    return nodeColor;
}

///////////////////////////////////////////////////////////////////////////////

JNodeSequence JMeshTopology :: getNodeSubComplex()
{
    JNodeSequence nodes;

    if( mesh == nullptr) return nodes;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->setVisitBit(0);
    }

    int topDim = getDimension();
    int numedges, numfaces, numcells;

    switch(topDim )
    {
    case 1:
        numedges = mesh->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            for( int j = 0; j > 2; j++) {
                const JNodePtr &vtx  = edge->getNodeAt(j);
                vtx->setVisitBit(1);
            }
        }
        break;
    case 2:
        numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            for( int j = 0; j < face->getSize(0); j++) {
                const JNodePtr &vtx = face->getNodeAt(j);
                vtx->setVisitBit(1);
            }
        }
        break;
    case 3:
        numcells = mesh->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            for( int j = 0; j < cell->getSize(0); j++) {
                const JNodePtr &vtx = cell->getNodeAt(j);
                vtx->setVisitBit(1);
            }
        }
        break;
    }

    for( size_t i = 0; i < numnodes; i++)  {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->getVisitBit() == 1)
            nodes.push_back(vtx);
    }
    return nodes;
}

///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JMeshTopology :: getEdgeSubComplex()
{
    JEdgeSequence edges;

    if( mesh == nullptr) return edges;

    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        edge->setVisitBit(0);
    }

    int topDim = getDimension();
    int  numfaces, numcells;

    switch(topDim )
    {
    case 2:
        numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            for( int j = 0; j < face->getSize(1); j++) {
                const JEdgePtr &edge = face->getEdgeAt(j);
                edge->setVisitBit(1);
            }
        }
        break;
    case 3:
        numcells = mesh->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            for( int j = 0; j < cell->getSize(1); j++) {
                const JEdgePtr &edge = cell->getEdgeAt(j);
                edge->setVisitBit(1);
            }
        }
        break;
    }

    for( size_t i = 0; i < numedges; i++)  {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->getVisitBit() == 1)
            edges.push_back(edge);
    }
    return edges;
}

///////////////////////////////////////////////////////////////////////////////

JFaceSequence JMeshTopology :: getFaceSubComplex()
{
    JFaceSequence faces;

    if( mesh == nullptr) return faces;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->setVisitBit(0);
    }


    size_t numcells = mesh->getSize(3);
    for( size_t i = 0; i < numcells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        for( int j = 0; j < cell->getSize(1); j++) {
            const JFacePtr &face = cell->getFaceAt(j);
            face->setVisitBit(1);
        }
    }

    for( size_t i = 0; i < numfaces; i++)  {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->getVisitBit() == 1)
            faces.push_back(face);
    }
    return faces;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getOrphaned( JNodeSequence &orphans)
{
    if( mesh == nullptr) return;

    orphans.clear();
    JNodeSet aset, bset;
    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) aset.insert(vtx);
    }

    size_t numCells = mesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < cell->getSize(0); j++)
                bset.insert(cell->getNodeAt(j));
        }
    }

    size_t numFaces = mesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            for( int j = 0; j < face->getSize(0); j++)
                bset.insert(face->getNodeAt(j));
        }
    }

    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            for( int j = 0; j < 2; j++)
                bset.insert(edge->getNodeAt(j));
        }
    }
    assert( aset.size() >= bset.size());

    boost::set_difference(aset, bset, std::inserter(orphans, orphans.begin()));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getOrphaned( JEdgeSequence &orphans)
{
    if( mesh == nullptr) return;

    orphans.clear();
    JEdgeSet aset, bset;
    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) aset.insert(edge);
    }

    size_t numCells = mesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < cell->getSize(1); j++)
                bset.insert(cell->getEdgeAt(j));
        }
    }

    size_t numFaces = mesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            for( int j = 0; j < face->getSize(1); j++)
                bset.insert(face->getEdgeAt(j));
        }
    }

    assert( aset.size() >= bset.size());

    boost::set_difference(aset, bset, std::inserter(orphans, orphans.begin()));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getOrphaned( JFaceSequence &orphans)
{
    if( mesh == nullptr) return;

    orphans.clear();
    JFaceSet aset, bset;
    size_t numFaces = mesh->getSize(3);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) aset.insert(face);
    }

    size_t numCells = mesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            for( int j = 0; j < cell->getSize(2); j++)
                bset.insert(cell->getFaceAt(j));
        }
    }

    assert( aset.size() >= bset.size());

    boost::set_difference(aset, bset, std::inserter(orphans, orphans.begin()));
}

///////////////////////////////////////////////////////////////////////////////
JFacePtr JMeshTopology :: getCenter( const JFaceSequence &faces)
{
    JEdgeSequence boundary;
    getSubmeshBoundary(  faces, boundary);
    if( boundary.empty() ) return nullptr;

    JFaceSet faceSet;

    int levelid = -1;
    for( const JFacePtr &f : faces) {
        f->setAttribute("Level", levelid);
        faceSet.insert(f);
        f->setVisitBit(0);
    }

    JFaceSequence faceneighs;
    deque<JFacePtr> faceQ;
    for( const JEdgePtr &edge : boundary) {
        JEdge::getRelations( edge, faceneighs);
        for( const JFacePtr &f : faceneighs) {
            if( faceSet.find(f) != faceSet.end() )
                faceQ.push_back(f);
        }
    }

    levelid = 0;
    for( const JFacePtr &f : faceQ)
        f->setAttribute("Level", levelid);

    while(!faceQ.empty() ) {
        JFacePtr currface = faceQ.front();
        faceQ.pop_front();
        if( currface->getVisitBit() == 0) {
            currface->getAttribute("Level", levelid);
            int nextid = levelid + 1;
            JFace::getRelations12(currface, faceneighs);
            for( const JFacePtr &f: faceneighs) {
                if( f->getVisitBit() == 0) {
                    f->setAttribute("Level", nextid);
                    faceQ.push_back(f);
                }
            }
        }
        currface->setVisitBit(1);
    }

    int maxid = -1;
    JFacePtr center;
    for( const JFacePtr &face: faces) {
        assert( face->getVisitBit() == 1);
        face->getAttribute("Level", levelid);
        if( levelid > maxid ) {
            center = face;
            maxid   = levelid;
        }
    }
    return center;
}
///////////////////////////////////////////////////////////////////////////////
size_t JMeshTopology :: countElementType(int etype) const
{
    size_t nCount = 0;

    size_t numfaces = mesh->getSize(2);
    size_t numcells = mesh->getSize(3);

    if( etype == JFace::TRIANGLE) {
        #pragma omp parallel for reduction(+:nCount)
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() )
                if( f->getSize(0) == 3) nCount++;
        }
        return nCount;
    }

    if( etype == JFace::QUADRILATERAL) {
        #pragma omp parallel for reduction(+:nCount)
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() )
                if( f->getSize(0) == 4) nCount++;
        }
        return nCount;
    }

    if( etype == JFace::POLYGON) {
        #pragma omp parallel for reduction(+:nCount)
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() && f->getSize(0) > 4) nCount++;
        }
        return nCount;
    }

    if( etype == JCell::TETRAHEDRON) {
        #pragma omp parallel for reduction(+:nCount)
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            if( c->isActive() && c->getSize(0) == 4) nCount++;
        }
        return nCount;
    }

    if( etype == JCell::HEXAHEDRON) {
        #pragma omp parallel for reduction(+:nCount)
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            if( c->isActive() && c->getSize(0) == 8) nCount++;
        }
        return nCount;
    }

    if( etype == JCell::POLYHEDRON) {
        #pragma omp parallel for reduction(+:nCount)
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            if( c->isActive() )
                if(c->getSize(0) !=4  && c->getSize(0) !=8) nCount++;
        }
        return nCount;
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

bool JMeshTopology :: isSameAs( const JMeshPtr &checkmesh) const
{
    if( mesh == nullptr || checkmesh == nullptr ) return 0;

    if( mesh->getSize(3) != checkmesh->getSize(3) ) {
        cout << "Info: number of cells are same in two mesh objects " << endl;
        return 0;
    }

    if( mesh->getSize(2) != checkmesh->getSize(2) ) {
        cout << "Info: number of cells are same in two mesh objects " << endl;
        return 0;
    }

    if( mesh->getSize(1) != checkmesh->getSize(1) ) {
        cout << "Info: number of cells are same in two mesh objects " << endl;
        return 0;
    }

    if( mesh->getSize(0) != checkmesh->getSize(0) ) {
        cout << "Info: number of vertices are same in two mesh objects " << endl;
        return 0;
    }

    size_t numCells = mesh->getSize(3);

    if( numCells ) {
        for( size_t i = 0; i < numCells; i++) {
            const JCellPtr &acell = mesh->getCellAt(i);
            const JCellPtr &bcell = mesh->getCellAt(i);
            if( !acell->hasSameNodesID(bcell) ) {
                cout << "Info: Two cell objects are not same" << endl;
                return 0;
            }
        }
        return 1;
    }

    size_t numFaces = mesh->getSize(2);
    if( numFaces ) {
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &aface = mesh->getFaceAt(i);
            const JFacePtr &bface = mesh->getFaceAt(i);
            if( !aface->hasSameNodesID(bface) ) {
                cout << "Info: Two face objects are not same " << endl;
                cout << "FaceID " << i << " Nodes ";
                for( int j = 0; j < aface->getSize(0); j++)
                    cout << aface->getNodeAt(j)->getID() << " ";
                cout << endl;
                cout << "FaceID " << i << " Nodes ";
                for( int j = 0; j < aface->getSize(0); j++)
                    cout << bface->getNodeAt(j)->getID() << " ";
                cout << endl;
                return 0;
            }
        }
        return 1;
    }

    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &aedge = mesh->getEdgeAt(i);
        const JEdgePtr &bedge = mesh->getEdgeAt(i);
        if( !aedge->hasSameNodesID(bedge) ) {
            cout << "Info: Two edge objects are not same " << endl;
            return 0;
        }
    }

    return 1;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getRim( const JNodePtr &vtx, JEdgeSequence &edges)
{
    edges.clear();
    JFaceSequence faceneighs;

    JNode::getRelations(vtx, faceneighs);
    if( faceneighs.empty() ) {
        cout << "Warning: Rim can not be identified: vertex-face relations missing" << endl;
        return;
    }

    JEdgeSet edgeSet;
    for( const JFacePtr &face : faceneighs) {
        int numedges = face->getSize(1);
        for( int i = 0; i < numedges; i++) {
            const JEdgePtr &edge  = face->getEdgeAt(i);
            const JNodePtr &v0 = edge->getNodeAt(0);
            const JNodePtr &v1 = edge->getNodeAt(1);
            if( (v0 != vtx) && (v1 != vtx) )
                edgeSet.insert(edge);
        }
    }
    if( edgeSet.empty() ) return;

    edges.resize( edgeSet.size() );
    boost::copy(edgeSet, edges.begin() );
}

////////////////////////////////////////////////////////////////////////////////

bool JMeshTopology :: isDisk()
{
    int eular  = getEulerCharacteristic();
    if(eular == 1) return 1;
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JMeshTopology :: getLexicographic( JEdgeSequence &seq)
{
    seq.clear();

    if( mesh == nullptr ) return 1;

    seq = mesh->getEdges();
    size_t nSize =  seq.size();
    for( size_t i = 0; i < nSize; i++) {
        int v0 = seq[i]->getNodeAt(0)->getID();
        int v1 = seq[i]->getNodeAt(1)->getID();
        if( v0 > v1 ) seq[i]->reverse();
    }

    std::sort( seq.begin(), seq.end(), JEdgeTopology::lexiCompare);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
int JMeshTopology :: getLexicographic( JFaceSequence &seq)
{
    seq.clear();
    if( mesh == nullptr ) return 1;
    seq = mesh->getFaces();

    std::sort( seq.begin(), seq.end(),
               []( const JFacePtr &a, const JFacePtr &b)
    {
        size_t amin = a->getNodeAt(0)->getID();
        for( int i = 1; i < a->getSize(0); i++)
            amin  = min(amin, a->getNodeAt(i)->getID());

        size_t bmin = b->getNodeAt(0)->getID();
        for( int i = 1; i < b->getSize(0); i++)
            bmin  = min(bmin, b->getNodeAt(i)->getID());
        return amin < bmin;
    }
             );
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
void JMeshTopology :: getEntitySet( const JEdgeSequence &edges, JNodeSequence &nodes)
{
    nodes.clear();
    JNodeSet vset;
    for( const JEdgePtr &edge : edges) {
        vset.insert( edge->getNodeAt(0) );
        vset.insert( edge->getNodeAt(1) );
    }
    if( vset.empty() ) return;
    nodes.resize( vset.size() );
    boost::copy(vset, nodes.begin() );
}
////////////////////////////////////////////////////////////////////////////////
void JMeshTopology :: getEntitySet( const JFaceSequence &faces, JNodeSequence &nodes)
{
    nodes.clear();
    JNodeSet vset;
    for( const JFacePtr &face : faces) {
        for( int i = 0; i < face->getSize(0); i++)
            vset.insert( face->getNodeAt(i) );
    }
    if( vset.empty() ) return;
    nodes.resize( vset.size() );
    boost::copy( vset, nodes.begin() );
}
////////////////////////////////////////////////////////////////////////////////
void JMeshTopology :: getEntitySet( const JCellSequence &cells, JNodeSequence &nodes)
{
    nodes.clear();
    JNodeSet vset;
    for( const JCellPtr &cell : cells) {
        for( int i = 0; i < cell->getSize(0); i++)
            vset.insert( cell->getNodeAt(i) );
    }
    if( vset.empty() ) return;
    nodes.resize( vset.size() );
    boost::copy( vset, nodes.begin() );
}
////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getEntitySet( const JFaceSequence &faces, JEdgeSequence &edges)
{
    edges.clear();
    JEdgeSet eset;
    for( const JFacePtr &face : faces) {
        for( int i = 0; i < face->getSize(1); i++)
            eset.insert( face->getEdgeAt(i) );
    }
    if( eset.empty() ) return;
    edges.resize( eset.size() );
    boost::copy( eset, edges.begin() );
}
////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getEntitySet( const JCellSequence &cells, JEdgeSequence &edges)
{
    edges.clear();
    JEdgeSet eset;
    for( const JCellPtr &cell : cells) {
        for( int i = 0; i < cell->getSize(1); i++)
            eset.insert( cell->getEdgeAt(i) );
    }
    if( eset.empty() ) return;
    edges.resize( eset.size() );
    boost::copy( eset, edges.begin() );
}
////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getEntitySet( const JCellSequence &cells, JFaceSequence &faces)
{
    faces.clear();
    JFaceSet fset;
    for( const JCellPtr &cell : cells) {
        for( int i = 0; i < cell->getSize(2); i++)
            fset.insert( cell->getFaceAt(i) );
    }
    if( fset.empty() ) return;
    faces.resize( fset.size() );
    boost::copy(fset, faces.begin() );
}
////////////////////////////////////////////////////////////////////////////////


void JMeshTopology :: reverseAll()
{
    if( mesh == nullptr ) return;

    logger->setInfo("Reversing the face-nodes order ");

    size_t numfaces = mesh->getSize(2);
    #pragma omp parallel for
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        f->reverse();
    }
}

////////////////////////////////////////////////////////////////////////////////

void
JMeshTopology::getNodesArray( vector<size_t> &nodearray, vector<int> &elmtopo)
{
    if( mesh == nullptr ) return;

    logger->setInfo("Collecting connectivity array");

    elmtopo.clear();
    nodearray.clear();
    int  topo, tDim  = getDimension();
    size_t numnodes = mesh->getSize(0);

    map<size_t,size_t> nodemap;
    size_t index = 0;
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) nodemap[v->getID()] = index++;
    }

    if( tDim == 3 ) {
        size_t numcells = mesh->getSize(3);
        topo = mesh->getTopology()->getElementsType( 3 );
        if (topo) nodearray.reserve(topo * numcells);

        elmtopo.reserve( numcells );

        for (size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                elmtopo.push_back( cell->getSize(0) );
                for (int j = 0; j < cell->getSize(0); j++) {
                    const JNodePtr &v = cell->getNodeAt(j);
                    int lid   = nodemap[v->getID()];
                    nodearray.push_back(lid);
                }
            }
        }
    }

    if( tDim == 2 ) {
        size_t numfaces = mesh->getSize(2);
        topo = mesh->getTopology()->getElementsType( 2 );
        if (topo) nodearray.reserve(topo * numfaces);

        elmtopo.reserve( numfaces );

        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                elmtopo.push_back( face->getSize(0) );
                for (int j = 0; j < face->getSize(0); j++) {
                    const JNodePtr &v = face->getNodeAt(j);
                    int lid   = nodemap[v->getID()];
                    nodearray.push_back(lid);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

size_t JMeshTopology::getNumIrregularNodes()
{
    if( mesh == nullptr ) return 0;

    logger->setInfo("Getting irregular nodes in the mesh");

    int degree_of_regular_node = 0;
    int topdim = getDimension();

    int elem_type;
    switch(topdim ) {
    case 2:
        elem_type = getElementsType(2);
        if( elem_type == JFace::TRIANGLE) degree_of_regular_node = 6;
        if( elem_type == JFace::QUADRILATERAL) degree_of_regular_node = 4;
        break;
    case 3:
        elem_type = getElementsType(3);
        if( elem_type == JCell::HEXAHEDRON) degree_of_regular_node = 6;
        break;
    }

    if( degree_of_regular_node == 0) {
        cout << "Warning: irregular nodes not counted " << endl;
        return 0;
    }

    mesh->buildRelations(0, topdim);
    searchBoundary();

    size_t numnodes = mesh->getSize(0);
    size_t ncount = 0;

    #pragma omp parallel for reduction(+:ncount)
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            int nSize = vertex->getNumRelations(topdim);
            if (!vertex->isBoundary()) {
                if ( nSize != degree_of_regular_node)
                    ncount++;
            }
        }
    }

    return ncount;
}

////////////////////////////////////////////////////////////////////////////////

int JMeshTopology :: collectEdges()
{
    if( mesh == nullptr ) return 1;

    size_t numcells = mesh->cells.size();
    size_t numfaces = mesh->faces.size();

    if( numfaces + numcells == 0) return 2;

    logger->setInfo("Collecting edge simplices");

    JEdgeSet  eset;
    JEdgeSequence faceedges, celledges;

    for (size_t icell = 0; icell < numcells; icell++) {
        const JCellPtr &cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            celledges = cell->getEdges();
            eset.insert( celledges.begin(), celledges.end() );
        }
    }

    if( numcells == 0) {
        for (size_t iface = 0; iface < numfaces; iface++) {
            const JFacePtr &face = mesh->getFaceAt(iface);
            if( face->isActive() ) {
                faceedges = face->getEdges();
                eset.insert( faceedges.begin(), faceedges.end() );
            }
        }
    }

    mesh->edges.clear();
    size_t numEdges = eset.size();
    if( numEdges) {
        mesh->edges.resize(numEdges);
        std::copy( eset.begin(), eset.end(), mesh->edges.begin() );
    }

    mesh->enumerate(1);

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getDerivedRelations(const JNodePtr &vtx, JEdgeSequence &edges)
{
    edges.clear();
    JFaceSequence faceneighs;
    JNode::getRelations(vtx, faceneighs);

    JEdgeSet   edgeSet;
    int numfaces = faceneighs.size();
    for( const JFacePtr &face : faceneighs ) {
        int numedges = face->getSize(1);
        for( int j = 0; j < numedges; j++) {
            const JEdgePtr &edge = face->getEdgeAt(j);
            const JNodePtr &v0   = edge->getNodeAt(0);
            const JNodePtr &v1   = edge->getNodeAt(1);
            if( v0 == vtx || v1 == vtx)
                edgeSet.insert(edge);
        }
    }

    edges.resize( edgeSet.size() );
    boost::copy(edgeSet, edges.begin());
}
////////////////////////////////////////////////////////////////////////////////

int JMeshTopology ::collectFaces()
{
    if( mesh == nullptr) return 1;

    size_t numcells = mesh->cells.size();
    if( numcells == 0) return 2;

    logger->setInfo("Collecting face simplices");

    JFaceSet  fset;
    JFaceSequence cellfaces;

    for (size_t icell = 0; icell < numcells; icell++) {
        const JCellPtr &cell = mesh->getCellAt(icell);
        if( cell->isActive() ) {
            cellfaces = cell->getFaces();
            fset.insert( cellfaces.begin(), cellfaces.end() );
        }
    }

    mesh->faces.clear();
    size_t numFaces = fset.size();
    if( numFaces ) {
        mesh->faces.resize(numFaces);
        std::copy( fset.begin(), fset.end(), mesh->faces.begin());
    }

    mesh->enumerate(2);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
void
JMeshTopology ::search_3d_boundary()
{
    size_t numfaces = mesh->getSize(2);
    int bmark, err;

    const string name = "Boundary";

    JCellSequence cellneighs;
    JEdgeSequence faceedges;
    #pragma omp parallel for private(faceedges,cellneighs)
    for( size_t iface = 0; iface < numfaces; iface++) {
        const JFacePtr &face = mesh->getFaceAt(iface);
        if( face->isActive() ) {
            JFace::getRelations( face, cellneighs);
            if( cellneighs.size() == 1) {
                err = face->getAttribute(name, bmark);
                if( err ) {
                    bmark = 1;
                    face->setAttribute( name, bmark);
                }
                faceedges = face->getEdges();
                for( size_t j = 0; j < faceedges.size(); j++) {
                    const JEdgePtr &edge = faceedges[j];
                    err = edge->getAttribute( name, bmark);
                    if( err ) {
                        bmark = 1;
                        edge->setAttribute( name, bmark);
                    }
                    const JNodePtr &v0 = edge->getNodeAt(0);
                    err = v0->getAttribute( name, bmark);
                    if( err ) {
                        bmark = 1;
                        v0->setAttribute( name, bmark);
                    }
                    const JNodePtr &v1 = edge->getNodeAt(1);
                    err = v1->getAttribute(name, bmark);
                    if( err ) {
                        bmark = 1;
                        v1->setAttribute( name, bmark);
                    }
                }
            } else
                face->deleteAttribute(name);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

void
JMeshTopology ::search_2d_boundary()
{
    size_t numnodes = mesh->getSize(0);
    size_t numedges = mesh->getSize(1);
    size_t numfaces = mesh->getSize(2);

    int bmark, err;

    const string name = "Boundary";

    bool  boundary_detected = 0;
    JNodeSet nodeSet;
    JFaceSequence faceneighs;

//  #pragma omp parallel for private(faceneighs)
    for( size_t iedge = 0; iedge < numedges; iedge++) {
        const JEdgePtr &edge = mesh->getEdgeAt(iedge);
        if( edge->isActive() ) {
            JEdge::getRelations(edge, faceneighs);
            assert( !faceneighs.empty() ) ;
            if( faceneighs.size() == 1) {
                boundary_detected = 1;
                err = edge->getAttribute(name, bmark);
                if( err ) {
                    bmark = 1;
                    edge->setAttribute( name, bmark);
                }
                const JNodePtr &v0 = edge->getNodeAt(0);
                nodeSet.insert(v0);
                err = v0->getAttribute(name, bmark);
                if( err ) {
                    bmark = 1;
                    v0->setAttribute( name, bmark);
                }
                const JNodePtr &v1 = edge->getNodeAt(1);
                nodeSet.insert(v1);
                err = v1->getAttribute(name, bmark);
                if( err ) {
                    bmark = 1;
                    v1->setAttribute( name, bmark);
                }
            } else
                edge->deleteAttribute(name);
        }
    }

    if( boundary_detected ) {
        size_t numnodes = mesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->hasAttribute(name) ) {
                if( nodeSet.find(vtx) == nodeSet.end() )
                    vtx->deleteAttribute(name);
            }
        }
    }

    if( !boundary_detected ) {
        boundary_detected = 1;
        for( size_t iedge = 0; iedge < numedges; iedge++) {
            const JEdgePtr &edge = mesh->getEdgeAt(iedge);
            if( edge->isActive() ) {
                JEdge::getRelations(edge, faceneighs);
                if( faceneighs.size() != 2) {
                    boundary_detected = 0;
                    break;
                }
            }
        }

        if( boundary_detected ) {
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &face = mesh->getFaceAt(i);
                int err = face->getAttribute(name, bmark);
                if( err) {
                    bmark = 1;
                    face->setAttribute(name, bmark);
                }
            }
            for( size_t i = 0; i < numedges; i++) {
                const JEdgePtr &edge = mesh->getEdgeAt(i);
                int err = edge->getAttribute(name, bmark);
                if( err) {
                    bmark = 1;
                    edge->setAttribute(name, bmark);
                }
            }
            for( size_t i = 0; i < numedges; i++) {
                const JEdgePtr &edge = mesh->getEdgeAt(i);
                int err = edge->getAttribute(name, bmark);
                if( err) {
                    bmark = 1;
                    edge->setAttribute(name, bmark);
                }
            }
            for( size_t i = 0; i < numnodes; i++) {
                const JNodePtr &vtx = mesh->getNodeAt(i);
                int err = vtx->getAttribute(name, bmark);
                if( err) {
                    bmark = 1;
                    vtx->setAttribute(name, bmark);
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

int JMeshTopology ::searchBoundary()
{
    logger->setInfo("Searching boundary in the mesh");

    int topDim = getDimension();
    if( topDim == 3) search_3d_boundary();
    if( topDim == 2) search_2d_boundary();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getDoublets( JNodeSequence &doublets)
{
    doublets.clear();

    int topDim = getDimension();
    if( topDim != 2 ) return;

    logger->setInfo("Searching doublets the mesh ");

    size_t numnodes = mesh->getSize(0);
    mesh->buildRelations(0,topDim);

    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            if( !vertex->isBoundary() ) {
                if( vertex->getNumRelations(topDim) == 2)
                    doublets.push_back(vertex);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

int JMeshTopology :: setMinColors()
{
    if( mesh == nullptr) return 0;

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        setMinColor( vtx );
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getDoublets( JEdgeSequence &doublets)
{
    doublets.clear();

    int topDim = getDimension();

    if( topDim != 3 ) return;
    size_t numedges = mesh->getSize(1);

    logger->setInfo("Searching doublets the mesh");

    mesh->buildRelations(1,topDim);

    for( size_t i = 0; i < numedges; i++) {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            if( !edge->isBoundary() ) {
                if( edge->getNumRelations(topDim) == 2)
                    doublets.push_back(edge);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JMeshTopology :: getSurfaceMesh()
{
    if( mesh == nullptr) return nullptr;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 3 ) {
        cout << "Warning: The given mesh has no 3D elements: Surface mesh not extracted" << endl;
        return nullptr;
    }

    logger->setInfo("Collecting surface mesh");

    searchBoundary();

    JMeshPtr smesh = JMesh::newObject();

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && vtx->isBoundary() ) smesh->addObject( vtx );
    }

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive()  && face->isBoundary() )
            smesh->addObject( face );
    }

    smesh->getTopology()->collectEdges();

    return smesh;
}

////////////////////////////////////////////////////////////////////////////////

int JMeshTopology ::getIrregularNodes(JNodeSequence &seq, int from_where)
{
    logger->setInfo("Collecting irregular nodes in the mesh");

    seq.clear();

    int internal_regular_count = 0;
    int bound_regular_count;

    int topDim = getDimension();
    int elem_type;

    size_t numnodes = mesh->getSize(0);
    int    nSize;

    if( topDim == 2 ) {
        elem_type = getElementsType(2);
        switch( elem_type) {
        case JFace::TRIANGLE:
            internal_regular_count = 6;
            bound_regular_count    = 3;
            break;
        case JFace::QUADRILATERAL:
            internal_regular_count = 4;
            bound_regular_count    = 2;
            break;
        default:
            logger->setError("Irregular nodes not supported for inhomogeneous mesh ");
            return 1;
        }

        mesh->buildRelations(0,2);
        searchBoundary();

        // Query from the boundary nodes ...
        if (from_where == 1) {
            for (size_t i = 0; i < numnodes; i++) {
                const JNodePtr &v = mesh->getNodeAt(i);
                if( v->isActive() ) {
                    nSize  = v->getNumRelations(2);
                    if (v->isBoundary() && nSize != bound_regular_count) seq.push_back(v);
                }
            }
        }

        // Query from the internal nodes ...
        if (from_where == 0) {
            for (size_t i = 0; i < numnodes; i++) {
                const JNodePtr &v = mesh->getNodeAt(i);
                if( v->isActive() ) {
                    nSize  = v->getNumRelations(2);
                    if (!v->isBoundary() &&  nSize != internal_regular_count) seq.push_back(v);
                }
            }
        }
    }

    if( topDim == 3 ) {
        elem_type = getElementsType(3);
        switch( elem_type) {
        case JCell::TETRAHEDRON:
            internal_regular_count = 6;
            break;
        case JCell::HEXAHEDRON:
            internal_regular_count = 6;
            break;
        default:
            logger->setWarn("Irregular nodes not supported for inhomogeneous mesh ");
            return 1;
        }
        mesh->buildRelations(0, 3);
    }

    searchBoundary();

    // Query from the internal nodes ...
    if (from_where == 0) {
        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() ) {
                nSize  = v->getNumRelations(topDim);
                if (!v->isBoundary() &&  (nSize != internal_regular_count)) seq.push_back(v);
            }
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshTopology :: getEulerCharacteristic()
{
    logger->setInfo("Calculating Euler Characterisistic of the mesh");

    int topDim = getDimension();

    size_t  V = mesh->getActiveSize(0);
    size_t  E = mesh->getActiveSize(1);
    size_t  F = mesh->getActiveSize(2);
    size_t  C = mesh->getActiveSize(3);

    if( topDim == 2 ) {
        if( E == 0 ) {
            collectEdges();
            E = mesh->getActiveSize(1);
        }
    }

    if( topDim == 3 ) {
        if( E == 0 ) {
            collectEdges();
            E = mesh->getActiveSize(1);
        }

        if( F == 0 ) {
            collectFaces();
            F = mesh->getActiveSize(2);
        }
    }

    return V - E + F - C;
}

///////////////////////////////////////////////////////////////////////////////
int JMeshTopology :: getEulerCharacteristic( const JFaceSequence &fseq)
{
    JEdgeSequence eseq;
    JMeshTopology::getEntitySet( fseq, eseq);

    JNodeSequence vseq;
    getEntitySet( fseq, vseq);

    size_t  V = vseq.size();
    size_t  E = eseq.size();
    size_t  F = fseq.size();
    return V - E + F;
}
///////////////////////////////////////////////////////////////////////////////

int JMeshTopology :: getBettiNumber( BettiMethod method)
{
    if( mesh == nullptr) return 0;
#ifdef USE_IGL

    JMeshEigenMatrix mat;
    mat.setMesh(mesh);
    Eigen::MatrixXi T = mat.getElementMatrix();

    Engine* engine;
    igl::matlab::mlinit(&engine);

    // Send Laplacian matrix to matlab
    igl::matlab::mlsetmatrix(&engine,"T",T);

//    double numloops = 0.0;
//    igl::matlab::mlsetscalar(&engine,"numloops",numloops);

    // Calling Yaron Lipman's Matlab code ....
    ostringstream oss;

    if( method == BettiMethod::SMITH_NORMAL_FORM )
        oss << "numloops = betti(T, 'Method', 'smith-normal-form')";

    if( method == BettiMethod::SPARSE_RANK )
        oss << "numloops = betti(T, 'Method', 'rank')";

    if( method == BettiMethod::FULL_RANK )
        oss << "numloops = betti(T, 'Method', 'full-rank')";

    if( method == BettiMethod::QR )
        oss << "numloops = betti(T, 'Method', 'qr')";

    if( method == BettiMethod::EIGENANALYSIS )
        oss << "numloops = betti(T,'Method', 'eigenanalysis')";

    string cmd = oss.str();
    igl::matlab::mleval(&engine, cmd);

    double numloops = igl::matlab::mlgetscalar(&engine,"numloops");
    igl::matlab::mlclose(&engine);
    return numloops;
    return 0;
#endif
}

///////////////////////////////////////////////////////////////////////////////

bool
JMeshTopology ::isConsistent() const
{
    logger->setInfo("Checking mesh consistency");

    int consistent = 1;

    int topDim = getDimension();

    if( topDim == 2 ) {
        mesh->buildRelations(1, 2);
        size_t numedges  = mesh->getSize(1);
        JFaceSequence neighs;
        for( size_t i = 0; i <  numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                JEdge::getRelations( edge, neighs );
                if (neighs.size() == 2) {
                    int dir1 = neighs[0]->getOrientation(edge);
                    int dir2 = neighs[1]->getOrientation(edge);
                    if (dir1 * dir2 == 1) {
                        logger->setWarn("Mesh is not consistently oriented" );
                        consistent = 0;
                        break;
                    }
                }
            }
        }
        return consistent;
    }

    size_t numfaces  = mesh->getSize(2);

    JCellSequence neighs;
    for (size_t j = 0; j < numfaces; j++) {
        const JFacePtr &face = mesh->getFaceAt(j);
        if( face->isActive() ) {
            JFace::getRelations(face, neighs );
            if (neighs.size() == 2) {
                int dir1 = neighs[0]->getOrientation(face);
                int dir2 = neighs[1]->getOrientation(face);
                if (dir1 * dir2 == 1) {
                    logger->setWarn("Mesh is not consistently oriented");
                    cout << "Debug: Face " << endl;
                    face->print();
                    cout << "Cell 1 " << endl;
                    neighs[0]->print();
                    cout << "Cell 2 " << endl;
                    neighs[1]->print();
                    consistent = 0;
                    break;
                }
            }
        }
    }
    return consistent;
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshTopology ::getConsistent()
{
    logger->setInfo("Building consistent mesh");

    int topDim = getDimension();
    JFacePtr face;
    JFaceSequence faceneighs;
    if( topDim == 2 ) {
        size_t numfaces = mesh->getSize(2);
        for (size_t iface = 0; iface < numfaces; iface++) {
            face = mesh->getFaceAt(iface);
            face->setVisitBit(0);
        }

        deque<JFacePtr> faceQ;
        for (size_t iface = 0; iface < numfaces; iface++) {
            face = mesh->getFaceAt( iface );
            if( face->isActive() ) {
                faceQ.push_back(face);
                break;
            }
        }

        while (!faceQ.empty()) {
            JFacePtr face = faceQ.front();
            faceQ.pop_front();
            if (!face->getVisitBit()) {
                face->setVisitBit(1);
                int nedges = face->getSize(0);
                for (int j = 0; j < nedges; j++) {
                    JEdgePtr edge = face->getEdgeAt(j);
                    JEdge::getRelations(edge, faceneighs);
                    if (faceneighs.size() == 2) {
                        int dir1 = faceneighs[0]->getOrientation(edge);
                        int dir2 = faceneighs[1]->getOrientation(edge);
                        if (dir1 * dir2 == 1) {
                            if (!faceneighs[0]->getVisitBit() && faceneighs[1]->getVisitBit())
                                faceneighs[0]->reverse();

                            if (!faceneighs[1]->getVisitBit() && faceneighs[0]->getVisitBit())
                                faceneighs[1]->reverse();
                        }
                        faceQ.push_back(faceneighs[0]);
                        faceQ.push_back(faceneighs[1]);
                    }
                }
            }
        }

        for (size_t iface = 0; iface < numfaces; iface++) {
            face = mesh->getFaceAt(iface);
            if( face->isActive() && !face->getVisitBit()) {
                cout << "Error: face not visited : " << face->getID() << " "
                     << face->getVisitBit() << endl;
            }
        }
    }


    JCellPtr cell;
    JCellSequence cellneighs;

    if( topDim == 3 ) {
        size_t numcells = mesh->getSize(3);
        for (size_t icell = 0; icell < numcells; icell++) {
            const JCellPtr &c = mesh->getCellAt(icell);
            c->setVisitBit(0);
        }

        deque<JCellPtr> cellQ;
        for (size_t icell = 0; icell < numcells; icell++) {
            const JCellPtr &c = mesh->getCellAt( icell );
            if( c->isActive() ) {
                cellQ.push_back(c);
                break;
            }
        }

        set<JCellPtr> cset;
        while (!cellQ.empty()) {
            JCellPtr cell = cellQ.front();
            cellQ.pop_front();
            if (cell->getVisitBit() == 0) {
                cell->setVisitBit(1);
                cset.insert(cell);
                int nfaces = cell->getSize(2);
                for (int j = 0; j < nfaces; j++) {
                    const JFacePtr &face = cell->getFaceAt(j);
                    JFace::getRelations(face, cellneighs);
                    if (cellneighs.size() == 2) {
                        int dir1 = cellneighs[0]->getOrientation(face);
                        int dir2 = cellneighs[1]->getOrientation(face);
                        if (dir1 * dir2 == 1) {
                            if (!cellneighs[0]->getVisitBit() && cellneighs[1]->getVisitBit())
                                cellneighs[0]->reverse();

                            if (!cellneighs[1]->getVisitBit() && cellneighs[0]->getVisitBit())
                                cellneighs[1]->reverse();
                        }
                        cellQ.push_back(cellneighs[0]);
                        cellQ.push_back(cellneighs[1]);
                    }
                }
            }
        }

        for (size_t icell = 0; icell < numcells; icell++) {
            cell = mesh->getCellAt(icell);
            if( cell->isActive() && cell->getVisitBit() == 0) {
                cout << "Error: cell not visited : " << cell->getID() << " "
                     << cell->getVisitBit() << endl;
            }
        }
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

bool JMeshTopology :: isClosed()
{
    logger->setInfo("Checking if the mesh is closed");

    int topDim = getDimension();
    if( topDim == 2 ) {
        size_t numedges = mesh->getSize(1);
        mesh->buildRelations(1,2);
        for( size_t iedge = 0; iedge < numedges; iedge++) {
            const JEdgePtr &edge = mesh->getEdgeAt(iedge);
            if( edge->isActive() ) {
                if( edge->getNumRelations(2) < 2 ) return 0;
            }
        }
        return 1;
    }

    logger->setWarn("JMeshTopology::isClosed() not implemented");
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getBoundary( JNodeSequence &bnodes)
{
    bnodes.clear();
    JEdgeSequence bedges;
    getBoundary(bedges);
    getEntitySet(bedges, bnodes);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getInternal( JNodeSequence &inodes)
{
    assert( mesh != nullptr);
    inodes.clear();

    // Collect all the nodes in the mesh
    JNodeSequence allnodes = mesh->getNodes();
    if( allnodes.empty() ) return;

    boost::sort( allnodes );

    // Collect all the boundary nodes in the mesh
    JNodeSequence bdnodes;
    getBoundary(bdnodes);
    boost::sort( bdnodes );

    // Internal nodes = allnodes - boundary nodes ...
    boost::set_difference(allnodes,bdnodes, std::back_inserter(inodes));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getBoundary( JEdgeSequence &bedges)
{
    assert( mesh != nullptr);

    bedges.clear();

    if( getDimension() == 1 ) {
        bedges = mesh->getEdges();
        return;
    }

    JFaceSequence faceneighs;

    if( getDimension() == 2 ) {
        size_t numEdges = mesh->getSize(1);
        // First check if the mesh is topological disk ...
        for( size_t i = 0; i < numEdges; i++) {
            const JEdgePtr &e = mesh->getEdgeAt(i);
            if( e->isActive() ) {
                JEdge::getRelations(e, faceneighs);
                if( faceneighs.size() == 1 ) bedges.push_back(e);
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getInternal( JEdgeSequence &iedges)
{
    assert( mesh != nullptr);
    iedges.clear();

    // Collect all the edges in the mesh
    JEdgeSequence alledges = mesh->getEdges();

    if( alledges.empty() ) return;
    boost::sort( alledges );

    // Collect all the boundary edges in the mesh
    JEdgeSequence bdedges;
    getBoundary(bdedges);
    boost::sort(bdedges );

    // Internal edges = alledges - boundary edges ...
    boost::set_difference(alledges, bdedges, std::back_inserter(iedges));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getInternal( JEdgeSequence &iedges, JNodeSequence &inodes)
{
    assert( mesh != nullptr);
    getInternal(iedges);
    getInternal(inodes);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getBoundary( vector<JEdgeSequence> &edgeloops)
{
    edgeloops.clear();

    JEdgeSequence edges;

    int dim = getDimension();

    if( dim == 1 ) edges = mesh->getEdges();
    if( dim == 2 ) getBoundary(edges);

    JEdgeTopology::getLoops( edges, edgeloops);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getSubmeshBoundary(  const JFaceSequence &faces, JEdgeSequence &boundary)
{
    boundary.clear();

    JFaceSet faceSet;
    for( const JFacePtr &f: faces)
        faceSet.insert(f);

    map<JEdgePtr, set<JFacePtr>>  edgefaces;

    JFaceSequence faceneighs;
    for( const JFacePtr &face: faces) {
        int numedges = face->getSize(0);
        for( int i = 0; i < numedges; i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            JEdge::getRelations(edge, faceneighs);
            for( const JFacePtr &fneigh : faceneighs) {
                if( faceSet.find(fneigh) != faceSet.end() )
                    edgefaces[edge].insert(fneigh);
            }
        }
    }

    for( const auto &keyVal : edgefaces)
        if( keyVal.second.size() == 1) boundary.push_back(keyVal.first);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getSubmeshInternal( const JFaceSequence &faces, JEdgeSequence &iedges, JNodeSequence &inodes)
{
    JFaceSet faceSet;
    for( const JFacePtr &f: faces)
        faceSet.insert(f);

    map<JEdgePtr, set<JFacePtr>>  edgefaces;

    JFaceSequence faceneighs;
    for( const JFacePtr &face: faces) {
        int numedges = face->getSize(0);
        for( int i = 0; i < numedges; i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            JEdge::getRelations(edge, faceneighs);
            for( const JFacePtr &fneigh : faceneighs) {
                if( faceSet.find(fneigh) != faceSet.end() )
                    edgefaces[edge].insert(fneigh);
            }
        }
    }

    iedges.clear();

    JNodeSet boundnodes;
    for( const auto &keyVal : edgefaces) {
        const JEdgePtr &edge = keyVal.first;
        if( keyVal.second.size() == 2)
            iedges.push_back(edge);
        else {
            boundnodes.insert( edge->getNodeAt(0));
            boundnodes.insert( edge->getNodeAt(1));
        }
    }

    JNodeSet inodeSet;

    for( const JEdgePtr &e : iedges) {
        const JNodePtr &v0 = e->getNodeAt(0);
        const JNodePtr &v1 = e->getNodeAt(1);
        if( boundnodes.find(v0) == boundnodes.end() )
            inodeSet.insert(v0);
        if( boundnodes.find(v1) == boundnodes.end() )
            inodeSet.insert(v1);
    }

    inodes.clear();
    for( const JNodePtr &v : inodeSet)
        inodes.push_back(v);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getSubmeshInternal( const JFaceSequence &faces, JNodeSequence &inodes)
{
    JEdgeSequence dummy;
    getSubmeshInternal( faces, dummy, inodes);
}

///////////////////////////////////////////////////////////////////////////////

JEdgeSequence JMeshTopology :: getNonManifoldEdges() const
{
    JEdgeSequence bedges;
    if( getDimension() != 2 ) {
        logger->setWarn("The input mesh must be a surface");
        return bedges;
    }

    mesh->buildRelations(1,2);

    size_t numEdges = mesh->getSize(1);

    JFaceSequence faceneighs;
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            if( edge->getNumRelations(2) > 2) bedges.push_back(edge);
        }
    }
    return bedges;
}
///////////////////////////////////////////////////////////////////////////////

JFaceSequence JMeshTopology :: getNonManifoldFaces() const
{
    JFaceSequence bfaces;
    if( getDimension() != 3 ) {
        logger->setWarn("The input mesh must be a volume");
        return bfaces;
    }

    mesh->buildRelations(2,3);

    size_t numFaces = mesh->getSize(2);

    JCellSequence cellneighs;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            if( face->getNumRelations(3) > 2) bfaces.push_back(face);
        }
    }
    return bfaces;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getBoundary( JFaceSequence &bfaces)
{
    bfaces.clear();

    if( getDimension() == 2 ) {
        JEdgeSequence edges;
        getBoundary( edges );
        if( !edges.empty() ) {
            JFaceSet faceSet;
            JFaceSequence faceneighs;
            for( const JEdgePtr &edge : edges) {
                JEdge::getRelations(edge, faceneighs);
                if( faceneighs.size() == 1)
                    faceSet.insert( faceneighs[0] );
            }
            for( const JFacePtr  &face : faceSet )
                bfaces.push_back(face);
            return;
        }
    }

    size_t numCells = mesh->getSize(3);
    if( numCells == 0) return;

    collectFaces();
    mesh->buildRelations(2,3);

    size_t numFaces = mesh->getSize(2);

    JCellSequence cellneighs;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            JFace::getRelations(face, cellneighs);
            if( cellneighs.size() == 1 ) {
                bfaces.push_back(face);
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getInternal( JFaceSequence &ifaces)
{
    ifaces.clear();
    if( mesh == nullptr) return;

    // Collect all the edges in the mesh
    JFaceSequence allfaces = mesh->getFaces();
    if( allfaces.empty() ) return;

    boost::sort(allfaces);

    // Collect all the boundary edges in the mesh
    JFaceSequence bdfaces;
    getBoundary(bdfaces);
    boost::sort( bdfaces );

    // Internal edges = alledges - boundary edges ...
    boost::set_difference(allfaces, bdfaces, std::back_inserter(ifaces));
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTopology :: getDuplicates( JEdgeSequence &dedges)
{
    dedges.clear();
    map<JNodePtr, JEdgeSequence> vmap;
    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        const JNodePtr &vhash = edge->getHashNode();
        if( vmap.find(vhash) == vmap.end() ) {
            vmap[vhash].push_back(edge);
        } else {
            if(find( vmap[vhash].begin(), vmap[vhash].end(), edge) == vmap[vhash].end())
                vmap[vhash].push_back(edge);
            else
                dedges.push_back( edge);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getDuplicates( JFaceSequence &dfaces)
{
    dfaces.clear();
    map<JNodePtr, JFaceSequence> vmap;
    size_t numFaces = mesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        const JNodePtr &vhash = face->getHashNode();
        if( vmap.find(vhash) == vmap.end() ) {
            vmap[vhash].push_back(face);
        } else {
            if(find( vmap[vhash].begin(), vmap[vhash].end(), face) == vmap[vhash].end())
                vmap[vhash].push_back(face);
            else
                dfaces.push_back( face );
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getDuplicates( JCellSequence &dcells)
{
    dcells.clear();
    map<JNodePtr, JCellSequence> vmap;
    size_t numCells = mesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell  = mesh->getCellAt(i);
        const JNodePtr &vhash = cell->getHashNode();
        if( vmap.find(vhash) == vmap.end() ) {
            vmap[vhash].push_back(cell);
        } else {
            if(find( vmap[vhash].begin(), vmap[vhash].end(), cell) == vmap[vhash].end())
                vmap[vhash].push_back(cell);
            else
                dcells.push_back( cell );
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

int
JMeshTopology ::getElementsType(int dim) const
{
    if( mesh == nullptr) return 0;

    int maxnodes = 0;
    int minnodes = std::numeric_limits<int>::max();

    size_t nSize;

    if( dim == 2 ) {
        nSize = mesh->getSize(2);
        #pragma omp parallel for reduction(max:maxnodes) reduction(min:minnodes)
        for (size_t i = 0; i < nSize; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if(f->isActive()) {
                maxnodes = max(maxnodes, f->getSize(0));
                minnodes = min(minnodes, f->getSize(0));
            }
        }
        if( minnodes != maxnodes ) return 0;
        if( minnodes == 3) return JFace::TRIANGLE;
        if( minnodes == 4) return JFace::QUADRILATERAL;
    }

    if( dim == 3 ) {
        size_t nSize = mesh->getSize(3);
        #pragma omp parallel for reduction(max:maxnodes) reduction(min:minnodes)
        for (size_t i = 0; i < nSize; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            if (c->isActive()) {
                maxnodes = max(maxnodes, c->getSize(0));
                minnodes = min(minnodes, c->getSize(0));
            }
        }
        if( minnodes != maxnodes ) return 0;
        if( minnodes == 4) return JCell::TETRAHEDRON;
        if( minnodes == 8) return JCell::HEXAHEDRON;
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

size_t JMeshTopology :: getBoundarySize( int d )  const
{
    size_t nSize, nCount = 0;

    if( d == 0) {
        nSize = mesh->getSize(0);
        for( size_t i = 0; i < nSize; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() && v->isBoundary() ) nCount++;
        }
        return nCount;
    }

    if( d == 1) {
        nSize = mesh->getSize(1);
        for( size_t i = 0; i < nSize; i++) {
            const JEdgePtr &e = mesh->getEdgeAt(i);
            if( e->isActive() && e->isBoundary() ) nCount++;
        }
        return nCount;
    }

    if( d == 2) {
        nSize = mesh->getSize(2);
        for( size_t i = 0; i < nSize; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            if( f->isActive() && f->isBoundary() ) nCount++;
        }
        return nCount;
    }

    if( d == 3) {
        nSize = mesh->getSize(3);
        for( size_t i = 0; i < nSize; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            if( c->isActive() && c->isBoundary() ) nCount++;
        }
        return nCount;
    }

    if( nCount == 0)
        cout << "Warning: No boundary elements found" << endl;

    return nCount;
}

////////////////////////////////////////////////////////////////////////////////

map<int,size_t> JMeshTopology :: getNodesDegreeDistribution(int where ) const
{
    map<int,size_t>  degree;

    int d;
    int topDim = getDimension();

    size_t  numnodes = mesh->getSize(0);
    mesh->buildRelations(0, topDim);

    size_t ncount = 0;

    if( where == JMeshEntity::ANY_ENTITY) {
        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() ) {
                d = v->getNumRelations(topDim );
                if( degree.find(d) == degree.end() )
                    degree[d] = 0;
                degree[d]++;
                ncount++;
            }
        }
    }

    if( where == JMeshEntity::BOUNDARY_ENTITY) {
        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() && v->isBoundary() ) {
                d = v->getNumRelations(topDim );
                if( degree.find(d) == degree.end() )
                    degree[d] = 0;
                degree[d]++;
                ncount++;
            }
        }
    }

    if( where == JMeshEntity::INTERNAL_ENTITY) {
        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( v->isActive() && !v->isBoundary() ) {
                d = v->getNumRelations(topDim );
                if( degree.find(d) == degree.end() )
                    degree[d] = 0;
                degree[d]++;
                ncount++;
            }
        }
    }

    return degree;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshTopology::removeOrphaned()
{
    /*
        if( mesh == nullptr) return;

        size_t numnodes = mesh->getSize(0);
        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            v->setVisitBit(0);
        }

        size_t numedges = mesh->getSize(1);
        for (size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            const JNodePtr &v0 = edge->getNodeAt(0);
            const JNodePtr &v1 = edge->getNodeAt(1);
            v0->setVisitBit(1);
            v1->setVisitBit(1);
        }

        size_t numfaces = mesh->getSize(2);
        for (size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            for (int j = 0; j < face->getSize(0); j++) {
                JNodePtr v = face->getNodeAt(j);
                v->setVisitBit(1);
            }
        }

        size_t numcells = mesh->getSize(3);
        for (size_t i = 0; i < numcells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            for (int j = 0; j < cell->getSize(0); j++) {
                JNodePtr v = cell->getNodeAt(j);
                v->setVisitBit(1);
            }
        }

        for (size_t i = 0; i < numnodes; i++) {
            const JNodePtr &v = mesh->getNodeAt(i);
            if( !v->getVisitBit() ) v->setStatus(JMeshEntity::REMOVE);
        }
        mesh->pruneNodes();
        JFaceSequence cellfaces;
        size_t numCells = mesh->getSize(3);
        for( size_t i = 0; i < numCells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            cell->remove_unattached_lower_entities();
        }

        JEdgeSequence faceedges;
        size_t numFaces = mesh->getSize(2);
        for( size_t i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            face->remove_unattached_lower_entities();
        }

        size_t numEdges = mesh->getSize(1);
        for( size_t i = 0; i < numEdges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            edge->remove_unattached_lower_entities();
        }
    */

    /*
       size_t numNodes = mesh->getSize(0);
       for( size_t i = 0; i < numNodes; i++) {
            Node *node = mesh->getNodeAt(i);
            node->remove_unattached_lower_entities();
       }
    */

    mesh->pruneAll();
}
////////////////////////////////////////////////////////////////////////////////

int
JMeshTopology::getNumComponents(bool )
{
    int  tDim  = getDimension();

    int numComponents = 0;

    if( tDim == 2 ) {
        JFacePtr face;
        deque<JFacePtr> faceQ;
        JFaceSequence fneighs;

        size_t numfaces = mesh->getSize(2);

        for (size_t iface = 0; iface < numfaces; iface++) {
            face = mesh->getFaceAt(iface);
            face->setID(iface);
            face->setVisitBit(0);
        }

        while (1) {
            face.reset();
            faceQ.clear();
            for (size_t iface = 0; iface < numfaces; iface++) {
                face = mesh->getFaceAt(iface);
                if (face->isActive() && !face->getVisitBit()) {
                    face->setAttribute("Component", numComponents);
                    faceQ.push_back(face);
                    break;
                }
            }

            if (faceQ.empty())
                break;

            while (!faceQ.empty()) {
                JFacePtr face = faceQ.front();
                faceQ.pop_front();
                if (!face->getVisitBit()) {
                    face->setAttribute( "Component", numComponents);
                    face->setVisitBit(1);
                    int nedges = face->getSize(1);
                    for (int j = 0; j < nedges; j++) {
                        JEdgePtr edge = face->getEdgeAt(j);

                        int proceed = 1;
                        /*
                        if (stop_at_interface) {
                             if (v0->isConstrained() && v1->isConstrained()) {
                                  Edge edge(v0, v1);
                                  if (hasFeatureEdge(edge)) proceed = 0;
                             }
                        }
                        */

                        if (proceed) {
                            JEdge::getRelations( edge, fneighs );
                            int numneighs = fneighs.size();
                            for( int k = 0; k < numneighs; k++) {
                                if( !fneighs[k]->getVisitBit() )
                                    faceQ.push_back(fneighs[k]);
                            }
                        }
                    }
                }
            } // Complete one Component
            numComponents++;
        }

        for (size_t iface = 0; iface < numfaces; iface++) {
            face = mesh->getFaceAt(iface);
            if (face->isActive() && !face->getVisitBit())
                logger->setError("Component search incomplete" );
        }
    }

    JCellPtr cell;
    JCellSequence cneighs;
    if( tDim == 3 ) {
        cell.reset();
        deque<JCellPtr> cellQ;

        size_t numcells = mesh->getSize(3);

        for (size_t icell = 0; icell < numcells; icell++) {
            cell = mesh->getCellAt(icell);
            cell->setID(icell);
            cell->setVisitBit(0);
        }


        while (1) {
            cell.reset();
            cellQ.clear();
            for (size_t icell = 0; icell < numcells; icell++) {
                cell = mesh->getCellAt(icell);
                if( cell->isActive() &&!cell->getVisitBit()) {
                    cell->setAttribute("Component", numComponents);
                    cellQ.push_back(cell);
                    break;
                }
            }

            if (cellQ.empty())
                break;

            while (!cellQ.empty()) {
                JCellPtr cell = cellQ.front();
                cellQ.pop_front();
                if (!cell->getVisitBit()) {
                    cell->setAttribute( "Component", numComponents);
                    cell->setVisitBit(1);
                    int nfaces = cell->getSize(2);
                    for (int j = 0; j < nfaces; j++) {
                        JFacePtr face = cell->getFaceAt(j);

                        int proceed = 1;
                        /*
                        if (stop_at_interface) {
                             if (v0->isConstrained() && v1->isConstrained()) {
                                  Edge edge(v0, v1);
                                  if (hasFeatureEdge(edge)) proceed = 0;
                             }
                        }
                        */

                        if (proceed) {
                            JFace::getRelations( face, cneighs );
                            int numneighs = cneighs.size();
                            for( int k = 0; k < numneighs; k++) {
                                if( !cneighs[k]->getVisitBit() )
                                    cellQ.push_back(cneighs[k]);
                            }
                        }

                    }
                }
            } // Complete one Component
            numComponents++;
        }

        for (size_t icell = 0; icell < numcells; icell++) {
            cell = mesh->getCellAt(icell);
            if( cell->isActive() && !cell->getVisitBit())
                logger->setError("Component search incomplete " );
        }
    }

    return numComponents;
}

////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getOrdered(const JFacePtr &seedface, JFaceSequence &bfs, int order, size_t maxSize)
{

    if( order == JMesh::BREADTH_FIRST_ORDER ) {
        mesh->buildRelations(1, 2);

        bfs.clear();

        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            f->setVisitBit(0);
        }

        list<JFacePtr> faceQ;
        faceQ.push_back(seedface);

        JFaceSequence neighs;

        while(1) {
            while(!faceQ.empty() ) {
                JFacePtr currface = faceQ.front();
                faceQ.pop_front();
                if( currface->isActive() ) {
                    JFace::getRelations12(currface, neighs);
                    for( size_t i = 0; i < neighs.size(); i++)
                        if( !neighs[i]->getVisitBit() ) faceQ.push_back(neighs[i] );
                    if( !currface->getVisitBit() ) {
                        bfs.push_back(currface);
                        if( bfs.size() == maxSize ) return;
                        currface->setVisitBit(1);
                    }
                }
            }

            for( size_t i = 0; i < numfaces; i++) {
                JFacePtr f = mesh->getFaceAt(i);
                if( !f->getVisitBit() && f->isActive() ) {
                    faceQ.push_back( f );
                    break;
                }
            }
            if( faceQ.empty() ) break;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JMeshTopology :: getOrdered(const JCellPtr &seedcell, JCellSequence &bfs, int order, size_t maxSize)
{
    mesh->buildRelations(2, 3);

    if( order == JMesh::BREADTH_FIRST_ORDER ) {
        bfs.clear();
        size_t numcells = mesh->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            const JCellPtr &c = mesh->getCellAt(i);
            c->setVisitBit(0);
        }

        list<JCellPtr> cellQ;
        cellQ.push_back(seedcell);

        JCellSequence neighs;

        while(1) {
            while(!cellQ.empty() ) {
                JCellPtr currcell = cellQ.front();
                cellQ.pop_front();
                if( currcell->isActive() ) {
                    JCell::getRelations23( currcell, neighs);
                    for( size_t i = 0; i < neighs.size(); i++)
                        if( !neighs[i]->getVisitBit() ) cellQ.push_back(neighs[i] );
                    if( !currcell->getVisitBit() ) {
                        bfs.push_back(currcell);
                        if( bfs.size() == maxSize ) return;
                        currcell->setVisitBit(1);
                    }
                }
            }

            for( size_t i = 0; i < numcells; i++) {
                JCellPtr c = mesh->getCellAt(i);
                if( !c->getVisitBit() && c->isActive() ) {
                    cellQ.push_back( c );
                    break;
                }
            }
            if( cellQ.empty() ) break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

int
JMeshTopology::set_nodes_wavefront(const JNodePtr &vsrc)
{
    int relexist = mesh->buildRelations(0, 0);

    size_t numNodes = mesh->getSize(0);
    for (size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->setAttribute("Layer", -1);
        vtx->setVisitBit(0);
    }

    JNodeSequence vertexQ;
    JNodeSet  nextQ;

    if( vsrc == nullptr) {
        searchBoundary();
        for (size_t i = 0; i < numNodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if (vtx->isBoundary() && !vtx->isRemoved() ) {
                vtx->setAttribute("Layer", 0);
                vtx->setVisitBit(1);
                vertexQ.push_back(vtx);
            }
        }
    } else
        vertexQ.push_back(vsrc);

    if (vertexQ.empty()) return 0;

    JNodeSequence vneighs;
    int layerid = 1;
    size_t nSize;
    while (!vertexQ.empty()) {
        nextQ.clear();
        nSize = vertexQ.size();
        for (size_t j = 0; j < nSize; j++) {
            JNodePtr currNode = vertexQ[j];
            JNode::getRelations( currNode, vneighs );
            for (size_t i = 0; i < vneighs.size(); i++) {
                JNodePtr vn = vneighs[i];
                if (!vn->getVisitBit()) nextQ.insert(vn);
            }
        }

        nSize = nextQ.size();
        if( nSize == 0) break;

        JNodeSet::const_iterator it;
        for( it = nextQ.begin(); it != nextQ.end(); ++it) {
            JNodePtr currNode = *it;
            currNode->setAttribute("Layer", layerid);
            currNode->setVisitBit(1);
        }

        vertexQ.clear();
        vertexQ.reserve( nSize );
        for( it = nextQ.begin(); it != nextQ.end(); ++it)
            vertexQ.push_back( *it );
        layerid++;
    }

    if (!relexist)
        mesh->clearRelations(0, 0);

    cout << "Info: Number of layers in node front : " << layerid - 1 << endl;

    return layerid - 1;
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshTopology::set_faces_wavefront(int start_from_boundary)
{
    if( mesh == nullptr) return;

    mesh->buildRelations(0, 2);

    if( start_from_boundary)
        searchBoundary();

    /*
             size_t numFaces = getSize(2);

             for (size_t i = 0; i < numFaces; i++) {
                  Face *f = getFaceAt(i);
                  f->setAttribute("Layer", 0);
                  f->setVisitBit(0);
             }

             JFaceSequence faceQ, nextQ, neighs;
             for (size_t i = 0; i < numFaces; i++) {
                  Face *f = getFaceAt(i);
                  if( !f->isRemoved() ) {
                       f->setAttribute("Layer", 1);
                       if (f->has_boundary_edge()) {
                            f->setAttribute("Layer", 0);
                            f->setVisitBit(1);
                            faceQ.push_back(f);
                       }
                  }
             }

             int layerid = 1;
             size_t nSize;

             while (!faceQ.empty()) {
                  nextQ.clear();
                  nextQ.reserve(faceQ.size());
                  nSize = faceQ.size();
                  for (size_t j = 0; j < nSize; j++) {
                       Face *currFace = faceQ[j];
                       currFace->getRelations12( neighs );
                       for (size_t i = 0; i < neighs.size(); i++) {
                            Face *nf = neighs[i];
                            if (!nf->getVisitBit()) nextQ.push_back(nf);
                       }
                  }

                  nSize = nextQ.size();
                  for (size_t i = 0; i < nSize; i++) {
                       Face *f = nextQ[i];
                       f->setAttribute("Layer", layerid);
                       f->setVisitBit(1);
                  }

                  layerid++;
                  faceQ = nextQ;
             }

             if (!relexist)
                  clearRelations(0, 2);

             cout << "Info: Number of layers in face front : " << layerid - 1 << endl;

             return layerid - 1;
        */
}
///////////////////////////////////////////////////////////////////////////////

void
JMeshTopology::set_cells_wavefront(int start_from_boundary)
{
    if( mesh == nullptr) return;

    mesh->buildRelations(2,3);

    size_t numCells = mesh->getSize(3);

    int nextid, layerid = 0;
    for (size_t i = 0; i < numCells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        if( cell->isActive() )
            cell->setAttribute("Layer", layerid);
    }

    JCellPtr cell;
    JCellSequence cellQ, neighs;
    size_t nSize;
    layerid = 1;
    if( start_from_boundary) {
        JFaceSequence boundfaces;
        this->getBoundary(boundfaces);
        nSize = boundfaces.size();
        for (size_t j = 0; j < nSize; j++) {
            JFacePtr face = boundfaces[j];
            JFace::getRelations(face, neighs);
            for (size_t i = 0; i < neighs.size(); i++) {
                cell = neighs[i];
                cell->setAttribute("Layer", layerid);
                cellQ.push_back(cell);
            }
        }
    } else {
        for (size_t i = 0; i < numCells; i++) {
            cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                cell->setAttribute("Layer", layerid);
                cellQ.push_back(cell);
                break;
            }
        }
    }

    JCellSet nextQ;

    layerid = 2;
    while (!cellQ.empty()) {
        nextQ.clear();
        nSize = cellQ.size();
        for (size_t j = 0; j < nSize; j++) {
            cell = cellQ[j];
            JCell::getRelations23(cell, neighs);
            for (size_t i = 0; i < neighs.size(); i++) {
                JCellPtr nc = neighs[i];
                if( nc->isActive()) {
                    nc->getAttribute("Layer", nextid);
                    if( nextid == 0) nextQ.insert(nc);
                }
            }
        }
        if( nextQ.empty() ) break;
        cellQ.clear();
        foreach_(cell, nextQ) {
            cell->setAttribute("Layer", layerid);
            cellQ.push_back(cell);
        }
        layerid++;
    }
}
///////////////////////////////////////////////////////////////////////////////

void
JMeshTopology::set_wavefronts( int start_from_boundary )
{
//  set_nodes_wavefront( start_from_boundary );
    set_faces_wavefront( start_from_boundary );
    set_cells_wavefront( start_from_boundary );
}

///////////////////////////////////////////////////////////////////////////////
bool JMeshTopology :: isManifold() const
{
    int topDim = getDimension();

    if( topDim == 2) {
        JEdgeSequence edges = getNonManifoldEdges();
        if( edges.empty() ) return 1;
        return 0;
    }

    if( topDim == 3) {
        JFaceSequence faces = getNonManifoldFaces();
        if( faces.empty() ) return 1;
        return 0;
    }
    return 0;
}

/*
int
Mesh::verify_front_ordering(int mentity)
{
     int l1, l2;
     size_t numfaces = getSize(2);
     if (mentity == 0) {
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               int nsize = face->getSize(0);
               for (int j = 0; j < nsize; j++) {
                    Node *v0 = face->getNodeAt( j + 0 );
                    Node *v1 = face->getNodeAt( j + 1 );
                    v0->getAttribute("Layer", l1 );
                    v1->getAttribute("Layer", l2 );
                    assert(l1 >= 0 && l2 >= 0);
                    if (abs(l2 - l1) > 1) return 1;
               }
          }
     }

     JFaceSequence neighs;

     if (mentity == 2) {
          int relexist2 = buildRelations(0, 2);
          for (size_t i = 0; i < numfaces; i++) {
               Face *face = getFaceAt(i);
               face->getRelations12( neighs );
               face->getAttribute("Layer", l1);
               int nf  = neighs.size();
               assert(l1 >= 0);
               for (int j = 0; j < nf; j++) {
                    neighs[j]->getAttribute("Layer", l2);
                    assert(l2 >= 0);
                    if (abs(l2 - l1) > 1) return 1;
               }
          }
          if (!relexist2) clearRelations(0, 2);
     }

     return 0;
}
*/

///////////////////////////////////////////////////////////////////////////////


/*
int Mesh::get_breadth_first_ordered_nodes(JNodeSequence &seq, Node *vstart, MeshFilter *filter)
{
     assert(vstart != NULL);

     seq.clear();

     int relexist0 = buildRelations(0, 0);

     size_t numnodes = getSize(0);

     if (numnodes == 0) return 1;

     for (size_t i = 0; i < numnodes; i++) {
          Node *v = getNodeAt(i);
          v->setVisitBit(0);
          v->setAttribute("Layer", 0);
     }

     if (vstart == 0) vstart = getNodeAt(0);

     seq.reserve(numnodes);

     list<Node*> vertexQ;
     vertexQ.push_back(vstart);
     JNodeSequence vneighs;

     int proceed = 1;
     int currlevel;
     while (!vertexQ.empty()) {
          Node *curr_vertex = vertexQ.front();
          vertexQ.pop_front();
          curr_vertex->getAttribute("Layer", currlevel);
          if (filter) {
               if (curr_vertex != vstart) proceed = filter->pass(curr_vertex);
          }
          if (!curr_vertex->getVisitBit()) {
               seq.push_back(curr_vertex);
               if (!proceed) break;
               curr_vertex->setVisitBit(1);
               curr_vertex->getRelations( vneighs );
               for (size_t i = 0; i < vneighs.size(); i++) {
                    if (!vneighs[i]->getVisitBit()) {
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

int Mesh::get_depth_first_ordered_nodes(JNodeSequence &seq, Node *vstart, MeshFilter *filter)
{
     int relexist0 = buildRelations(0, 0);

     size_t numnodes = getSize(0);
     list<Node*> vertexQ;
     for (size_t i = 0; i < numnodes; i++) {
          Node *v = getNodeAt(i);
          v->setVisitBit(0);
     }

     seq.clear();

     if (vstart == 0) vstart = getNodeAt(0);
     vertexQ.push_back(vstart);
     JNodeSequence vneighs;

     while (!vertexQ.empty()) {
          Node *curr_vertex = vertexQ.front();
          vertexQ.pop_front();
          if (!curr_vertex->getVisitBit()) {
               seq.push_back(curr_vertex);
               curr_vertex->setVisitBit(1);
               curr_vertex->getRelations( vneighs );
               for (size_t i = 0; i < vneighs.size(); i++) {
                    if (!vneighs[i]->getVisitBit())
                         vertexQ.push_front(vneighs[i]);
               }
          }
     }

     if (!relexist0)
          clearRelations(0, 0);

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

Jaal::JNodeSequence
Mesh::boundary_chain_nodes(Node *v0, Node *v1)
{
     JNodeSequence bndnodes;

     JFaceSequence neighs;
     Mesh::getRelations112(v0, v1, neighs);

     if (neighs.size() != 2) return bndnodes;

     vector<Edge> bndedges;

     bndedges.reserve(6);
     Edge sharededge(v0, v1);

     for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 4; j++) {
               Node *vf0 = neighs[i]->getNodeAt(j + 0);
               Node *vf1 = neighs[i]->getNodeAt(j + 1);
               Edge edge(vf0, vf1);
               if (!edge.isSameAs(sharededge))
                    bndedges.push_back(edge);
          }
     }

     int err = Mesh::make_chain(bndedges);
     if( err ) return bndnodes;

     bndnodes.reserve(6);
     bndnodes.push_back(bndedges[0].getNodeAt(0));
     bndnodes.push_back(bndedges[0].getNodeAt(1));

     size_t nSize = bndedges.size();
     for (size_t i = 1; i < nSize - 1; i++) {
          Node *v0 = bndedges[i].getNodeAt(0);
          Node *v1 = bndedges[i].getNodeAt(1);
          if (v0 == bndnodes[i]) {
               bndnodes.push_back(v1);
          } else if (v1 == bndnodes[i])
               bndnodes.push_back(v0);
          else {
               cout << "Error in bound edges : " << endl;
               bndnodes.clear();
               return bndnodes;
          }
     }
     return bndnodes;
}

/////////////////////////////////////////////////////////////////////////////
Node* MeshDualGrapher :: getDualNode( const Face *face)
{
    Point3D xyz;
    Node *dvertex = NULL;

    int numnodes = face->getSize(0);
    if( numnodes < 3) return dvertex;

    // For general Polygon, we return the average node.
    if( numnodes > 3 ) {
        dvertex = Node::newObject();
        face->getAvgXYZ(xyz);
        dvertex->setXYZCoords(xyz);
        return dvertex;
    }

    Node *v0 = face->getNodeAt(0);
    Node *v1 = face->getNodeAt(1);
    Node *v2 = face->getNodeAt(2);

    Point3D pa = v0->getXYZCoords();
    Point3D pb = v1->getXYZCoords();
    Point3D pc = v2->getXYZCoords();

    // For any obtuse triangle, the circumcenter lies outside the triangle.
    // but return the mid node of the longest edge...

    double angle;
    int pos = JMath::getMaxTriAngle(pa, pb, pc, angle);
    if( angle > 90.0) {
        dvertex = Node::newObject();
        face->getAvgXYZ(xyz);
        dvertex->setXYZCoords(xyz);
        return dvertex;
    }
    //
    // If the triangle is acute the circumcenter will lies within the
    // triangle: Bring any 3D triangle to a canonical oientation (z = 0)
    //
    double center[3], param[3];
    TriCircumCenter3D( &pa[0], &pb[0], &pc[0], center, param);

    unused_parameter(param);

    xyz[0] = center[0];
    xyz[1] = center[1];
    xyz[2] = center[2];

    dvertex = Node::newObject();
    dvertex->setXYZCoords(xyz);
    return dvertex;
}

////////////////////////////////////////////////////////////////////////////////

Mesh* MeshDualGrapher :: getGraph(Mesh *m)
{
    mesh = m;

    if( mesh == NULL ) return NULL;

    logger->setInfo("Building Dual graph");

    Mesh *dualGraph = Mesh::newObject();
    dualGraph->setName("DualGraph");

    int bmark;
    Point3D xyz;
    Node *dv0, *dv1, *vmid;
    JFaceSequence faceneighs;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim == 2 ) {
        mesh->buildRelations(1,2);
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            Face *face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                if( !face->hasAttribute("DualNode") ) {
                    Node *dvertex = getDualNode(face);
                    dualGraph->addObject(dvertex);
                    face->setAttribute("DualNode", dvertex);
                    dvertex->setAttribute("PrimalFace", face);
                }
            }
        }

        size_t numedges = mesh->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            Edge *edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                edge->getRelations(faceneighs);
                if( faceneighs.size() == 2 ) {
                    faceneighs[0]->getAttribute("DualNode", dv0);
                    faceneighs[1]->getAttribute("DualNode", dv1);
                    Edge *dedge = Edge::newObject(dv0,dv1);
                    if( midnodes ) {
                        vmid = Node::mid_node(edge->getNodeAt(0), edge->getNodeAt(1));
                        dedge->setAttribute("Steiner", vmid);
                    }
                    dualGraph->addObject(dedge);
                } else {
                    if( !edge->hasAttribute("DualNode") ) {
                        dv0 = Node::mid_node(edge->getNodeAt(0), edge->getNodeAt(1));
                        edge->setAttribute("DualNode", dv0);
                    }
                    edge->getAttribute("DualNode", dv0);
                    dualGraph->addObject(dv0);
                    bmark = max( 1, edge->getBoundaryMark() );
                    dv0->setBoundaryMark(bmark);
                    faceneighs[0]->getAttribute("DualNode", dv1);
                    Edge *dedge = Edge::newObject(dv0,dv1);
                    dualGraph->addObject(dedge);
                }
            }
        }
    }

    JCellSequence cellneighs;
    if( topDim == 3 ) {
        mesh->buildRelations(2,3);
        size_t numcells = mesh->getSize(3);
        for( size_t i = 0; i < numcells; i++) {
            Cell *cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                if(!cell->hasAttribute("DualNode") ) {
                    cell->getAvgXYZ(xyz);
                    Node *dvertex = Node::newObject();
                    dvertex->setXYZCoords(xyz);
                    dualGraph->addObject(dvertex);
                    cell->setAttribute("DualNode", dvertex);
                    dvertex->setAttribute("DualCell", cell);
                }
            }
        }

        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            Face *face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getRelations(cellneighs);
                if( cellneighs.size() == 2 ) {
                    cellneighs[0]->getAttribute("DualNode", dv0);
                    cellneighs[1]->getAttribute("DualNode", dv1);
                    Edge *dedge = Edge::newObject(dv0,dv1);
                    dualGraph->addObject(dedge);
                } else {
                    if(!face->hasAttribute("DualNode") ) {
                        dv0 = Node::newObject();
                        face->setAttribute("DualNode", dv0);
                    }
                    face->getAttribute("DualNode", dv0);
                    face->getAvgXYZ(xyz);
                    dv0->setXYZCoords(xyz);
                    dualGraph->addObject(dv0);
                    bmark = max( 1, face->getBoundaryMark() );
                    dv0->setBoundaryMark(bmark);
                    cellneighs[0]->getAttribute("DualNode", dv1);
                    Edge *dedge = Edge::newObject(dv0,dv1);
                    dualGraph->addObject(dedge);
                }
            }
        }
    }

    return dualGraph;
}
*/

void JMeshTopology :: setMinColor( const JNodePtr &vtx)
{
    /*
        if( !vtx->isActive() ) return;

        JNodeSequence neighs;
        vector<int>  assigned;

        vtx->getRelations(neighs);
        int numneighs = neighs.size();
        for( int i = 0; i < numneighs; i++)
             int err = neighs[i]->getAttribute("MinColor", cid);
             if( !err) assigned.push_back(cid);
        }
        sort( assigned.begin(), assigned.end() );

        int nsize;
        for( int i = 0; i < nsize; i++) {
             if( assigned[i] != i ) {
                 vtx->setAttribute("MinColor", i);
                 return;
             }
        }
        int maxval = *max_element( assigned.begin(), assigned.end() );
        vtx->setAttribute("MinColor", maxval+1);
    */
}


#ifdef CSV
int MeshUtil::rcm_ordering(SharedFaceMesh &facemesh, int entity)
{
    cout << " Need to relook " << endl;
    exit(0);

    /*
        MeshAdjacencyMatrix mesh_adjacency;
        JMath::SymmetricSparseMatrix<short> symmat = mesh_adjacency(facemesh, entity);

        JMath::GeneralSparseMatrix<short> genmat;
        genmat = JMath::general_sparse_matrix(symmat, 1);

        vector<MeshEntity::KeyType> perm = MeshUtil::rcm_ordering(genmat);

        if (entity == 0)
        {
            SharedJNodeSet nodeset = facemesh->getNodes(1);
            nodeset->resetGlobalIDs();
        }

        if (entity == 2)
        {
            facemesh->resetGlobalIDs(perm);
        }
    */
}


vector<int>
JMeshTopology::getStatistics(int entity, bool sorted)
{
    mesh->buildRelations(0, 2);

    assert(mesh->getAdjTable(0, 2));

    int numnodes = mesh->getSize(0);

    vector<int> degree;
    degree.reserve(numnodes);
    for (int i = 0; i < numnodes; i++) {
        Node *v = mesh->getNodeAt(i);
        if( v->isActive() )
            degree.push_back(v->getNumRelations(2) );
    }

    int mindegree = *min_element(degree );
    int maxdegree = *max_element(degree );

    cout << " Mesh Topological Quality : " << endl;

    cout << " ************************ " << endl;
    cout << " Degree         FaceCount " << endl;
    cout << " ************************ " << endl;

    for (int i = mindegree; i <= maxdegree; i++) {
        int ncount = 0;
        for (size_t j = 0; j < degree.size(); j++)
            if (degree[j] == i)
                ncount++;
        cout << setw(5) << i << setw(15) << ncount << endl;
    }

    return degree;
}
////////////////////////////////////////////////////////////////////////////////
bool
JMeshTopology ::isSimple()
{
    if( mesh == nullptr) return 0;

    logger->setInfo("Checking is the mesh is simple");

    /*
    int topoDim = getDimension();
        JCellSequence cellfaces;
        if( topoDim == 3) {
            mesh->buildRelations(2,3);
            size_t numfaces = mesh->getSize(2);
            for (size_t iface = 0; iface < numfaces; iface++) {
                JFacePtr face = mesh->getFaceAt(iface);
                if( face->isActive() ) {
                    JFace::getRelations(face, cellfaces);
                    if( cellfaces.size() > 2 ) return 0;
                }
            }
        }

        JEdgeSequence faceedges;
        if( topoDim == 2) {
            mesh->buildRelations(1,2);
            size_t numedges = mesh->getSize(1);
            for (size_t iedge = 0; iedge < numedges; iedge++) {
                JEdgePtr edge = mesh->getEdgeAt(iedge);
                if( edge->isActive() ) {
                    edge->getRelations(faceedges);
                    if( faceedges.size() > 2 ) return 0;
                }
            }
        }
    */

    return 1;
}

#endif
