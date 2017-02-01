#pragma once

#include "MeshEntity.hpp"
#include <Eigen/Core>

////////////////////////////////////////////////////////////////////////////////

class JMeshGeometry;
typedef boost::shared_ptr<JMeshGeometry> JMeshGeometryPtr;

class JMeshTopology;
typedef boost::shared_ptr<JMeshTopology> JMeshTopologyPtr;

struct JMeshFilter {
    virtual bool passThrough(const JNodePtr &) const
    {
        cout << "Warning: passing the base " << endl;
        return 1;
    }

    virtual bool passThrough(const JFacePtr &) const
    {
        return 1;
    };
};
typedef boost::shared_ptr<JMeshFilter>  JMeshFilterPtr;

////////////////////////////////////////////////////////////////////////////////

template<class T>
inline void shrink_to_zero( vector<T> &v)
{
    v.clear();
    v.shrink_to_fit();
}

////////////////////////////////////////////////////////////////////////////////


class JMeshRelationManager {
    static JLogger *logger;
public:

    explicit JMeshRelationManager(JMesh *m)
    {
        mesh = m;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) adjTable[i][j] = 0;
        }
    }

    int setAdjTable(int i, int j, int val)
    {
        adjTable[i][j] = val;
        return 0;
    }

    bool getAdjTable(int i, int j) const
    {
        return adjTable[i][j];
    }

    void printAdjTable() const;

    void updateRelations();

    void buildRelations( const JEdgePtr  &e);
    void buildRelations( const JFacePtr  &f);
    void buildRelations( const JCellPtr  &c);

    int  removeRelations( const JNodePtr &v);
    int  removeRelations( const JEdgePtr &e);
    int  removeRelations( const JFacePtr &f);
    int  removeRelations( const JCellPtr &c);

    int buildRelations(int src, int dst);
    int clearRelations(int src, int dst);
private:
    JMesh  *mesh;   // Do not keep the shared ptr ...
    volatile char adjTable[4][4];

    int buildRelations00();
    int buildRelations01();
    int buildRelations02();
    int buildRelations03();

    int buildRelations12();
    int buildRelations13();

    int buildRelations21();
    int buildRelations23();

    int buildRelations31();
    int buildRelations32();
};

typedef boost::shared_ptr<JMeshRelationManager> JMeshRelationManagerPtr;

////////////////////////////////////////////////////////////////////////////////

class JMesh
{
    static size_t nCounter;
    friend class  JMeshRelationManager;
    friend class  JMeshTopology;
    friend class  JMeshGeometry;
    static JLogger *logger;
public:

    static  JLogger*  getLogger();
    typedef JNodeSequence::iterator  NodeIterator;
    typedef JEdgeSequence::iterator  EdgeIterator;
    typedef JFaceSequence::iterator  FaceIterator;
    typedef JCellSequence::iterator  CellIterator;

    typedef JNodeSequence::const_iterator  NodeConstIterator;
    typedef JEdgeSequence::const_iterator  EdgeConstIterator;
    typedef JFaceSequence::const_iterator  FaceConstIterator;
    typedef JCellSequence::const_iterator  CellConstIterator;

    static const int REMOVE  = 0;
    static const int ACTIVE  = 1;

    static const int RCM_ORDER = 0;           // Recursive-Cuthill Ordering.
    static const int BREADTH_FIRST_ORDER = 1;
    static const int DEPTH_FIRST_ORDER   = 2;

    static JMeshPtr newObject();
    static JMeshPtr newObject(JEdgeSequence &s);
    static JMeshPtr newObject(JFaceSequence &s);
    static JMeshPtr newObject(JCellSequence &s);

    static void generate_objects(size_t n, JNodeSequence &objects);
    static void generate_objects(size_t n, JEdgeSequence &objects);
    static void generate_objects(size_t n, JFaceSequence &objects, int type);
    static void generate_objects(size_t n, JCellSequence &objects, int type);

    JMesh()
    {
        initialize();
    }

    ~JMesh()
    {
        finalize();
    }

    JMeshPtr deepCopy();

    /*
        std::pair<NodeIterator,NodeIterator> node_iterators();
        std::pair<NodeConstIterator,NodeConstIterator> node_iterators() const;

        std::pair<EdgeIterator,EdgeIterator> edge_iterators();
        std::pair<EdgeConstIterator,EdgeConstIterator> edge_iterators() const;

        std::pair<FaceIterator,FaceIterator> face_iterators();
        std::pair<FaceConstIterator,FaceConstIterator> face_iterators() const;

        std::pair<CellIterator,CellIterator> cell_iterators();
        std::pair<CellConstIterator,CellConstIterator> cell_iterators() const;

        const JNodePtr &getNext( NodeConstIterator &it);
        const JEdgePtr &getNext( EdgeConstIterator &it);
        const JFacePtr &getNext( FaceConstIterator &it);
        const JCellPtr &getNext( CellConstIterator &it);
    */

    void setActiveBit( bool b) {
        activeBit = b;
    }

    bool isActive() const {
        return activeBit;
    }

    // Set one good name for the mesh...
    void setName(const string &s)
    {
        name = s;
    }

    const string getName() const
    {
        return name;
    }

    void setFileName( const string &s) {
        filename = s;
    }
    string getFileName() const {
        return filename;
    }

    void setStatus(int s)
    {
        status = s;
    }

    int  getStatus() const
    {
        return status;
    }

    size_t isEmpty()
    {
        return nodes.empty();
    }

    void progressiveSort(int etype);

    // Get the Current Container Size ...
    size_t getSize(int d) ;

    // Does this mesh have associated Model ?
    bool hasGeometricModel() const
    {
        return geomModel;
    }

    // Get information about the attribute associated with the mesh entity...
    int  getAttributeType(const string &s, int entity) const;

    //  What attributes are attached to the given entity..
    int  getAttributeNames(vector<string> &s, int entity) const;

    // How many attributes are attached to the given entity ( both sparse
    // and dense attributes) are identified..
    int  getNumAttributes(int e) const;

    // Identify if the given attribute is attached to any entity...
    bool hasAttribute( const string &name, int e) const;

    size_t getNumAttributes( const string &name, int e) const;

    // Assign the same Attribute value to each edge in the given sequence ...
    template<class T>
    void setAttribute( JEdgeSequence &eseq, const string &name, const T &val)
    {
        for(JEdgePtr edge: eseq) edge->setAttribute(name, val);
    }

    template<class T>
    void getEntities( const string &name, const T &val, JNodeSequence &vseq)
    {
        vseq.clear();
        T thisval;
        for(JNodePtr vtx : nodes) {
            if( vtx->isActive() ) {
                int err = vtx->getAttribute(name, thisval);
                if( err == 0 && thisval == val) vseq.push_back(vtx);
            }
        }
    }

    // Collect all the edges with the given attribute value ...
    template<class T>
    void getEntities( const string &name, const T &val, JEdgeSequence &eseq)
    {
        eseq.clear();
        T thisval;
        for(JEdgePtr edge : edges) {
            if( edge->isActive() ) {
                int err = edge->getAttribute(name, thisval);
                if( err == 0 && thisval == val) eseq.push_back(edge);
            }
        }
    }

    void getEntities( const string &name, JEdgeSequence &eseq)
    {
        eseq.clear();
        for(JEdgePtr edge : edges) {
            if( edge->isActive() && edge->hasAttribute(name) ) eseq.push_back(edge);
        }
    }

    // Collect all the edges within th range of attribute values (both inclusive)
    template<class T>
    void getEntities( const string &name, const T &lowVal, const T &highVal,
                      JEdgeSequence &eseq)
    {
        eseq.clear();
        T thisval;
        for(JEdgePtr edge : edges) {
            int err = edge->getAttribute(name, thisval);
            if( err == 0 && thisval >= lowVal && thisval <= highVal)
                eseq.push_back(edge);
        }
    }

    void setGeometry( const JMeshGeometryPtr &g) { geometry = g; }
    const JMeshGeometryPtr getGeometry() const { return  geometry; }

    void setTopology( const JMeshTopologyPtr &t) { topology = t; }
    const JMeshTopologyPtr getTopology() const { return topology; }

    void setVisitBits(int e, bool val);

    // Number of elements that are active now.
    size_t getActiveSize( int d )  const;

    // What is the highest dimension of the entity ?
    int getDimension() const
    {
        if( !cells.empty() ) return 3;
        if( !faces.empty() ) return 2;
        if( !edges.empty() ) return 1;
        if( !nodes.empty() ) return 0;
        return -1;
    }

    // What is the capacity of each mesh entity container...
    size_t getCapacity(int d)
    {
        if (d == 0) return nodes.capacity();
        if (d == 1) return edges.capacity();
        if (d == 2) return faces.capacity();
        if (d == 3) return cells.capacity();
        return 0;
    }

    // You may reserve the memory, if you know in advance...
    void reserve(size_t nSize, int entity)
    {
        if( nSize < 1 ) return;
        if (entity == 0) nodes.reserve(nSize);
        if (entity == 1) edges.reserve(nSize);
        if (entity == 2) faces.reserve(nSize);
        if (entity == 3) cells.reserve(nSize);
    }

    // Change the Adjacency table ....
    int setAdjTable(int i, int j, int val)
    {
        return relManager->setAdjTable(i,j,val);
    }

    bool getAdjTable(int i, int j) const
    {
        return relManager->getAdjTable(i,j);
    }

    // Which mesh entities, you would like to persist in the memory...
    // Nodes are always persist. For 2D simplices, faces are always
    // known, only the edges are derived. For 3D simplicaes, both the
    // mesh edges and faces be derived: therefore. only cells and
    // nodes are alaways represnt.
    void keepEntity( int e, bool v)
    {
        if( e < 0 || e > 3 ) return;
        keep_mesh_entity[e] = v;
    }

    // Query, if you could if the user requested for persistence of
    // some meshenttity "e";
    //
    bool hasEntity(int e ) const
    {
        return keep_mesh_entity[e];
    }

    void addObject( const JMeshPtr &m )
    {
        const JNodeSequence &nlist = m->getNodes();
        const JEdgeSequence &elist = m->getEdges();
        const JFaceSequence &flist = m->getFaces();
        const JCellSequence &clist = m->getCells();

        addObjects(nlist);
        addObjects(elist);
        addObjects(flist);
        addObjects(clist);
    }

    // Appends a node in the mesh.
    void addObject(const JNodePtr &vtx, bool embedded = 0)
    {
        vtx->setStatus(JMeshEntity::ACTIVE);
        if( embedded )
            embedded_nodes.push_back(vtx);
        else
            nodes.push_back(vtx);
    }

    // Appends a bulk of nodes in the mesh.
    void addObjects(const JNodeSequence &vnodes, bool embedded = 0);
    int  remove(const JNodePtr &v);

    // Get the node at the specified position in the container.
    const JNodePtr &getNodeAt(size_t id) const
    {
        assert(id < nodes.size());
        return nodes[id];
    }

    int sort(const JEntityCompare *compare, int e = 0);

    int getPosOf(const JNodePtr &v, size_t &pos) const;

    // Get All the nodes in the mesh after lazy prunning ( Garbage collection ).
    JNodeSequence  getNodes() const ;
    JNodeSequence  getEmbeddedNodes() const ;
    JNodeSequence  getNodes(const string &attribname) const;
    JNodeSequence  getElementsNodes() const ;

    // Check if the node is present in the mesh..
    bool contains(const JNodePtr &v) const
    {
        if (find(nodes.begin(), nodes.end(), v) == nodes.end()) return 0;
        return 1;
    }
    // Swap the position of v1 with v2
    int swap( const JNodePtr &v1, const JNodePtr &v2);

    // Clear all the nodes, but do not delete them. Perhaps they are
    // associated with the other mesh...
    void deleteNodes(int type = JMeshEntity::ANY_ENTITY);

    // Delete the attribute ( if present) from the node.
    void deleteNodeAttribute(const string &s);

    ///////////////////////////////////////////////////////////////////////////
    // Edge functions ....
    ///////////////////////////////////////////////////////////////////////////

    int  addObject(const JEdgePtr &e, bool embedded = 0 );
    void addObjects(const JEdgeSequence  &seq, bool embedded = 0 );

    JEdgePtr getEdgeAt(size_t id) const
    {
        assert(id < edges.size());
        return edges[id];
    }

    int getPosOf(const JEdgePtr &e, size_t &pos) const;

    // Return all the edges of the primal mesh ...
    JEdgeSequence getEdges();
    JEdgeSequence getEdges(const string &s);
    JEdgeSequence getEmbeddedEdges() const ;
    JEdgeSequence getElementsEdges() const ;

    int  remove( const JEdgePtr &v);
    int  swap( const JEdgePtr &e1, const JEdgePtr &e2);

    void  deleteEdges( int type = JMeshEntity::ANY_ENTITY);
    void  deleteEdgeAttribute(const string &s);

    // Check if the edge is present in the mesh..
    bool contains(const JEdgePtr &e) const
    {
        if (find(edges.begin(), edges.end(), e) == edges.end()) return 0;
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Face Functions ....
    ///////////////////////////////////////////////////////////////////////////

    // Add a face in the mesh. No duplication is checked..
    int   addObject(const JFacePtr &f, bool embedded = 0);
    void  addObjects(const JFaceSequence &vfaces);

    // Get the face at the specified position in the container.
    JFacePtr getFaceAt(size_t id) const
    {
        if( id >= faces.size() ) {
            cout << "Warning: Requesting invalid face " << id << endl;
            return nullptr;
        }
        return faces[id];
    }

    int getPosOf(const JFacePtr &f, size_t &pos) const;

    JFaceSequence getFaces() const;
    JFaceSequence getFaces( const string &s) const;
    JFaceSequence getEmbeddedFaces() const ;
    JFaceSequence getElementsFaces() const ;

    size_t getNumFaces(int facetype);

    bool contains(const JFacePtr &f) const
    {
        if (find(faces.begin(), faces.end(), f) == faces.end()) return 0;
        return 1;
    }

    int  swap( const JFacePtr &f1, const JFacePtr &f2);
    int  remove(const JFacePtr &f);
    void deleteFaces(int type = JMeshEntity::ANY_ENTITY);
    void deleteFaceAttribute(const string &s);

    void invert_faces();

    //////////////////////////////////////////////////////////////////////////
    // Cell functions ....
    //////////////////////////////////////////////////////////////////////////

    const JCellPtr &getCellAt(size_t id) const
    {
        assert(id < cells.size());
        return cells[id];
    }

    int getPosOf(const JCellPtr &c, size_t &pos) const;

    int  addObject(const JCellPtr &c);
    void addObjects( const JCellSequence &vcells);

    size_t getNumCells(int facetype);

    JCellSequence getCells() const;
    JCellSequence getCells(const string &s) const;
    JCellSequence getEmbeddedCells() const ;

    int  remove( const JCellPtr &c);
    int  swap( const JCellPtr &c1, const JCellPtr &c2);

    void deleteCells();
    void deleteCellAttribute(const string &s);

    // Check if the cell is present in the mesh..
    bool contains(const JCellPtr &c) const
    {
        if (find(cells.begin(), cells.end(), c) == cells.end()) return 0;
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////

    int buildRelations(int src, int dst)
    {
        return relManager->buildRelations(src,dst);
    }

    void updateRelations()
    {
        relManager->updateRelations();
    }

    // clean specified entity-entity relations.
    int clearRelations(int src, int dst)
    {
        return relManager->clearRelations(src,dst);
    }
    int  clearAllRelations()
    {
        for( int src = 0; src < 4; src++)
            for( int dst = 0; dst < 4; dst++)
                return relManager->clearRelations(src,dst);
        return 0;
    }

    int swapNodePositions( size_t node1, size_t node2 )
    {
        if( node1 >= nodes.size() || node2 >= nodes.size() ) return 1;
        if( nodes[node1]->getID() != node1 ) return 2;
        if( nodes[node2]->getID() != node2 ) return 2;
        std::swap( nodes[node1], nodes[node2] );
        nodes[node1]->setID( node1 );
        nodes[node2]->setID( node2 );
        return 0;
    }

    // Check if the lazy garbage collection is performed..
    bool isPruned() const;
    // Get rid of entities which are marked "removed" from the mesh.
    // a la Lazy garbage collection.
    //
    void pruneNodes();
    void pruneEdges();
    void pruneFaces();
    void pruneCells();
    void pruneAll();

    void clearAll()
    {
        cells.clear();
        faces.clear();
        edges.clear();
        nodes.clear();
    }

    // Update the mesh after the insertion, deletion, and modification..
    void update();

    void deleteVolumeMesh();
    void deleteSurfaceMesh();

    // Empty every thing in the Mesh, and also deallocate all the objects.
    void deleteAll();

    // Renumber mesh entities starting from index = 0
    void enumerate(int etype);

//  void progressiveSort(int etype);

    template<class T>
    int setAttribute(const string &s, const T &val)
    {
        if( attribManager == nullptr ) attribManager.reset(new JAttributeManager);
        return attribManager->setAttribute(s,val);
    }

    template<class T>
    int getAttribute(const string &s, T &val) const
    {
        if( attribManager )
            return attribManager->getAttribute(s,val);
        return 1;
    }

    int getNumAttributes() const
    {
        int ncount = 0;
        if( attribManager ) {
            ncount = attribManager->getNumAttributes();
        }
        return ncount;
    }

    int getAttributeNames( vector<string> &names) const
    {
        names.clear();
        if( attribManager )
            return attribManager->getAttributeNames( names );
        return 0;
    }

    int hasAttribute(const string &s) const
    {
        if( attribManager )
            return attribManager->hasAttribute(s);
        return 0;
    }

    void  deleteAttribute(const string &s)
    {
        if( attribManager ) attribManager->deleteAttribute(s);
    }

    const Array4I getFVector() const {
        Array4I fvec;
        fvec[0] =  getActiveSize(0);
        fvec[1] =  getActiveSize(1);
        fvec[2] =  getActiveSize(2);
        fvec[3] =  getActiveSize(3);
        return fvec;
    }

    // Save the mesh and all its attributes ( in the simple format ).
    int saveAs(const string &s);

private:
    volatile int topoChanged[4], geomChanged;
    int status;
    bool activeBit;

    string name, filename;

    bool  geomModel;
    int   keep_mesh_entity[4];

    boost::shared_ptr<JMeshGeometry>  geometry;
    boost::shared_ptr<JMeshTopology>  topology;

    boost::scoped_ptr<JAttributeManager>  attribManager;
    boost::scoped_ptr<JMeshRelationManager>  relManager;

    // Contains all the mesh entities.
    JNodeSequence nodes, embedded_nodes;
    JEdgeSequence edges, embedded_edges;
    JFaceSequence faces, embedded_faces;
    JCellSequence cells;

    void initialize();
    void finalize();

};

////////////////////////////////////////////////////////////////////////////////

void set_tfi_coords(int i, int j, int nx, int ny, JNodeSequence &qnodes);
///////////////////////////////////////////////////////////////////////////////

struct ParametricSurface {
    static JMeshPtr getCrossCapSurface(int nx, int ny);
    static JMeshPtr getBoySurface(int nx, int ny);
    static JMeshPtr getAperyBoySurface(int nx, int ny);
};

struct JMeshIO
{
    static JMeshPtr readFile( const string &f);
    static int  saveAs( const JMeshPtr &m, const string &f);
};

} // namespace Jaal

#include "MeshGeometry.hpp"
#include "MeshTopology.hpp"
#include "MeshImporter.hpp"
#include "MeshExporter.hpp"
