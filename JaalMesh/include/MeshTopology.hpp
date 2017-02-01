#pragma once

#include "Mesh.hpp"
#include "basic_math.hpp"
#include "circumcenter.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#ifdef USE_IGL
#include <igl/matlab/matlabinterface.h>
#endif

namespace Jaal
{
enum class BettiMethod {
    SMITH_NORMAL_FORM, SPARSE_RANK, FULL_RANK, QR, EIGENANALYSIS
};

class JMeshComponents {
public:
    void setMesh(const JMeshPtr &m);
    void   searchComponents();
    int    getNumComponents() const ;
    JMeshPtr getComponent( int n) const;
private:
    JMeshPtr mesh;
    int      dim;
    vector<JMeshPtr> components;
    void searchFaceComponents();
    void searchCellComponents();
    JMeshPtr getFaceComponent(int id) const;
    JMeshPtr getCellComponent(int id) const;
};

////////////////////////////////////////////////////////////////////////////////

class JMeshTopology {
    static JLogger *logger;
public:

    static const int  PRIMAL_GRAPH = 0;
    static const int  DUAL_GRAPH   = 0;

    static void getEntitySet( const JEdgeSequence &edges, JNodeSequence &nodes);
    static void getEntitySet( const JFaceSequence &faces, JNodeSequence &nodes);
    static void getEntitySet( const JCellSequence &faces, JNodeSequence &nodes);

    static void getEntitySet( const JFaceSequence &faces, JEdgeSequence &edges);
    static void getEntitySet( const JCellSequence &cells, JEdgeSequence &edges);
    static void getEntitySet( const JCellSequence &cells, JFaceSequence &faces);
    static int  getEulerCharacteristic(const JFaceSequence &f);
    static int  getEulerCharacteristic(const JCellSequence &f);

    explicit JMeshTopology( const JMeshPtr &m )
    {
        mesh = m;
        boundary_known = 0;
        dualGraph = nullptr;
    }

    int getDimension() const
    {
        if( mesh->getSize(3) ) return 3;
        if( mesh->getSize(2) ) return 2;
        if( mesh->getSize(1) ) return 1;
        if( mesh->getSize(0) ) return 0;
        return 0;
    }

    bool isSameAs(const JMeshPtr &m) const;
    size_t countElementType(int e) const;

    // Collect faces from the cells
    int  collectFaces();
    // Collect all the active edges from the hash Container ( i.e) Vertex and put
    // into Mesh Object. Active edges are not removed from the hash container.
    //
    int  collectEdges();

    int get_nearest_neighbours( const JNodePtr &vtx, int k, JNodeSequence &seq);

    int getLexicographic(JEdgeSequence &e);
    int getLexicographic(JFaceSequence &e);
    int getLexicographic(JCellSequence &e);

    // Number of elements on the boundary..
    size_t getBoundarySize( int d )  const;

    //  Is the mesh elements homogeneous, if yes return the shape type..
    int  getElementsType( int dim) const;

    // A mesh is simple, when each edge is shared by at the most two faces which is
    // topological simple. A geometrically simple mesh will not have any crossing,
    // or overlapping edges, but this has not been checked right now.
    bool isSimple();

    // Is the mesh closed. Both edge and faces can be queried.
    bool isClosed();
    bool isDisk();

    int  getEulerCharacteristic();
    int  getBettiNumber(BettiMethod m );
    int  getNumComponents( bool stop_at_interface = 0);
    bool isConsistent() const;
    int  getConsistent();
    bool isManifold() const;

    /////////////////////////////////////////////////////////////////////////
    //  Check if the vertex is ideal.
    //  Complex     : Ideal degree
    //  Triangular    :  6
    //  Quadrilateral :  4
    //  Hexahedral    :  6
    ////////////////////////////////////////////////////////////////////////
    int isIdeal( JNodePtr vtx, int &idealsign) const;

    int    getIrregularNodes(JNodeSequence &seq, int where = JMeshEntity::INTERNAL_ENTITY);
    size_t getNumIrregularNodes();

    void getDoublets(JNodeSequence &nodes);
    void getDoublets(JEdgeSequence &nodes);

    // Returns the maximum degree of a node in the mesh...
    size_t getMaxNodeDegree();

    map<int,size_t>  getNodesDegreeDistribution( int where = JMeshEntity::ANY_ENTITY) const;

    JNodeSequence getNodeSubComplex();
    JEdgeSequence getEdgeSubComplex();
    JFaceSequence getFaceSubComplex();

    // Search the boundary. If the mesh is 2D, search for edges and if 3D
    //  search for faces. If an entity is on the surface then its lower entities
    //  are also marked on boundaries...
    int  searchBoundary();

    // Collect all the boundary nodes in the mesh. In case of 2D mesh. boundary
    // edges have two neighboring faces, and in 3D all the edges are part of
    // boundary faces are boundary edges..
    void getBoundary( JEdgeSequence &seq);
    void getBoundary( vector<JEdgeSequence> &seq);
    void getRim( const JNodePtr &vtx, JEdgeSequence &e);

    // If an edge is shared by more than two faces, it is called non-manifold.
    JEdgeSequence  getNonManifoldEdges()const;
    JFaceSequence  getNonManifoldFaces()const;

    // Collect all the boundary nodes in the mesh..
    void getBoundary( JNodeSequence &seq);
    void getBoundary( vector<JNodeSequence> &seq);

    // Collect all the boundary faces in the mesh (in 3D)
    void getBoundary( JFaceSequence &seq);
    void getBoundary( vector<JFaceSequence> &seq);

    void getSubmeshBoundary( const JFaceSequence &f, JEdgeSequence &e);
    void getSubmeshInternal( const JFaceSequence &f, JNodeSequence &);
    void getSubmeshInternal( const JFaceSequence &f, JEdgeSequence &e, JNodeSequence &);

    void getInternal( JEdgeSequence &e, JNodeSequence &v);

    void getInternal(JNodeSequence &seq);
    void getInternal(JEdgeSequence &seq);
    void getInternal(JFaceSequence &seq);
    void getInternal(JCellSequence &seq);
    void removeSubmesh(const JMeshPtr &m);

    JFacePtr  getCenter(const JFaceSequence &seq);

    // Collect all duplicated faces ...
    void getDuplicates( JEdgeSequence &seq);

    // Collect all duplicated faces ...
    void getDuplicates( JFaceSequence &seq);

    // Collect all duplicated faces ...
    void getDuplicates( JCellSequence &seq);

    // If the topology of the mesh is changed, make sure to setBoundaryStatus = 0
    // so that it is again determined.

    void setBoundaryKnown(bool v)
    {
        boundary_known = v;
    }

    // If the boundary nodes and faces are knowm, return the value = 1. otherwise 0

    bool isBoundaryKnown() const
    {
        return boundary_known;
    }

    void getOrphaned( JNodeSequence &v);
    void getOrphaned( JEdgeSequence &e);
    void getOrphaned( JFaceSequence &f);
    void removeOrphaned();

    void getDerivedRelations( const JNodePtr &vtx, JEdgeSequence &edges);

    JMeshPtr getSurfaceMesh();

    vector<int> getStatistics(int entity, bool sorted);

    // Backup the connectivity information...
    void backup();

    void fromBackup();
    void clearBackup();

    // Get the node connectivity as a stream ( Used for interface with other
    // software. example Mesquite )
    void getNodesArray( vector<size_t> &a, vector<int> &topo);

    // Although the nodes are appended, yet an user can request the nodes
    // in some typical sequence
    void getOrdered( const JNodePtr &v, JNodeSequence &seq, int ordered,
                               size_t maxSize = INT_MAX);

    void getOrdered( const JFacePtr &f, JFaceSequence &seq, int ordered,
                               size_t maxSize = INT_MAX );

    void getOrdered( const JCellPtr &c, JCellSequence &seq, int ordered,
                               size_t maxSize = INT_MAX );

    int  set_nodes_wavefront(const JNodePtr &src);
    void set_faces_wavefront(int start_from_boundary);
    void set_cells_wavefront(int start_from_boundary);
    void set_wavefronts(int start_from_boundary);

    vector<int> getMinNodesColor() const;
    void  setRCMOrdering();  // Reorder nodes in RCM order ...

    void reverseAll();
    int  setMinColors();
private:
    JMeshPtr mesh, dualGraph;
    bool boundary_known;
    void  search_3d_boundary();
    void  search_2d_boundary();

    void setMinColor( const JNodePtr &p);
};
};

