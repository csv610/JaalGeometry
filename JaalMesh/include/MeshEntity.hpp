#pragma once

#include <iostream>
#include <cassert>
#include <fstream>
#include <math.h>
#include <string>
#include <values.h>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <queue>
#include <algorithm>
#include <unordered_map>

#ifdef HAVE_VERDICT
#include <verdict.h>
#endif

#include <boost/any.hpp>
#include <boost/cast.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>
#include <boost/range/numeric.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/fusion/include/push_back.hpp>
#include <boost/algorithm/cxx11/is_sorted.hpp>

#define foreach_    BOOST_FOREACH
#define foreach_r_  BOOST_REVERSE_FOREACH
typedef boost::tokenizer<> StringTokenizer;

#include "basic_math.hpp"
#include "tfiblend.hpp"
#include "BoundingBox.hpp"
#include "circumcenter.hpp"
#include "Geometry.hpp"
#include "Logger.hpp"
#include "Attributes.hpp"

# include <omp.h>

struct MutexType {
    MutexType()
    {
        omp_init_lock(&thislock);
    }
    ~MutexType()
    {
        omp_destroy_lock(&thislock);
    }
    void lock()
    {
        omp_set_lock(&thislock);
    }
    void unlock()
    {
        omp_unset_lock(&thislock);
    }

    MutexType(const MutexType& )
    {
        omp_init_lock(&thislock);
    }
    MutexType& operator= (const MutexType& )
    {
        return *this;
    }
public:
    omp_lock_t thislock;
};
/* A dummy mutex that doesn't actually exclude anything,
 * but as there is no parallelism either, no worries. */
/*
struct MutexType {
    void lock() {}
    void unlock() {}
};
*/

/* An exception-safe scoped lock-keeper. */
struct ScopedLock {
    explicit ScopedLock(MutexType& m) : mut(m), locked(true)
    {
        mut.lock();
    }
    ~ScopedLock()
    {
        unlock();
    }

    void unlock()
    {
        if(!locked) return;
        locked=false;
        mut.unlock();
    }

    void lockAgain()
    {
        if(locked) return;
        mut.lock();
        locked=true;
    }
private:
    MutexType& mut;
    bool locked;
private: // prevent copying the scoped lock.
    void operator=(const ScopedLock&);
    ScopedLock(const ScopedLock&);
};

#include <unordered_set>
#include <unordered_map>

using namespace std;
using namespace boost;

template<typename T>
void unused_parameter( T const &) {}

template <typename T>
inline string get_pod_name( T )
{
    return "Unknown";
}

template<>
inline string get_pod_name( int )         {
    return "int";
}

template<>
inline string get_pod_name( short int)    {
    return "sint";
}

template<>
inline string get_pod_name( float )       {
    return "float";
}

template<>
inline string get_pod_name( double )      {
    return "double";
}

template<>
inline string get_pod_name( char  )       {
    return "char";
}

template<>
inline string get_pod_name( unsigned char ) {
    return "uchar";
}

template<>
inline string get_pod_name( string ) {
    return "string";
}

namespace Jaal {
#define GEOMETRIC_ERROR   2
#define TOPOLOGICAL_ERROR 3


class JMeshEntity;
class JNode;
class JEdge;
class JFace;
class JCell;

class JTriangle;
class JQuadrilateral;
class JPolygon;

class JTetrahedron;
class JHexahedron;
class JPolyhedron;
class JTriangularPrism;

class JMesh;

typedef boost::shared_ptr<JMeshEntity> JMeshEntityPtr;
typedef boost::shared_ptr<JNode> JNodePtr;
typedef boost::shared_ptr<JEdge> JEdgePtr;
typedef boost::shared_ptr<JFace> JFacePtr;
typedef boost::shared_ptr<JCell> JCellPtr;

typedef boost::shared_ptr<JTriangle> JTrianglePtr;
typedef boost::shared_ptr<JPolygon>  JPolygonPtr;
typedef boost::shared_ptr<JQuadrilateral> JQuadrilateralPtr;
typedef boost::shared_ptr<JTetrahedron> JTetrahedronPtr;
typedef boost::shared_ptr<JHexahedron> JHexahedronPtr;
typedef boost::shared_ptr<JPolyhedron> JPolyhedronPtr;
typedef boost::shared_ptr<JTriangularPrism> JTriangularPrismPtr;
typedef boost::shared_ptr<JMesh>    JMeshPtr;

typedef std::vector<JNodePtr> JNodeSequence;
typedef std::vector<JEdgePtr> JEdgeSequence;
typedef std::vector<JFacePtr> JFaceSequence;
typedef std::vector<JCellPtr> JCellSequence;

typedef std::list<JNodePtr> JNodeList;
typedef std::list<JEdgePtr> JEdgeList;
typedef std::list<JFacePtr> JFaceList;
typedef std::list<JCellPtr> JCellList;

typedef std::set<JNodePtr> JNodeSet;
typedef std::set<JEdgePtr> JEdgeSet;
typedef std::set<JFacePtr> JFaceSet;
typedef std::set<JCellPtr> JCellSet;

typedef std::pair<JFacePtr, JFacePtr> JFacePair;

#define TOPOLOGICAL_DISTANCE  0
#define EUCLIDEAN_DISTANCE    1
#define EUCLIDEAN_DISTANCE2   2
#define CITY_BLOCK_DISTANCE   3

////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

template <class T>
int split_stl_vector(const std::vector<T> &a, size_t pos, std::vector<T> &b, std::vector<T> &c)
{
    size_t nSize = a.size();
    if (pos >= nSize ) return 1;

    b.resize(pos);

    size_t index = 0;
    for (size_t i = 0; i < pos; i++)
        b[index++] = a[i];

    c.resize(nSize - pos + 1);
    index = 0;

    for (size_t i = pos - 1; i < nSize; i++)
        c[index++] = a[i];

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

template <class T>
int split_stl_vector(const std::deque<T> &a, size_t pos, std::deque<T> &b, std::deque<T> &c)
{
    size_t nSize = a.size();

    if (pos >= nSize) return 1;

    b.resize(pos);

    size_t index = 0;
    for (size_t i = 0; i < pos; i++)
        b[index++] = a[i];

    c.resize(nSize - pos + 1);
    index = 0;
    for (size_t i = pos - 1; i < nSize; i++)
        c[index++] = a[i];

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

class JRelationManager {
    static JLogger* logger;
public:

    static const int SORTED_RELATIONS   = 0;
    static const int ORDERED_RELATIONS  = 1;

    ~JRelationManager()
    {
        clearAll();
    }

    // Add a vertex in the relation ( duplication not allowed )
    void addRelation(const JNodePtr &vertex, int type = SORTED_RELATIONS);

    // Add an edge in the relation ( duplication not allowed )
    void addRelation(const JEdgePtr &edge, int type = SORTED_RELATIONS);

    // Add a face in the relation ( duplication not allowed )
    void addRelation(const JFacePtr &face, int type = SORTED_RELATIONS);

    // Add a cell in the relation ( duplication not allowed )
    void addRelation(const JCellPtr &cell, int type = SORTED_RELATIONS);

    // Get the number of relations of given rank stored.
    int  getNumRelations( int e ) const;

    // Check if the given vertex is in the node relation list ...
    bool hasRelation(const JNodePtr &vertex) const;

    // Check if the given edge is in the edge relation list ...
    bool hasRelation(const JEdgePtr &edge) const;

    // Check if the given edge is in the face relation list ...
    bool hasRelation(const JFacePtr &face) const;

    // Check if the given edge is in the cell relation list ...
    bool hasRelation(const JCellPtr &cell) const;

    // Edge ...
    JEdgePtr getEdgeOf( const JEdgePtr &dummy)  const;

    // Edge ...
    JEdgePtr getEdgeOf( const JNodePtr &v0, const JNodePtr &v1)  const;

    // Triangle Face ...
    JFacePtr getFaceOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2) const;

    // Quadrilateral Face ...
    JFacePtr getFaceOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3)  const;

    JFacePtr getFaceOf( const JFacePtr &f) const;
    // General Polygon ....
    JFacePtr getFaceOf( const JNodeSequence &n) const;

    // Tetrahedron cell...
    JCellPtr getCellOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3)  const;
    // Hexahedron cell...
    JCellPtr getCellOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3,
                        const JNodePtr &v4, const JNodePtr &v5, const JNodePtr &v6, const JNodePtr &v7)  const;

    JCellPtr getCellOf( const JCellPtr &c) const;
    // General Polyhedra ...
    JCellPtr getCellOf( const JNodeSequence &n) const;

    // removes the given vertex from the relation list
    void removeRelation(const JNodePtr &vertex);

    // removes the given edge from the relation list
    void removeRelation(const JEdgePtr &edge);

    // removes the given face from the relation list
    void removeRelation(const JFacePtr &face);

    // removes the given cell from the relation list
    void removeRelation(const JCellPtr &cell);

    void prune( int e );
    void pruneAll();

    // Clear relations of the given rank...
    void clearRelations(int t);

    void clearAll()
    {
        clearRelations(0);
        clearRelations(1);
        clearRelations(2);
        clearRelations(3);
    }

    int getRelations( JNodeSequence &seq, bool cyclic_ordered = 0) const;
    int getRelations( JEdgeSequence &seq, bool cyclic_ordered = 0) const;
    int getRelations( JFaceSequence &seq, bool cyclic_ordered = 0) const;
    int getRelations( JCellSequence &seq, bool cyclic_ordered = 0) const;

    void setRelations( const JNodeSequence &seq)
    {
        relations0 = seq;
    }
    void setRelations( const JEdgeSequence &seq)
    {
        relations1 = seq;
    }
    void setRelations( const JFaceSequence &seq)
    {
        relations2 = seq;
    }
    void setRelations( const JCellSequence &seq)
    {
        relations3 = seq;
    }

    JNodeSequence relations0; // nodes
    JEdgeSequence relations1; // edges
    JFaceSequence relations2; // faces
    JCellSequence relations3; // cells
};

typedef boost::shared_ptr<JRelationManager> JRelationManagerPtr;

/////////////////////////////////////////////////////////////////////////////////////

class JMeshEntity {
public:
    static const int INTERNAL_ENTITY      = 0;
    static const int BOUNDARY_ENTITY      = 1;
    static const int INTERFACE_ENTITY     = 2;
    static const int FEATURE_ENTITY       = 3;
    static const int CONSTRAINED_ENTITY   = 4;
    static const int UNATTACHED_ENTITY    = 5;
    static const int EMBEDDED_ENTITY      = 6;
    static const int ANY_ENTITY           = 7;

    typedef size_t idtype;

    static const int REMOVE    = 0;
    static const int ACTIVE    = 1;
    static const int INACTIVE  = 2;

    static void init_random_number();

    JMeshEntity()
    {
        visitBit     = 0; // Default: Not yet visited.
        statusMark   = ALIVE; // Default: Active< not removable.
        attribManager.reset(new JAttributeManager);
        relationManager.reset(new JRelationManager);
    }

    virtual ~JMeshEntity()
    {
        if( attribManager )   attribManager->clearAll();
        if( relationManager ) relationManager->clearAll();
    }

    void reset();

    virtual JNodePtr getHashNode() const = 0;

    virtual int getDimension() const = 0;

    void setVisitBit(bool r)
    {
        visitBit = r;
    }

    bool getVisitBit() const
    {
        return visitBit;
    }

    int getStatus() const
    {
        return statusMark;
    }

    bool isRemoved() const
    {
        if (statusMark == REMOVE) return 1;
        return 0;
    }

    virtual bool isActive() const
    {
        if (statusMark == ACTIVE) return 1;
        return 0;
    }

    bool isBoundary() const
    {
        return this->hasAttribute("Boundary");
    }

    void setID(idtype id)
    {
        gid = id;
    }

    const idtype &getID() const
    {
        return gid;
    }

    template<class T>
    int setAttribute(const string &s, const T &val)
    {
        return attribManager->setAttribute(s,val);
    }

    template<class T>
    int getAttribute(const string &s, T &val) const
    {
        return attribManager->getAttribute(s,val);
    }

    int getNumAttributes() const
    {
        return attribManager->getNumAttributes();
    }

    int getAttributeNames( vector<string> &names) const
    {
        names.clear();
        return attribManager->getAttributeNames( names );
    }

    int hasAttribute(const string &s) const
    {
        return attribManager->hasAttribute(s);
    }

    void  deleteAttribute(const string &s)
    {
        attribManager->deleteAttribute(s);
    }

    void clearAllRelations()
    {
        relationManager->clearRelations(0);
        relationManager->clearRelations(1);
        relationManager->clearRelations(2);
        relationManager->clearRelations(3);
    }

    void clearRelations(int e)
    {
        relationManager->clearRelations(e);
    }

    void addRelation(const JNodePtr &vertex, int type = JRelationManager::SORTED_RELATIONS)
    {
        relationManager->addRelation(vertex, type);
    }

    void setRelations( const JNodeSequence &seq)
    {
        relationManager->setRelations(seq);
    }

    void addRelation(const JEdgePtr  &edge, int type =  JRelationManager::SORTED_RELATIONS)
    {
        relationManager->addRelation(edge, type);
    }

    void setRelations( const JEdgeSequence &seq)
    {
        relationManager->setRelations(seq);
    }

    void addRelation(const JFacePtr &face, int type = JRelationManager::SORTED_RELATIONS)
    {
        relationManager->addRelation(face, type);
    }

    void setRelations( const JFaceSequence &seq)
    {
        relationManager->setRelations(seq);
    }

    void addRelation(const JCellPtr &cell, int type = JRelationManager::SORTED_RELATIONS)
    {
        relationManager->addRelation(cell, type);
    }

    void setRelations( const JCellSequence &seq)
    {
        relationManager->setRelations(seq);
    }

    int  getNumRelations( int e ) const
    {
        return relationManager->getNumRelations(e);
    }

    int getNumHigherRelations() const
    {
        int myDimension = getDimension();
        int nCount = 0;
        for( int i = myDimension+1; i < 4; i++)
            nCount += getNumRelations(i);
        return nCount;
    }

    int getNumLowerRelations() const
    {
        int myDimension = getDimension();
        int nCount = 0;
        for( int i = 0; i < myDimension-1; i++)
            nCount += getNumRelations(i);
        return nCount;
    }

    bool hasRelation(const JNodePtr &vertex) const
    {
        return relationManager->hasRelation( vertex );
    }

    bool hasRelation(const JEdgePtr &edge) const
    {
        return relationManager->hasRelation( edge );
    }

    bool hasRelation(const JFacePtr &face) const
    {
        return relationManager->hasRelation( face );
    }

    bool hasRelation(const JCellPtr &cell) const
    {
        return relationManager->hasRelation( cell );
    }

    void removeRelation(const JNodePtr &vertex)
    {
        relationManager->removeRelation( vertex );
    }

    void removeRelation(const JEdgePtr &edge)
    {
        relationManager->removeRelation( edge );
    }

    void removeRelation(const JFacePtr &face)
    {
        relationManager->removeRelation( face );
    }

    void removeRelation(const JCellPtr &cell)
    {
        relationManager->removeRelation( cell );
    }

    // Since you can not have the same name as static member, I have put "_" after the name...

    int getRelations_( JNodeSequence &seq, bool cyclic_ordered = 0) const
    {
        seq.clear();
        return relationManager->getRelations(seq, cyclic_ordered);
    }

    int getRelations_( JEdgeSequence &seq, bool cyclic_ordered = 0) const
    {
        seq.clear();
        return relationManager->getRelations(seq, cyclic_ordered);
    }

    int getRelations_( JFaceSequence &seq, bool cyclic_ordered = 0) const
    {
        seq.clear();
        return relationManager->getRelations(seq, cyclic_ordered);
    }

    int getRelations_( JCellSequence &seq, bool cyclic_ordered = 0) const
    {
        seq.clear();
        return relationManager->getRelations(seq, cyclic_ordered);
    }

    void pruneRelations()
    {
        relationManager->pruneAll();
    }

    virtual void setStatus(int a)
    {
        statusMark = a;
    }

    virtual int unitTest()
    {
        return 0;
    }

protected:
    idtype gid;
    volatile bool visitBit;
    volatile short int statusMark;
    boost::scoped_ptr<JRelationManager>  relationManager;
    boost::scoped_ptr<JAttributeManager> attribManager;
};

///////////////////////////////////////////////////////////////////////////////

struct EntityRemovedPred {

    bool operator() (const JMeshEntityPtr entity) const
    {
        if (entity) return entity->isRemoved();
        return 0;
    }
};

///////////////////////////////////////////////////////////////////////////////

struct BaseNode : public JMeshEntity {
    BaseNode()
    {
    }

    ~BaseNode()
    {
    }

    virtual int getDimension() const
    {
        return 0;
    }

    void attach( const JEdgePtr &edge)
    {
        if( hashRel == nullptr) hashRel.reset(new JRelationManager);
        hashRel->addRelation(edge);
    }

    void attach( const JFacePtr &face)
    {
        if( hashRel == nullptr) hashRel.reset(new JRelationManager);
        hashRel->addRelation(face);
    }

    void attach( const JCellPtr &cell)
    {
        if( hashRel == nullptr) hashRel.reset(new JRelationManager);
        hashRel->addRelation(cell);
    }

    void detach( const JEdgePtr &edge)
    {
        if( hashRel == nullptr) return;
        hashRel->removeRelation(edge);
    }

    void detach( const JFacePtr &face)
    {
        if( hashRel == nullptr) return;
        hashRel->removeRelation(face);
    }

    void detach( const JCellPtr &cell)
    {
        if( hashRel == nullptr) return;
        hashRel->removeRelation(cell);
    }

    void detachAll()
    {
        if( hashRel == nullptr) return;
        hashRel->clearAll();
    }

    int getNumEntitiesHashed() const
    {
        if( hashRel == nullptr) return 0;
        int nCount = 0;
        nCount += hashRel->getNumRelations(0);
        nCount += hashRel->getNumRelations(1);
        nCount += hashRel->getNumRelations(2);
        nCount += hashRel->getNumRelations(3);
        return nCount;
    }

    void setStatus( int s )
    {
        if( s == JMeshEntity::REMOVE) detachAll();
        JMeshEntity::setStatus(s);
    }

    JEdgePtr getEdgeOf( const JNodePtr &v0, const JNodePtr &v1)  const
    {
        if( hashRel) return hashRel->getEdgeOf(v0,v1);
        JEdgePtr nullPtr;
        return nullPtr;
    }

    // General Polygon...
    JFacePtr getFaceOf( const JNodeSequence &n) const
    {
        if( hashRel) return hashRel->getFaceOf(n);
        JFacePtr nullPtr;
        return nullPtr;
    }

    // Triangle Face ...
    JFacePtr getFaceOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2)  const
    {
        if( hashRel) return hashRel->getFaceOf(v0,v1,v2);
        JFacePtr nullPtr;
        return nullPtr;
    }

    // Quadrilateral Face ...
    JFacePtr getFaceOf( const JNodePtr &v0, const JNodePtr &v1,
                        const JNodePtr &v2, const JNodePtr &v3)  const
    {
        if( hashRel) return hashRel->getFaceOf(v0,v1,v2,v3);
        JFacePtr nullPtr;
        return nullPtr;

    }

    // Tetrahedron Face ...
    JCellPtr getCellOf( const JNodePtr &v0, const JNodePtr &v1,
                        const JNodePtr &v2, const JNodePtr &v3)  const
    {

        if( hashRel) return hashRel->getCellOf(v0,v1,v2,v3);
        JCellPtr nullPtr;
        return nullPtr;
    }

    // Hexahedron Face ...
    JCellPtr getCellOf( const JNodePtr &v0, const JNodePtr &v1,
                        const JNodePtr &v2, const JNodePtr &v3,
                        const JNodePtr &v4, const JNodePtr &v5,
                        const JNodePtr &v6, const JNodePtr &v7)  const
    {
        if( hashRel) return hashRel->getCellOf(v0,v1,v2,v3,v4,v5,v6,v7);
        JCellPtr nullPtr;
        return nullPtr;
    }

    // General polyhedra...
    JCellPtr getCellOf( const JNodeSequence &n) const
    {
        if( hashRel) return hashRel->getCellOf(n);
        JCellPtr nullPtr;
        return nullPtr;
    }

// Some vertices keep the information of other higher order entities.
// This is different from relationship, as relationship are created
// on-demand and can be destroyed any time, but not the information
// contain with the vertex..

    boost::scoped_ptr<JRelationManager> hashRel;
};

///////////////////////////////////////////////////////////////////////////////
struct EntityFeatureSet {
public:
    void addAttribute( const string &s)
    {
        attributes.insert( s );
    }
    void deleteAttribute( const string &s)
    {
        attributes.erase( s );
    }
    bool hasAttribute( const string &s);

    vector<string>  getFeatureAttributes() const
    {
        vector<string> vs;
        if( attributes.empty() ) return vs;
        vs.resize( attributes.size() );
        std::copy( attributes.begin(), attributes.end(), vs.begin() );
        return vs;
    }
    set<string> attributes;
};

///////////////////////////////////////////////////////////////////////////////

class JNode: public BaseNode {
    static size_t GlobalID;
    static std::map<string,string> attribInfo;
public:
    static size_t NumObjectsCreated;

    static JNodePtr down_cast( const JMeshEntityPtr &m);
    static JNodePtr newObject();
    static JNodeSequence newObjects(size_t n);

    static void  getRelations( const JNodePtr &e, JNodeSequence &seq);
    static void  getRelations( const JNodePtr &e, JEdgeSequence &seq);
    static void  getRelations( const JNodePtr &e, JFaceSequence &seq);
    static void  getRelations( const JNodePtr &e, JCellSequence &seq);

    static int  registerAttribute(const string &s, const string &type);
    static string getAttributeTypeName( const string &s);

    JNode() {}

    ~JNode() {}

    void setXYZCoords( double x, double y, double z)
    {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
    }

    void setXYZCoords(const Point3D &p)
    {
        xyz[0] = p[0];
        xyz[1] = p[1];
        xyz[2] = p[2];
    }

    JNodePtr getHashNode() const
    {
        JNodePtr nullPtr;
        return nullPtr;
    }

    const Point3D &getXYZCoords() const
    {
        return xyz;
    }

    Point2D getXYCoords() const
    {
        Point2D p2d;
        p2d[0] = xyz[0];
        p2d[1] = xyz[1];
        return p2d;
    }

    Point3D &getXYZCoords()
    {
        return xyz;
    }

    JNodePtr getClone() const;

    bool isFeature() const;
    bool isInternal() const;
    int    get_ideal_face_degree( int n ) const;
private:
    Point3D xyz;
};

struct JNodeGeometry {
    static double  getAngleDefect( const JNodePtr &vertex);
    static double  getArea( const JNodePtr &vtx);
    static double  getFeatureLength(const JNodePtr &vertex);
    static double  getLength(const JNodePtr &v0, const JNodePtr &v1);
    static double  getLength2(const JNodePtr &v0, const JNodePtr &v1);
    static Point3D  getMidPoint(const JNodePtr &v0, const JNodePtr &v1, double alpha = 0.5);
    static JNodePtr getMidNode(const JNodePtr &v0, const JNodePtr &v1,  double alpha = 0.5);
    static double   getOrientation( const Point3D &p0, const Point3D &p1, const Point3D &qpoint);
    static Point3D getCenter( const JNodeSequence &nodes);
    static JBoundingBox   getBoundingBox(const JNodeSequence &nodes);
    static double getSpanAngleAt( const JNodePtr &v, int measure);
    static double getSolidAngleAt( const JNodePtr &v);
};


///////////////////////////////////////////////////////////////////////////////

inline JNodePtr JNode::down_cast( const JMeshEntityPtr &entity)
{
    JNodePtr vtx = boost::dynamic_pointer_cast<JNode>(entity);
    assert(vtx);
    return vtx;
}

struct GeomProperty {
};

struct NodeGeomProperty : public GeomProperty {
    double getSpanAngle() const;
};

struct EdgeGeomProperty : public GeomProperty {
    static double getLength( const JEdgePtr &f );
};

struct FaceGeomProperty : public GeomProperty {
    static double getArea(const JFacePtr &f );
};

struct CellGeomProperty : public GeomProperty {
    static double getVolume( const JCellPtr &f );
};

////////////////////////////////////////////////////////////////////////////////

struct JEntityCompare
{
    virtual ~JEntityCompare() {}

    virtual bool isLess( const  JNodePtr &v1, const JNodePtr &v2)  const {
        return v1 < v2;
    }

    virtual bool isLess( const  JEdgePtr &e1, const JEdgePtr &e2)  const {
        return e1 < e2;
    }

    virtual bool isLess( const  JFacePtr &e1, const JFacePtr &e2)  const {
        return e1 < e2;
    }

    virtual bool isLess( const  JCellPtr &e1, const JCellPtr &e2)  const {
        return e1 < e2;
    }
};

////////////////////////////////////////////////////////////////////////////////

struct JEntityIDCompare : public JEntityCompare {

    ~JEntityIDCompare() {}

    template<class T>
    bool operator() (const T &e1, const T &e2) const
    {
        if( e1->isActive() && e2->isActive() ) {
            size_t d1 = e1->getID();
            size_t d2 = e2->getID();
            return d1 < d2;
        }
        return 0;
    }


    bool isLess(const JNodePtr &e1, const JNodePtr &e2) const
    {
        if( e1->isActive() && e2->isActive() ) {
            size_t d1 = e1->getID();
            size_t d2 = e2->getID();
            return d1 < d2;
        }
        return 0;
    }
};

///////////////////////////////////////////////////////////////////////////////
inline Point3D JNodeGeometry::getCenter( const JNodeSequence &nodes)
{
    Point3D p3d;
    p3d[0] = 0.0;
    p3d[1] = 0.0;
    p3d[2] = 0.0;
    for( size_t i = 0; i < nodes.size(); i++) {
        const Point3D &xyz = nodes[i]->getXYZCoords();
        p3d[0] += xyz[0];
        p3d[1] += xyz[1];
        p3d[2] += xyz[2];
    }
    p3d[0] /= ( double) nodes.size();
    p3d[1] /= ( double) nodes.size();
    p3d[2] /= ( double) nodes.size();
    return p3d;
}
///////////////////////////////////////////////////////////////////////////////

class JSimplex : public JMeshEntity
{
public:

    static double linear_interpolation(const vector<double> &x, const vector<double> &w);

    static JEdgePtr getEdgeOf( const JNodePtr &v0, const JNodePtr &v1, bool build = 0);
    // Triangular face ...
    static JFacePtr getFaceOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, bool build = 0);
    // Quadrangular face ..
    static JFacePtr getFaceOf( const JNodePtr &v0, const JNodePtr &v1,
                               const JNodePtr &v2, const JNodePtr &v3, bool build = 0);
    // Tetrahedral face..
    static JCellPtr getCellOf( const JNodePtr &v0, const JNodePtr &v1,
                               const JNodePtr &v2, const JNodePtr &v3, bool build = 0);
    // Hexahedral face ...
    static JCellPtr getCellOf( const JNodePtr &v0, const JNodePtr &v1,
                               const JNodePtr &v2, const JNodePtr &v3,
                               const JNodePtr &v4, const JNodePtr &v5,
                               const JNodePtr &v6, const JNodePtr &v7, bool build = 0);

    // Presently very simple, but shouldn't matter.
    JNodePtr getHashNode() const
    {
        if( !nodes.empty() ) return *min_element(nodes);
        JNodePtr nullPtr;
        return nullPtr;
    }

// Simplex() {}
    ~JSimplex() {}

    bool hasActiveNodes() const
    {
        int nSize = nodes.size();
        for (int i = 0; i < nSize; i++)
            if( !nodes[i]->isActive() ) return 0;
        return 1;
    }

    virtual int getDimension() const = 0;

    virtual int build_lower_entities(int dim) = 0;
    virtual int remove_unattached_lower_entities() = 0;

    int setNodes(const JNodeSequence &v)
    {
        int nSize = v.size();
        bool err = 0;
        for (int i = 0; i < nSize; i++) {
            for (int j = i+1; j < nSize; j++)
                if(v[i] == v[j]) err = 1;
        }

        if( !err ) {
            nodes = v;
            setStatus( JMeshEntity :: ACTIVE );
        }
        return err;
    }

    bool hasUniqueNodes() const;

    bool hasNode(const JNodePtr  &vertex) const
    {
        if (find(nodes.begin(), nodes.end(), vertex) != nodes.end()) return 1;
        return 0;
    }

    bool hasNodeID( size_t id) const
    {
        for( const JNodePtr &vtx: nodes)
            if( vtx->getID() == id ) return 1;
        return 0;
    }

    int getPosOf(const JNodePtr &vertex) const
    {
        int nSize = nodes.size();
        for (int i = 0; i < nSize; i++)
            if (nodes[i] == vertex) return i;

        return -1;
    }

    virtual int getPosOf(const JEdgePtr &) const
    {
        return -1;
    }

    virtual int getPosOf(const JFacePtr &) const
    {
        return -1;
    }

    int replace(const JNodePtr &oldvertex, const JNodePtr &newvertex)
    {

        int  nSize = nodes.size();
        for (int i = 0; i < nSize; i++) {
            if (nodes[i] == oldvertex) {
                nodes[i] = newvertex;
                return 0;
            }
        }
        return 1;
    }

    const JNodeSequence &getNodes() const
    {
        return nodes;
    }

    void getNodes( vector<size_t> &v) const
    {
        int nSize = nodes.size();
        v.resize( nSize );
        for(int i = 0; i < nSize; i++)
            v[i] = nodes[i]->getID();
    }

    bool hasBoundaryNode() const
    {
        int nSize = nodes.size();
        for (int i = 0; i < nSize; i++)
            if (nodes[i]->isBoundary()) return 1;
        return 0;
    }

    void getXYZCoords( vector<Point3D> &p)
    {
        p.clear();
        int nsize = nodes.size();
        if( nsize < 2) return;
        p.resize(nsize);
        for( int i = 0; i < nsize; i++)
            p[i] = nodes[i]->getXYZCoords();
    }

    int getAvgXYZ( Point3D &p) const;
    int getAvgUV( Point2D &p) const;

    void print() const
    {
        for( size_t i = 0; i < nodes.size() ; i++)
            cout << nodes[i]->getID() << " ";
        cout << endl;
    }

    JBoundingBox getBoundingBox() const;

    JNodeSequence nodes;
};

///////////////////////////////////////////////////////////////////////////////

class JEdge : public JSimplex {
    static std::map<string,string> attribInfo;
public:
    static size_t NumObjectsCreated;
    static size_t AnyObjectModified;

    static EntityFeatureSet *featureSet;
    static const JEdgePtr down_cast( const JMeshEntityPtr &m);

    static int   registerAttribute(const string &s, const string &type);
    static string getAttributeTypeName( const string &s);

    static JEdgePtr newObject();
    static JEdgePtr newObject( const JNodeSequence &vn);
    static JEdgePtr newObject( const JNodePtr &v0, const JNodePtr &v1);
    static JEdgeSequence newObjects(size_t n);

    static void  getRelations( const JEdgePtr &e, JFaceSequence &seq);
    static void  getRelations( const JEdgePtr &e, JCellSequence &seq);

    JEdge()
    {
        nodes.resize(2);
    }

    JEdge(const JNodePtr &n1, const JNodePtr &n2)
    {
        setNodes(n1, n2);
    }

    ~JEdge() {}

    virtual JEdgePtr getClone() const
    {
        JEdgePtr newedge = newObject(nodes[0], nodes[1]);
        return newedge;
    }

    int getDimension() const
    {
        return 1;
    }

    int build_lower_entities(int )
    {
        cout <<"Warning: There are no derived lower entities for an edge" << endl;
        return 1;
    }
    int remove_unattached_lower_entities();

    void setNodes(const JNodeSequence &v)
    {
        assert(v.size() == 2);
        setNodes(v[0], v[1]);
        setStatus( JMeshEntity :: ACTIVE );
    }

    void setNodes(const JNodePtr &v1, const JNodePtr &v2)
    {
        assert(v1 != v2);
        nodes.resize(2);
        nodes[0] = v1;
        nodes[1] = v2;
        setStatus( JMeshEntity :: ACTIVE );
    }

    const JNodePtr &getNodeAt(int id) const
    {
        return nodes[id];
    }

    JNodePtr getOtherNode(const JNodePtr &v) const
    {
        if( nodes[0] == v) return nodes[1];
        if( nodes[1] == v) return nodes[0];
        JNodePtr nullPtr;
        return nullPtr;
    }

    bool isSameAs(const JEdgePtr &rhs) const
    {
        if ((nodes[0] == rhs->nodes[0]) && (nodes[1] == rhs->nodes[1])) return 1;
        if ((nodes[0] == rhs->nodes[1]) && (nodes[1] == rhs->nodes[0])) return 1;
        return 0;
    }

    bool hasSameNodesID( const JEdgePtr &rhs) const
    {
        if ((nodes[0]->getID() == rhs->nodes[0]->getID()) && (nodes[1]->getID() == rhs->nodes[1]->getID())) return 1;
        if ((nodes[0]->getID() == rhs->nodes[1]->getID()) && (nodes[1]->getID() == rhs->nodes[0]->getID())) return 1;
        return 0;
    }

    bool isFeature() const;

    int getOrientation( const JEdgePtr &rhs) const
    {
        if( rhs->getNodeAt(0) == nodes[0] && rhs->getNodeAt(1) == nodes[1] ) return  1;
        if( rhs->getNodeAt(0) == nodes[1] && rhs->getNodeAt(1) == nodes[0] ) return -1;
        return 0;
    }

    int getOrientation( const JNodePtr &v0, const JNodePtr &v1) const
    {
        if( nodes[0] == v0 && nodes[1] == v1 ) return  1;
        if( nodes[0] == v1 && nodes[1] == v0 ) return -1;
        return 0;
    }

    int getMidNodes( int n, JNodeSequence &newnodes) const
    {
        if( n < 1 ) return 1;
        newnodes.resize(n);
        double dl = 1.0/(double)(n+1);
        for( int i = 0; i < n; i++)
            newnodes[i] = JNodeGeometry::getMidNode( nodes[0], nodes[1], (i+1)*dl );
        return 0;
    }

    void reverse()
    {
        std::swap( nodes[0], nodes[1] );
    }

    bool isManifold() const {
        if( this->isBoundary()  && this->getNumRelations(2) == 1) return 1;
        if( !this->isBoundary() && this->getNumRelations(2) == 2) return 1;
        return 0;
    }
};

/*
inline JEdgePtr JEdge::down_cast( const JMeshEntityPtr &entity)
{
    JEdgePtr edge = dynamic_pointer_cast<JEdge>(entity);
    assert(edge);
    return edge;
}
*/

////////////////////////////////////////////////////////////////////////////////

class JEdgeGeometry {
public:
    static bool isWithin( const double *p0, const double *p1, const double *pq);

    // Check, if the two segments intersects, It does not returns the intersecting point coordiantes.
    static int intersectPredicate2d( const double *p0, const double *p1,
                                     const double *p2, const double *p4);

    static int intersectPredicate2d( const JEdgePtr &edge, const JEdgePtr &edge2);

    // If we know two segments intersect, we can query for their intersection point. The return value is in
    // the range of 0.0, 1.0). if the range is outside, then the lines do not interesects.
    static int intersectPos2d( const double *p1, const double *p2,
                               const double *p3, const double *p4,
                               double *xy = nullptr, double *uv = nullptr);

    static int  intersectPos2d( const JEdgePtr &edge, const JEdgePtr &edge2,
                                double *xy = nullptr, double *pu = nullptr);

    static double getU( const double *p0, const double *p1, const double *xy);
    static void   smooth( JEdgeSequence &seq, int numIterations = 1);

    static double getLength( const JEdgePtr &edge)
    {
        const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
        return JMath::length( p0, p1);
    }

    static double getLength2( const JEdgePtr &edge)
    {
        const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
        return JMath::length2( p0, p1);
    }

    static double getLength( const JEdgeSequence &edges)
    {
        double sum = 0.0;
        for( const JEdgePtr &edge : edges)  sum += JEdgeGeometry::getLength(edge);
        return sum;
    }
    static int  makeUniform(const JEdgeSequence &e);
    static int  makeUniformRefine(const JEdgeSequence &e, int numnodes);
 
    static double   getCreaseAngle( const JEdgePtr &p);
    static int      getOrientation(const JEdgeSequence &edge);
    static Point3D  getMidPoint( const JEdgePtr &p, double alpha = 0.5);
    static JNodePtr getMidNode( const JEdgePtr &p, double alpha = 0.5);

    static void generateLinearNodes( const JNodePtr &v0, const JNodePtr &v1, int n, JNodeSequence &s);
    static void generateLinearNodes( const JEdgePtr &edge, int n, JNodeSequence &s);
    static void generateLinearNodes( const JEdgePtr &edge, int n, const JMeshPtr &mesh);
};

struct JEdgeTopology
{
    static bool isTopologicalSimple( const JEdgeSequence &edges);
    static int  getChain( JEdgeSequence &edges);
    static int  getChain( JEdgeSequence &edges, const JNodePtr &start);
    static bool isChain(const JEdgeSequence &edges);
    static int  getChain( JEdgeSequence &edges, const JEdgePtr &start);
    static int  isClosed(const JEdgeSequence &edges);
    static int  isCloseable(const JEdgeSequence &edges);
    static int  getChainNodes(const JEdgeSequence &e, JNodeSequence &nodes);
    static void reverse(JEdgeSequence &edges);
    static bool lexiCompare(const JEdgePtr &edge1, const JEdgePtr &edge2);
    static int  getLoops( const JEdgeSequence &seq, vector<JEdgeSequence> &l);
};


/////////////////////////////////////////////////////////////////////////////////
class JFace : public JSimplex {
    friend class JFaceGeometry;
    static std::map<string,string> attribInfo;
public:
    static size_t NumObjectsCreated;
    static size_t AnyObjectModified;
    static EntityFeatureSet *featureSet;

    static const int POLYGON       = 20;
    static const int TRIANGLE      = 23;
    static const int QUADRILATERAL = 24;

    static JFacePtr down_cast( const JMeshEntityPtr &m);
    static int    registerAttribute(const string &s, const string &type);
    static string getAttributeTypeName( const string &s);

    static JFacePtr getProduct(int type );
    static JFacePtr newObject( const JNodeSequence &seq);
    static JFacePtr getRandomPolygon( int n );


    /////////////////////////////////////////////////////////////////////////////
    // Use modified Heron formula for find out the area of a triangle
    /////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////////
    // Calculate Area of Quadrilateral:
    // Example : Convex Case
    //           Coorindates  (0.0,0.0)  ( 1.0, 1.0), (2.0, 0.0), (1,0, -1.0)
    //           Result :  2*( 0.5*2.0*1.0 ) = 2.0;
    //           Concave case:
    //           Coodinates: ( 0.0, 0.0), 10, 2), (2.0, 0.0), (10, -2)
    //           Result : 2 *( 0.5*2.0* 2.0) = 4.0
    /////////////////////////////////////////////////////////////////////////////

    static JNodePtr getDiagonalNode(const JFacePtr &quad, const JNodePtr &n1);
    static void getOppositeNodes(const JFacePtr &quad, JNodePtr &n1, JNodePtr &n2,
                                 JNodePtr &n3, JNodePtr &n4);
    static int check_on_boundary(const JFacePtr tri);

    static void  getRelations02( const JFacePtr &e, JFaceSequence &seq);
    static void  getRelations12( const JFacePtr &e, JFaceSequence &seq);

    static void  getRelations( const JFacePtr &e, JCellSequence &seq);
    static void  getSharedEntities( const JFacePtr &t1, const JFacePtr &t2, JNodeSequence &nodes);
    static void  getSharedEntities( const JFacePtr &t1, const JFacePtr &t2, JEdgeSequence &edges);

    static bool lexiCompare(const JFacePtr f1, const JFacePtr f2);

//  JFace() { }

    virtual ~JFace() {}

    virtual JFacePtr getClone() const = 0;

    JFacePtr explode( double alpha = 0.99) const;

    virtual int  getTypeID()   const = 0;
    virtual string getTypeName() const = 0;

    int getDimension() const
    {
        return 2;
    }

    int getType() const
    {
        if (nodes.size() == 3) return JFace::TRIANGLE;
        if (nodes.size() == 4) return JFace::QUADRILATERAL;
        return JFace::POLYGON;
    }

    int build_lower_entities( int dim ) ;
    int remove_unattached_lower_entities();

    void reverse()
    {
        std::reverse(nodes.begin(), nodes.end());
    }

    bool isSameAs(const JFacePtr rhs) const
    {
        size_t nSize = rhs->nodes.size();
        if (nSize != nodes.size() ) return 0;

        for (size_t i = 0; i < nSize; i++)
            if (!rhs->hasNode(nodes[i])) return 0;

        return 1;
    }

    bool hasSameNodesID(const JFacePtr rhs) const
    {
        size_t nSize = rhs->nodes.size();
        if (nSize != nodes.size() ) return 0;

        for (size_t i = 0; i < nSize; i++)
            if (!rhs->hasNodeID(nodes[i]->getID())) return 0;
        return 1;
    }

    int getOrientation(const JEdgePtr &rhs) const;
    int getOrientation( const JNodePtr &n0, const JNodePtr &n1, const JNodePtr &n2) const;
    int getOrientation( const JNodePtr &n0, const JNodePtr &n1, const JNodePtr &n2,
                        const JNodePtr &n3) const;

    int getSize(int etype) const
    {
        if (etype == 0) return nodes.size();
        if (etype == 1) return nodes.size();
        return 0;
    }

    const JNodePtr &getNodeAt(int  id) const
    {
        return nodes[id%nodes.size()];
    }

    int getPosOf( const JNodePtr &v) const
    {
        return JSimplex::getPosOf(v);
    }

    int getPosOf( const JEdgePtr &edge) const
    {
        int nSize = getSize(0);
        for( int i = 0; i < nSize; i++)
            if( getEdgeAt(i) == edge ) return i;
        return -1;
    }

    JEdgePtr getEdgeAt( int pos ) const
    {
        const JNodePtr &v0 = getNodeAt(pos);
        const JNodePtr &v1 = getNodeAt(pos+1);
        return JSimplex::getEdgeOf(v0,v1,1);
    }

    bool hasEdgeAt( int pos) const
    {
        const JNodePtr &v0 = getNodeAt(pos);
        const JNodePtr &v1 = getNodeAt(pos+1);
        JEdgePtr e    =  JSimplex::getEdgeOf(v0,v1,0);
        if( e) return 1;
        return 0;
    }

    JEdgePtr getEdgeOf( const JNodePtr &v0, const JNodePtr &v1);

    JEdgeSequence getEdges() const
    {
        JEdgeSequence localedges;
        int  nn = getSize(0);
        if( nn > 2) {
        localedges.resize(nn);
        for( int i = 0; i < nn; i++)
            localedges[i] = getEdgeAt(i);
        }
        return localedges;
    }

    bool has_all_bound_nodes() const;
    bool has_boundary_edge() const;

    int triangulate(JFaceSequence &newfaces, JNodeSequence &newnodes, int n = 2)const;

    virtual int refine( JNodeSequence &newnodes, JFaceSequence &newfaces)
    {
        newnodes.clear();
        newfaces.clear();
        return 1;
    }

    int setStartNode( const JNodePtr &p);

    bool isManifold() const {
        if( this->isBoundary()  && this->getNumRelations(3) == 1) return 1;
        if( !this->isBoundary() && this->getNumRelations(3) == 2) return 1;
        return 0;
    }
};

class JFaceGeometry {
public:

    static int is_3_sided_convex_loop_quad_meshable(const int *s, int *sdiv);
    static int is_5_sided_convex_loop_quad_meshable(const int *s, int *sdiv);

    static int getDualPosition(const JFacePtr f, Point3D &xyz);

//  static void getCentroid( const double *x, const double *y, int n, double *c);

    // Centroid of Triangle element ...
    static JNodePtr  getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2);
    static JNodePtr  getCentroid(const JFacePtr &f);
    static void   getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, Point3D &pd);
    static void   getCentroid(const JFacePtr &f, Point3D &p);

    // Centroid of Quad  element ...
    static JNodePtr getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2,
                                const JNodePtr &v3);
    static void    getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2,
                               const JNodePtr &v3, Point3D &p3d);

    static JNodePtr getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2,
                                const JNodePtr &v3, const JNodePtr &v4);
    static void    getCentroid(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2,
                               const JNodePtr &v3, const JNodePtr &v4, Point3D &p3d);

    static Vec3D  getNormal( const JFacePtr &f);
    static double getMaxAngle( const JFacePtr &f, int measure);
    static double getMinAngle( const JFacePtr &f, int measure);
    static double getAngleAt( const JFacePtr &f, int pos, int measure);
    static double getAngleAt( const JFacePtr &f, const JNodePtr &v, int measure);
    static void   getAngles( const JFacePtr &f, vector<double> &angles, int measure);

    // Are two faces congruent, i.e. do they have same interior angles...
    static int  areCongruent( const JFacePtr &f1, const JFacePtr &f2, double epsilon = 1.0E-03);

    // is the given face congrunent to the given angles...
    static bool   isSimple(const JFacePtr &f);
    static int    isCongruent( const JFacePtr& f1, const vector<double> &angles, double epsilon = 1.0E-03);
    static bool   isConvex(const JFacePtr &f);
    static int    reflexAngleAt( const JFacePtr &f);
    static double getAspectRatio( const JFacePtr &f);

    static double getArea( const JFacePtr &f);
    static double getSignedArea( const JFacePtr &f);

    /*
        static double getArea2D( const double *x, const double *y, int n);
        static double getArea2D( const vector<Point2D> &p);
    */
    static double getArea(   const JFaceSequence &fs);

    static bool   intersectPredicate2d( const JFacePtr &f1, const JFacePtr &f2);

//  static int    getBoundedSide(const double *p0, const double *p1, const double *p2, const double *ptest);
    static int    getBoundedSide(const JFacePtr &f1, const Point3D &p);
    static int    getOrientation2D( const JFacePtr &f);
    static bool   isInverted( const JFacePtr &f);
    static bool   isPlanar( const JFacePtr &f);
};

///////////////////////////////////////////////////////////////////////////////

inline JFacePtr JFace::down_cast( const JMeshEntityPtr &entity)
{
    JFacePtr face = dynamic_pointer_cast<JFace>(entity);
    return face;
}

////////////////////////////////////////////////////////////////////////////////

class JTriangle : public JFace {
public:

    static JTrianglePtr newObject();
    static JTrianglePtr newObject( const JNodeSequence &s);
    static JTrianglePtr newObject( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2);

    static JTrianglePtr getCanonical( double len = 1.0);
    static JTrianglePtr down_cast( const JFacePtr &f);

    static JNodePtr getOppositeNode(const JFacePtr &f, const JNodePtr &n1, const JNodePtr &n2);
    static JEdgePtr getOppositeEdge(const JFacePtr &f, const JNodePtr &e);

    static JFaceSequence newObjects(size_t n);

    ~JTriangle() {}

    JFacePtr getClone() const
    {
        JFacePtr newtri = newObject();
        newtri->setNodes( nodes );
        return newtri;
    }

    string getTypeName() const
    {
        return "Triangle";
    }

    int  getTypeID() const
    {
        return JFace::TRIANGLE;
    }

    int setNodes( const JNodeSequence &seq)
    {
        assert(seq.size() == 3 );
        int err = JFace::setNodes(seq);
        return err;
    }

    int setNodes(const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3)
    {
        JNodeSequence seq(3);
        seq[0] = (v1);
        seq[1] = (v2);
        seq[2] = (v3);
        JFace::setNodes(seq);
        return 0;
    }

//   int getPosOf( const JNodePtr &v0, const JNodePtr &v1) const;

};

class JTriGeometry {
public:

    static double getArea(const Point3D &p0, const Point3D &p1, const Point3D &p2);
    static double getSignedArea(const double *p0, const double *p1, const double *p2);
    static double getSignedArea(const Point3D &p0, const Point3D &p1, const Point3D &p2);

    static void   getMaxAngle( const JFacePtr &f, double &angle, int &pos);
    static double getMaxAngle(const JFacePtr &f);

    static Vec3D   getNormal(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2);
    static Vec3D   getNormal(const Point3D &p0, const Point3D &p1, const Point3D &v2);
    static Point3D getConformal(const Point3D &pa, const Point3D &pb, const Point3D &pc,
                                const Point3D &angles);
    static Point3D getIdealPosition(const JFacePtr &f, const JNodePtr &v, double tipangle);
    static double  getDistortion(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2);

    static int  getBaryCoords( const double *p0, const double *p1, const double *p2,
                               const double *xy, double *uv);
    static int  getUVCoords( const double *p0, const double *p1, const double *p2,
                             const double *xy, double *uv);
    static int  getXYZCoordinates( const JFacePtr &f, const Point3D &baryCoord, Point3D &xyz);
    static int  getXYZCoordinates( const double *p0,  const double *p1, const double *p2,
                                   const double *uvw, double *xyz);

    static bool isInside(const JFacePtr &c, const Point3D &p,  bool include_boundary = 0);
    static bool isOutside(const JFacePtr &c, const Point3D &p, bool include_boundary = 0);
    static bool isOnBoundary(const JFacePtr &c, const Point3D &p);
    static int  getBaryCoordinates( const JFacePtr &f, const Point3D &xyz, Point3D &baryCoords);
    static Point3D  getRandomPoint( const JFacePtr &p);
};

////////////////////////////////////////////////////////////////////////////////

class JQuadrilateral : public JFace
{
public:
    static JQuadrilateralPtr newObject();
    static JQuadrilateralPtr newObject( const JNodeSequence &s);
    static JQuadrilateralPtr newObject( const JNodePtr &v0, const JNodePtr &v1,
                                        const JNodePtr &v2, const JNodePtr &v3);
    static JQuadrilateralPtr newObject(const JFacePtr &t1, const JFacePtr &t2); // From two triangles.

    static JQuadrilateralPtr getCanonical(double l = 1.0);
    static JQuadrilateralPtr getCanonical(const Point3D &center, double l = 1.0);
    static JQuadrilateralPtr getCanonical(double lx, double ly);
    static JQuadrilateralPtr getCanonical(const Point3D &normal, double dist, double len);

    static JNodePtr getDiagonalNode( const JFacePtr &f, const JNodePtr &v) ;
    static JEdgePtr getOppositeEdge( const JFacePtr &f, const JEdgePtr &e) ;

    static JQuadrilateralPtr down_cast( const JFacePtr &f);

    static JFaceSequence newObjects(size_t n);

    // This function is for concave quads. Arrange vertices so that OpenGL can render them
    // correctly. When you save the mesh, probably the nodesivity may change.
    static int quad_tessalate(const JNodeSequence &orgNodes, JNodeSequence &rotatedNodes);
    static int hexagon_2_quads(const JNodeSequence &hexnodes, JFaceSequence &quads, int start_from);

    ~JQuadrilateral() {}

    JFacePtr getClone() const
    {
        JFacePtr newquad = JQuadrilateral::newObject();
        newquad->setNodes( nodes );
        return newquad;
    }

    string getTypeName() const
    {
        return "Quad";
    }

    int  getTypeID() const
    {
        return JFace::QUADRILATERAL;
    }

    int setNodes( const JNodeSequence &seq)
    {
        assert(seq.size() == 4 );
        int err = JFace::setNodes(seq);
        return err;
    }

    int setNodes(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3)
    {
        JNodeSequence seq(4);
        seq[0] = (v0);
        seq[1] = (v1);
        seq[2] = (v2);
        seq[3] = (v3);
        JFace::setNodes(seq);
        return 0;
    }

    int triangulate( JFaceSequence &newfaces, JNodeSequence &newnodes, int n = 2) const;

    bool isPlanar() const;

};

///////////////////////////////////////////////////////////////////////////////

struct JQuadGeometry {
    static void bilinear_weights(double xi, double eta, vector<double> &weight);
    static double getArea(const Point3D &p0, const Point3D &p1,
                          const Point3D &p2, const Point3D &p3);
    static double distortion(const JNodePtr v0, const JNodePtr v1, const JNodePtr v2,
                             const JNodePtr v3, const JNodePtr v4);
    static int isCyclic(const Point3D &p0, const Point3D &p1, const Point3D &p2, const Point3D &p3);

    static int getBilinearCoords( const JFacePtr face, const Point2D &param, Point3D &xyz);
    static int getUVCoords( const JFacePtr face, const Point3D &xy, Point2D &uv);

    static bool isConvex(const Point3D &p0, const Point3D &p1,
                         const Point3D &p2, const Point3D &p3);
    static bool isConvex(const JFacePtr f);
};
///////////////////////////////////////////////////////////////////////////////

class  JPolygon : public JFace {
public:
    static JPolygonPtr newObject();
    static JPolygonPtr newObject( const JNodeSequence &s);
    static JPolygonPtr down_cast( const JFacePtr &f);
    static JPolygonPtr getCanonical( int n, double r = 1.0);
    static JEdgeSequence getCircle( const Point2D &center, double r, int np);
    static JEdgeSequence getSquare( const Point2D &center, double lx, double ly, int np);
    static JFaceSequence newObjects(size_t n);

    ~JPolygon() {}

    string getTypeName() const
    {
        return "Polygon";
    }

    int  getTypeID() const
    {
        return JFace::POLYGON;
    }

    int setNodes( const JNodeSequence &seq)
    {
        int err = JFace::setNodes(seq);
        return err;
    }

    JFacePtr getClone() const
    {
        JFacePtr newpoly = JPolygon::newObject();
        newpoly->setNodes( nodes );
        return newpoly;
    }
};

///////////////////////////////////////////////////////////////////////////////
class JCell : public JSimplex {
    static std::map<string,string> attribInfo;
public:
    static size_t NumObjectsCreated;
    static size_t AnyObjectModified;

    static const int POLYHEDRON  = 30;
    static const int TETRAHEDRON = 34;
    static const int TRIPRISM    = 36;
    static const int HEXAHEDRON  = 38;

    static JCellPtr newObject( const JNodeSequence &seq);

    static void getRelations03( const JCellPtr &e, JCellSequence &seq);
    static void getRelations23( const JCellPtr &e, JCellSequence &seq);
    static int  get_shared_entities( const JCellPtr &c1, const JCellPtr &c2, JNodeSequence  &nodes);
    static int  get_shared_entities( const JCellPtr &c1, const JCellPtr &c2, JEdgeSequence  &edges);
    static int  get_shared_entities( const JCellPtr &c1, const JCellPtr &c2, JFaceSequence  &faces);
    static JNodePtr getCentroid(const JCellPtr f);

    static JCellPtr down_cast( const JMeshEntityPtr &m);
    static int    registerAttribute(const string &s, const string &type);
    static string getAttributeTypeName( const string &s);

    virtual ~JCell() {}

    virtual string getTypeName() const = 0;
    virtual JCellPtr  getClone() const = 0;

    int getDimension() const
    {
        return 3;
    }

    virtual int getTypeID() const = 0;

    bool isSameAs(const JCellPtr rhs) const
    {
        size_t nSize = rhs->nodes.size();
        if (nSize != nodes.size() ) return 0;

        for (size_t i = 0; i < nSize; i++)
            if (!rhs->hasNode(nodes[i])) return 0;

        return 1;
    }

    bool hasSameNodesID(const JCellPtr rhs) const
    {
        size_t nSize = rhs->nodes.size();
        if (nSize != nodes.size() ) return 0;

        for (size_t i = 0; i < nSize; i++)
            if (!rhs->hasNodeID(nodes[i]->getID() )) return 0;

        return 1;
    }

    virtual int getSize(int e) const = 0;

    virtual int build_lower_entities(int dim ) = 0;
    int remove_unattached_lower_entities();

    // Get the global edges of the simplex. these objects are shared by other objects. In
    // order to get these values, edges and faces must be prebuild from the entire mesh...
    // This is very fragile code, because if the edges and the faces are not build correctly.
    // then the result may be incorrect. Unfortunately, at present, there is dilemma, that
    // who should get control of creating entities, and presently, I have given this
    // responsibility to the mesh object and not the cell object. This is because that
    // very often, we need "All" the unique edges and faces. And if this responsibiliy
    // is given to the cell/faces class, then there are chances of duplication. but it may
    // not be the best way to do things, and I will come to this topic in future.
    // Vertices which are the lowest entities, keep all the higher order entities in
    // there relation.

    const JNodePtr &getNodeAt(int  id) const
    {
        return nodes[id];
    }

    virtual JEdgePtr getEdgeAt(int id) const = 0;
    virtual JEdgePtr getEdgeAt(int id, int &ori) const = 0;

    virtual bool  hasFaceAt(int id) const = 0;
    virtual JFacePtr getFaceAt(int id) const = 0;
    virtual JFacePtr getFaceAt(int id, int &ori) const = 0;

    virtual JEdgeSequence  getEdges() const = 0;
    virtual JFaceSequence  getFaces() const = 0;

    virtual void reverse() {};

    int   setNodes( const JNodeSequence &v)
    {
        return JSimplex::setNodes(v);
    }

    virtual int getOrientation(const JEdgePtr &f) const = 0;
    virtual int getOrientation(const JFacePtr &f) const = 0;

    virtual int refine( JNodeSequence &newnodes, JCellSequence &newcells )
    {
        newnodes.clear();
        newcells.clear();
        return 1;
    }


    JCellPtr explode( double alpha = 0.99) const;
};

struct JCellGeometry {
    static bool isInverted( const JCellPtr &c);
};

///////////////////////////////////////////////////////////////////////////////

inline JCellPtr JCell::down_cast( const JMeshEntityPtr &entity)
{
    JCellPtr cell = boost::dynamic_pointer_cast<JCell>(entity);
    return cell;
}

///////////////////////////////////////////////////////////////////////////////

class JTetrahedron: public JCell {

public:
    static const int NumNodes = 4;
    static const int NumEdges = 6;
    static const int NumFaces = 4;
    static const int NumCells = 1;

    static int getEdgeTopology( int id,  int &v0, int &v2);
    static int getFaceTopology( int id,  int &v0, int &v2, int &v3);

    static JTetrahedronPtr  getCanonical();

    static JTetrahedronPtr  newObject();
    static JTetrahedronPtr  newObject( const JNodeSequence &n);
    static JTetrahedronPtr  newObject( const JNodePtr &v0, const JNodePtr &v1,
                                       const JNodePtr &v2, const JNodePtr &v3);
    static JTetrahedronPtr  down_cast(const JCellPtr &c);

    static JCellSequence  newObjects(size_t n);

    static JMeshPtr  getSchonhardt();

    ~JTetrahedron() {}

    JCellPtr getClone() const
    {
        JTetrahedronPtr tet = JTetrahedron::newObject();
        tet->setNodes( this->getNodes() );
        return tet;
    }

    int getTypeID() const
    {
        return JCell::TETRAHEDRON;
    }

    string getTypeName() const
    {
        return "Tetrahedron";
    }

    int getSize(int e )  const
    {
        if( e == 0) return NumNodes;
        if( e == 1) return NumEdges;
        if( e == 2) return NumFaces;
        if( e == 3) return NumCells;
        return 0;
    }

    int setNodes( const JNodeSequence &v)
    {
        assert( v.size() == 4 );
        return JCell::setNodes(v);
    }

    int setNodes( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3)
    {
        JNodeSequence v(4);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        return setNodes( v );
    }

    int getHexCells( vector<JHexahedron> &h) const;

    int generate( int n );

    int getPosOf( const JEdgePtr &e) const;
    int getPosOf( const JNodePtr &v0, const JNodePtr &v1) const;

    int getPosOf( const JFacePtr &f) const;
    int getPosOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2) const;

    JEdgePtr getEdgeAt(int id) const;
    JEdgePtr getEdgeAt(int id, int &ori) const;
    JEdgePtr getEdgeOf(const JNodePtr &v0, const JNodePtr &v1) const;
    JEdgeSequence  getEdges() const;

    bool  hasFaceAt(int id) const;
    JFacePtr getFaceAt(int id) const;
    JFacePtr getFaceAt(int id, int &ori) const;
    JFacePtr getFaceOf(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2) const;
    JFaceSequence getFaces() const;

    int build_lower_entities(int dim );

    void reverse()
    {
        JNodePtr v2 = nodes[2];
        JNodePtr v3 = nodes[3];
        nodes[2] = v3;
        nodes[3] = v2;
    }

    int getOrientation( const JEdgePtr &e) const;
    int getOrientation( const JFacePtr &f) const;
};

///////////////////////////////////////////////////////////////////////////
class TetGeometry
{
public:
    static double getRegularTetrahedronArea( double a);
    static double getRegularTetrahedronVolume( double a);
    static double getRegularTetrahedronHeight( double a);

    static bool isDegenerate(const JCellPtr &c);
    static bool isInside(const JCellPtr &c, const Point3D &p,  bool include_boundary = 0);
    static bool isOutside(const JCellPtr &c, const Point3D &p, bool include_boundary = 0);
    static bool isOnBoundary(const JCellPtr &c, const Point3D &p);

    static int  getOrientation(const JCellPtr &c);
    static int  getOrientedSide(const JCellPtr &c, const Point3D &p);
    static int  getCircumSphere( const JCellPtr &c, Point3D &center, double &radius);
    static int  getInSphere( const JCellPtr &c, Point3D &center, double &radius);

    static int  getBaryCoordinates( const JCellPtr &cell, const Point3D &xyz, Point4D &baryCoords);
    static int  getXYZCoordinates( const JCellPtr &cell, const Point4D &baryCoord, Point3D &xyz);

    static int  getBaryCoordinates( const double *p0, const double *p1, const double *p2,
                                    const double *xyz, double *bcoords);
    static int  getXYZfromBaryCoordinates( const double *p0, const double *p1, const double *p2,
                                           const double *bccoords, double *zyz);
    static Point3D  getRandomPoint( const JCellPtr &p);

    static void     getDihedralAngles( const JCellPtr &p, vector<double> &a);
    static double   getMinDihedralAngle( const JCellPtr &p);
    static double   getMaxDihedralAngle( const JCellPtr &p);
};

class JHexahedron : public JCell {

public:
    static const int LEFT_SIDE    = 0;
    static const int RIGHT_SIDE   = 1;
    static const int BOTTOM_SIDE  = 2;
    static const int TOP_SIDE     = 3;
    static const int BACK_SIDE    = 4;
    static const int FRONT_SIDE   = 5;

    static const int NumNodes = 8;
    static const int NumEdges = 12;
    static const int NumFaces = 6;
    static const int NumCells = 1;

    static int getEdgeTopology( int id,  int &v0, int &v2);
    static int getFaceTopology( int id,  int &v0, int &v2, int &v3, int &v4);
    static int getFaceID( const string &side);
    static int getCommonEdgeID(int face1, int face2);

    static int rotate_common_face( const JMeshPtr &m, const JHexahedronPtr &hex1, const JHexahedronPtr &hex2, int dir);
    static int rotate_common_edge( const JMeshPtr &m, const JEdgePtr &edge, int dir);

    static JHexahedronPtr  newObject();
    static JHexahedronPtr  newObject( const JNodeSequence &n);
    static JHexahedronPtr  newObject( const JFacePtr &quad1, const JFacePtr &quad2);
    static JHexahedronPtr  newObject( const JFacePtr &quad1, const JFacePtr &quad2,
                                      const JFacePtr &quad3, const JFacePtr &quad4,
                                      const JFacePtr &quad5, const JFacePtr &quad6);
    static JHexahedronPtr  newObject( const JNodePtr &v0,  const JNodePtr &v1,
                                      const JNodePtr &v2,  const JNodePtr &v3,
                                      const JNodePtr &v4,  const JNodePtr &v5,
                                      const JNodePtr &v6,  const JNodePtr &v7);
    static JHexahedronPtr  down_cast(const JCellPtr &c);
    static JCellSequence  newObjects(size_t n);

    ~JHexahedron() {}

    JCellPtr getClone() const
    {
        JHexahedronPtr hex = JHexahedron::newObject();
        hex->setNodes( this->getNodes() );
        return hex;
    }

    static JHexahedronPtr  getCanonical( double l = 1.0);

    int getTypeID() const
    {
        return JCell::HEXAHEDRON;
    }

    string getTypeName() const
    {
        return "Hexahedron";
    }

    int getSize(int e )  const
    {
        if( e == 0) return NumNodes;
        if( e == 1) return NumEdges;
        if( e == 2) return NumFaces;
        if( e == 3) return NumCells;
        return 0;
    }

    int setNodes( const JNodeSequence &v)
    {
        assert( v.size() == 8);
        return JCell::setNodes(v);
    }

    int setNodes( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3,
                  const JNodePtr &v4, const JNodePtr &v5, const JNodePtr &v6, const JNodePtr &v7)
    {
        JNodeSequence v(8);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        v[5] = v5;
        v[6] = v6;
        v[7] = v7;
        return this->setNodes(v);
    }

    void reverse()
    {
        std::swap( nodes[1], nodes[3] );
        std::swap( nodes[5], nodes[7] );
    }

    // Given two nodes of the element, it will return the local id of an edge.
    // (With respect to the element ).
    int getPosOf( const JEdgePtr &e ) const;
    int getPosOf( const JNodePtr &v0, const JNodePtr &v1) const;

    // Given four nodes of the element, it will return the local id of the face.
    // (With respect to the element )
    int getPosOf( const JFacePtr &f ) const;
    int getPosOf( const JNodePtr &v0, const JNodePtr &v1,
                  const JNodePtr &v2, const JNodePtr &v3) const;

    int build_lower_entities(int dim );

    JEdgePtr getEdgeOf(const JNodePtr &v0, const JNodePtr &v1) const;
    JFacePtr getFaceOf(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3) const;

    int getTetrahedra( vector<JCellPtr> &h) const;

    JEdgePtr getEdgeAt(int id) const;
    JEdgePtr getEdgeAt(int id, int &ori) const;
    JEdgeSequence  getEdges() const;

    bool  hasFaceAt(int id) const;
    JFacePtr getFaceAt(int id) const;
    JFacePtr getFaceAt(int id, int &ori) const;
    JFaceSequence getFaces() const;

    // Useful functions for creating dual sheets and chords ..
    JFacePtr getOppositeFace(const JFacePtr &f) const;
    JEdgePtr getDiagonalEdge(const JEdgePtr &e) const;
    JNodePtr getDiagonalNode(const JNodePtr &v) const;

    int  getEdgesAt( const JNodePtr &v , JEdgeSequence &s);
    int  getFacesAt( const JNodePtr &v , JFaceSequence &s);
    int  getFacesAt( const JEdgePtr &e , JFaceSequence &s);

    int   getOppositeFaceNodes( const JNodeSequence &this_side, JNodeSequence &oppo_side )const;

    // Every edge has three parallel edge
    int get_topological_parallel_edges( const JEdgePtr &edge, JEdgeSequence &paredges) const;

    // Every edge four cyclic face neighbors
    int getCyclicFaces(const JEdgePtr &edge, JFaceSequence &efaces ) const;
//   int getCyclicFaces(const Edge *edge, JNodeSequence &vnodes ) const;

    // Exactly three faces per vertex are neighbors.
    int getNeighbors( const JNodePtr &v,   JFaceSequence &faces);
    int get_local_neighbors(  const JNodePtr &v,   JFaceSequence &faces);

    // Exactly two faces per edge are  neighbors.
    int getNeighbors( const JEdgePtr &edge, JFaceSequence &faces);
    int get_local_neighbors(  const JEdgePtr &edge, JFaceSequence &faces);

    // Exactly three edges per vertex are neighbors.
    int getNeighbors( const JNodePtr &v,  JEdgeSequence &edges);
    int get_lobal_neighbors ( const JNodePtr &v,  JEdgeSequence &edges);

    int dice( JEdgePtr &edge, int npieces, JNodeSequence &newnodes, JCellSequence &newcells);

    int getOrientation( const JEdgePtr &edge) const;
    int getOrientation( const JFacePtr &face) const;

    int getSignature(const JFacePtr &face) const;
    int getSignature(const JEdgePtr &edge) const;

    void get_parallel_face_diagonal( const JNodePtr &v0, const JNodePtr &v1,
                                     JNodePtr &v2, JNodePtr &v3);


    void tesseract( JNodeSequence &newnodes, JCellSequence &newcells);
};

////////////////////////////////////////////////////////////////////////////////

class JTriangularPrism : public JCell {
public:
    static const int NumNodes = 6;
    static const int NumEdges = 9;
    static const int NumFaces = 5;
    static const int NumCells = 1;

    static JTriangularPrismPtr newObject();
    static JTriangularPrismPtr getCanonical();

    static JTriangularPrismPtr  down_cast( const JCellPtr &c);

    ~JTriangularPrism() {}

    JCellPtr getClone() const
    {
        JTriangularPrismPtr tp = JTriangularPrism::newObject();
        tp->setNodes( this->getNodes() );
        return tp;
    }

    int getTypeID() const
    {
        return JCell::TRIPRISM;
    }

    string getTypeName() const
    {
        return "TriangularPrism";
    }

    int getSize(int e )  const
    {
        if( e == 0) return NumNodes;
        if( e == 1) return NumEdges;
        if( e == 2) return NumFaces;
        if( e == 3) return NumCells;
        return 0;
    }

    JEdgePtr getEdgeAt(int ) const
    {
        JNoImpl();
    }

    JEdgePtr getEdgeAt(int, int &) const
    {
        JNoImpl();
        JEdgePtr nullPtr;
        return nullPtr;
    }

    JEdgeSequence  getEdges() const;

    bool  hasFaceAt(int ) const
    {
        JNoImpl();
    }

    JFacePtr getFaceAt(int ) const
    {
        JNoImpl();
        JFacePtr nullPtr;
        return nullPtr;
    }

    JFacePtr getFaceAt(int, int &) const
    {
        JNoImpl();
        JFacePtr nullPtr;
        return nullPtr;
    }

    int getOrientation( const JEdgePtr &) const
    {
        JNoImpl();
        return 0;
    }

    int getOrientation( const JFacePtr &) const
    {
        JNoImpl();
        return 0;
    }

    int build_lower_entities(int dim );

    JFaceSequence getFaces() const;
};

class JPolyhedron : public JCell {
public:
    static JPolyhedronPtr newObject();
    static JPolyhedronPtr newObject(const JCellPtr &c);

    JCellPtr getClone() const
    {
        JPolyhedronPtr tp = JPolyhedron ::newObject();
        tp->setNodes( this->getNodes() );
        return tp;
    }

    int getSize(int )  const
    {
        return 0;
    }

    int getTypeID() const
    {
        return JCell::POLYHEDRON;
    }

    string getTypeName() const
    {
        return "Polyhedron";
    }

    int identify() const;

    void addFace( JFacePtr f)
    {
        faces.push_back(f);
    }

    JEdgePtr getEdgeAt(int ) const
    {
        JNoImpl();
    }

    JEdgePtr getEdgeAt(int, int &) const
    {
        JNoImpl();
    }

    bool  hasFaceAt(int ) const
    {
        JNoImpl();
        return 0;
    }

    JFacePtr getFaceAt(int ) const
    {
        JNoImpl();
        JFacePtr nullPtr;
        return nullPtr;
    }

    JFacePtr getFaceAt(int, int &) const
    {
        JNoImpl();
        JFacePtr nullPtr;
        return nullPtr;
    }

    int build_lower_entities(int dim );

    JEdgeSequence getEdges() const
    {
        JEdgeSet eset;
        JEdgeSequence faceedges;
        for(  size_t i = 0; i < faces.size(); i++) {
            faceedges = faces[i]->getEdges();
            for( size_t j = 0; j < faceedges.size(); j++)
                eset.insert( faceedges[j] );
        }

        JEdgeSequence eseq;
        JEdgeSet::const_iterator it;
        for( it = eset.begin(); it != eset.end(); ++it)
            eseq.push_back(*it);

        return eseq;
    };

    JFaceSequence getFaces()  const
    {
        return faces;
    }

    int getOrientation( const JEdgePtr &) const
    {
        JNoImpl();
        return 0;
    }
    int getOrientation( const JFacePtr &) const
    {
        JNoImpl();
        return 0;
    }
private:
    JFaceSequence faces;
};

///////////////////////////////////////////////////////////////////////////////

inline
JNodePtr JNode ::newObject()
{
    JNodePtr v(new JNode);
    assert(v);
    v->setID(GlobalID);
    GlobalID++;
    NumObjectsCreated++;
    return v;
}

///////////////////////////////////////////////////////////////////////////////

inline JNodePtr JNode ::getClone() const
{
    JNodePtr v = newObject();
    assert(v);
    v->setXYZCoords(xyz);
    return v;
}

///////////////////////////////////////////////////////////////////////////////

inline JEdgePtr JEdge::newObject()
{
    JEdgePtr e(new JEdge);
    assert(e);
    NumObjectsCreated++;
    return e;
}

///////////////////////////////////////////////////////////////////////////////

inline JEdgePtr JEdge::newObject( const JNodeSequence &vn)
{
    assert( vn.size() == 2);
    JEdgePtr e = newObject();
    assert(e);
    e->setNodes(vn[0],vn[1] );
    return e;
}
///////////////////////////////////////////////////////////////////////////////

inline JEdgePtr JEdge::newObject( const JNodePtr &v0, const JNodePtr &v1)
{
    JEdgePtr e = newObject();
    assert(e);
    e->setNodes(v0,v1);
    return e;
}

///////////////////////////////////////////////////////////////////////////////

inline JTrianglePtr JTriangle::newObject()
{
    JTrianglePtr f(new JTriangle);
    assert(f);
    NumObjectsCreated++;
    return f;
}
///////////////////////////////////////////////////////////////////////////////

inline JTrianglePtr JTriangle::newObject(const JNodeSequence &vn)
{
    assert( vn.size() == 3 );
    JTrianglePtr f(new JTriangle);
    assert(f);
    f->setNodes(vn);
    NumObjectsCreated++;
    return f;
}
///////////////////////////////////////////////////////////////////////////////

inline JTrianglePtr JTriangle::newObject( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2)
{
    JTrianglePtr f(new JTriangle);
    assert(f);
    f->setNodes( v0, v1, v2);
    assert(f);
    NumObjectsCreated++;
    return f;
}

///////////////////////////////////////////////////////////////////////////////

inline
JTrianglePtr JTriangle::down_cast( const JFacePtr &f)
{
    JTrianglePtr tri;
    tri = boost::dynamic_pointer_cast<JTriangle>(f);
    return tri;
}

inline
JQuadrilateralPtr JQuadrilateral::down_cast( const JFacePtr &f)
{
    JQuadrilateralPtr quad;
    quad = boost::dynamic_pointer_cast<JQuadrilateral>(f);
    return quad;
}

///////////////////////////////////////////////////////////////////////////////
inline JQuadrilateralPtr JQuadrilateral::newObject()
{
    JQuadrilateralPtr f(new JQuadrilateral);
    assert(f);
    NumObjectsCreated++;
    return f;
}
///////////////////////////////////////////////////////////////////////////////

inline JQuadrilateralPtr JQuadrilateral::newObject(const JNodeSequence &s)
{
    JQuadrilateralPtr f(new JQuadrilateral);
    assert(f);
    f->setNodes(s);
    assert(f);
    NumObjectsCreated++;
    return f;
}
///////////////////////////////////////////////////////////////////////////////

inline JQuadrilateralPtr JQuadrilateral::newObject( const JNodePtr &v0, const JNodePtr &v1,
        const JNodePtr &v2, const JNodePtr &v3)
{
    JQuadrilateralPtr f(new JQuadrilateral);
    assert(f);
    f->setNodes( v0, v1, v2, v3);
    assert(f);
    NumObjectsCreated++;
    return f;
}
///////////////////////////////////////////////////////////////////////////////
inline JPolygonPtr JPolygon::newObject(const JNodeSequence &s)
{
    JPolygonPtr f(new JPolygon );
    assert(f);
    f->setNodes(s);
    assert(f);
    NumObjectsCreated++;
    return f;
}

inline JTetrahedronPtr JTetrahedron::newObject( const JNodePtr &v0, const JNodePtr &v1,
        const JNodePtr &v2, const JNodePtr &v3)
{
    JTetrahedronPtr c(new JTetrahedron);
    JNodeSequence qn(4);
    qn[0] = v0;
    qn[1] = v1;
    qn[2] = v2;
    qn[3] = v3;
    c->setNodes( qn);
    NumObjectsCreated++;

    return c;
}
///////////////////////////////////////////////////////////////////////////////

inline
JPolygonPtr JPolygon ::newObject()
{
    JPolygonPtr f(new JPolygon);
    assert(f);
    NumObjectsCreated++;
    return f;
}
///////////////////////////////////////////////////////////////////////////////

/*
inline
JPolygonPtr Polygon::newObject(const JNodeSequence &s)
{
    JPolygonPtr p(new Polygon);
    assert(p);
    p->setNodes(s);
    NumObjectsCreated++;
    return p;
}
*/

///////////////////////////////////////////////////////////////////////////////
inline JPolygonPtr JPolygon::down_cast( const JFacePtr &)
{
    JPolygonPtr poly;
//  poly = dynamic_pointer_cast<Polygon>(f);
    return poly;
}
///////////////////////////////////////////////////////////////////////////////

inline
JTetrahedronPtr JTetrahedron ::newObject()
{
    JTetrahedronPtr c(new JTetrahedron);
    assert(c);
    NumObjectsCreated++;
    return c;
}

///////////////////////////////////////////////////////////////////////////////
inline
JTetrahedronPtr JTetrahedron ::newObject( const JNodeSequence &vn)
{
    JTetrahedronPtr c(new JTetrahedron);
    assert(c);
    assert(vn.size() == 4 );
    c->setNodes(vn);
    NumObjectsCreated++;
    return c;
}

///////////////////////////////////////////////////////////////////////////////

inline
JTriangularPrismPtr JTriangularPrism ::newObject()
{
    JTriangularPrismPtr c(new JTriangularPrism);
    assert(c);
    NumObjectsCreated++;
    return c;
}
///////////////////////////////////////////////////////////////////////////////

inline
JTetrahedronPtr JTetrahedron::down_cast( const JCellPtr &c)
{
    JTetrahedronPtr tet = boost::dynamic_pointer_cast<JTetrahedron>(c);
    return tet;
}

inline
JHexahedronPtr JHexahedron::down_cast( const JCellPtr &c)
{
    JHexahedronPtr hex = boost::dynamic_pointer_cast<JHexahedron>(c);
    return hex;
}

///////////////////////////////////////////////////////////////////////////////

inline
JHexahedronPtr JHexahedron ::newObject()
{
    JHexahedronPtr c(new JHexahedron);
    assert(c);
    NumObjectsCreated++;
    return c;
}

///////////////////////////////////////////////////////////////////////////////
inline
JHexahedronPtr JHexahedron ::newObject( const JNodeSequence &vn )
{
    assert( vn.size() == 8);
    JHexahedronPtr c(new JHexahedron);
    assert(c);
    c->setNodes( vn );
    NumObjectsCreated++;
    return c;
}
///////////////////////////////////////////////////////////////////////////////

inline
JPolyhedronPtr JPolyhedron ::newObject()
{
    JPolyhedronPtr c(new JPolyhedron);
    assert(c);
    NumObjectsCreated++;
    return c;
}

///////////////////////////////////////////////////////////////////////////////

inline
JPolyhedronPtr JPolyhedron ::newObject( const JCellPtr &knownCell)
{
    JPolyhedronPtr c(new JPolyhedron);

    JFaceSequence cellfaces = knownCell->getFaces();
    for( size_t i = 0; i < cellfaces.size(); i++)
        c->addFace( cellfaces[i] );   // Can be unsorted....
    assert(c);
    NumObjectsCreated++;
    return c;
}

///////////////////////////////////////////////////////////////////////////////

inline
JHexahedronPtr JHexahedron ::newObject( const JFacePtr &quad1, const JFacePtr &quad2)
{
    JHexahedronPtr c(new JHexahedron);
    assert(c);

    JNodeSequence nodes(8);
    nodes[0] = quad1->getNodeAt(0);
    nodes[1] = quad1->getNodeAt(1);
    nodes[2] = quad1->getNodeAt(2);
    nodes[3] = quad1->getNodeAt(3);

    nodes[4] = quad2->getNodeAt(0);
    nodes[5] = quad2->getNodeAt(1);
    nodes[6] = quad2->getNodeAt(2);
    nodes[7] = quad2->getNodeAt(3);
    c->setNodes(nodes);
    NumObjectsCreated++;
    return c;
}
///////////////////////////////////////////////////////////////////////////////

/*
inline
JHexahedronPtr Hexahedron::down_cast( JCellPtr &c)
{
    JHexahedronPtr hex = dynamic_pointer_cast<Hexahedron>(c);
    assert(hex);
    return hex;
}
*/

///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager::addRelation(const JNodePtr &vertex, int type)
{
    if (!hasRelation(vertex)) relations0.push_back(vertex);
    if( type == SORTED_RELATIONS ) boost::sort(relations0 );
}
///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager::addRelation(const JEdgePtr &edge, int type)
{
    if (!hasRelation(edge)) relations1.push_back(edge);
    if( type == SORTED_RELATIONS ) boost::sort( relations1 );
}
///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager ::addRelation(const JFacePtr &face, int type)
{
    if (!hasRelation(face)) relations2.push_back(face);
    if( type == SORTED_RELATIONS ) boost::sort(relations2);
}
///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager ::addRelation(const JCellPtr &cell, int type)
{
    if (!hasRelation(cell)) relations3.push_back(cell);

    if( type == SORTED_RELATIONS ) boost::sort(relations3);
}
///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager::removeRelation(const JNodePtr &vertex)
{
    if (hasRelation(vertex))
        boost::remove_erase(relations0, vertex);
}
///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager::removeRelation(const JEdgePtr &edge)
{
    if (hasRelation(edge))
        boost::remove_erase(relations1, edge);
}
///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager ::removeRelation(const JFacePtr &face)
{
    if (hasRelation(face))
        boost::remove_erase(relations2, face);
}
///////////////////////////////////////////////////////////////////////////////

inline
void JRelationManager ::removeRelation(const JCellPtr &cell)
{
    if (hasRelation(cell))
        boost::remove_erase(relations3, cell);
}

///////////////////////////////////////////////////////////////////////////////
inline
void JRelationManager ::clearRelations(int t)
{
    switch( t ) {
    case 0:
        relations0.clear();
        break;
    case 1:
        relations1.clear();
        break;
    case 2:
        relations2.clear();
        break;
    case 3:
        relations3.clear();
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////

inline
bool JRelationManager ::hasRelation(const JNodePtr &vertex) const
{
    if (relations0.empty()) return 0;

    auto it = boost::find(relations0, vertex);
    if( it == boost::end(relations0)) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

inline bool JRelationManager ::hasRelation(const JEdgePtr &edge) const
{
    if (relations1.empty()) return 0;

    auto it = boost::find(relations1, edge);
    if(it == boost::end(relations1)) return 0;

    return 1;
}
///////////////////////////////////////////////////////////////////////////////
inline JEdgePtr JRelationManager :: getEdgeOf( const JEdgePtr &edge) const
{
    return getEdgeOf(edge->getNodeAt(0), edge->getNodeAt(1) );
}

///////////////////////////////////////////////////////////////////////////////
inline JEdgePtr JRelationManager :: getEdgeOf( const JNodePtr &v0, const JNodePtr &v1) const
{

    for( const JEdgePtr &edge: relations1) {
        if( edge->isActive() ) {
            if( !edge->hasNode(v0) ) continue;
            if( !edge->hasNode(v1) ) continue;
            return edge;
        }
    }
    JEdgePtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////
inline JFacePtr JRelationManager :: getFaceOf( const JFacePtr &dummyface) const
{
    return getFaceOf(dummyface->getNodes());
}

inline JFacePtr JRelationManager :: getFaceOf( const JNodeSequence &vnodes ) const
{
    for( const JFacePtr &face: relations2) {
        if( face->isActive() ) {
            int nn = vnodes.size();
            if( face->getSize(0) != nn) continue;
            int found = 1;
            for( int i = 0; i < nn; i++) {
                if( !face->hasNode(vnodes[i]) ) {
                    found = 0;
                    break;
                }
            }
            if( found ) return face;
        }
    }

    JFacePtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////

inline JFacePtr JRelationManager :: getFaceOf( const JNodePtr &v0, const JNodePtr &v1,
        const JNodePtr &v2) const
{
    for( JFacePtr face: relations2) {
        if( face->isActive() ) {
            if( face->getSize(0) != 3) continue;
            if( !face->hasNode(v0) ) continue;
            if( !face->hasNode(v1) ) continue;
            if( !face->hasNode(v2) ) continue;
            return face;
        }
    }
    JFacePtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////

inline JFacePtr JRelationManager :: getFaceOf( const JNodePtr &v0, const JNodePtr &v1,
        const JNodePtr &v2, const JNodePtr &v3) const
{
    for(const JFacePtr &face : relations2) {
        if( face->isActive() ) {
            if( face->getSize(0) != 4) continue;
            if( !face->hasNode(v0) ) continue;
            if( !face->hasNode(v1) ) continue;
            if( !face->hasNode(v2) ) continue;
            if( !face->hasNode(v3) ) continue;
            return face;
        }
    }

    JFacePtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////

inline bool JRelationManager ::hasRelation(const JFacePtr &face) const
{
    if (relations2.empty()) return 0;

    auto it = boost::find(relations2, face);
    if( it == boost::end(relations2)) return 0;

    return 1;
}

///////////////////////////////////////////////////////////////////////////////

inline JCellPtr JRelationManager :: getCellOf( const JCellPtr &dummycell) const
{
    return getCellOf(dummycell->getNodes());
}

///////////////////////////////////////////////////////////////////////////////
inline JCellPtr JRelationManager :: getCellOf( const JNodeSequence &vnodes ) const
{
    for( const JCellPtr &cell: relations3) {
        if( cell->isActive() ) {
            int nn = vnodes.size();
            if( cell->getSize(0) != nn) continue;
            int found = 1;
            for( int i = 0; i < nn; i++) {
                if( !cell->hasNode(vnodes[i]) ) {
                    found = 0;
                    break;
                }
            }
            if( found ) return cell;
        }
    }

    JCellPtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////

inline JCellPtr JRelationManager ::getCellOf( const JNodePtr &v0, const JNodePtr &v1,
        const JNodePtr &v2, const JNodePtr &v3) const
{
    for( const JCellPtr &cell: relations3) {
        if( cell->getSize(0) != 4 ) continue;
        if( !cell->hasNode(v0) ) continue;
        if( !cell->hasNode(v1) ) continue;
        if( !cell->hasNode(v2) ) continue;
        if( !cell->hasNode(v3) ) continue;
        return cell;
    }

    JCellPtr nullPtr;
    return  nullPtr;
}
///////////////////////////////////////////////////////////////////////////////

inline JCellPtr JRelationManager ::getCellOf( const JNodePtr &v0, const JNodePtr &v1,
        const JNodePtr &v2, const JNodePtr &v3,
        const JNodePtr &v4, const JNodePtr &v5,
        const JNodePtr &v6, const JNodePtr &v7 )const
{
    for( JCellPtr cell: relations3) {
        if( cell->getSize(0) != 8 ) continue;
        if( !cell->hasNode(v0) ) continue;
        if( !cell->hasNode(v1) ) continue;
        if( !cell->hasNode(v2) ) continue;
        if( !cell->hasNode(v3) ) continue;
        if( !cell->hasNode(v4) ) continue;
        if( !cell->hasNode(v5) ) continue;
        if( !cell->hasNode(v6) ) continue;
        if( !cell->hasNode(v7) ) continue;
        return cell;
    }

    JCellPtr nullPtr;
    return nullPtr;
}
///////////////////////////////////////////////////////////////////////////////

inline bool JRelationManager ::hasRelation(const JCellPtr &cell) const
{
    if (relations3.empty()) return 0;

    auto it = boost::find(relations3, cell);
    if( it == boost::end(relations3)) return 0;

    return 1;
}
///////////////////////////////////////////////////////////////////////////////

inline int JRelationManager :: getNumRelations( int e ) const
{
    size_t nCount = 0;

    if( e ==  0) {
        nCount = 0;
        for( const JNodePtr &vtx: relations0)
            if(vtx->isActive() ) nCount++;
        return nCount;
    }

    if( e ==  1) {
        nCount = 0;
        for( const JEdgePtr &edge : relations1)
            if( edge->isActive() ) nCount++;
        return nCount;
    }

    if( e ==  2) {
        nCount = 0;
        for( const JFacePtr &face : relations2)
            if( face->isActive() ) nCount++;
        return nCount;
    }

    if( e ==  3) {
        nCount = 0;
        for( const JCellPtr &cell : relations3)
            if( cell->isActive() ) nCount++;
        return nCount;
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

inline
int JRelationManager :: getRelations( JNodeSequence &seq, bool ) const
{
    seq.clear();
    size_t nSize = relations0.size();
    if( nSize == 0) return 1;
    seq.reserve( nSize );
    for( size_t i = 0; i < nSize; i++) {
        if( relations0[i]->isActive() ) seq.push_back( relations0[i] );
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

inline
int JRelationManager :: getRelations( JEdgeSequence &seq, bool ) const
{
    seq.clear();
    size_t nSize = relations1.size();
    if( nSize == 0) return 1;
    seq.reserve( nSize );
    for( size_t i = 0; i< nSize; i++) {
        if( relations1[i]->isActive() ) seq.push_back(relations1[i]);
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

inline
int JRelationManager :: getRelations( JFaceSequence &seq, bool ) const
{
    seq.clear();
    size_t nSize = relations2.size();
    if( nSize == 0) return 1;
    seq.reserve( nSize );
    for( size_t i = 0; i< nSize; i++) {
        if( relations2[i]->isActive()) seq.push_back( relations2[i] );
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

inline
int JRelationManager :: getRelations( JCellSequence &seq, bool ) const
{
    seq.clear();
    size_t nSize = relations3.size();
    if( nSize == 0) return 1;
    seq.reserve( nSize );
    for( size_t i = 0; i< nSize; i++) {
        if( relations3[i]->isActive() ) seq.push_back(relations3[i]);
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
inline void JRelationManager :: prune( int e )
{
    if( e == 0) {
        JNodeSequence::iterator vend;
        vend = remove_if(relations0.begin(), relations0.end(), EntityRemovedPred());
        relations0.erase(vend, relations0.end());
        return;
    }

    if( e == 1) {
        JEdgeSequence::iterator eend;
        eend = remove_if(relations1.begin(), relations1.end(), EntityRemovedPred());
        relations1.erase(eend, relations1.end());
        return;
    }

    if( e == 2) {
        JFaceSequence::iterator fend;
        fend = remove_if(relations2.begin(), relations2.end(), EntityRemovedPred());
        relations2.erase(fend, relations2.end());
        return;
    }

    if( e == 3) {
        JCellSequence::iterator cend;
        cend = remove_if(relations3.begin(), relations3.end(), EntityRemovedPred());
        relations3.erase(cend, relations3.end());
    }

}

///////////////////////////////////////////////////////////////////////////////

inline
void  JRelationManager :: pruneAll()
{
    prune(0);
    prune(1);
    prune(2);
    prune(3);
}

///////////////////////////////////////////////////////////////////////////////

inline
JEdgePtr JSimplex :: getEdgeOf( const JNodePtr &v0, const JNodePtr &v1, bool build)
{
    JNodePtr vhash = std::min(v0, v1);
    JEdgePtr e =  vhash->getEdgeOf(v0,v1);
    if( e ) return e;

    if( build ) {
        e = JEdge::newObject(v0,v1);
        vhash->attach(e);
        return e;
    }

    JEdgePtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////
inline
JFacePtr JSimplex :: getFaceOf( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2,
                               bool build)
{
    JNodePtr vhash = JMath::min_value(v0, v1, v2);
    JFacePtr f = vhash->getFaceOf(v0,v1,v2);
    if( f ) return f;

    if( build ) {
        f = JTriangle::newObject(v0,v1,v2);
        vhash->attach(f);
        return f;
    }

    JFacePtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////
inline
JFacePtr JSimplex :: getFaceOf( const JNodePtr &v0, const JNodePtr &v1,
                               const JNodePtr &v2, const JNodePtr &v3, bool build)
{
    JNodePtr vhash = JMath::min_value(v0, v1, v2, v3);
    JFacePtr f = vhash->getFaceOf(v0,v1,v2,v3);
    if( f ) return f;

    if( build ) {
        f = JQuadrilateral::newObject(v0,v1,v2,v3);
        vhash->attach(f);
        return f;
    }

    JFacePtr nullPtr;
    return nullPtr;
}
///////////////////////////////////////////////////////////////////////////////

inline
JCellPtr JSimplex :: getCellOf( const JNodePtr &v0, const JNodePtr &v1,
                               const JNodePtr &v2, const JNodePtr &v3, bool build)
{
    JNodePtr vhash = JMath::min_value(v0, v1, v2, v3);
    JCellPtr c =  vhash->getCellOf(v0,v1,v2,v3);
    if( c ) return c;

    if( build ) {
        c = JTetrahedron::newObject(v0,v1,v2,v3);
        vhash->attach(c);
        return c;
    }

    JCellPtr nullPtr;
    return nullPtr;
}

///////////////////////////////////////////////////////////////////////////////
inline
JCellPtr JSimplex::getCellOf( const JNodePtr &v0, const JNodePtr &v1,
                             const JNodePtr &v2, const JNodePtr &v3,
                             const JNodePtr &v4, const JNodePtr &v5,
                             const JNodePtr &v6, const JNodePtr &v7,
                             bool build)
{
    JNodePtr vhash1 = JMath::min_value(v0, v1, v2, v3);
    JNodePtr vhash2 = JMath::min_value(v4, v5, v6, v7);
    JNodePtr vhash  = min( vhash1, vhash2);
    JCellPtr c =  vhash->getCellOf(v0,v1,v2,v3);
    if( c ) return c;

    if( build ) {
        c = JHexahedron::newObject(v0,v1,v2,v3,v4,v5,v6,v7);
        vhash->attach(c);
        return c;
    }

    JCellPtr nullPtr;
    return nullPtr;
}

/*

////////////////////////////////////////////////////////////////////////////////

inline
Edge* Hexahedron :: getEdgeAt(int id) const
{
    assert( id >=0 && id < 12);
    int p0, p1;
    getEdgeTopology(id, p0, p1);
    Vertex *v0 = getNodeAt(p0);
    Vertex *v1 = getNodeAt(p1);
    return Simplex::getEdgeOf( v0, v1, 1);
}
////////////////////////////////////////////////////////////////////////////////

inline
Edge* Hexahedron :: getEdgeAt(int id, int &ori) const
{
    assert( id >=0 && id < 12);
    int p0, p1;
    getEdgeTopology(id, p0, p1);
    Vertex *v0 = getNodeAt(p0);
    Vertex *v1 = getNodeAt(p1);
    Edge *edge = Simplex::getEdgeOf( v0, v1, 1);
    ori = edge->getOrientation( v0, v1);
    return edge;
}
////////////////////////////////////////////////////////////////////////////////

inline
Face* Hexahedron :: getFaceAt(int id) const
{
    assert( id >=0 && id < 6);

    int p0, p1, p2, p3;
    getFaceTopology(id, p0, p1, p2, p3);

    Vertex *v0 = getNodeAt(p0);
    Vertex *v1 = getNodeAt(p1);
    Vertex *v2 = getNodeAt(p2);
    Vertex *v3 = getNodeAt(p3);

    return Simplex::getFaceOf( v0, v1, v2, v3, 1);
}
////////////////////////////////////////////////////////////////////////////////

inline
Face* Hexahedron :: getFaceAt(int id, int &ori) const
{
    assert( id >=0 && id < 6);

    int p0, p1, p2, p3;
    getFaceTopology(id, p0, p1, p2, p3);

    Vertex *v0 = getNodeAt(p0);
    Vertex *v1 = getNodeAt(p1);
    Vertex *v2 = getNodeAt(p2);
    Vertex *v3 = getNodeAt(p3);

    Face *fs = Simplex::getFaceOf( v0, v1, v2, v3, 1);
    ori  = fs->getOrientation(v0,v1,v2,v3);
    return fs;
}

////////////////////////////////////////////////////////////////////////////////
struct LowDegreeCompare : public JEntityCompare {
    bool lessCompare(const JNodePtr &v1, const JNodePtr &v2) const
    {
        if( v1->isActive() && v2->isActive() )  {
            size_t d1 = v1->getNumRelations(0);
            size_t d2 = v2->getNumRelations(0);
            return d1 < d2;
        }
        if( v1->isActive()  ) return 1;
        return 0;
    }
};
////////////////////////////////////////////////////////////////////////////////

struct HighVertexDegreeCompare : public JEntityCompare {
    bool lessCompare(const JNodePtr &v1, const JNodePtr &v2) const
    {
        if( v1->isActive() && v2->isActive() ) {
            size_t d1 = v1->getNumRelations(0);
            size_t d2 = v2->getNumRelations(0);
            return d1 > d2;
        }
        if( v1->isActive()  ) return 1;
        return 0;
    }
};
////////////////////////////////////////////////////////////////////////////////

struct LowLayerCompare : public JEntityCompare
{
    bool lessCompare(const JNodePtr &v1, const JNodePtr &v2) const
    {
        int val1, val2;
        if( v1->isActive() && v2->isActive() ) {
            v1->getAttribute("Layer", val1);
            v2->getAttribute("Layer", val2);
            return val1 < val2;
        }
        if( v1->isActive()  ) return 1;
        return 0;
    }
};
*/


