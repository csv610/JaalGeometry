#pragma once

#include "Mesh.hpp"
#include "basic_math.hpp"
#include "AllTriMeshGenerator.hpp"
#include <geodesic.h>

////////////////////////////////////////////////////////////////////////////////

class  JMeshGeodesics
{
public:

    static const int  APPROXIMATE_DIJKSTRA = 0;
    static const int  EXACT_DIJKSTRA       = 1;
    static const int  SUBDIVISION_DIJKSTRA = 2;
    static const int  KEENAN_HEAT_FLOW     = 3;

    void setMesh( const JMeshPtr &m);

    void setAlgorithm( int a)
    {
        algorithm = a;
        algo_initialized = 0;
    }

    void setFilter(const JMeshFilterPtr &f)  {
        meshFilterPtr = f;
    }

    JNodeSequence getPath( const JNodePtr &src, const JNodePtr &dst);

    JNodeSequence  getNearestSource( const JNodePtr &v);

    double getDistance( const JNodePtr &src, const JNodePtr &dst);
    int setDistanceField(const JNodeSequence &source);

    JNodeSequence getFarthestNode() const;

protected:
    JMeshPtr mesh;
    bool mesh_initialized  = 0;
    bool algo_initialized  = 0;
    int  algorithm    = 0;
    JMeshFilterPtr  meshFilterPtr;
};
////////////////////////////////////////////////////////////////////////////////

class  JTriMeshGeodesics : public JMeshGeodesics
{
public:
    double getDistance( const JNodePtr &src, const JNodePtr &dst);

    JNodeSequence  getPath( const JNodePtr &src, const JNodePtr &dst);

    int setDistanceField( const JNodeSequence &b);

private:
    std::vector<double>   points;
    std::vector<unsigned> faces;

    boost::scoped_ptr<geodesic::Mesh> googleTriMesh;
    boost::scoped_ptr<geodesic::GeodesicAlgorithmBase> geodesicAlgorithm;;
    int  subdivision_level = 3;
    void initGoogleTriMesh(const JMeshPtr &m);
    void google_geodesics();
    void initMesh();
    void initAlgo();
};

////////////////////////////////////////////////////////////////////////////////

class  JGraphGeodesics : public JMeshGeodesics
{
public:

    JEdgeSequence getPath(const JNodePtr &src, const JNodePtr &dst);
    JEdgeSequence getPath(const JNodePtr &src) {
        return getPath(src, nullptr);
    }

    void   initialize();
    int    setDistanceField(const JNodeSequence &src);
    double getDistance(const JNodePtr &src, const JNodePtr &dst);
    void   clear();

private:

    struct DistInfo {
        DistInfo() {
            distance  = 0.0;
            prevNode  = nullptr;
            thisNode  = nullptr;
        }

        size_t getID() const {
            return thisNode->getID();
        }

        bool operator > ( const DistInfo &rhs) const {
            return distance > rhs.distance;
        }

        bool operator < ( const DistInfo &rhs) const {
            return distance < rhs.distance;
        }

        double   distance;   // Shortest distance from the source to this point
        JNodePtr thisNode;   // Current Node
        JNodePtr prevNode;  // Previous vertex
    };

    typedef std::priority_queue<DistInfo, vector<DistInfo>, greater<DistInfo> > PriorityQ;

    void init();
    void initMesh();
    int  atomicOp( DistInfo &node);
    int  traceback(const JNodePtr &src, const JNodePtr &dst, JNodeSequence &nodeSeq);
    int  getApproxPath(const JNodePtr &src, const JNodePtr &dst, JEdgeSequence &edgeSeq);
    void fastmarching( PriorityQ  &nodeQ, const JNodePtr &dst);
};
