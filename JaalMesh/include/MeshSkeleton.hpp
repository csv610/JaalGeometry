#pragma once

#include <algorithm>
#include <fstream>
#include <Eigen/Core>
#include <map>

#include <boost/foreach.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/property_map/property_map.hpp>

#ifdef USE_CGAL
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Surface_mesh.h>
#endif

#include "Mesh.hpp"
#include "AllTriMeshGenerator.hpp"
#include "UniqueCoordinates.hpp"
#include "NearestNeighbours.hpp"
#include "MeshRefine.hpp"

struct JMeshSkeletonBranch
{
    int id;
    int type;   // {00,11,13,33}

    JEdgeSequence  skelEdges;
    JNodeSequence  surfNodes;
    JFaceSequence  surfFaces;
//  vector<JEdgeSequence>  getContours();

    bool isTopologicalCylinder();
    void recenter();
    void reparameterize() {
         JEdgeGeometry::makeUniform(skelEdges);
    }

    JFaceSequence getCap(int i = 0);
    JNodePtr  getCapCenter(int i = 0);

    void refine(int numseg);
    void fitCylinder();

    double getMeanRadius();
    double getLength();
};

typedef boost::shared_ptr<JMeshSkeletonBranch> JMeshSkeletonBranchPtr;

class JMeshSkeleton
{
#ifdef USE_CGAL
    typedef CGAL::Simple_cartesian<double>   Kernel;
    typedef Kernel::Point_3                  Point;
    typedef CGAL::Surface_mesh<Point>  TriMesh;
    typedef boost::graph_traits<TriMesh>::vertex_descriptor    vertex_descriptor;
    typedef CGAL::Mean_curvature_flow_skeletonization<TriMesh> Skeletonization;
    typedef Skeletonization::Skeleton        Skeleton;
    typedef Skeleton::vertex_descriptor      Skeleton_vertex;
    typedef Skeleton::edge_descriptor        Skeleton_edge;

    struct SkelPolylines {
        const Skeleton& skeleton;
        std::ofstream& out;
        int polyline_size;
        std::stringstream sstr;

        SkelPolylines(const Skeleton& skeleton, std::ofstream& out)
            : skeleton(skeleton), out(out)
        {}
        void start_new_polyline() {
            polyline_size=0;
            sstr.str("");
            sstr.clear();
        }
        void add_node(Skeleton_vertex v) {
            ++polyline_size;
            sstr << " " << skeleton[v].point << endl;
        }
        void end_polyline()
        {
            out << polyline_size << endl;
            out << sstr.str() << endl;
        }
    };
    TriMesh  tmesh;
    boost::scoped_ptr<Skeletonization> mcs;
#endif

public:
    void setMesh(const JMeshPtr &m) {
        inMesh = m;
//      mcs.reset();
    }
    // Return the input mesh for which the skeleton was required. Locally
    // we a differentt mesh and update it. The input mesh remains unchanged.
    JMeshPtr getMesh() const { return inMesh; }

    void clear();
       
    void setRefinement(int n) {
        numRefine = n;
    }

    void contractGeometry();

    // We will be working with a duplicate mesh which may be refined or modified..
    JMeshPtr getWorkingMesh();

    void setSpeedQuality( double q) {
        speedQuality = max(0.01, q);
    }
    void setCenterQuality( double q)
    {
        centerQuality = max(0.01,q);
    }

    JMeshPtr getSkeleton() {
        genSkeleton();
        return skelGraph;
    }

    JMeshPtr getTopoGraph();

    JNodeSequence getLeafNodes()  const {
        return leafNodes;
    }

    JNodeSequence getJunctionNodes() const {
        return junctionNodes;
    }

    JNodeSequence getExtremeNodes() const;

    size_t getNumBranches() const {
        return branches.size();
    }

    JMeshSkeletonBranchPtr  getBranch(int id){
         if( id < (int)branches.size() ) return branches[id];
         return nullptr;
    }
    vector<JMeshSkeletonBranchPtr>  getBranches() const { return branches;} 

    void splitBranch(int id);
    void removeBranch(int id);

private:
    JMeshPtr inMesh;
    JMeshPtr triMesh;
    JMeshPtr skelGraph;
    JMeshPtr topoGraph;

    double speedQuality = 0.1;
    double centerQuality = 0.2;
    int    numRefine   = 0;

    vector<JMeshSkeletonBranchPtr> branches;
    vector<Point3D> faceCentroid;
    vector<double>  nodeSDF;


    JNodeSequence leafNodes;
    JNodeSequence junctionNodes;

    void initMesh();
    void storeTriMesh();
    void genSkeleton();
    void segmentBranches();
    void segmentSurface();
    void processGraph();

    void      setNodeRadius( const JNodePtr &v);
    JNodePtr  getStartNode( JEdgeSequence &s);
    JMeshSkeletonBranchPtr  getNewBranch(int id);
    JNodePtr  getNearestSurfaceNode( const Point3D &p);
    void classify( JMeshSkeletonBranchPtr &b);
};

typedef boost::shared_ptr<JMeshSkeleton> JMeshSkeletonPtr;
