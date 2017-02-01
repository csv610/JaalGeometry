#pragma once

#include "Mesh.hpp"
#include "MeshDual.hpp"
#include "MeshLaplacian.hpp"
#include "AllTriMeshGenerator.hpp"
#include "BarycentricCoords.hpp"
#include "MeshContour.hpp"
#include "DelaunayMesh.hpp"
#include "PolyPartitioner.hpp"
#include "QuadDefectPatch.hpp"
#include "MeshRefine.hpp"
#include "LloydOptimizer.hpp"
#include "Geometry.hpp"
#include "MeshUntangle.hpp"
#include "MeshAffineTransforms.hpp"

class JMSTQuadMesher
{
    typedef boost::shared_ptr<QDefectivePatch>  QDefectivePatchPtr;
    static const int   RANDOM_REMOVAL   = 0;
    static const int   BOUNDARY_REMOVAL = 1;
    static const int   INTERIOR_REMOVAL = 2;

public:
    JMSTQuadMesher()
    {
        cornerAngle    = 15.0;
        geomDilation   = 2.0;
        select_global_patch = 0;
        handleInvertedElements  = 0;
        geodesic_path  = 0;
        minSingularNodes = 3;
    }

    void setBoundarySplit( int v );
    JMeshPtr getTemplate( const vector<int> &segments);
    JMeshPtr getTemplate( const JMeshPtr &m, int side);
    JMeshPtr getQuadTemplate( const JFacePtr &f, double edgelen);

    void setCornerAngle(double a) {
        cornerAngle = a;
    }
    void setGeometricDilation( double d) {
        geomDilation = d;
    }
    void setOptimize( bool opt) {
        optimize = opt;
    }

    void      setMesh( const JMeshPtr &m);
    JMeshPtr  getBoundarySingularGraph();
    JMeshPtr  getInteriorSingularGraph();
//  JNodeSequence getInteriorSingularNodes();
    JNodeSequence getSingularNodes();

    JMeshPtr  getQuadMesh();

    void setMinSingularNodes( int v ) {
        minSingularNodes = v;
    }
    void selectGlobalPatch(int r) {
        select_global_patch = r;
    }
    void setGeodesicPath(int v) {
        geodesic_path = v;
    }

    void setSingularityRemovalPolicy( int p ) {
        singularity_removal_policy = p;
        reorderSingularities();
    }

    QDefectivePatchPtr getDefectivePatch( const JNodePtr &v);
    QDefectivePatchPtr getDefectivePatch( const JFacePtr &f);
    QDefectivePatchPtr getDefectivePatch( const JNodeSequence &v);
    QDefectivePatchPtr getAnyDefectivePatch();

    int  remesh(const QDefectivePatchPtr &df);

    int  setHandleInvertedElements( bool v) {
        handleInvertedElements = v;
    }
    JNodeSequence getCorners() const { return patchCorners; }

    JMeshPtr getPatch( const JNodeSequence &srcnodes, const vector<int> &segments);
    JMeshPtr getAdaptivePatch( const JNodeSequence &srcnodes, const vector<int> &segments, double factor);

private:

    double cornerAngle, geomDilation, avg_edgelen;
    bool   planarMesh = 0;  // A planar mesh will have boundary nodes...
    bool   select_global_patch;
    bool   handleInvertedElements;
    bool   geodesic_path;
    bool   optimize = 0;
    int    minSingularNodes;
    int    singularity_removal_policy = RANDOM_REMOVAL;

    JMeshPtr mesh, singularMesh, quadmesh;
    std::multimap<double,JNodePtr>  distMap;

    std::deque<JNodePtr> singularNodesQ;
    JNodeSequence  defect3nodes, defect5nodes, higherDefectNodes;
    JNodeSequence  patchCorners;

    JNodeSequence domainBoundaryNodes, qnodes;
    JNodeSequence newnodes;
    JFaceSequence newfaces;

    void clear()
    {
        newfaces.clear();
        newnodes.clear();
    }

    void updateDistanceMap( const JNodePtr &v);
    int splitPolygon( const JFacePtr &face);
    int nsidedQuads( const JFacePtr &face);
    int stitchBoundary( const JMeshPtr &mtemplate, const JNodeSequence &boundary);

    void pillowPatch( const JMeshPtr &patch, const vector<int> &segments, int minQuads = 0);

    JMeshPtr getPatch( const JNodeSequence &srcnodes, const JNodeSequence  &corners);
    JMeshPtr getLoftPatch( const JNodeSequence &v1, const JNodeSequence &bridge1,
                           const JNodeSequence &v2, const JNodeSequence &bridge2);

    int applyLIM(const JMeshPtr &m);

    int  segmentBoundary();
    int  getConvexPatches();
    int  genEdgeNodes();
    int  getBoundarySingularGraph( JEdgeSequence &);
    int  make_even_segments(JEdgeSequence &edgeloop, const JEdgePtr &segment);
    void refineSegment( const JNodePtr &v0, const JNodePtr &v1);
    void straighten( const JEdgePtr &edge);
    void getVoxelMesh();
    void reorderSingularities();
    void UV2XYZMap(const JMeshPtr &mesh, const JNodeSequence &v, const vector<int> &segments);
    QDefectivePatchPtr buildCircularPatch( const JNodePtr &v);
    QDefectivePatchPtr getGlobalPatch();
};
///////////////////////////////////////////////////////////////////////////////
