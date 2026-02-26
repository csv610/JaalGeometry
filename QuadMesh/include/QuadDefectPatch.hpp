#pragma once

#include "Mesh.hpp"
#include "MeshGeodesics.hpp"

namespace Jaal
{
class QDefectivePatch
{
public:
    QDefectivePatch() {
        validPatch = 0;
        specifiedCorners = 0;
        minSingularNodes = 5;
    }

    // Input ...
    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    void setSeed( const JNodePtr &v) {
        apexNode = v;
    }
    void setSeed( const JFacePtr &f) {
        apexFace = f;
    }

    JNodePtr   getSeedNode() const {
        return apexNode;
    }
    JFacePtr   getSeedFace() const {
        return apexFace;
    }

    void setInitialPath( const JNodeSequence &sq) {
        nodePath = sq;
    }

    void setCorners(int n) {
        specifiedCorners = n;
    }

    void setMaximumFaces(int m );

    bool isValid() const;

    int  build();
    int  build( const JFaceSequence &f);

    // Output ...
    JFaceSequence getFaces() const;
    JEdgeSequence getBoundEdges() const {
        return boundEdges;
    }

    void  setMinSingularNodes( int v ) {
        minSingularNodes = v;
    }

    size_t getNumSingularNodes(int where) const;

    JNodeSequence getCorners()  const {
        return cornerNodes;
    }

    vector<int>  getSegments() const {
        return segSize;
    }

    void  setTemplateMesh( const JMeshPtr &m) {
        newSubmesh = m;
    }

    JMeshPtr  getTemplateMesh() {
        return newSubmesh;
    }

    Point3D  getCenter() const;

    void clear();

private:
    JMeshPtr  mesh, newSubmesh, dualGraph;
    JMeshPtr  meshTemplate;
    JNodePtr  apexNode;               // Seed: Irregular vertex to start from.
    JFacePtr  apexFace;               // Seed: Irregular vertex to start from.
    JMeshGeodesics  djkPath;

    struct NonQuadFilter : public JMeshFilter
    {
        bool pass( const JNodePtr &vtx) const
        {
            JFacePtr face;
            int err = vtx->getAttribute("PrimalFace", face);
            assert( !err);
            if( face->getSize(0) == 4) return 1;
            return 0;
        }
    };

    boost::shared_ptr<JMeshFilter> nonQuadFilter;

    JNodeSequence  nodePath;     // Initial joining two irregular nodes..
    bool  validPatch;
    int   specifiedCorners = 0;
    int   minSingularNodes = 3;
    int   maxSearchedFaces = 100;
    bool  stopAtBoundary   = 1;

    JFaceSet faceSet, newFaceSet;
    JEdgeSet edgeSet;
    std::map<JNodePtr, JFaceSet> relations02;

    // Local data ...
    JNodeSet nodes;                   // All the nodes (inner + boundary)
    JNodeSequence cornerNodes;                 // Corners of the blob
    JNodeSequence innerNodes;        // Inner nodes (not on the boundary ) of the blob
    JNodeSequence boundNodes;        // Boundary nodes
    JEdgeSequence boundEdges;          // boundary of the blob.
    vector<int>  segSize;

    // Get the position on the boundary ...
    int getPosOf( const JNodePtr &v) const {
        size_t nSize = boundNodes.size();
        for (size_t i = 0; i <  nSize; i++)
            if (boundNodes[i] == v) return i;
        return -1;
    }
    // Return nodes within the range (src, dst)
    void getBoundNodes( const JNodePtr &src, const JNodePtr &dst, JNodeSequence &s);

    // re-orient boundary nodes so that it starts from a given vertex.
    void start_boundary_loop_from (const JNodePtr &v);

    // Patch creation functions...
    void  buildDualGraph();
    int   initBlob();
    int   updateBoundary();
    int   finalizeBoundary();
    int   expandBlob(const JNodePtr &v);
    int   expandBlob();
    int   getTopologicalAngle( const JNodePtr &v);
    bool  isSimple();
    void  getTemplate();
    bool  isReplacementProfitable();

    void set_boundary_segments();

    double getArea() const {
        double a = 0.0;
        for( const JFacePtr &f : faceSet)
            a += fabs(JFaceGeometry::getArea(f));
        return a;
    }

    double getPerimeter() const {
        double l = 0.0;
        for( const JEdgePtr &e : boundEdges) {
            const JNodePtr &v0 = e->getNodeAt(0);
            const JNodePtr &v1 = e->getNodeAt(0);
            l += JNodeGeometry::getLength(v0, v1);
        }
        return l;
    }

    double getIsoperimeticQuotient() const;

    size_t getSize(int e) const {
        if( e == 0) return innerNodes.size() + boundNodes.size();
        if( e == 2) return faceSet.size();
        return 0;
    }
    int  buildSingularityPatch();
    int  buildTrianglesPatch();
};
}
typedef boost::shared_ptr<Jaal::QDefectivePatch> QDefectivePatchPtr;

