#pragma once

#include "Mesh.hpp"
#include <bitset>
#include "MeshLaplacian.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
//   MeshRefine2D:
//   		LongestEdgeRefine2D: ConsistencyRefine2D
//   		ObtuseRefine2D: ConsistencyRefine2D
//   		GradeRefine2D: ConsistencyRefine2D
//   		Refine2D14:  ConsistencyRefine2D
//   		CentroidRefine2D
//   		DelaunayRefine2D
//   		Square3SubDivision
////////////////////////////////////////////////////////////////////////////////
//Class Refine, refines an existing triangulation. This class works for 2D cases
//or the parametric domain of 3D surface triangulation. The user need to provide
//which cells need to be refined. There could be many criteria for deviding the
//cell, so this is not part of the class.
//
//Input :   Connection :  Connectivity of the triangulation ( 1D array)
//      :   ParamCoords:  parametric coordinates of the vertices.
//      :   markFace   :  Faces which are marked for refinement.
//
//Output:   Connection :
//          ParamCoords:
//
//If the user is working with 3D triangulated surface, it is assumed that he
//has some representation ( example NURBS ) from which he can calculate the
//physical Coordinates from the paramCoords. Again this is not the responsiblity
//of the class.
//
//------------------------------------------------------------------------------
//Number of Division     :   New Vertices           :  New Faces
//------------------------------------------------------------------------------
//    2                  :     1                    :  1
//    3                  :     0                    :  2
//    4                  :     3                    :  3
//------------------------------------------------------------------------------
// Default Number of Subdivisions is 2, which is conservative approach, but
// good at equalizing the aspect ratio with respect to neighbouring triangles.
//
// Programmer : Chaman Singh Verma
// Place      : Argonne National Lab.
//              Argonne, IL, USA
//
//
////////////////////////////////////////////////////////////////////////////////
//
/**
 * REFINE AREA         : Increase the density where grid cells have high area/volume
 * REFINE ASPECT_RATIO : Increase the aspect Ratio
 * REFINE CURVATURE    : Create high density mesh near high curvature.
 */

namespace Jaal
{
enum RefinePolicy { CENTROID_PLACEMENT, CIRCUMCENTER_PLACEMENT, LONGEST_EDGE_BISECTION};
enum RefineObjective {REFINE_AREA, REFINE_ASPECT_RATIO, REFINE_CURVATURE};

///////////////////////////////////////////////////////////////////////////////
class JMeshRefiner
{
public:
    JMeshRefiner() {
        uvCoords = 0;
        numIters = 1;
        attribs  = 0;
        selective_refinement = 1;
        make_consistent = 1;
    }

    void setMesh( const JMeshPtr &m ) {
        mesh  = m;
        outmesh = m;
    }

    void updateMesh( const JMeshPtr &m)
    {
        outmesh = m;
    }

    void makeConsistent( bool b )
    {
        make_consistent = b;
    }

    // Calculate the UV-Coords of the new nodes, if they are available at old nodes.
    void setUVCoords( bool v) {
        uvCoords = v;
    }

    // Calculate the attributes of the new nodes, if they are available at old nodes.
    //
    void setAttributes( bool v) {
        attribs = v;
    }

    JNodeSequence getNewNodes() const {
        return newnodes;
    }

    JEdgeSequence getNewEdges() const {
        return newedges;
    }

    JFaceSequence getNewFaces() const {
        return newfaces;
    }

    JCellSequence getNewCells() const {
        return newcells;
    }

    void setNumOfIterations( int i ) {
        numIters = i;
    }

protected:
    JMeshPtr mesh, outmesh;
    int  numIters;
    bool uvCoords;
    bool attribs;
    bool selective_refinement;
    bool make_consistent;

    JNodeSequence newnodes;
    JEdgeSequence newedges;
    JFaceSequence newfaces;
    JCellSequence newcells;

};

////////////////////////////////////////////////////////////////////////////
class JEdgeRefiner : public JMeshRefiner
{
public:
    JEdgeRefiner() {}
    ~JEdgeRefiner() {}

    void refineAll();
    void refine( const JEdgeSequence &e, int n = 2);
};
////////////////////////////////////////////////////////////////////////////

//! \brief 2D Mesh Refinement class.
class JMeshRefiner2D : public JMeshRefiner
{
public:

    JMeshRefiner2D() {
        boundary_split_flag = 0;
    }
    virtual ~JMeshRefiner2D() {}

// void setGeometry(  const iGeom_Instance &g ) { geom = g; }

    void setBoundarySplitFlag( bool f ) {
        boundary_split_flag = f;
    }

    size_t  getNumFacesRefined() const {
        return numfacesRefined;
    }

    // Set Desired Mesh Quality
    void setAspectRatio( double a ) {
        desiredAspectRatio = a;
    }
    void setDesiredArea( double a ) {
        desiredArea     = a;
    }
    void setMinimumAngle( double a ) {
        desiredMinAngle = a;
    }

    void setMaximumAngle( double a ) {
        desiredMaxAngle = a;
    }

    void setFeatureAngle( double a ) {
        featureAngle    = a;
    }
    void setMaximumCells( size_t a ) {
        maxAllowedCells = a;
    }

protected:
    JFaceSet      inConsistentSet;  // Neigbors may have hanging nodes.
    JNodeSequence newNodes;
    JEdgeSequence newEdges;
    JFaceSequence newFaces;

    bool    boundary_split_flag;   // Should we allow boundary edges to be refinment ?
    size_t  numfacesRefined;       // How many faces are refined...

    // Some desiderta .....
    double  desiredAspectRatio, desiredArea;
    double  desiredMinAngle, desiredMaxAngle;
    double  featureAngle;
    size_t  maxAllowedCells;

    JNodePtr setNodeOnEdge( const JEdgePtr &edge ) const;
    JNodePtr getNodeOnEdge( const JEdgePtr &edge ) const;

    bool allow_edge_refinement( const JEdgePtr &edge) const;

    void  append( const JNodePtr &v0 );
    JFacePtr append(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2);     // Triangle
    JFacePtr append(const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3); // Quad

    int initialize();
    int finalize();

    void makeConsistent();
    void upgradeInconsistent();
};

///////////////////////////////////////////////////////////////////////////////

class JTriRefiner: public JMeshRefiner2D
{
public:
    static const int CENTROID_REFINE     = 0;
    static const int BARYCENTRIC_REFINE  = 1;
    static const int ISOTROPIC_REFINE    = 2;
    static const int TRI_2_QUADS_REFINE  = 3;

    JTriRefiner() { }

    int  getQuads(const JFacePtr &f);
    int  refine(const JFacePtr &f, const JEdgePtr &e);
    int  refine(const JFacePtr &f, int type);
    int  refine( JFaceSequence &fs, int type);
    int  refineAll(int type);

    // Refine each triangle into three faces by inserting a new vertex at the centroid...
    int  refine3(const JFacePtr &f);
    int  refine3(JFaceSequence &seq);

    // Refine each triangle into four faces..
    int  refine4(const JFacePtr &f);
    int  refine4(JFaceSequence &seq);

    // Refine each triangle into six  faces..
    int  refine6(const JFacePtr &f);
    int  refine6(JFaceSequence &seq);

    int  refine_vertex_degree_below_5(const JNodePtr &v);
    int  refine_vertex_degree_above_7(const JNodePtr &v);

    int  get567Degree();
    int  refineObtuse(const JFacePtr &f);
    int  refineObtuseTriangles();

    JNodeSequence getNewNodes() const {
        return newNodes;
    }

    JFaceSequence getNewFaces() const {
        return newFaces;
    }

    int  getDeletedNodes(JNodeSequence &n);
    int  getDeletedEdges(JEdgeSequence &n);
    int  getDeletedFaces(JFaceSequence &n);

    int  getConstrainedDelaunay( const JFacePtr &f);
    int  getQualityDelaunay( const JFacePtr &f);

private:
    JNodeSequence nodeneighs;
    JFaceSequence faceneighs;

    int  refineNode3( const JNodePtr &v, const JEdgePtr &edge1, const JEdgePtr &edge2);
    int  refineNode4( const JNodePtr &v, const JEdgePtr &edge1);
    int  refineNode3( const JNodePtr &v);
    int  refineNode4( const JNodePtr &v);
};
////////////////////////////////////////////////////////////////////////////////

class JQuadRefiner: public JMeshRefiner2D
{
public:
    static const int QUAD14 = 14;
    static const int QUAD15 = 15;

    JQuadRefiner() {}

    int  mstRefine( const JFacePtr &f);
    int  refine4(const JFacePtr &f);
    int  refine5(const JFacePtr &f);
    int  refine5(const JFacePtr &f, JNodeSequence &vnew, JFaceSequence &fnew);

    int  refineAll(int scheme);
    int  refineAll(const JFaceSequence &v, int scheme);

    int  refineAt( const JNodePtr &vertex);

    int  refine2( const JFacePtr &f);
    int  refine3( const JFacePtr &f);
    int  refine2( JFaceSequence &s);
    int  refine3( JFaceSequence &s);

    int  refineEdges( const JFacePtr &f);
    int  refineEdges( const JFaceSequence &face, int n);

    int  refine( const JFacePtr &f);
    int  refine( const JFacePtr &f, int *dim);
    int  refineAll( int *dim);
    int  refineAll( const JFacePtr &f, int *dim);

    int  insert_boundary_pillows();
    int  getConstrainedDelaunay( const JFacePtr &f);
    int  getQualityDelaunay( const JFacePtr &f);
    int  getAllQuads(const JFacePtr &f);

private:
    void refine4Consistency( const JFacePtr &f);
    int refine15_convex(const JFacePtr &face);
    int refine15_concave(const JFacePtr &face);
};

struct JMeshRefiner3D : public JMeshRefiner {
};

struct JTetRefiner : JMeshRefiner3D {
};

struct JHexRefiner : JMeshRefiner3D {
    static int  refine(const JCellPtr &c, const JMeshPtr &m);
    static int  refine(const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells);
    static int  refine(const JCellPtr &c, int *nodesdim, const JMeshPtr &m, bool nodes_shared = 0);
    static int  refine(const JCellPtr &c, int *nodesdim, JNodeSequence &newnodes, JCellSequence &newcells, bool nodes_shared = 0);

    static int  refineAt( const JCellPtr &c, const JEdgePtr &e,   JNodeSequence &newnodes, JCellSequence &newcells );

    static int  refine18( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refine18(const JMeshPtr &hexmesh);

    static int  refine17( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refine17(const JMeshPtr &hexmesh);

    static int  refineNode0( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineNode1( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineNode2( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineNode3( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineNode4( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineNode5( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineNode6( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineNode7( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );

    static int  refineEdge0(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge1(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge2(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge3(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge4(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge5(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge6(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge7(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge8(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge9(  const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge10( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
    static int  refineEdge11( const JCellPtr &c, JNodeSequence &newnodes, JCellSequence &newcells );
};
///////////////////////////////////////////////////////////////////////////////

class JCentroidRefiner2D : public JMeshRefiner2D
{
public:
    void refineAll( int n = 1);

    void refine( const JFacePtr &f ) {
        atomicOp( f );
    }

    void operator() (const JFacePtr &f ) {
        refine(f);
    }

private:
    int atomicOp(  const JFacePtr &f);
    int refine_tri( const JFacePtr &f);
    int refine_quad( const JFacePtr &f);
};

///////////////////////////////////////////////////////////////////////////////

class JLongestEdgeRefiner2D : public JMeshRefiner2D
{
public:
    JLongestEdgeRefiner2D() {
        cutOffAspectRatio = 0.50;
    }

    ~JLongestEdgeRefiner2D() {}

    void setCutOffAspectRatio(double asp) {
        cutOffAspectRatio = asp;
    }

    int  refineAll();

private:
    double cutOffAspectRatio;
    int  atomicOp( const JFacePtr &face);
};

///////////////////////////////////////////////////////////////////////////////

class JConsistencyRefiner2D : public JMeshRefiner2D
{
public:
    JConsistencyRefiner2D() { }

    void refineAll();
    void refine( JFaceSequence &f);
private:
    JNodePtr edgeNode[3];

    int   atomicOp( const JFacePtr &f);
    void  quad2tris( const JNodePtr &v0, const JNodePtr &v1, const JNodePtr &v2, const JNodePtr &v3);
    void  makeConsistent1( const JFacePtr &f );
    void  makeConsistent2( const JFacePtr &f );
    void  makeConsistent3( const JFacePtr &f );
    void  makeConsistent();
};

///////////////////////////////////////////////////////////////////////////////

class JRefiner2D14 : public JMeshRefiner2D
{
public:
    JRefiner2D14() {}

    int  refineAll();
    int  refine( JFaceSequence &face);
private:
    int  atomicOp( const JFacePtr &f);
    int  refine_tri( const JFacePtr &f);
    int  refine_quad( const JFacePtr &f);
};

///////////////////////////////////////////////////////////////////////////////

struct JDelaunayRefiner2D : public JMeshRefiner2D
{
    ~JDelaunayRefiner2D() {}

    int  initialize();
    int  finalize();

    int  refineAll() {

        return 0;
    }
};

///////////////////////////////////////////////////////////////////////////////

class JObtuseRefiner2D : public JMeshRefiner2D
{
public:
    JObtuseRefiner2D( ) {
        cutoffAngle = 90.0;
    }

    void setCutOffAngle( double a ) {
        cutoffAngle = std::max(90.0, a);
    }

    int  refineAll();

private:
    int  initialize();
    double cutoffAngle;
    int   atomicOp(const JFacePtr &f);
};

///////////////////////////////////////////////////////////////////////////////

class JGradeRefine2D : public JMeshRefiner2D
{
public:
    int  initialize();
    int  finalize();
    int  execute();

private:
    int atomicOp( const JNodePtr &v);
};

///////////////////////////////////////////////////////////////////////////////

class JAfterRefinement
{
};

}

