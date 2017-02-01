#pragma once

#ifndef ALLHEX_H
#define ALLHEX_H

#include "Mesh.hpp"
#include "Curve.hpp"
#include "MeshTopology.hpp"
#include "DynamicEulerCharacteristic.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////

class JStructuredMesh3D
{
public:
    JStructuredMesh3D();

    void  setCellDimensions( const vector<int> &dim);
    vector<int>  getCellDimensions() const {
        vector<int> d(3);
        d[0] = cellDim[0];
        d[1] = cellDim[1];
        d[2] = cellDim[2];
        return d;
    }

    JMeshPtr getMesh() const {
        return mesh;
    }

    JCellPtr  getCellAt( int i, int j, int k) const;
    JNodePtr  getNodeAt( int i, int j, int k) const;
    int        getIJK( const JCellPtr &c, int &i, int &j, int &k) const;
    int        getIJK( const JNodePtr &n, int &i, int &j, int &k) const;

    size_t     getSize(int e) const;

    int   getRelations(  const JNodePtr &, JNodeSequence &) const;
    int   getRelations(  const JNodePtr &, JCellSequence &) const;
    int   getRelations6( const JCellPtr &, JCellSequence &) const;
    int   getRelations8( const JCellPtr &, JCellSequence &) const;

    void  setActiveCells(const string &);
    int   getActiveCells( JCellSequence &cell);

    bool  isBoundaryNode( const JNodePtr &) const;
    bool  isBoundaryCell( const JCellPtr &) const;

    int   getBoundary( JNodeSequence  &seq) const;
    int   getBoundary( JEdgeSequence  &seq) const;
    int   getBoundary( JFaceSequence  &seq) const;
    int   getBoundary( JCellSequence  &seq) const;

    int   getNumSingularPoints() const;
    bool  isSingular( const JNodePtr &p) const;

    int   searchComponents();

    int   getEulerCharacteristic() const;
    bool  isActive( int i, int j, int k) const;
    int   fillSpace();

    JMeshPtr mesh;
    boost::tuple<int,int> getRange(int dir) const;
    int  cellDim[3], nodeDim[3];
    size_t numNodes, numCells;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////

class JVoxelMesh : public JStructuredMesh3D
{
public:
    // A background mesh is a complete structured mesh in which we will switch on/off the cells....


    void   setBackgroundMesh( const JMeshPtr &m, int *gridim);
    JMeshPtr getBackgroundMesh(int *gridim, double *length = nullptr, double *origin = nullptr);

    JMeshPtr getBackgroundMesh() const {
        return mesh;
    }

    int  getNumSingularPoints() const;
    int   makeSolid();

    JMeshPtr getModelMesh();
    string   getBitString() const;

    int      getCellsInPlane(int id, int dir, JCellSequence &cells);

    boost::tuple<int,int> getRange(int dir) const;

    double getVolume() const;
    double getSurfaceArea() const;
    int    saveAs( const string &s) const;
    int    readFrom( const string &s);
    int    searchComponents();
private:
    JMeshPtr modelmesh;
    double  origin[3], length[3];
    boost::tuple<int,int>  getXRange() const;
    boost::tuple<int,int>  getYRange() const;
    boost::tuple<int,int>  getZRange() const;
    bool isBoundary(int i, int j, int k) const;
};

typedef boost::shared_ptr<JVoxelMesh>  JVoxelMeshPtr;

///////////////////////////////////////////////////////////////////////////////////////////////////////

struct JMonotoneVoxelizer
{
    void  setVoxelMesh( const JVoxelMeshPtr &vmesh) {
        voxmesh = vmesh;
    }
    void  optimize();

private:
    JVoxelMeshPtr voxmesh;

    JDynamicEulerCharacteristic   deuler;
    vector<float> distance;
    deque<JCellPtr> addGroup, delGroup;

    void makeGroups();
    bool isSingular(const JNodePtr &v) const;
    int  getNumSingularPoints( const JCellPtr &c) const;
    int  getChangeInSingularPointsAfterAddition( const JCellPtr &c);
    int  getChangeInSingularPointsAfterDeletion( const JCellPtr &c);
    void setDistance();

};

//////////////////////////////////////////////////////////////////////////////////////////////////////

class AllHexMeshGenerator {
public:
    static JMeshPtr SchneiderPyramid();
    static JMeshPtr getStructuredMesh(int *gridim, double *length = nullptr, double *origin = nullptr);
    static JMeshPtr getGeodeTemplate();
    static JMeshPtr getSojiTemplate();
    static JMeshPtr getSierpinski( int nlevel ) ;

    // Stack templates....
    static int stackQuadMesh( const JMeshPtr &mesh1, const JMeshPtr &mesh2, JCellSequence &newcells);
    static JMeshPtr stackQuadMesh(vector<JMeshPtr> &quadmesh, bool closed = 0);

    // Tetrahedral to Hexhedral mesh ...
    static JMeshPtr Tetrahedral2Hexahedral(const JMeshPtr &m);

    static bool  isDoublet(const JNodePtr &v);
    static bool  isSinglet(const JNodePtr &v);
    static int   removeDoublet(const JMeshPtr &m, const JNodePtr &v);
    static void  searchDoublets( const JMeshPtr &m, JNodeSequence &v);
    static void  collapse( const JMeshPtr &m, const JNodePtr &v0, const JNodePtr &v1);
    static int   rotate_boundary_edge( const JMeshPtr &mesh, const JEdgePtr &oldedge, JEdgeSequence &affected);
    static int   rotate_common_face(const JMeshPtr &mesh, const JHexahedronPtr &hex1, const JHexahedronPtr &hex2, int dir);
    static int   rotate3hex( const JEdgePtr &edge, JNodeSequence &newnodes, JCellSequence &newcells);
    static int   rotate3hex( const JEdgePtr &edge, const JMeshPtr &m);

    int setSurfaceQuadMesh(const JMeshPtr &m);

    JMeshPtr  getTopologicalHexMesh();
private:
    JMeshPtr surfQuadmesh;
    JMeshPtr tetmesh, hexmesh;
    void genTetVolMesh();
};

#endif

