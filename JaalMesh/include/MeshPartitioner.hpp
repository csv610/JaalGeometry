#pragma once 

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

#include<vector>

#include "Mesh.hpp"
#include "MeshTopology.hpp"
#include "basic_math.hpp"

using namespace std;
using namespace Jaal;

class JMeshPartitioner {
    static JLogger *logger;
public:
    JMeshPartitioner() {}

    JMeshPartitioner( const JMeshPtr &m) {
        mesh = m;
    }

    void setMesh(const JMeshPtr &m)
    {
        mesh = m;
    }
    void clear();

    virtual int getPartitions();

    int setElementsPartitionFromNodes();

    int searchComponents();

    // If you are given interfaces, find the regions using flood-filling ...
    int searchRegions();

    // If you are given partitions of faces/cells, find the interfaces ...
    int searchInterfaces();

    // If you are given partitions of faces/cells, find the interfaces ...
    int searchCorners();

    int getNumPartitions() const;
    int getNumInterfaces() const;
    int getNumCorners() const;

    JMeshPtr getSubMesh(int id);

    void   getInterface(int id, JEdgeSequence &edges);
    void   getInterface(int id, JFaceSequence &faces);
    void   getPartition(int id, JFaceSequence &faces);
    void   getCorners( JNodeSequence &nodes);

    void   removeInterface(int id);

    JMeshPtr  getGraph();
    JMeshPtr  getToplogyGraph();

    int getRegionBoundary(int id, JEdgeSequence &edges);

    int enumerate();   // Renumber all partitions starting from zero...
    void savePartitions(int index);
protected:
    JMeshPtr  mesh;

    int   numParts =2;
    vector<int> elemPart, nodePart;
    std::map<int, JEdgeSequence>  edgeGroup;

    // When the faces/cells have been partitioned, create interfaces.
    int setInterface();
    int setCorners();

    int searchEdgesInterface();
    int searchFacesInterface();

    int searchFaceComponents();

    // For 2D mesh, interface consists of  edges..
    void getInterface( JEdgeSequence &edges);

    // For 3D mesh, interface consists of  faces..
    void  getInterface( JFaceSequence &edges);

    int searchRegion( const JFacePtr &f);
};

///////////////////////////////////////////////////////////////////////////////

extern "C" {
#include <metis.h>
}
class JMetisPartitioner : public JMeshPartitioner 
{
public:
    int getPartitions(int nparts);

    vector<JFaceSequence> getPartitions( const JFaceSequence &fseq, int npart);
    vector<JCellSequence> getPartitions( const JCellSequence &fseq, int npart);
    int getTopologicalDisks();
private:
    void init();
/*
    int tri_or_quad_mesh(int etype);
    int tet_or_hex_mesh(int etype);
*/
    int getFacePartitions(int n);
    int getCellPartitions(int n);
    int general_graph();
};

///////////////////////////////////////////////////////////////////////////////

class JMeanShiftSegmentation : public JMeshPartitioner
{
public:
    void segment(const JMeshPtr &m,  const string &s, int nsamples = 3);

private:
    int numSamples;
    string attribname;
    void getNeighbours( const JFacePtr &f, JFaceSequence &seq);
    int  getMeanValue( const JFacePtr &f);
};

