#pragma once

#ifndef POLYSQUARES_H
#define POLYSQUARES_H

#include "Mesh.hpp"
#include "MeshPartitioner.hpp"
#include "LocallyInjectiveMap.hpp"
#include "AllTetMeshGenerator.hpp"
#include "MeshLaplacian.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "PointLocation.hpp"
#include "DelaunayMesh.hpp"
#include "PointLocation.hpp"

using namespace Jaal;

class JPolySquares
{
public:
    static const int MODEL_TO_POLYSQUARE  =  1;
    static const int POLYSQUARE_TO_MODEL  = -1;

    JPolySquares();

    // Input :: A model for which we want to generate a hex mesh..
    // Presently, we assume that the input model consists of
    // triangle mesh. Later, we plan to support quad mesh also.
    void setMesh( const JMeshPtr &m);

    void clear();

    int  setDeformDirection( int d ) {
        deformDirection = d;
    }

    // A Second hop from the tight polycube to integer polycube..
    const JMeshPtr & getPolySquares();

    // Create the hex mesh for the input model...
    const JMeshPtr& getPolyQuadMesh();
    const JMeshPtr& getQuadMesh();
    const JMeshPtr& getSingularGraph();

private:
    vector<Vec3F> edgeNormal;

    JMeshPtr inmesh, polymesh, polyQuadmesh, trimesh, voxmesh, modelQuadmesh;
    int  deformDirection;

    JMeshPtr singularGraph;
    JPointLocation  pointLoc;

    boost::scoped_ptr<JLocallyInjectiveMap>  limDeformer;

    map<double,int> xmap, ymap;

    bool  valid_topology;

    // Step 9: Using locally injective mapping, assign new coordinates to the surface
    //  So that they are perfectly assined to the cube side and in the end optimize
    // the tetmesh...
    int  deformSource( int dir = 1);

    // Step :: Using Barycentric coordinates, project the interior locations to the
    // Original domain ...
    int  projectBack();

    // Generate cubical ( voxelized ) mesh from the integer map..
    JMeshPtr generateQuadMesh();

    void boundarySegmentation();

    int  reorient( const JEdgePtr &f);

    void integerMap( const vector<double> &v, map<double,int> &imap) const;

    void smoothSingularGraph(JEdgeSequence &edges);
    void cornerSingularGraph(JEdgeSequence &edges);
    void getSingularGraph(JEdgeSequence &edges);

    void integerSnap();
    int  setOriginalPosition();
    int  setOriginalPosition( const JNodePtr &v);

    bool isIntegerMapped( const JMeshPtr &m) const;

    // Step 4: The above step may generate many small components which may be due
    //         noise in the surface. Mean Shift may reduce the components and
    // give better segmentation...
    void searchComponents();         // Step 6
    void removeSmallComponents();    // Step 7

    int  incrementDeform();
    void getNodesInBetween( const JNodeSequence &nin, const JNodePtr &vstart, JNodePtr &vend,
                            JNodeSequence &nout);
    void setLinePoints( const JNodeSequence &nodes);
    void removeOutsideVoxels();

};

#endif
