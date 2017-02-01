#pragma once

#ifndef POLYCUBES_H
#define POLYCUBES_H

#include "Mesh.hpp"
#include "MeshPartitioner.hpp"
#include "LocallyInjectiveMap.hpp"
#include "AllTetMeshGenerator.hpp"
#include "MeshLaplacian.hpp"
#include "AllHexMeshGenerator.hpp"
#include "PointLocation.hpp"

using namespace Jaal;

class JPolyCubes
{
public:
    static const int MODEL_TO_POLYCUBE  =  1;
    static const int POLYCUBE_TO_MODEL  = -1;

    static int alignAlong( JEdgeSequence &seq, int axis);

    JPolyCubes();

    // Input :: A model for which we want to generate a hex mesh..
    // Presently, we assume that the input model consists of
    // triangle mesh. Later, we plan to support quad mesh also.
    void setModelMesh( const JMeshPtr &m);

    // Input:  If the base polycube is coming from other sources, then
    // we do not have to do much hardwork of creating basic polycubes.
    // which is often the most time consuming and the weakest part of the
    // algorithm.
    //
    // We are assuming that the basic polycube is topological equaivalent
    // and there is one-to-one matching with the model. It has same number
    // of vertices, same number of triangles.
    //
    // If the user has provided the tetrahedral mesh, we will just discard
    // it and use our own tet mesh, which is just an intermediate mesh
    // to create the hex mesh.
    int  setPolycubes(  const JMeshPtr &m);

    void setDeformDirection( int d ) {
        deformDirection = d;
    }

    // If not, we have our simple way to create basic polycubes, which at
    // present, is not automatic ....
    void  initialSegmentation();
    void  regionGrow();
    int   buildCubicalTopology();
    void  genTetMesh();
    int   alignAlongXYZPlanes();

    // Step :: All the corners of the cube are assigned to a unique integer location.
    JMeshPtr integerSnap();

    // A tetreahedral mesh using the Tight polymesh. We use injective mapping
    // deforme the model mesh into the polybubes mesh and make every to ensure
    // that the mapping is bijective and there are no flip-overs...
    JMeshPtr getPolycubeMesh1();

    // A Second hop from the tight polycube to integer polycube..
    JMeshPtr getPolycubeMesh2();

    // Generate cubical ( voxelized ) mesh from the integer map..
    JMeshPtr generateCubicalMesh();

    // Create the hex mesh for the input model...
    JMeshPtr getPolyHexMesh();
    JMeshPtr getModelHexMesh();

    bool isIntegerMapped( const JMeshPtr &m) const;
    int  getCubeNodes() const {
        return singularNodes.size();
    }
    int  getCubeEdges() const {
        return numInterfaces;
    }
    int  getCubeFaces() const {
        return numPatches;
    }

    // Step 9: Using locally injective mapping, assign new coordinates to the surface
    //  So that they are perfectly assined to the cube side and in the end optimize
    // the tetmesh...
    int  deformSource( int dir = 1);

    // Step :: Using Barycentric coordinates, project the interior locations to the
    // Original domain ...
    int  projectBack();

    void saveForMatlab();

private:
    vector<Vec3F> cubeNormal;

    JMeshPtr insurfmesh, polysurfmesh, tetmesh, voxmesh, hexmesh;
    JMeshPtr sourcemesh, targetmesh;
    int      deformDirection;

    boost::scoped_ptr<JLocallyInjectiveMap>  limDeformer;

    bool  valid_topology;

    int flatten( const JFaceSequence &);

    vector<JEdgeSequence> interfaces;
    vector<JFaceSequence> patches;

    JNodeSequence singularNodes;
    int  numInterfaces, numPatches;

    void seedClusters( const JMeshPtr &m);
    int  buildPatches( const JMeshPtr &m);

    int  reorient( const JFacePtr &f);
    void searchComponent( const JFacePtr &f, int val);

    // Predict the position of nodes on the interface ...
    int targetPositions( const JEdgeSequence &e);

    // Predict the position of the nodes on the patch ...
    int targetPositions( const JFaceSequence &e);

    int setTarget1( const JFaceSequence &s);
    int setTarget2( const JFaceSequence &s);

    void integerMap( const vector<double> &v, map<double,int> &imap) const;

    void scaleModel();

    // Step 1: All the surface normal should be consistently oriented. If not
    //         make them consistent.
    void getConsistentNormals();

    // Step 2: Based on the normals generated in the step 1, Identify faces
    //         which are strongly facing toward the six sides of a cube. If
    // the initial surface mesh is noisy, the segmentation will be noisy and
    // next four steps will cleanup the segmentation ....
    void generateSeeds();           // Step 2

    // Step 3: Faces which have not been identified in the step 2 are marked
    //         "Unknown". Now we use greedy BFS search to fill the faces. which
    // cluster touches the face first, it include the face inside....
    void floodFillRegions();        // Step 3

    // Step 4: The above step may generate many small components which may be due
    //         noise in the surface. Mean Shift may reduce the components and
    // give better segmentation...
    void searchComponents();         // Step 6
    void removeSmallComponents();    // Step 7
    void meanShiftCleanup();         // Step 4;

    // Step 5: After a good segmentation, we create the  network of interfaces and
    // patches....
    void generateCurveNetwork();    // Step 5;

    // Step 7:  Now we assign cube topology to the interfaces and patches. i.e.
    // each patch is between 1-6 and each interface between 1-12..
    void assignFaces();

    int  incrementDeform();
    void getIntegerPoints( const JCellPtr &tet, vector<Point3D> &iPoints);

};

#endif
