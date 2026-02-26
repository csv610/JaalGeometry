#pragma once

#include "Mesh.hpp"
#include "NearestNeighbours.hpp"
#include "MSTQuadMesher.hpp"
#include "MeshRefine.hpp"
#include "Doublet.hpp"
#include "BasicShapes.hpp"

class JAlphaMSTQuadMesh
{
public:
    void setMesh(const JMeshPtr &m);

    void clear();

    bool isEmpty() { return patchFaces.empty(); } 

    // Extract patch from the circle...
    void setCircle( const JCircle &m);
    void setCenter( const Point3D &c);
    JCircle  getCircle() const { return circle; }
    // Create the submesh where the template will be applied.
    void buildPatch();

    // When we build the submesh, we identify faces, nodes and edges..
    JFaceSequence getFaces()   { return patchFaces; }
    JEdgeSequence getBoundaryEdges(){ return patchEdges; }
    JNodeSequence getBoundaryNodes(){ return patchNodes; }

    int  getNumSingularities() const { return numSingularities; }
    int  getNumBoundaryNodes() const { return patchEdges.size(); }
    
    void setAdaptationFactor( double a) { adaptFactor = a; }
    void setExpectedEdgeLength( double a) { expectedEdgeLength = a; }

    // Only the patch which is identlfied is remeshed.
    void remeshPatch();

    // Entire mesh can be remesh. i.e. No need to identlfy a patch...
    // In this case, we keep the number faces same.
    void remeshAll();

    JFaceSequence getNewFaces() const { return newFaces; }

private:
    JMeshPtr mesh;

    JCircle  circle;

    JNodeSequence patchNodes, newNodes;
    JEdgeSequence patchEdges, newEdges;
    JFaceSequence patchFaces, newFaces;

    double meanArea;
    int    numSingularities = 0;
    double adaptFactor = 1.0;
    double expectedEdgeLength = 0.0;

    boost::scoped_ptr<JNearestNeighbours> nearSearch;
};
