#pragma once

#include "Mesh.hpp"
#include "MeshOptimization.hpp"
#include "MeshLaplacian.hpp"
#include "LocallyInjectiveMap.hpp"
#include "LloydOptimizer.hpp"
#include "AllTriMeshGenerator.hpp"
#include "AllTetMeshGenerator.hpp"

class JMeshUntangle
{
public:
    JMeshUntangle();

    ~JMeshUntangle();

    // Input mesh containing tangled elements...
    void setMesh( const JMeshPtr &m);

    // How many elements are inverted in the tri/tet mesh ..
    size_t countInverted( const JMeshPtr &m);

    // Back projection energy ....
    void setEnergyType( int e) {
        energyType = e;
    }

    // settting minimum number of inflation iterations ?
    void setMaxInflationSteps( int n) {
        maxInflationSteps = n;
    }
    void   setMaxBackProjectionSteps( int n) {
        maxBackProjectionSteps = n;
    }

    // Inflate all the boundary nodes...
    bool   setInflateAll( bool b) {
        inflateAll = b;
    }

    // Return the inflated mesh, which sould be tangle free ...
    JMeshPtr getInflatedMesh() const { return inflatedMesh; }

    int  startInflation();

    // Inflated mesh to be bring back to the original shape along with the
    // uninverted elements...
    int  startBackProjection();

    int  optimize();

    int execute() {
         int err = 0;
         err = startInflation();
         if( err ) return err;

         err = startBackProjection();
         if( err ) return err;
         return 0;
    }

    vector<double>  getOffset() const {
        return offset;
    }

    double  getMaxDistance() const;

    JFaceSequence getInvertedFaces(); // get Inverted faces in the inflated mesh.
    double  getAreaChanged() const {
        return 100*(area[1] - area[0] )/area[0];
    }

    JCellSequence getInvertedCells(); // get Inverted cells in the inflated mesh.
    double  getVolumeChanged() const
    {
        return 100*(volume[1] - volume[0])/volume[0];
    }

private:
    JMeshPtr  srcMesh;      // The input mesh; Could be any mesh..
    JMeshPtr  inflatedMesh; // Always a simplicial mesh...
    JMeshPtr  triMesh, tetMesh; // Simplicial mesh ..

    vector<double> offset;   // Offset of the boundary nodes from the inflated model.
    double area[2], volume[2]; // Keep track of area and volume of the mesh. "0" indicates
                                // the value of the original mesh and "1" indicates the 
                                // value of the completeley inflated mesh which doesn't
                                // contain any inverted element.

    int  entityDim;             // Are we solving  2D or a 3D problem.
    int  energyType;
    int  maxInflationSteps  = 100;
    int  maxBackProjectionSteps = 100;
    bool org_simplicial = 0;  // Is the source mesh is simplicial.

    bool inflateAll    = 0;  // Do we want to updated all the boundary nodes or subset of nodes where 
                             // inversion occurs ?

    JLaplaceMeshSmoother  lap;   // Laplacian Smoothing ...
    JLloydMeshOptimizer   lloyd; // Lloyd smoothing ...
    JLocallyInjectiveMap  lim;   // Used for back projection ....
    JMeshNonlinearOptimization mopt; // Once the mesh is inversion-free, apply aggressive optimization..

    JNodeSet concaveCorners;     // Boundary nodes where inversion occurs. Usually there are concave corners.
    int smoothCorner( const JNodePtr &v); // Smooth a given node..
    int smoothBoundary();

    JMeshPtr getSimplicialMesh();
    JMeshPtr getInflated2D();
    JMeshPtr getInflated3D();

    int  updateSource();
    int  clear();  // remove inflated mesh,if it wasn't derived from the srcmesh and other local dataset.
};
