#pragma once

#include "Mesh.hpp"
#include "MeshPartitioner.hpp"
#include "MeshGeodesics.hpp"
#include "MeshMatrix.hpp"

#ifdef USE_IGL
#include <igl/matlab/matlabinterface.h>
#endif

class JMeshSamples
{
public:
    static const int  RANDOM     = 0;
    static const int  REGULAR    = 1;
    static const int  OCTREE     = 2;
    static const int  GEODESIC   = 3;
    static const int  BIHARMONIC = 4;
    static const int  METIS      = 5;
    static const int  POISSON    = 6;
    static const int  LLOYDS     = 7;

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    void setAlgorithm( int a ) {
        algorithm = a;
    }

    JNodeSequence getNodeSamples( int n);

    JNodePtr getRandomNode( bool fromWhere = 0) const;
    JEdgePtr getRandomEdge( bool fromWhere = 0) const;
    JFacePtr getRandomFace( bool fromWhere = 0) const;
    JCellPtr getRandomCell() const;

private:
    JMeshPtr mesh;
    int  algorithm = RANDOM;

    JNodePtr  getRandomNode( const JMeshPtr &m);

    void getRandomSamples(int n, JNodeSequence &nodes);
    void getMetisSamples(int n, JNodeSequence &nodes);
    void getGeodesicSamples(int n, JNodeSequence &nodes);
    void getIGLSamples(int n, JNodeSequence &nodes);
    void getRegularSamples(int n, JNodeSequence &nodes);
};
