#pragma once

#ifndef LLOYD_RELAX_H
#define LLOYD_RELAX_H

#include "Mesh.hpp"
#include "MeshTopology.hpp"
#include "MeshDualGraph.hpp"

using namespace Jaal;

class JLloydMeshOptimizer
{
public:
    JLloydMeshOptimizer() {
        numIters = 1;
        tolerance = 1.0E-05;
    }

    void setMesh( const JMeshPtr &m);
    JMeshPtr getDualGraph() const {
        return dualGraph;
    }

    void setNumIterations( int n ) {
        numIters = n;
    }

    void   setPreserveBoundary( bool p) {
        preserve_boundary = p;
    }
    int    getNumIterations() const {
        return numIters;
    }
    void   setTolerance(double t ) {
        tolerance = t;
    }
    double getMaxResidue() const {
        return maxResidue;
    }

    int smoothAll();
    int smooth( const JNodeSequence &v);

private:
    JMeshPtr mesh, dualGraph;

    int  numIters;
    double tolerance, maxResidue;
    JEdgeSequence neighedges;
    bool  preserve_boundary = 1;

    int updateDual();
    int updatePrimal();
    int setDualEdge(const JNodePtr &vertex);
    int getVoroCentroid(const JNodePtr &v, Point3D &xyz);
};
#endif
