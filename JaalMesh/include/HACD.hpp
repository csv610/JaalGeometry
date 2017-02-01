#pragma once

#ifdef CSV

#define _CRT_SECURE_NO_WARNINGS

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <hacdHACD.h>
#include <hacdMicroAllocator.h>

#include "Mesh.hpp"

class JApproxConvexDecomposition
{
public:
    JApproxConvexDecomposition();

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    void setNumClusters( int n ) {
        nClusters = n;
    }
    void setMaxConcavity( double c);
    int  getPartitions();

private:
    JMeshPtr mesh;

    int   nClusters;
    bool  invert, addExtraDistPoints, addFacesPoints;
    double maxConcavity, ccConnectDist;
    size_t targetNTrianglesDecimatedMesh;
    std::vector< HACD::Vec3<HACD::Real> > points;
    std::vector< HACD::Vec3<long> > triangles;
    int preprocess();

};

#endif

