#pragma once

#include "Mesh.hpp"

class JMeshGeodesicVoronoiRegions
{
public:
    static const int  APPROXIMATE_BOUNDARY = 0;
    static const int  EXACT_BOUNDARY       = 1;

    void setMesh(JMeshPtr &m);
    void setNumRegions( int n);
    void extractRegions();
private:
    JMeshPtr mesh;
    int  boundary_type;
    boost::scoped_ptr<JMeshGeodesic>  meshGeodesic;

    struct DistInfo
    {
        double  minDistance;
        size_t  source;
    }
    JNodeSequence  centers;
};

