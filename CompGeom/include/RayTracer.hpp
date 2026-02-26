#pragma once
#include "Mesh.hpp"

#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>

struct RTCNode  {
    float x,y,z,r;
};
struct RTCTriangle {
    int v0, v1, v2;
};

class JRayTracer
{
public:
    JRayTracer();
    ~JRayTracer();

    void setMesh( const JMeshPtr &m);
    JFacePtr getFirstHit( const JFacePtr &f);
private:
    JMeshPtr mesh;
    RTCScene sceneID;
    RTCRay   ray;
    unsigned int rtcmesh;
};
