#include "MeshRayTracer.hpp"


void JMeshRayTracer :: initialize()
{
    rtcInit ();
    RTCScene scene = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);

 

    rtcCommit (scene);
    rtcExit ();
}

void JMeshRayTracer :: addGeometry( RTCScene scene)
{
      size_t numFaces = mesh->getSize(0);
      size_t numNodes = mesh->getSize(1);


}
