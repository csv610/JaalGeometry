#include "RayTracer.hpp"

JRayTracer :: JRayTracer()
{
    rtcInit();
}

JRayTracer ::~JRayTracer()
{
    rtcDeleteScene(sceneID);
}

/////////////////////////////////////////////////////////////////////////////////
void JRayTracer :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    sceneID = rtcNewScene(RTC_SCENE_STATIC,RTC_INTERSECT1);

    int numNodes = m->getSize(0);
    int numFaces = m->getSize(2);

    rtcmesh = rtcNewTriangleMesh (sceneID, RTC_GEOMETRY_STATIC, numNodes, numFaces);

    RTCNode* vertices = (RTCNode*) rtcMapBuffer(sceneID, rtcmesh, RTC_VERTEX_BUFFER);
    for(size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        const Point3D &p3d = vtx->getXYZCoords();
        vertices[i].x = p3d[0];
        vertices[i].y = p3d[1];
        vertices[i].z = p3d[2];
    }
    rtcUnmapBuffer(sceneID, rtcmesh, RTC_VERTEX_BUFFER);

    RTCTriangle* triangles = (RTCTriangle*) rtcMapBuffer(sceneID, rtcmesh, RTC_INDEX_BUFFER);

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        triangles[i].v0 = face->getNodeAt(0)->getID();
        triangles[i].v1 = face->getNodeAt(0)->getID();
        triangles[i].v2 = face->getNodeAt(0)->getID();
    }
    rtcUnmapBuffer(sceneID, rtcmesh, RTC_INDEX_BUFFER);
    rtcCommit(sceneID);
}

/////////////////////////////////////////////////////////////////////////////////

JFacePtr JRayTracer :: getFirstHit(const JFacePtr &fsrc)
{
    RTCRay ray;

    Point3D  pcenter;
    fsrc->getAvgXYZ(pcenter);

    ray.org[0] = pcenter[0];
    ray.org[1] = pcenter[1];
    ray.org[2] = pcenter[2];

    Vec3F normal;
    fsrc->getAttribute("Normal", normal);
    ray.dir[0] = -normal[0];
    ray.dir[1] = -normal[1];
    ray.dir[2] = -normal[2];

    rtcIntersect( sceneID, ray);
    if( ray.primID < 0) return nullptr;
    const JFacePtr  &fdst = mesh->getFaceAt( ray.primID );
    return fdst;
}

