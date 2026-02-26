#pragma once

#include "Mesh.hpp"
#include "PointLocation.hpp"
#include "AlphaMSTQuadMesh.hpp"

class JUVQuadMesher
{
    //  Map Cycle
    //  Patch (XYZ) -> Patch(UV) -> ImprovedPatch(UV)->ImprovedPatch(XYZ)
    //
public:
    void setSubmesh( const JMeshPtr &s);
    void remesh();
private:
    JMeshPtr mesh;   // A patch lying on 3D surface..
    JMeshPtr uvMesh; // A parameteric UV mesh corresponding to the patch.
    JMeshPtr uvTriMesh; // A parameteric UV mesh corresponding to the patch.
    JMeshPtr newUVMesh; // A new improved UV mesh.
    JMeshPtr newMesh; // A new improved UV mesh.
};

