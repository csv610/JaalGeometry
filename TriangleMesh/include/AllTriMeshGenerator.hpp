#pragma once

#include "Mesh.hpp"

struct AllTriMeshGenerator {

    static JMeshPtr  getStructuredMesh(int nx, int ny);
    static JMeshPtr  getCylinder( const Point3D &p0, const Point3D &p1, double radius, int nR);
    static JMeshPtr  getFromQuad2Tri(int *dim, double *len = nullptr, double *org= nullptr, bool texCoord = 0);
    static JMeshPtr  getFromQuad4Tri(int *dim, double *len = nullptr, double *org= nullptr, bool texCoord = 0);
    static JMeshPtr  getFromQuadMesh(const JMeshPtr &quadmesh, JNodeSequence &steiner, int type = 2, bool randomize=0); // or type == 4
    static JMeshPtr  getFromQuadMesh(const JMeshPtr &quadmesh, int type = 2, bool randomize = 0);
    static JMeshPtr  getFromPolyMesh(const JMeshPtr &polymesh, JNodeSequence &steiner);
    static JMeshPtr  getFromPolyMesh(const JMeshPtr &polymesh);
    static JMeshPtr  getSierpinski(int nlevel ) ;
    static int       getIsotropicMesh( const JMeshPtr &m) ;

    static void  getCongruentMesh(const JMeshPtr &mesh);
};

