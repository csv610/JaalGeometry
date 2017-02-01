#pragma once


#include "Mesh.hpp"

class JMeshRayTracer
{
    public:
	void addMesh( const JMeshPtr &m);
	void addPlane( const JQuadrilateral &m);
        JMeshPtr  getShadows();
    private:
        vector<JMeshPtr>   meshes;
        vector<JQuadrilateral> planes;

};
