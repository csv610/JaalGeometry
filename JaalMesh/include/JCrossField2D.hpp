#pragma once

#include "Mesh.hpp"


class JCrossField2D
{
   public:
	void setMesh(const JMeshPtr &m);
        vector<double>     getAngles();
        vector<Vec3D>      getVecField(int n, double len);
   private:
        JMeshPtr mesh;
};
      
