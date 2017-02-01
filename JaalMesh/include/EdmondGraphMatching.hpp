#pragma once

#include "Mesh.hpp"

class JEdmondGraphMatching
{
   public:
         void setMesh( const JMeshPtr &m) { mesh = m; }
         JMeshPtr getEdgeMatching();
   private:
         JMeshPtr mesh;
};
