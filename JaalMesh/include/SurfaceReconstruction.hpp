#pragma once

#include "Mesh.hpp"

class JSurfaceReconstruction
{
   public:
        void setMesh( const JMeshPtr &m) { mesh = m; }
        virtual int generate() = 0;
        JMeshPtr getMesh() { return newMesh; }
   
   protected:
        JMeshPtr mesh, newMesh;
};

class JPoissonSurfaceReconstruction : public JSurfaceReconstruction
{
  public:
       int generate();

  private:
       int boundaryType = 3;
       int depth = 8;
       double scale = 1.1;
       double samplesPerNode = 1.5;
       double pointWeight = 4.0;

       void write_ply_file();
};
