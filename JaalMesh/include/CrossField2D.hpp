#pragma once

#include "Mesh.hpp"
#include "HarmonicField.hpp"

class JCrossField2D 
{
   public:
        void setMesh(const JMeshPtr &m);

        void genVecField();

        const vector<double> &getCos4tField() const
              { return cos4t; }
        const vector<double> &getSin4tField() const
              { return sin4t; }

   private:
        JMeshPtr mesh;
        vector<double> cos4t, sin4t, theta;
        void  getBoundaryAngles();
        void  solveSystem();
};
   
