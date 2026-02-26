#pragma once

#include "Mesh.hpp"
#include "tetgen.h"

class AllTetMeshGenerator {
public:
    AllTetMeshGenerator()
    {
//       maxVolume = 0.99*std::numeric_limits<double>::max();
    }

    void     setOptions( const string &s) {
        options = s;
    }

    bool     isDelaunay( const JMeshPtr &p);
    JMeshPtr getConvexHull( const JNodeSequence &m);
    JMeshPtr getConvexHull( const JMeshPtr &m);
    JMeshPtr getReMesh(const JMeshPtr &m, const vector<Point3D> &addPoints);
    JMeshPtr fromHexMesh(const JMeshPtr &m);
    JMeshPtr getSierpinski( int nlevel ) ;
    JMeshPtr genMesh(const JMeshPtr &mesh, const string &cmd);
    JMeshPtr getIsotropicMesh( const JMeshPtr &m, bool modifySurface = 0);
    JMeshPtr getQualityMesh( const JMeshPtr &m, bool modifySurface = 0);
    JMeshPtr getQualityMesh( const JMeshPtr &m, const vector<Point3D> &addPoints, bool modifySurface = 0);
    JMeshPtr getConstrainedMesh( const JMeshPtr &m);

private:
    JMeshPtr mesh, newtetmesh;
    tetgenio  inmesh, outmesh;
    string    options;

    void    preprocess();
    void    postprocess();
    void    replaceNodes( const JNodeSequence &oldnodes, JNodeSequence &newnodes);
};

