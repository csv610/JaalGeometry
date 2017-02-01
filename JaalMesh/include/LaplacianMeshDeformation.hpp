#ifndef LAPMESHDEFORM_H
#define LAPMESHDEFORM_H

#include "MeshLaplacian.hpp"

class LaplacianMeshDeformation
{
public:
    LaplacianMeshDeformation( JMeshPtr m )
    {   mesh = m;
        algorithm = MeshLaplacian::FLOATER_MEAN_VALUE;
        lap  = NULL;
    }

    ~LaplacianMeshDeformation()
    {
        if( lap ) delete lap;
    }

    void setAlgorithm ( int a ) {
        algorithm = a;
    }
    void setConstraints( const vector<Vertex*>  &c) {
        constraints = c;
    }

    int  execute();
private:
    JMeshPtr mesh;
    int   algorithm;
    MeshLaplacian    *lap;
    vector<Vertex*>  constraints;
    int checkPreconditions();
};

#endif



