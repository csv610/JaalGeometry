#include "LaplacianMeshDeformation.hpp"

int LaplacianMeshDeformation ::checkPreconditions()
{
    if( constraints.empty() )  {
        cout << "Warning: At least one constraint must be specified " << endl;
        return 1;
    }

    return 0;
}

int LaplacianMeshDeformation :: execute()
{
    if( checkPreconditions() ) return 1;

    if( lap == NULL ) lap = MeshLaplacian::getProduct( algorithm );
    if( lap == NULL ) return 1;

    return 1;
}



