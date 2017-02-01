#include "Mesh.hpp"
#include "MeshDual.hpp"

using namespace Jaal;

int main(int argc, char **argv)
{
    assert( argc == 3);
    MeshImporter mimp;
    Mesh *tetmesh = mimp.readFrom( argv[1] );
    assert( mesh );

    return 0;
}

