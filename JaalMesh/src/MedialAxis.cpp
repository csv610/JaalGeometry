#include <iostream>
#include <string>

#include "Mesh.hpp"

using namespace Jaal;


int main(int argc, char **argv)
{
    if( argc != 3) {
        cout << "Usage: Executable infile outfile" << endl;
        return 1;
    }

    Mesh *trimesh = new Mesh;
    trimesh->readFromFile( argv[1] );

    delete trimesh;
}

///////////////////////////////////////////////////////////////////////////////
