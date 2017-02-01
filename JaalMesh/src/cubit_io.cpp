#include "MeshImporter.hpp"

#include "JaalMoabConverter.hpp"

using namespace Jaal;

JMeshPtr JMeshCubitImporter ::readFile(const string &fname)
{
    JMeshPtr mesh = JMesh::newObject();

#ifdef HAVE_MOAB
    const char *file_name = fname.c_str();

    int ierr;
    iMesh_Instance imesh = 0;
    iMesh_newMesh( NULL, &imesh, &ierr, 0 );
    if (iBase_SUCCESS != ierr) {
        return 0;
    }

    iBase_EntitySetHandle root_set;
    iMesh_getRootSet( imesh, &root_set, &ierr );
    if (iBase_SUCCESS != ierr) {
        iMesh_dtor( imesh, &ierr );
        return 0;
    }

    iMesh_load( imesh, root_set, file_name, 0, &ierr, strlen(file_name), 0 );
    if (iBase_SUCCESS != ierr) {
        iMesh_dtor( imesh, &ierr );
        return 0;
    }

    JaalMoabConverter jm;
    mesh = jm.fromMOAB(imesh, mesh, root_set);
    return 0;
#endif

    return mesh;
}

//##############################################################################
