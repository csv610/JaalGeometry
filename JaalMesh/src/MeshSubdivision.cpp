#include "MeshSubdivision.hpp"

JMeshPtr  JMeshSubdivision :: getSubdivided()
{
    if( inMesh == nullptr) return nullptr;

    JMeshPtr newmesh;
#ifdef USE_CGAL
    JMeshIO   meshio;
    meshio.saveAs(inMesh, "tmp.off");

    CGALPolyhedron P;

    ifstream ifile("tmp.off", ios::in);
    ifile  >> P; // read the .off

    switch(algo)
    {
    case CATMULL_CLARK:
        CGAL::Subdivision_method_3::CatmullClark_subdivision(P,numIters);
        break;
    case LOOP:
        CGAL::Subdivision_method_3::Loop_subdivision(P,numIters);
        break;
    case SQRT3:
        CGAL::Subdivision_method_3::Sqrt3_subdivision(P,numIters);
        break;
    case DOOSABIN:
        CGAL::Subdivision_method_3::DooSabin_subdivision(P,numIters);
        break;
    }

    ofstream ofile("tmp.off", ios::out);
    ofile << P; // write the .off
    newmesh = meshio.readFile("tmp.off");
#endif

    return newmesh;
}
