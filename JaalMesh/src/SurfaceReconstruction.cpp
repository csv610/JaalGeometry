#include "SurfaceReconstruction.hpp"

////////////////////////////////////////////////////////////////////////////////

int JPoissonSurfaceReconstruction :: generate()
{
    newMesh.reset();
    if( mesh == nullptr) return 1;
    write_ply_file();

    ostringstream oss;
    oss << "PoissonRecon ";
    oss << " --in tmp1.ply ";
    oss << " --out tmp2.ply ";
    oss << " --bType " << boundaryType;
    oss << " --depth " << depth;
    oss << " --scale " << scale;
    oss << " --samplesPerNode " << samplesPerNode;
    oss << " --pointWeight " << pointWeight;
    string cmd = oss.str();
    system( cmd.c_str() );
    system ("mesh_filter tmp2.ply -scale 1.0 tmp2.off");
    JMeshIO mio;
    newMesh = mio.readFile("tmp2.off");
}

////////////////////////////////////////////////////////////////////////////////

void JPoissonSurfaceReconstruction :: write_ply_file()
{
    if( mesh == nullptr) return;
    ofstream ofile("tmp1.ply");
    ofile << "ply" << endl;
    ofile << "format ascii 1.0" << endl;
    ofile << "comment Jaal generated" << endl;
    ofile << "element vertex " << mesh->getSize(0) << endl;
    ofile << "property float x" << endl;
    ofile << "property float y" << endl;
    ofile << "property float z" << endl;
    ofile << "property float nx" << endl;
    ofile << "property float ny" << endl;
    ofile << "property float nz" << endl;
    ofile << "element face 0" << endl;
    ofile << "property list uchar int vertex_indices" << endl;
    ofile << "end_header" << endl;

    mesh->getGeometry()->setNodesNormal();

    size_t numNodes = mesh->getSize(0);
    Vec3F normal;
    Point3D xyz;

    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        xyz = vtx->getXYZCoords();
        vtx->getAttribute("Normal", normal);
        ofile << xyz[0] << " " << xyz[1] << " " << xyz[2] << " "
              << normal[0] << " " << normal[1] << " " << normal[2] << endl;
    }

}
////////////////////////////////////////////////////////////////////////////////

