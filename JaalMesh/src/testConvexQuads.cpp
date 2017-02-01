#include "Mesh.hpp"
#include"MeshOptimization.hpp"
#include "MeshExporter.hpp"

using namespace Jaal;

int meshtype     = 4;
int randomtangle = 1;


///////////////////////////////////////////////////////////////////////////////
int CountConcave( Mesh *mesh)
{
    int numfaces = mesh->getSize(2);
    int nCount = 0;
    for( int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int  pos =  FaceGeometry::concaveAt(face);
        if( pos >= 0) nCount++;
    }
    return nCount;
}
///////////////////////////////////////////////////////////////////////////////

void ConvexQuads( Mesh *mesh)
{
    int numfaces = mesh->getSize(2);

    Point3D pmid;
    for( int i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        int  pos =  FaceGeometry::concaveAt(face);
        if( pos >= 0) {
            const Point3D  &p1 = face->getNodeAt(pos+1)->getXYZCoords();
            const Point3D  &p2 = face->getNodeAt(pos+3)->getXYZCoords();
            pmid[0] = 0.5*(p1[0]+p2[0]);
            pmid[1] = 0.5*(p1[1]+p2[1]);
            face->getNodeAt(pos)->setXYZCoords(pmid);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void fromTriangles()
{
    int  xlength = 1.1;
    int  ylength = 1.1;

    ofstream ofile( "test.poly", ios::out);

    int numPoints = 100;

    ofile << numPoints+4 << " 2  0  0" << endl;
    ofile << " 0  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << " 1  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << " 2  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << " 3  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << "4  1 " << endl;
    ofile << " 0  0  1  1 "  << endl;
    ofile << " 1  1  2  2 "  << endl;
    ofile << " 2  2  3  3 "  << endl;
    ofile << " 3  3  0  4 "  << endl;
    ofile << " 0 " << endl;
    ofile.close();

    system( "triangle -peq30a0.1 test.poly");

    Mesh *mesh = Mesh::newObject();
    mesh->readFromFile( "test.1.ele");

    if( meshtype == 4) {
        Mesh *quadmesh= AllQuadMeshGenerator::SimpleTris2Quads(mesh);
        mesh = quadmesh;
    }
    mesh->saveAs("tmp.vtk");
}

void fromQuads()
{
//  CSV
    int grid_dim[]   = {4, 4};
    double length[]  = {1, 1};
    double origin[]  = {0.0, 0.0};
    origin[0] = -0.5*length[0];
    origin[1] = -0.5*length[1];

    Mesh *quadmesh= getStructuredMesh(2, grid_dim, length, origin);
    Mesh *mesh = quadmesh;

    if( meshtype == 3) {
        Mesh *trimesh = AllTriMeshGenerator::getFromQuadMesh(quadmesh,4);
        mesh = trimesh;
    }

    JMeshNonlinearOptimization mopt;
    mopt.setMesh(mesh);
    mopt.setBoundaryPreservation(1);

    int numnodes = mesh->getSize(0);

    int val = 2;
    int index = 0;
    if( randomtangle) {
        double xc = 0.0;
        double yc = 0.0;
        double rcut = 0.1;
        int nCount = 0;
        for( int i = 0; i < numnodes; i++) {
            Vertex *v = mesh->getNodeAt(i);
            Point3D xyz = v->getXYZCoords();
            double dx = xyz[0] - xc;
            double dy = xyz[1] - yc;
            double r = sqrt(dx*dx + dy*dy);
            if( r < rcut)  {
                xyz[0] = xc + (-rcut + 2.0*rcut*drand48());
                xyz[1] = yc + (-rcut + 2.0*rcut*drand48());
                v->setXYZCoords(xyz);
                nCount++;
            }
        }
    }

    cout << "Info: Storing mesh in tmp.vtk" << endl;
    mesh->saveAs("tmp.vtk");

    /*
        for( int i = 0; i < 10; i++) {
        ConvexQuads(mesh);
        cout << "#Concave elements " << CountConcave( mesh) << endl;
        }
        mesh->saveAs("tmp1.vtk");
        exit(0);
    */


//  Assign random Z value to the coordinates;
    for( int i = 0; i < numnodes; i++) {
        Vertex *v = mesh->getNodeAt(i);
        Point3D &xyz = v->getXYZCoords();
        xyz[2] = drand48();
        v->setXYZCoords(xyz);
    }

    JMeshTRIExporter mexp;
    mexp.setDimension(3);
    mexp.writeNodes( mesh, "rand.node");

    system("tetgen rand.node");
    JMeshTRIImporter mimp;

    Mesh *tetmesh = Mesh::newObject();
    mimp.readNodes("rand.1.node", tetmesh);
    mimp.readCells("rand.1.ele",  tetmesh);

    cout << "Delaunay Tets in tmp1.vtk" << endl;
    tetmesh->saveAs("tmp1.vtk");

    Mesh *projmesh = Mesh::newObject();

    assert( tetmesh->getSize(0) == numnodes);
    for( int i = 0; i < numnodes; i++) {
        Vertex *v = tetmesh->getNodeAt(i);
        Point3D &xyz = v->getXYZCoords();
        xyz[2] = 0.0;
        v->setXYZCoords(xyz);
        projmesh->addObject(v);
    }
    cout << "Projection in tmp2.vtk" << endl;

    int nCells = tetmesh->getSize(3);
    JNodeSequence tetnodes;
    for( int i = 0; i < numnodes; i++) {
        Cell *c = tetmesh->getCellAt(i);
        tetnodes = c->getNodes();
        assert( tetnodes.size() == 4 );

        Face *f = Quadrilateral::newObject(tetnodes[0], tetnodes[1], tetnodes[3], tetnodes[2]);
        projmesh->addObject(f);
    }
    projmesh->saveAs("tmp2.vtk");
}

int main()
{

    fromQuads();
    return 0;
}
