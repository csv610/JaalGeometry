#include <iostream>
#include <string>

#include "Mesh.hpp"
#include "HarmonicField.hpp"

#include "MeshExporter.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "AllHexMeshGenerator.hpp"
#include "AffineTransforms.hpp"
#include <stdlib.h>
#include <boost/range/algorithm.hpp>

using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    vector<int> a(100), b(50);

    for( int i = 0; i < 100; i++) {
        a[i] = 100-i;
    }

    boost::sort(a);

    for( int i = 0; i < 50; i++) {
        b[i] = 75+i;
    }
    boost::sort(b);

    vector<int> c;
    boost::set_difference(a, b, inserter(c, c.end()) );

    for( int i = 0; i < c.size(); i++)
        cout << i << " " << c[i] << endl;
    exit(0);


    JMeshIO mexp;
//  JMeshPtr mesh = mexp.readFile("sph.off");

    int dim[3] = {10,10,10};
//  JMeshPtr mesh = AllHexMeshGenerator::getStructuredMesh(dim);
    JMeshPtr mesh = AllQuadMeshGenerator::getStructuredMesh(dim);
    /*
       JAffineTransform affine;
       srand(time(0));
       JMeshPtr catmesh = JMesh::newObject();
       for( int i = 0; i < 67; i++) {
            JMeshPtr dupm = mesh->deepCopy();
            affine.setMesh(dupm);
            affine.translate(i,i, 0);
            catmesh->addObject(dupm);
       }
       mexp.saveAs(catmesh, "B.xml");
    */
    mexp.saveAs(mesh, "A.xml");


    /*
        mesh->getTopology()->searchBoundary();

        Point3D xyz;
        size_t numnodes = mesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
             const JNodePtr &vtx = mesh->getNodeAt(i);
             if( !vtx->isBoundary() ) {
                  xyz = vtx->getXYZCoords();
                  xyz[0] = xyz[0] + 1.5*drand48();
                  xyz[1] = xyz[1] + 1.5*drand48();
                  xyz[2] = xyz[2] + 1.5*drand48();
                  vtx->setXYZCoords(xyz);
             }
        }
        mexp.saveAs(mesh, "A.xml");

        cout << "Solver " << endl;
        int solver;
        cin >> solver;

        JHarmonicField hf;
        hf.setSolver( solver );
        hf.setMesh(mesh);
        hf.solveSystem();
        mexp.saveAs(mesh, "B.xml");
    */
}

///////////////////////////////////////////////////////////////////////////////
