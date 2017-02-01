#include <stdio.h>

#include "Mesh.hpp"
#include "AllTriMeshGenerator.hpp"
#include "AffineTransforms.hpp"




int generate()
{
    JMeshOFFExporter mexp;

    JMeshPtr  mesh = JMesh::newObject();
    ifstream infile("Bracket3Hole_data.txt", ios::in);
    if( infile.fail() ) return 1;

    int numnodes;
    infile >> numnodes;

    system("mesh_make ssphere 4 2 sph.off");


    JMeshOFFImporter mimp;
    JMeshPtr sphtemp = mimp.readFile("sph.off");
    JAffineTransform affine;
    affine.setMesh(sphtemp);
    affine.scale(0.0005, 0.0005, 0.0005);

    JMeshPtr sphmesh = JMesh::newObject();

    JNodeSequence nodes;
    Point3D p3d;
    int id;
    for( int i = 0; i < numnodes; i++) {
        cout << "I = " << i << endl;
        infile >> id >> p3d[0] >> p3d[1] >> p3d[2];
        JNodePtr vtx = JVertex::newObject();
        vtx->setXYZCoords(p3d);
        nodes.push_back(vtx);
        JMeshPtr sp = sphtemp->deepCopy();
        affine.setMesh(sp);
        affine.translate(p3d[0], p3d[1], p3d[2]);
        sphmesh->addObject(sp);
    }
    mesh->addObjects(nodes);
    mexp.writeFile(sphmesh, "ball.off");

    /*
       int numEdges;
       infile >> numEdges;
       double    radius;
       int       id0, id1;
       for( int i = 0; i < numEdges; i++){
            infile >> id0 >> id1 >> radius;
            JNodePtr &v0 = nodes[id0];
            JNodePtr &v1 = nodes[id1];
            JEdgePtr e   = JEdge::newObject(v0,v1);
            mesh->addObject(e);
       }

       for( int i = 0; i < numEdges; i++){
            cout << "Edge Mesh " << i << endl;
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            const JNodePtr &v0 = edge->getNodeAt(0);
            const JNodePtr &v1 = edge->getNodeAt(1);
            JMeshPtr edgemesh = AllTriMeshGenerator::getCylinder(v0->getXYZCoords(), v1->getXYZCoords(), radius, 10);
            mesh->addObject(edgemesh);
       }
       mexp.writeFile(mesh, "microstruct.off");
    */
    return 0;
}

int main()
{
    JHexahedronPtr   hex = Jaal::Hexahedron::getCanonical();
    cout << "HELLO " << endl;

    JMeshQuality mq;
    double val;
    mq.getJacobian(hex,val);
    cout << val << endl;
    vector<JCellPtr> tets;
    hex->getTetrahedra(tets);
    for( int i = 0; i < tets.size(); i++) {
        mq.getJacobian(tets[i],val);
        cout << val << endl;
    }
    JNodeSequence nodes;
    nodes = hex->getNodes();
    JMeshPtr m = JMesh::newObject();
    m->addObjects(nodes);
    m->addObjects(tets);
    JMeshXMLExporter mex;
    mex.writeFile(m, "tets.xml");




    /*
        int N = 100000;
        vector<double>  v(N);
        for( int i = 0; i < N; i++) v[i] = drand48();


        double maxval = -1.0;
        for( int i = 0; i < N; i++)
             if( v[i] > maxval ) maxval = max( maxval, v[i]);
        cout << maxval << endl;


        maxval = -1.0;
    #pragma omp parallel for reduction(max:maxval)
        for( int i = 0; i < N; i++) {
             if( v[i] > maxval ) maxval = max( maxval, v[i]);
        }
        cout << maxval << endl;
        cout << omp_get_max_threads() << endl;
    */

//   generate();
}


