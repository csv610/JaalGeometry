#include "AllHexMeshGenerator.hpp"

using namespace std;
using namespace Jaal;

JMeshPtr AllHexMeshGenerator :: getGeodeTemplate()
{

    ifstream infile( "geode.txt", ios::in);

    int id;

    JMeshPtr mesh = JMesh::newObject();
    Point3D p3d;
    for( int i = 0; i < 48; i++) {
        JNodePtr vtx = JNode::newObject();
        infile >> id >> p3d[0] >>  p3d[1] >> p3d[2];
        vtx->setXYZCoords(p3d);
        mesh->addObject(vtx);
    }

    JNodeSequence hexnodes(8);
    for( int i = 0; i < 26; i++) {
        infile >> id;
        for( int j = 0; j < 8; j++) {
            infile >> id;
            hexnodes[j] = mesh->getNodeAt(id-1);
        }
        JCellPtr c = JHexahedron::newObject();
        c->setNodes( hexnodes );
        mesh->addObject(c);
    }
    return mesh;
}
