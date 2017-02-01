#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

JMeshPtr AllQuadMeshGenerator::SchneiderPyramid()
{
    cout << "HELLO " << endl;
    JMeshPtr mesh = JMesh::newObject();

    JNodeSequence nodes(100);

    Point3D p3d;

    p3d[0] = -1.0;
    p3d[1] = -1.0;
    p3d[2] = 0.0;
    JNodePtr v0 = JNode::newObject();
    v0->setXYZCoords(p3d);
    mesh->addObject(v0);

    p3d[0] =  1.0;
    p3d[1] = -1.0;
    p3d[2] = 0.0;
    JNodePtr v1 = JNode::newObject();
    v1->setXYZCoords(p3d);
    mesh->addObject(v1);

    p3d[0] = 1.0;
    p3d[1] = 1.0;
    p3d[2] = 0.0;
    JNodePtr v2  = JNode::newObject();
    v2->setXYZCoords(p3d);
    mesh->addObject(v2);

    p3d[0] = -1.0;
    p3d[1] =  1.0;
    p3d[2] =  0.0;
    JNodePtr v3 = JNode::newObject();
    v3->setXYZCoords(p3d);
    mesh->addObject(v3);

    p3d[0] =  0.0;
    p3d[1] =  0.0;
    p3d[2] =  2.0;
    JNodePtr v4 = JNode::newObject();
    v4->setXYZCoords(p3d);
    mesh->addObject(v4);

    JFacePtr face;
    /*
       TriRefiner  triRefiner(mesh);
       QuadRefiner quadRefiner(mesh);

       // Base Quad
       face  = Quadrilateral::newObject( v0, v1, v2, v3);
       quadRefiner.refine14( face);

       // 1st triangle ...
       face  = Triangle::newObject( v0, v1, v4);
       triRefiner.tri2quads( face);
       delete face;

       // 2nd triangle ...
       face  = Triangle::newObject( v1, v2, v4);
       triRefiner.tri2quads( face);
       delete face;

       // 3rd triangle ...
       face  = Triangle::newObject( v2, v3, v4);
       triRefiner.tri2quads( face);
       delete face;

       // 4th triangle ...
       face  = Triangle::newObject( v3, v0, v4);
       triRefiner.tri2quads( face);

       mesh->delete_edge_attribute( "Steiner");
       mesh->pruneAll();
    */
    return mesh;
}


