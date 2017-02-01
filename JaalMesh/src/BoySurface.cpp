#include "BoySurface.hpp"

///////////////////////////////////////////////////////////////////////////////
void JBoySurface:: getQuadMesh(int nx, int ny)
{
    for( int j = 0; j < ny-1; j++) {
        for( int i = 0; i < nx-1;  i++) {
            JNodePtr v0 = newmesh->getNodeAt( j*nx + i );
            JNodePtr v1 = newmesh->getNodeAt( j*nx + i + 1);
            JNodePtr v2 = newmesh->getNodeAt( (j+1)*nx + i + 1);
            JNodePtr v3 = newmesh->getNodeAt( (j+1)*nx + i );
            JFacePtr face = JQuadrilateral::newObject( v0, v1, v2, v3 );
            newmesh->addObject( face );
        }
    }

    /*
         for( int j = 0; j < ny-1; j++) {
              Vertex *v0 = newmesh->getNodeAt( j*nx + (nx-1) );
              Vertex *v1 = newmesh->getNodeAt( j*nx );
              Vertex *v2 = newmesh->getNodeAt( (j+1)*nx);
              Vertex *v3 = newmesh->getNodeAt( (j+1)*nx + (nx-1) );
              Face *face = Quadrilateral::newObject( v0, v1, v2, v3 );
              newmesh->addObject( face );
         }

         for( int i = 0; i < nx-1;  i++) {
              Vertex *v0 = newmesh->getNodeAt( (ny-1)*nx + i );
              Vertex *v1 = newmesh->getNodeAt( (ny-1)*nx + i + 1);
              Vertex *v2 = newmesh->getNodeAt( i + 1);
              Vertex *v3 = newmesh->getNodeAt( i );
              Face *face = Quadrilateral::newObject( v0, v1, v2, v3 );
              newmesh->addObject( face );
         }
    */
}
///////////////////////////////////////////////////////////////////////////////

JMeshPtr JBoySurface:: getApery(int nx, int ny)
{
    newmesh =  JMesh::newObject();

    double du = M_PI/(double)(nx-1);
    double dv = M_PI/(double)(ny-1);
    double alpha = 1.0;

    Point3D xyz;
    for( int j = 0; j < ny; j++) {
        for( int i = 0; i < nx;  i++) {
            JNodePtr vtx  = JNode::newObject();
            double   u = i*du;
            double   v = j*dv;
            xyz[0] = (2/3.0)*cos(u)*(cos(u)*cos(2*v) + sqrt(2)*sin(u)*cos(v))/
                     (sqrt(2) - sin(2*u)*sin(3*v));
            xyz[1] = (2/3.0)*cos(u)*(cos(u)*sin(2*v) - sqrt(2)*sin(u)*sin(v))/
                     (sqrt(2)  - sin(2*u)*cos(u)*sin(3*v));
            xyz[2] = sqrt(2)*(cos(u)*cos(u))/(sqrt(2) - alpha*sin(2*u)*sin(3*v));
            vtx->setXYZCoords( xyz );
            newmesh->addObject(vtx);
        }
    }
    getQuadMesh(nx,ny);

    return newmesh;
}

///////////////////////////////////////////////////////////////////////////////

JMeshPtr JBoySurface :: getBryant(int nx, int ny)
{
    newmesh = JMesh::newObject();

    double du = 1.0/(double)(nx-1);
    double dv = 1.0/(double)(ny-1);

    Point3D xyz;
    double g1, g2, g3, ga;
    std::complex<double> f1, f2, f3, z,d, one(1.0, 0.0);

    for( int j = 0; j < ny; j++) {
        for( int i = 0; i < nx;  i++) {
            JNodePtr vtx  = JNode::newObject();
            double   u = i*du;
            double   v = j*dv;
            z = std::complex<double>(u, v);
            if( abs(z) < 1.0 ) {
                d = pow(z,6) + sqrt(5.0)*pow(z,3) - one ;
                f1 = z*( one - pow(z,4) )/d;
                f2 = z*( one + pow(z,4) )/d;
                f3 = (one + pow(z,6))/d;
                g1 = -1.5*f1.imag();
                g2 = -1.5*f2.real();
                g3 = f3.imag() - 0.5;
                ga = g1*g1 + g2*g2 + g3*g3;
                xyz[0] =  g1/ga;
                xyz[1] =  g2/ga;
                xyz[2] =  g3/ga;
                vtx->setXYZCoords( xyz );
                newmesh->addObject(vtx);
            }
        }
    }

    getQuadMesh(nx,ny);

    return newmesh;
}

///////////////////////////////////////////////////////////////////////////////


void JBoySurface :: createNodes()
{
    nodes.resize(9);
    Point3D xyz;

    xyz[0] =  -2.0;
    xyz[1] =  0.0;
    xyz[2] =  0.0;
    nodes[0] = JNode::newObject();
    nodes[0]->setXYZCoords(  xyz );

    xyz[0] =   0.0;
    xyz[1] = -2.0;
    xyz[2] =  0.0;
    nodes[1] = JNode::newObject();
    nodes[1]->setXYZCoords(  xyz );

    xyz[0] =  0.0;
    xyz[1] = 0.0;
    xyz[2] = -2.0;
    nodes[2] = JNode::newObject();
    nodes[2]->setXYZCoords(  xyz );

    xyz[0] = -1.0;
    xyz[1] = 2.0;
    xyz[2] =  1.0;
    nodes[3] = JNode::newObject();
    nodes[3]->setXYZCoords(  xyz );

    xyz[0] =  1.0;
    xyz[1] = -1.0;
    xyz[2] =  2.0;
    nodes[4] = JNode::newObject();
    nodes[4]->setXYZCoords(  xyz );

    xyz[0] =  2.0;
    xyz[1] =  1.0;
    xyz[2] = -1.0;
    nodes[5] = JNode::newObject();
    nodes[5]->setXYZCoords(  xyz );

    xyz[0] = -1.0;
    xyz[1] =  1.0;
    xyz[2] = 0.0;
    nodes[6] = JNode::newObject();
    nodes[6]->setXYZCoords(  xyz );

    xyz[0] =  0.0;
    xyz[1] = -1.0;
    xyz[2] = 1.0;
    nodes[7] = JNode::newObject();
    nodes[7]->setXYZCoords(  xyz );

    xyz[0] =  1.0;
    xyz[1] =  0.0;
    xyz[2] =-1.0;
    nodes[8] = JNode::newObject();
    nodes[8]->setXYZCoords(  xyz );
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JBoySurface :: getMin1()
{
    createNodes();
    JMeshPtr mesh = JMesh::newObject();
    mesh->addObjects( nodes );

    // a : 0  b : 1  c : 2  d : 3
    // e : 4  f : 5  g : 6  h : 7 i : 8

    // [agfi]
    JFacePtr q1 = JQuadrilateral::newObject( nodes[0], nodes[6], nodes[5], nodes[8] );

    // [cieh]
    JFacePtr q2 = JQuadrilateral::newObject( nodes[2], nodes[8], nodes[4], nodes[7] );

    // [bhdg]
    JFacePtr q3 = JQuadrilateral::newObject( nodes[1], nodes[7], nodes[3], nodes[6] );

    // [abhe]
    JFacePtr q4 = JQuadrilateral::newObject( nodes[0], nodes[1], nodes[7], nodes[4] );

    // [bfic]
    JFacePtr q5 = JQuadrilateral::newObject( nodes[1], nodes[5], nodes[8], nodes[2] );

    // [acdg]
    JFacePtr q6 = JQuadrilateral::newObject( nodes[0], nodes[2], nodes[3], nodes[6] );

    // [aie]
    JFacePtr t1 = JTriangle::newObject( nodes[0], nodes[8], nodes[4] );

    // [bgf]
    JFacePtr t2 = JTriangle::newObject( nodes[1], nodes[6], nodes[5] );

    // [abc]
    JFacePtr t3 = JTriangle::newObject( nodes[0], nodes[1], nodes[2] );

    // [chd]
    JFacePtr t4 = JTriangle::newObject( nodes[2], nodes[7], nodes[3] );

    mesh->addObject(q1);
    mesh->addObject(q2);
    mesh->addObject(q3);
    mesh->addObject(q4);
    mesh->addObject(q5);
    mesh->addObject(q6);
    mesh->addObject(t1);
    mesh->addObject(t2);
    mesh->addObject(t3);
    mesh->addObject(t4);

    return mesh;
}

////////////////////////////////////////////////////////////////////////////////

JMeshPtr JBoySurface :: getMin2()
{
    createNodes();
    JMeshPtr mesh = JMesh::newObject();
    mesh->addObjects( nodes );

    // a : 0  b : 1  c : 2  d : 3
    // e : 4  f : 5  g : 6  h : 7 i : 8

    // [aifgd]
    JFacePtr q1 = JQuadrilateral::newObject( nodes[0], nodes[6], nodes[5], nodes[8] );

    // [chieif]
    JFacePtr q2 = JQuadrilateral::newObject( nodes[2], nodes[8], nodes[4], nodes[7] );

    // [bgdhe]
    JFacePtr q3 = JQuadrilateral::newObject( nodes[1], nodes[7], nodes[3], nodes[6] );

    // [aei]
    JFacePtr q4 = JQuadrilateral::newObject( nodes[0], nodes[1], nodes[7], nodes[4] );

    // [abe]
    JFacePtr q5 = JQuadrilateral::newObject( nodes[1], nodes[5], nodes[8], nodes[2] );

    // [adc]
    JFacePtr q6 = JQuadrilateral::newObject( nodes[0], nodes[2], nodes[3], nodes[6] );

    // [acb]
    JFacePtr t1 = JTriangle::newObject( nodes[0], nodes[8], nodes[4] );

    // [bgf]
    JFacePtr t2 = JTriangle::newObject( nodes[1], nodes[6], nodes[5] );

    // [abc]
    JFacePtr t3 = JTriangle::newObject( nodes[0], nodes[1], nodes[2] );

    // [chd]
    JFacePtr t4 = JTriangle::newObject( nodes[2], nodes[7], nodes[3] );

    mesh->addObject(q1);
    mesh->addObject(q2);
    mesh->addObject(q3);
    mesh->addObject(q4);
    mesh->addObject(q5);
    mesh->addObject(q6);
    mesh->addObject(t1);
    mesh->addObject(t2);
    mesh->addObject(t3);
    mesh->addObject(t4);

    return mesh;
}
