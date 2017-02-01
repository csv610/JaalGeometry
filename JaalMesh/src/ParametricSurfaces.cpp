#include "Mesh.hpp"

///////////////////////////////////////////////////////////////////////////////

Mesh *ParametricSurface :: getCrossCapSurface(int nx, int ny)
{
    Mesh *mesh = Mesh::newObject();

    double du = M_PI/(double)nx;
    double dv = M_PI/(double)ny;

    Point3D xyz;
    for( int j = 0; j < ny; j++) {
        for( int i = 0; i < nx;  i++) {
            Vertex *vtx  = Vertex::newObject();
            double   u = i*du;
            double   v = j*dv;
            xyz[0]   = sin(u)*sin(2*v)/2;
            xyz[1]   = sin(2*u)*cos(v)*cos(v);
            xyz[2]   = cos(2*u)*cos(v)*cos(v);
            vtx->setXYZCoords( xyz );
            mesh->addObject(vtx);
        }
    }
    return mesh;
}

///////////////////////////////////////////////////////////////////////////////

Mesh *ParametricSurface :: getAperyBoySurface(int nx, int ny)
{
    Mesh *mesh =  Mesh::newObject();

    double du = M_PI/(double)nx;
    double dv = M_PI/(double)ny;
    double alpha = 1.0;

    Point3D xyz;
    for( int j = 0; j < ny; j++) {
        for( int i = 0; i < nx;  i++) {
            Vertex *vtx  = Vertex::newObject();
            double   u = i*du;
            double   v = j*dv;
            xyz[0] = (2/3.0)*cos(u)*(cos(u)*cos(2*v) + sqrt(2)*sin(u)*cos(v))/
                     (sqrt(2) - sin(2*u)*sin(3*v));
            xyz[1] = (2/3.0)*cos(u)*(cos(u)*sin(2*v) - sqrt(2)*sin(u)*sin(v))/
                     (sqrt(2)  - sin(2*u)*cos(u)*sin(3*v));
            xyz[2] = sqrt(2)*(cos(u)*cos(u))/(sqrt(2) - alpha*sin(2*u)*sin(3*v));
            vtx->setXYZCoords( xyz );
            mesh->addObject(vtx);
        }
    }
    return mesh;
}

///////////////////////////////////////////////////////////////////////////////

Mesh *ParametricSurface :: getBoySurface(int nx, int ny)
{
    Mesh *mesh = Mesh::newObject();

    double du = 1.0/(double)nx;
    double dv = 1.0/(double)ny;

    Point3D xyz;
    double g1, g2, g3, ga;
    std::complex<double> f1, f2, f3, z,d, one(1.0, 0.0);

    for( int j = 0; j < ny; j++) {
        for( int i = 0; i < nx;  i++) {
            Vertex *vtx  = Vertex::newObject();
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
                mesh->addObject(vtx);
            }
        }
    }
    return mesh;
}

