#pragma once

#include "Mesh.hpp"
#include <boost/math/quaternion.hpp>

class JMeshAffineTransform {
public:

    static void xRotateMatrix( double theta, double m[3][3] );
    static void yRotateMatrix( double theta, double m[3][3] );
    static void zRotateMatrix( double theta, double m[3][3] );

    static void xRotatePoint(Point3D &p, double theta);
    static void yRotatePoint(Point3D &p, double theta);
    static void zRotatePoint(Point3D &p, double theta);

    static void getRotationMatrix( double A[3][3], double B[3][3], double C[3][3]);
    static void getRotationMatrix( double A[3][3], double B[3][3], double C[3][3], double D[3][3]);
    static int  getConformalTriangle(const JFacePtr &face, Point3D &pa, Point3D &pb, Point3D &pc);

    JMeshAffineTransform()
    {
        mesh = nullptr;
    }

    JMeshAffineTransform( const JMeshPtr &m)
    {
        mesh = m;
    }

    void setMesh( const JMeshPtr &m) {
        mesh = m;
    }

    void toCenter();

    void rotate( const JNodeSequence &nodes, double angle, int axis );
    void scale( const JNodeSequence &nodes, double x, double y, double z );
    void translate( const JNodeSequence &nodes, double x, double y, double z );
    void alignAlong( JNodeSequence &nodes, const Point3D &vsrc, const Point3D &vdst);

    void rotate( double angle, int axis );
    void scale( double x, double y, double z );
    void translate( double x, double y, double z );
    void alignAlong( const Point3D &vsrc, const Point3D &vdst);
    void normalize();
    void setCenter( double xc, double yc, double zc );

    void apply( double m[3][3] );
//  void alignAlong( const Point3D &vsrc, const Point3D &vdst);

private:
    JMeshPtr mesh;

    void xrotate( double t );
    void yrotate( double t );
    void zrotate( double t );
    Point3D  getCenter( const JNodeSequence &n) const;
};
