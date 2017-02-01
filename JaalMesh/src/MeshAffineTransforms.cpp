#include "MeshAffineTransforms.hpp"

////////////////////////////////////////////////////////////////////////////////

Point3D JMeshAffineTransform :: getCenter( const JNodeSequence &nodes) const
{
    Point3D  center;
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;
    for( const JNodePtr &vtx: nodes) {
        const Point3D &xyz = vtx->getXYZCoords();
        center[0] += xyz[0];
        center[1] += xyz[1];
        center[2] += xyz[2];
    }
    center[0] /= (double) nodes.size();
    center[1] /= (double) nodes.size();
    center[2] /= (double) nodes.size();
    return center;
}
////////////////////////////////////////////////////////////////////////////////

int JMeshAffineTransform :: getConformalTriangle( const JFacePtr &face, Point3D &pa, Point3D &pb, Point3D &pc)
{
    if( face->getSize(0) != 3 ) return 1;

    const Point3D &p0  = face->getNodeAt(0)->getXYZCoords();
    const Point3D &p1  = face->getNodeAt(1)->getXYZCoords();
    const Point3D &p2  = face->getNodeAt(2)->getXYZCoords();

    pa[0] = 0.0;
    pa[1] = 0.0;
    pa[2] = 0.0;

    double len  = JMath::length(p0, p1);
    pb[0] = len;
    pb[1] = 0.0;
    pb[2] = 0.0;

    len  = JMath::length(p0, p2);
    double angle = JMath::getTriAngle(p0, p1, p2, ANGLE_IN_RADIANS);

    pc[0] = len*cos(angle);
    pc[1] = len*sin(angle);
    pc[2] = 0.0;

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: toCenter()
{
    if( mesh == nullptr ) return;
    mesh->getLogger()->setInfo("Bringing mesh to center ");

    Point3D p3d = mesh->getGeometry()->getCenter();
    translate( -p3d[0], -p3d[1], -p3d[2] );
}

void JMeshAffineTransform :: setCenter( double xc, double yc, double zc)
{
   if( mesh == nullptr) return;

   Point3D p3d = mesh->getGeometry()->getCenter();
   translate( xc-p3d[0], yc-p3d[1], zc-p3d[2] );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: translate( const JNodeSequence &nodes, double x, double y, double z )
{
    for(const JNodePtr &v: nodes) {
        Point3D xyz = v->getXYZCoords();
        xyz[0] += x;
        xyz[1] += y;
        xyz[2] += z;
        v->setXYZCoords( xyz );
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: translate( double x, double y, double z )
{
    if( mesh == nullptr) return;
    mesh->getLogger()->setInfo("Translating mesh");

    JNodeSequence nodes = mesh->getNodes();
    translate( nodes, x, y, z );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: scale( const JNodeSequence &nodes, double x, double y, double z )
{
    Point3D center = getCenter( nodes);
    translate( nodes, -center[0], -center[1], -center[2] );

    for( const JNodePtr &v: nodes) {
        Point3D xyz = v->getXYZCoords();
        xyz[0] *= x;
        xyz[1] *= y;
        xyz[2] *= z;
        v->setXYZCoords( xyz );
    }

    translate( nodes, center[0], center[1], center[2] );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: scale( double x, double y, double z )
{
    if( mesh == nullptr) return;
    mesh->getLogger()->setInfo("Scaling mesh");

    JNodeSequence nodes = mesh->getNodes();
    scale( nodes, x, y, z );
}
///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: xRotateMatrix( double theta, double mat[3][3])
{
    mat[0][0]  =  1.0;
    mat[0][1]  =  0.0;
    mat[0][2]  =  0.0;

    mat[1][0]  =  0.0;
    mat[1][1]  =  cos( theta );
    mat[1][2]  = -sin( theta );

    mat[2][0]  = 0.0;
    mat[2][1]  = sin( theta );
    mat[2][2]  = cos( theta );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: xRotatePoint( Point3D &p, double theta)
{
    double mat[3][3];
    xRotateMatrix( theta, mat);
    double x = mat[0][0]*p[0] + mat[0][1]*p[1] + mat[0][2]*p[2];
    double y = mat[1][0]*p[0] + mat[1][1]*p[1] + mat[1][2]*p[2];
    double z = mat[2][0]*p[0] + mat[2][1]*p[1] + mat[2][2]*p[2];
    p[0] = x;
    p[1] = y;
    p[2] = z;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: yRotateMatrix( double theta, double mat[3][3])
{
    mat[0][0]  =  cos(theta);
    mat[0][1]  =  0.0;
    mat[0][2]  =  sin(theta);

    mat[1][0]  = 0.0;
    mat[1][1]  = 1.0;
    mat[1][2]  = 0.0;

    mat[2][0]  = -sin(theta);
    mat[2][1]  =  0.0;
    mat[2][2]  =  cos( theta );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: yRotatePoint( Point3D &p, double theta)
{
    double mat[3][3];
    yRotateMatrix( theta, mat);
    double x = mat[0][0]*p[0] + mat[0][1]*p[1] + mat[0][2]*p[2];
    double y = mat[1][0]*p[0] + mat[1][1]*p[1] + mat[1][2]*p[2];
    double z = mat[2][0]*p[0] + mat[2][1]*p[1] + mat[2][2]*p[2];
    p[0] = x;
    p[1] = y;
    p[2] = z;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: zRotateMatrix( double theta, double mat[3][3])
{
    mat[0][0]  =  cos(theta);
    mat[0][1]  = -sin(theta);
    mat[0][2]  =  0.0;

    mat[1][0]  =  sin(theta);
    mat[1][1]  =  cos(theta);
    mat[1][2]  =  0.0;

    mat[2][0]  =  0.0;
    mat[2][1]  =  0.0;
    mat[2][2]  =  1.0;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: zRotatePoint(Point3D &p, double theta)
{
    double mat[3][3];
    zRotateMatrix( theta, mat);
    double x = mat[0][0]*p[0] + mat[0][1]*p[1] + mat[0][2]*p[2];
    double y = mat[1][0]*p[0] + mat[1][1]*p[1] + mat[1][2]*p[2];
    double z = mat[2][0]*p[0] + mat[2][1]*p[1] + mat[2][2]*p[2];

    p[0] = x;
    p[1] = y;
    p[2] = z;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: rotate( const JNodeSequence &nodes, double theta, int axis)
{
    double m[3][3];
    switch(axis) {
    case 0:
        xRotateMatrix( theta, m);
        break;
    case 1:
        yRotateMatrix( theta, m);
        break;
    case 2:
        zRotateMatrix( theta, m);
        break;
    }

    Point3D p3d;
    size_t nSize = nodes.size();
    for( size_t i = 0; i < nSize; i++) {
        const JNodePtr &v = nodes[i];
        const Point3D xyz = v->getXYZCoords();
        p3d[0] = m[0][0]*xyz[0] + m[0][1]*xyz[1] + m[0][2]*xyz[2];
        p3d[1] = m[1][0]*xyz[0] + m[1][1]*xyz[1] + m[1][2]*xyz[2];
        p3d[2] = m[2][0]*xyz[0] + m[2][1]*xyz[1] + m[2][2]*xyz[2];
        v->setXYZCoords( p3d );
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: rotate( double theta, int axis)
{
    if( mesh == nullptr) return;
    mesh->getLogger()->setInfo("Scaling mesh");

    JNodeSequence nodes = mesh->getNodes();
    rotate( nodes, theta, axis);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: getRotationMatrix( double A[3][3], double B[3][3], double C[3][3] )
{
    for( int i = 0; i < 3; i++) {
        C[i][0] = 0.0;
        C[i][1] = 0.0;
        C[i][2] = 0.0;
    }

    for( int i = 0; i < 3; i++) {
        for( int j = 0; j < 3; j++) {
            double sum = 0.0;
            for( int k = 0; k <3; k++)
                sum = sum + A[i][k]*B[k][j];
            C[i][j] = sum;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: getRotationMatrix( double A[3][3], double B[3][3], double C[3][3], double D[3][3] )
{
    double AB[3][3];
    getRotationMatrix(A, B, AB );
    getRotationMatrix(AB, C, D );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform :: alignAlong(JNodeSequence &nodes, const Vec3D &srcVec, const Vec3D &dstVec)
{
    Vec3D  perpAxis;
    JMath::cross_product( srcVec, dstVec, perpAxis);
    double dl = JMath::magnitude(perpAxis);
    perpAxis[0] /= dl;
    perpAxis[1] /= dl;
    perpAxis[2] /= dl;
    double angle = JMath::getVecAngle(srcVec, dstVec, ANGLE_IN_RADIANS);
    if( fabs(angle) < 1.0E-15) return;

    double qcos = cos(0.5*angle);
    double qsin = sin(0.5*angle);

    boost::math::quaternion<double> q(qcos, qsin*perpAxis[0], qsin*perpAxis[1], qsin*perpAxis[2]);
    boost::math::quaternion<double> q1 = boost::math::conj(q);

    boost::math::quaternion<double> result;
    size_t numNodes = mesh->getSize(0);
    for( const JNodePtr &vtx : nodes) {
        Point3D p3d = vtx->getXYZCoords();
        boost::math::quaternion<double> v(0.0, p3d[0], p3d[1], p3d[2]);
        result = q*v*q1;
        p3d[0]  = result.R_component_2();
        p3d[1]  = result.R_component_3();
        p3d[2]  = result.R_component_4();
        vtx->setXYZCoords(p3d);
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransform :: alignAlong(const Vec3D &srcVec, const Vec3D &dstVec)
{
    if( mesh == nullptr) return;
    JNodeSequence nodes = mesh->getNodes();
    alignAlong(nodes, srcVec, dstVec);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransform  :: normalize()
{
    if( mesh == nullptr ) return;

    size_t numnodes = mesh->getSize(0);

    if( numnodes == 0) return;

    mesh->getLogger()->setInfo("Normalizing the mesh" );

    double xmin, xmax, ymin, ymax, zmin, zmax;
    Point3D xyz;
    xyz = mesh->getNodeAt(0)->getXYZCoords();

    xmin = xyz[0];
    xmax = xyz[0];
    ymin = xyz[1];
    ymax = xyz[1];
    zmin = xyz[2];
    zmax = xyz[2];

    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            xyz = v->getXYZCoords();
            xmin = min(xmin, xyz[0]);
            xmax = max(xmax, xyz[0]);
            ymin = min(ymin, xyz[1]);
            ymax = max(ymax, xyz[1]);
            zmin = min(zmin, xyz[2]);
            zmax = max(zmax, xyz[2]);
        }
    }

    double xlen = fabs(xmax - xmin);
    double ylen = fabs(ymax - ymin);
    double zlen = fabs(zmax - zmin);
    double scale = max(max(xlen, ylen), zlen);

    #pragma omp parallel for private(xyz)
    for (size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        if( v->isActive() ) {
            xyz = v->getXYZCoords();
            xyz[0] /= scale;
            xyz[1] /= scale;
            xyz[2] /= scale;
            v->setXYZCoords(xyz);
        }
    }
}

