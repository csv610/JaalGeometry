#include "Geometry.hpp"
#include "Mesh.hpp"

int JGeometry::getBoundedSide( const double *p0, const double *p1, const double *p2, const double *ptest)
{
    /*
        // Given three points ( forming the triangle), this function detemines if a test point "ptest" is within
        // the triangle, on the boundary or outside the triangle ....
        double area = TriGeometry::getArea2D(p0,p1,p2);

        double *pa, *pb, *pc;
        if( area > 0.0) {
            pa = const_cast<double*>(p0);
            pb = const_cast<double*>(p1);
            pc = const_cast<double*>(p2);
        } else {
            pa = const_cast<double*>(p0);
            pb = const_cast<double*>(p2);
            pc = const_cast<double*>(p1);
        }

        int ori;
        ori = getPointOrientation(pa,pb,ptest);
        if( ori == GeomOrient::CLOCKWISE) return GeomOrient::OUTSIDE;
        if( ori == GeomOrient::BOUNDARY ) return GeomOrient::BOUNDARY;

        ori = getPointOrientation(pb,pc,ptest);
        if( ori == GeomOrient::CLOCKWISE) return GeomOrient::OUTSIDE;
        if( ori == GeomOrient::BOUNDARY ) return GeomOrient::BOUNDARY;

        ori = getPointOrientation(pc,pa,ptest);
        if( ori == GeomOrient::CLOCKWISE) return GeomOrient::OUTSIDE;
        if( ori == GeomOrient::BOUNDARY ) return GeomOrient::BOUNDARY;

        return GeomOrient::INSIDE;
    */
    cout << "Exit" << __LINE__ << endl;
    exit(0);
}

////////////////////////////////////////////////////////////////////////////////
int JGeometry :: isInside( const vector<Point2D> &polyPoints,  Point2D &queryPoint)
{
    Point2D  extremePoint;
    extremePoint[0] = 0.90*numeric_limits<double>::max();
    extremePoint[1] = queryPoint[1];

    int nCount = 0;

    int n = polyPoints.size();
    for( int i = 0; i < n; i++) {
        const Point2D  &pa = polyPoints[i];
        const Point2D  &pb = polyPoints[(i+1)%n];
        int stat = JEdgeGeometry::intersectPredicate2d( &pa[0], &pb[0], &queryPoint[0], &extremePoint[0]);
        if( stat > 0) nCount++;
    }
    return nCount%2 == 1;
}
////////////////////////////////////////////////////////////////////////////////

double JGeometry :: getSignedArea(const double *x, const double *y, int n)
{
    // Given "n" points ( 2D ) on a polygon, (P(n+1) = p(0)) this function
    // return the signed area..
    assert( n >= 3);

    double sum = 0.0;
    for( int i  = 0; i < n; i++) {
        double x0 = x[i];
        double y0 = y[i];
        double x1 = x[(i+1)%n];
        double y1 = y[(i+1)%n];
        sum += x0*y1 - x1*y0;
    }
    return 0.5*sum;
}
////////////////////////////////////////////////////////////////////////////////

double JGeometry ::  getSignedArea( const vector<Point2D> &poly)
{
    int np = poly.size();
    vector<double> x(np), y(np);
    for( int i = 0; i < np; i++) {
        x[i] = poly[i][0];
        y[i] = poly[i][1];
    }

    return getSignedArea( &x[0], &y[0], np);
}
////////////////////////////////////////////////////////////////////////////////

void JGeometry :: getCentroid( const double *x, const double *y, int n, double *center )
{
    assert( n >= 3);

    // Formula from wikipedia. Note that centroid does not depend on the
    // orientation of the polygon. The division by Area will take care
    // of correct value.

    // For convex bodies, centroid is always inside the region.

    double cx = 0.0;
    double cy = 0.0;
    double cf;

    for( int i  = 0; i < n; i++) {
        cf  = x[i]*y[(i+1)%n] - x[(i+1)%n]*y[i];
        cx +=  (x[i]+x[(i+1)%n])*cf;
        cy +=  (y[i]+y[(i+1)%n])*cf;
    }

    double A = getSignedArea(x, y, n);

    center[0] = cx/(6.0*A);
    center[1] = cy/(6.0*A);
}

