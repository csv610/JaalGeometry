/*
 * FastArea.c++
 *
 * From the paper:
 *
 *      Daniel Sunday
 *      "Fast Polygon Area and Newell Normal Computation"
 *      journal of graphics tools, 7(2):9-13, 2002
 *
 */

// assume vertex coordinates are in arrays x[], y[], and z[]

// with room to duplicate the first two vertices at the end


// return the signed area of a 2D polygon

#include <math.h>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

inline double
PolygonArea(int n, double *x, double *y)             // 2D polygon

{
    // guarantee the first two vertices are also at array end

    x[n] = x[0];
    y[n] = y[0];
    x[n+1] = x[1];
    y[n+1] = y[1];

    double sum = 0.0;
    double *xptr = x+1, *ylow = y, *yhigh = y+2;
    for (int i=1; i <= n; i++) {
        sum += (*xptr++) * ( (*yhigh++) - (*ylow++) );
    }

    return 0.5*sum ;
}

///////////////////////////////////////////////////////////////////////////////

// return the signed area of a 3D planar polygon (given normal vector)

inline double
PolygonArea3D(int n, double *x, double *y, double *z,  // 3D planar polygon

              double nx, double ny, double nz)         // and plane normal

{
    // select largest normal coordinate to ignore for projection
    // CSV: Original Code seeems to have bug: Need to verify it.

    double ax = (nx>0 ? nx : -nx);	// abs nx

    double ay = (ny>0 ? ny : -ny);	// abs ny

    double az = (nz>0 ? nz : -nz);	// abs nz

    double len = sqrt(nx*nx + ny*ny + nz*nz); // length of normal


    if (ax > ay) {
        if (ax > az)			       // ignore x-coord

// CSV   return PolygonArea(n, y, z) * (len / nx);
            return PolygonArea(n, y, z) * (len / ax);
    }
    else if (ay > az)			       // ignore y-coord

// CSV   return PolygonArea(n, z, x) * (len / ny);
        return PolygonArea(n, z, x) * (len / ay);

//CSV   return PolygonArea(n, x, y) * (len / nz); // ignore z-coord
    return PolygonArea(n, x, y) * (len / az); // ignore z-coord

}

///////////////////////////////////////////////////////////////////////////////

// output the approximate unit normal of a 3D nearly planar polygon
// return the area of the polygon

inline double
PolygonNormal3D(int n, double *x, double *y, double *z,
                double *nx, double *ny, double *nz)

{
    // get the Newell normal
    double nwx = PolygonArea(n, y, z);
    double nwy = PolygonArea(n, z, x);
    double nwz = PolygonArea(n, x, y);

    // get length of the Newell normal

    double nlen = sqrt( nwx*nwx + nwy*nwy + nwz*nwz );
    // compute the unit normal

    double multby = 1.0/( double)nlen;

    *nx = multby*nwx;
    *ny = multby*nwy;
    *nz = multby*nwz;

    return nlen;    // area of polygon = length of Newell normal
}

///////////////////////////////////////////////////////////////////////////////

inline
double PolygonArea3D( int n, double *x, double *y, double *z)
{
    double nx, ny, nz;
    PolygonNormal3D( n, x, y, z, &nx, &ny, &nz );

    return PolygonArea3D( n, x, y, z, nx, ny, nz );
}

///////////////////////////////////////////////////////////////////////////////
