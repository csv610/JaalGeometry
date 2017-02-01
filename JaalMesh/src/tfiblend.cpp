//****************************************************************************//
// Acknowlegement:  The original Code "Blend" from John Burkardt ( Under GNU
//                  lincense) have been significantly modified. Anyhow full
//                  credit is given to original developers.
//
//                  Document is also copied from the original code.
//
// Modification Date :  7th July 2009;
// Modified By       :  Chaman Singh Verma
//                      Argonne National Lab, Argonne, IL, USA.
// Major Changes:
//                 (1) Parameters range (-1,1) instead of (0,1).
//                 (2) Probably better function names
//                 (3) Better modularity and reuse
//                 (4) Changed ordering (k,j,i) instead of (i,j,k)
//                 (5) Guass point introduction.
//                 (6) Probably much easier to understand than the original code.
//                 (7) Probably better input parameter passing order.
//                 (8) Error checking using assertions introduced.
// Downside:
//                 (1) May be little more expensive ( than the original code) because of
//                     function calling overheads.
//
//                 (2) Redundant calculations to make code easier to understand.
//
////////////////////////////////////////////////////////////////////////////////

#include <assert.h>
#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

namespace TFI {

double gauss_node(int i, int N)
{
    double du = 2.0 / (double) (N - 1);
    return -1.0 + i*du;
}

//****************************************************************************80

double linear_interpolation(double r, double x0, double x1)
{

    //************************************************************************80
    //
    //  Purpose: extends scalar data at endpoints to a line.
    //
    //  Diagram:
    //
    //    -1-----r-----1
    //
    //  Author:
    //
    //    John Burkardt: Modified by Chaman Singh Verma
    //
    //  Parameters:
    //
    //    Input, double R, the coordinate where an interpolated value is desired.
    //
    //    Input, double X0, X1, the data values at the ends of the line.
    //
    //    Output, double *X, the interpolated data value at (R).
    //
    ///////////////////////////////////////////////////////////////////////////

    assert(r >= -1.0 && r <= 1.0);

    double val = (1.0 - r) * x0 + (1.0 + r) * x1;

    val *= 0.5;

    return val;
}

//****************************************************************************80

double bilinear_interpolation(double r, double s, double *valCorners)
{
    assert(r >= -1.0 && r <= 1.0);
    assert(s >= -1.0 && s <= 1.0);

    double val = (1.0 - r) * (1.0 - s) * valCorners[0] +
                 (1.0 + r) * (1.0 - s) * valCorners[1] +
                 (1.0 + r) * (1.0 + s) * valCorners[2] +
                 (1.0 - r) * (1.0 + s) * valCorners[3];

    val *= 0.25;
    return val;
}

//****************************************************************************80

double bilinear_interpolation(double r, double s,
                              double x00, double x10, double x11, double x01)
{
    double valCorners[4];

    valCorners[0] = x00;
    valCorners[1] = x10;
    valCorners[2] = x11;
    valCorners[3] = x01;

    double val = bilinear_interpolation(r, s, valCorners);

    return val;
}

//****************************************************************************80

double trilinear_interpolation(double r, double s, double t, double *valCorners)
{
    assert(r >= -1.0 && r <= 1.0);
    assert(s >= -1.0 && s <= 1.0);
    assert(t >= -1.0 && t <= 1.0);

    double val = (1.0 - r) * (1.0 - s) * (1.0 - t) * valCorners[0] +
                 (1.0 + r) * (1.0 - s) * (1.0 - t) * valCorners[1] +
                 (1.0 + r) * (1.0 + s) * (1.0 - t) * valCorners[2] +
                 (1.0 - r) * (1.0 + s) * (1.0 - t) * valCorners[3] +
                 (1.0 - r) * (1.0 - s) * (1.0 + t) * valCorners[4] +
                 (1.0 + r) * (1.0 - s) * (1.0 + t) * valCorners[5] +
                 (1.0 + r) * (1.0 + s) * (1.0 + t) * valCorners[6] +
                 (1.0 - r) * (1.0 + s) * (1.0 + t) * valCorners[7];

    val *= 0.125;

    return val;
}

//****************************************************************************80

double trilinear_interpolation(double r, double s, double t,
                               double x000, double x100, double x110, double x010,
                               double x001, double x101, double x111, double x011)
{
    double valCorners[8];

    valCorners[0] = x000;
    valCorners[1] = x100;
    valCorners[2] = x110;
    valCorners[3] = x010;

    valCorners[4] = x001;
    valCorners[5] = x101;
    valCorners[6] = x111;
    valCorners[7] = x011;

    double val = trilinear_interpolation(r, s, t, valCorners);

    return val;
}

//****************************************************************************80

double transfinite_blend(double r, double s,
                         double x00, double x10, double x11, double x01,
                         double xr0, double x1s, double xr1, double x0s)
{

    //****************************************************************************80
    //
    //  Purpose:
    //
    //    BLEND_1D1 extends scalar data along the boundary into a square.
    //
    //  Diagram:
    //
    //    01-----r1-----11
    //     |      .      |
    //     |      .      |
    //    0s.....rs.....1s
    //     |      .      |
    //     |      .      |
    //    00-----r0-----10
    //
    //  Formula:
    //
    //    Written as a polynomial in R and S, the interpolation map has the form
    //
    //      X(R,S) =
    //           1     * ( x0s + xr0 - x00 )
    //         + r     * ( x00 + x1s - x0s - x10 )
    //         + s     * ( x00 + xr0 - x01 - xr1 )
    //         + r * s * ( x01 + x10 - x00 - x11 )
    //
    //    The nonlinear term ( r * s ) has an important role:
    //
    //    If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in a plane,
    //    and the mapping is affine.  All the interpolated data will lie
    //    on the plane defined by the four corner values.  In particular,
    //    on any line through the square, data values at intermediate points
    //    will lie between the values at the endpoints.
    //
    //    If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
    //    not lie in a plane, and the interpolation map is nonlinear.  On
    //    any line through the square, data values at intermediate points
    //    may lie above or below the data values at the endpoints.  The
    //    size of the coefficient of r * s will determine how severe this
    //    effect is.
    //
    //  Author:
    //
    //    John Burkardt: Modified by Chaman Singh Verma
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input, double R, S, the coordinates where an interpolated value is desired.
    //
    //    Input, double X00, X01, X10, X11, the data values at the corners.
    //
    //    Input, double XR0, XR1, X0S, X1S, the data values at points along the edges.
    //
    //    Output, double *X, the interpolated data value at (R,S).
    //
    ////////////////////////////////////////////////////////////////////////////

    double u = linear_interpolation(r, x0s, x1s);
    double v = linear_interpolation(s, xr0, xr1);
    double uv = bilinear_interpolation(r, s, x00, x10, x11, x01);

    double val = u + v - uv;

    return val;
}
//****************************************************************************80

double transfinite_blend(double r, double s, double t,
                         double x000, double xr00, double x100,
                         double x0s0, double xrs0, double x1s0,
                         double x010, double xr10, double x110,
                         double x00t, double xr0t, double x10t,
                         double x0st, double x1st,
                         double x01t, double xr1t, double x11t,
                         double x001, double xr01, double x101,
                         double x0s1, double xrs1, double x1s1,
                         double x011, double xr11, double x111)
{

    //************************************************************************80
    //
    //  Purpose:  extends scalar data along the surface into a cube.
    //
    //  Diagram:
    //
    //    010-----r10-----110    011-----r11-----111
    //      |       .       |      |       .       |
    //      |       .       |      |       .       |
    //    0s0.....rs0.....1s0    0s1.....rs1.....1s1     S
    //      |       .       |      |       .       |     |
    //      |       .       |      |       .       |     |
    //    000-----r00-----100    001-----r01-----101     +----R
    //       BACK                     FRONT
    //
    //    001-----0s1-----011    101-----1s1-----111
    //      |       .       |      |       .       |
    //      |       .       |      |       .       |
    //    00t.....0st.....01t    10t.....1st.....11t     T
    //      |       .       |      |       .       |     |
    //      |       .       |      |       .       |     |
    //    000-----0s0-----010    100-----1s0-----110     +----S
    //       LEFT                       RIGHT

    //    001-----r01-----101    011-----r11-----111
    //      |       .       |      |       .       |
    //      |       .       |      |       .       |
    //    00t.....r0t.....10t    01t.....r1t.....11t     T
    //      |       .       |      |       .       |     |
    //      |       .       |      |       .       |     |
    //    000-----r0t-----100    010-----r10-----110     +----R
    //       BOTTOM                       TOP
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input, double R, S, T, the coordinates where an interpolated value is desired.
    //
    //    Input, double X000, X001, X010, X011, X100, X101, X110, X111, the data
    //    values at the corners.
    //
    //    Input, double XR00, XR01, XR10, XR11, X0S0, X0S1, X1S0, X1S1, X00T, X01T,
    //    X10T, X11T, the data values at points along the edges.
    //
    //    Input, double X0ST, X1ST, XR0T, XR1T, XRS0, XRS1, the data values
    //    at points on the faces.
    //
    //    Output, double *X, the interpolated data value at (R,S,T).
    //
    ////////////////////////////////////////////////////////////////////////////

    double u = linear_interpolation(r, x0st, x1st);
    double v = linear_interpolation(s, xr0t, xr1t);
    double w = linear_interpolation(t, xrs0, xrs1);

    double uv = bilinear_interpolation(r, s, x00t, x10t, x11t, x01t);
    double uw = bilinear_interpolation(r, t, x0s0, x1s0, x1s1, x0s1);
    double vw = bilinear_interpolation(s, t, xr00, xr10, xr11, xr01);

    double uvw = trilinear_interpolation(r, s, t,
                                         x000, x100, x110, x010,
                                         x001, x101, x111, x011);

    double val = u + v + w - uw - uv - vw + uvw;

    return val;
}

//****************************************************************************80

void blend_from_corners(double *x, int n)
{

//****************************************************************************80
//
//  Purpose: extends indexed scalar data at endpoints along a line.
//
//  Diagram:
//
//    ( X0, ..., ..., ..., ..., ..., X6 )
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gordon, Charles A Hall,
//    Construction of Curvilinear Coordinate Systems and Application to
//    Mesh Generation,
//    International Journal of Numerical Methods in Engineering,
//    Volume 7, pages 461-477, 1973.
//
//    Joe Thompson, Bharat Soni, Nigel Weatherill,
//    Handbook of Grid Generation,
//    CRC Press,
//    1999.
//
//  Parameters:
//
//    Input/output, double X[M].
//
//    On input, X[0] and X[M-1] contain scalar values which are to be
//    interpolated through the entries X[1] through X[M-2].  It is assumed
//    that the dependence of the data is linear in the vector index I.
//
//    On output, X[1] through X[M-2] have been assigned interpolated values.
//
//    Input, int M, the number of entries in X.
//
    int i;
    double r;

    for (i = 1; i < n - 1; i++) {
        r = gauss_node(i, n);
        x[i] = linear_interpolation(r, x[0], x[n - 1]);
    }

    return;
}
//****************************************************************************80

void blend_from_edges(double *x, int nx, int ny)
{

    //****************************************************************************80
    //
    //  Purpose: extends indexed scalar data along edges into a table.
    //
    //  Diagram:
    //
    //    ( X00,  X01,  X02,  X03,  X04,  X05,  X06 )
    //    ( X10,  ...,  ...,  ...,  ...,  ...,  X16 )
    //    ( X20,  ...,  ...,  ...,  ...,  ...,  X26 )
    //    ( X30,  ...,  ...,  ...,  ...,  ...,  X36 )
    //    ( X40,  X41,  X42,  X43,  X44,  X45,  X46 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    7th July December 2009
    //
    //  Author:
    //
    //    John Burkardt: Modified by Chaman Singh Verma
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY], a singly dimensioned array that
    //    is "really" doubly dimensioned.  The double dimension index [I][J]
    //    corresponds to the single dimension index J * NX + I.
    //
    //    On input, data is contained in the "edge entries"
    //    X[0][0], X[I][0], X[0][NY-1] and X[NX-1][NY-1],
    //    for I = 0 to NX-1, and J = 0 to NY-1.
    //
    //    On output, all entries in X have been assigned a value.
    //
    //    Input, int NX, NY  the number of rows and columns in X.
    //
    ///////////////////////////////////////////////////////////////////////////////

    double x0 = x[0];
    double x1 = x[nx - 1];
    double x2 = x[nx * ny - 1];
    double x3 = x[nx * ny - nx];

    for (int j = 1; j < ny - 1; j++) {
        double s = gauss_node(j, ny);

        for (int i = 1; i < nx - 1; i++) {
            double r = gauss_node(i, nx);
            int il = j*nx;
            int ir = j * nx + nx - 1;
            int it = (ny - 1) * nx + i;
            int ib = i;
            int offset = j * nx + i;

            double xr0 = x[ib];
            double x1s = x[ir];
            double xr1 = x[it];
            double x0s = x[il];

            x[offset] = transfinite_blend(r, s, x0, x1, x2, x3, xr0, x1s, xr1, x0s);
        }
    }
    return;
}

////////////////////////////////////////////////////////////////////////////////

void blend_from_corners(double *x, int nx, int ny)
{

    //****************************************************************************80
    //
    //  Purpose: extends indexed scalar data at corners into a table.
    //
    //  Diagram:
    //
    //    ( X00,  ..., ..., ..., ..., ..., X06 )
    //    ( ...,  ..., ..., ..., ..., ..., ... )
    //    ( ...,  ..., ..., ..., ..., ..., ... )
    //    ( ...,  ..., ..., ..., ..., ..., ... )
    //    ( X40,  ..., ..., ..., ..., ..., X46 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    19 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY], a singly dimensioned array that
    //    is "really" doubly dimensioned.  The double dimension index [I][J]
    //    corresponds to the single dimension index J * NX + I.
    //
    //    On input, data values have been stored in the entries
    //    [0], [NX-1], [NX*NY-NX] and [NX*NY-1], which correspond to the double
    //    dimension entries [0][0], [1][NX-1], [0][NY-1] and [NX-1][NY-1].
    //
    //    On output, all entries in X have been assigned a value.
    //
    //    Input, int NX, NY, the number of rows and columns in the doubly
    //    dimensioned data.
    //
    ///////////////////////////////////////////////////////////////////////////////

    int offset;
    double r, s, x0, x1;

    int nxy = nx*ny;

    for (int i = 1; i < nx - 1; i++) {
        r = gauss_node(i, nx);

        x0 = x[0];
        x1 = x[nx - 1];
        offset = i;
        x[offset] = linear_interpolation(r, x0, x1);


        x0 = x[nxy - nx];
        x1 = x[nxy - 1 ];
        offset = (ny - 1) * nx + i;
        x[offset] = linear_interpolation(r, x0, x1);
    }


    for (int j = 1; j < ny - 1; j++) {
        s = gauss_node(j, ny);

        x0 = x[0];
        x1 = x[nx * ny - nx];
        offset = j*nx;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[nx - 1];
        x1 = x[nxy - 1];
        offset = j * nx + (nx - 1);
        x[offset] = linear_interpolation(s, x0, x1);
    }

    blend_from_edges(x, nx, ny);

    return;
}


//****************************************************************************80

void blend_from_faces(double *x, int nx, int ny, int nz)
{

    //************************************************************************80
    //
    //  Purpose:
    //
    //    BLEND_FROM_FACES extends indexed scalar data along faces into a cubic table.
    //
    //  Diagram:
    //
    //    ( X000    X010    X020    X030    X040    X050 )
    //    ( X100    X110    X120    X130    X140    X150 )
    //    ( X200    X210    X220    X230    X240    X250 )   Layer 1
    //    ( X300    X310    X320    X330    X340    X350 )
    //    ( X400    X410    X420    X430    X440    X450 )
    //
    //    ( X001    X011    X021    X031    X041    X051 )
    //    ( X101    ...     ....    ....    ....    X151 )
    //    ( X201    ...     ....    ....    ....    X251 )   Layer K
    //    ( X301    ...     ....    ....    ....    X351 )   1 < K < M3
    //    ( X401    X411    X421    X431    X441    X451 )
    //
    //    ( X002    X012    X022    X032    X042    X052 )
    //    ( X102    X112    X122    X132    X142    X152 )
    //    ( X202    X212    X222    X232    X242    X252 )   Layer M3
    //    ( X302    X312    X322    X332    X342    X352 )
    //    ( X402    X412    X422    X432    X442    X452 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    22 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY*NZ], a singly dimensioned array that
    //    is "really" triply dimensioned.  The triple dimension index
    //    [I][J][K] corresponds to the single dimension index
    //    K * NX*NY + J * NX + I
    //
    //    On input, there is already scalar data in the entries X[I][J][K]
    //    corresponding to "faces" of the table, that is, entries for which
    //    at least one of the three indices I, J and K is equal to their
    //    minimum or maximum possible values.
    //
    //    On output, all entries in X have been assigned a value, using the
    //    table indices as independent variables.
    //
    //    Input, int NX, NY, NY, the number of rows, columns, and layers in X.
    //
    ////////////////////////////////////////////////////////////////////////////

    double r, s, t;

    int offset, nxy = nx*ny;

    for (int k = 1; k < nz - 1; k++) {
        t = gauss_node(k, nz);
        for (int j = 1; j < ny - 1; j++) {
            s = gauss_node(j, ny);
            for (int i = 1; i < nx - 1; i++) {
                r = gauss_node(i, nx);

                // Points on Back Plane
                double x000 = x[0];
                double xr00 = x[i];
                double x100 = x[(nx - 1)];

                double x0s0 = x[j * nx];
                double xrs0 = x[i + j * nx];
                double x1s0 = x[(nx - 1) + j * nx];

                double x010 = x[(ny - 1) * nx];
                double xr10 = x[i + (ny - 1) * nx];
                double x110 = x[(ny - 1) * nx + (nx - 1)];

                // Intermediate Plane

                double x00t = x[k * nxy];
                double xr0t = x[i + k * nxy];
                double x10t = x[(nx - 1) + k * nxy];

                double x0st = x[j * nx + k * nxy];
                double x1st = x[(nx - 1) + j * nx + k * nxy];

                double x01t = x[ (ny - 1) * nx + k * nxy];
                double xr1t = x[i + (ny - 1) * nx + k * nxy];
                double x11t = x[(nx - 1) + (ny - 1) * nx + k * nxy];

                // Front Plane
                double x001 = x[(nz - 1) * nxy];
                double xr01 = x[ i + (nz - 1) * nxy];
                double x101 = x[(nz - 1) * nxy + (nx - 1)];

                double x0s1 = x[ j * nx + (nz - 1) * nxy];
                double xrs1 = x[ i + j * nx + (nz - 1) * nxy];
                double x1s1 = x[ (nx - 1) + j * nx + (nz - 1) * nxy];

                double x011 = x[(ny - 1) * nx + (nz - 1) * nxy];
                double xr11 = x[ i + (ny - 1) * nx + (nz - 1) * nxy];
                double x111 = x[(nz - 1) * nxy + (ny - 1) * nx + (nx - 1)];

                offset = k * nxy + j * nx + i;
                x[offset] = transfinite_blend(r, s, t,
                                              x000, xr00, x100,
                                              x0s0, xrs0, x1s0,
                                              x010, xr10, x110,
                                              x00t, xr0t, x10t,
                                              x0st, x1st,
                                              x01t, xr1t, x11t,
                                              x001, xr01, x101,
                                              x0s1, xrs1, x1s1,
                                              x011, xr11, x111);
            }

        }

    }

    return;
}

////////////////////////////////////////////////////////////////////////////////

void blend_from_edges(double *x, int nx, int ny, int nz)
{

    //************************************************************************80
    //
    //  Purpose: extends indexed scalar data along "edges" into a cubic table.
    //
    //  Diagram:
    //
    //    ( X000,   X010,   X020,   X030,   X040,   X050 )
    //    ( X100,   ...,    ...,    ...,    ...,    X150 )
    //    ( X200,   ...,    ...,    ...,    ...,    X250 )   Layer 1
    //    ( X300,   ...,    ...,    ...,    ...,    X350 )
    //    ( X400,   X410,   X420,   X430,   X440,   X450 )
    //
    //    ( X001,   ...,    ...,    ...,    ...,    X051 )
    //    ( ....,   ...,    ...,    ...,    ...,    ...  )
    //    ( ....,   ...,    ...,    ...,    ...,    ...  )   Layer K
    //    ( ....,   ...,    ...,    ...,    ...,    ...  )   1 < K < M3
    //    ( X401,   ...,    ...,    ...,    ...,    X451 )
    //
    //    ( X002,   X012,   X022,   X032,   X042,   X052 )
    //    ( X102,   ...,    ...,    ...,    ...,    X152 )
    //    ( X202,   ...,    ...,    ...,    ...,    X252 )   Layer M3
    //    ( X302    ...,    ...,    ...,    ...,    X352 )
    //    ( X402,   X412,   X422,   X432,   X442,   X452 )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    22 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY*NZ], a singly dimensioned array that
    //    is "really" triply dimensioned.  The triple dimension index
    //    [I][J][K] corresponds to the single dimension index
    //    K * NX*NY + J * NX + I
    //
    //    On input, there is already scalar data in the entries X[I][J][K]
    //    corresponding to "edges" of the table, that is, entries for which
    //    at least two of the three indices I, J and K are equal to their
    //    minimum or maximum possible values.
    //
    //    Input, int NX, NY, NZ, the number of rows, columns, and layers in X.
    //
    ////////////////////////////////////////////////////////////////////////////

    int offset;
    int nxy = nx*ny;

    double r, s, t;
    double x00, x10, x11, x01;
    double xs0, x1t, xs1, x0t;

    for (int k = 1; k < nz - 1; k++) {
        t = gauss_node(k, nz);
        for (int j = 1; j < ny - 1; j++) {
            s = gauss_node(j, ny);
            //Left Face ...
            x00 = x[0];
            x10 = x[(ny - 1) * nx];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx];
            x01 = x[(nz - 1) * nxy];

            xs0 = x[j * nx];
            x1t = x[k * nxy + (ny - 1) * nx];
            xs1 = x[(nz - 1) * nxy + j * nx];
            x0t = x[k * nxy];

            offset = k * nxy + j*nx;

            x[offset] = transfinite_blend(s, t, x00, x10, x11, x01, xs0, x1t, xs1, x0t);

            // Right Face ...
            x00 = x[nx - 1];
            x10 = x[(ny - 1) * nx + (nx - 1)];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx + (nx - 1)];
            x01 = x[(nz - 1) * nxy + (nx - 1)];

            xs0 = x[j * nx + (nx - 1)];
            x1t = x[k * nxy + (ny - 1) * nx + (nx - 1)];
            xs1 = x[(nz - 1) * nxy + j * nx + (nx - 1)];
            x0t = x[k * nxy + (nx - 1)];

            offset = k * nxy + j * nx + nx - 1;

            x[offset] = transfinite_blend(s, t, x00, x10, x11, x01, xs0, x1t, xs1, x0t);
        }
    }

    double xr0, xr1;

    for (int k = 1; k < nz - 1; k++) {
        t = gauss_node(k, nz);
        for (int i = 1; i < nx - 1; i++) {
            r = gauss_node(i, nx);

            // Bottom Face ...
            x00 = x[0];
            x10 = x[nx - 1];
            x11 = x[(nz - 1) * nxy + (nx - 1)];
            x01 = x[(nz - 1) * nxy];

            xr0 = x[i];
            x1t = x[k * nxy + (nx - 1)];

            xr1 = x[(nz - 1) * nxy + i];
            x0t = x[k * nxy];

            offset = k * nxy + i;

            x[offset] = transfinite_blend(r, t, x00, x10, x11, x01, xr0, x1t, xr1, x0t);

            // Top Face ...
            x00 = x[(ny - 1) * nx];
            x10 = x[(ny - 1) * nx + nx - 1];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx + nx - 1];
            x01 = x[(nz - 1) * nxy + (ny - 1) * nx];

            xr0 = x[(ny - 1) * nx + i];
            x1t = x[k * nxy + (ny - 1) * nx + nx - 1];
            xr1 = x[(nz - 1) * nxy + (ny - 1) * nx + i];
            x0t = x[k * nxy + (ny - 1) * nx];

            offset = k * nxy + (ny - 1) * nx + i;

            x[offset] = transfinite_blend(r, t, x00, x10, x11, x01, xr0, x1t, xr1, x0t);

        }
    }

    double x0s, x1s;
    for (int j = 1; j < ny - 1; j++) {
        s = gauss_node(j, ny);
        for (int i = 1; i < nx - 1; i++) {
            r = gauss_node(i, nx);

            // Back Face ...
            x00 = x[0];
            x10 = x[nx - 1];
            x11 = x[(ny - 1) * nx + (nx - 1)];
            x01 = x[(ny - 1) * nx];

            xr0 = x[i];
            x1s = x[j * nx + nx - 1];
            xr1 = x[(ny - 1) * nx + i];
            x0s = x[j * nx];

            offset = j * nx + i;

            x[offset] = transfinite_blend(r, s, x00, x10, x11, x01, xr0, x1s, xr1, x0s);

            // Front Face ...
            x00 = x[(nz - 1) * nxy];
            x10 = x[(nz - 1) * nxy + nx - 1];
            x11 = x[(nz - 1) * nxy + (ny - 1) * nx + (nx - 1)];
            x01 = x[(nz - 1) * nxy + (ny - 1) * nx];

            xr0 = x[(nz - 1) * nxy + i];
            x1s = x[(nz - 1) * nxy + j * nx + nx - 1];
            xr1 = x[(nz - 1) * nxy + (ny - 1) * nx + i];
            x0s = x[(nz - 1) * nxy + j * nx];

            offset = (nz - 1) * nx * ny + j * nx + i;

            x[offset] = transfinite_blend(r, s, x00, x10, x11, x01, xr0, x1s, xr1, x0s);
        }
    }

    blend_from_faces(x, nx, ny, nz);

    return;
}

////////////////////////////////////////////////////////////////////////////////

void blend_from_corners(double *x, int nx, int ny, int nz)
{

    //************************************************************************80
    //
    //  Purpose:
    //
    //    BLEND_IJK_0D1 extends indexed scalar data along corners into a cubic table.
    //
    //  Diagram:
    //
    //    ( X000,   ...,  ...,  ...,  ...,  ...,  X060 )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   First "layer"
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( X400,   ...,  ...,  ...,  ...,  ...,  X460 )
    //
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   Middle "layers"
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //
    //    ( X003,  ...,  ...,  ...,  ...,  ...,  X063  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )   Last "layer"
    //    ( ....,   ...,  ...,  ...,  ...,  ...,  ...  )
    //    ( X403,  ...,  ...,  ...,  ...,  ...,  X463  )
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    22 December 1998
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    William Gordon, Charles A Hall,
    //    Construction of Curvilinear Coordinate Systems and Application to
    //    Mesh Generation,
    //    International Journal of Numerical Methods in Engineering,
    //    Volume 7, pages 461-477, 1973.
    //
    //    Joe Thompson, Bharat Soni, Nigel Weatherill,
    //    Handbook of Grid Generation,
    //    CRC Press,
    //    1999.
    //
    //  Parameters:
    //
    //    Input/output, double X[NX*NY*NZ], a singly dimensioned array that
    //    is "really" triply dimensioned.  The triple dimension index
    //    [I][J][K] corresponds to the single dimension index
    //    K * NX*NY + J * NX + I
    //
    //    On input, there is already scalar data in the entries X[I][J][K]
    //    corresponding to "cornders" of the table, that is, entries for which
    //    each of the three indices I, J and K is equal to their
    //    minimum or maximum possible values.
    //
    //    Input, int NX, NY, NZ, the number of rows, columns, and layers in X.
    //
    ///////////////////////////////////////////////////////////////////////////

    int offset;
    double r, s, t, x0, x1;

    int nxy = nx*ny;

    for (int i = 1; i < nx - 1; i++) {
        r = gauss_node(i, nx);

        x0 = x[0];
        x1 = x[nx - 1];
        offset = i;
        x[offset] = linear_interpolation(r, x0, x1);

        x0 = x[(ny - 1) * nx];
        x1 = x[(ny - 1) * nx + nx - 1];
        offset = (ny - 1) * nx + + i;
        x[offset] = linear_interpolation(r, x0, x1);

        x0 = x[(nz - 1) * nxy];
        x1 = x[(nz - 1) * nxy + nx - 1];
        offset = (nz - 1) * nxy + i;
        x[offset] = linear_interpolation(r, x0, x1);

        x0 = x[(nz - 1) * nxy + (ny - 1) * nx ];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx + nx - 1];
        offset = (nz - 1) * nxy + (ny - 1) * nx + i;
        x[offset] = linear_interpolation(r, x0, x1);

    }

    for (int j = 1; j < ny - 1; j++) {
        s = gauss_node(j, ny);

        x0 = x[0];
        x1 = x[ (ny - 1) * nx ];
        offset = j*nx;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[nx - 1];
        x1 = x[(ny - 1) * nx + (nx - 1)];
        offset = j * nx + nx - 1;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[(nz - 1) * nxy];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx];
        offset = (nz - 1) * nxy + j*nx;
        x[offset] = linear_interpolation(s, x0, x1);

        x0 = x[ (nz - 1) * nxy + nx - 1];
        x1 = x[ (nz - 1) * nxy + (ny - 1) * nx + nx - 1];
        offset = (nz - 1) * nxy + j * nx + nx - 1;
        x[offset] = linear_interpolation(s, x0, x1);
    }

    for (int k = 1; k < nz - 1; k++) {
        t = gauss_node(k, nz);

        x0 = x[0];
        x1 = x[(nz - 1) * nxy];
        offset = k*nxy;
        x[offset] = linear_interpolation(t, x0, x1);

        x0 = x[nx - 1];
        x1 = x[(nz - 1) * nxy + nx - 1];
        offset = k * nxy + nx - 1;
        x[offset] = linear_interpolation(t, x0, x1);

        x0 = x[(ny - 1) * nx];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx ];
        offset = k * nxy + (ny - 1) * nx;
        x[offset] = linear_interpolation(t, x0, x1);

        x0 = x[(ny - 1) * nx + nx - 1];
        x1 = x[(nz - 1) * nxy + (ny - 1) * nx + nx - 1];
        offset = k * nx * ny + (ny - 1) * nx + nx - 1;
        x[offset] = linear_interpolation(t, x0, x1);
    }

    blend_from_edges(x, nx, ny, nz);

    return;
}

}

///////////////////////////////////////////////////////////////////////////////

