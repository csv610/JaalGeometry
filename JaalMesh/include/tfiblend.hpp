#ifndef TFIBLEND_H
#define TFIBLEND_H

namespace TFI
{
double linear_interpolation(double r, double x0, double x1);
double bilinear_interpolation(double r, double s, double *valCorners);
double trilinear_interpolation(double r, double s, double t, double *valCorners);

void blend_from_corners ( double *x, int m );
void blend_from_corners ( double *x, int nx, int ny );
void blend_from_corners ( double *x, int nx, int ny, int nz );

void blend_from_edges   ( double *x, int nx, int ny );
void blend_from_edges   ( double *x, int nx, int ny, int nz );

void blend_from_faces   ( double *x, int nx, int ny, int nz );

double transfinite_blend(double r, double s,
                         double x00, double x10, double x11, double x01,
                         double xr0, double x1s, double xr1, double x0s);
}

#endif
