#ifndef CCENTER_H
#define CCENTER_H

#include "GeomPredicates.hpp"
#include "basic_math.hpp"

#ifdef EXACT
#define EXACT 1
#endif

void TriCircumCenter2D( const double *a, const double *b, const double *c, double *center);
void TriCircumCenter2D( const double *a, const double *b, const double *c, double *center, double *param);

void TriCircumCenter3D( const double *a, const double *b, const double *c, double *center);
void TriCircumCenter3D( const double *a, const double *b, const double *c, double *center, double *param);

void TetCircumCenter( const double *a, const double *b, const double *c, const double *d, double *r, double *p);

namespace UnitTest
{
int test_tet_circumcenter();
int test_tri2d_circumcenter();
int test_tri3d_circumcenter();
}

#endif
