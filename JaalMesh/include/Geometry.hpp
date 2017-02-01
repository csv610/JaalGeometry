#pragma once

#ifndef BASICGEOM_H
#define BASICGEOM_H

#include "basic_math.hpp"

struct JGeometry
{
    static int isInside( const vector<Point2D> &polyPoints,  Point2D &queryPoint);
    static double getSignedArea(const double *x, const double *y, int n);
    static double getSignedArea( const vector<Point2D> &poly);
    static void getCentroid( const double *x, const double *y, int n, double *center );
    static int  getBoundedSide( const double *p0, const double *p1, const double *p2, const double *ptest);
};

#endif





