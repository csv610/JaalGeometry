#pragma once

#ifndef KNOTS_H
#define KNOTS_H

#include "Mesh.hpp"
#include "Curve.hpp"

struct Knot
{
    vector<JCurve*> loops;
    vector<int>     getGaussCode();
    vector<int>     getDTCode();
    vector<int>     getConway();
};

class ClassicalKnots
{
public:
    Knot*  getKnot( int component, int crossing, int index );
    Knot*  getKnot( int id );
private:
    void readKnotFile();
};

#endif

