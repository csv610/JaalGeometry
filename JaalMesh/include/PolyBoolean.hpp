#pragma once

#ifndef JPOLYBOOL_H
#define JPOLYBOOL_H

#include <vector>
#include "Mesh.hpp"

using namespace std;
using namespace Jaal;

struct JPolyBoolean {
    typedef vector<Point2D>  PolyType;

    int getUnionOf( PolyType &a, PolyType &b, PolyType &c);
    int getDifferencef( PolyType &a, PolyType &b, PolyType &c);
    int getIntersectionOf(PolyType &a, PolyType &b, PolyType &c);

    int getOrientation(const PolyType &p);
    int intersect(const PolyType &a, const PolyType &b);

    int unitTests();

private:
    struct PolyPoint
    {
        Point2D xy;
        int     sign;
        int     id;
        double  u;
        bool operator < ( const PolyPoint &rhs) const {
            return u < rhs.u;
        }
    };

    vector<PolyPoint> points;
    vector< pair<int,int> >  edges;

    void extract_segments( const vector<int> &p, const vector< vector<int> > &edgepoints);
    void create_chain();
};

#endif

