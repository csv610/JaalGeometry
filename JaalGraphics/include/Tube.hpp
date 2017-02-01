#pragma once

#include "MeshEntity.hpp"
#include "EntityColor.hpp"

#include <GL/glu.h>
#include <GL/gle.h>

class Tube {
    typedef gleDouble  RowArray[3];

public:
    Tube() {
        numSides = 20;
        numEdges = 0;
        numPoints = 0;
        points = NULL;
        radius = 0.01;
        color[0] = 1.0;
        color[1] = 0.5;
        color[2] = 0.2;
        color[3] = 1.0;
    }

    void clear() {
        if( points ) delete [] points;
        points = NULL;
        numEdges = 0;
    }

    void setRadius( double r) {
        radius = r;
    }

    void setColor( const JColor &c)  {
        color  = c;
    }

    void setNumSides(int n) {
        numSides = n;
    }

    void setSegments( const JEdgeSequence &eseq);

    void draw();
private:
    size_t numEdges, numPoints;
    int    numSides;
    double radius;
    JColor  color;
    RowArray *points;
};
