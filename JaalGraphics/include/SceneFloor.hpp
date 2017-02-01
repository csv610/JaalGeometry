#pragma once

#include "Texture.hpp"
#include "EntityColor.hpp"

class JSceneFloor
{
public:
    static const int  GRID_LINES   = 0;
    static const int  CHECKERBOARD = 1;
    static const int  CUSTOM_IMAGE = 2;

    JSceneFloor();

    void setPattern(int p );

    void setStdPlane( int d )
    {
        plane = d;
    }

    void setStatus( bool v )
    {
        status = v;
    }

    void setDistance( double d )
    {
        floorDistance = d;
    }

    void setLength( double d )
    {
        floorLength = d;
        genTexture();
    }

    void setNumLines( int  d )
    {
        numFloorLines = d;
        genTexture();
    }

    void setColor( const JColor &c) {
        color = c;
    }

    void draw();


private:
    int    status   = 0;
    int    plane    = 1;
    int    pattern  = 1;
    int    numFloorLines = 10;
    double floorDistance = 0.0;
    double floorLength   = 5.0;
    JColor color;
    boost::scoped_ptr<JTexture> texturePtr;

    void genTexture();
    void texturedFloor();
    void drawLines();
};
///////////////////////////////////////////////////////////////////////////////
typedef boost::shared_ptr<JSceneFloor>  JSceneFloorPtr;
