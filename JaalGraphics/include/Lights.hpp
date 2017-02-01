#pragma once

#include <GL/gl.h>

#include "EntityColor.hpp"
#include "JaalViewer.hpp"

using namespace Jaal;

class JaalViewer;

class JLights
{
public:
    static const int POSITIONAL  = 0;
    static const int DIRECTIONAL = 1;
    static const int SPOT        = 2;

    JLights() {
        init();
    }

    void setViewManager( JaalViewer *v) {
        viewManager = v;
    }

    int getNumLights() const;

    bool isDirectional( int lnum ) {
        if(lightPos[lnum][3] == 0)
            return 1;
    }

    void Switch( bool v );

    void Switch( int id, bool v) {
        status[id] = v;
    }

    /*
        void setLightPlacementPolicy( int p )
        {
            lights_placement_policy = p;
            init();
        }
    */

    void setLocalViewer( bool v );
    void setTwoSided( bool v );
    void setAmbientModel( const JColor &color);

    void setPosition( int light,  const Point4F &pos);
    void setAmbientColor( int light, const JColor &color);
    void setDiffuseColor( int light, const JColor &color);
    void setSpecularColor( int light, const JColor &color);

    void setSpotDirection( int light,  const Point3F &dir);
    void setSpotExponent( int  light,  double val);
    void setSpotCutoff(  int light,  double angle);

    Point4F getPosition( int lnum );

    void setConstantAttenuation( int light,  double val);
    void setLinearAttenuation( int light,  double angle);
    void setQuadraticAttenuation( int light,  double val);

    void setLightBulbs( bool v ) {
        light_bulbs = v;
    }

    void setBulbRadius( double r ) {
        bulbRadius = r;
    }

    void draw();
private:
    JaalViewer *viewManager;
    vector<GLenum> lightID;

    Point4F lightPos[8];
    int   status[8];

    bool   light_bulbs;
    double bulbRadius;
    void   setPositions();
    void   draw_bulb( int light);
    void   init();
};

typedef boost::shared_ptr<JLights> JLightsPtr;
