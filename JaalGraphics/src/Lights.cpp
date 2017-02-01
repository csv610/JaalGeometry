#include "MeshViewer.hpp"
#include "Lights.hpp"

void JLights :: init()
{
    viewManager = nullptr;
    light_bulbs = 0;
    bulbRadius = 0.1;

    lightID.resize(8);
    lightID[0] = GL_LIGHT0;
    lightID[1] = GL_LIGHT1;
    lightID[2] = GL_LIGHT2;
    lightID[3] = GL_LIGHT3;
    lightID[4] = GL_LIGHT4;
    lightID[5] = GL_LIGHT5;
    lightID[6] = GL_LIGHT6;
    lightID[7] = GL_LIGHT7;

    Point4F lp;
    lp[0] = 1.0;
    lp[1] = 1.0;
    lp[2] = 1.0;
    lp[3] = 0.0;
    for( int i = 0; i < 8; i++) {
        lightPos[i] = lp;
        status[i]   = 0;
    }

    lp[0] = 0.0;
    lp[1] = 0.0;
    lp[2] = 1.0;
    setPosition(0, lp);

    lp[0] =  0.0;
    lp[1] =  0.0;
    lp[2] = -1.0;
    setPosition(1, lp);

    lp[0] =  1.0;
    lp[1] =  0.0;
    lp[2] =  0.0;
    setPosition(2, lp);

    lp[0] = -1.0;
    lp[1] =  0.0;
    lp[2] =  0.0;
    setPosition(3, lp);

    status[0] = 1;
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );

}

////////////////////////////////////////////////////////////////////////////////

Point4F JLights :: getPosition( int lnum )
{
    GLfloat pos[4];
    GLenum light = lightID[lnum];
    glGetLightfv( light, GL_POSITION, pos);

    Point4F p;
    p[0] = pos[0];
    p[1] = pos[1];
    p[2] = pos[2];
    p[3] = pos[3];
    return p;
}

////////////////////////////////////////////////////////////////////////////////

int JLights :: getNumLights() const
{
    int nCount = 0;
    for( int i = 0; i < 8; i++) {
        if(status[i]) nCount++;
    }
    return nCount;
}

///////////////////////////////////////////////////////////////////////////////
void JLights::setPositions()
{
    for( int i = 0; i < 8; i++) {
        GLenum light = lightID[i];
        if( glIsEnabled( light ) )
            glLightfv( light, GL_POSITION, &lightPos[i][0] ) ;
    }
}
///////////////////////////////////////////////////////////////////////////////
void JLights :: Switch(bool v)
{
    if(v)
        glEnable(GL_LIGHTING);
    else
        glDisable(GL_LIGHTING);

    for( int inum = 0; inum < 8; inum++) {
        GLenum light = lightID[inum];
        if( status[inum] && v )
            glEnable( light );
        else
            glDisable( light );
    }
}
///////////////////////////////////////////////////////////////////////////////

void JLights :: setLocalViewer( bool v )
{
    glLightModelf( GL_LIGHT_MODEL_LOCAL_VIEWER, v);
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setTwoSided( bool v )
{
    glLightModelf( GL_LIGHT_MODEL_TWO_SIDE, v);
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setAmbientModel( const JColor &color)
{
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, &color[0] );
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setAmbientColor( int lnum, const JColor &color)
{
    GLenum light = lightID[lnum];
    glLightfv( light, GL_AMBIENT, &color[0] );
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setDiffuseColor( int lnum, const JColor &color)
{
    GLenum light = lightID[lnum];
    glLightfv( light, GL_DIFFUSE, &color[0] );
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setSpecularColor( int lnum, const JColor &color)
{
    GLenum light = lightID[lnum];
    glLightfv( light, GL_SPECULAR, &color[0] );
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setPosition( int lnum,  const Point4F &pos)
{
    GLenum light = lightID[lnum];
    glLightfv( light, GL_POSITION, &pos[0] );

    lightPos[lnum] = pos;
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setSpotDirection( int lnum,  const Point3F &dir)
{
    GLenum light = lightID[lnum];
    glLightfv( light, GL_SPOT_DIRECTION, &dir[0] );
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setSpotExponent( int lnum,  double val)
{
    GLenum light = lightID[lnum];
    glLightf( light, GL_SPOT_EXPONENT, val);
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setSpotCutoff( int lnum,  double angle)
{
    GLenum light = lightID[lnum];
    glLightf( light, GL_SPOT_CUTOFF, angle);
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setConstantAttenuation( int lnum,  double val)
{
    GLenum light = lightID[lnum];
    glLightf( light, GL_CONSTANT_ATTENUATION, val);
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setLinearAttenuation( int lnum,  double angle)
{
    GLenum light = lightID[lnum];
    glLightf( light, GL_LINEAR_ATTENUATION, angle);
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: setQuadraticAttenuation( int lnum,  double val)
{
    GLenum light = lightID[lnum];
    glLightf( light, GL_QUADRATIC_ATTENUATION, val);
    if( viewManager ) viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLights :: draw()
{
    for( int i = 0; i < 8; i++)
        draw_bulb( i );
}

///////////////////////////////////////////////////////////////////////////////

void JLights :: draw_bulb( int lnum)
{
    if( status[lnum] ) {
        GLenum light = lightID[lnum];

        glDisable( GL_LIGHTING ) ;

        GLfloat color[4];
        glGetLightfv( light, GL_DIFFUSE, color );
        glColor3fv(color);

        GLfloat pos[4];
        glGetLightfv( light, GL_POSITION, pos);

        glPushMatrix();
        glTranslatef( pos[0], pos[1], pos[2] );
//        gluSphere( JNodeDraw::sphereObj, bulbRadius, 20, 20);
        glPopMatrix();

        glEnable( GL_LIGHTING ) ;
    }
}

///////////////////////////////////////////////////////////////////////////////
