#include "SceneFloor.hpp"

/////////////////////////////////////////////////////////////////////////////////
JSceneFloor :: JSceneFloor()
{
}
/////////////////////////////////////////////////////////////////////////////////
void JSceneFloor :: draw()
{
    pattern = 1;
    if( texturePtr == nullptr) genTexture();

    glColor4fv( &color[0] );

    if( status ) {
        glDisable(GL_LIGHTING);
        glDisable(GL_CULL_FACE);

        if( pattern == 0)
            drawLines();
        else
            texturedFloor();
    }
}

/////////////////////////////////////////////////////////////////////////////////
void JSceneFloor :: setPattern(int p )
{
    pattern = p;
    if( pattern == 1 ) {
        genTexture();
    }

}

/////////////////////////////////////////////////////////////////////////////////
void JSceneFloor :: drawLines()
{
    double dl = floorLength/(double)(numFloorLines-1);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glColor3f( 1.0, 1.0, 0.0);
    glEnable (GL_LINE_STIPPLE);
    glLineStipple (1, 0x0101);   /*  dotted   */

    if( plane == 0) {
        glBegin( GL_LINES);
        for( int j = 0; j < numFloorLines; j++) {
            double d = -0.5*floorLength + j*dl;
            glVertex3f( floorDistance, -0.5*floorLength, d);
            glVertex3f( floorDistance, +0.5*floorLength, d);
            glVertex3f( floorDistance, d, -0.5*floorLength);
            glVertex3f( floorDistance, d, +0.5*floorLength);
        }
        glEnd();
    }

    if( plane == 1) {
        glBegin( GL_LINES);
        for( int j = 0; j < numFloorLines; j++) {
            double d = -0.5*floorLength + j*dl;
            glVertex3f(  -0.5*floorLength, floorDistance, d);
            glVertex3f(  +0.5*floorLength, floorDistance, d);
            glVertex3f(  d, floorDistance, -0.5*floorLength);
            glVertex3f(  d, floorDistance, +0.5*floorLength);
        }
        glEnd();
    }

    if( plane == 2) {
        glBegin( GL_LINES);
        for( int j = 0; j < numFloorLines; j++) {
            double d = -0.5*floorLength + j*dl;
            glVertex3f(  -0.5*floorLength, d, floorDistance);
            glVertex3f(  +0.5*floorLength, d, floorDistance);
            glVertex3f(  d, -0.5*floorLength, floorDistance);
            glVertex3f(  d, +0.5*floorLength, floorDistance);
        }
        glEnd();
    }

    glDisable (GL_LINE_STIPPLE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

///////////////////////////////////////////////////////////////////////////////
void JSceneFloor :: genTexture()
{
    texturePtr.reset( new JTexture );
    GLuint texID = texturePtr->getCheckerBoard(floorLength, floorLength, 8, 8);
}
///////////////////////////////////////////////////////////////////////////////

void JSceneFloor :: texturedFloor()
{
    double dl = floorLength/(double)(numFloorLines-1);
    double du = 1.0/(double)(numFloorLines-1);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if( plane == 0) {
        glBegin( GL_QUADS);
        for( int j = 0; j < numFloorLines; j++) {
            for( int i = 0; i < numFloorLines; i++) {
                double u = i*du;
                double v = j*du;
                double y = -0.5*floorLength + i*dl;
                double z = -0.5*floorLength + j*dl;
                glTexCoord2f( u, v);
                glVertex3f( floorDistance, y,    z);

                glTexCoord2f( u+du, v);
                glVertex3f( floorDistance, y+dl, z);

                glTexCoord2f( u+du, v+du);
                glVertex3f( floorDistance, y+dl, z+dl);

                glTexCoord2f( u, v+du);
                glVertex3f( floorDistance, y,    z+dl);
            }
        }
        glEnd();
    }

    if( plane == 1) {
        glBegin( GL_QUADS);
        for( int j = 0; j < numFloorLines; j++) {
            for( int i = 0; i < numFloorLines; i++) {
                double u = i*du;
                double v = j*du;
                double x = -0.5*floorLength + i*dl;
                double z = -0.5*floorLength + j*dl;
                glTexCoord2f( u, v);
                glVertex3f( x, floorDistance,  z);

                glTexCoord2f( u+du, v);
                glVertex3f( x + dl, floorDistance, z);

                glTexCoord2f( u+du, v+du);
                glVertex3f( x+ dl, floorDistance, z+dl);

                glTexCoord2f( u, v+du);
                glVertex3f( x, floorDistance,  z+dl);
            }
        }
        glEnd();
    }

    /*
        if( plane == 2) {
            glBegin( GL_LINES);
            for( int j = 0; j < numFloorLines; j++) {
                double d = -0.5*floorLength + j*dl;
                glVertex3f(  -0.5*floorLength, d, floorDistance);
                glVertex3f(  +0.5*floorLength, d, floorDistance);
                glVertex3f(  d, -0.5*floorLength, floorDistance);
                glVertex3f(  d, +0.5*floorLength, floorDistance);
            }
            glEnd();
        }
    */

}
///////////////////////////////////////////////////////////////////////////////
