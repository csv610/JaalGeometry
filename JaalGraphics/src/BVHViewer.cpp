#include "BVHViewer.h"
#include <assert.h>

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>

static int stencilReflection = 1, stencilShadow = 1, offsetShadow = 1;
static int renderShadow = 1, renderReflection = 1;
static int linearFiltering = 0, useMipmaps = 0, floorTextured = 1, makeGround = 1;
static int directionalLight = 1;
static int useColor = 0;

/* Time varying or user-controled variables. */
static float lightAngle = 0.0, lightHeight = 50;

static GLfloat lightPosition[4];
static GLfloat lightColor[] = {0.8, 1.0, 0.8, 1.0};

enum
{
    X, Y, Z, W
};

enum
{
    A, B, C, D
};

static GLfloat floorVertices[4][3] = {
    { -500.0, 0.0, 500.0},
    { 500.0, 0.0, 500.0},
    { 500.0, 0.0, -500.0},
    { -500.0, 0.0, -500.0},
};

static GLfloat floorPlane[4];
static GLfloat floorShadow[4][4];

int texWidth = 16;
int texHeight = 16;

/////////////////////////////////////////////////////////////////////////////////

void
renderCylinder (float x1, float y1, float z1,
                float x2, float y2, float z2,
                float radius, int subdivisions, GLUquadricObj *quadric)
{
    float vx = x2 - x1;
    float vy = y2 - y1;
    float vz = z2 - z1;

    //handle the degenerate case of z1 == z2 with an approximation
    if (vz == 0)
        vz = .0001;

    float v = sqrt (vx * vx + vy * vy + vz * vz);
    float ax = 57.2957795 * acos (vz / v);
    if (vz < 0.0)
        ax = -ax;
    float rx = -vy*vz;
    float ry = vx*vz;
    glPushMatrix ();

    //draw the cylinder body
    glTranslatef (x1, y1, z1);
    glRotatef (ax, rx, ry, 0.0);
    gluQuadricOrientation (quadric, GLU_OUTSIDE);
    gluCylinder (quadric, radius, radius, v, subdivisions, 1);

    //draw the first cap
    gluQuadricOrientation (quadric, GLU_INSIDE);
    gluDisk (quadric, 0.0, radius, subdivisions, 1);
    glTranslatef (0, 0, v);

    //draw the second cap
    gluQuadricOrientation (quadric, GLU_OUTSIDE);
    gluDisk (quadric, 0.0, radius, subdivisions, 1);
    glPopMatrix ();
}

///////////////////////////////////////////////////////////////////////////////

void
renderCylinder (float *src, float *dst, float radius, int subdivisions)
{
    GLUquadricObj *quadric = gluNewQuadric ();
    gluQuadricNormals (quadric, GLU_SMOOTH);

    float x1 = src[0];
    float y1 = src[1];
    float z1 = src[2];
    float x2 = dst[0];
    float y2 = dst[1];
    float z2 = dst[2];

    renderCylinder (x1, y1, z1, x2, y2, z2, radius, subdivisions, quadric);

    gluDeleteQuadric (quadric);
}

///////////////////////////////////////////////////////////////////////////////

void
drawTranslateAgent (float *pos, float len)
{
    GLUquadricObj *quadric = gluNewQuadric ();
    /*
      glPointSize (5.0);
      glColor3f (1.0, 1.0, 1.0);
      glBegin (GL_POINTS);
      glVertex3f (pos[0], pos[1], pos[2]);
      glEnd ();
    */

    float radius = 1.0, top[3];
    glEnable( GL_LIGHTING);
    glColor3f (1.0, 0.0, 0.0);
    top[0] = pos[0] + len;
    top[1] = pos[1] ;
    top[2] = pos[2];
    renderCylinder( pos, top, radius, 20);
    glDisable( GL_LIGHTING);

    /*
      glLineWidth (2.0);

      glColor3f (0.0, 1.0, 0.0);
      glBegin (GL_LINES);
      glVertex3f (pos[0], pos[1], pos[2]);
      glVertex3f (pos[0], pos[1] + len, pos[2]);
      glEnd ();

      glColor3f (0.0, 0.0, 1.0);
      glBegin (GL_LINES);
      glVertex3f (pos[0], pos[1], pos[2]);
      glVertex3f (pos[0], pos[1], pos[2] + len);
      glEnd ();
    */


    glLineWidth (1.0);
    glPointSize (1.0);

    gluDeleteQuadric (quadric);
}


static
void
makeCheckImage (void)
{
    GLubyte checkImage[texWidth][texHeight][3];
    int i, j, c;

    for (i = 0; i < texHeight; i++)
    {
        for (j = 0; j < texWidth; j++)
        {
            if (floorTextured)
                c = ((((i & 0x8) == 0)^((j & 0x8)) == 0))*255;
            else
                c = 200;

            checkImage[i][j][0] = (GLubyte) c;
            checkImage[i][j][1] = (GLubyte) c;
            checkImage[i][j][2] = (GLubyte) c;
        }
    }

    glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    if (useMipmaps)
    {
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                         GL_LINEAR_MIPMAP_LINEAR);
        gluBuild2DMipmaps (GL_TEXTURE_2D, 3, texWidth, texHeight,
                           GL_RGB, GL_UNSIGNED_BYTE, checkImage);
        return;
    }

    if (linearFiltering)
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    else
        glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glTexImage2D (GL_TEXTURE_2D, 0, 3, texWidth, texHeight, 0,
                  GL_RGB, GL_UNSIGNED_BYTE, checkImage);
}

/* Create a matrix that will project the desired shadow. */
void
shadowMatrix (GLfloat shadowMat[4][4],
              GLfloat groundplane[4],
              GLfloat lightpos[4])
{
    GLfloat dot;

    /* Find dot product between light position vector and ground plane normal. */
    dot = groundplane[X] * lightpos[X] +
          groundplane[Y] * lightpos[Y] +
          groundplane[Z] * lightpos[Z] +
          groundplane[W] * lightpos[W];

    shadowMat[0][0] = dot - lightpos[X] * groundplane[X];
    shadowMat[1][0] = 0.f - lightpos[X] * groundplane[Y];
    shadowMat[2][0] = 0.f - lightpos[X] * groundplane[Z];
    shadowMat[3][0] = 0.f - lightpos[X] * groundplane[W];

    shadowMat[X][1] = 0.f - lightpos[Y] * groundplane[X];
    shadowMat[1][1] = dot - lightpos[Y] * groundplane[Y];
    shadowMat[2][1] = 0.f - lightpos[Y] * groundplane[Z];
    shadowMat[3][1] = 0.f - lightpos[Y] * groundplane[W];

    shadowMat[X][2] = 0.f - lightpos[Z] * groundplane[X];
    shadowMat[1][2] = 0.f - lightpos[Z] * groundplane[Y];
    shadowMat[2][2] = dot - lightpos[Z] * groundplane[Z];
    shadowMat[3][2] = 0.f - lightpos[Z] * groundplane[W];

    shadowMat[X][3] = 0.f - lightpos[W] * groundplane[X];
    shadowMat[1][3] = 0.f - lightpos[W] * groundplane[Y];
    shadowMat[2][3] = 0.f - lightpos[W] * groundplane[Z];
    shadowMat[3][3] = dot - lightpos[W] * groundplane[W];

}

/* Find the plane equation given 3 points. */
void
findPlane (GLfloat plane[4],
           GLfloat v0[3], GLfloat v1[3], GLfloat v2[3])
{
    GLfloat vec0[3], vec1[3];

    /* Need 2 vectors to find cross product. */
    vec0[X] = v1[X] - v0[X];
    vec0[Y] = v1[Y] - v0[Y];
    vec0[Z] = v1[Z] - v0[Z];

    vec1[X] = v2[X] - v0[X];
    vec1[Y] = v2[Y] - v0[Y];
    vec1[Z] = v2[Z] - v0[Z];

    /* find cross product to get A, B, and C of plane equation */
    plane[A] = vec0[Y] * vec1[Z] - vec0[Z] * vec1[Y];
    plane[B] = -(vec0[X] * vec1[Z] - vec0[Z] * vec1[X]);
    plane[C] = vec0[X] * vec1[Y] - vec0[Y] * vec1[X];

    plane[D] = -(plane[A] * v0[X] + plane[B] * v0[Y] + plane[C] * v0[Z]);
}

void
BVHViewer::animate ()
{
    animbody->nextFrame ();
}

void
BVHViewer::drawObjects (void)
{
    glDisable( GL_LIGHTING );
    float pos[3];
    pos[0] = 0.0;
    pos[1] = 2.0;
    pos[2] = 0.0;
    drawTranslateAgent( pos, 5.0);
    glEnable( GL_LIGHTING );

    /*
      if (animbody == nullptr) return;

      int id = selectedName ();
      animbody->setSelectedFrame(id);
      animbody->displayTrace (1);
      animbody->draw ();

      if (useColor)
        glColor3f (1.0, 0.0, 0.0);
    */
}

static void
drawFloor (void)
{
    static bool hasDisplayList  = 0;
    static GLuint displayListID = 0;

    if (!makeGround) return;

    glDisable (GL_LIGHTING);

    if (floorTextured) glEnable (GL_TEXTURE_2D);

    if( !hasDisplayList )  {
        displayListID = glGenLists(1);
        glNewList( displayListID, GL_COMPILE_AND_EXECUTE);
        glBegin (GL_QUADS);
        glTexCoord2f (0.0, 0.0);
        glVertex3fv (floorVertices[0]);
        glTexCoord2f (0.0, texHeight);
        glVertex3fv (floorVertices[1]);
        glTexCoord2f (texWidth, texHeight);
        glVertex3fv (floorVertices[2]);
        glTexCoord2f (texWidth, 0.0);
        glVertex3fv (floorVertices[3]);
        glEnd ();
        glEndList();
        hasDisplayList = 1;
    } else {
        glCallList( displayListID );
    }

    if (floorTextured) glDisable (GL_TEXTURE_2D);

    glEnable (GL_LIGHTING);
}

void
BVHViewer::draw ()
{
    useColor = 0;
    if ((stencilReflection && renderReflection) || (stencilShadow && renderShadow))
        glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    else
        glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /* Reposition the light source. */
    lightPosition[0] = 0.0* cos (lightAngle);
    lightPosition[1] = lightHeight;
    lightPosition[2] = 0.0 * sin (lightAngle);
    if (directionalLight)
    {
        lightPosition[3] = 0.0;
    }
    else
    {
        lightPosition[3] = 1.0;
    }

    shadowMatrix (floorShadow, floorPlane, lightPosition);

    glPushMatrix ();
    /* Tell GL new light source position. */
    glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);

    if (renderReflection)
    {
        if (stencilReflection)
        {
            /* We can eliminate the visual "artifact" of seeing the "flipped"
               dinosaur underneath the floor by using stencil.  The idea is
               draw the floor without color or depth update but so that
               a stencil value of one is where the floor will be.  Later when
               rendering the dinosaur reflection, we will only update pixels
               with a stencil value of 1 to make sure the reflection only
               lives on the floor, not below the floor. */

            /* Don't update color or depth. */
            glDisable (GL_DEPTH_TEST);
            glColorMask (GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

            /* Draw 1 into the stencil buffer. */
            glEnable (GL_STENCIL_TEST);
            glStencilOp (GL_REPLACE, GL_REPLACE, GL_REPLACE);
            glStencilFunc (GL_ALWAYS, 1, 0xffffffff);

            /* Now render floor; floor pixels just get their stencil set to 1. */
            drawFloor ();

            /* Re-enable update of color and depth. */
            glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
            glEnable (GL_DEPTH_TEST);

            /* Now, only render where stencil is set to 1. */
            glStencilFunc (GL_EQUAL, 1, 0xffffffff); /* draw if ==1 */
            glStencilOp (GL_KEEP, GL_KEEP, GL_KEEP);
        }

        glPushMatrix ();

        /* The critical reflection step: Reflect dinosaur through the floor
           (the Y=0 plane) to make a relection. */
        glScalef (1.0, -1.0, 1.0);

        /* Reflect the light position. */
        glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);

        /* To avoid our normals getting reversed and hence botched lighting
           on the reflection, turn on normalize.  */
        glEnable (GL_NORMALIZE);
        glCullFace (GL_FRONT);

        /* Draw the reflected dinosaur. */

        drawObjects ();

        /* Disable noramlize again and re-enable back face culling. */
        glDisable (GL_NORMALIZE);
        glCullFace (GL_BACK);

        glPopMatrix ();

        /* Switch back to the unreflected light position. */
        glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);

        if (stencilReflection)
        {
            glDisable (GL_STENCIL_TEST);
        }
    }

    /* Back face culling will get used to only draw either the top or the
       bottom floor.  This let's us get a floor with two distinct
       appearances.  The top floor surface is reflective and kind of red.
       The bottom floor surface is not reflective and blue. */

    /* Draw "bottom" of floor in blue. */
    glFrontFace (GL_CW); /* Switch face orientation. */
    glColor4f (0.1, 0.1, 0.7, 1.0);
    drawFloor ();
    glFrontFace (GL_CCW);

    if (renderShadow)
    {
        if (stencilShadow)
        {
            /* Draw the floor with stencil value 3.  This helps us only
               draw the shadow once per floor pixel (and only on the
               floor pixels). */
            glEnable (GL_STENCIL_TEST);
            glStencilFunc (GL_ALWAYS, 3, 0xffffffff);
            glStencilOp (GL_KEEP, GL_KEEP, GL_REPLACE);
        }
    }

    /* Draw "top" of floor.  Use blending to blend in reflection. */
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f (0.7, 0.0, 0.0, 0.3);
    glColor4f (1.0, 1.0, 1.0, 0.3);
    drawFloor ();
    glDisable (GL_BLEND);

    if (renderShadow)
    {

        /* Render the projected shadow. */
        if (stencilShadow)
        {

            /* Now, only render where stencil is set above 2 (ie, 3 where
               the top floor is).  Update stencil with 2 where the shadow
               gets drawn so we don't redraw (and accidently reblend) the
               shadow). */
            glStencilFunc (GL_LESS, 2, 0xffffffff); /* draw if ==1 */
            glStencilOp (GL_REPLACE, GL_REPLACE, GL_REPLACE);
        }

        /* To eliminate depth buffer artifacts, we use polygon offset
           to raise the depth of the projected shadow slightly so
           that it does not depth buffer alias with the floor. */
        if (offsetShadow)
            glEnable (GL_POLYGON_OFFSET_FILL);

        /* Render 50% black shadow color on top of whatever the
           floor appareance is. */
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDisable (GL_LIGHTING); /* Force the 50% black. */
        glColor4f (0.0, 0.0, 0.0, 0.5);

        glPushMatrix ();
        /* Project the shadow. */
        glMultMatrixf ((GLfloat *) floorShadow);
        drawObjects ();
        glPopMatrix ();

        glDisable (GL_BLEND);
        glEnable (GL_LIGHTING);

        if (offsetShadow)
            glDisable (GL_POLYGON_OFFSET_FILL);

        if (stencilShadow)
            glDisable (GL_STENCIL_TEST);
    }

#ifdef CSV
    glPushMatrix ();
    glDisable (GL_LIGHTING);
    glColor3f (1.0, 1.0, 0.0);
    if (directionalLight)
    {
        /* Draw an arrowhead. */
        glDisable (GL_CULL_FACE);
        glTranslatef (lightPosition[0], lightPosition[1], lightPosition[2]);
        glRotatef (lightAngle * -180.0 / M_PI, 0, 1, 0);
        glRotatef (atan (lightHeight / 12) * 180.0 / M_PI, 0, 0, 1);
        glBegin (GL_TRIANGLE_FAN);
        glVertex3f (0, 0, 0);
        glVertex3f (2, 1, 1);
        glVertex3f (2, -1, 1);
        glVertex3f (2, -1, -1);
        glVertex3f (2, 1, -1);
        glVertex3f (2, 1, 1);
        glEnd ();
        /* Draw a white line from light direction. */
        glColor3f (1.0, 1.0, 1.0);
        glBegin (GL_LINES);
        glVertex3f (0, 0, 0);
        glVertex3f (5, 0, 0);
        glEnd ();
        glEnable (GL_CULL_FACE);
    }
    else
    {
        /* Draw a yellow ball at the light source. */
        glTranslatef (lightPosition[0], lightPosition[1], lightPosition[2]);
    }
    glEnable (GL_LIGHTING);
    glPopMatrix ();
#endif


    glPushMatrix ();
    useColor = 1;
    drawObjects ();
    glPopMatrix ();

    glPopMatrix ();

    glFinish ();
}

void
BVHViewer::init ()
{
    // Restore previous viewer state.
    restoreStateFromFile ();

    camera ()->setZClippingCoefficient (500);
    camera ()->setZNearCoefficient (0.00001);

    glPolygonOffset (-2.0, -1.0);

    glEnable (GL_CULL_FACE);
    glEnable (GL_DEPTH_TEST);
    glEnable (GL_TEXTURE_2D);
    glLineWidth (3.0);

    glLightModeli (GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
    glLightfv (GL_LIGHT0, GL_DIFFUSE, lightColor);
    glLightf (GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0.1);
    glLightf (GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.05);
    glEnable (GL_LIGHT0);
    glEnable (GL_LIGHTING);

    makeCheckImage ();

    /* Setup floor plane for projected shadow calculations. */
    findPlane (floorPlane, floorVertices[1], floorVertices[2], floorVertices[3]);

    //setAnimationPeriod( 10 );

    // Opens help window
    help ();
}

QString
BVHViewer::helpString () const
{
    QString text ("<h2>S i m p l e V i e w e r</h2>");
    text += "Use the mouse to move the camera around the object. ";
    text += "You can respectively revolve around, zoom and translate with the three mouse buttons. ";
    text += "Left and middle buttons pressed together rotate around the camera view direction axis<br><br>";
    text += "Pressing <b>Alt</b> and one of the function keys (<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
    text += "Simply press the function key again to restore it. Several keyFrames define a ";
    text += "camera path. Paths are saved when you quit the application and restored at next start.<br><br>";
    text += "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
    text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save a snapshot. ";
    text += "See the <b>Keyboard</b> tab in this window for a complete shortcut list.<br><br>";
    text += "Double clicks automates single click actions: A left button double click aligns the closer axis with the camera (if close enough). ";
    text += "A middle button double click fits the zoom of the camera and the right button re-centers the scene.<br><br>";
    text += "A left button double click while holding right button pressed defines the camera <i>Revolve Around Point</i>. ";
    text += "See the <b>Mouse</b> tab and the documentation web pages for details.<br><br>";
    text += "Press <b>Escape</b> to exit the viewer.";
    return text;
}

#include <qapplication.h>

int
main (int argc, char** argv)
{
    assert (argc == 2);
    QApplication application (argc, argv);

    BVHViewer viewer;
    viewer.setWindowTitle ("BVHViewer");
    viewer.loadFile (argv[1]);

    viewer.show ();
    return application.exec ();
}
