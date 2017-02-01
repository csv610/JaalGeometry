#include <assert.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sstream>
#include <iostream>

#include "MeshViewer3.hpp"
#include "JaalViewer.hpp"

using namespace Jaal;
using namespace std;

#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/gle.h>

///////////////////////////////////////////////////////////////////////////////
// Static members
///////////////////////////////////////////////////////////////////////////////

log4cxx::LoggerPtr JMeshViewer ::logger = 0;
GLUquadricObj* DrawNode :: sphereObj = gluNewQuadric();

///////////////////////////////////////////////////////////////////////////////

Color MeshEntityColor :: getRandomColor()
{
    Color clr;
    double r = drand48();
    r  = max(0.2, r);
    r  = min(0.8, r);

    double g = drand48();
    g  = max(0.2, g);
    g  = min(0.8, g);

    double b = drand48();
    b  = max(0.2, b);
    b  = min(0.8, b);
    clr[0] = r;
    clr[1] = g;
    clr[2] = b;
    return clr;
}

///////////////////////////////////////////////////////////////////////////////

DrawNode ::  DrawNode()
{
    offset        = 0;
    offset_factor = 1.0;
    offset_unit   = 1.0;
    antiAlias    = 0;
    display_ids  = 0;
    display_flag = 1;
    display_normals = 0;
    glyphScale = 1.0;

    normalSign   = 1;
    normalLength = 1.0;
    normalsColor[0] = 0.6;
    normalsColor[1] = 0.6;
    normalsColor[2] = 0.6;
    normalsColor[3] = 0.0;

    pointSize = 1.0;
    sphRadius = 0.1;
    numSlices = 16;
    numStacks = 16 ;
    salientRadius = 0.2;
    glyph  = NODE_AS_POINT;
    use_display_list = 0;
    colorMethod  = nullptr;
    defaultColorMethod  = nullptr;

    gluQuadricDrawStyle( sphereObj, GLU_FILL);
}

///////////////////////////////////////////////////////////////////////////////

void DrawNode ::  preRender()
{
    lights = JLights::getInstance().isAnyOn();

    if( offset ) {
        glEnable( GL_POLYGON_OFFSET_POINT);
        glPolygonOffset( offset_factor, offset_unit);
    } else
        glDisable( GL_POLYGON_OFFSET_POINT);

    if( antiAlias ) {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glHint( GL_POINT_SMOOTH_HINT, GL_NICEST);
        glEnable( GL_POINT_SMOOTH );
    } else {
        glDisable (GL_BLEND);
        glDisable( GL_POINT_SMOOTH );
    }

    switch( glyph ) {
    case NODE_AS_POINT:
        glDisable(GL_COLOR_MATERIAL );
        JLights::getInstance().Switch(0);
        break;
    case NODE_AS_SPHERE:
        JLights::getInstance().Switch(1);
        glEnable(GL_COLOR_MATERIAL );
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glShadeModel(GL_SMOOTH );
        break;
    }
}
///////////////////////////////////////////////////////////////////////////////


void DrawNode ::  postRender()
{
    glDisable(GL_BLEND );
    glDisable(GL_POINT_SMOOTH );
    glDisable( GL_COLOR_MATERIAL );
    JLights::getInstance().Switch(lights);
    glDisable( GL_POLYGON_OFFSET_POINT);
    this->setGlyphScale(1.0);
}

///////////////////////////////////////////////////////////////////////////////

int DrawNode :: withName( const Vertex *vertex, size_t vid) const
{
    if( preConditions(vertex) == FAILURE ) return 1;

    glPushMatrix();
    glPushName(vid);
    draw(vertex);
    glPopName();
    glPopMatrix();
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool DrawNode :: preConditions( const Vertex *vertex) const
{
    if( vertex == nullptr ) return FAILURE;

    if( !vertex->isActive() ) return FAILURE;

    bool val = 0;
    vertex->getAttribute("Display", val);
    if( val == 0) return FAILURE;

    return SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////

int DrawNode :: draw( const Vertex *vertex ) const
{
    if( preConditions(vertex) == FAILURE ) return 1;

    if( setColor(vertex) == SUCCESS ) {
        if( glyph ) {
            const Point3D &xyz = vertex->getXYZCoords();
            glPushMatrix();
            glTranslatef( xyz[0], xyz[1], xyz[2] );
            gluSphere( sphereObj, glyphScale*sphRadius, numSlices, numStacks);
            glPopMatrix();
            return 0;
        }

        glPointSize( glyphScale*pointSize );
        glBegin(GL_POINTS);
        const Point3D &xyz = vertex->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2] );
        glEnd();
        return 0;
    }
    return 2;
}
///////////////////////////////////////////////////////////////////////////////

int DrawNode :: setColor(const Vertex *vertex) const
{
    int err = 0;
    if( colorMethod ) {
        Color clr;
        err = colorMethod->getColor(vertex, clr);
        glColor4fv( &clr[0] );
    }
    return err;
}

////////////////////////////////////////////////////////////////////////////////
int DrawNode :: draw_normal(const Vertex *vertex) const
{
    if( preConditions(vertex) == FAILURE ) return 1;
    assert( normalLength > 0.0);

    Vec3F normal;
    vertex->getAttribute("Normal", normal);

    Point3D xyz;
    double x,y,z;
    glBegin(GL_LINES);
    xyz = vertex->getXYZCoords();

    x = xyz[0];
    y = xyz[1];
    z = xyz[2];
    glVertex3f(x, y, z);

    x += normalLength * normalSign * normal[0];
    y += normalLength * normalSign * normal[1];
    z += normalLength * normalSign * normal[2];
    glVertex3f(x, y, z);
    glEnd();

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

DrawEdge :: DrawEdge()
{
    offset        = 0;
    offset_factor = 1.0;
    offset_unit   = 1.0;
    display_flag = 1;
    antiAlias   = 0;
    cylRadius  = 0.1;
    drawNode   = nullptr;
    glyphScale = 1.0;
    lineWidth  = 1.0;
    glyph  = EDGE_AS_LINE;
    use_display_list = 0;
    display_lower_nodes = 0;
    numCylSides = 20;

    colorMethod  = nullptr;
    defaultColorMethod  = nullptr;

    gleSetJoinStyle (TUBE_NORM_EDGE | TUBE_JN_ANGLE | TUBE_JN_CAP);
    gleSetJoinStyle( TUBE_JN_ROUND );
}
////////////////////////////////////////////////////////////////////////////////

void DrawEdge :: display_lower_entity( int e, bool v )
{
    if( e == 0) display_lower_nodes = v;
}

////////////////////////////////////////////////////////////////////////////////
bool DrawEdge :: display_lower_entity( int e) const
{
    if( e == 0) return display_lower_nodes;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void DrawEdge :: preRender()
{
    lights = JLights::getInstance().isAnyOn();

    if( offset ) {
        glEnable( GL_POLYGON_OFFSET_LINE);
        glPolygonOffset( offset_factor, offset_unit);
    } else
        glDisable( GL_POLYGON_OFFSET_LINE);

    if( antiAlias ) {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glHint( GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable( GL_LINE_SMOOTH );
    } else {
        glDisable (GL_BLEND);
        glDisable( GL_LINE_SMOOTH );
    }

    switch(glyph) {
    case EDGE_AS_LINE:
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable( GL_COLOR_MATERIAL );
        JLights::getInstance().Switch(0);
        break;
    case EDGE_AS_CYLINDER:
        glEnable( GL_COLOR_MATERIAL );
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glShadeModel(GL_SMOOTH );
        JLights::getInstance().Switch(1);
        break;
    }

}
////////////////////////////////////////////////////////////////////////////////

void DrawEdge :: postRender()
{
    JLights::getInstance().Switch(lights);

    glDisable( GL_BLEND );
    glDisable( GL_COLOR_MATERIAL );
    glDisable( GL_LINE_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable( GL_POLYGON_OFFSET_LINE);
    this->setGlyphScale(1.0);
}
////////////////////////////////////////////////////////////////////////////////

void DrawEdge :: drawCylinder( const Point3D &p0, const Point3D &p1)
{
    Point3D xyz;

    glPushMatrix();
    gleDouble endPoints[4][3];

    endPoints[0][0] = p0[0];
    endPoints[0][1] = p0[1];
    endPoints[0][2] = p0[2];

    endPoints[3][0] = p1[0];
    endPoints[3][1] = p1[1];
    endPoints[3][2] = p1[2];

    JMath::interpolate( p0, p1, xyz, 0.0001);
    endPoints[1][0] = xyz[0];
    endPoints[1][1] = xyz[1];
    endPoints[1][2] = xyz[2];

    JMath::interpolate( p0, p1, xyz, 0.9999);
    endPoints[2][0] = xyz[0];
    endPoints[2][1] = xyz[1];
    endPoints[2][2] = xyz[2];

    glePolyCylinder (4, endPoints, nullptr, cylRadius);

    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////

void DrawEdge :: drawCylinder( const Edge *edge, double radius ) const
{
    if( !edge->isActive() ) return;

    Point3D xyz;

    glPushMatrix();
    gleDouble endPoints[4][3];

    Vertex *n0 = edge->getNodeAt(0);
    xyz = n0->getXYZCoords();
    endPoints[0][0] = xyz[0];
    endPoints[0][1] = xyz[1];
    endPoints[0][2] = xyz[2];

    Vertex *n1 = edge->getNodeAt(1);
    xyz = n1->getXYZCoords();
    endPoints[3][0] = xyz[0];
    endPoints[3][1] = xyz[1];
    endPoints[3][2] = xyz[2];

    Vertex::mid_point( n0, n1, xyz, 0.0001);
    endPoints[1][0] = xyz[0];
    endPoints[1][1] = xyz[1];
    endPoints[1][2] = xyz[2];

    Vertex::mid_point( n0, n1, xyz, 0.9999);
    endPoints[2][0] = xyz[0];
    endPoints[2][1] = xyz[1];
    endPoints[2][2] = xyz[2];

    gleSetNumSides(numCylSides);
    glePolyCylinder (4, endPoints, nullptr, radius);

    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
int DrawEdge :: draw( const Edge *edge) const
{
    if( preConditions(edge) == FAILURE ) return 1;

    if( setColor(edge) == SUCCESS ) {
        if( glyph ) {
            JLights::getInstance().Switch(1);
            drawCylinder(edge, glyphScale*cylRadius);
        } else {
            Point3D xyz;
            glLineWidth(glyphScale*lineWidth);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glBegin(GL_LINES);
            Vertex *n0 = edge->getNodeAt(0);
            xyz = n0->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);

            Vertex *n1 = edge->getNodeAt(1);
            xyz = n1->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
            glEnd();
        }
    }

    if( display_lower_nodes ) {
        drawNode->preRender();
        drawNode->draw( edge->getNodeAt(0) );
        drawNode->draw( edge->getNodeAt(1) );
        drawNode->postRender();
    }

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int DrawEdge :: withName(const Edge *edge, size_t eid) const
{
    if( preConditions(edge) == FAILURE ) return 1;

    glPushMatrix();
    glPushName(eid);
    draw(edge);
    glPopName();
    glPopMatrix();
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
int  DrawEdge :: draw_id( const Edge *edge, size_t eid ) const
{
    if( preConditions(edge) == FAILURE ) return 1;

    Point3D xyz;
    char number[128];
    qglviewer::Vec wc;
    Vertex::mid_point( edge->getNodeAt(0), edge->getNodeAt(1), xyz );
    wc.x = xyz[0];
    wc.y = xyz[1];
    wc.z = xyz[2];
    glPushMatrix();
    glTranslatef(wc.x, wc.y, wc.z);
    glRotatef( rotFontAngles[0], 1.0, 0.0, 0.0);
    glRotatef( rotFontAngles[1], 0.0, 1.0, 0.0);
    glRotatef( rotFontAngles[2], 0.0, 0.0, 1.0);
    glScalef(fontScale, fontScale, fontScale);
    sprintf(number, "%ld", eid);
    font->Render(number);
    glPopMatrix();
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

bool DrawEdge :: preConditions( const Edge *edge) const
{
    if( edge  == nullptr ) return 0;
    if( !edge->isActive() ) return 0;
    bool val = 0;
    edge->getAttribute("Display", val);
    if( val == 0) return FAILURE;
    return ALL_PASSED;
}

////////////////////////////////////////////////////////////////////////////////

int DrawEdge :: setColor(const Edge *edge) const
{
    int err = 0;

    if(  colorMethod ) {
        Color clr;
        err = colorMethod->getColor( edge, clr);
        glColor4fv( &clr[0] );
    }

    return err;
}
////////////////////////////////////////////////////////////////////////////////

DrawFace :: DrawFace()
{
    offset        = 1;
    offset_factor = 1.0;
    offset_unit   = 1.0;
    backfaceCull = 0;
    display_lower_nodes = 0;
    display_lower_edges = 0;
    display_normals = 0;
    drawEdge = nullptr;
    drawNode = nullptr;
    normalSign  = 1;
    display_flag    = 0;
    normalsColor[0] = 0.6;
    normalsColor[1] = 0.6;
    normalsColor[2] = 0.6;
    use_normal  = 0;
    use_display_list = 0;
    style  = FACE_FILL;
    screendoorPattern = 8;
    transparencyMethod = OPAQUE;
    use_material   = 0;

    colorMethod   = nullptr;
    defaultColorMethod   = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void DrawFace :: display_lower_entity( int e, bool v )
{
    switch(e) {
    case 0:
        display_lower_nodes = v;
        break;
    case 1:
        display_lower_edges = v;
        break;
    }
}
////////////////////////////////////////////////////////////////////////////////

bool DrawFace :: display_lower_entity( int e) const
{
    switch(e) {
    case 0:
        return display_lower_nodes;
    case 1:
        return display_lower_edges;
        break;
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void DrawFace :: preRender()
{
    if( offset ) {
        glEnable( GL_POLYGON_OFFSET_FILL);
        glPolygonOffset( offset_factor, offset_unit);
    } else
        glDisable( GL_POLYGON_OFFSET_FILL);

    glEnable( GL_POLYGON_OFFSET_FILL);
    glEnable( GL_POLYGON_OFFSET_LINE);

    if( backfaceCull )
        glEnable( GL_CULL_FACE );
    else
        glDisable( GL_CULL_FACE );

    if( style == 0)
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);

    glEnable( GL_COLOR_MATERIAL );

    switch( transparencyMethod ) {
    case SCREENDOOR:
        glDisable( GL_BLEND );
        glEnable(GL_POLYGON_STIPPLE);
        glPolygonStipple(stippleMask[screendoorPattern]);  /* 0% opaqueness */
        break;
    case BLENDING:
        glDisable(GL_POLYGON_STIPPLE);
        glEnable( GL_BLEND );
        break;
    default:
        glDisable( GL_BLEND );
    }
}
////////////////////////////////////////////////////////////////////////////////
void DrawFace :: postRender()
{
    glDisable( GL_BLEND );
    glDisable( GL_CULL_FACE );
    glDisable( GL_COLOR_MATERIAL );
    glDisable(GL_POLYGON_STIPPLE);
    glDisable( GL_POLYGON_OFFSET_FILL);
}

////////////////////////////////////////////////////////////////////////////////

bool DrawFace :: preConditions( const Face *face) const
{
    if( face  == nullptr ) return FAILURE;
    if( !face->isActive() ) return FAILURE;

    bool val = 0;
    face->getAttribute("Display", val);
    if( val == 0) return FAILURE;
    return ALL_PASSED;
}

////////////////////////////////////////////////////////////////////////////////
int DrawFace :: setColor( const Face *face )const
{
    int err = 0;
    if( colorMethod ) {
        Color clr;
        err = colorMethod->getColor(face, clr);
        glColor4fv( &clr[0] );
    }
    return err;
}
////////////////////////////////////////////////////////////////////////////////

int DrawFace :: draw(const Face *face, int ori ) const
{
    if( preConditions(face) == FAILURE) return 1;

    if( setColor(face) != SUCCESS ) return 1;

    int nnodes = face->getSize(0);
    switch( nnodes ) {
    case 3:
        draw_tri(face, ori);
        break;
    case 4:
        draw_quad(face, ori);
        break;
    default:
        draw_poly(face, ori);
        break;
    }


    if( display_lower_edges ) {
        if( drawEdge) {
            drawEdge->preRender();
            int numedges = face->getSize(1);
            bool val = 1;
            for( int i = 0; i < numedges; i++)  {
                Edge *edge = face->getEdgeAt(i);
                edge->setAttribute("Display", val);
                drawEdge->draw(edge);
            }
            drawEdge->postRender();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int DrawFace :: withName(const Face *face, size_t fid) const
{
    if( preConditions(face) == FAILURE ) return 1;

    glPushName(fid);
    draw(face);
    glPopName();
    glPopMatrix();

    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void DrawFace :: draw_tri( const Face *face, int ori) const
{
    if( ori < 0) {
        glBegin(GL_TRIANGLES);
        for (int j = 0; j < 3; j++) {
            Vertex *vtx = face->getNodeAt(2-j);
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }
    glBegin(GL_TRIANGLES);
    for (int j = 0; j < 3; j++) {
        Vertex *vtx = face->getNodeAt(j);
        const Point3D &xyz = vtx->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
}
////////////////////////////////////////////////////////////////////////////////

void DrawFace :: draw_quad( const Face *face, int ori) const
{
    if( ori < 0) {
        glBegin(GL_QUADS);
        for (int j = 0; j < 4; j++) {
            Vertex *vtx = face->getNodeAt(3-j);
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }

    glBegin(GL_QUADS);
    for (int j = 0; j < 4; j++) {
        Vertex *vtx = face->getNodeAt(j);
        const Point3D &xyz = vtx->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
}

////////////////////////////////////////////////////////////////////////////////
void DrawFace :: draw_poly( const Face *face, int ori) const
{
    int nnodes = face->getSize(0);

    if( ori < 0) {
        glBegin(GL_POLYGON);
        for (int j = 0; j < nnodes; j++) {
            Vertex *vtx = face->getNodeAt(nnodes-1-j);
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }

    glBegin(GL_POLYGON);
    for (int j = 0; j < nnodes; j++) {
        Vertex *vtx = face->getNodeAt(j);
        const Point3D &xyz = vtx->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
}
////////////////////////////////////////////////////////////////////////////////


void DrawFace :: setNormal( const Face *face ) const
{
    Vec3F normal;
    face->getAttribute( "Normal", normal );
    float nx = 1.0 * normalSign * normal[0];
    float ny = 1.0 * normalSign * normal[1];
    float nz = 1.0 * normalSign * normal[2];
    glNormal3f(nx, ny, nz);
}

////////////////////////////////////////////////////////////////////////////////

int DrawFace :: draw_normal( const Face *face ) const
{
    if( preConditions(face) == FAILURE ) return 1;

    assert( normalLength > 0.0);

    Vec3F   normal;
    face->getAttribute("Normal", normal);

    Point3D centroid;
    face->getAvgPos( centroid );

    double x,y,z;
    glBegin( GL_LINES );
    x = centroid[0];
    y = centroid[1];
    z = centroid[2];
    glVertex3f(x, y, z);

    x += normalLength * normalSign * normal[0];
    y += normalLength * normalSign * normal[1];
    z += normalLength * normalSign * normal[2];
    glVertex3f(x, y, z);
    glEnd();

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

DrawCell :: DrawCell()
{
    display_flag        = 0;
    display_lower_nodes = 0;
    display_lower_edges = 0;
    display_lower_faces = 1;
    drawNode = nullptr;
    drawEdge = nullptr;
    drawFace = nullptr;
    colorMethod = nullptr;
    defaultColorMethod = nullptr;
}

/////////////////////////////////////////////////////////////////////////////

bool DrawCell :: preConditions( const Cell *cell) const
{
    if( cell  == nullptr ) return 0;
    if( !cell->isActive() ) return 0;

    bool val = 0;
    cell->getAttribute("Display", val);
    if( val == 0) return 0;

    return 1;
}
/////////////////////////////////////////////////////////////////////////////

void DrawCell :: display_lower_entity( int e, bool v )
{
    switch(e) {
    case 0:
        display_lower_nodes = v;
        break;
    case 1:
        display_lower_edges = v;
        break;
    case 2:
        display_lower_faces = v;
        break;
    }
}
////////////////////////////////////////////////////////////////////////////////

bool DrawCell :: display_lower_entity( int e) const
{
    switch(e) {
    case 0:
        return display_lower_nodes;
    case 1:
        return display_lower_edges;
    case 2:
        return display_lower_faces;
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////

int DrawCell :: setColor( const Cell *cell )const
{
    int err = 0;
    if( colorMethod ) {
        Color clr;
        err = colorMethod->getColor(cell, clr);
        glColor4fv( &clr[0] );
    }
    return err;
}

/////////////////////////////////////////////////////////////////////////////

int DrawCell :: draw(const Cell *cell)
{
    if( preConditions(cell) == 0) return 1;

    CellColor *clr = getColorMethod();

    bool val = 1;
    int ori;

    if( display_lower_faces  && drawFace != nullptr ) {
        drawFace->preRender();
        setColor(cell);
        cell->getFaces( faces );
        int numfaces = faces.size();
        for( int iface = 0; iface < numfaces; iface++)  {
            Face *face = cell->getFaceAt(iface, ori);
            if( face) {
                face->setAttribute("Display", val);
                drawFace->draw(face, ori);
            }
        }
        drawFace->postRender();
    }

    /*
         if( display_lower_edges ) {
              if( drawEdge) {
                   drawEdge->preRender();
                   int numedges = cell->getSize(1);
                   for( int i = 0; i < numedges; i++)  {
                        Edge *edge = cell->getEdgeAt(i);
                        edge->setAttribute("Display", val);
                        drawEdge->draw(edge);
                   }
                   drawEdge->postRender();
              }
         }

         if( display_lower_nodes ) {
              if( drawNode ) {
                   drawNode->preRender();
                   int numnodes = cell->getSize(0);
                   for( int i = 0; i < numnodes; i++)  {
                        Vertex *vertex= cell->getNodeAt(i);
                        vertex->setAttribute("Display", val);
                        drawNode->draw(vertex);
                   }
                   drawNode->postRender();
              }
         }
    */
    return 0;
}
/////////////////////////////////////////////////////////////////////////////

int DrawCell :: withName(const Cell *cell, int cid)
{
    if( preConditions(cell) == 0) return 1;

    glPushMatrix();
    glPushName(cid);
    draw(cell);
    glEnd();
    glPopName();
    glPopMatrix();
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
void JMeshViewer:: alignAlong( Mesh *mesh, const Vec3D &srcVec, const Vec3D &dstVec)
{
    Vec3D  perpAxis;
    JMath::cross_product( dstVec, srcVec, perpAxis);

    double angle = JMath::getVecAngle(srcVec, dstVec, ANGLE_IN_RADIANS);

    qglviewer::Vec rotaxis(perpAxis[0], perpAxis[1], perpAxis[2] );

    qglviewer::Quaternion quaternion(rotaxis, -1.0*angle);

    qglviewer::Vec prot;

    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        Vertex *v = mesh->getNodeAt(i);
        Point3D p3d = v->getXYZCoords();
        prot[0] = p3d[0];
        prot[1] = p3d[1];
        prot[2] = p3d[2];
        prot    = quaternion.rotate(prot);
        p3d[0]  = prot[0];
        p3d[1]  = prot[1];
        p3d[2]  = prot[2];
        v->setXYZCoords(p3d);
    }
}

/////////////////////////////////////////////////////////////////////////////////


void JMeshViewer :: alignAlong( Edge *currEdge, int along, bool refresh )
{
    if( currEdge == nullptr ) return;

    assert( mesh->getStatus() == Mesh::ACTIVE );

    const Point3D &p1 = currEdge->getNodeAt(0)->getXYZCoords();
    const Point3D &p2 = currEdge->getNodeAt(1)->getXYZCoords();

    AffineTransform af(mesh);
    af.translate(-p1[0], -p1[1], -p1[2] );

    // Specify the destination vector ...
    Vec3D dstVec;
    dstVec[0] = 0.0;
    dstVec[0] = 0.0;
    dstVec[2] = 0.0;
    dstVec[along] = 1.0;

    // Where is the vector now ...
    Point3D currVec;
    currVec[0] = p2[0] - p1[0];
    currVec[1] = p2[1] - p1[1];
    currVec[2] = p2[2] - p1[2];

    // Let the Quaternion rotate the "Current Vector" to "Destination Vector"
    alignAlong(mesh, currVec, dstVec);

    if( refresh ) refreshDisplay();
}
/////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: alignAlong( Face *currFace, int along, bool refresh )
{
    if( currFace == nullptr ) return;

    const Point3D &p1 = currFace->getNodeAt(0)->getXYZCoords();
    const Point3D &p2 = currFace->getNodeAt(1)->getXYZCoords();
    const Point3D &p3 = currFace->getNodeAt(2)->getXYZCoords();

    AffineTransform af(mesh);
    af.translate(-p1[0], -p1[1], -p1[2] );

    // Specify the destination vector ...
    Vec3D dstVec;
    dstVec[0] = 0.0;
    dstVec[1] = 0.0;
    dstVec[2] = 0.0;
    dstVec[along] = 1.0;

    // Where is the vector now ...
    Point3D vec1, vec2, currVec;

    vec1[0] = p2[0] - p1[0];
    vec1[1] = p2[1] - p1[1];
    vec1[2] = p2[2] - p1[2];

    vec2[0] = p3[0] - p1[0];
    vec2[1] = p3[1] - p1[1];
    vec2[2] = p3[2] - p1[2];

    JMath::cross_product( vec1, vec2, currVec);

    // Let the Quaternion rotate the "Current Vector" to "Destination Vector"
    alignAlong( mesh, currVec, dstVec);

    if( refresh ) refreshDisplay();
}

////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: resetAll( int entity, bool val)
{
    if( mesh == nullptr ) return;

    assert( mesh->getStatus() == Mesh::ACTIVE );

    size_t nSize;
    switch(entity) {
    case 0:
        nSize = mesh->getSize(0);
        for( size_t i = 0; i < nSize; i++) {
            Vertex *vtx = mesh->getNodeAt(i);
            vtx->setAttribute("Display", val);
        }
        break;
    case 1:
        nSize = mesh->getSize(1);
        for( size_t i = 0; i < nSize; i++) {
            Edge *edge = mesh->getEdgeAt(i);
            edge->setAttribute("Display", val);
        }
    case 2:
        nSize = mesh->getSize(2);
        for( size_t i = 0; i < nSize; i++) {
            Face *face = mesh->getFaceAt(i);
            face->setAttribute("Display", val);
        }
        break;
    case 3:
        nSize = mesh->getSize(3);
        for( size_t i = 0; i < nSize; i++) {
            Cell *cell= mesh->getCellAt(i);
            cell->setAttribute("Display", val);
        }
        break;
    }
}
/////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: resetAll( bool val)
{
    if( mesh == nullptr) return;
    resetAll(0, val);
    resetAll(1, val);
    resetAll(2, val);
    resetAll(3, val);
}

/////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: init_logger()
{
    BasicConfigurator::configure();
    JMeshViewer::logger = log4cxx::Logger::getLogger("JMeshViewer");
}

///////////////////////////////////////////////////////////////////////////////

JMeshViewer :: JMeshViewer( JaalViewer *vm)
{
    viewManager = vm;
    init();
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::init()
{
    mesh = nullptr;
    dualGraph = nullptr;
    minBox    = nullptr;

    drawNode = new DrawNode;
    drawEdge = new DrawEdge;
    drawFace = new DrawFace;
    drawCell = new DrawCell;

    meshNodeColor = new MeshNodeColor;
    meshEdgeColor = new MeshEdgeColor;
    meshFaceColor = new MeshFaceColor;
    meshCellColor = new MeshCellColor;

    drawEdge->setDrawNode( drawNode );

    drawFace->setDrawNode( drawNode );
    drawFace->setDrawEdge( drawEdge );

    drawCell->setDrawNode( drawNode );
    drawCell->setDrawEdge( drawEdge );
    drawCell->setDrawFace( drawFace );

    drawNode->setDefaultColorMethod( meshNodeColor );
    drawEdge->setDefaultColorMethod( meshEdgeColor );
    drawFace->setDefaultColorMethod( meshFaceColor );
    drawCell->setDefaultColorMethod( meshCellColor );

    edgemeshStyle = MESH_EDGES_LINE;
    facemeshStyle = DrawFace::FACE_FILL;
    hiddenlinesMethod = HIDDENLINE_WITH_DEPTH_TEST;

    display_entity[0] = 0;
    display_entity[1] = 0;
    display_entity[2] = 0;
    display_entity[3] = 0;

    display_boundary  = 0;
    display_enclosure = 0;
    currCounter  = 0;
    display_dual_graph  = 0;

    entityPicker =  new JMeshEntityPicker;
}

///////////////////////////////////////////////////////////////////////////////
JMeshViewer :: ~JMeshViewer()
{
    if( mesh != nullptr ) {
        mesh->deleteAll();
        delete mesh;
    }

    if( drawNode ) delete drawNode;
    if( drawEdge ) delete drawEdge;
    if( drawFace ) delete drawFace;
    if( drawCell ) delete drawCell;

    if( meshNodeColor ) delete meshNodeColor;
    if( meshEdgeColor ) delete meshEdgeColor;
    if( meshFaceColor ) delete meshFaceColor;
    if( meshCellColor ) delete meshCellColor;

    if( entityPicker  ) delete entityPicker;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: refreshDisplay()
{
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: loadNewMesh( const string &s)
{
    // Make sure that all the OpenGL events have finished.
    glFlush();

    meshFileName = s;
    if (!meshFileName.empty()) readData(meshFileName);
}
///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::saveMesh( const string &s)
{
    if (!s.empty()) {
        if( mesh ) mesh->saveAs( s );
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: animate()
{
    currCounter++;
//   QGLViewer :: animate();
}

///////////////////////////////////////////////////////////////////////////////
void
JMeshViewer::draw_bounding_box()
{
    JLights::getInstance().Switch(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glLineWidth(2.0);

    Point3D lower = box.getLower();
    Point3D upper = box.getUpper();

    float xmin = lower[0];
    float ymin = lower[1];
    float zmin = lower[2];

    float xmax = upper[0];
    float ymax = upper[1];
    float zmax = upper[2];

    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmax, ymin, zmin);
    glEnd();

    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmin, ymax, zmin);
    glEnd();

    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmin, ymin, zmax);
    glEnd();

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymax, zmin);

    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmin, ymax, zmin);

    glVertex3f(xmin, ymin, zmax);
    glVertex3f(xmax, ymin, zmax);

    glVertex3f(xmax, ymin, zmax);
    glVertex3f(xmax, ymax, zmax);

    glVertex3f(xmax, ymax, zmax);
    glVertex3f(xmin, ymax, zmax);

    glVertex3f(xmin, ymax, zmax);
    glVertex3f(xmin, ymin, zmax);

    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymin, zmax);

    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmin, ymax, zmax);

    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmax, ymax, zmax);
    glEnd();
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw_minimum_box()
{
    if( minBox == nullptr ) return;
    JLights::getInstance().Switch(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glLineWidth(2.0);
    vector<Point3D> hexpoints(8);
    for( int i = 0; i< 8; i++)
        hexpoints[i] = minBox->getNodeAt(i)->getXYZCoords();

    int v0, v1;
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    for( int i = 0; i < 12; i++) {
        Hexahedron::getEdgeTopology(i, v0, v1);
        glVertex3f(hexpoints[v0][0], hexpoints[v0][1], hexpoints[v0][2]);
        glVertex3f(hexpoints[v1][0], hexpoints[v1][1], hexpoints[v1][2]);
    }
    glEnd();
}


void
JMeshViewer::readData(const string &fname)
{
    if( mesh != nullptr ) {
        mesh->deleteAll();
        mesh->setStatus( Mesh::REMOVED );
        deletedmesh.push_back(mesh);
    }

    JNodeSequence nodes;
    JEdgeSequence hedges;

    mesh = new Jaal::Mesh();
    assert( mesh );
    mesh->readFromFile(fname);

    init_mesh();

}

///////////////////////////////////////////////////////////////////////////////

void Floor :: draw()
{
    if( active  == 0) return;

    JLights::getInstance().Switch(0);

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

/////////////////////////////////////////////////////////////////////////////////
size_t JMeshViewer :: getNumVisible( int entity )
{
    if( mesh->getStatus() != Mesh::ACTIVE ) return 1;

    bool   val = 0;
    size_t nCount = 0;
    size_t nSize;

    switch( entity ) {
    case 0:
        nSize = mesh->getSize(0);
        for( size_t i = 0; i < nSize; i++) {
            Vertex *vertex = mesh->getNodeAt(i);
            if( vertex->isActive() ) {
                vertex->getAttribute("Display", val) ;
                if( val ) nCount++;
            }
        }
        break;
    case 1:
        nSize = mesh->getSize(1);
        for( size_t i = 0; i < nSize; i++) {
            Edge *edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                edge->getAttribute("Display", val) ;
                if( val ) nCount++;
            }
        }
        break;
    case 2:
        nSize = mesh->getSize(2);
        for( size_t i = 0; i < nSize; i++) {
            Face *face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Display", val) ;
                if( val ) nCount++;
            }
        }
        break;
    case 3:
        nSize = mesh->getSize(3);
        for( size_t i = 0; i < nSize; i++) {
            Cell *cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                cell->getAttribute("Display", val) ;
                if( val ) nCount++;
            }
        }
        break;
    }
    return nCount;
}



///////////////////////////////////////////////////////////////////////////////
void
JMeshViewer::draw_nodes()
{
    if( display_entity[0] == 0 ) return;

    drawNode->preRender();

    static GLuint dlistID = 0;
    if( drawNode->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        }
        dlistID = glGenLists(1) + 1;
        glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
    }

    size_t numNodes = mesh->getSize(0);
    for (size_t j = 0; j < numNodes; j++) {
        Vertex *vtx = mesh->getNodeAt(j);
        drawNode->draw(vtx);
    }

    if( drawNode->displayList() ) {
        if( dlistID > 0) glEndList();
    }

    drawNode->postRender();
}

////////////////////////////////////////////////////////////////////////////////
void
JMeshViewer::draw_graph_nodes()
{
    if( dualGraph == nullptr ) return;

    drawNode->preRender();

    static GLuint dlistID = 0;
    if( drawNode->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        }
        dlistID = glGenLists(1) + 1;
        glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
    }

    size_t numNodes = dualGraph->getSize(0);
    for (size_t j = 0; j < numNodes; j++) {
        Vertex *vtx = dualGraph->getNodeAt(j);
        drawNode->draw(vtx);
    }

    if( drawNode->displayList() ) {
        if( dlistID > 0) glEndList();
    }

    drawNode->postRender();
}

////////////////////////////////////////////////////////////////////////////////



void JMeshViewer :: attach( JCurve *c )
{
    curves.push_back(c);
}

////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: detach( JCurve *c )
{
    vector<JCurve*>::iterator vend;
    vend = remove(curves.begin(), curves.end(), c);
    curves.erase(vend, curves.end());
}
////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw_curves()
{
//   if( display_entity[0] == 0 ) return;

    drawNode->preRender();
    JNodeSequence nodes;
    for( size_t i = 0; i < curves.size(); i++) {
        curves[i]->getNodes(nodes);
        size_t numNodes = nodes.size();
        for (size_t j = 0; j < numNodes; j++)
            drawNode->draw( nodes[j] );
    }
    drawNode->postRender();

    drawEdge->preRender();
    JEdgeSequence edges;
    for( size_t i = 0; i < curves.size(); i++) {
        curves[i]->getEdges(edges);
        size_t numEdges = edges.size();
        for (size_t j = 0; j < numEdges; j++) {
            drawEdge->draw( edges[j] );
        }
    }
    drawEdge->postRender();
}


////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw_edges()
{
    if( mesh == nullptr ) return;

    if( display_entity[1] == 0 ) return;

    if( edgemeshStyle == MESH_EDGES_HIDDENLINES_REMOVED) {
        draw_hidden_lines();
        return;
    }

    drawEdge->preRender();

    static GLuint dlistID = 0;
    if( drawEdge->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        }
        dlistID = glGenLists(1) + 1;
        glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
    }

    size_t numedges = mesh->getSize(1);

    if( numedges == 0 ) mesh->getTopology()->collect_edges();

    numedges = mesh->getSize(1);
    for (size_t i = 0; i < numedges; i++) {
        Edge *edge = mesh->getEdgeAt(i);
        drawEdge->draw(edge);
    }

    if( drawEdge->displayList() ) {
        if( dlistID > 0) glEndList();
    }

    drawEdge->postRender();
}

///////////////////////////////////////////////////////////////////////////////
void
JMeshViewer::draw_graph_edges()
{
    if( dualGraph == nullptr ) return;

    drawEdge->preRender();

    static GLuint dlistID = 0;
    if( drawEdge->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        }
        dlistID = glGenLists(1) + 1;
        glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
    }

    size_t numedges = dualGraph->getSize(1);

    numedges = dualGraph->getSize(1);
    for (size_t i = 0; i < numedges; i++) {
        Edge *edge = dualGraph->getEdgeAt(i);
        drawEdge->draw(edge);
    }

    if( drawEdge->displayList() ) {
        if( dlistID > 0) glEndList();
    }

    drawEdge->postRender();
}

///////////////////////////////////////////////////////////////////////////////


void
JMeshViewer::draw_nodes_normal()
{
    if( display_entity[0] == 0 ) return;

    if( !drawNode->displayNormals() ) return;

    JLights::getInstance().Switch(0);

    glDisable( GL_LINE_SMOOTH);

    const Color &nColor = drawNode->getNormalsColor();
    glColor3fv( &nColor[0] );

    glLineWidth(2.0);
    static GLuint dlistID = 0;

    if( drawNode->displayList()  ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        } else {
            dlistID = glGenLists(1) + 1;
            glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
        }
    }

    size_t numNodes = mesh->getSize(0);
    for (size_t j = 0; j < numNodes; j++) {
        Vertex *vtx = mesh->getNodeAt(j);
        drawNode->draw_normal(vtx);
    }

    if( drawNode->displayList()  )  {
        if( dlistID ) glEndList();
    }
}

////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw_faces_normal()
{
    if( display_entity[2] == 0 ) return;

    if( !drawFace->displayNormals()  ) return;

    JLights::getInstance().Switch(0);

    const Color &nColor = drawFace->getNormalsColor();

    glColor3fv( &nColor[0] );
    glDisable( GL_LINE_SMOOTH);

    glLineWidth(2.0);

    static GLuint dlistID = 0;

    if( drawFace->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        } else {
            dlistID = glGenLists(1) + 1;
            glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
        }
    }

    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        drawFace->draw_normal(face);
    }

    if( drawFace->displayList() )  {
        if( dlistID ) glEndList();
    }
}

////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: draw_all_faces()
{
    size_t numfaces = mesh->getSize(2);
    for (size_t i = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        drawFace->draw(face);
    }
}

////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: draw_filled_faces()
{
    /*
         float ambient[] = {0.19125, 0.0735, 0.0225, 1.0};
         float diffuse[] = {0.7038, 0.27048, 0.0828, 1.0};
         float specular[] = {1.0, 1.0, 1.0, 1.0};
         float shininess = 120;

         glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  ambient);
         glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  diffuse);
         glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
         glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
         JLights::getInstance().Switch(1);
    */

    glShadeModel(GL_FLAT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_COLOR_MATERIAL);

    static GLuint dlistID = 0;
    if( drawFace->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        } else {
            dlistID = glGenLists(1) + 1;
            glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
        }
    }

    draw_all_faces();

    glDisable(GL_COLOR_MATERIAL);

    if( dlistID > 0) glEndList();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: draw_wired_faces()
{
    JLights::getInstance().Switch(0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    static GLuint dlistID = 0;
    if( drawFace->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            return;
        } else {
            dlistID = glGenLists(1) + 1;
            glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
        }
    }

    draw_all_faces();

    if( dlistID > 0) glEndList();
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw_hidden_lines()
{
    if( display_entity[1] == 0 ) return;

    Color clr;
    clr  = meshEdgeColor->getInternalColor();
    glColor3f( clr[0], clr[1], clr[2] );

    int fstatus = display_entity[2];
    int filltype = drawFace->getStyle();

    display_entity[2] = 1;

    switch( hiddenlinesMethod ) {
    case HIDDENLINE_WITH_BACKFACES_CULL:
        drawFace->setBackfaceCull(1);
        drawFace->setStyle(DrawFace::FACE_LINES);
        draw_faces();
        drawFace->setBackfaceCull(0);
        break;
    case HIDDENLINE_WITH_FRONTLINES:
        clr  = viewManager->getBackgroundColor();
        glColor3f( clr[0], clr[1], clr[2] );
        glPolygonMode(GL_FRONT, GL_FILL);
        drawFace->setColorMethod( nullptr );
        drawFace->setStyle( DrawFace::FACE_FILL);
        draw_faces();

        clr  = meshEdgeColor->getInternalColor();
        glColor3f( clr[0], clr[1], clr[2] );
        drawFace->setStyle(DrawFace::FACE_LINES);
        draw_faces();
        break;
    case HIDDENLINE_WITH_DEPTH_TEST:
        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
        drawFace->setStyle( DrawFace::FACE_FILL);
        draw_faces();

        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
        glDepthFunc(GL_LEQUAL);
        drawFace->setStyle(DrawFace::FACE_LINES);
        draw_faces();
        glDepthFunc(GL_LESS);
        break;
    }

    display_entity[2] = fstatus;
    drawFace->setStyle( filltype );
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw_faces()
{
    if( display_entity[2] == 0 ) return;

    drawFace->preRender();

    size_t numfaces = mesh->getSize(2);

    if( numfaces == 0) {
        mesh->getTopology()->collect_faces();
    }

    numfaces = mesh->getSize(2);

    int faceStyle = drawFace->getStyle();

    switch( faceStyle ) {
    case DrawFace::FACE_LINES:
        draw_wired_faces();
        break;
    case DrawFace::FACE_FILL:
        draw_filled_faces();
        break;
    }

    drawFace->postRender();
}

///////////////////////////////////////////////////////////////////////////////
void
JMeshViewer::draw_cells()
{
    if( display_entity[3] == 0 ) return;

    size_t numcells = mesh->getSize(3);
    if( numcells < 1) return;

    drawCell->preRender();

    static GLuint dlistID = 0;

    if( drawCell->displayList() ) {
        if( dlistID > 0) {
            glCallList(dlistID);
            drawCell->postRender();
            return;
        } else {
            dlistID = glGenLists(1) + 1;
            glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
        }
    }

    for (size_t i = 0; i < numcells; i++) {
        Cell *cell = mesh->getCellAt(i);
        drawCell->draw(cell);
    }

    if( drawCell->displayList() ) {
        if( dlistID > 0) glEndList();
    }

    if( dlistID > 0) glEndList();

    drawCell->postRender();
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::drawWithNames()
{
    if( mesh == nullptr ) return;

    int pick_entity = entityPicker->getPickableEntity();

    if (pick_entity == 0) {
        size_t numnodes = mesh->getSize(0);
        for (size_t i = 0; i < numnodes; i++) {
            Vertex *vtx = mesh->getNodeAt(i);
            drawNode->withName( vtx, i );
        }
    }

    if (pick_entity == 1) {
        size_t numedges = mesh->getSize(1);
        for (size_t i = 0; i < numedges; i++) {
            Edge *edge = mesh->getEdgeAt(i);
            drawEdge->withName(edge,i);
        }
    }

    if (pick_entity == 2) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        size_t numfaces = mesh->getSize(2);
        for ( size_t i = 0; i < numfaces; i++) {
            Face *face = mesh->getFaceAt(i);
            drawFace->withName(face,i);
        }
    }

    if (pick_entity == 3) {
        size_t numcells = mesh->getSize(3);
        for (size_t i = 0; i < numcells; i++) {
            Cell *cell = mesh->getCellAt(i);
            drawCell->withName(cell,i);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: updateMouseReleasedEvent()
{
    if( entityPicker ) entityPicker->select_entity();
}

///////////////////////////////////////////////////////////////////////////////
int DrawNode :: draw_id( const Vertex *vertex, size_t vid ) const
{
    if( preConditions(vertex) == FAILURE ) return 1;

    qglviewer::Vec wc;
    char number[128];

    const Point3D  &xyz = vertex->getXYZCoords();
    wc.x = xyz[0];
    wc.y = xyz[1];
    wc.z = xyz[2];
    glPushMatrix();
    glTranslatef(wc.x, wc.y, wc.z);
    glRotatef( rotFontAngles[0], 1.0, 0.0, 0.0);
    glRotatef( rotFontAngles[1], 0.0, 1.0, 0.0);
    glRotatef( rotFontAngles[2], 0.0, 0.0, 1.0);
    glScalef(fontScale, fontScale, fontScale);
    sprintf(number, "%ld", vid);
    font->Render(number);
    glPopMatrix();
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int DrawFace :: draw_id( const Face *face, size_t fid ) const
{
    if( preConditions(face) == FAILURE) return 1;

    Point3D xyz;
    char number[128];
    qglviewer::Vec wc;
    wc.x = 0.0;
    wc.y = 0.0;
    wc.z = 0.0;
    int nsize = face->getSize(0);
    for (int j = 0; j < nsize; j++) {
        Vertex *n0 = face->getNodeAt(j);
        xyz = n0->getXYZCoords();
        wc.x += xyz[0];
        wc.y += xyz[1];
        wc.z += xyz[2];
    }
    wc.x /= (double) nsize;
    wc.y /= (double) nsize;
    wc.z /= (double) nsize;
    glPushMatrix();
    glTranslatef(wc.x, wc.y, wc.z + 0.001);
    glScalef(fontScale, fontScale, fontScale);
    sprintf(number, "%ld", fid);
    font->Render(number);
    glPopMatrix();

    return 0;

}

///////////////////////////////////////////////////////////////////////////////

void  JMeshViewer :: draw_ids()
{
    JLights::getInstance().Switch(0);

    glEnable( GL_CULL_FACE );

    float fontScale = 0.001*FontsManager::Instance().getFontScale();
    FTFont *font    = FontsManager::Instance().getFont(nullptr);
    Color color     = FontsManager::Instance().getColor();
    Point3F angles  = FontsManager::Instance().getRotateAngles();
    if( font == nullptr ) return;

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glColor3fv( &color[0] );
    if( drawNode->displayIDs() )  {
        drawNode->setFont( font, color, fontScale);
        drawNode->setRotateFonts( angles );
        int numnodes = mesh->getSize(0);
        for (int id = 0; id < numnodes; id++) {
            Vertex *vtx = mesh->getNodeAt(id);
            drawNode->draw_id(vtx, id);
        }
    }

    if( drawEdge->displayIDs() )  {
        drawEdge->setFont( font, color, fontScale);
        drawEdge->setRotateFonts( angles );
        int numedges = mesh->getSize(1);
        for (int id = 0; id < numedges; id++) {
            Edge *edge= mesh->getEdgeAt(id);
            drawEdge->draw_id(edge, id);
        }
    }

    if( drawFace->displayIDs() ) {
        drawFace->setFont( font, color, fontScale);
        drawFace->setRotateFonts( angles );
        for (size_t id = 0; id < mesh->getSize(2); id++) {
            Face *face = mesh->getFaceAt(id);
            drawFace->draw_id(face, id);
        }
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw_primal_mesh()
{
    draw_nodes();
    draw_edges();
    draw_faces();
    draw_cells();
}

///////////////////////////////////////////////////////////////////////////////
void
JMeshViewer::draw_dual_graph()
{
    if( display_dual_graph == 0) return;

    Jaal::DualGraph dgrapher;
    dualGraph = dgrapher.getGraph(mesh);

    draw_graph_nodes();
    draw_graph_edges();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshViewer :: draw_enclosure()
{
    if( !display_enclosure) return;

    switch( enclosure_type )
    {
    case AXIS_ALIGNED_BOUNDING_BOX:
        draw_bounding_box();
        break;
    case MINIMUM_BOUNDING_BOX:
        draw_minimum_box();
        break;
    }
}
///////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::draw()
{
    draw_curves();

    if(  mesh  ) {

        draw_enclosure();
        draw_ids();
        if( entityPicker) entityPicker->draw();
        draw_primal_mesh();
        draw_dual_graph();
        draw_nodes_normal();
        draw_faces_normal();
    }

}
////////////////////////////////////////////////////////////////////////////////
void JMeshViewer :: updateGeometry()
{
    if( mesh == nullptr ) return;

    if( enclosure_type == AXIS_ALIGNED_BOUNDING_BOX ) {
        box = mesh->getGeometry()->getBoundingBox();
    }

    if( enclosure_type == MINIMUM_BOUNDING_BOX ) {
        if( minBox ) delete minBox;
        minBox = mesh->getGeometry()->getMinimumBox();
    }

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();

    double minlen = mesh->getGeometry()->getMinEdgeLength();
    double maxlen = mesh->getGeometry()->getMaxEdgeLength();
    double avglen = 0.5*(maxlen + minlen );
    drawNode->setSphereRadius( 0.10*avglen);
    drawEdge->setCylinderRadius( 0.05*avglen);
    viewManager->refreshDisplay();
}

////////////////////////////////////////////////////////////////////////////////

void
JMeshViewer::init_mesh()
{
    if( mesh == nullptr ) return;

    drawNode->setColorMethod( meshNodeColor );
    drawEdge->setColorMethod( meshEdgeColor );
    drawFace->setColorMethod( meshFaceColor );
    drawCell->setColorMethod( meshCellColor );

    entityPicker->setMeshViewer(this);

    int topDim = mesh->getTopology()->getDimension();
    switch( topDim ) {
    case 3:
        display_entity[3] = 1;
        display_entity[1] = 1;
        display_entity[0] = 1;
        break;
    case 2:
        display_entity[2] = 1;
        display_entity[1] = 1;
        display_entity[0] = 1;
        break;
    case 1:
        display_entity[1] = 1;
        display_entity[0] = 1;
        break;
    case 0:
        display_entity[0] = 1;
        break;
    }
    updateGeometry();

    bool val = 1;
    size_t nSize;

    nSize = mesh->getSize(0);
    for( size_t i = 0; i < nSize; i++) {
        Vertex *v = mesh->getNodeAt(i);
        v->setAttribute("Display", val);
    }

    nSize = mesh->getSize(1);
    for( size_t i = 0; i < nSize; i++) {
        Edge *e = mesh->getEdgeAt(i);
        e->setAttribute("Display", val);
    }

    nSize = mesh->getSize(2);
    for( size_t i = 0; i < nSize; i++) {
        Face *f = mesh->getFaceAt(i);
        f->setAttribute("Display", val);
    }

    nSize = mesh->getSize(3);
    for( size_t i = 0; i < nSize; i++) {
        Cell *c = mesh->getCellAt(i);
        c->setAttribute("Display", val);
    }

    double minlen = mesh->getGeometry()->getMinEdgeLength();
    double maxlen = mesh->getGeometry()->getMaxEdgeLength();
    double avglen = 0.5*(maxlen + minlen );
    drawNode->setSphereRadius( 0.10*avglen );
    drawEdge->setCylinderRadius( 0.05*avglen );

    if( mesh->getTopology()->getDimension() == 2 ) {
        bool val = mesh->getTopology()->isClosed();   // Let there be atleast one light.
        JLights::getInstance().Switch(val);
    }
}

////////////////////////////////////////////////////////////////////////////////

