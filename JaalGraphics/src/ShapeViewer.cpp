#include "ShapeViewer.hpp"


JShapeViewer :: JShapeViewer( JaalViewer *vm)
{
    viewManager = vm;
    /*    entityPicker.reset( new JEntityPicker()); */
    qObj = gluNewQuadric();
}

/////////////////////////////////////////////////////////////////////////////////

JShapeViewer :: ~JShapeViewer()
{
}

///////////////////////////////////////////////////////////////////////////////

JViewComponentPtr JShapeViewer :: registerComponent(JaalViewer *viewer)
{
    boost::shared_ptr<JShapeViewer> obj;
    obj.reset(new JShapeViewer(viewer));
    obj->setName("ShapeViewer");
    viewer->attach(obj);
    return obj;
}

///////////////////////////////////////////////////////////////////////////////
void
JShapeViewer::refreshDisplay()
{

}
///////////////////////////////////////////////////////////////////////////////

void
JShapeViewer::drawWithNames()
{
    if( !isActive()  ) return;

/*
    drawSpheresWithNames();
    drawBoxesWithNames();
    drawCylindersWithNames();
    drawTubesWithNames();
*/
}

//////////////////////////////////////////////////////////////////////////////////

void JShapeViewer::drawObject( const JSpherePtr &obj)
{
    Point3D xyz      = obj->getCenter();
    double radius    = obj->getRadius();
    JSphereRenderPtr attrib;
    obj->getAttribute("Render", attrib);
    int    numSlices = attrib->numSlices;
    int    numStacks = attrib->numStacks;
    JColor  color    = attrib->color;

/*
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColor3fv(  &color[0] );
    glPushMatrix();
    glTranslatef(xyz[0], xyz[1], xyz[2]);
    gluSphere( qObj, radius, numSlices, numStacks);
    glPopMatrix();
*/
}

///////////////////////////////////////////////////////////////////////////////
void JShapeViewer :: drawCylinder( const float *p0, const float *p1, float radius) const
{
/*
    glPushMatrix();
    gleDouble endPoints[4][3];

    endPoints[0][0] = p0[0];
    endPoints[0][1] = p0[1];
    endPoints[0][2] = p0[2];

    endPoints[1][0] = p1[0];
    endPoints[1][1] = p1[1];
    endPoints[1][2] = p1[2];

    endPoints[2][0] = p0[0];
    endPoints[2][1] = p0[1];
    endPoints[2][2] = p0[2];

    endPoints[3][0] = p1[0];
    endPoints[3][1] = p1[1];
    endPoints[3][2] = p1[2];

    glePolyCylinder(4, endPoints, nullptr, radius);
    glPopMatrix();
*/
}

///////////////////////////////////////////////////////////////////////////////

void JShapeViewer::drawObject( const JCylinderPtr  &obj)
{
/*
        Point3D xyz      = obj->getCenter();
        double radius    = obj->getRadius();
        double height    = obj->getHeight();
        int    numSlices = obj->getNumSlices();
        Vec3D  axis      = obj->getAxis();

        Point3F p0;
        p0[0] = xyz[0];
        p0[1] = xyz[1];
        p0[2] = xyz[2];

        Point3F p1;
        p1[0] = p0[0] + height*axis[0];
        p1[1] = p0[1] + height*axis[1];
        p1[2] = p0[2] + height*axis[2];

        Color  color;
        obj->getAttribute("Color", color);

        gleNumSides( numSlices);

        drawCylinder( &p0[0], &p1[0], radius);
    */
}


void
JShapeViewer::fastDraw()
{
}

///////////////////////////////////////////////////////////////////////////////
void
JShapeViewer::draw()
{
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glPushMatrix();

    /*
        for( const auto obj : spheres)   drawObject(obj);
        for( const auto obj : cylinders) drawObject(obj);
        for( const auto obj : cones)     drawObject(obj);
    */

    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////

