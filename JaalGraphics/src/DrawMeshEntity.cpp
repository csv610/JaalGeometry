
#include "MeshViewer.hpp"
#include "JaalViewer.hpp"
#include "Lights.hpp"

using namespace Jaal;
using namespace std;

GLUquadricObj* JNodeDraw :: sphereObj = gluNewQuadric();
GLUquadricObj* JNodeDraw :: diskObj   = gluNewQuadric();

JNodeRenderPtr JNodeRender :: newObject() {
    return  boost::shared_ptr<JNodeRender>( new JNodeRender );
}

JEdgeRenderPtr JEdgeRender :: newObject() {
    return  boost::shared_ptr<JEdgeRender>( new JEdgeRender );
}

JFaceRenderPtr JFaceRender :: newObject() {
    return  boost::shared_ptr<JFaceRender>( new JFaceRender );
}

JCellRenderPtr JCellRender :: newObject() {
    return  boost::shared_ptr<JCellRender>( new JCellRender );
}

int JMeshRender :: setNodeScalarFieldColor( const JMeshPtr &mesh, const vector<double> &darray)
{
    double minVal = *boost::min_element(darray );
    double maxVal = *boost::max_element(darray );

    size_t numnodes = mesh->getSize(0);

    JColor clr;
    size_t index = 0;
    double  r, g, b, val, p;
    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( !vertex->isActive() )  continue;
        val = darray[index++];
        p = (val-minVal)/(maxVal-minVal);
        if( p < 0.5) {
            r = 1.0 - 2.0*p;
            g = 2.0*p;
            b = 0.0;
        } else {
            p = p - 0.5;
            r = 0.0;
            g = 1.0 - 2.0*p;
            b = 2.0*p;
        }
        clr[0] = r;
        clr[1] = g;
        clr[2] = b;
        clr[3] = 1.0;
        vertex->getAttribute("Render", nAttrib);
        nAttrib->color = clr;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JMeshRender :: setNodeScalarFieldColor( const JMeshPtr &mesh, const string &attribname)
{
    if( mesh == nullptr ) return 1;

    size_t numnodes = mesh->getSize(0);
    vector<double> scalar(numnodes);
    int index = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            vertex->getAttribute(attribname, scalar[index++] );
        }
    }
    JMeshRender :: setNodeScalarFieldColor( mesh, scalar);
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

void JMeshEntityRender :: drawTextAt( char *text, const float *xyz) const
{
    glPushMatrix();
    glTranslatef(xyz[0], xyz[1], xyz[2] );
    glScalef(fontScale, fontScale, fontScale);
    font->Render(text);
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////

JNodeDraw :: JNodeDraw()
{
    offset        = 0;
    offset_factor = 1.0;
    offset_unit   = 1.0;
    antiAlias    = 0;
    alpha        = 0.8;
    glyph        =  JNodeRender::NODE_AS_POINT;

    normalLength = 1.0;
    normalsColor[0] = 0.5;
    normalsColor[1] = 0.5;
    normalsColor[2] = 0.5;
    normalsColor[3] = 1.0;

    color = JEntityColor::getColor("Pink");

    pointSize =  1.0;
    ballRadius = 0.1;
    numSlices =  16;
    numStacks =  16 ;

    gluQuadricDrawStyle( sphereObj, GLU_FILL);
    gluQuadricNormals( sphereObj,   GLU_SMOOTH);
    gluQuadricDrawStyle( diskObj,   GLU_FILL);
}

///////////////////////////////////////////////////////////////////////////////

JNodeRenderPtr JNodeDraw :: initRenderAttribute(const JNodePtr &vtx)
{
    JNodeRenderPtr attrib;
    int err = vtx->getAttribute("Render", attrib);
    if( err ) {
        attrib.reset(new JNodeRender);
        vtx->setAttribute("Render", attrib);
    }

    assert(attrib);

    attrib->glyph    = getGlyph();
    attrib->scale    = 1.0;
    attrib->display  = 1;
    attrib->color       = color;
    attrib->ballRadius  = getBallRadius();
    attrib->pointSize   = getPointSize();
    return attrib;
}

///////////////////////////////////////////////////////////////////////////////
void JNodeDraw :: initRenderAttributes(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() && !vtx->hasAttribute("Render")  )
            initRenderAttribute(vtx);
    }
}
///////////////////////////////////////////////////////////////////////////////

void JNodeDraw :: updateBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;
    mesh->enumerate(0);

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JNodeBuffers> vbo = mrender->nodeBuffers;

    vbo->clearAll();

    size_t numNodes = mesh->getSize(0);
    if( numNodes == 0) return;

    JNodeRenderPtr nattrib;
    Vec3F normal;
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            int err = vtx->getAttribute("Render", nattrib);
            if( err ) nattrib = initRenderAttribute(vtx);
            nattrib->color[3] = alpha;
            if( nattrib->display) {
                const Point3D &xyz = vtx->getXYZCoords();
                vbo->rootGroup.push_back(vtx->getID() );
                if( nattrib->glyph == JNodeRender::NODE_AS_SPHERE) {
                    double radius  = nattrib->ballRadius*nattrib->scale;
                    vbo->sphereGroup[radius].push_back( vtx->getID() );
                }
                if( nattrib->glyph == JNodeRender::NODE_AS_POINT) {
                    double pointsize = nattrib->pointSize*nattrib->scale;
                    vbo->pointGroup[pointsize].push_back( vtx->getID() );
                }
                if( nattrib->glyph == JNodeRender::NODE_AS_SPLAT) {
                    double radius  = nattrib->ballRadius*nattrib->scale;
                    vbo->splatGroup[radius].push_back(vtx->getID() );
                }
                vbo->coordsArray.push_back(xyz[0] );
                vbo->coordsArray.push_back(xyz[1] );
                vbo->coordsArray.push_back(xyz[2] );
                vbo->colorArray.push_back( nattrib->color[0] );
                vbo->colorArray.push_back( nattrib->color[1] );
                vbo->colorArray.push_back( nattrib->color[2] );
                vbo->colorArray.push_back( nattrib->color[3] );
                err = vtx->getAttribute("Normal", normal);
                if( !err) {
                    vbo->normalArray.push_back( normal[0] );
                    vbo->normalArray.push_back( normal[1] );
                    vbo->normalArray.push_back( normal[2] );
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JNodeDraw :: preRender()
{
    glEnable( GL_POLYGON_OFFSET_POINT);
    glEnable( GL_POLYGON_OFFSET_LINE);

    if( glyph == JNodeRender::NODE_AS_POINT) {
        glDisable(GL_LIGHTING);
        glEnable( GL_DEPTH_TEST );
        if( antiAlias) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable (GL_BLEND);
            glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glHint( GL_POINT_SMOOTH_HINT, GL_NICEST);
            glEnable( GL_POINT_SMOOTH );
        } else {
            glDisable (GL_BLEND);
            glDisable( GL_POINT_SMOOTH );
        }
        glDisable(GL_COLOR_MATERIAL);
    }

    if( glyph == JNodeRender::NODE_AS_SPHERE) {
        glEnable( GL_LIGHTING );
        glEnable( GL_LIGHT0);
        glEnable( GL_DEPTH_TEST );
        glEnable( GL_CULL_FACE );
        glMaterialfv(GL_FRONT, GL_SPECULAR, White);
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT, GL_SHININESS, 100);
        glEnable(GL_COLOR_MATERIAL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
}
///////////////////////////////////////////////////////////////////////////////
void JNodeDraw ::  postRender()
{
    glDisable(GL_BLEND );
    glDisable(GL_POINT_SMOOTH );
    glDisable( GL_COLOR_MATERIAL );
    glDisable( GL_POLYGON_OFFSET_POINT);
}

///////////////////////////////////////////////////////////////////////////////
void
JNodeDraw::draw(const JNodePtr &vtx)
{
    JNodeRenderPtr attrib;
    vtx->getAttribute("Render", attrib);
    const Point3D xyz = vtx->getXYZCoords();

    if( attrib->glyph == JNodeRender::NODE_AS_POINT) {
        double pointSize = (attrib->pointSize)*(attrib->scale);
        glPointSize(pointSize);
        glColor3fv(  &attrib->color[0] );
        glBegin(GL_POINTS);
        glVertex3dv( &xyz[0] );
        glEnd();
    }

    if( glyph == JNodeRender::NODE_AS_SPHERE) {
        glEnable(GL_LIGHTING);
        glEnable(GL_COLOR_MATERIAL);
        double ballRadius = (attrib->ballRadius)*(attrib->scale);
        glColor3fv(  &attrib->color[0] );
        glPushMatrix();
        glTranslatef(xyz[0], xyz[1], xyz[2]);
        gluSphere( sphereObj, ballRadius, numSlices, numStacks);
        glPopMatrix();
    }

    /*
    double innerRadius = ballRadius*(1.0-borderThickness);
            if( glyph == NODE_AS_CIRCLE || glyph == NODE_AS_SQUARE) {
                glDisable(GL_LIGHTING);
                glEnable( GL_CULL_FACE);
                Vec3F srcVec, dstVec, perpAxis;
                glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
                glMaterialf(GL_FRONT, GL_SHININESS, 100);
                GLdouble mat[16];
                for( size_t i = 0; i < numNodes; i++) {
                    srcVec[0]  = vbo->normalHeadArray[3*i];
                    srcVec[1]  = vbo->normalHeadArray[3*i+1];
                    srcVec[0]  = vbo->normalHeadArray[3*i+2];
                    srcVec[0] = 0.0;
                    srcVec[1] = 0.0;
                    srcVec[2] = 1.0;
                    JMath::cross_product( dstVec, srcVec, perpAxis);
                    double angle = JMath::getVecAngle(srcVec, dstVec, ANGLE_IN_RADIANS);
                    qglviewer::Vec rotaxis(perpAxis[0], perpAxis[1], perpAxis[2] );
                    qglviewer::Quaternion quaternion(rotaxis, -1.0*angle);
                    glNormal3f( dstVec[0], dstVec[1], dstVec[2] );
                    glPushMatrix();
                    float x = vbo->coordsArray[3*i];
                    float y = vbo->coordsArray[3*i+1];
                    float z = vbo->coordsArray[3*i+2];
                    glTranslatef( x, y, z);
                    quaternion.getMatrix(mat);
                    glMultMatrixd( mat );

                    glDisable (GL_BLEND);
                    glColor3fv( &borderColor[0] );
                    if( glyph == NODE_AS_SQUARE)
                        gluDisk( diskObj, innerRadius, ballRadius, 4, 1);
                    else
                        gluDisk( diskObj, innerRadius, ballRadius, numSlices, 1);

                    glEnable (GL_BLEND);
                    glColor4fv(&vbo->colorArray[4*i] );
                    if( glyph == NODE_AS_SQUARE)
                        gluDisk( diskObj, 0.0, innerRadius, 4, 1);
                    else
                        gluDisk( diskObj, 0.0, innerRadius, numSlices, 1);


                    glPopMatrix();
                }
            }
    */
}
///////////////////////////////////////////////////////////////////////////////

void
JNodeDraw::draw(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);

    boost::shared_ptr<JNodeBuffers> vbo = mrender->nodeBuffers;

    glPushMatrix();

    preRender();
    glDisable( GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    for( auto &keyval : vbo->pointGroup) {
        double pointSize = keyval.first;
        const vector<int> &group = keyval.second;
        glPointSize(pointSize);
        size_t numNodes = group.size();
        glBegin(GL_POINTS);
        for( size_t i = 0; i < numNodes; i++) {
            int id =  group[i];
            glColor3fv(  &vbo->colorArray[4*id]);
            glVertex3fv( &vbo->coordsArray[3*id]);
        }
        glEnd();
    }

    for( auto &keyval : vbo->sphereGroup) {
        glEnable( GL_LIGHTING);
        glEnable( GL_LIGHT0);
        glEnable(GL_COLOR_MATERIAL);
        double radius = keyval.first;
        const vector<int> &group  = keyval.second;
        size_t numNodes = group.size();
        for( size_t i = 0; i < numNodes; i++) {
            int id = group[i];
            glColor4fv(  &vbo->colorArray[4*id]);
            glPushMatrix();
            float x = vbo->coordsArray[3*id];
            float y = vbo->coordsArray[3*id+1];
            float z = vbo->coordsArray[3*id+2];
            glTranslatef(x, y, z);
            gluSphere( sphereObj, radius, numSlices, numStacks);
            glPopMatrix();
        }
    }

    Vec3F srcVec, dstVec, perpAxis;
    GLdouble mat[16];
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glMaterialf(GL_FRONT, GL_SHININESS, 100);
    for( auto &keyval : vbo->splatGroup) {
        double radius = keyval.first;
        double innerRadius = radius*(1.0-borderThickness);
        const vector<int> &group  = keyval.second;
        size_t numNodes = group.size();
        for( size_t i = 0; i < numNodes; i++) {
            int id = group[i];
            dstVec[0]  = vbo->normalArray[3*id];
            dstVec[1]  = vbo->normalArray[3*id+1];
            dstVec[2]  = vbo->normalArray[3*id+2];
            srcVec[0] = 0.0;
            srcVec[1] = 0.0;
            srcVec[2] = 1.0;
            JMath::cross_product( dstVec, srcVec, perpAxis);
            double angle = JMath::getVecAngle(srcVec, dstVec, ANGLE_IN_RADIANS);
            qglviewer::Vec rotaxis(perpAxis[0], perpAxis[1], perpAxis[2] );
            qglviewer::Quaternion quaternion(rotaxis, -1.0*angle);
            glNormal3f( dstVec[0], dstVec[1], dstVec[2] );
            glPushMatrix();
            float x = vbo->coordsArray[3*id+0] + 0.001*srcVec[0];
            float y = vbo->coordsArray[3*id+1] + 0.001*srcVec[1];
            float z = vbo->coordsArray[3*id+2] + 0.001*srcVec[2];
            glTranslatef( x, y, z);
            quaternion.getMatrix(mat);
            glMultMatrixd( mat );
            glDisable (GL_BLEND);
            glColor3fv( &borderColor[0] );
            gluDisk( diskObj, innerRadius, radius, numSlices, 1);
            glEnable (GL_BLEND);
            glColor4fv(&vbo->colorArray[4*id] );
            gluDisk( diskObj, 0.0, innerRadius, numSlices, 1);
            glPopMatrix();
        }
    }
    postRender();
    glPopMatrix();

    drawNormals( mesh );
    glPushMatrix();
    drawIDs( mesh );
    glPopMatrix();
}
////////////////////////////////////////////////////////////////////////////////
void
JNodeDraw::drawNormals(const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->displayNormals[0] == 0) return;

    boost::shared_ptr<JNodeBuffers> vbo = mrender->nodeBuffers;

    glDisable( GL_LINE_SMOOTH);
    glColor3fv( &normalsColor[0] );
    glLineWidth(1.0);

    glDisable( GL_LIGHTING );

    glBegin(GL_LINES);
    for( auto id : vbo->rootGroup) {
        float x = vbo->coordsArray[3*id+0];
        float y = vbo->coordsArray[3*id+1];
        float z = vbo->coordsArray[3*id+2];
        glVertex3f( x, y, z);
        x += normalLength*vbo->normalArray[3*id+0];
        y += normalLength*vbo->normalArray[3*id+1];
        z += normalLength*vbo->normalArray[3*id+2];
        glVertex3f( x, y, z);
    }
    glEnd();
}

////////////////////////////////////////////////////////////////////////////////
void
JNodeDraw::withName(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->pickableEntity != 0) return;

    boost::shared_ptr<JNodeBuffers> vbo = mrender->nodeBuffers;

    preRender();

    glPushMatrix();
    for( auto &keyval : vbo->pointGroup) {
        glDisable( GL_LIGHTING);
        double pointSize = keyval.first;
        const vector<int> &group = keyval.second;
        glPointSize(pointSize);
        size_t numNodes = group.size();
        for( size_t i = 0; i < numNodes; i++) {
            size_t id =  group[i];
            glPushName(id);
            glBegin(GL_POINTS);
            glColor4fv(  &vbo->colorArray[4*id]);
            glVertex3fv( &vbo->coordsArray[3*id]);
            glEnd();
            glPopName();
        }
    }

    for( auto &keyval : vbo->sphereGroup) {
        glEnable( GL_LIGHTING);
        double radius = keyval.first;
        const vector<int> &group  = keyval.second;
        size_t  numNodes = group.size();
        for( size_t i = 0; i < numNodes; i++) {
            int id = group[i];
            glPushName(id);
            glPushMatrix();
            float x = vbo->coordsArray[3*id];
            float y = vbo->coordsArray[3*id+1];
            float z = vbo->coordsArray[3*id+2];
            glTranslatef(x, y, z);
            glColor4fv(  &vbo->colorArray[4*id]);
            gluSphere( sphereObj, radius, numSlices, numStacks);
            glPopMatrix();
            glPopName();
        }
    }
    glPopMatrix();

    double innerRadius = ballRadius*(1.0-borderThickness);

    Vec3F srcVec, dstVec, perpAxis;
    GLdouble mat[16];
    glDisable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glMaterialf(GL_FRONT, GL_SHININESS, 100);
    for( auto &keyval : vbo->splatGroup) {
//      double radius = keyval.first;
        const vector<int> &group  = keyval.second;
        size_t numNodes = group.size();
        for( size_t i = 0; i < numNodes; i++) {
            int id = group[i];
            glPushName(id);
            dstVec[0]  = vbo->normalArray[3*id];
            dstVec[1]  = vbo->normalArray[3*id+1];
            dstVec[2]  = vbo->normalArray[3*id+2];
            srcVec[0] = 0.0;
            srcVec[1] = 0.0;
            srcVec[2] = 1.0;
            JMath::cross_product( dstVec, srcVec, perpAxis);
            double angle = JMath::getVecAngle(srcVec, dstVec, ANGLE_IN_RADIANS);
            qglviewer::Vec rotaxis(perpAxis[0], perpAxis[1], perpAxis[2] );
            qglviewer::Quaternion quaternion(rotaxis, -1.0*angle);
            glNormal3f( dstVec[0], dstVec[1], dstVec[2] );
            glPushMatrix();
            float x = vbo->coordsArray[3*id+0] + 0.0001*srcVec[0];
            float y = vbo->coordsArray[3*id+1] + 0.0001*srcVec[1];
            float z = vbo->coordsArray[3*id+2] + 0.0001*srcVec[2];
            glTranslatef( x, y, z);
            quaternion.getMatrix(mat);
            glMultMatrixd( mat );
            glDisable (GL_BLEND);
            glColor4fv( &borderColor[0] );
            gluDisk( diskObj, innerRadius, ballRadius, numSlices, 1);
            glEnable (GL_BLEND);
            glColor4fv(&vbo->colorArray[4*id] );
            gluDisk( diskObj, 0.0, innerRadius, numSlices, 1);
            glPopMatrix();
            glPopName();
        }
    }
    postRender();
}
////////////////////////////////////////////////////////////////////////////////

void
JNodeDraw::drawIDs( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->displayIDs[0] == 0) return;
    boost::shared_ptr<JNodeBuffers> vbo = mrender->nodeBuffers;

    glDisable( GL_LIGHTING );
    glDisable (GL_BLEND);
    glColor3f(0.1, 0.1, 0.1);

    int height = viewManager->height();

    glEnable( GL_DEPTH_TEST);

    glPushMatrix();

    QFont newFont("Times", 12);
    Vec pvec, normal;
    viewManager->startScreenCoordinatesSystem();

    /*
       Vec camdir = -viewManager->camera()->viewDirection();
        for( size_t id  : vbo->rootGroup) {
            normal[0] = vbo->normalArray[3*id+0];
            normal[1] = vbo->normalArray[3*id+1];
            normal[2] = vbo->normalArray[3*id+2];
            if( normal[0]*camdir[0] + normal[1]*camdir[1] + normal[2]*camdir[2] >= 0.1)  {
                pvec[0]  = vbo->coordsArray[3*id+0];
                pvec[1]  = vbo->coordsArray[3*id+1];
                pvec[2]  = vbo->coordsArray[3*id+2];
                pvec     = viewManager->camera()->projectedCoordinatesOf(pvec);
                viewManager->drawText( pvec[0], pvec[1], QString::number(id), newFont );
            }
        }
    */


    for( size_t id  : vbo->rootGroup) {
        pvec[0]  = vbo->coordsArray[3*id+0];
        pvec[1]  = vbo->coordsArray[3*id+1];
        pvec[2]  = vbo->coordsArray[3*id+2];
        pvec     = viewManager->camera()->projectedCoordinatesOf(pvec);
        viewManager->drawText( pvec[0], pvec[1], QString::number(id), newFont );
    }


    glPopMatrix();
}
////////////////////////////////////////////////////////////////////////////////

JEdgeDraw :: JEdgeDraw()
{
    offset        = 0;
    offset_factor = 1.0;
    offset_unit   = 2.0;
    antiAlias     = 0;
    cylRadius     = 0.1;
    lineWidth     = 1.0;
    withSteiner   = 0;
    numCylSides = 20;
    glyph  = EDGE_AS_LINE;

    gleSetJoinStyle (TUBE_NORM_EDGE | TUBE_JN_ANGLE | TUBE_JN_CAP);
    gleSetJoinStyle( TUBE_JN_ROUND );
}

////////////////////////////////////////////////////////////////////////////////

JEdgeRenderPtr JEdgeDraw :: initRenderAttribute(const JEdgePtr &edge)
{
    JEdgeRenderPtr attrib;
    int err = edge->getAttribute("Render", attrib);
    if( err ) {
        attrib.reset(new JEdgeRender);
        assert(attrib);
        edge->setAttribute("Render", attrib);
    }
    attrib->display = 1;

    if( edge->isBoundary() ) {
        attrib->glyph   = 0;
        attrib->scale   = 1;
        color = JEntityColor::getColor("Red");
    } else {
        attrib->glyph   = 0;
        attrib->scale   = 1;
        color[0]        = 0.3;
        color[1]        = 0.3;
        color[2]        = 0.3;
        color[3]        = 1.0;
    }

    attrib->lineWidth  = getLineWidth();
    attrib->cylinderRadius   = getCylinderRadius();
    attrib->color       = color;
    return attrib;
}

////////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: initRenderAttributes(const JMeshPtr &mesh)
{
    size_t numEdges = mesh->getSize(1);
    JEdgeRenderPtr  attrib;
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() && !edge->hasAttribute("Render")  ) {
            initRenderAttribute(edge);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: activateLowerEntities( const JEdgePtr &edge)
{
    if( edge == nullptr ) return;

    JEdgeRenderPtr  eattrib;
    int err = edge->getAttribute("Render", eattrib);
    if( err) return;

    JNodeRenderPtr   vattrib;
    for( int i = 0; i < 2; i++) {
        const JNodePtr &vtx = edge->getNodeAt(i);
        int err = vtx->getAttribute("Render", vattrib);
        if( !err) vattrib->display = eattrib->display;
    }
}
////////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: updateGeometryBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JEdgeBuffers> ebo =  mrender->edgeBuffers;

    ebo->coordsArray.clear();
    size_t numNodes = mesh->getSize(0);
    ebo->coordsArray.reserve(3*numNodes);

    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            const Point3D &p3d = vtx->getXYZCoords();
            ebo->coordsArray.push_back(p3d[0]);
            ebo->coordsArray.push_back(p3d[1]);
            ebo->coordsArray.push_back(p3d[2]);
        }
    }

    /*
        size_t numEdges  = mesh->getSize(1);
        ebo->normalArray.reserve(3*numEdges);

        Vec3F n0, n1;
        for( size_t i = 0; i < numEdges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                const JNodePtr &v0 = edge->getNodeAt(0);
                err  = v0->getAttribute("Normal", n0);
                assert( !err);
                const JNodePtr &v1 = edge->getNodeAt(1);
                err  = v1->getAttribute("Normal", n1);
                assert( !err);
                float nx = 0.5*( n0[0] + n1[0]);
                float ny = 0.5*( n0[1] + n1[1]);
                float nz = 0.5*( n0[2] + n1[2]);
                ebo->normalArray.push_back(nx);
                ebo->normalArray.push_back(ny);
                ebo->normalArray.push_back(nz);
            }
        }
    */
}
////////////////////////////////////////////////////////////////////////////////

void JEdgeDraw :: updateTopologyBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JEdgeBuffers> ebo =  mrender->edgeBuffers;

    ebo->lineGroup.clear();
    ebo->cylGroup.clear();
    ebo->nodesArray.clear();
    ebo->colorArray.clear();
    ebo->edgePos.clear();

    JEdgeRenderPtr attrib;

    size_t numEdges = mesh->getSize(1);
    ebo->colorArray.reserve(4*numEdges);
    ebo->nodesArray.reserve(2*numEdges);

    JNodePtr vertex;
    double val;
    size_t pos = 0;
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            int err = edge->getAttribute("Render", attrib);
            if( err ) attrib = initRenderAttribute(edge);
            if( attrib->display) {
                ebo->edgePos[edge->getID()] = pos++;
                if( attrib->glyph == EDGE_AS_LINE) {
                    val = attrib->lineWidth*attrib->scale;
                    ebo->lineGroup[val].push_back(edge->getID());
                }
                if( attrib->glyph == EDGE_AS_CYLINDER) {
                    val  = attrib->cylinderRadius*attrib->scale;
                    ebo->cylGroup[val].push_back(edge->getID());
                }
                vertex = edge->getNodeAt(0);
                int v1 = vertex->getID();
                vertex = edge->getNodeAt(1);
                int v2 = vertex->getID();
                ebo->nodesArray.push_back(v1);
                ebo->nodesArray.push_back(v2);
                ebo->colorArray.push_back( attrib->color[0] );
                ebo->colorArray.push_back( attrib->color[1] );
                ebo->colorArray.push_back( attrib->color[2] );
                ebo->colorArray.push_back( attrib->color[3] );
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: updateBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;
    mesh->enumerate(1);
    updateGeometryBuffers(mesh);
    updateTopologyBuffers(mesh);
}
////////////////////////////////////////////////////////////////////////////////

void JEdgeDraw :: preRender()
{
// glEnable( GL_POLYGON_OFFSET_POINT);
    glEnable( GL_POLYGON_OFFSET_LINE);
// glPolygonOffset( offset_factor, offset_unit);
    glPolygonOffset( 1.9, 1.0);

    glDisable (GL_BLEND);
    glDisable( GL_LINE_SMOOTH );

    if( antiAlias ) {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glHint( GL_LINE_SMOOTH_HINT, GL_NICEST);
        glEnable( GL_LINE_SMOOTH );
    }

    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

}
////////////////////////////////////////////////////////////////////////////////

void JEdgeDraw :: postRender()
{
    glDisable( GL_BLEND );
    glDisable( GL_COLOR_MATERIAL );
    glDisable( GL_LINE_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable( GL_POLYGON_OFFSET_LINE);
}
////////////////////////////////////////////////////////////////////////////////

void JEdgeDraw :: drawCylinder( const float *p0, const float *p1, float radius) const
{
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
}
///////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: drawCylinder( const JEdgePtr &edge, double radius ) const
{
    if( !edge->isActive() ) return;

    Point3D xyz;

    glPushMatrix();
//  gleDouble endPoints[4][3];

    xyz = edge->getNodeAt(0)->getXYZCoords();
    float p0[3];
    p0[0] = xyz[0];
    p0[1] = xyz[1];
    p0[2] = xyz[2];

    xyz = edge->getNodeAt(1)->getXYZCoords();
    float p1[3];
    p1[0] = xyz[0];
    p1[1] = xyz[1];
    p1[2] = xyz[2];
    drawCylinder( p0, p1, (float)radius);
}


///////////////////////////////////////////////////////////////////////////////

void
JEdgeDraw::draw(const JEdgePtr &edge)
{
    JEdgeRenderPtr attrib;
    int err = edge->getAttribute("Render", attrib);
    if( err) return;

    JNodePtr vtx;
    if( attrib->glyph == JEdgeDraw::EDGE_AS_LINE) {
        double lineWidth = attrib->lineWidth*attrib->scale;
        glDisable( GL_LIGHTING);
        glLineWidth(lineWidth);
        glColor3fv( &attrib->color[0] );
        glBegin(GL_LINES);
        for( int i = 0; i < 2; i++) {
            vtx = edge->getNodeAt(i);
            const Point3D xyz = vtx->getXYZCoords();
            glVertex3dv( &xyz[0] );
        }
        glEnd();
    }

    if( attrib->glyph == JEdgeDraw::EDGE_AS_CYLINDER) {
        glEnable( GL_LIGHTING );
        glEnable( GL_LIGHT0);
        glEnable( GL_CULL_FACE );
        glMaterialfv(GL_FRONT, GL_SPECULAR, White);
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT, GL_SHININESS, 100);
        glEnable(GL_COLOR_MATERIAL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3fv( &attrib->color[0] );
        double cylRadius = attrib->cylinderRadius*attrib->scale;
        drawCylinder(edge, cylRadius);
    }
}
///////////////////////////////////////////////////////////////////////////////

void
JEdgeDraw::draw(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JEdgeBuffers> ebo =  mrender->edgeBuffers;

    size_t numEdges = ebo->nodesArray.size()/2;
    if( numEdges == 0) return;

    preRender();

    for( auto &keyval : ebo->lineGroup) {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDisable( GL_LIGHTING);
        double linewidth = keyval.first;
        glLineWidth(linewidth);
        const vector<int> &group = keyval.second;
        int numlines = group.size();
        glBegin(GL_LINES);
        for( int i = 0; i < numlines; i++) {
            int eid = group[i];
            int pos = ebo->edgePos[eid];
            int v0  = ebo->nodesArray[2*pos];
            int v1  = ebo->nodesArray[2*pos+1];
            glColor4fv(  &ebo->colorArray[4*eid] );
            glVertex3fv( &ebo->coordsArray[3*v0]);
            glVertex3fv( &ebo->coordsArray[3*v1]);
        }
        glEnd();
        glDisable (GL_BLEND);
    }

    for( auto &keyval : ebo->cylGroup) {
        glDisable( GL_LIGHTING);
        double cylRadius = keyval.first;
        const vector<int> &group = keyval.second;
        int numCylinders = group.size();
        glEnable( GL_LIGHTING );
        glEnable( GL_LIGHT0);
        glEnable( GL_CULL_FACE );
        glMaterialfv(GL_FRONT, GL_SPECULAR, White);
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT, GL_SHININESS, 100);
        glEnable(GL_COLOR_MATERIAL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        for( int i = 0; i < numCylinders; i++) {
            int eid = group[i];
            int pos = ebo->edgePos[eid];
            int v0 = ebo->nodesArray[2*pos];
            int v1 = ebo->nodesArray[2*pos+1];
            glColor3fv( &ebo->colorArray[4*eid] );
            drawCylinder( &ebo->coordsArray[3*v0], &ebo->coordsArray[3*v1], cylRadius);
        }
    }
    postRender();

    drawIDs( mesh );
}
////////////////////////////////////////////////////////////////////////////////
void
JEdgeDraw::withName(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->pickableEntity != 1 ) return;

    boost::shared_ptr<JEdgeBuffers> ebo =  mrender->edgeBuffers;

    size_t numEdges = ebo->nodesArray.size()/2;
    if( numEdges == 0) return;

    glPushMatrix();
    preRender();
    for( auto &keyval : ebo->lineGroup) {
        glDisable( GL_LIGHTING);
        double linewidth = keyval.first;
        glLineWidth(linewidth);
        const vector<int> &group = keyval.second;
        int numlines = group.size();
        for( int i = 0; i < numlines; i++) {
            int eid = group[i];
            int pos = ebo->edgePos[eid];
            glPushName( eid );
            glBegin(GL_LINES);
            int v0 = ebo->nodesArray[2*pos];
            int v1 = ebo->nodesArray[2*pos+1];
            glColor3fv( &ebo->colorArray[4*eid] );
            glVertex3fv( &ebo->coordsArray[3*v0]);
            glVertex3fv( &ebo->coordsArray[3*v1]);
            glEnd();
            glPopName();
        }
    }

    for( auto &keyval : ebo->cylGroup) {
        glDisable( GL_LIGHTING);
        double cylRadius = keyval.first;
        const vector<int> &group = keyval.second;
        int numCylinders = group.size();
        glEnable( GL_LIGHTING );
        glEnable( GL_LIGHT0);
        glEnable( GL_CULL_FACE );
        glMaterialfv(GL_FRONT, GL_SPECULAR, White);
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT, GL_SHININESS, 100);
        glEnable(GL_COLOR_MATERIAL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        for( int i = 0; i < numCylinders; i++) {
            int  eid = group[i];
            int  pos = ebo->edgePos[eid];
            glPushName( eid );
            int v0 = ebo->nodesArray[2*pos];
            int v1 = ebo->nodesArray[2*pos+1];
            glColor3fv( &ebo->colorArray[4*eid] );
            drawCylinder( &ebo->coordsArray[3*v0], &ebo->coordsArray[3*v1], cylRadius);
            glPopName();
        }
    }

    postRender();
    glPopMatrix();
}
////////////////////////////////////////////////////////////////////////////////

void
JEdgeDraw::drawIDs( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->displayIDs[1] == 0) return;

    glDisable( GL_LIGHTING );
    glDisable (GL_BLEND);

    glColor3f(0.0, 0.0, 0.0);

    //  Vec3D  srcVec, dstVec, perpAxis;
    //  GLdouble mat[16];

    Vec camdir = -viewManager->camera()->viewDirection();

    glEnable( GL_DEPTH_TEST);

    //   char number[128];
    //    float pos[3];
    glPushMatrix();

    QFont newFont("Times", 20);

    Vec    pvec;
    Vec3F normal;

    /*
        viewManager->startScreenCoordinatesSystem();
        for( const auto &keyVal : ebo->edgePos) {
        int  id   = keyVal.first;
        int  vpos = keyVal.second;

                normal[0] = ebo->normalArray[3*id+0];
                normal[1] = ebo->normalArray[3*id+1];
                normal[2] = ebo->normalArray[3*id+2];
                if( normal[0]*camdir[0] + normal[1]*camdir[1] + normal[2]*camdir[2] >= 0.0)  {
                    pvec[0]  = 0.0;
                    pvec[1]  = 0.0;
                    pvec[2]  = 0.0;

                    vid      =  ebo->nodesArray[2*vpos];
                    pvec[0]  += ebo->coordsArray[3*vid+0];
                    pvec[1]  += ebo->coordsArray[3*vid+1];
                    pvec[2]  += ebo->coordsArray[3*vid+2];

                    vid      =  ebo->nodesArray[2*vpos+1];
                    pvec[0]  += ebo->coordsArray[3*vid+0];
                    pvec[1]  += ebo->coordsArray[3*vid+1];
                    pvec[2]  += ebo->coordsArray[3*vid+2];

                    pvec[0]  *= 0.5;
                    pvec[1]  *= 0.5;
                    pvec[2]  *= 0.5;
                    pvec     = viewManager->camera()->projectedCoordinatesOf(pvec);
                    viewManager->drawText( pvec[0], pvec[1], QString::number(id), newFont );
                }
        }
        glPopMatrix();
    */
}

////////////////////////////////////////////////////////////////////////////////

JFaceDraw :: JFaceDraw()
{
    offset        = 1;
    offset_factor = 1.0;
    offset_unit   = 1.0;
    culling       = 0;
    normalsColor[0] = 0.6;
    normalsColor[1] = 0.6;
    normalsColor[2] = 0.6;
    use_normal  = 0;
    screendoorPattern = 8;
    color          = JEntityColor::getColor("SnowWhite");
    color[0]  = 213.0/255.0;
    color[1]  = 229.0/255.0;
    color[2]  = 1.0;
    color[3]  = 0.5;
    lights         = 1;
    frontface      = GL_CCW;
}

////////////////////////////////////////////////////////////////////////////////
JFaceRenderPtr JFaceDraw :: initRenderAttribute(const JFacePtr &face)
{
    JFaceRenderPtr attrib;
    int err = face->getAttribute("Render", attrib);
    if( err ) {
        attrib.reset(new JFaceRender);
        face->setAttribute("Render", attrib);
    }
    assert(attrib);
    attrib->display = 1;
    attrib->color   = color;
    return attrib;
}
////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: initRenderAttributes(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    size_t numFaces = mesh->getSize(2);
    JFaceRenderPtr  attrib;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() && !face->hasAttribute("Render")  )
            initRenderAttribute(face);
    }
}
////////////////////////////////////////////////////////////////////////////////
void JFaceDraw :: activateLowerEntities( const JFacePtr &face)
{
    if( face == nullptr ) return;

    JFaceRenderPtr   fattrib;
    int err = face->getAttribute("Render", fattrib);
    if( err ) return;

    JEdgeRenderPtr   eattrib;
    for( int i = 0; i < face->getSize(1); i++) {
        const JEdgePtr &edge = face->getEdgeAt(i);
        int err = edge->getAttribute("Render", eattrib);
        if( !err) eattrib->display = fattrib->display;
    }

    JNodeRenderPtr   vattrib;
    for( int i = 0; i < face->getSize(0); i++) {
        const JNodePtr &vtx = face->getNodeAt(i);
        int err = vtx->getAttribute("Render", vattrib);
        if( !err) vattrib->display = fattrib->display;
    }
}

////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: updateGeometryBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    fbo->nodeCoords.clear();
    fbo->nodeNormal.clear();
    fbo->nodeColor.clear();

    Vec3F normal;
    JNodeRenderPtr vattrib;

    size_t numNodes = mesh->getSize(0);
    fbo->nodeCoords.reserve(3*numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Render", vattrib);
        if( vtx->isActive() ) {
            const Point3D &p3d = vtx->getXYZCoords();
            fbo->nodeCoords.push_back(p3d[0]);
            fbo->nodeCoords.push_back(p3d[1]);
            fbo->nodeCoords.push_back(p3d[2]);
            vtx->getAttribute("Normal", normal);
            fbo->nodeNormal.push_back(normal[0]);
            fbo->nodeNormal.push_back(normal[1]);
            fbo->nodeNormal.push_back(normal[2]);
            fbo->nodeColor.push_back(vattrib->color[0]);
            fbo->nodeColor.push_back(vattrib->color[1]);
            fbo->nodeColor.push_back(vattrib->color[2]);
            fbo->nodeColor.push_back(vattrib->color[3]);
        }
    }

    size_t numFaces = mesh->getSize(2);
    fbo->faceNormal.clear();
    fbo->faceNormalHead.clear();
    fbo->faceNormalTail.clear();
    JFaceRenderPtr fattrib;

    float x,y,z;
    Point3D centroid;

    size_t nCount = 0;
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int err = face->getAttribute("Render", fattrib);
            if( err ) fattrib = initRenderAttribute(face);
            if( fattrib->display) {
                err = face->getAttribute("Normal", normal);
                if( err )  {
                    mesh->getGeometry()->setFacesNormal();
                    face->getAttribute("Normal", normal);
                }
                fbo->faceNormal.push_back(normal[0]);
                fbo->faceNormal.push_back(normal[1]);
                fbo->faceNormal.push_back(normal[2]);
                JFaceGeometry::getCentroid(face,centroid);
                x = centroid[0];
                y = centroid[1];
                z = centroid[2];
                fbo->faceNormalTail.push_back(x);
                fbo->faceNormalTail.push_back(y);
                fbo->faceNormalTail.push_back(z);
                x += normalLength * normal[0];
                y += normalLength * normal[1];
                z += normalLength * normal[2];
                fbo->faceNormalHead.push_back(x);
                fbo->faceNormalHead.push_back(y);
                fbo->faceNormalHead.push_back(z);
                nCount++;
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
void JFaceDraw :: updateTopologyBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    size_t numFaces = mesh->getSize(2);

    fbo->textID.clear();
    fbo->idArray.clear();
    fbo->triArray.clear();
    fbo->quadArray.clear();
    fbo->polyArray.clear();
    fbo->faceColor.clear();

    fbo->textID.reserve(numFaces);
    fbo->idArray.reserve(2*numFaces);
    fbo->faceColor.reserve(4*numFaces);

    JFaceRenderPtr fattrib;
    JNodePtr vertex;

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int nsize = face->getSize(0);
            if( nsize == 3 ) {
                int err = face->getAttribute("Render", fattrib);
                if( err ) fattrib = initRenderAttribute(face);
                if( fattrib->display) {
                    fbo->idArray.push_back(face->getID());
                    for( int j = 0; j < 3; j++) {
                        vertex = face->getNodeAt(j);
                        fbo->triArray.push_back(vertex->getID());
                    }
                    fbo->faceColor.push_back( fattrib->color[0] );
                    fbo->faceColor.push_back( fattrib->color[1] );
                    fbo->faceColor.push_back( fattrib->color[2] );
                    fbo->faceColor.push_back( fattrib->color[3] );
                }
            }
        }
    }

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int nsize = face->getSize(0);
            if( nsize == 4) {
                int err = face->getAttribute("Render", fattrib);
                if( err ) fattrib = initRenderAttribute(face);
                if( fattrib->display) {
                    fbo->idArray.push_back(face->getID());
                    for( int j = 0; j < 4; j++) {
                        vertex = face->getNodeAt(j);
                        fbo->quadArray.push_back(vertex->getID());
                    }
                    fbo->faceColor.push_back( fattrib->color[0] );
                    fbo->faceColor.push_back( fattrib->color[1] );
                    fbo->faceColor.push_back( fattrib->color[2] );
                    fbo->faceColor.push_back( fattrib->color[3] );
                }
            }
        }
    }

    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int nsize = face->getSize(0);
            if( nsize > 4) {
                int err = face->getAttribute("Render", fattrib);
                if( err ) fattrib = initRenderAttribute(face);
                if( fattrib->display) {
                    fbo->idArray.push_back(face->getID());
                    vector<int> pnodes(nsize);
                    for( int j = 0; j < nsize; j++) {
                        vertex = face->getNodeAt(j);
                        pnodes[j] = vertex->getID();
                    }
                    fbo->polyArray.push_back(pnodes);
                    fbo->faceColor.push_back( fattrib->color[0] );
                    fbo->faceColor.push_back( fattrib->color[1] );
                    fbo->faceColor.push_back( fattrib->color[2] );
                    fbo->faceColor.push_back( fattrib->color[3] );
                }
            }
        }
    }

}
////////////////////////////////////////////////////////////////////////////////
void JFaceDraw :: updateBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;
    mesh->enumerate(2);

    updateGeometryBuffers( mesh );
    updateTopologyBuffers( mesh );
}
////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: preRender()
{
    if( offset ) {
        glEnable( GL_POLYGON_OFFSET_FILL);
        glPolygonOffset( offset_factor, offset_unit);
    } else
        glDisable( GL_POLYGON_OFFSET_FILL);

    glEnable( GL_POLYGON_OFFSET_FILL);
    glEnable( GL_POLYGON_OFFSET_LINE);

    if( culling )
        glEnable( GL_CULL_FACE );
    else
        glDisable( GL_CULL_FACE );

    glFrontFace( frontface );

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    JLightsPtr lightPtr = viewManager->getLights();
    lightPtr->Switch(lights);
}
////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: postRender()
{
    glShadeModel(GL_FLAT);
    glDisable( GL_BLEND );
    glDisable( GL_CULL_FACE );
    glDisable( GL_COLOR_MATERIAL );
    glDisable(GL_POLYGON_STIPPLE);
//  glDisable( GL_POLYGON_OFFSET_FILL);
}

////////////////////////////////////////////////////////////////////////////////
void
JFaceDraw::drawNormals(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->displayNormals[2] == 0) return;

    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    glDisable( GL_LINE_SMOOTH);
    glColor3fv( &normalsColor[0] );
    glLineWidth(1.0);

    size_t numFaces = fbo->getSize();

    glDisable( GL_LIGHTING );

    glBegin(GL_LINES);
    for (size_t i = 0; i < numFaces; i++) {
        glVertex3fv( &fbo->faceNormalTail[3*i] );
        glVertex3fv( &fbo->faceNormalHead[3*i] );
    }
    glEnd();
}

////////////////////////////////////////////////////////////////////////////////

void JFaceDraw::drawPolySurface(const JMeshPtr &mesh)
{
    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    if( mrender->transparent ) {
        switch( mrender->transparencyMethod ) {
        case JRender::TRANSPARENT_SCREENDOOR:
            glDisable( GL_BLEND );
            glEnable(GL_POLYGON_STIPPLE);
            glPolygonStipple(stippleMask[screendoorPattern]);
            break;
        case JRender::TRANSPARENT_BLEND:
            glDisable(GL_POLYGON_STIPPLE);
            glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable( GL_BLEND );
            break;
        default:
            glDisable( GL_BLEND );
        }
    }

    int numTris = fbo->triArray.size()/3;

    size_t offset  = 0;
    if( numTris ) {
        for( int i = 0; i < numTris; i++) {
            int v0 = fbo->triArray[3*i];
            int v1 = fbo->triArray[3*i+1];
            int v2 = fbo->triArray[3*i+2];
            glBegin(GL_TRIANGLES);
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glEnd();
        }
    }

    offset += numTris;
    int numQuads = fbo->quadArray.size()/4;
    if( numQuads ) {
        for( int i = 0; i < numQuads; i++) {
            int v0 = fbo->quadArray[4*i];
            int v1 = fbo->quadArray[4*i+1];
            int v2 = fbo->quadArray[4*i+2];
            int v3 = fbo->quadArray[4*i+3];
            glBegin(GL_QUADS);
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glVertex3fv( &fbo->nodeCoords[3*v3]);
            glEnd();
        }
    }

    offset += numQuads;
    int numPoly = fbo->polyArray.size();
    if( numPoly ) {
        for( int i = 0; i < numPoly; i++) {
            glBegin(GL_POLYGON);
            int nn = fbo->polyArray[i].size();
            for( int j = 0; j < nn; j++) {
                int v0 = fbo->polyArray[i][j];
                glVertex3fv( &fbo->nodeCoords[3*v0]);
            }
            glEnd();
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void JFaceDraw::drawFlatSurface(const JMeshPtr &mesh)
{
    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    if( mrender->transparent ) {
        switch( mrender->transparencyMethod ) {
        case JRender::TRANSPARENT_SCREENDOOR:
            glDisable( GL_BLEND );
            glEnable(GL_POLYGON_STIPPLE);
            glPolygonStipple(stippleMask[screendoorPattern]);
            break;
        case JRender::TRANSPARENT_BLEND:
            glDisable(GL_POLYGON_STIPPLE);
            glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable( GL_BLEND );
            break;
        default:
            glDisable( GL_BLEND );
        }
    } else
        glDisable( GL_BLEND );

    if( mrender->useMaterial ) {
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, White);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
        glEnable(GL_COLOR_MATERIAL);
    }

    glShadeModel(GL_FLAT);

    int numTris = fbo->triArray.size()/3;

    size_t offset  = 0;
    if( numTris ) {
        for( int i = 0; i < numTris; i++) {
            int v0 = fbo->triArray[3*i];
            int v1 = fbo->triArray[3*i+1];
            int v2 = fbo->triArray[3*i+2];
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor4fv(  &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_TRIANGLES);
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glEnd();
        }
    }

    offset += numTris;
    int numQuads = fbo->quadArray.size()/4;
    if( numQuads ) {
        for( int i = 0; i < numQuads; i++) {
            int v0 = fbo->quadArray[4*i];
            int v1 = fbo->quadArray[4*i+1];
            int v2 = fbo->quadArray[4*i+2];
            int v3 = fbo->quadArray[4*i+3];
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor4fv(  &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_QUADS);
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glVertex3fv( &fbo->nodeCoords[3*v3]);
            glEnd();
        }
    }

    offset += numQuads;
    int numPoly = fbo->polyArray.size();
    if( numPoly ) {
        for( int i = 0; i < numPoly; i++) {
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor4fv( &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_POLYGON);
            for( size_t j = 0; j < fbo->polyArray[i].size(); j++) {
                int v0 = fbo->polyArray[i][j];
                glVertex3fv( &fbo->nodeCoords[3*v0]);
            }
            glEnd();
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

void JFaceDraw::drawSmoothSurface(const JMeshPtr &mesh)
{
    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    if( mrender->transparent ) {
        switch( mrender->transparencyMethod ) {
        case JRender::TRANSPARENT_SCREENDOOR:
            glDisable( GL_BLEND );
            glEnable(GL_POLYGON_STIPPLE);
            glPolygonStipple(stippleMask[screendoorPattern]);
            break;
        case JRender::TRANSPARENT_BLEND:
            glDisable(GL_POLYGON_STIPPLE);
            glEnable( GL_BLEND );
            break;
        default:
            glDisable( GL_BLEND );
        }
    }

    if( mrender->useMaterial ) {
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, White);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100);
        glEnable(GL_COLOR_MATERIAL);
    }

    glShadeModel(GL_SMOOTH);

    int numTris = fbo->triArray.size()/3;

    size_t offset  = 0;
    if( numTris ) {
        for( int i = 0; i < numTris; i++) {
            int v0 = fbo->triArray[3*i];
            int v1 = fbo->triArray[3*i+1];
            int v2 = fbo->triArray[3*i+2];
            glBegin(GL_TRIANGLES);
            glNormal3fv( &fbo->nodeNormal[3*v0] );
            glColor3fv(  &fbo->nodeColor[4*v0] );
            glVertex3fv( &fbo->nodeCoords[3*v0]);

            glNormal3fv( &fbo->nodeNormal[3*v1] );
            glColor3fv(  &fbo->nodeColor[4*v1] );
            glVertex3fv( &fbo->nodeCoords[3*v1]);

            glNormal3fv( &fbo->nodeNormal[3*v2] );
            glColor3fv(  &fbo->nodeColor[4*v2] );
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glEnd();
        }
    }

    offset += numTris;
    int numQuads = fbo->quadArray.size()/4;
    if( numQuads ) {
        for( int i = 0; i < numQuads; i++) {
            int v0 = fbo->quadArray[4*i];
            int v1 = fbo->quadArray[4*i+1];
            int v2 = fbo->quadArray[4*i+2];
            int v3 = fbo->quadArray[4*i+3];
            glBegin(GL_QUADS);
            glNormal3fv( &fbo->nodeNormal[3*v0] );
            glColor3fv(  &fbo->nodeColor[4*v0] );
            glVertex3fv( &fbo->nodeCoords[3*v0]);

            glNormal3fv( &fbo->nodeNormal[3*v1] );
            glColor3fv(  &fbo->nodeColor[4*v1] );
            glVertex3fv( &fbo->nodeCoords[3*v1]);

            glNormal3fv( &fbo->nodeNormal[3*v2] );
            glColor3fv(  &fbo->nodeColor[4*v2] );
            glVertex3fv( &fbo->nodeCoords[3*v2]);

            glNormal3fv( &fbo->nodeNormal[3*v3] );
            glColor3fv(  &fbo->nodeColor[4*v3] );
            glVertex3fv( &fbo->nodeCoords[3*v3]);
            glEnd();
        }
    }

    offset += numQuads;
    int numPoly = fbo->polyArray.size();
    if( numPoly ) {
        for( int i = 0; i < numPoly; i++) {
            glBegin(GL_POLYGON);
            int nn = fbo->polyArray[i].size();
            for( int j = 0; j < nn; j++) {
                int v0 = fbo->polyArray[i][j];
                glNormal3fv( &fbo->nodeNormal[3*v0] );
                glColor3fv(  &fbo->nodeColor[4*v0] );
                glVertex3fv( &fbo->nodeCoords[3*v0]);
            }
            glEnd();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void JFaceDraw::drawTexturedSurface(const JMeshPtr &mesh)
{
    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    GLuint texid;
    int err = mesh->getAttribute("TextureID", texid);
    if( err ) {
        cout << "Warning: No texture is bound with the mesh" << endl;
        return;
    }

    glBindTexture( GL_TEXTURE_2D, texid);
    glEnable(GL_TEXTURE_2D);

    int numTris = fbo->triArray.size()/3;

    size_t offset  = 0;
    if( numTris ) {
        for( int i = 0; i < numTris; i++) {
            int v0 = fbo->triArray[3*i];
            int v1 = fbo->triArray[3*i+1];
            int v2 = fbo->triArray[3*i+2];
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor3fv(  &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_TRIANGLES);
            glTexCoord2fv(&fbo->uvCoords[2*v0] );
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glTexCoord2fv(&fbo->uvCoords[2*v1] );
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glTexCoord2fv(&fbo->uvCoords[2*v2] );
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glEnd();
        }
    }

    offset += numTris;
    int numQuads = fbo->quadArray.size()/4;
    if( numQuads ) {
        for( int i = 0; i < numQuads; i++) {
            int v0 = fbo->quadArray[4*i];
            int v1 = fbo->quadArray[4*i+1];
            int v2 = fbo->quadArray[4*i+2];
            int v3 = fbo->quadArray[4*i+3];
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor4fv(  &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_QUADS);
            glTexCoord2fv(&fbo->uvCoords[2*v0] );
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glTexCoord2fv(&fbo->uvCoords[2*v1] );
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glTexCoord2fv(&fbo->uvCoords[2*v2] );
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glTexCoord2fv(&fbo->uvCoords[2*v3] );
            glVertex3fv( &fbo->nodeCoords[3*v3]);
            glEnd();
        }
    }

    offset += numQuads;
    int numPoly = fbo->polyArray.size();
    if( numPoly ) {
        for( int i = 0; i < numPoly; i++) {
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor4fv( &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_POLYGON);
            for( size_t j = 0; j < fbo->polyArray[i].size(); j++) {
                int v0 = fbo->polyArray[i][j];
                glTexCoord2fv(&fbo->uvCoords[2*v0] );
                glVertex3fv( &fbo->nodeCoords[3*v0]);
            }
            glEnd();
        }
    }
    glDisable(GL_TEXTURE_2D);
}

////////////////////////////////////////////////////////////////////////////////

void JFaceDraw::drawHiddenlines(const JMeshPtr &mesh)
{
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f( 1.0, 1.0, 1.0);
    drawPolySurface( mesh );

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(2.0, 2.0);
    drawFlatSurface( mesh );
    glDisable(GL_POLYGON_OFFSET_FILL);

    glEnable( GL_POLYGON_OFFSET_POINT);
    glEnable( GL_POLYGON_OFFSET_LINE);

    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    drawFlatSurface( mesh );
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDepthFunc(GL_LEQUAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawFlatSurface( mesh );
    glDepthFunc(GL_LESS);
}

////////////////////////////////////////////////////////////////////////////////

void
JFaceDraw::draw(const JFacePtr &face)
{
    int numnodes = face->getSize(0);
    JFaceRenderPtr attrib;
    int err = face->getAttribute("Render", attrib);
    if( err ) return;

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    Point3F normal;
    if( numnodes == 3) {
        face->getAttribute("Normal", normal);
        glNormal3fv( &normal[0]);
        glColor3fv( &attrib->color[0] );
        glBegin(GL_TRIANGLES);
        for( int i = 0; i < 3; i++) {
            JNodePtr vtx = face->getNodeAt(i);
            const Point3D &p3d = vtx->getXYZCoords();
            glVertex3dv( &p3d[0]);
        }
        glEnd();
        return;
    }

    if( numnodes == 4) {
        face->getAttribute("Normal", normal);
        glNormal3fv( &normal[0]);
        glColor3fv( &attrib->color[0] );
        glBegin(GL_QUADS);
        for( int i = 0; i < 4; i++) {
            JNodePtr vtx = face->getNodeAt(i);
            const Point3D &p3d = vtx->getXYZCoords();
            glVertex3dv( &p3d[0]);
        }
        glEnd();
        return;
    }
    if( numnodes == 4) {
        face->getAttribute("Normal", normal);
        glNormal3fv( &normal[0]);
        glColor3fv( &attrib->color[0] );
        glBegin(GL_POLYGON);
        for( int i = 0; i < numnodes; i++) {
            JNodePtr vtx = face->getNodeAt(i);
            const Point3D &p3d = vtx->getXYZCoords();
            glVertex3dv( &p3d[0]);
        }
        glEnd();
        return;
    }
}
///////////////////////////////////////////////////////////////////////////////
void
JFaceDraw::draw(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    preRender();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);

    glEnable( GL_BLEND);

    int renderStyle = mrender->getRenderStyle();

    if( renderStyle == JRender::HIDDENLINES) {
        glDisable(GL_LIGHTING);
        drawHiddenlines(mesh);
    } else {
        if( renderStyle == JRender::WIREFRAME)
            glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
        else
            glPolygonMode( GL_FRONT_AND_BACK, GL_FILL);

        int shade = mrender->getSurfaceShade();

        if( shade == JRender::FLAT_SHADE)
            drawFlatSurface(mesh);
        else
            drawSmoothSurface(mesh);
    }

    postRender();
    drawNormals(mesh);
    drawIDs(mesh);
}
/////////////////////////////////////////////////////////////////////////////

void
JFaceDraw::drawIDs( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->displayIDs[2] == 0) return;

    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    glDisable( GL_LIGHTING );
    glDisable (GL_BLEND);
    glColor3f( 0.1, 0.1, 0.1);

    Vec camdir = -viewManager->camera()->viewDirection();

    glEnable( GL_DEPTH_TEST);

    glPushMatrix();
    QFont newFont("Times", 12);

    Vec pvec, normal;
    viewManager->startScreenCoordinatesSystem();
    for( size_t id  : fbo->idArray) {
        normal[0] = fbo->faceNormal[3*id+0];
        normal[1] = fbo->faceNormal[3*id+1];
        normal[2] = fbo->faceNormal[3*id+2];
        if( normal[0]*camdir[0] + normal[1]*camdir[1] + normal[2]*camdir[2] >= 0.0)  {
            pvec[0]  = fbo->faceNormalTail[3*id+0];
            pvec[1]  = fbo->faceNormalTail[3*id+1];
            pvec[2]  = fbo->faceNormalTail[3*id+2];
            pvec     = viewManager->camera()->projectedCoordinatesOf(pvec);
            viewManager->drawText( pvec[0], pvec[1], QString::number(id), newFont );
        }
    }
    glPopMatrix();
    viewManager->stopScreenCoordinatesSystem();

}
/////////////////////////////////////////////////////////////////////////////

void
JFaceDraw::withName( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender->pickableEntity != 2) return;

    boost::shared_ptr<JFaceBuffers> fbo = mrender->faceBuffers;

    preRender();

    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glDisable( GL_BLEND );

    int numTris = fbo->triArray.size()/3;

    size_t offset  = 0;
    if( numTris ) {
        for( int i = 0; i < numTris; i++) {
            glPushName(fbo->idArray[i]);
            int v0 = fbo->triArray[3*i];
            int v1 = fbo->triArray[3*i+1];
            int v2 = fbo->triArray[3*i+2];
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor3fv(  &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_TRIANGLES);
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glEnd();
            glPopName();
        }
    }

    offset += numTris;
    int numQuads = fbo->quadArray.size()/4;
    if( numQuads ) {
        for( int i = 0; i < numQuads; i++) {
            glPushName(fbo->idArray[i+offset]);
            int v0 = fbo->quadArray[4*i];
            int v1 = fbo->quadArray[4*i+1];
            int v2 = fbo->quadArray[4*i+2];
            int v3 = fbo->quadArray[4*i+3];
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor4fv(  &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_QUADS);
            glVertex3fv( &fbo->nodeCoords[3*v0]);
            glVertex3fv( &fbo->nodeCoords[3*v1]);
            glVertex3fv( &fbo->nodeCoords[3*v2]);
            glVertex3fv( &fbo->nodeCoords[3*v3]);
            glEnd();
            glPopName();
        }
    }

    offset += numQuads;
    int numPoly = fbo->polyArray.size();
    if( numPoly ) {
        for( int i = 0; i < numPoly; i++) {
            glPushName(fbo->idArray[i+offset]);
            glNormal3fv( &fbo->faceNormal[3*(i+offset)] );
            glColor4fv( &fbo->faceColor[4*(i+offset)] );
            glBegin(GL_POLYGON);
            for( size_t j = 0; j < fbo->polyArray[i].size(); j++) {
                int v0 = fbo->polyArray[i][j];
                glVertex3fv( &fbo->nodeCoords[3*v0]);
            }
            glEnd();
            glPopName();
        }
    }
    postRender();
}
/////////////////////////////////////////////////////////////////////////////

JCellDraw :: JCellDraw()
{
}

/////////////////////////////////////////////////////////////////////////////
void JCellDraw :: initRenderAttributes(const JMeshPtr &mesh)
{
    size_t numCells = mesh->getSize(3);
    JCellRenderPtr attrib;
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() && !cell->hasAttribute("Render")  ) {
            attrib = JCellRender::newObject();
            cell->setAttribute("Render", attrib);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////

void JCellDraw :: activateLowerEntities( const JCellPtr &cell)
{
    JCellRenderPtr   cattrib;
    cell->getAttribute("Render", cattrib);
    bool val = cattrib->display;

    JFaceRenderPtr   fattrib;
    for( int i = 0; i < cell->getSize(2); i++) {
        const JFacePtr &face = cell->getFaceAt(i);
        int err = face->getAttribute("Render", fattrib);
        if( !err &&  val ) fattrib->display = val;
    }

    JEdgeRenderPtr   eattrib;
    for( int i = 0; i < cell->getSize(1); i++) {
        const JEdgePtr &edge = cell->getEdgeAt(i);
        int err = edge->getAttribute("Render", eattrib);
        if( !err && val ) eattrib->display = val;
    }

    JNodeRenderPtr   vattrib;
    for( int i = 0; i < cell->getSize(0); i++) {
        const JNodePtr &vtx = cell->getNodeAt(i);
        int err = vtx->getAttribute("Render", vattrib);
        if( !err && val) vattrib->display = val;
    }
}

/////////////////////////////////////////////////////////////////////////////
void JCellDraw :: updateBuffers( const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    size_t numCells = mesh->getSize(3);
    if( numCells == 0) return;

    JCellRenderPtr  cattrib;
    JFaceRenderPtr  fattrib;
    JEdgeRenderPtr  eattrib;
    JNodeRenderPtr  vattrib;

    /*
        FaceSet faceSet;
        for( size_t i = 0; i < numCells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                for( int i = 0; i < cell->getSize(2); i++) {
                    const JFacePtr &face = cell->getFaceAt(i);
                    int err = face->getAttribute("Render", fattrib);
                    if( !err) fattrib->display = 0;
                }

                for( int i = 0; i < cell->getSize(1); i++) {
                    const JEdgePtr &edge = cell->getEdgeAt(i);
                    int err = edge->getAttribute("Render", eattrib);
                    if(!err) eattrib->display = 0;
                }

                for( int i = 0; i < cell->getSize(0); i++) {
                    const JNodePtr &vtx = cell->getNodeAt(i);
                    int err = vtx->getAttribute("Render", vattrib);
                    if( !err) vattrib->display = 0;
                }
            }
        }
    */

    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            int err = cell->getAttribute("Render", cattrib);
            if( err) {
                cattrib = JCellRender::newObject();
                cell->setAttribute("Render", cattrib);
            }
            activateLowerEntities(cell);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
void
JCellDraw::draw( const JCellPtr &)
{
}

void
JCellDraw::draw( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;
}
/////////////////////////////////////////////////////////////////////////////

#ifdef DEACTIVATE
void
JCellDraw::withName( const JMeshPtr &mesh)
{
   if( mesh == nullptr ) return;
    JMeshRenderPtr mrender;
    int err = mesh->getAttribute("Render", mrender);
    if( mrender->pickableEntity == 3) {
        size_t numcells = mesh->getSize(3);
        for (size_t i = 0; i < numcells; i++) {
            JCellPtr cell = mesh->getCellAt(i);
            withName(cell, i);
        }
     }
}
/////////////////////////////////////////////////////////////////////////////
void
JCellDraw::drawIDs( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

        size_t nSize = mesh->getSize(3);
        size_t index = 0;
        for( size_t i = 0; i < nSize; i++) {
                const JCellPtr &cell = mesh->getCellAt(i);
                if( cell->isActive() ) {
                    cell->setID(index++);
                    draw_id(cell);
                }
            }
         }
}

/////////////////////////////////////////////////////////////////////////////

void JCellDraw :: draw(const JCellPtr &cell)
{
    if( !isDrawable(cell) ) return;

    bool val = 1;
    int ori;

    if( display_lower_faces ) {
        cell->getFaces( faces );
        int numfaces = faces.size();
        for( int iface = 0; iface < numfaces; iface++) {
            JFacePtr face = cell->getFaceAt(iface, ori);
            face->setAttribute("Display", val);
        }
    }

    if( display_lower_edges ) {
        int numedges = cell->getSize(1);
        for( int i = 0; i < numedges; i++) {
            JEdgePtr edge = cell->getEdgeAt(i);
            edge->setAttribute("Display", val);
        }
    }

    if( display_lower_nodes ) {
        int numnodes = cell->getSize(0);
        for( int i = 0; i < numnodes; i++) {
            JNodePtr vertex= cell->getNodeAt(i);
            vertex->setAttribute("Display", val);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////

void JCellDraw :: withName(const JCellPtr &cell, size_t pos)
{
    glPushMatrix();
    glPushName(pos);
    draw(cell);
    glEnd();
    glPopName();
    glPopMatrix();
}

/////////////////////////////////////////////////////////////////////////////
void JCellDraw :: drawID( const JCellPtr &cell)
{
    if( !isDrawable(cell) ) return;

    char number[128];
    size_t id = cell->getID();
    sprintf(number, "%ld", id);

    Point3D xyz;
    cell->getAvgXYZ(xyz);
//  drawTextAt(number,xyz);
}

/////////////////////////////////////////////////////////////////////////////

void
JCellDraw::draw( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    size_t numcells = mesh->getSize(3);
    for (size_t i = 0; i < numcells; i++) {
        Cell *cell = mesh->getCellAt(i);
        draw(cell);
    }
}
/////////////////////////////////////////////////////////////////////////////
void
JCellDraw::withName( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    size_t numcells = mesh->getSize(3);
    for (size_t i = 0; i < numcells; i++) {
        Cell *cell = mesh->getCellAt(i);
        withName(cell, i);
    }
}
/////////////////////////////////////////////////////////////////////////////
void
JCellDraw::drawIDs( const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    size_t nSize = mesh->getSize(3);
    size_t index = 0;
    for( size_t i = 0; i < nSize; i++) {
        Cell *cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            cell->setID(index++);
            drawID(cell);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////

void JNodeDraw :: draw(const JNodePtr &vertex ) const
{
    if( vertex == nullptr ) return;
    if( !vertex->isActive() ) return;

    JNodeRenderPtr nAttrib;

    vertex->getAttribute("Render", nAttrib);
    if( nAttrib == nullptr ) return;

    if(  nAttrib->display ==  0) return;

    glColor4fv( &nAttrib->color[0] );

    int glyph    =  nAttrib->glyph;
    float scale  =  nAttrib->scale;
    float radius;

    if( glyph == JNodeRender::NODE_AS_SPHERE ) {
        JLights::getInstance().Switch(1);
        glEnable( GL_DEPTH_TEST );
        glMaterialfv(GL_FRONT, GL_SPECULAR, White);
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT, GL_SHININESS, 100);
        glEnable(GL_COLOR_MATERIAL);
        const Point3D &xyz = vertex->getXYZCoords();
        glPushMatrix();
        glTranslatef( xyz[0], xyz[1], xyz[2] );
        radius =  nAttrib->ballRadius;
        gluSphere( sphereObj, radius*scale, numSlices, numStacks);
        glPopMatrix();
        return;
    }

    if( glyph == JNodeRender::NODE_AS_CIRCLE ) {
        JLights::getInstance().Switch(0);
        Vec3F srcVec, dstVec, perpAxis;
        vertex->getAttribute("Normal", dstVec);
        srcVec[0] = 0.0;
        srcVec[1] = 0.0;
        srcVec[2] = 1.0;
        JMath::cross_product( dstVec, srcVec, perpAxis);
        double angle = JMath::getVecAngle(srcVec, dstVec, ANGLE_IN_RADIANS);
        qglviewer::Vec rotaxis(perpAxis[0], perpAxis[1], perpAxis[2] );
        qglviewer::Quaternion quaternion(rotaxis, -1.0*angle);
        const Point3D &xyz = vertex->getXYZCoords();
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glMaterialf(GL_FRONT, GL_SHININESS, 100);
        glNormal3f( dstVec[0], dstVec[1], dstVec[2] );
        glPushMatrix();
        glTranslatef( xyz[0], xyz[1], xyz[2] );
        GLdouble mat[16];
        quaternion.getMatrix(mat);
        glMultMatrixd( mat );
        radius =  nAttrib->ballRadius;
        gluDisk( diskObj, 0.0, radius*scale, numSlices, 1);
        glPopMatrix();
        return;
    }

    JLights::getInstance().Switch(0);
    radius =  nAttrib->pointSize;
    glPointSize( radius*scale );

    glColor3f( 1.0, 0.0, 0.0);

    glBegin(GL_POINTS);
    const Point3D &xyz = vertex->getXYZCoords();
    glVertex3f(xyz[0], xyz[1], xyz[2] );
    glEnd();
}
////////////////////////////////////////////////////////////////////////////////
void JNodeDraw :: drawNormal(const JNodePtr &vertex) const
{
    if( vertex == nullptr ) return;
    if( !vertex->isActive() )  return;

    JNodeRenderPtr nodeAttrib;
    vertex->getAttribute("Render", nodeAttrib);
    if( nodeAttrib->display == 0)  return;

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
}

////////////////////////////////////////////////////////////////////////////////
void JNodeDraw :: drawID(const JNodePtr &vertex) const
{
    if( vertex == nullptr ) return;
    if( !vertex->isActive() )  return;

    JNodeRenderPtr nodeAttrib;
    vertex->getAttribute("Render", nodeAttrib);
    if( nodeAttrib->display == 0)  return;

    char number[128];
    const  Point3D  &xyz = vertex->getXYZCoords();
    size_t vid = vertex->getID();
    sprintf(number, "%ld", vid);
//  drawTextAt(number, xyz);
}
///////////////////////////////////////////////////////////////////////////////

void JNodeDraw :: withName( const JNodePtr &vertex, size_t pos) const
{
    if( !vertex->isActive() ) return;

    glPushMatrix();
    glPushName(pos);
//  draw(vertex);
    glPopName();
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: draw( const JEdgePtr &edge) const
{
    if( edge  == nullptr ) return ;
    if( !edge->isActive() ) return;

    JEdgeRenderPtr edgeAttrib;
    edge->getAttribute("Render", edgeAttrib);
    if( edgeAttrib == nullptr ) return;
    if( edgeAttrib->display == 0) return;

    JNodePtr v0, vm, v1;

    int glyph = edgeAttrib->glyph;

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
        break;
    }

    glColor4fv( &edgeAttrib->color[0] );

    if( glyph ) {
        JLights::getInstance().Switch(1);
        double radius = edgeAttrib->cylinderRadius*edgeAttrib->scale;
        drawCylinder(edge, radius);
    } else {
        float radius = edgeAttrib->lineWidth*edgeAttrib->scale;
        JLights::getInstance().Switch(0);
        Point3D xyz;
        glLineWidth(radius);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glBegin(GL_LINES);
        v0 = edge->getNodeAt(0);
        xyz = v0->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
        if( withSteiner) {
            if( edge->hasAttribute("Steiner") ) {
                edge->getAttribute("Steiner", vm);
                xyz = vm->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
        }
        v1 = edge->getNodeAt(1);
        xyz = v1->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
        glEnd();
    }

    if( display_lower_nodes ) {
        bool val = 1;
        JNodeRenderPtr nAttrib;
        edge->getNodeAt(0)->getAttribute("Render", nAttrib);
        nAttrib->display = val;
        edge->getNodeAt(1)->getAttribute("Render", nAttrib);
        nAttrib->display = val;
    }
}
////////////////////////////////////////////////////////////////////////////////

void JEdgeDraw :: withName(const JEdgePtr &edge, size_t pos) const
{
    if( !edge->isActive() ) return;

    glPushMatrix();
    glPushName(pos);
    draw(edge);
    glPopName();
    glPopMatrix();
}
////////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: drawID( const JEdgePtr &edge) const
{
    if( edge  == nullptr ) return ;
    if( !edge->isActive() ) return;
    JEdgeRenderPtr edgeAttrib;
    edge->getAttribute("Render", edgeAttrib);
    if( edgeAttrib == nullptr ) return;
    if( edgeAttrib->display == 0) return;

    char number[128];
    size_t eid = edge->getID();
    sprintf(number, "%ld", eid);

    Point3D xyz;
    Vertex::getMidPoint( edge->getNodeAt(0), edge->getNodeAt(1), xyz );
//  drawTextAt(number, xyz);
}
////////////////////////////////////////////////////////////////////////////////

void
JFaceDraw::draw_hidden_lines( const JMeshPtr &)
{
    /*
         if( display_entity[1] == 0 ) return;

         clr  = meshEdgeColor->getInternalColor();
         glColor3f( clr[0], clr[1], clr[2] );

         int fstatus = display_entity[2];
         int filltype = drawFace->getStyle();

         display_entity[2] = 1;

         switch( hiddenlinesMethod ) {
         case HIDDENLINE_WITH_BACKFACES_CULL:
              drawFace->setBackfaceCull(1);
              drawFace->setStyle(JFaceDraw::FACE_LINES);
              draw_faces();
              drawFace->setBackfaceCull(0);
              break;
         case HIDDENLINE_WITH_FRONTLINES:
              clr  = viewManager->getBackgroundColor();
              glColor3f( clr[0], clr[1], clr[2] );
              glPolygonMode(GL_FRONT, GL_FILL);
              drawFace->setColorMethod( nullptr );
              drawFace->setStyle( JFaceDraw::FACE_FILL);
              draw_faces();

              clr  = meshEdgeColor->getInternalColor();
              glColor3f( clr[0], clr[1], clr[2] );
              drawFace->setStyle(JFaceDraw::FACE_LINES);
              draw_faces();
              break;
         case HIDDENLINE_WITH_DEPTH_TEST:
              glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
              drawFace->setStyle( JFaceDraw::FACE_FILL);
              draw_faces();

              glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
              glDepthFunc(GL_LEQUAL);
              drawFace->setStyle(JFaceDraw::FACE_LINES);
              draw_faces();
              glDepthFunc(GL_LESS);
              break;
         }
         display_entity[2] = fstatus;
         drawFace->setStyle( filltype );
    */
}

///////////////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: draw(const JFacePtr &face, int ori ) const
{
    if( face  == nullptr ) return;
    if( !face->isActive() ) return;

    JFaceRender *fAttrib = nullptr;
    face->getAttribute("Render", fAttrib);
    if( fAttrib == nullptr ) return;
    if( fAttrib->display == 0) return;

    glColor4fv( &fAttrib->color[0]);

    if( use_material ) {
        if( face->hasAttribute("Normal") ) {
            Vec3F normal;
            face->getAttribute( "Normal", normal );
            float nx = 1.0 * normalSign * normal[0];
            float ny = 1.0 * normalSign * normal[1];
            float nz = 1.0 * normalSign * normal[2];
            glNormal3f(nx, ny, nz);
        }
    }

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

    if( display_lower_nodes ) {
        int numnodes = face->getSize(0);
        bool val = 1;
        for( int i = 0; i < numnodes; i++)  {
            JNodePtr vtx = face->getNodeAt(i);
            vtx->setAttribute("Display", val);
        }
    }

    if( display_lower_edges ) {
        int numedges = face->getSize(1);
        bool val = 1;
        for( int i = 0; i < numedges; i++)  {
            JEdgePtr edge = face->getEdgeAt(i);
            edge->setAttribute("Display", val);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: withName(const JFacePtr &face, size_t pos) const
{
    if( !face->isActive() ) return;

    glPushName(pos);
    draw(face);
    glPopName();
    glPopMatrix();
}
////////////////////////////////////////////////////////////////////////////////

void
JFaceDraw::withName(const JMeshPtr &mesh)
{
    if( mesh == nullptr ) return;

    preRender();

    dfdsfsdfs
    size_t nSize = mesh->getSize(2);
    for( size_t i = 0; i < nSize; i++) {
        Face *face = mesh->getFaceAt(i);
        withName(face, i);
    }

    postRender();
}

////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: draw_tri( const JFacePtr &face, int ori) const
{
    if( !face->isActive() ) return;

    if( shade == FLAT_SHADE ) {
        if( ori < 0) {
            glBegin(GL_TRIANGLES);
            for (int j = 0; j < 3; j++) {
                JNodePtr vtx = face->getNodeAt(2-j);
                const Point3D &xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
            glEnd();
            return;
        }

        glBegin(GL_TRIANGLES);
        for (int j = 0; j < 3; j++) {
            JNodePtr vtx = face->getNodeAt(j);
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }

    JNodeRenderPtr nAttrib;

    if( ori < 0) {
        glBegin(GL_TRIANGLES);
        for (int j = 0; j < 3; j++) {
            JNodePtr vtx = face->getNodeAt(2-j);
            vtx->getAttribute("Render", nAttrib);
            glColor4fv( &nAttrib->color[0] );
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }

    glBegin(GL_TRIANGLES);
    for (int j = 0; j < 3; j++) {
        JNodePtr vtx = face->getNodeAt(j);
        vtx->getAttribute("Render", nAttrib);
        glColor4fv( &nAttrib->color[0] );
        const Point3D &xyz = vtx->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
    return;
}
////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: draw_quad( const JFacePtr &face, int ori) const
{
    if( !face->isActive() ) return;

    if( shade == FLAT_SHADE ) {
        if( ori < 0) {
            glBegin(GL_QUADS);
            for (int j = 0; j < 4; j++) {
                JNodePtr vtx = face->getNodeAt(3-j);
                const Point3D &xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], xyz[2]);
            }
            glEnd();
            return;
        }

        glBegin(GL_QUADS);
        for (int j = 0; j < 4; j++) {
            JNodePtr vtx = face->getNodeAt(j);
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }

    JNodeRenderPtr nAttrib;

    if( ori < 0) {
        glBegin(GL_QUADS);
        for (int j = 0; j < 4; j++) {
            JNodePtr vtx = face->getNodeAt(3-j);
            vtx->getAttribute("Render", nAttrib);
            glColor4fv( &nAttrib->color[0] );
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }

    glBegin(GL_QUADS);
    for (int j = 0; j < 4; j++) {
        JNodePtr vtx = face->getNodeAt(j);
        vtx->getAttribute("Render", nAttrib);
        glColor4fv( &nAttrib->color[0] );
        const Point3D &xyz = vtx->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();

}

////////////////////////////////////////////////////////////////////////////////
void JFaceDraw :: draw_poly( const JFacePtr &face, int ori) const
{
    if( !face->isActive() ) return;

    int nnodes = face->getSize(0);

    if( ori < 0) {
        glBegin(GL_POLYGON);
        for (int j = 0; j < nnodes; j++) {
            JNodePtr vtx = face->getNodeAt(nnodes-1-j);
            const Point3D &xyz = vtx->getXYZCoords();
            glVertex3f(xyz[0], xyz[1], xyz[2]);
        }
        glEnd();
        return;
    }

    glBegin(GL_POLYGON);
    for (int j = 0; j < nnodes; j++) {
        JNodePtr vtx = face->getNodeAt(j);
        const Point3D &xyz = vtx->getXYZCoords();
        glVertex3f(xyz[0], xyz[1], xyz[2]);
    }
    glEnd();
}
////////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: drawNormal( const JFacePtr &face ) const
{
    if( face  == nullptr ) return;
    if( !face->isActive() ) return;

    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);
    if( fAttrib == nullptr ) return;
    if( fAttrib->display == 0) return;

    assert( normalLength > 0.0);

    Vec3F   normal;
    face->getAttribute("Normal", normal);

    Point3D centroid;
    face->getAvgXYZ( centroid );

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
}

/////////////////////////////////////////////////////////////////////////////

void JFaceDraw :: drawID( const JFacePtr &face) const
{
    if( !face->isActive() ) return;

    char number[128];
    size_t fid = face->getID();
    sprintf(number, "%ld", fid);

    Point3D xyz;
    face->getAvgXYZ(xyz);
//  drawTextAt(number,xyz);
}


/////////////////////////////////////////////////////////////////////////////
void JEdgeDraw :: activateLowerEntities(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JNodeRender  attrib;
    size_t numEdges = mesh->getSize(1);
    for( int i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            for( int j = 0; j < 2; j++) {
                JNodePtr v0 = edge->getNodeAt(0);
                v0->getAttribute("Render", attrib);
                attrib->display = 1;
            }
        }
    }
}
/////////////////////////////////////////////////////////////////////////////
void JFaceDraw :: activateLowerEntities(const JMeshPtr &mesh)
{
    if( mesh == nullptr) return;

    JEdgeRenderPtr attrib;
    size_t numFaces = mesh->getSize(2);
    for( int i = 0; i < numFaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int numEdges = face->getSize(1);
            for( int j = 0; j < numEdges; j++) {
                Edge *edge = face->getEdgeAt(j);
                edge->getAttribute("Render", attrib);
                attrib->display = 1;
            }
        }
    }
}

#endif

