#include "JaalViewer.hpp"

void JViewComponent :: displayList( GLuint &dlistID )
{
    if( dlistID > 0) {
        glCallList(dlistID);
        return;
    } else {
        dlistID = glGenLists(1) + 1;
        glNewList( dlistID, GL_COMPILE_AND_EXECUTE);
    }
}

////////////////////////////////////////////////////////////////////////////////

JaalViewer :: JaalViewer(QWidget *w) : QGLViewer(w)
{
    parent = w;
    sceneFloor.reset( new JSceneFloor);

    lights.reset( new JLights);
//  lights_placement_policy = MOVE_LIGHTS_WITH_AXIS;

    scaleFactor    = 1.0;
    saveAnimation  = 1;
    numAnimFrames  = 360;
    currAnimFrame = 0;
    mousetracking = 0;
    scene_radius  = 1.0;

    displayAxis   = 0;
    displayBox    = 0;
    displayFloor  = 0;
    frozenView    = 0;

    Point3D p0;
    p0[0] = -1.0;
    p0[1] = -1.0;
    p0[2] = -1.0;
    Point3D p1;
    p1[0] = 1.0;
    p1[1] = 1.0;
    p1[2] = 1.0;
    box.setPoints(p0, p1);
    worldConstraint = new WorldConstraint();
    qtWaitSpinner.reset( new QtWaitingSpinner(this));
    init();
}

///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::init()
{
    backgroundColor[0] = 1.0;
    backgroundColor[1] = 1.0;
    backgroundColor[2] = 1.0;

    perspective_view = 1;
    projectionType = PERSPECTIVE_PROJECTION;
    viewer_dim  = 3;
    faceCulling = 0;

    if( projectionType == ORTHOGRAPHIC_PROJECTION )
        camera()->setType(Camera::ORTHOGRAPHIC);
    else
        camera()->setType(Camera::PERSPECTIVE);

    setRotationAxis(0);
    int w = parent->width();
    int h = parent->height();
    camera()->setScreenWidthAndHeight(w,h);
}

///////////////////////////////////////////////////////////////////////////////

void  JaalViewer :: freezeView( bool f ) {
    frozenView = f;
    if( f ) {
        setRotationAxis(NO_ROTATION);
        setTranslationAxis(NO_TRANSLATION);
    } else {
        setRotationAxis(FREE_ROTATION);
        setTranslationAxis(FREE_TRANSLATION);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JaalViewer:: paintEvent(QPaintEvent *e)
{
    for( size_t i = 0; i < listeners.size(); i++)
        QApplication::sendEvent(listeners[i], e);

    QGLViewer::paintEvent(e);
}

bool
JaalViewer::notify( QObject *, QEvent *)
{
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
void
JaalViewer::keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_A) {
        QGLViewer::drawAxis( sceneRadius() );
    }

    if( e->key() == Qt::Key_Enter) {
        return;
    }
    QGLViewer::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////
int JaalViewer :: getPixelPosition( const Point2I &pixelPos, Point3D &xyz) const
{
    QPoint qp(pixelPos[0], pixelPos[1] );
    bool found;
    Vec p = camera()->pointUnderPixel( qp, found);
    if( found ) {
        xyz[0] = p.x;
        xyz[1] = p.y;
        xyz[2] = p.z;
        return 0;
    }
    cout << "Warning: Point under pixel not found " << endl;
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
void
JaalViewer::wheelEvent( QWheelEvent *e)
{
    for( size_t i = 0; i < listeners.size(); i++)
        QApplication::sendEvent(listeners[i], e);
    QGLViewer::wheelEvent(e);
}

void
JaalViewer::mousePressEvent( QMouseEvent *e)
{
    notifyMouseEvent(1);

    int x = e->x();
    int y = e->y();

    x = max(0,x);
    y = max(0,y);
    x = min(x, camera()->screenWidth() );
    y = min(y, camera()->screenHeight() );

    mouseStartPixelPos[0] = x;
    mouseStartPixelPos[1] = y;

    mouseCurrPixelPos[0]  = x;
    mouseCurrPixelPos[1]  = y;

    if( frozenView ) return;

    for( size_t i = 0; i < listeners.size(); i++)
        QApplication::sendEvent(listeners[i], e);

    QGLViewer::mousePressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::mouseReleaseEvent( QMouseEvent *e)
{
    notifyMouseEvent(3);

    int x = e->x();
    int y = e->y();

    x = max(0,x);
    y = max(0,y);
    x = min(x, camera()->screenWidth() );
    y = min(y, camera()->screenHeight() );

    mouseEndPixelPos[0]  = x;
    mouseEndPixelPos[1]  = y;

    mouseCurrPixelPos[0] = x;
    mouseCurrPixelPos[1] = y;

    if( frozenView ) return;

    for( size_t i = 0; i < listeners.size(); i++)
        QApplication::sendEvent(listeners[i], e);
    QGLViewer::mouseReleaseEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::mouseMoveEvent( QMouseEvent *e)
{
    if( frozenView ) return;

    notifyMouseEvent(2);

    int x = e->x();
    int y = e->y();

    x = max(0,x);
    y = max(0,y);
    x = min(x, camera()->screenWidth() );
    y = min(y, camera()->screenHeight() );

    mouseCurrPixelPos[0] = x;
    mouseCurrPixelPos[1] = y;

    for( size_t i = 0; i < listeners.size(); i++)
        QApplication::sendEvent(listeners[i], e);

    QGLViewer::mouseMoveEvent(e);
}

/////////////////////////////////////////////////////////////////////////////////

void JaalViewer:: showEvent( QShowEvent *)
{
    int w = parent->width();
    int h = parent->height();
    this->resize(w,h);
    updateGL();
}

///////////////////////////////////////////////////////////////////////////////
void JaalViewer::endSelection(const QPoint&)
{
    picked_entities.clear();
    // Flush GL buffers
    glFlush();

    // Get the number of objects that were seen through the pick matrix frustum. Reset GL_RENDER mode.
    GLint nbHits = glRenderMode(GL_RENDER);

    if (nbHits <= 0)
        setSelectedName(-1);
    else
    {
        // Interpret results: each object created 4 values in the selectBuffer().
        // selectBuffer[4*i+1] is the object minimum depth value, while selectBuffer[4*i+3] is the id pushed on the stack.
        // Of all the objects that were projected in the pick region, we select the closest one (zMin comparison).
        // This code needs to be modified if you use several stack levels. See glSelectBuffer() man page.
        GLuint zMin = (selectBuffer())[1];
        setSelectedName((selectBuffer())[3]);
        for (int i=1; i<nbHits; ++i) {
//            size_t id = selectBuffer()[4*i+3];
            if ((selectBuffer())[4*i+1] < zMin)
                zMin = (selectBuffer())[4*i+1];
        }


        for (int i=1; i<nbHits; ++i) {
            size_t id = selectBuffer()[4*i+3];
            if ((selectBuffer())[4*i+1] == zMin) {
                picked_entities.push_back(id);
                setSelectedName((selectBuffer())[4*i+3]);
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::initializeGL()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glEnable( GL_LIGHTING);

    glPolygonOffset(2.0, 3.0);
    glEnable( GL_DEPTH_TEST);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_POLYGON_OFFSET_LINE);

    glEnable(GL_NORMALIZE);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    camera()->setSceneRadius(scene_radius);
    camera()->lookAt(sceneCenter() );
    GLfloat  pos[] = { 0.0, 1.0, 1.0, 0.0};
    glLightfv( GL_LIGHT0,  GL_POSITION, pos);

//  restoreStateFromFile ();

    QGLViewer::initializeGL();
    showEntireScene();
}

///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::animate()
{
    int nSize = components.size();
    for( int i = 0; i < nSize; i++)
        components[i]->animate();

    if( saveAnimation ) {
        if( !QDir("Animation").exists() ) QDir().mkdir("Animation");
        ostringstream oss;
        oss << "./Animation/snapshot";
        this->setSnapshotFileName( QString::fromStdString( oss.str() ) );
//      this->setSnapshotFormat( QString::fromStdString( ".PNG" ) );
        this->setSnapshotQuality(100);
        this->saveSnapshot();
        currAnimFrame++;
    }

    if( currAnimFrame == numAnimFrames ) {
        saveAnimation = 0;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JaalViewer::draw()
{
    glClearColor( 1.0, 1.0, 1.0, 1.0);
    glClearColor(backgroundColor[0], backgroundColor[1], backgroundColor[2], 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
    drawScene();

    boost::remove_erase_if(zoomLens, [](JMagnifyingLensPtr &l ) {
        return l->deleted;
    });

    for( size_t i = 0; i < zoomLens.size(); i++)  {
        if( zoomLens[i]->active ) zoomLens[i]->draw();
    }
}

///////////////////////////////////////////////////////////////////////////////
void JaalViewer::postDraw()
{
    QGLViewer::postDraw();
}
///////////////////////////////////////////////////////////////////////////////

void JaalViewer::drawScene()
{
    if(displayAxis) {
        glPushMatrix();
        drawAxis();
        glPopMatrix();
    }

    if( sceneFloor) {
        glPushMatrix();
        sceneFloor->draw();
        glPopMatrix();
    }

    glPushMatrix();
    int nSize = components.size();
    for( int i = 0; i < nSize; i++) {
        if( components[i]->isActive() ) components[i]->draw();
    }
    glPopMatrix();

    if( displayBox ) {
        glPushMatrix();
        drawBox(&box);
        glPopMatrix();
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////
void
JaalViewer::resetSceneBox()
{
    /*
        box[0] = 0.0;
        box[1] = 0.0;
        box[2] = 0.0;
        int nSize = components.size();
        for( int i = 0; i < nSize; i++) {
            if( components[i]->isActive() ) {
                JBoundingBox bx = component[i]->getAxixAlignedBox();
                box.getUnion(bx);
            }
        }
    */
}
///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::drawWithNames()
{
    int nSize = components.size();
    for( int i = 0; i < nSize; i++)
        components[i]->drawWithNames();
}
///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::fastDraw()
{
    glClearColor(backgroundColor[0], backgroundColor[1], backgroundColor[2], 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_ACCUM_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );

    /*
        if( lights_placement_policy  == MOVE_LIGHTS_WITH_AXIS ) {
            glPushMatrix();
            glLoadIdentity();
            JLights::getInstance().setPositions();
            glPopMatrix();
        }

        if( lights_placement_policy == MOVE_LIGHTS_WITH_CAMERA ) {
            glPushMatrix();
            JLights::getInstance().setPositions();
            glPopMatrix();
        }

    */

    int nSize = components.size();
    for( int i = 0; i < nSize; i++)
        if( components[i]->isActive() ) components[i]->fastDraw();

    QGLViewer::fastDraw();
}

///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::drawBox( const JBoundingBox *box)
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glLineWidth(2.0);

    Point3D lower = box->getLower();
    Point3D upper = box->getUpper();

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

    glColor3f(0.5, 0.5, 0.5);
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
JaalViewer::drawBox( const JHexahedronPtr &minBox)
{
    if( minBox == nullptr ) return;
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glLineWidth(2.0);
    vector<Point3D> hexpoints(8);
    for( int i = 0; i< 8; i++)
        hexpoints[i] = minBox->getNodeAt(i)->getXYZCoords();

    int v0, v1;
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    for( int i = 0; i < 12; i++) {
        JHexahedron::getEdgeTopology(i, v0, v1);
        glVertex3f(hexpoints[v0][0], hexpoints[v0][1], hexpoints[v0][2]);
        glVertex3f(hexpoints[v1][0], hexpoints[v1][1], hexpoints[v1][2]);
    }
    glEnd();
}

///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::setRotationAxis(int p)
{
    if( p == FREE_ROTATION) {
        worldConstraint->setRotationConstraintType(AxisPlaneConstraint::FREE);
        this->camera()->frame()->setConstraint(worldConstraint);
        return ;
    }

    if( p == NO_ROTATION) {
        worldConstraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
        camera()->frame()->setConstraint(worldConstraint);
        return ;
    }
    double x = 0, y = 0, z = 0;

    switch( p ) {
    case XAXIS_ROTATION:
        x = 1.0;
        break;
    case YAXIS_ROTATION:
        y = 1.0;
        break;
    case ZAXIS_ROTATION:
        z = 1.0;
        break;
    }

    Vec vc(x, y, z);
    worldConstraint->setRotationConstraintType(AxisPlaneConstraint::AXIS);
    worldConstraint->setRotationConstraintDirection(vc);
    camera()->frame()->setConstraint(worldConstraint);
}

////////////////////////////////////////////////////////////////////////////////

void JaalViewer::setTranslationAxis(int p)
{
    if( p == FREE_TRANSLATION) {
        worldConstraint->setTranslationConstraintType(AxisPlaneConstraint::FREE);
        this->camera()->frame()->setConstraint(worldConstraint);
        return ;
    }

    if( p == NO_TRANSLATION) {
        worldConstraint->setTranslationConstraintType(AxisPlaneConstraint::FORBIDDEN);
        camera()->frame()->setConstraint(worldConstraint);
        return ;
    }
    double x = 0, y = 0, z = 0;

    switch( p ) {
    case XAXIS_TRANSLATION:
        x = 1.0;
        break;
    case YAXIS_TRANSLATION:
        y = 1.0;
        break;
    case ZAXIS_TRANSLATION:
        z = 1.0;
        break;
    }

    Vec vc(x, y, z);
    worldConstraint->setTranslationConstraintType(AxisPlaneConstraint::AXIS);
    worldConstraint->setTranslationConstraintDirection(vc);
    camera()->frame()->setConstraint(worldConstraint);
}

////////////////////////////////////////////////////////////////////////////////

void
JaalViewer::resetView( const JViewDirection &v)
{
    Point3D pcenter = box.getCenter();
    Vec vc( pcenter[0], pcenter[1], pcenter[2] );
    camera()->lookAt(vc);
    scene_radius = box.getMaxLength();

    Vec pos = camera()->position();
    double dist = sqrt( pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2] );

    Vec upvec;

    if( v == JViewDirection::FRONT_VIEW) {
        vc.x = 0.0;
        vc.y = 0.0;
        vc.z = dist;
        camera()->setPosition(vc);
        upvec.x = 0.0;
        upvec.y = 1.0;
        upvec.z = 0.0;
        camera()->setUpVector(upvec);
    }

    if( v == JViewDirection::BACK_VIEW) {
        vc.x =  0.0;
        vc.y =  0.0;
        vc.z = -dist;
        camera()->setPosition(vc);
        upvec.x = 0.0;
        upvec.y = 1.0;
        upvec.z = 0.0;
        camera()->setUpVector(upvec);
    }

    if( v == JViewDirection::RIGHT_VIEW) {
        vc.x = dist;
        vc.y = 0.0;
        vc.z = 0.0;
        camera()->setPosition(vc);
        upvec.x = 0.0;
        upvec.y = 0.0;
        upvec.z = 1.0;
        camera()->setUpVector(upvec);
    }

    if( v == JViewDirection::LEFT_VIEW) {
        vc.x = -dist;
        vc.y = 0.0;
        vc.z = 0.0;
        camera()->setPosition(vc);
        upvec.x = 0.0;
        upvec.y = 0.0;
        upvec.z = 1.0;
        camera()->setUpVector(upvec);
    }

    if( v == JViewDirection::TOP_VIEW) {
        vc.x = 0.0;
        vc.y = dist;
        vc.z = 0.0;
        camera()->setPosition(vc);
        upvec.x = 0.0;
        upvec.y = 0.0;
        upvec.z = 1.0;
        camera()->setUpVector(upvec);
    }

    if( v == JViewDirection::BOTTOM_VIEW) {
        vc.x = 0.0;
        vc.y = -dist;
        vc.z = 0.0;
        camera()->setPosition(vc);
        upvec.x = 0.0;
        upvec.y = 0.0;
        upvec.z = 1.0;
        camera()->setUpVector(upvec);
    }

    scene_radius = box.getMaxLength();
    camera()->setSceneRadius(scene_radius);
    camera()->lookAt(sceneCenter() );

    updateGL();

    showEntireScene();
}

///////////////////////////////////////////////////////////////////////////////

void
JaalViewer::drawAxis()
{
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//  QGLViewer::drawAxis( sceneRadius() );
    double  d = 10*camera()->sceneRadius();
    glLineWidth(2.0);

    /*
        qglviewer::Vec from = sceneCenter();
        glEnable( GL_LINE_STIPPLE);
        glLineStipple(2, 0XAAAA);
        qglviewer::Vec to;
        to[0] = from[0] + 10*r;
        to[1] = from[1] + 0.0;
        to[2] = from[2] + 0.0;
        drawArrow(from, to);

        to[0] = from[0] + 0;
        to[1] = from[1] + 10.0*r;
        to[2] = from[2] + 0.0;
        drawArrow(from, to);

        to[0] = from[0] + 0;
        to[1] = from[1] + 0;
        to[2] = from[2] + 10*r;
        drawArrow(from, to);
    */
    glDisable(GL_LIGHTING);

    qglviewer::Vec vc = sceneCenter();
    glColor3f( 1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f( vc[0], vc[1], vc[2]);
    glVertex3f(vc[0] + d, 0.0, 0.0);
    glEnd();

    glColor3f( 0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3f( vc[0], vc[1], vc[2]);
    glVertex3f(0.0, vc[1] + d, 0.0);
    glEnd();

    glColor3f( 0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3f( vc[0], vc[1], vc[2]);
    glVertex3f(0.0, 0.0, vc[2] + d);
    glEnd();
    glDisable( GL_LINE_STIPPLE);

    glDisable(GL_LIGHTING);
}

///////////////////////////////////////////////////////////////////////////////

