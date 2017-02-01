#include "MeshTangleDialog.hpp"
///////////////////////////////////////////////////////////////////////////////

void JMeshTangleViewer :: draw()
{
    glDisable( GL_LIGHTING );
    glDisable( GL_BLEND );
    glEnable( GL_DEPTH_TEST);

    glPointSize(5.0);
    glColor3f( 0.0, 1.0, 1.0);

    int numNodes = intersectPoints.size();
    glBegin(GL_POINTS);
    for( int i = 0; i < numNodes; i++) {
        const Point2D &p0  = intersectPoints[i];
        glVertex3f( p0[0], p0[1], 0.001);
    }
    glEnd();
    glPointSize(1.0);


    int nfaces = negativeFaces.size();
    if( nfaces ) {
        glColor3f( 1.0, 0.2, 0.2);
        glBegin(GL_TRIANGLES);
        for( int i = 0; i < nfaces; i++) {
            const Point3D &p0  = negativeFaces[i]->getNodeAt(0)->getXYZCoords();
            const Point3D &p1  = negativeFaces[i]->getNodeAt(1)->getXYZCoords();
            const Point3D &p2  = negativeFaces[i]->getNodeAt(2)->getXYZCoords();
            glVertex3f( p0[0], p0[1], 0.0000);
            glVertex3f( p2[0], p2[1], 0.0000);
            glVertex3f( p1[0], p1[1], 0.0000);
        }
        glEnd();
    }

    glEnable( GL_LIGHTING );
}

/////////////////////////////////////////////////////////////////////////////////

JMeshTangleDialog :: JMeshTangleDialog( QWidget *parent) : QDialog(parent)
{
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    setupUi(this);
    makeConnections();

    numNodes2MoveLineEdit->setText( QString::number(1) );
    maxShiftLineEdit->setText( QString::number(0.1) );
    numNegativesLineEdit->setText( QString::number(0) );
    numTangledLineEdit->setText( QString::number(0 ) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshTangleDialog :: ~JMeshTangleDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

//  mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;

    mesh->getGeometry()->getCoordsArray(orgCoords,l2g);
    mesh->getTopology()->searchBoundary();
    meshtangle.setMesh(mesh);

    if( tangleViewer == nullptr ) {
        tangleViewer.reset( new JMeshTangleViewer() );
        tangleViewer->setName("MeshTangleViewer");
    }

    tangleViewer->setViewManager(viewManager);
    viewManager->attach( tangleViewer );

    viewManager->camera()->setType(Camera::ORTHOGRAPHIC);
    getTangle();
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: getOrgMesh()
{
    if( meshViewer == nullptr || mesh == nullptr) return;

    mesh->getGeometry()->setCoordsArray(orgCoords, l2g);

    negativeFaces.clear();
    tangledFaces.clear();

    clearColors();
    setNegativeColor();

    numNegativesLineEdit->setText( QString::number(0) );
    numTangledLineEdit->setText( QString::number(0 ) );

    tangleViewer->clear();
    meshtangle.clear();

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTangleDialog :: clearColors()
{
    size_t numfaces = mesh->getSize(2);
    JColor greenColor;
    greenColor[0] = 0.0;
    greenColor[1] = 1.0;
    greenColor[2] = 0.0;
    greenColor[2] = 1.0;

    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 0;
        fAttrib->color = greenColor;
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: setEdgesColor()
{
    size_t numedges = mesh->getSize(1);
    JEdgeRenderPtr eAttrib;
    JColor greenColor;
    greenColor[0] = 0.0;
    greenColor[1] = 1.0;
    greenColor[2] = 0.0;
    greenColor[3] = 0.5;

    for( size_t i = 0; i < numedges; i++) {
        JEdgePtr edge = mesh->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        eAttrib->display = 1;
        eAttrib->scale   = 1;
        eAttrib->color = greenColor;
    }

    JColor redColor;
    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 0.5;
    for( JEdgePtr edge: intersectEdges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->display = 1;
        eAttrib->scale   = 2;
        eAttrib->color = redColor;
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshTangleDialog :: getNegativeFaces()
{
    negativeFaces = meshtangle.getNegativeFaces();
    setNegativeColor();

    numNegativesLineEdit->setText( QString::number(negativeFaces.size() ) );
    tangleViewer->setNegativeFaces( negativeFaces );
}

////////////////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: setNegativeColor()
{
    if( !showNegativeElemCheckBox->isChecked() ) return;

    /*
        RenderFaceAttribute *fAttrib = nullptr;
        for( Face *face : negativeFaces) {
            face->getAttribute("Render", fAttrib);
            fAttrib->display = 1;
            fAttrib->color = MeshEntityColor::getRandomColor();
        }
    */
}
////////////////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: setOverlapColor()
{
    if( !showOverlapElemCheckBox->isChecked() ) return;

    JColor yellowColor;
    yellowColor[0] = 1.0;
    yellowColor[1] = 0.5;
    yellowColor[2] = 0.0;
    yellowColor[3] = 0.5;

    JFaceRenderPtr fAttrib;
    for( JFacePtr face : tangledFaces) {
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
        fAttrib->color = yellowColor;
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTangleDialog :: getOverlappedFaces()
{
    tangledFaces = meshtangle.getOverlapFaces();
    numTangledLineEdit->setText( QString::number(tangledFaces.size() ) );

    intersectEdges = meshtangle.getIntersectEdges();
    numCrossingEdgesLineEdit->setText( QString::number(intersectEdges.size() ) );

    vector<Point2D> points = meshtangle.getIntersectPoints();
    numEdgePointsLineEdit->setText( QString::number(points.size() ) );

    if( tangleViewer )
        tangleViewer->setIntersectPoints( points );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: moveNodes()
{
    if( mesh == nullptr ) return;

    QString qstr;
    qstr = numNodes2MoveLineEdit->text() ;
    int  numNodesMove = qstr.toInt();

    qstr = maxShiftLineEdit->text() ;
    double maxShift = qstr.toDouble();
    maxShift = max(1.0E-10, maxShift);

    double xshift, yshift;
    for( int i = 0; i < numNodesMove; i++) {
//       Vertex *vtx = mesh->getRandomNode();
        JNodePtr vtx = mesh->getNodeAt(10);
        if( vtx ) {
            Point3D xyz = vtx->getXYZCoords();
//            xshift = JMath::random_value(-maxShift, maxShift);
//            yshift = JMath::random_value(-maxShift, maxShift);
            xshift  =0.28;
            yshift  =0.0;
            xyz[0] += xshift;
            xyz[1] += yshift;
            vtx->setXYZCoords(xyz);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: getTangle()
{
    clearColors();

    meshtangle.searchOverlap();

    // Classify faces into positive or negative ...
    getNegativeFaces();

    // Check the faces which overlap with the negative elements...
    getOverlappedFaces();

//    setColor();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTangleDialog :: setColor()
{
    clearColors();
    setOverlapColor();
    setEdgesColor();
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTangleDialog :: makeConnections()
{
    PushButton( getOrgMeshPushButton,  [=] {getOrgMesh(); });
    PushButton( moveNodePushButton,  [=] {moveNodes();});
    PushButton( getTangledPushButton,  [=] {getTangle();});

    CheckBox( showNegativeElemCheckBox, [=] { setColor();});
    CheckBox( showOverlapElemCheckBox,  [=] { setColor();});

    PushButton( closePushButton,  [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
