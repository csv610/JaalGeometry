#include "CurveGenDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JCurveGenDialog :: JCurveGenDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    nCount    = 0;
    lastVertex = nullptr;

    discretizeLineEdit->setText( QString::number(1) );
    xCoordLineEdit->setText( QString::number(0.0) );
    yCoordLineEdit->setText( QString::number(0.0) );
    zCoordLineEdit->setText( QString::number(0.0) );
    mergeDistLineEdit->setText( QString::number(0.0) );
    nodeIDLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////

JCurveGenDialog :: ~JCurveGenDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    viewManager->attach(this);
    viewManager->setRotationAxis(JaalViewer::NO_ROTATION);
}

///////////////////////////////////////////////////////////////////////////////
void JCurveGenDialog :: startCurve()
{
    nCount = 0;
    mesh = JMesh::newObject();
    meshViewer->addObject(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: mouseMoveEvent( QMouseEvent *event)
{
    if( canvas == nullptr) return;

    if(!enableSketchCheckBox->isChecked() )  {
        viewManager->freezeView(0);
        return;
    }

    viewManager->setRotationAxis(JaalViewer::NO_ROTATION);

    if( mesh == nullptr) startCurve();

    JColor color = JEntityColor::getColor("Blue");

    Point3D xyz;
    int err = viewManager->getMouseCurrentXYZPosition(xyz);
    if( !err) {
        JNodePtr vnew = JNode::newObject();
        vnew->setID(nCount++);
        xyz[2] = 0.0;   // Always on the XY plane ...
        vnew->setXYZCoords(xyz);
        mesh->addObject(vnew);
        if( nCount > 1) {
            const JNodePtr &vlast = mesh->getNodeAt(nCount-2);
            JEdgePtr newedge = JEdge::newObject( vlast, vnew);
            mesh->addObject(newedge);
        }
        JNodeRenderPtr attrib = JNodeRender::newObject();
        attrib->color     = color;
        attrib->pointSize = 5.0;
        vnew->setAttribute("Render", attrib);
        meshViewer->updateBuffers(mesh);
        numNodesLineEdit->setText( QString::number(nCount));
    }
}

///////////////////////////////////////////////////////////////////////////////
void JCurveGenDialog :: setCanvasColor()
{
    if( canvas == nullptr ) return;

    JColor rgb;
    QColor color = QColorDialog::getColor();
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    rgb[3] = 1.0;

    size_t numfaces = canvas->getSize(2);
    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = canvas->getFaceAt(i);
        f->getAttribute("Render", fAttrib);
        fAttrib->color = rgb;
    }
    meshViewer->updateBuffers(canvas);
}

///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: startSketching()
{
    nCount = 0;
    if(enableSketchCheckBox->isChecked()  && canvas == nullptr)  {
        createCanvas();
        return;
    }

    if(!enableSketchCheckBox->isChecked() )  {
        meshViewer->removeObject( canvas );
        canvas.reset();
        meshViewer->refreshDisplay();
        return;
    }
}
///////////////////////////////////////////////////////////////////////////////
void JCurveGenDialog :: createCanvas()
{
    if( canvas ) meshViewer->removeObject(canvas);

    int dim[2];
    int ngrid = gridSizeSpinBox->value();
    dim[0] = ngrid;
    dim[1] = ngrid;

    double len[2];
    int  clen  = canvasSizeSpinBox->value();
    len[0] = clen;
    len[1] = clen;

    double  org[2];
    org[0] = -0.5*len[0];
    org[1] = -0.5*len[1];
    canvas = AllQuadMeshGenerator::getStructuredMesh(dim, len, org);
    JMeshAffineTransform maffine(canvas);
    maffine.translate(0.0, 0.0, -0.00001);
    meshViewer->addObject(canvas);
}

///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: openKnotsDialog()
{
    if( knotsDialog == nullptr )
        knotsDialog.reset(new JKnotsDialog(this));
    knotsDialog->setViewManager( viewManager );
    newCurve = knotsDialog->getKnot();
    knotsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: getKnots()
{
    double dist;
    if( knotsDialog )
    {
        newCurve = knotsDialog->getKnot();
        dist = newCurve->getMinimumLocalDistance();
        minLocalDistLineEdit->setText( QString::number(dist) );

        dist = newCurve->getMinimumGlobalDistance();
        minGlobalDistLineEdit->setText( QString::number(dist) );

        dist = newCurve->getLength();
        lengthLineEdit->setText( QString::number(dist) );

        int ncount = newCurve->getSize(0);
        nodeIDLineEdit->setText( QString::number(nCount) );
    }
}
///////////////////////////////////////////////////////////////////////////////
void JCurveGenDialog :: closeCurve()
{
    size_t numnodes = mesh->getSize(0);
    if( numnodes < 2) return;

    JNodePtr v0 = mesh->getNodeAt(0);
    JNodePtr v1 = mesh->getNodeAt(numnodes-1);
    JEdgePtr edge = JEdge::newObject(v1,v0);
    mesh->addObject( edge);
    meshViewer->updateBuffers(mesh);

    // You are also done with this curve... Start with a new one.
    mesh.reset();
}
///////////////////////////////////////////////////////////////////////////////
void JCurveGenDialog :: deleteCurve()
{
    if( mesh == nullptr) return;

    mesh->deleteAll();
    meshViewer->updateBuffers(mesh);
    mesh.reset();
    nCount = 0;
    numNodesLineEdit->setText( QString::number(nCount));
}

///////////////////////////////////////////////////////////////////////////////
void JCurveGenDialog :: deletePoint()
{
    if( mesh == nullptr) return;

    size_t numedges = mesh->getSize(1);
    if( numedges ) {
        JEdgePtr edge = mesh->getEdgeAt(numedges-1);
        edge->setStatus(JMeshEntity::REMOVE);
    }

    size_t numnodes = mesh->getSize(0);
    if( numnodes) {
        JNodePtr node = mesh->getNodeAt(numnodes-1);
        node->setStatus(JMeshEntity::REMOVE);
    }
    mesh->pruneAll();
    meshViewer->updateBuffers(mesh);

    nCount = numnodes-1;
    numNodesLineEdit->setText( QString::number(nCount));
}
///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: closeDialog()
{
    if( canvas) {
        meshViewer->removeObject(canvas);
        canvas.reset();
    }

    viewManager->setRotationAxis(JaalViewer::FREE_ROTATION);
    viewManager->detach(this);
    this->parentWidget()->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: refineSegments()
{
    if( newCurve )
    {
        QString qstr = zCoordLineEdit->text();
        int n =  qstr.toInt();
        if( n > 2) newCurve->refineSegments(n);
    }
}
///////////////////////////////////////////////////////////////////////////////

void JCurveGenDialog :: makeConnections()
{
    PushButton( knotsPushButton,  [=] { openKnotsDialog();});
    PushButton( discretizeSegmentPushButton,  [=] { refineSegments();});

    CheckBox( enableSketchCheckBox, [=] {startSketching(); });

    PushButton( newCurvePushButton, [=] { startCurve();});
    PushButton( closeCurvePushButton, [=] {closeCurve();});
    PushButton( deletePointPushButton, [=] { deletePoint(); });
    PushButton( deleteCurvePushButton, [=] { deleteCurve(); });
    PushButton( canvasColorPushButton, [=] { setCanvasColor(); });
    SpinBoxi(gridSizeSpinBox,   [=] {createCanvas(); });
    SpinBoxi(canvasSizeSpinBox, [=] {createCanvas(); });
    connect( canvasColorPushButton, SIGNAL( clicked() ), this, SLOT( setCanvasColor() ));

    PushButton( closePushButton, [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
/*
void JCurveGenDialog :: addNode()
{
     double  dist;
     Point3D p3d;

     QString qstr;
     qstr = xCoordLineEdit->text();
     p3d[0] = qstr.toDouble();

     qstr = yCoordLineEdit->text();
     p3d[1] = qstr.toDouble();

     qstr = zCoordLineEdit->text();
     p3d[2] = qstr.toDouble();

     if( lastVertex ) {
          dist = JMath::length(lastVertex->getXYZCoords(), p3d );
          if( dist < 1.0E-10)  {
               QMessageBox msg;
               msg.setIcon(QMessageBox::Warning);
               msg.setText("Distance too small with previous node");
               msg.setInformativeText("Do you still want to add it ? ");
               msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel);
               int ret = msg.exec();
               if( ret == QMessageBox::Cancel ) return;
          }
     }

     Vertex *vtx = Vertex::newObject();
     vtx->setXYZCoords( p3d );

     bool val = 1;
     vtx->setAttribute( "Display", val );
     newCurve->addControlNode( vtx );
     nCount++;

     lastVertex = vtx;

     JEdgeSequence edges;
     newCurve->getEdges(edges);
     for( size_t i = 0; i < edges.size(); i++)
          edges[i]->setAttribute( "Display", val );

     dist = newCurve->getMinimumLocalDistance();
     minLocalDistLineEdit->setText( QString::number(dist) );

     dist = newCurve->getMinimumGlobalDistance();
     minGlobalDistLineEdit->setText( QString::number(dist) );

     nodeIDLineEdit->setText( QString::number(nCount) );
     meshViewer->refreshDisplay();
}
*/
