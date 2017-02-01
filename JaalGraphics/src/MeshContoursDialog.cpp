#include "MeshContoursDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshContoursDialog :: JMeshContoursDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    cornerAngleSpinBox->setValue(90.0);
}

///////////////////////////////////////////////////////////////////////////////

JMeshContoursDialog :: ~JMeshContoursDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshContoursDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");

    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    setMesh( meshViewer->getCurrentMesh() );

    viewManager->attach( this );
    picker = meshViewer->getEntityPicker();
    if( picker) picker->setMode(1);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshContoursDialog :: setMesh( const JMeshPtr &m)
{
    boundedges.clear();
    mesh = m;
    if( mesh == nullptr ) return ;

    mesh->getGeometry()->getCoordsArray(orgCoords,l2g);

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    if( !mesh->getTopology()->isClosed() ) {
        mesh->getTopology()->getBoundary(boundedges);
    }
    numContoursLineEdit->setText( QString::number( boundedges.size() ) );
    currContourSpinBox->setMaximum(boundedges.size()-1);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickable = 0;

    displayContour();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshContoursDialog :: keyPressEvent( QKeyEvent *e)
{
    // After selecting the vertex, if the user presses "P" Key then the
    // vertex becomes constraints vertex and changes color to red..
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshContoursDialog :: mousePressEvent(QMouseEvent *e)
{
    if( e->button() == Qt::LeftButton && e->modifiers() == Qt::ShiftModifier )
        if( meshViewer) meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshContoursDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( picker == nullptr) return;

    JNodeSequence nodeSeq = picker->getPickedNodes();
    if( nodeSeq.empty() ) return;


    int id = 1;
    if( e->button() == Qt::LeftButton && e->modifiers() == Qt::ShiftModifier ) {
        JColor red  = JEntityColor::getColor("Red");
        JColor blue = JEntityColor::getColor("Blue");
        JNodeRenderPtr attrib;
        JNodePtr picked = nodeSeq[0];
        picked->getAttribute("Render", attrib);
        if( picked->hasAttribute("Constraint") ) {
            picked->deleteAttribute("Constraint");
            attrib->color = blue;
        }  else {
            picked->setAttribute("Constraint", id);
            attrib->color = red;
        }
        meshViewer->updateBuffers(mesh);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshContoursDialog :: displayContour()
{
    int nSize = boundedges.size();
    if( nSize < 1) return;

    int currID = currContourSpinBox->value();

    JColor red = JEntityColor::getColor("Red");
    JColor blue = JEntityColor::getColor("Blue");

    JEdgeRenderPtr eAttrib;
    JNodeRenderPtr vAttrib;
    JNodeSequence   bnodes;

    for( size_t i = 0; i < nSize; i++) {
        for( const JEdgePtr &edge : boundedges[i]) {
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = 1;
            eAttrib->scale   = 2;
            eAttrib->color   = red;
        }
        JMeshTopology::getEntitySet( boundedges[i], bnodes);
        for( const JNodePtr &vtx : bnodes) {
            vtx->getAttribute("Render", vAttrib);
            vAttrib->display = 1;
            vAttrib->glyph   = 0;
        }
    }

    JMeshTopology::getEntitySet( boundedges[currID], bnodes);
    for( const JNodePtr &vtx : bnodes) {
        vtx->getAttribute("Render", vAttrib);
        vAttrib->display = 1;
        vAttrib->glyph   = 1;
        vAttrib->color   = blue;
    }

    JMeshContour mc;
    mc.setSegments( boundedges[currID] );

    if( cornersCheckBox->isChecked() ) {
        double theta = cornerAngleSpinBox->value();
        bnodes = mc.getCorners( theta);
        for( const JNodePtr &vtx : bnodes) {
            vtx->getAttribute("Render", vAttrib);
            vAttrib->display = 1;
            vAttrib->glyph   = 1;
            vAttrib->color   = red;
        }
        numCornersLineEdit->setText(QString::number(bnodes.size()));
    }
    meshViewer->updateBuffers(mesh);

    double len = mc.getLength();
    lengthLineEdit->setText(QString::number(len));

    numSegmentsLineEdit->setText(QString::number(boundedges[currID].size()));
    reparamEdgesLineEdit->setText(QString::number(boundedges[currID].size()));
}

///////////////////////////////////////////////////////////////////////////////
void JMeshContoursDialog :: getSmooth()
{
    int id = currContourSpinBox->value();
    if( id >= boundedges.size() ) return;

    JMeshContour polyline;
    polyline.setSegments( boundedges[id] );
    polyline.smooth();
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshContoursDialog :: getReparam()
{
    int id = currContourSpinBox->value();
    if( id >= boundedges.size() ) return;

    JMeshContour polyline;
    polyline.setSegments( boundedges[id] );

    QString qstr = reparamEdgesLineEdit->text();
    int N = qstr.toInt();
    JEdgeSequence newedges = polyline.getDiscretized(N);

    mesh->setActiveBit(0);
    JMeshPtr newmesh = JMesh::newObject(newedges);
    meshViewer->addObject(newmesh);
    boundedges[id] = newedges;
    displayContour();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshContoursDialog :: getOriginal()
{
    if( mesh == nullptr) return;
    mesh->getGeometry()->setCoordsArray(orgCoords,l2g);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshContoursDialog :: getSimplify()
{
    if( mesh == nullptr) return;

    if( mesh->getTopology()->getDimension() != 1) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: Mesh has internal elements: Do you want to delete them ");
        msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel);
        int ret = msg.exec();
        if( ret == QMessageBox::Cancel) return;
    }
    mesh->deleteCells();
    mesh->deleteFaces( JMeshEntity::INTERNAL_ENTITY);
    mesh->deleteEdges( JMeshEntity::INTERNAL_ENTITY);

    int id = currContourSpinBox->value();
    if( id >= boundedges.size() ) return;

    /*
        QString qstr = simplifyComboBox->currentText();
        string  str = StdString(qstr);
        int     algo = JMeshContour::getSimplifyAlgorithmID(str);

        qstr = toleranceLineEdit->text();
        double tol = qstr.toDouble();

        JMeshContour polyline;
        polyline.setSegments( boundedges[id] );

        JEdgeSequence newedges =  polyline.getSimplify(algo, tol);
    */

}
///////////////////////////////////////////////////////////////////////////////
void JMeshContoursDialog :: closeDialog()
{
    if( viewManager ) viewManager->detach( this );
    this->close();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshContoursDialog :: makeConnections()
{
    connect( cornersCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( displayContour() ));
    connect( cornerAngleSpinBox, SIGNAL( valueChanged(double) ), this, SLOT( displayContour() ));
    connect( reparamPushButton, SIGNAL( clicked() ), this, SLOT( getReparam() ));
    connect( simplifyPushButton, SIGNAL( clicked() ), this, SLOT( getSimplify() ));

    connect( originalPushButton, SIGNAL( clicked() ), this, SLOT( getOriginal() ));
    connect( smoothPushButton, SIGNAL( clicked() ), this, SLOT( getSmooth() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
