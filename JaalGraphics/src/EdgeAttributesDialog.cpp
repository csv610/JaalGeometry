#include "EdgeAttributesDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JEdgeAttributesDialog :: JEdgeAttributesDialog( QWidget *parent) : QDialog(parent)
{
    viewManager = nullptr;
    meshViewer  = nullptr;
    picker      = nullptr;
    update_default_values = 0;

    setupUi(this);

    widthLineEdit->setText( QString::number(1.0) );
    numCylSidesLineEdit->setText( QString::number(20) );

    boundaryTagLineEdit->setText( QString::number(1) );

    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JEdgeAttributesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

    JEdgeDraw *drawEdge = meshViewer->getEdgeDraw();
    if( drawEdge == nullptr ) return;

    float radius = drawEdge->getCylinderRadius();
    cylRadiusLineEdit->setText( QString::number(radius) );

    viewManager->attach( this );
    viewManager->setMouseTracking(1);

    /*
        picker = meshViewer->getEntityPicker();
        if( picker ) {
            picker->setPickableEntity(1);
            picker->setMode(1);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString::fromStdString(name));
    JEdgeSequence medges = mesh->getEdges();
    setEdges(medges);
}
///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: setEdges(const JEdgeSequence &seq)
{
    edges = seq;

    int numVis = 0;
    JEdgeRenderPtr attrib;
    for( size_t i = 0; i < edges.size(); i++) {
        edges[i]->getAttribute("Render", attrib);
        if(attrib->display) numVis++;
    }

    numEdgesLineEdit->setText( QString::number(edges.size()) );
    numVisibleLineEdit->setText( QString::number(numVis) );

    bool val = numVis == 0 ? 0:1;
    displayCheckBox->setChecked(val);
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: setColor()
{
    if( meshViewer == nullptr) return;

    QColor color = QColorDialog::getColor();
    float rgb[3];
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    rgb[3] = 1.0;

    JEdgeRenderPtr eAttrib;
    for( size_t i = 0; i < edges.size(); i++) {
        edges[i]->getAttribute("Render", eAttrib);
        eAttrib->color[0] = rgb[0];
        eAttrib->color[1] = rgb[1];
        eAttrib->color[2] = rgb[2];
        eAttrib->display  = 1;
    }
    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////
void JEdgeAttributesDialog :: setAlpha()
{
    if( meshViewer == nullptr) return;
    double alpha = transparencySpinBox->value();

    JEdgeRenderPtr eAttrib;
    for( size_t i = 0; i < edges.size(); i++)  {
        edges[i]->getAttribute("Render", eAttrib);
        eAttrib->color[3] = alpha;
    }
    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: setGlyph()
{
    if( meshViewer ==  nullptr) return;

    int glyph = 0;
    if( lineRadioButton->isChecked() ) {
        glyph = JEdgeDraw::EDGE_AS_LINE;
    }

    if( cylinderRadioButton->isChecked() )
        glyph = JEdgeDraw::EDGE_AS_CYLINDER;

    size_t numedges = edges.size();
    JEdgeRenderPtr  eAttrib;
    for( size_t i = 0; i < numedges; i++) {
        edges[i]->getAttribute("Render", eAttrib);
        if( eAttrib) eAttrib->glyph  = glyph;
    }

    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: setLineWidth()
{
    if( meshViewer ==  nullptr) return;
    JEdgeDraw *drawEdge = meshViewer->getEdgeDraw();

    QString str = widthLineEdit->text() ;
    float  radius = str.toDouble();

    if( update_default_values) drawEdge->setLineWidth(radius);

    JEdgeRenderPtr eAttrib;
    if( lineRadioButton->isChecked() ) {
        size_t numedges = edges.size();
        for( size_t i = 0; i < numedges; i++) {
            int err = edges[i]->getAttribute("Render", eAttrib);
            if( !err) {
                eAttrib->lineWidth = radius;
                eAttrib->scale     = 1.0;
            }
        }
        meshViewer->updateBuffers( mesh );
    }
}
///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: setCylinderRadius()
{
    if( meshViewer ==  nullptr ) return;
    JEdgeDraw *drawEdge = meshViewer->getEdgeDraw();
    if( drawEdge == nullptr ) return;

    QString str = cylRadiusLineEdit->text() ;
    double radius  = str.toDouble();

    if( update_default_values ) drawEdge->setCylinderRadius( radius );

    JEdgeRenderPtr eAttrib;
    if( cylinderRadioButton->isChecked() ) {
        size_t numedges = edges.size();
        for( size_t i = 0; i < numedges; i++) {
            int err = edges[i]->getAttribute("Render", eAttrib);
            if( !err ) eAttrib->cylinderRadius = radius;
        }
        meshViewer->updateBuffers( mesh );
    }
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: setCylinderSides()
{
    if( meshViewer ==  nullptr ) return;
    JEdgeDraw *drawEdge = meshViewer->getEdgeDraw();

    if( drawEdge == nullptr ) return;

    QString str = numCylSidesLineEdit->text() ;
    int n  = str.toInt();

    drawEdge->setNumCylinderSides( n );

    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////
void JEdgeAttributesDialog :: checkDisplay()
{
    if( meshViewer ==  nullptr ) return;

    int numVis = 0;
    bool val = displayCheckBox->isChecked();

    JEdgeRenderPtr nAttrib;
    for( size_t i = 0; i < edges.size(); i++) {
        edges[i]->getAttribute("Render", nAttrib);
        nAttrib->display = val;
        if( val ) numVis++;
    }

    numVisibleLineEdit->setText( QString::number(numVis) );
    displayCheckBox->setChecked(val);
    meshViewer->updateBuffers( mesh );

    /*
        cout << "Do not display these edges " << endl;
        for( size_t i = 0; i < mesh->getSize(1); i++) {
            const JEdgePtr &e = mesh->getEdgeAt(i);
            if( e->isActive() ) {
            e->getAttribute("Render", nAttrib);
            if( !nAttrib->display)
            cout << e->getID() << " NodeID " << e->getNodeAt(0)->getID() << "  " << e->getNodeAt(1)->getID() << endl;
            }
        }
    */
}
///////////////////////////////////////////////////////////////////////////////
void JEdgeAttributesDialog :: mousePressEvent(QMouseEvent *e)
{
    if( meshViewer ==  nullptr ) return;
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JEdgeAttributesDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( meshViewer ==  nullptr ) return;

    if( boundaryTagCheckBox->isChecked() ) {
        if( picker ) {
            QString str = boundaryTagLineEdit->text() ;
            int bmark  = str.toInt();
            JEdgeSequence edgeSeq = picker->getPickedEdges();
            if( !edgeSeq.empty() ) {
                edgeSeq[0]->setAttribute("Boundary", bmark);
            }
            meshViewer->refreshDisplay();
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

void JEdgeAttributesDialog :: closeDialog()
{
    edges.clear();
    this->parentWidget()->show();
    this->hide();

    /*
        if( picker ) {
            picker->clearAll();
            picker->setPickableEntity(-1);
        }
    */

    if( viewManager) {
        viewManager->detach( this );
        viewManager->setMouseTracking(0);
    }

    picker = nullptr;
}

///////////////////////////////////////////////////////////////////////////////
void JEdgeAttributesDialog :: makeConnections()
{
    connect(colorPushButton, SIGNAL( clicked() ), this, SLOT( setColor() ));
    SpinBoxd( transparencySpinBox, [=] {setAlpha();});
    RadioButton( lineRadioButton,  [=] {setGlyph();});
    RadioButton( cylinderRadioButton, [=] {setGlyph();});

    LineEdit( widthLineEdit,  [=] {setLineWidth();});
    LineEdit( cylRadiusLineEdit,  [=] {setCylinderRadius();});
    LineEdit( numCylSidesLineEdit,  [=] {setCylinderSides();});
    CheckBox( displayCheckBox, [=] {checkDisplay();});
    PushButton( closePushButton,  [=] {closeDialog();});
}
///////////////////////////////////////////////////////////////////////////////

/*
void JEdgeAttributesDialog :: setOffset()
{
    if( meshViewer ==  nullptr ) return;
    DrawEdge *drawEdge = meshViewer->getDrawEdge();

    if( drawEdge == nullptr ) return;

    bool val = offsetCheckBox->isChecked();
    drawEdge->setOffset(val);

    QString str;
    if( val )
    {
	str = factorLineEdit->text() ;
	float factor  = str.toDouble();

	str = unitLineEdit->text() ;
	float unit  = str.toDouble();
	drawEdge->setOffset(factor, unit);
    }

    meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////
