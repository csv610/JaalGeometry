#include "NodeAttributesDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int JNodeDegreeColor :: assign(const JNodePtr &vertex)
{
    if( vertex == nullptr ) return 1;

    JNodeRenderPtr nAttrib;
    vertex->getAttribute("Render", nAttrib);

    color[3] = alpha;

    int nSize = vertex->getNumRelations(0);
    if( nSize < lowDegree ) {
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        nAttrib->color = color;
        return 0;
    }

    if( nSize == equalDegree ) {
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
        nAttrib->color = color;
        return 0;
    }

    if( nSize >  highDegree ) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        nAttrib->color = color;
        return 0;
    }

    return 1;
}
///////////////////////////////////////////////////////////////////////////////
JNodeAttributesDialog :: JNodeAttributesDialog( QWidget *parent) : QDialog(parent)
{
    meshViewer = nullptr;
    viewManager = nullptr;

    setupUi(this);
    makeConnections();

    numSlicesLineEdit->setText( QString::number(16) );
    numStacksLineEdit->setText( QString::number(16) );
}

///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JNodeAttributesDialog :: setMesh(const JMeshPtr &m)
{
    mesh = m;
    string name;
    if( mesh) name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}
///////////////////////////////////////////////////////////////////////////////
void JNodeAttributesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    JNodeDraw *drawNode = meshViewer->getNodeDraw();
    if( drawNode == nullptr ) return;

    double radius;
    radius = drawNode->getPointSize();
    pointSizeSpinBox->setValue( radius );

    radius = drawNode->getBallRadius();
    sphereRadiusLineEdit->setText( QString::number(radius) );
    update_default_values = 0;
}

///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: setNodes(JNodeSequence &seq)
{
    nodes = seq;

    size_t numVis   = 0;
    size_t numNodes = nodes.size();

    JNodeRenderPtr attrib;
    for( size_t i = 0; i < numNodes; i++) {
        nodes[i]->getAttribute("Render", attrib);
        if( attrib->display) numVis++;
    }

    numNodesLineEdit->setText( QString::number(numNodes) );
    numVisibleLineEdit->setText( QString::number(numVis) );

    bool val = 0;
    if( numVis ) val = 1;
    displayCheckBox->setChecked(val);
    meshViewer->updateBuffers();
}

///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: setColor()
{
    if( meshViewer == nullptr ) return;

    if( nodes.empty() ) return;

    QColor color = QColorDialog::getColor();

    float rgb[3];
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;

    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < nodes.size(); i++) {
        nodes[i]->getAttribute("Render", nAttrib);
        nAttrib->color[0] = rgb[0];
        nAttrib->color[1] = rgb[1];
        nAttrib->color[2] = rgb[2];
    }

    meshViewer->updateBuffers( mesh );
}
///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: setBorderColor()
{
    if( meshViewer == nullptr ) return;

    QColor color = QColorDialog::getColor();
    JColor rgba;
    rgba[0] = color.red()/255.0;
    rgba[1] = color.green()/255.0;
    rgba[2] = color.blue()/255.0;
    rgba[3] = 1.0;
    meshViewer->getNodeDraw()->setBorderColor(rgba);
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: setGlyph()
{
    if( meshViewer ==  nullptr ) return;

    if( nodes.empty() ) return;

    JNodeRenderPtr nAttrib;

    int glyph;
    double radius;
    if( pointRadioButton->isChecked() ) {
        glyph  = JNodeRender::NODE_AS_POINT;
        radius = meshViewer->getNodeDraw()->getPointSize();
        if( update_default_values)
            meshViewer->getNodeDraw()->setGlyph(glyph);
        for( size_t i = 0; i < nodes.size(); i++) {
            nodes[i]->getAttribute("Render", nAttrib);
            nAttrib->glyph  = glyph;
            nAttrib->pointSize = radius;
        }
    }

    if( sphereRadioButton->isChecked() ) {
        glyph = JNodeRender::NODE_AS_SPHERE;
        radius = meshViewer->getNodeDraw()->getBallRadius();
        if( update_default_values)
            meshViewer->getNodeDraw()->setGlyph(glyph);
        for( size_t i = 0; i < nodes.size(); i++) {
            nodes[i]->getAttribute("Render", nAttrib);
            nAttrib->glyph  = glyph;
            nAttrib->ballRadius = radius;
        }
    }

    if( splatRadioButton->isChecked() ) {
        glyph = JNodeRender::NODE_AS_SPLAT;
        radius = meshViewer->getNodeDraw()->getBallRadius();
        if( update_default_values)
            meshViewer->getNodeDraw()->setGlyph(glyph);
        for( size_t i = 0; i < nodes.size(); i++) {
            nodes[i]->getAttribute("Render", nAttrib);
            nAttrib->glyph  = glyph;
        }
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: setPointSize()
{
    if( meshViewer ==  nullptr ) return;

    if( nodes.empty() ) return;

    float psize = (float) pointSizeSpinBox->value();

    meshViewer->getNodeDraw()->setPointSize( psize );

    if( pointRadioButton->isChecked() ) {
        JNodeRenderPtr nAttrib;
        for( size_t i = 0; i < nodes.size(); i++) {
            nodes[i]->getAttribute("Render", nAttrib);
            if(nAttrib) nAttrib->pointSize = psize;
        }
    }

    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////
void JNodeAttributesDialog :: setBorderThickness()
{
    float thickness = borderThickSpinBox->value();
    meshViewer->getNodeDraw()->setBorderThickness(thickness );
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: setSphereRadius()
{
    if( meshViewer ==  nullptr ) return;

    if( nodes.empty() ) return;

    QString str = sphereRadiusLineEdit->text() ;
    double radius  = str.toDouble();

    meshViewer->getNodeDraw()->setBallRadius( radius );

    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < nodes.size(); i++) {
        int err = nodes[i]->getAttribute("Render", nAttrib);
        if( !err ) nAttrib->ballRadius = radius;
    }
    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////
void JNodeAttributesDialog :: closeDialog()
{
    nodes.clear();
    this->parentWidget()->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JNodeAttributesDialog :: setSphResolution()
{
    if( meshViewer ==  nullptr ) return;

    JNodeDraw *drawNode = meshViewer->getNodeDraw();
    if( drawNode == nullptr ) return;

    QString str = numSlicesLineEdit->text() ;
    int nSlices  = str.toInt();

    str = numStacksLineEdit->text() ;
    int nStacks  = str.toInt();

    drawNode->setSphereResolution( nSlices, nStacks);
}

///////////////////////////////////////////////////////////////////////////////
void JNodeAttributesDialog :: checkDisplay()
{
    bool val = displayCheckBox->isChecked();
    int numVis = 0;
    JNodeRenderPtr nAttrib;

    for( size_t i = 0; i < nodes.size(); i++) {
        nodes[i]->getAttribute("Render", nAttrib);
        nAttrib->display = val;
        if(val) numVis++;
    }

    numVisibleLineEdit->setText( QString::number(numVis) );
    meshViewer->updateBuffers();
    emit nodesChanged();
}
///////////////////////////////////////////////////////////////////////////////

void JNodeAttributesDialog :: makeConnections()
{
    PushButton( pointColorPushButton,  [=] {setColor();});
    PushButton( borderColorPushButton, [=] {setBorderColor();});

    LineEdit( sphereRadiusLineEdit,  [=] {setSphereRadius(); });
    LineEdit( numSlicesLineEdit,     [=] {setSphResolution();});
    LineEdit( numStacksLineEdit,     [=] {setSphResolution();});

    RadioButton( pointRadioButton,   [=] { setGlyph();});
    RadioButton( splatRadioButton,   [=] { setGlyph();});
    RadioButton( sphereRadioButton,  [=] { setGlyph();});

    SpinBoxi( pointSizeSpinBox,   [=] {setPointSize();});
    SpinBoxd( borderThickSpinBox, [=] {setBorderThickness();});

    CheckBox( displayCheckBox,  [=] { checkDisplay();});
    PushButton( closePushButton,   [=] { closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////

/*
void JNodeAttributesDialog :: setOffset()
{
    if( meshViewer ==  nullptr ) return;
    NodeDraw *drawNode = meshViewer->getNodeDraw();

    if( drawNode == nullptr ) return;

    bool val = offsetCheckBox->isChecked();
    drawNode->setOffset(val);
    QString str;
    if( val )
    {
        str = factorLineEdit->text() ;
        float factor  = str.toDouble();

        str = unitLineEdit->text() ;
        float unit  = str.toDouble();
        drawNode->setOffset(factor, unit);
    }
}
*/

