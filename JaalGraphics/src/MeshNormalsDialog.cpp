#include "MeshNormalsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshNormalsDialog :: JMeshNormalsDialog( QWidget *parent) : QDialog(parent)
{
    viewManager = nullptr;
    meshViewer  = nullptr;

    setupUi(this);

    makeConnections();

    mesh = nullptr;
    lengthLineEdit->setText( QString::number(0.1) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshNormalsDialog :: ~JMeshNormalsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNormalsDialog :: setMesh( const JMeshPtr &m, int pos)
{
    mesh = m;
    if( mesh == nullptr) return;
    normalAt = pos;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ));

}
///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: setLength()
{
    if( meshViewer == nullptr) return;
    QString str = lengthLineEdit->text();

    double len = fabs(str.toDouble())+ 1.0E-06;

    if( normalAt == 0)
        meshViewer->getNodeDraw()->setNormalsLength(len);

    if( normalAt == 2)
        meshViewer->getFaceDraw()->setNormalsLength(len);

    meshViewer->updateBuffers();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: setColor()
{
    QColor color = QColorDialog::getColor();

    JColor rgba;
    rgba[0] = color.red()/255.0;
    rgba[1] = color.green()/255.0;
    rgba[2] = color.blue()/255.0;
    rgba[3] = 0.0;

    if( meshViewer == nullptr) return;

    if( normalAt == 0) meshViewer->getNodeDraw()->setNormalsColor( rgba );
    if( normalAt == 2) meshViewer->getFaceDraw()->setNormalsColor( rgba );

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: init()
{
    if( viewManager == nullptr) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr) return;
    if( mesh == nullptr) return ;

    double len = 0.0;
    if( normalAt == 0)
        len = meshViewer->getNodeDraw()->getNormalsLength();

    if( normalAt == 2)
        len = meshViewer->getFaceDraw()->getNormalsLength();

    lengthLineEdit->setText( QString::number(len) );

    checkDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: checkDisplay()
{
    if( meshViewer == nullptr) return;

    bool val = displayCheckBox->isChecked();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayNormals[normalAt] = val;

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: setScale()
{
    if( meshViewer == nullptr) return;
    QString str = lengthLineEdit->text();
    double len = str.toDouble();

    if( normalAt == 0)
        meshViewer->getNodeDraw()->setNormalsLength(len);

    if( normalAt == 2)
        meshViewer->getFaceDraw()->setNormalsLength(len);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: getConsistent()
{
    if( mesh == nullptr) return;
    mesh->getTopology()->getConsistent();

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNormalsDialog :: recalculate()
{
    if( mesh == nullptr) return;

    cout << "NORMAL AT " << normalAt << endl;

    if( normalAt == 0)
        mesh->getGeometry()->setNodesNormal();

    if( normalAt == 2)
        mesh->getGeometry()->setFacesNormal();

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshNormalsDialog :: reverseAll()
{
    if( mesh == nullptr) return;

    mesh->getTopology()->reverseAll();
    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNormalsDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->close();

}
///////////////////////////////////////////////////////////////////////////////
void JMeshNormalsDialog :: makeConnections()
{
    LineEdit( lengthLineEdit,  [=] {setLength();});
    CheckBox( displayCheckBox, [=] {checkDisplay();});

    PushButton( colorPushButton, [=] {setColor();});
    PushButton( reverseAllPushButton, [=] {reverseAll();});
    PushButton( recalculatePushButton,  [=] {recalculate(); });
    PushButton( getConsistentPushButton, [=] {getConsistent(); });

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
