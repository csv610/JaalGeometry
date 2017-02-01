#include "MeshExtrudeDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshExtrudeDialog :: JMeshExtrudeDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
    distanceLineEdit->setText(QString::number(1.0));
}

///////////////////////////////////////////////////////////////////////////////

JMeshExtrudeDialog :: ~JMeshExtrudeDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshExtrudeDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshExtrudeDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JMeshExtrudeDialog :: genMesh()
{
    if( meshExtrude == nullptr)
        meshExtrude.reset( new JMeshExtrude);
    meshExtrude->setMesh(mesh);


    QString qstr = distanceLineEdit->text();
    double d = qstr.toDouble();

    JMeshPtr newmesh = meshExtrude->getMesh(d);
    meshViewer->addObject( newmesh );

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;
    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////
void JMeshExtrudeDialog :: closeDialog()
{
    close();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshExtrudeDialog :: makeConnections()
{
    PushButton( extrudePushButton,   [=] {genMesh(); });
    PushButton( closePushButton,     [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
