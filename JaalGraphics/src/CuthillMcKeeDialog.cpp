#include "CuthillMcKeeDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JCuthillMcKeeDialog :: JCuthillMcKeeDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JCuthillMcKeeDialog :: ~JCuthillMcKeeDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JCuthillMcKeeDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JCuthillMcKeeDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JCuthillMcKeeDialog :: applyAlgorithm()
{
}
///////////////////////////////////////////////////////////////////////////////

void JCuthillMcKeeDialog :: closeDialog()
{
}
///////////////////////////////////////////////////////////////////////////////

void JCuthillMcKeeDialog :: makeConnections()
{
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
