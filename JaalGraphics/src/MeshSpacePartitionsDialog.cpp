#include "MeshSpacePartitionsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshSpacePartitionsDialog :: JMeshSpacePartitionsDialog( QWidget *parent) : QDialog(parent)
{
    /*
        setupUi(this);
        makeConnections();

        mesh = nullptr;
        viewManager = nullptr;
       meshViewer  = nullptr;
    i*/

}

///////////////////////////////////////////////////////////////////////////////

JMeshSpacePartitionsDialog :: ~JMeshSpacePartitionsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSpacePartitionsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSpacePartitionsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSpacePartitionsDialog :: genPartitions()
{

}


void JMeshSpacePartitionsDialog :: makeConnections()
{
//    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
