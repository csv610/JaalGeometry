#include "MeshSkeletonShapesDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshSkeletonShapesDialog :: JMeshSkeletonShapesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshSkeletonShapesDialog :: ~JMeshSkeletonShapesDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonShapesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonShapesDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonShapesDialog :: getContours()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonShapesDialog :: closeDialog()
{
    parentWidget()->show();
    this->close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonShapesDialog :: makeConnections()
{
    PushButton( closePushButton,     [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
