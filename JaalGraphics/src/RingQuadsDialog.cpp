#include "RingQuadsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JRingQuadsDialog :: JRingQuadsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JRingQuadsDialog :: ~JRingQuadsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JRingQuadsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JRingQuadsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JRingQuadsDialog :: openMeshSkeletonDialog()
{
    if( skeletonDialog == nullptr)
        skeletonDialog.reset( new JMeshSkeletonDialog(this) );

    skeletonDialog->setViewManager(viewManager);
    skeletonDialog->setParentAfterClose(1);
    skeletonDialog->setMesh(mesh);
    skeletonDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JRingQuadsDialog :: getDisks()
{
}
///////////////////////////////////////////////////////////////////////////////
void JRingQuadsDialog :: getQuadMesh()
{
}
///////////////////////////////////////////////////////////////////////////////
void JRingQuadsDialog :: closeDialog()
{
}

void JRingQuadsDialog :: makeConnections()
{
    PushButton( skeletonPushButton,     [=] {openMeshSkeletonDialog(); });
    PushButton( closePushButton,     [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
