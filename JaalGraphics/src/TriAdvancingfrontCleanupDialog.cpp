#include "TriAdvancingfrontCleanupDialog.hpp"
///////////////////////////////////////////////////////////////////////////////

JTriAdvancingfrontCleanupDialog :: JTriAdvancingfrontCleanupDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();
    mesh = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JTriAdvancingfrontCleanupDialog :: ~JTriAdvancingfrontCleanupDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JTriAdvancingfrontCleanupDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;
//   mesh = meshViewer->getMesh();
}

///////////////////////////////////////////////////////////////////////////////

void JTriAdvancingfrontCleanupDialog :: makeConnections()
{
    PushButton( closePushButton,  [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
