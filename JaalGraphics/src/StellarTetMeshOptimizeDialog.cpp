#include "StellarTetMeshOptimizeDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JStellarTetMeshOptimizeDialog :: JStellarTetMeshOptimizeDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JStellarTetMeshOptimizeDialog :: ~JStellarTetMeshOptimizeDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JStellarTetMeshOptimizeDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

//  mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;
}

///////////////////////////////////////////////////////////////////////////////
void JStellarTetMeshOptimizeDialog :: genConfigFile()
{

}
///////////////////////////////////////////////////////////////////////////////
void JStellarTetMeshOptimizeDialog :: optmesh()
{

}

///////////////////////////////////////////////////////////////////////////////
void JStellarTetMeshOptimizeDialog :: makeConnections()
{
    PushButton( optPushButton, [=] {optmesh(); });
}

///////////////////////////////////////////////////////////////////////////////
