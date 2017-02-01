#include "SketchQuadsDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
JSketchQuadsDialog :: JSketchQuadsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JSketchQuadsDialog :: ~JSketchQuadsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JSketchQuadsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
};
///////////////////////////////////////////////////////////////////////////////

void JSketchQuadsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;
    rayTracer.reset( new JRayTracer);
    rayTracer->setMesh(m);
}

///////////////////////////////////////////////////////////////////////////////

void JSketchQuadsDialog :: makeConnections()
{
//  connect(closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
