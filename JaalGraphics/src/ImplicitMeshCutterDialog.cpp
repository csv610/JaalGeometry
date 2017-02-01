#include "ImplicitMeshCutterDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JImplicitMeshCutterDialog :: JImplicitMeshCutterDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JImplicitMeshCutterDialog :: ~JImplicitMeshCutterDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JImplicitMeshCutterDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}
///////////////////////////////////////////////////////////////////////////////
void JImplicitMeshCutterDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    double len =  box.getDiameter();
}
///////////////////////////////////////////////////////////////////////////////
void JImplicitMeshCutterDialog :: setPlane()
{
    double t;
    t = nxScrollBar->value()/100.0;
    double nx = -1.0 + 2.0*t;

    t = nyScrollBar->value()/100.0;
    double ny = -1.0 + 2.0*t;

    t = nzScrollBar->value()/100.0;
    double nz = -1.0 + 2.0*t;

    cout << nx << "  " << nz << "  " << endl;
}
///////////////////////////////////////////////////////////////////////////////
void JImplicitMeshCutterDialog :: makeConnections()
{
    /*
        connect( nxScrollBar, SIGNAL( valueChanged(int) ), this, SLOT( setPlane() ));
        connect( nyScrollBar, SIGNAL( valueChanged(int) ), this, SLOT( setPlane() ));
        connect( nzScrollBar, SIGNAL( valueChanged(int) ), this, SLOT( setPlane() ));
    */

    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
