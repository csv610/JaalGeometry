#include "MeshRefine3DDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshRefine3DDialog :: JMeshRefine3DDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    meshViewer  = nullptr;
    viewManager = nullptr;

    int nx = 2;
    nxLineEdit->setText( QString::number(nx ) );
    nyLineEdit->setText( QString::number(nx ) );
    nzLineEdit->setText( QString::number(nx ) );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRefine3DDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine3DDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 3;

    string name = mesh->getName();
//  objectNameLineEdit->setText( QString(name.c_str() ) );
//  entityPicker->setMesh(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshRefine3DDialog :: refineHex17()
{
    if( meshViewer == nullptr ) return;

//     hexmesh = meshViewer->getMesh();
    if( hexmesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JHexRefiner::refine17(hexmesh);

//   meshViewer->setNewMesh(hexmesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshRefine3DDialog :: refineHex18()
{
    if( meshViewer == nullptr ) return;

//   hexmesh = meshViewer->getMesh();
    if( hexmesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JHexRefiner::refine18(hexmesh);

//     meshViewer->setNewMesh(hexmesh);
}


///////////////////////////////////////////////////////////////////////////////
void JMeshRefine3DDialog :: insertPillows()
{
    /*
         if( meshViewer == nullptr ) return;
         quadmesh  = meshViewer->getMesh();
         if( quadmesh == nullptr ) return;

         QuadCleanUp qClean(quadmesh);
         qClean.insert_boundary_pillows();

         meshViewer->setNewMesh(quadmesh);
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshRefine3DDialog :: genHexBlocks()
{
    /*
         QString str;
         str = nxLineEdit->text();
         int nx = str.toInt();

         str = nyLineEdit->text();
         int ny = str.toInt();

         str = nzLineEdit->text();
         int nz = str.toInt();
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshRefine3DDialog :: makeConnections()
{
    PushButton( hex17PushButton,  [=] {refineHex17();});
    PushButton( hex18PushButton,  [=] {refineHex18();});
    PushButton( blocksPushButton, [=] {genHexBlocks();});
    PushButton( closePushButton,  [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
