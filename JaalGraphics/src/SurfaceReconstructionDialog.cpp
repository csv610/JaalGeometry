#include "SurfaceReconstructionDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JSurfaceReconstructionDialog :: JSurfaceReconstructionDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JSurfaceReconstructionDialog :: ~JSurfaceReconstructionDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceReconstructionDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JSurfaceReconstructionDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    size_t numNodes = mesh->getSize(0);
    inNodesLineEdit->setText(QString::number(numNodes));
    outNodesLineEdit->setText(QString::number(0));
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceReconstructionDialog :: genMesh()
{
    JWaitCursor waitCursor;
    waitCursor.start();

    if( newMesh ) meshViewer->removeObject(newMesh);

    JPoissonSurfaceReconstruction srecon;
    srecon.setMesh(mesh);
    srecon.generate();
    newMesh = srecon.getMesh();

    meshViewer->addObject(newMesh);

    size_t numNodes = newMesh->getSize(0);
    outNodesLineEdit->setText(QString::number(numNodes));

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;

    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceReconstructionDialog :: rejectMesh()
{
    if( meshViewer == nullptr) return;

    if( newMesh ) meshViewer->removeObject(newMesh);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 1;
    mrender->displayEntity[2] = 1;
    mrender->displayEntity[3] = 0;
    meshViewer->refreshDisplay();
    outNodesLineEdit->setText(QString::number(0));

}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceReconstructionDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceReconstructionDialog :: makeConnections()
{
    PushButton( applyPushButton,     [=] {genMesh(); });
    PushButton( rejectPushButton,    [=] {rejectMesh(); });
    PushButton( closePushButton,     [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
