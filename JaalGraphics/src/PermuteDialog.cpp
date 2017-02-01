#include "PermuteDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JPermuteDialog :: JPermuteDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

void JPermuteDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

//   mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return;

    size_t numNodes = mesh->getSize(0);

    node1SpinBox->setMinimum(0);
    node1SpinBox->setMaximum(numNodes-1);

    node2SpinBox->setMinimum(0);
    node2SpinBox->setMaximum(numNodes-1);

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JPermuteDialog :: permute()
{
    if( mesh == nullptr ) return;

    int node1 = node1SpinBox->value();
    int node2 = node2SpinBox->value();

    if( node1 == node2 ) return;

    mesh->swapNodePositions(node1,node2);

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JPermuteDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////
void JPermuteDialog :: makeConnections()
{
    PushButton( permutePushButton, [=] {permute();});
    PushButton( closePushButton,   [=] {close();});
}
///////////////////////////////////////////////////////////////////////////////
