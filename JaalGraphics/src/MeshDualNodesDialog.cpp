#include "MeshDualNodesDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JMeshDualNodesDialog :: JMeshDualNodesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    dualGraph  = nullptr;

    numNodesLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualNodesDialog :: init()
{
    /*
        if( dualViewer ) {
            bool val = dualViewer->getMeshViewer()->isEnabled(0);
            displayNodesCheckBox->setChecked(val);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualNodesDialog :: openAttribDialog()
{
    /*
         if( attribDialog == nullptr ) {
              attribDialog.reset(new JNodeAttributesDialog( this ));
              attribDialog->setMeshViewer(dualViewer->getMeshViewer() );
         }
         this->hide();
         attribDialog->show();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualNodesDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
//       if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualNodesDialog :: checkNodes()
{
    /*
         if( dualViewer == nullptr ) return;
         bool val;

         Mesh *primalMesh = dualViewer->getPrimalMesh();
         if( primalMesh == nullptr ) return;

         val = displayNodesCheckBox->isChecked();
         JMeshViewer *mviewer = dualViewer->getMeshViewer();
         mviewer->displayAll(0, val);

         if(val) {
              dualGraph = primalMesh->getDualGraph(1);
         }

         if( dualGraph) {
              dualViewer->setDualMesh(dualGraph);
              size_t numNodes = dualGraph->getSize(0);
              numNodesLineEdit->setText( QString::number(numNodes) );
              for( size_t i = 0; i < numNodes; i++) {
                   Vertex *node = dualGraph->getNodeAt(i);
                   node->setAttribute("Display", val);
              }
              dualViewer->refreshDisplay();
         }
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualNodesDialog :: deleteNodes()
{
    if( dualViewer == nullptr) return;

    /*
         if( dualGraph ) {
              int err = dualGraph->deleteNodes();
              if( !err) {
                   displayNodesCheckBox->setChecked( false );
                   dualViewer->refreshDisplay();
                   numNodesLineEdit->setText( QString::number(0) );
              }
         }
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualNodesDialog :: makeConnections()
{
    PushButton( displayNodesCheckBox, [=] { checkNodes();});
    PushButton( attributesPushButton, [=] { openAttribDialog(); });
    PushButton( deleteAllnodesPushButton, [=] { deleteNodes(); });
    PushButton( closePushButton, [=] {close(); });
}
///////////////////////////////////////////////////////////////////////////////
