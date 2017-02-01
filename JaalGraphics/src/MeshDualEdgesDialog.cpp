#include "MeshDualEdgesDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshDualEdgesDialog :: JMeshDualEdgesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    dualGraph    = nullptr;
    numEdgesLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualEdgesDialog :: init()
{
    if( dualViewer ) {
//      bool val = dualViewer->getMeshViewer()->isEnabled(0);
//        displayEdgesCheckBox->setChecked(val);
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshDualEdgesDialog :: openAttribDialog()
{
    /*
         if( attribDialog == nullptr ) {
              attribDialog.reset(new JEdgeAttributesDialog( this ));
              attribDialog->setMeshViewer(dualViewer->getMeshViewer() );
         }
         this->hide();
         attribDialog->show();
    */

}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualEdgesDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
//       if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualEdgesDialog :: checkEdges()
{
    /*
         if( dualViewer == nullptr ) return;
         bool val;

         Mesh *primalMesh = dualViewer->getPrimalMesh();
         if( primalMesh == nullptr ) return;

         val = displayEdgesCheckBox->isChecked();
         JMeshViewer *mviewer = dualViewer->getMeshViewer();
         if( mviewer ) mviewer->displayAll(1, val);

         if(val) {
              dualGraph = primalMesh->getDualGraph(1);
         }

         if( dualGraph ) {
              dualViewer->setDualMesh(dualGraph);
              size_t numEdges = dualGraph->getSize(1);
              numEdgesLineEdit->setText( QString::number(numEdges) );
              for( size_t i = 0; i < numEdges; i++) {
                   Edge *edge = dualGraph->getEdgeAt(i);
                   edge->setAttribute("Display", val);
              }
              dualViewer->refreshDisplay();
         }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualEdgesDialog :: deleteEdges()
{
    if( dualViewer == nullptr) return;

    /*
         if( dualGraph ) {
              int err = dualGraph->deleteEdges();
              if( !err) {
                   displayEdgesCheckBox->setChecked( false );
                   dualViewer->refreshDisplay();
                   numEdgesLineEdit->setText( QString::number(0) );
              }
         }
    */

}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualEdgesDialog :: makeConnections()
{
    connect( displayEdgesCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( checkEdges() ));
    connect( attributesPushButton, SIGNAL( clicked( bool ) ), this, SLOT( openAttribDialog() ));
    connect( closePushButton, SIGNAL( clicked( bool ) ), this, SLOT( close() ));
    connect( deleteAllEdgesPushButton, SIGNAL( clicked( bool ) ), this, SLOT( deleteEdges() ));
}
///////////////////////////////////////////////////////////////////////////////
