#include "GPOperatorsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JGPOperatorsDialog :: JGPOperatorsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    mesh = nullptr;
    currEdge = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

}

///////////////////////////////////////////////////////////////////////////////

JGPOperatorsDialog :: ~JGPOperatorsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
/*
void JGPOperatorsDialog :: setColor()
{
     float rgb[3];
     QColor color = QColorDialog::getColor();
     rgb[0] = color.red()/255.0;
     rgb[1] = color.green()/255.0;
     rgb[2] = color.blue()/255.0;
     if( singletColor) singletColor->setColor(rgb);

     if( meshViewer == nullptr ) return;
     meshViewer->getViewManager()->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////

void JGPOperatorsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

//  mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;
    /*
         qdual = new QuadDual( mesh );

         if( !isQuadMesh() ) {
              QMessageBox msg;
              msg.setIcon(QMessageBox::Warning);
              msg.setText("At present singlet operations are only for all quad mesh ");
              msg.setStandardButtons( QMessageBox::Ok);
              int ret = msg.exec();
              if( ret == QMessageBox::Ok ) {
                   return;
              }
         }

         mesh->buildRelations(0,2);
         meshViewer->setNodeColorMethod( singletColor );

         int nSize = GPOperators::getSize(mesh);
         numGPOperatorssLineEdit->setText( QString::number(nSize) );
    */
}

///////////////////////////////////////////////////////////////////////////////
void JGPOperatorsDialog :: getSeeds()
{
    if( meshViewer == nullptr ) return;
    JMeshEntityPickerPtr entityPicker = meshViewer->getEntityPicker();

    JEdgeSequence edges = entityPicker->getPickedEdges();
    if( edges.size() != 1 ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Only one edge must be selected ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    JFaceSequence faces = entityPicker->getPickedFaces();
    if( faces.size() != 1 ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Only one face must be selected ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    JFacePtr gpFace = faces[0];
    if( gpFace->getType() != JFace::QUADRILATERAL ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The selected face must be quadrilateral");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    /*
      Edge *gpEdge = edges[0];
      meshViewer->setGPOD( gpEdge, gpFace );
    */
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JGPOperatorsDialog :: shiftLeftOp()
{
    if( mesh == nullptr ) return ;
    if( currEdge == nullptr ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("No initial edge picked: Pick it now ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
}

///////////////////////////////////////////////////////////////////////////////
void JGPOperatorsDialog :: shiftRightOp()
{
    if( mesh == nullptr ) return ;
    if( currEdge == nullptr ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("No initial edge picked: Pick it now ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
}
///////////////////////////////////////////////////////////////////////////////

void JGPOperatorsDialog :: collapseOp()
{
    if( mesh == nullptr ) return ;
    if( currEdge == nullptr ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("No initial edge picked: Pick it now ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
}

void JGPOperatorsDialog :: makeConnections()
{
    connect( shiftLeftPushButton,  SIGNAL( clicked() ), this, SLOT( shiftLeftOp() ));
    connect( shiftRightPushButton, SIGNAL( clicked() ), this, SLOT( shiftRightOp() ));
    connect( collapsePushButton,  SIGNAL( clicked() ),  this, SLOT( collapseOp() ));
    connect( seedsPushButton,  SIGNAL( clicked() ),  this, SLOT( getSeeds() ));
    connect( closePushButton,  SIGNAL( clicked() ),  this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
