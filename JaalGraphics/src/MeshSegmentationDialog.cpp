#include "MeshSegmentationDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshSegmentationDialog :: JMeshSegmentationDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    meshViewer  = nullptr;
    viewManager = nullptr;

    creaseAngleLineEdit->setText( QString::number(30) );
    normalAngleLineEdit->setText( QString::number(5) );
    numEdgesLineEdit->setText( QString::number(0) );
    numPlanesLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshSegmentationDialog :: ~JMeshSegmentationDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSegmentationDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
}

void JMeshSegmentationDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSegmentationDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSegmentationDialog :: openMeshGeodesicsDialog()
{
    if( meshGeodesicsDialog.get()  == nullptr )
        meshGeodesicsDialog.reset(new JMeshGeodesicsDialog( this ));

    meshGeodesicsDialog->setViewManager(viewManager);
    meshGeodesicsDialog->show();
    this->hide();
}


///////////////////////////////////////////////////////////////////////////////
void JMeshSegmentationDialog :: segment()
{
    if( mesh == nullptr ) return;

    QString str = creaseAngleLineEdit->text() ;
    double angle  = str.toDouble();

    int ncount = mesh->getGeometry()->setSharpEdges( angle );
    numEdgesLineEdit->setText( QString::number(ncount) );

    // Register Sharp edges as features...
// JEdge::featureSet->addAttribute( "SharpEdge");
//  meshViewer->getEdgeDraw()->setColorMethod( edgeColor );

    /*
         MeshSegmentation mseg(mesh);
         mseg.addStopper("SharpEdge");
         mseg.getPartitions();
    //   partColor->assignColors( mesh );
         meshViewer->getDrawFace()->setColorMethod(partColor);
    */
    meshViewer->getViewManager()->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSegmentationDialog :: getSDFSegments()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JMeshSDF  sdf;
    sdf.setMesh(mesh);
    sdf.setNumOfClusters(16);
    sdf.execute();

    JFaceColorPtr faceColor(new JFacePartitionColor);
    faceColor->setMesh(mesh);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 1;
    mrender->displayEntity[3] = 0;

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSegmentationDialog :: openMeshPartitionDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSegmentationDialog :: closeDialog()
{
    this->close();
     parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSegmentationDialog :: makeConnections()
{
    PushButton( applyPushButton,  [=] {segment();});
    PushButton( sdfPushButton,    [=] {getSDFSegments();});
    PushButton( geodesicsPushButton, [=] {openMeshGeodesicsDialog();});
    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
