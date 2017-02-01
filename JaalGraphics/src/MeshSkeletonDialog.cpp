#include "MeshSkeletonDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshSkeletonDialog :: JMeshSkeletonDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
    speedLineEdit->setText( QString::number(0.1));
    centralityLineEdit->setText( QString::number(0.2));
}

///////////////////////////////////////////////////////////////////////////////

JMeshSkeletonDialog :: ~JMeshSkeletonDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c1 = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c1);

    JViewComponentPtr c2 = viewManager->getComponent("ShapeViewer");
    if( c2 == nullptr)
        c2 = JShapeViewer::registerComponent(viewManager);

    shapeViewer = dynamic_pointer_cast<JShapeViewer>(c2);

/*
    if( meshViewer )
        setMesh( meshViewer->getCurrentMesh() );
*/
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: setMesh( const JMeshPtr &amesh)
{
/*
    surfMesh.reset();
    if( amesh == nullptr ) return ;

    QMessageBox msg;
    if( !amesh->getTopology()->isClosed() ) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The input mesh is not closed: The skeleton will not be extraded");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    JEdgeSequence edges = amesh->getTopology()->getNonManifoldEdges();

    if( !amesh->getTopology()->isManifold() ) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The input mesh non-manifold : The skeleton will not be extraded");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    inMesh = amesh;
    string name = inMesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    meshSkel.reset( new JMeshSkeleton);
    surfMesh.reset();

    meanCurvature = 1;
    lapSkel.reset( new JSimpleLaplacianSkeleton);
*/
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: initMesh()
{
/*
    if(surfMesh) return;

    if(meanCurvature ) {
        meshSkel->setMesh(inMesh);
        int nrefine = refineMeshSpinBox->value();
        meshSkel->setRefinement(nrefine);
        surfMesh  = meshSkel->getWorkingMesh();
    } else {
        lapSkel->setMesh(inMesh);
        surfMesh  = lapSkel->getWorkingMesh();
        meshViewer->addObject(surfMesh);
    }

    if( surfMesh == nullptr) return;

    meshViewer->addObject(surfMesh);
    string name = surfMesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    // Keep the original mesh on the top...
    meshViewer->setCurrentMesh(surfMesh);

    JMeshRenderPtr mrender;
    inMesh->getAttribute("Render", mrender);
    mrender->display = 1;
    mrender->setTransparency(1);
    mrender->setAlpha(0.5);
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 1;
    meshViewer->updateBuffers(inMesh);
*/
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: displaySkeleton()
{
/*
    if( skelGraph == nullptr) return;

    int  numBranches = meshSkel->getNumBranches();
    if( numBranches < 1) return;

    int  numNodes  = skelGraph->getSize(0);
    JNodeRenderPtr vattrib;
    double radius;
    if( displaySpheresCheckBox->isChecked() ) {
        for( int i = 0; i < numNodes; i++) {
            JNodePtr vtx = skelGraph->getNodeAt(i);
            vtx->getAttribute("Radius", radius);
            vtx->getAttribute("Render", vattrib);
            vattrib->glyph  = JNodeRender::NODE_AS_SPHERE;
            vattrib->ballRadius = radius;
        }
    } else {
        for( int i = 0; i < numNodes; i++) {
            JNodePtr vtx = skelGraph->getNodeAt(i);
            vtx->getAttribute("Render", vattrib);
            vattrib->glyph  = JNodeRender::NODE_AS_POINT;
        }
    }

    JEdgeRenderPtr attrib;
    for( int  i = 0; i < numBranches; i++) {
        JColor color = JEntityColor::getRandomColor();
        JMeshSkeletonBranchPtr branch = meshSkel->getBranch(i);
        for( const JEdgePtr &e : branch->skelEdges) {
            e->getAttribute("Render", attrib);
            attrib->lineWidth = 2;
            attrib->color = color;
        }
    }
    meshViewer->updateBuffers(skelGraph);

    // Now the working mesh is collapsed, do not show it ..
    JMeshRenderPtr mrender;
    surfMesh->getAttribute("Render", mrender);
    mrender->display = 0;
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 0;
    meshViewer->updateBuffers(surfMesh);
*/

}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: getSkeleton()
{
/*
    if( inMesh == nullptr) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("No valid input mesh present");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    JWaitCursor wcursor;
    wcursor.start();

    initMesh();
    if( surfMesh == nullptr) return;

    QString qstr;
    if( skelGraph == nullptr) {
        qstr = speedLineEdit->text();
        double speed = qstr.toDouble();
        meshSkel->setSpeedQuality(speed);

        qstr = centralityLineEdit->text();
        double center= qstr.toDouble();
        meshSkel->setCenterQuality(center);

        skelGraph = meshSkel->getSkeleton();
        meshViewer->addObject(skelGraph);
    }

    displaySkeleton();
*/
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: openNodesDialog()
{
    if( skelGraph == nullptr ) return;

    if( nodeAttribDialog == nullptr )
        nodeAttribDialog.reset(new JNodeAttributesDialog( this ));

    nodeAttribDialog->setViewManager( viewManager );
    nodeAttribDialog->setMesh(skelGraph);

    JNodeSequence nodes = skelGraph->getNodes();

    nodeAttribDialog->setNodes(nodes);
    nodeAttribDialog->show();
    this->hide();

}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: openLeafNodesDialog()
{
/*
    if( skelGraph == nullptr ) return;

    if( nodeAttribDialog == nullptr )
        nodeAttribDialog.reset(new JNodeAttributesDialog( this ));

    nodeAttribDialog->setViewManager( viewManager );
    nodeAttribDialog->setMesh(skelGraph);

    JNodeSequence nodes = meshSkel->getLeafNodes();

    nodeAttribDialog->setNodes(nodes);
    nodeAttribDialog->show();
    this->hide();
*/

}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: openJunctionNodesDialog()
{
/*
    if( skelGraph == nullptr ) return;

    if( nodeAttribDialog == nullptr )
        nodeAttribDialog.reset(new JNodeAttributesDialog( this ));

    nodeAttribDialog->setViewManager( viewManager );
    nodeAttribDialog->setMesh(skelGraph);

    JNodeSequence nodes = meshSkel->getJunctionNodes();

    nodeAttribDialog->setNodes(nodes);
    nodeAttribDialog->show();
    this->hide();
*/
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: openEdgesDialog()
{
    QMessageBox msg;
    if( skelGraph == nullptr ) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Skeleton has not been extracted for the mode");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    if( edgeAttribDialog.get() == nullptr )
        edgeAttribDialog.reset(new JEdgeAttributesDialog( this ));

    edgeAttribDialog->setViewManager( viewManager );
    edgeAttribDialog->setMesh(skelGraph);

    JEdgeSequence edges = skelGraph->getEdges();

    edgeAttribDialog->setEdges(edges);
    edgeAttribDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: createSurfaceParts()
{
    JWaitCursor wcursor;
    wcursor.start();

    if( skelGraph == nullptr) getSkeleton();

    displaySurfaceParts();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: displaySurfaceParts()
{
    JMeshRenderPtr mrender;
    inMesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 1;
    mrender->displayEntity[3] = 0;
    mrender->setTransparency(0);
    mrender->setAlpha(1.0);

    JFacePartitionColor faceColor;
    faceColor.setMesh(inMesh);

    JMeshPartitioner mpart;
    mpart.setMesh(inMesh);

    int nparts  = mpart.getNumPartitions();
    JFaceSequence ambiguous;
    mpart.getPartition( nparts, ambiguous);
    JFaceRenderPtr fAttrib;
    JColor black = JEntityColor::getColor("Black");
    for( const JFacePtr &f: ambiguous) {
         f->getAttribute("Render", fAttrib);
         fAttrib->color = black;
    }


/*
    JNodeRenderPtr vAttrib;
    size_t numNodes = surfMesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &node = surfMesh->getNodeAt(i);
        node->getAttribute("Render", vAttrib);
        vAttrib->display = 1;
    }

    int partID = 0;
    size_t numFaces = surfMesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = surfMesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 0;
        face->getAttribute("Partition", partID);
        if( partID ) {
            fAttrib->display = 1;
            for( int j = 0; j < face->getSize(0); j++) {
                const JNodePtr &node = face->getNodeAt(j);
                node->getAttribute("Render", vAttrib);
                vAttrib->display = 0;
            }
        }
    }
*/
    meshViewer->updateBuffers(inMesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: closeDialog()
{
    if( skelGraph ) {
        meshViewer->removeObject(skelGraph);
        skelGraph.reset();
    }

    if( inMesh ) {
        JMeshRenderPtr mrender;
        inMesh->getAttribute("Render", mrender);
        mrender->display = 1;
        mrender->setTransparency(0);
        mrender->setAlpha(1.0);
        mrender->displayEntity[0] = 1;
        mrender->displayEntity[1] = 1;
        mrender->displayEntity[2] = 1;
        meshViewer->updateBuffers(inMesh);
    }

    if( surfMesh ) {
        meshViewer->removeObject(surfMesh);
        surfMesh.reset();
    }

    if( showParentAfterClose)
        this->parentWidget()->show();
    this->close();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: openSkelEditingDialog()
{
/*
    if( skelEditingDialog == nullptr )
        skelEditingDialog.reset(new JMeshSkeletonEditingDialog( this ));

    skelEditingDialog->setViewManager( viewManager );
    skelEditingDialog->setSkeleton(meshSkel);

    skelEditingDialog->show();
    this->hide();
*/
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: openSkelShapesDialog()
{

/*
    if( skelShapesDialog == nullptr )
        skelShapesDialog.reset(new JMeshSkeletonShapesDialog( this ));

    skelShapesDialog->setViewManager( viewManager );
    skelShapesDialog->setMesh(surfMesh);
    skelShapesDialog->setSkeleton(meshSkel);

    skelShapesDialog->show();
    this->hide();
*/
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: contractGeometry()
{
/*
    if( meshSkel == nullptr) return;

    JWaitCursor wcursor;
    wcursor.start();

    meanCurvature = meanCurvatureRadioButton->isChecked();

    initMesh();

    QString qstr;
    if( meanCurvature ) {
        qstr = speedLineEdit->text();
        double speed = qstr.toDouble();
        meshSkel->setSpeedQuality(speed);

        qstr = centralityLineEdit->text();
        double center= qstr.toDouble();
        meshSkel->setCenterQuality(center);

        meshSkel->contractGeometry();

        // Remove the Old mesh ..
        meshViewer->removeObject(surfMesh);

        // get the new mesh ...
        surfMesh  = meshSkel->getWorkingMesh();
        meshViewer->addObject(surfMesh);
    } else {
        for( int i = 0; i < 10; i++)
            lapSkel->applyOneStep();
    }

    // Keep the original mesh on the top...
    meshViewer->setCurrentMesh(surfMesh);

    // See the Skeleton through the wireframe by default ...
    JEdgeRenderPtr eAttrib;
    size_t numEdges = surfMesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = surfMesh->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        eAttrib->display = 0;
        eAttrib->lineWidth = 1;
    }

    JColor red = JEntityColor::getColor("Red");
    JFaceRenderPtr fAttrib;
    size_t numFaces = surfMesh->getSize(2);
    for( size_t i = 0; i < numFaces; i++) {
        const JFacePtr &face = surfMesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
        fAttrib->color   = red;
    }
    meshViewer->updateBuffers(surfMesh);
*/
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: openMeshRenderDialog()
{
    if( meshRenderDialog == nullptr )
        meshRenderDialog.reset(new JMeshRenderDialog( this ));

    meshRenderDialog->setViewManager( viewManager );
    meshRenderDialog->setMesh(inMesh);

    meshRenderDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSkeletonDialog :: loadSkeleton()
{
    /*
        QString qstr = QFileDialog::getOpenFileName(this,
                       *new QString("Select Mesh File "),
                       lastSelectedDirectory,
                       *new QString( "Mesh Format (*.off"));

        string filename = qstr.toUtf8().constData();
        JMeshIO mio;
        skelGraph = mio.readFile(filename);
        meshviewer->addObject(skelGraph);
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: saveSkeleton()
{
    /*
        QString qstr  = QFileDialog::getSaveFileName(this,
                        *new QString("Select Mesh File "),
                        QString::fromStdString( name ),
                        *new QString( "Mesh Format (*.off)"));
        string filename = qstr.toUtf8().constData();
        JMeshIO mio;
        mio.saveAs(skelGraph, filename);
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: rejectSkeleton()
{
    if( skelGraph ) {
        meshViewer->removeObject(skelGraph);
        skelGraph.reset();
    }
    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: displayDisks()
{
    if( skelGraph == nullptr) return;

    JWaitCursor wcursor;
    wcursor.start();

/*
    JMeshSkeletonContours skCont;
    skCont.setSkeleton( meshSkel );
    contourDisks = skCont.getDisks();
    meshViewer->addObject( contourDisks);
*/

    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSkeletonDialog :: displaySpheres()
{
/*
    if( skelGraph == nullptr) return;
    JNodeRenderPtr vAttrib;

    size_t numnodes = skelGraph->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr  &vtx = skelGraph->getNodeAt(i);
        vtx->getAttribute("Render", vAttrib);
        vAttrib->glyph   = 0;
        vAttrib->display = 0;
    }
    if( displaySpheresCheckBox->isChecked() ) {

        double radius;
        JNodeSequence junctionNodes = meshSkel->getJunctionNodes();
        for( const JNodePtr &vtx : junctionNodes) {
            vtx->getAttribute("Render", vAttrib);
            vAttrib->glyph = 1;
            vtx->getAttribute("Radius", radius);
            vAttrib->ballRadius = radius;
            vAttrib->display = 1;
        }

        JNodeSequence leafNodes = meshSkel->getLeafNodes();

        for( const JNodePtr &vtx : leafNodes) {
            vtx->getAttribute("Render", vAttrib);
            vAttrib->glyph = 1;
            vtx->getAttribute("Radius", radius);
            vAttrib->ballRadius = radius;
            vAttrib->display = 1;
        }
    }
    meshViewer->updateBuffers(skelGraph);
*/
}
///////////////////////////////////////////////////////////////////////////////


void JMeshSkeletonDialog :: makeConnections()
{
    CheckBox( displayDisksCheckBox,  [=] {displayDisks(); });
    CheckBox( displaySpheresCheckBox,[=] {displaySpheres(); });
    CheckBox( displaySurfaceCheckBox,[=] {displaySurfaceParts(); });

    PushButton( meshRenderPushButton,  [=] {openMeshRenderDialog(); });

    PushButton( contractGeomPushButton,  [=] {contractGeometry(); });
    PushButton( nodesAttribPushButton,   [=] {openNodesDialog(); });
    PushButton( leafNodesPushButton,     [=] {openLeafNodesDialog(); });
    PushButton( junctionNodesPushButton, [=] {openJunctionNodesDialog(); });

    PushButton( edgeAttribPushButton,    [=] {openEdgesDialog(); });
    PushButton( getSkeletonPushButton,   [=] {getSkeleton(); });
    PushButton( skelEditingPushButton,   [=] {openSkelEditingDialog(); });
    PushButton( shapesFittingPushButton,   [=] {openSkelShapesDialog(); });

    PushButton( rejectPushButton,    [=] {rejectSkeleton(); });
    PushButton( closePushButton,     [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
