#include "MeshPartitionDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
int JNodePartitionColor :: assign( const JNodePtr &vertex)
{
    if( vertex == nullptr ) return 1;
    if( !vertex->isActive() ) return 1;

    int pid = 0;
    int err = vertex->getAttribute("Partition", pid);
    if( err ) return 1;

    if( colormap.find(pid) == colormap.end() )
        colormap[pid] = JEntityColor::getRandomColor();

    JNodeRenderPtr nAttrib;
    vertex->getAttribute("Render", nAttrib);
    nAttrib->color = colormap[pid];

    if( vertex->hasAttribute("PartitionCorner") ) {
        JColor clr;
        clr[0] = 1.0;
        clr[1] = 1.0;
        clr[2] = 1.0;
        clr[3] = 1.0;
        nAttrib->color = clr;
        nAttrib->scale = 1.5;
        nAttrib->glyph = JNodeRender::NODE_AS_SPHERE;
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JEdgePartitionColor :: assign( const JEdgePtr &edge)
{
    if( edge  == nullptr ) return 1;
    if( !edge->isActive() ) return 1;

    int pid = 0;

    JEdgeRenderPtr eAttrib;
    edge->getAttribute("Render", eAttrib);
    eAttrib->scale  = 1.0; 
    eAttrib->lineWidth = 1; 

    JColor clr;
    clr[0] = 0.0;
    clr[1] = 0.0;
    clr[2] = 0.0;
    clr[3] = 0.0;
    float scale = 1.0;

    int err = edge->getAttribute("Partition", pid);
    if( err ) return 1;
    if( colormap.find(pid) == colormap.end() )
        colormap[pid] = JEntityColor::getRandomColor();
    clr = colormap[pid];

    if( edge->hasAttribute("Interface") ) {
        clr[0] = 1.0;
        clr[1] = 1.8;
        clr[2] = 1.0;
        clr[3] = 1.0;
        scale  = 1.5;
    }
    eAttrib->color = clr;
    eAttrib->scale = scale;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JFacePartitionColor :: assign( const JFacePtr &face)
{
    if( face  == nullptr ) return 1;
    if( !face->isActive() ) return 1;

    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);
    assert(fAttrib);

    int err, pid = 0;
    if( minColors ) {
        err = face->getAttribute("ColorIndex", pid);
        if( err) return 1;
        fAttrib->color = JEntityColor::getMinColor( pid );
    } else {
        err = face->getAttribute("Partition", pid);
        if( err) return 1;
        if( colormap.find(pid) == colormap.end() )
            colormap[pid] = JEntityColor::getRandomColor();
        fAttrib->color = colormap[pid];
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JCellPartitionColor :: assign( const JCellPtr &cell)
{
    int pid = 0;
    int err = cell->getAttribute("Partition", pid);
    if( err ) return 1;

    if( colormap.find(pid) == colormap.end() )
        colormap[pid] = JEntityColor::getRandomColor();

    cell->setAttribute("Color", colormap[pid]);
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

JMeshPartitionDialog :: JMeshPartitionDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;

    numPartsLineEdit->setText(QString::number(2));
}

///////////////////////////////////////////////////////////////////////////////
JMeshPartitionDialog :: ~JMeshPartitionDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str()));

    faceColor.reset(new JFacePartitionColor);
    edgeColor.reset(new JEdgePartitionColor);
    nodeColor.reset(new JNodePartitionColor);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->getViewManager()->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: assignColors()
{
    if( mesh == nullptr ) return;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim == 2 ) {
        faceColor->setMesh(mesh);
        edgeColor->setMesh(mesh);
    }

    /*
        if( topDim == 3 ) {
            JCellPartitionColor cellPartColor;
    //      CellColor::assign(mesh, &cellPartColor);
    //      meshViewer->displayAll(2,0);
    //      meshViewer->displayAll(3,1);
        }
    */
    meshViewer->updateBuffers( mesh );
}
///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: metisPartition()
{
    if( meshViewer == nullptr ) return;

    if( mesh == nullptr ) return;
    QString qstr = numPartsLineEdit->text();
    int nparts  = qstr.toInt();

    if( nparts < 2 ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JMetisPartitioner mp;
    mp.setMesh(mesh);
    mp.getPartitions(nparts);
    assignColors();

}

///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: applyAlgorithm()
{
    if( mesh == nullptr) return;

    size_t numnodes = mesh->getSize(0);
    JNodeRenderPtr nattrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        vertex->getAttribute("Render", nattrib);
        nattrib->glyph = 0;
    }

    JEdgeRenderPtr eattrib;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        edge->getAttribute("Render", eattrib);
        eattrib->glyph = 0;
        eattrib->color[0] = 0.1;
        eattrib->color[1] = 0.1;
        eattrib->color[2] = 0.1;
    }

    QString qs = algorithmComboBox->currentText();
    string str = qs.toUtf8().constData();

    if( str == "Metis") metisPartition();

    JMeshPartitioner mp;
    mp.setMesh(mesh);
    mp.searchInterfaces();
    mp.searchCorners();
    assignColors();
}

///////////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: clearAll()
{
    if( mesh == nullptr ) return;
    mesh->deleteCellAttribute("Partition");
    mesh->deleteFaceAttribute("Partition");
    mesh->deleteFaceAttribute("Interface");
    mesh->deleteEdgeAttribute("Interface");
    mesh->deleteNodeAttribute("PartitionCorner");
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: openCellsPartitionsDialog()
{
    /*
        if( editRegionDialog == nullptr)
            editRegionDialog.reset( new JEditPartitionRegionDialog(this));

        editRegionDialog->setViewManager( viewManager );
        editRegionDialog->setMesh(mesh);
        editRegionDialog->show();
        this->hide();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: openFacesPartitionsDialog()
{
    if( facesPartitionsDialog == nullptr)
        facesPartitionsDialog.reset( new JMeshFacesPartitionsDialog(this));

    facesPartitionsDialog->setViewManager( viewManager );
    facesPartitionsDialog->setMesh(mesh);
    facesPartitionsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: openEdgesPartitionsDialog()
{
    if( edgesPartitionsDialog == nullptr)
        edgesPartitionsDialog.reset( new JMeshEdgesPartitionsDialog(this));


    edgesPartitionsDialog->setViewManager( viewManager );
    edgesPartitionsDialog->setMesh(mesh);
    edgesPartitionsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: openNodesPartitionsDialog()
{
/*
    if( nodesPartitionsDialog == nullptr)
        nodePartitionsDialog.reset( new JMeshNodesPartitionsDialog(this));

    nodesPartitionsDialog->setViewManager( viewManager );
    nodesPartitionsDialog->setMesh(mesh);
    nodesPartitionsDialog->show();
    this->hide();
*/
}
///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: openNormalClustersDialog()
{
    if( normalClustersDialog == nullptr)
        normalClustersDialog.reset( new JMeshNormalClustersDialog(this));

    normalClustersDialog->setViewManager( viewManager );
    normalClustersDialog->setMesh(mesh);
    normalClustersDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: getTopologicalDisks()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JMetisPartitioner mpart;
    mpart.setMesh(mesh);
    mpart.getTopologicalDisks();

    int nParts = mpart.getNumPartitions();

    numPartsLineEdit->setText(QString::number(nParts));

    assignColors();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: closeDialog()
{
    this->close();
    this->parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: savePartitions()
{
    JMetisPartitioner mpart;
    mpart.setMesh(mesh);
    mpart.savePartitions(0);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: getRegionGrowingClusters()
{
    JWaitCursor wCursor;
    wCursor.start();

    JMeshFacesClustering fclusters;
    fclusters.setMesh(mesh);
    fclusters.expand();
 // fclusters.getSpectralClusters(10);

    JMeshPartitioner mp;
    mp.setMesh(mesh);
    mp.searchInterfaces();
    mp.searchCorners();
    assignColors();

    JFaceSequence seeds = fclusters.getSeeds();
    JColor white = JEntityColor::getColor("White");
    JFaceRenderPtr fAttrib;
    for( const JFacePtr &f : seeds) {
         f->getAttribute("Render", fAttrib);
         fAttrib->color = white;
    }

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: removeZigZagInterfaces()
{
    JWaitCursor wCursor;
    wCursor.start();

    JMeshFacesClustering fclusters;
    fclusters.setMesh(mesh);
    fclusters.removeZigZagInterfaces();

    JMeshPartitioner mp;
    mp.setMesh(mesh);
    mp.searchInterfaces();
    mp.searchCorners();
    assignColors();



}
///////////////////////////////////////////////////////////////////////////////
void JMeshPartitionDialog :: makeConnections()
{
    PushButton( applyPushButton, [=] {applyAlgorithm();});
    PushButton( normalClustersPushButton, [=] { openNormalClustersDialog();});
    PushButton( topoDiskPartsPushButton, [=] { getTopologicalDisks();});

    PushButton( cellsPartitionsPushButton, [=] {openCellsPartitionsDialog();});
    PushButton( facesPartitionsPushButton, [=] {openFacesPartitionsDialog();});
    PushButton( edgesPartitionsPushButton, [=] {openEdgesPartitionsDialog();});
    PushButton( nodesPartitionsPushButton, [=] {openNodesPartitionsDialog();});
    PushButton( clearAllPushButton, [=] {clearAll();});
    PushButton( regionGrowingClustersPushButton, [=] {getRegionGrowingClusters();});
    PushButton( removeZigZagInterfacesPushButton, [=] {removeZigZagInterfaces();});

    PushButton( savePartitionsPushButton, [=] {savePartitions();});

    PushButton( closePushButton,  [=] {closeDialog();});

}

///////////////////////////////////////////////////////////////////////////////


void JMeshPartitionDialog :: convexPartition()
{
    if( meshViewer == nullptr ) return;
    QString qstr = objectNameLineEdit->text();

    /*
        if( mesh == nullptr ) return;
        QString str = numPartsLineEdit->text() ;
        int nparts  = str.toInt();

        if( nparts < 2 ) return;

        cout  << "HACD algorithm is pathetic "  << endl;
        JWaitCursor waitCursor;
        waitCursor.start();

        JApproxConvexDecomposition mp;
        mp.setMesh(mesh);
        mp.setNumClusters(nparts);
        mp.getPartitions();
        assignColors();
    */
}


/*
void JMeshPartitionDialog :: searchNodeCorners()
{
    if( mesh == nullptr || meshViewer == nullptr ) return;

    if( nodeAttribDialog == nullptr )
        nodeAttribDialog.reset(new JNodeAttributesDialog( this ));

    nodeAttribDialog->setViewManager( viewManager );
    nodeAttribDialog->setMesh( mesh );

    JNodeSequence nodes;
    size_t numNodes = mesh->getSize(0);
    for( size_t i = 0; i < numNodes; i++)  {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->hasAttribute("PartitionCorner") ) {
            nodes.push_back(vertex);
        }
    }

    nodeAttribDialog->setNodes(nodes);
    nodeAttribDialog->show();
}
*/

////////////////////////////////////////////////////////////////////////////////////

/*
void JMeshPartitionDialog :: searchEdgeInterfaces()
{
    if(mesh == nullptr) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning No mesh was provided");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    JMeshPartitioner mp;
    mp.setMesh(mesh);
    int err = mp.searchInterfaces();
    if( err ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: Elements have not been partitioned: interface search failed ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    assignColors();
    int numInterfaces = mp.getNumInterfaces();
    numInterfacesLineEdit->setText( QString::number(numInterfaces) );
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshPartitionDialog :: searchFaceInterface()
{
    JMeshPartitioner mp;
    mp.setMesh(mesh);
    int err = mp.searchRegions();

    if( err ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: No interface found: region search failed");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    assignColors();

    int numParts = mp.getNumPartitions();
//  numRegionsLineEdit->setText( QString::number(numParts) );
}
*/

///////////////////////////////////////////////////////////////////////////////

void JMeshPartitionDialog :: optPartition()
{
}

///////////////////////////////////////////////////////////////////////////////
/*
void JMeshPartitionDialog :: displayPart()
{
    if( mesh == nullptr || meshViewer == nullptr ) return;

    int partID = displayPartSpinBox->value();

    size_t numnodes = mesh->getSize(0);

    JNodeRenderPtr nAttrib;
    for( size_t i  = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            vertex->getAttribute("Render", nAttrib);
            nAttrib->display = 1;
            nAttrib->glyph   = JNodeDraw::NODE_AS_POINT;
        }
    }

    size_t numfaces = mesh->getSize(2);
    int pid = 0;
    JFaceRenderPtr fAttrib;
    JFacePtr lastface;
    if( displayPartCheckBox->isChecked() ) {
        for( size_t i  = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Partition", pid );
                face->getAttribute("Render", fAttrib);
                fAttrib->display = 0;
                if( pid == partID ) {
                    fAttrib->display = 1;
                    lastface = face;
                }
            }
        }
//        if( lastface ) meshViewer->lookAt(lastface->getNodeAt(0));
    } else {
        for( size_t i  = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            int err = face->getAttribute("Render", fAttrib);
            if( !err ) fAttrib->display = 1;
        }
    }

    meshViewer->updateBuffers(mesh);
}
*/

