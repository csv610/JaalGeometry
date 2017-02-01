#include "MeshEdgesDialog.hpp"

JMeshEdgesDialog :: JMeshEdgesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    drawEdge = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    renderStyle = WIREFRAME;

    numEdgesLineEdit->setText( QString::number(0) );
    numBoundaryLineEdit->setText( QString::fromStdString("Unknown"));
    numVisibleLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    size_t numEdges = mesh->getActiveSize(1);
    numEdgesLineEdit->setText( QString::number(numEdges) );

    lookAtSpinBox->setMinimum( 0);
    lookAtSpinBox->setMaximum( numEdges);

    countBoundEdges();

    if( drawEdge ) drawEdge->setMesh(mesh);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    bool val = mrender->displayEntity[1];
    displayCheckBox->setChecked(val);
   
    setNumVisible();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    drawEdge = meshViewer->getEdgeDraw();
//  checkDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: showEvent( QShowEvent *e)
{
    if( meshViewer == nullptr || mesh == nullptr) return;

    size_t numedges = mesh->getSize(1);
    size_t numVis   = 0;
    JEdgeRenderPtr attrib;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", attrib);
            if( attrib->display) numVis++;
        }
    }
    numVisibleLineEdit->setText( QString::number(numVis) );

    QDialog::showEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: openAttribListDialog()
{
    if( attribListDialog.get() == nullptr )
        attribListDialog.reset(new JMeshEntityAttribListDialog( this ));

    attribListDialog->setMesh( mesh, 1);
    attribListDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: setNumVisible()
{
    if( mesh == nullptr ) return;
    size_t nCount = meshViewer->getNumVisible(1);
    numVisibleLineEdit->setText( QString::number(nCount) );
    bool val  = nCount == 0 ? 0:1;
    displayCheckBox->setChecked(val);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: keyPressEvent( QKeyEvent *e)
{
    // This is ugly requirement as the "Return key" will close the dialog box.
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->getViewManager()->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesDialog :: checkDisplay()
{
    if( mesh == nullptr) return;

    bool val;
    val  = displayCheckBox->isChecked();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[1] = val;

    size_t numedges = mesh->getSize(1);
    size_t numVis   = 0;
    JEdgeRenderPtr attrib;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", attrib);
            attrib->display = val;
            if( val ) numVis++;
        }
    }
    meshViewer->updateTopologyBuffers(mesh);
    numVisibleLineEdit->setText( QString::number(numVis) );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: countBoundEdges()
{
    if( mesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();


    mesh->getTopology()->searchBoundary();
    size_t numedges = mesh->getTopology()->getBoundarySize(1);
    numBoundaryLineEdit->setText( QString::number(numedges) );
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: lookAt()
{
    if( viewManager == nullptr || mesh == nullptr ) return;

    Vec3F srcVec, dstVec, perpAxis;
    GLdouble mat[16];

    size_t id = lookAtSpinBox->value();
    const JEdgePtr &edge = mesh->getEdgeAt(id);
    const JNodePtr &vtx  = edge->getNodeAt(0);

    Vec3F normal;
    vtx->getAttribute("Normal", normal);

    Vec vec;
    vec = viewManager->camera()->position();
    double dist = sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
    vec[0] = dist*normal[0];
    vec[1] = dist*normal[1];
    vec[2] = dist*normal[2];

    viewManager->camera()->setPosition(vec);

    const Point3D &p3d = vtx->getXYZCoords();
    vec[0] = p3d[0];
    vec[1] = p3d[1];
    vec[2] = p3d[2];

    viewManager->camera()->lookAt( vec );

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesDialog :: setInternal()
{
    if( mesh == nullptr) return;

    if( attribDialog.get() == nullptr )
        attribDialog.reset(new JEdgeAttributesDialog( this ));

    attribDialog->setViewManager(viewManager);
    attribDialog->setMesh(mesh);

    JEdgeSequence eseq;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = mesh->getEdgeAt(i);
        if( e->isActive() && !e->isBoundary()  ) eseq.push_back(e);
    }
    attribDialog->setEdges(eseq);

    attribDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: setBoundary()
{
    if( mesh == nullptr) return;

    if( attribDialog.get() == nullptr )
        attribDialog.reset(new JEdgeAttributesDialog( this ));

    attribDialog->setViewManager(viewManager);
    attribDialog->setMesh( mesh );

    JEdgeSequence eseq;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = mesh->getEdgeAt(i);
        if( e->isActive() &&  e->isBoundary()  ) eseq.push_back(e);
    }

    attribDialog->setEdges(eseq);
    attribDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: checkState()
{
    if( meshViewer == nullptr) return;

    bool val;
    JMeshRenderPtr attrib;
    mesh->getAttribute("Render", attrib);
    attrib->displayEntity[1] = displayCheckBox->isChecked();

    if( drawEdge ) {
        val = antiAliasCheckBox->isChecked();
        drawEdge->setAntiAliasing(val);
    }

    attrib->displayIDs[1] =enumCheckBox->isChecked();

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: setDefaultColorMethod()
{
/*
    if( meshViewer == nullptr ) return;
//  meshViewer->displayAll(1,1);
    meshViewer->refreshDisplay();
*/
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesDialog :: openEditDialog()
{
    if( editDialog.get() == nullptr )
        editDialog.reset(new JMeshEdgeEditDialog( this ));

    editDialog->setViewManager( viewManager );
    editDialog->setMesh( mesh);
    editDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesDialog :: saveBoundary()
{
    static QString lastSelectedDirectory;

    QString qstr  = QFileDialog::getSaveFileName(this,
                    *new QString("Select File Name "),
                    lastSelectedDirectory,
                    *new QString( "Mesh Format (*.off)"));

    string fileName = StdString(qstr);
    if (fileName.empty()) return;

    JEdgeSequence bedges;
    mesh->getTopology()->getBoundary(bedges);

    JNodeSequence bnodes;
    mesh->getTopology()->getBoundary(bnodes);

    ofstream ofile(fileName.c_str(), ios::out);
    if( ofile.fail() ) return;
    ofile << "OFF" << endl;
    ofile << bnodes.size() << "  0  " << bedges.size() << endl;
    int index = 0;
    map<JNodePtr, int> localID;
    for( const JNodePtr &v : bnodes) {
        const Point3D &xyz = v->getXYZCoords();
        ofile << xyz[0] << "  " << xyz[1] << " " << xyz[2] << endl;
        localID[v] = index++;
    }

    for( const JEdgePtr &e : bedges) {
        const JNodePtr &v0 = e->getNodeAt(0);
        const JNodePtr &v1 = e->getNodeAt(1);
        ofile << "2 " << localID[v0] << " " << localID[v1] << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgesDialog :: getNonManifoldEdges()
{
    if( mesh == nullptr) return;

    JEdgeSequence edges= mesh->getTopology()->getNonManifoldEdges();
    if( edges.empty() ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Information);
        msg.setText("There are no non-manifold edges in the mesh");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    size_t numedges = mesh->getSize(1);

    JEdgeRenderPtr eattrib;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        edge->getAttribute("Render", eattrib);
        eattrib->display = 0;
    }

    size_t numfaces = mesh->getSize(2);
    JFaceRenderPtr fattrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAttribute("Render", fattrib);
        fattrib->display = 0;
    }

    JFaceSequence faceneighs;
    for( const JEdgePtr edge: edges) {
        edge->getAttribute("Render", eattrib);
        eattrib->display = 1;
        JEdge::getRelations(edge, faceneighs);
        for( const JFacePtr &face : faceneighs) {
            face->getAttribute("Render", fattrib);
            fattrib->display = 1;
        }
    }

    meshViewer->updateBuffers(mesh);

    if( attribDialog.get() == nullptr )
        attribDialog.reset(new JEdgeAttributesDialog( this ));

    attribDialog->setViewManager(viewManager);
    attribDialog->setMesh( mesh );

    attribDialog->setEdges(edges);
    attribDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: makeConnections()
{
    PushButton( lookAtPushButton,      [=] {lookAt();});
    PushButton( internalPushButton,    [=] {setInternal();});
    PushButton( boundaryPushButton,    [=] {setBoundary();});
    PushButton( attribListPushButton,  [=] {openAttribListDialog();});
    PushButton( numBoundaryPushButton, [=] {countBoundEdges();});
    PushButton( numBoundaryPushButton, [=] {countBoundEdges();});
    PushButton( editPushButton,        [=] {openEditDialog();});
    PushButton( nonManifoldEdgesPushButton, [=] {getNonManifoldEdges();});

    CheckBox( displayCheckBox, [=] { checkDisplay(); });
    CheckBox( antiAliasCheckBox, [=] { checkState(); });
    CheckBox( enumCheckBox, [=] {checkState(); });

    PushButton( saveBoundaryPushButton, [=] {saveBoundary();} );
    PushButton( closePushButton, [=] {closeDialog();} );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgesDialog :: getHiddenlines()
{
    if( meshViewer == nullptr ) return;

    /*
    // How the edges have to rendered ....
    if( wireframeRadioButton->isChecked() )
         meshViewer->setEdgeMeshStyle( JMeshViewer::MESH_EDGES_LINE );
    else
         meshViewer->setEdgeMeshStyle( JMeshViewer::MESH_EDGES_HIDDENLINES_REMOVED );

    // Set the hiddenline removal method ...
    if( backfacescullRadioButton->isChecked() )
         meshViewer->setHiddenLinesMethod( JMeshViewer::HIDDENLINE_WITH_BACKFACES_CULL);

    if( frontlinesRadioButton->isChecked() )
         meshViewer->setHiddenLinesMethod( JMeshViewer::HIDDENLINE_WITH_FRONTLINES);

    // Default ...
    if( depthtestRadioButton->isChecked() )
         meshViewer->setHiddenLinesMethod( JMeshViewer::HIDDENLINE_WITH_DEPTH_TEST );

    meshViewer->refreshDisplay();
    */
}

