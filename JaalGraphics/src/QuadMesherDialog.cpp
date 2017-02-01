#include "QuadMesherDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JQuadMesherDialog :: JQuadMesherDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JQuadMesherDialog :: ~JQuadMesherDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");

    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    setMesh( meshViewer->getCurrentMesh());
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;

    if( mesh == nullptr) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    QMessageBox msg;
    int topDim = mesh->getTopology()->getDimension();
    if( topDim > 2) {
        cout << "Warning: Topological Dimension " << topDim << endl;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Only surface mesh can be converted into a quad mesh");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    /*
        double elen = mesh->getGeometry()->getMeanBoundaryEdgeLength();
        desiredEdgeLengthLineEdit->setText( QString::number(elen));
    */
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: showEvent( QShowEvent *event)
{
    if( meshViewer ) setMesh( meshViewer->getCurrentMesh() );
    QDialog::showEvent(event);
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: refineOneBoundEdge( const JEdgePtr &edge)
{
    JFaceSequence neighs;
    JEdge::getRelations(edge, neighs);

    if( neighs.size() != 1) return;

    int elemType = neighs[0]->getSize(0);

    JTriRefiner trefiner;
    switch(elemType)
    {
    case 3:
        trefiner.setMesh(mesh);
        trefiner.refine( neighs[0], edge);
        break;
    case 4:
//     AllQuadMeshGenerator::SplitBoundQuad2QuadTriangle(mesh, edge);
        break;
    }

    mesh->pruneAll();

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: displayNonQuads()
{
    JColor  highlightColor;
    highlightColor[0] = 0.5;
    highlightColor[1] = 0.1;
    highlightColor[2] = 0.1;
    highlightColor[3] = 0.0;

    size_t numfaces = mesh->getSize(2);
    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face= mesh->getFaceAt(i);
        int err = face->getAttribute("Render", fAttrib);
        if( !err) {
            fAttrib->display   =  0;
            if( face->getSize(0) != 4 ) {
                fAttrib->color     =  highlightColor;
                fAttrib->display   =  1;
            }
        }
    }

    meshViewer->updateBuffers(mesh);
}

////////////////////////////////////////////////////////////////////////////////k

void JQuadMesherDialog :: displayStructure()
{
    if( mesh == nullptr ) return ;

    cout << "There is some bug here " << endl;
    exit(0);

/*
    JMotorcycleGraph  mg;
    mg.setMesh(mesh);
    mg.getPartitions();

    JNodeSequence nodes = mg.getJunctions();

    JColor  highlightColor;
    highlightColor[0] = 1.0;
    highlightColor[1] = 0.0;
    highlightColor[2] = 0.0;
    highlightColor[3] = 0.0;

    JNodeRenderPtr nAttrib;
    for( const JNodePtr &vtx : nodes) {
        vtx->getAttribute("Render", nAttrib);
        nAttrib->glyph   =  1;
        nAttrib->display =  1;
        nAttrib->color   =  highlightColor;
    }

    JEdgeSequence edges = mg.getAllPaths();
    highlightColor[0] = 0.0;
    highlightColor[1] = 0.0;
    highlightColor[2] = 1.0;
    highlightColor[3] = 0.0;

    JEdgeRenderPtr eAttrib;
    for( const JEdgePtr &edge : edges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->glyph   =  0;
        eAttrib->display =  1;
        eAttrib->color   =  highlightColor;
    }

    meshViewer->updateBuffers(mesh);
*/
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: genmesh()
{
    QString qstr = tri2QuadsComboBox->currentText();
    string str = StdString(qstr);

    if( str == "Simple3Quads")     simpleTri2Quad();
    if( str == "HamiltonQuads")    getHamiltonQuads();
    if( str == "EdmondMatching")   getEdmondMatching();
    if( str == "GreedyMatching")   greedyMatching();
    if( str == "DualTreeMatching") treeMatching();
}
///////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: simpleTri2Quad()
{
    if( meshViewer == nullptr) return;

    string name = StdString(objectNameLineEdit->text());
    JMeshPtr trimesh = meshViewer->getMesh(name);
    if( trimesh == nullptr ) return;

    int nelm = trimesh->getTopology()->getElementsType(2);
    if( nelm != JFace::TRIANGLE ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("All triangle mesh required for generating quad mesh ");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    AllQuadMeshGenerator qmesher;
    qmesher.setMesh(trimesh);
    JMeshPtr quadmesh = qmesher.getSimpleTris2Quads();
    if( quadmesh ) {
        meshViewer->addObject(quadmesh);
        trimesh->clearAll();
        meshViewer->removeObject(trimesh);
        mesh = quadmesh;
    }

    waitCursor.stop();
}

/////////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: greedyMatching()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    AllQuadMeshGenerator  qmesher;
    qmesher.setMesh(mesh);

    cout << "Gereedy " << endl;

    JMeshPtr quadmesh = qmesher.getTrianglesMatching(AllQuadMeshGenerator::GREEDY_MATCHING);
    meshViewer->removeObject(mesh);
    meshViewer->addObject( quadmesh );
    mesh = quadmesh;
    cout << "Gereedy " << endl;
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: getHamiltonQuads()
{
    if( meshViewer == nullptr) return;

    string name = StdString(objectNameLineEdit->text());
    JMeshPtr trimesh = meshViewer->getMesh(name);
    if( trimesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    AllQuadMeshGenerator  qmesher;
    qmesher.setMesh(trimesh);

    JMeshPtr quadmesh = qmesher.getHamiltonianQuads();
    meshViewer->removeObject(trimesh);
    meshViewer->addObject( quadmesh );
    mesh = quadmesh;
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: treeMatching()
{
    string name = StdString(objectNameLineEdit->text());
    JMeshPtr trimesh = meshViewer->getMesh(name);
    if( trimesh == nullptr ) return;
    if( trimesh == nullptr ) return;

    int topDim = trimesh->getTopology()->getDimension();
    if( topDim != 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Only surface mesh can be converted into a quad mesh");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    int nelm = trimesh->getTopology()->getElementsType(2);
    if( nelm != JFace::TRIANGLE ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("All triangle mesh required for generating quad mesh ");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    /*
        Mesh *quadmesh = Jaal::AllQuadMeshGenerator::BinaryTreeMatching(trimesh, 0);
        if( quadmesh ) meshViewer->setNewMesh(quadmesh);
    */
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: getEdmondMatching()
{
    string name = StdString(objectNameLineEdit->text());
    JMeshPtr trimesh = meshViewer->getMesh(name);
    if( trimesh == nullptr ) return;

    int topDim = trimesh->getTopology()->getDimension();
    if( topDim != 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Only surface mesh can be converted into a quad mesh");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    int nelm = trimesh->getTopology()->getElementsType(2);
    if( nelm != JFace::TRIANGLE ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("All triangle mesh required for generating quad mesh ");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();
    AllQuadMeshGenerator qmesher;
    qmesher.setMesh(trimesh);
    JMeshPtr quadmesh = qmesher.getTrianglesMatching(AllQuadMeshGenerator::EDMONDS_MATCHING);
    if( quadmesh ) {
        trimesh->setActiveBit(0);
        meshViewer->removeObject(trimesh);
        meshViewer->addObject(quadmesh);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: blossumMatching()
{
    QMessageBox msg;
    msg.setIcon(QMessageBox::Warning);
    msg.setText("Not yet integrated ");
    msg.setStandardButtons( QMessageBox::Ok);
    int ret = msg.exec();
    if( ret == QMessageBox::Ok ) return;
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openMSTQuadMesherDialog()
{
    if( mesh == nullptr) return;


    int elemType  = mesh->getTopology()->getElementsType(2);
    if( elemType != JFace::QUADRILATERAL) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("This module requires pure quad mesh");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        displayNonQuads();
        return;
    }

    if( mstQuadMesherDialog == nullptr)
        mstQuadMesherDialog.reset( new JMSTQuadMesherDialog(this));

    mstQuadMesherDialog->setViewManager( viewManager);
    mstQuadMesherDialog->setMesh( mesh);
    mstQuadMesherDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openAlphaMSTDialog()
{
    if( alphaMSTDialog == nullptr)
        alphaMSTDialog.reset( new JAlphaMSTQuadMeshDialog(this));

    alphaMSTDialog->setViewManager( viewManager);
    alphaMSTDialog->setMesh( mesh);
    alphaMSTDialog->show();
    this->hide();

}
///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openPolygonSimplifyDialog()
{
    /*
        if( polySimplifyDialog == nullptr)
            polySimplifyDialog.reset( new JPolygonSimplifyDialog(this));

        polySimplifyDialog->setViewManager( viewManager);
        polySimplifyDialog->setMesh( mesh);
        this->hide();
        polySimplifyDialog->show();
    */
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openGmsh2DDialog()
{
    if( gmsh2DDialog == nullptr)
        gmsh2DDialog.reset( new JGmsh2DDialog(this));

    gmsh2DDialog->setViewManager( viewManager);
    gmsh2DDialog->setMesh( mesh);
    gmsh2DDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: getCheckerBoardPattern()
{
    if( mesh == nullptr) return;
    JQuadCheckerBoard checker;
    checker.setMesh(mesh);
    checker.genPattern();

    JColor  whiteColor, blackColor;
    whiteColor[0] = 1.0;
    whiteColor[1] = 1.0;
    whiteColor[2] = 1.0;
    whiteColor[3] = 1.0;

    blackColor[0] = 0.0;
    blackColor[1] = 0.0;
    blackColor[2] = 0.0;
    blackColor[3] = 1.0;

    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    char val;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int err = face->getAttribute("CheckerBoard", val);
        if( !err ) {
            face->getAttribute("Render", fAttrib);
            if( val == 'W') fAttrib->color = whiteColor;
            if( val == 'B') fAttrib->color = blackColor;
        }
    }
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openCleanupDialog()
{
    if( cleanupDialog == nullptr)
        cleanupDialog.reset( new JQuadMeshCleanupDialog(this));

    cleanupDialog->setViewManager( viewManager);
    cleanupDialog->setMesh( mesh);
    cleanupDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openPureQuadsDialog()
{
    if( pureQuadsDialog == nullptr)
        pureQuadsDialog.reset( new JQuadDominant2PureQuadsDialog(this));

    pureQuadsDialog->setViewManager( viewManager);
    pureQuadsDialog->setMesh( mesh);
    pureQuadsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openQualityDialog()
{
    if( quadVerdictDialog == nullptr)
        quadVerdictDialog.reset( new JQuadVerdictDialog(this));

    quadVerdictDialog->setViewManager( viewManager);
    quadVerdictDialog->setMesh( mesh);
    quadVerdictDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: openCrossFieldDialog()
{
    if( crossFieldDialog == nullptr)
        crossFieldDialog.reset( new JCrossField2DDialog(this));

    crossFieldDialog->setViewManager( viewManager);
    crossFieldDialog->setMesh( mesh);
    crossFieldDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: openSingularityGraphDialog()
{
    if( meshSingularityGraphDialog == nullptr)
        meshSingularityGraphDialog.reset( new JMeshSingularityGraphDialog(this));

    meshSingularityGraphDialog->setViewManager( viewManager);
    meshSingularityGraphDialog->setMesh( mesh);
    meshSingularityGraphDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: displaySingularNodes()
{
    if( mesh == nullptr) return;

    mesh->buildRelations(0,2);

    JNodeRenderPtr vAttrib;

    JColor redColor  = JEntityColor::getColor("Red");
    JColor blueColor = JEntityColor::getColor("Blue");

    size_t numnodes = mesh->getSize(0);
    for(size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", vAttrib);
            vAttrib->display = 1;
            vAttrib->glyph   = 0;
            int nd = vtx->getNumRelations(2);
            if( nd != 4) {
                vAttrib->glyph   = 1;
                if( nd < 4)
                    vAttrib->color = redColor;
                else
                    vAttrib->color = blueColor;
            }
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: getCyclicQuads()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JCyclicQuadMeshOptimizer  cyclicQuadOpt;
    cyclicQuadOpt.setMesh(mesh);

    bool val = preserveBoundariesCheckBox->isChecked();
    cyclicQuadOpt.setBoundaryPreserve(val);
    cyclicQuadOpt.smoothAll();

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: getBaseQuadMesh()
{
    QString qstr = baseQuadsComboBox->currentText();
    string  str  = StdString(qstr);

    AllQuadMeshGenerator allquads;
    allquads.setMesh(mesh);

    JPolyPartitioner polyparts;
    polyparts.setMesh(mesh);

    JMeshPtr basemesh;

    JWaitCursor waitCursor;
    waitCursor.start();

    if( str == "Tri2Quads") basemesh = allquads.getBaseQuadMesh(0);
    if( str == "TriMatch")  basemesh = allquads.getBaseQuadMesh(1);
    if( str == "Skeleton")  basemesh = allquads.getBaseQuadMesh(2);
    if( str == "ConvexPartitions")  basemesh =  polyparts.getPartitions();
    if( str == "EdgeMatch")  basemesh =  allquads.getBaseQuadMesh(4);

    if( basemesh == nullptr) return;

    mesh->clearAll();
    mesh->setActiveBit(0);
    meshViewer->addObject( basemesh );
    mesh = basemesh;
}

///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: refineAll()
{
/*
    JWaitCursor waitCursor;
    waitCursor.start();

    JMeshMinSingularity refine;
    refine.setMesh(mesh);

    QString qstr = desiredEdgeLengthLineEdit->text();
    double  elen = qstr.toDouble();
    refine.setEdgeLength( elen );

    JMeshPtr quadmesh = refine.refineAll();

    if( quadmesh == nullptr) return;
    meshViewer->addObject( quadmesh );
    mesh = quadmesh;
*/
}
///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openStructMeshDialog()
{
    if( structmeshDialog == nullptr)
        structmeshDialog.reset( new JStructuredMeshDialog(this));

    structmeshDialog->setViewManager( viewManager);
    this->hide();
    structmeshDialog->show();

}
///////////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: openRingQuadsDialog()
{
    if( ringQuadsDialog == nullptr)
        ringQuadsDialog.reset( new JRingQuadsDialog(this));

    ringQuadsDialog->setViewManager( viewManager);
    this->hide();
    ringQuadsDialog->show();
}
///////////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openMeshExtrudeDialog()
{
    if( meshExtrudeDialog == nullptr)
        meshExtrudeDialog.reset( new JMeshExtrudeDialog(this));

    meshExtrudeDialog->setViewManager( viewManager);
    meshExtrudeDialog->setMesh( mesh );

    meshExtrudeDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: openInstantMeshDialog()
{
    if( instantMeshDialog == nullptr)
        instantMeshDialog.reset( new JInstantMeshDialog(this));

    instantMeshDialog->setViewManager( viewManager);
    instantMeshDialog->setMeshType(4);
    instantMeshDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////////

void JQuadMesherDialog :: closeDialog()
{
    this->close();
    parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////////


void JQuadMesherDialog :: makeConnections()
{
//    PushButton( quadEditingPushButton, [=] {openQuadEditingDialog();});

    PushButton( genStructMeshPushButton, [=] {openStructMeshDialog();});
    PushButton( gmshPushButton,        [=] {openGmsh2DDialog();});
    PushButton( qualityPushButton,     [=] {openQualityDialog();});
    PushButton( meshcleanupPushButton, [=] {openCleanupDialog();});
    PushButton( pureQuadsPushButton,   [=] {openPureQuadsDialog();});
    PushButton( simplifyPushButton,    [=] {openPolygonSimplifyDialog();});
    PushButton( alphaMSTPushButton, [=] {openAlphaMSTDialog();});
    PushButton( crossFieldPushButton,  [=] {openCrossFieldDialog();});
    PushButton( mstQuadMesherPushButton,[=] {openMSTQuadMesherDialog();});
    PushButton( ringQuadsPushButton,   [=] {openRingQuadsDialog();});
    PushButton( instantMeshPushButton, [=] {openInstantMeshDialog();});
    PushButton( extrudePushButton,     [=] {openMeshExtrudeDialog();});
    PushButton( getBaseQuadsPushButton,  [=] {getBaseQuadMesh();});
//    PushButton( mstRefinePushButton,     [=] {refineAll();});
    PushButton( checkerboardPushButton,   [=] {getCheckerBoardPattern();});
    PushButton( singularGraphPushButton,  [=] {openSingularityGraphDialog();});
    PushButton( cyclicQuadsPushButton,   [=] {getCyclicQuads();});
    PushButton( applyPushButton,     [=] {genmesh();});
    PushButton( closePushButton,    [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
/*
void JQuadMesherDialog :: extractQuads()
{
    if( mesh == nullptr) return;

    JMeshPtr qmesh;
    mesh->getAttribute("QuadMesh", qmesh);
    if( qmesh ) meshViewer->removeObject(qmesh);

    JWaitCursor waitCursor;
    waitCursor.start();

    int scale = uvScaleSpinBox->value();
    JQuadExtractor  allquads;
    allquads.setMesh(mesh);
    allquads.setUVScale(scale);

    qmesh = allquads.getQuadMesh();
    int nCount = allquads.getNumIntegerUV();

    integerUVLineEdit->setText( QString::number(nCount) );

    if( qmesh ) {
        mesh->setAttribute("QuadMesh", qmesh);
        mesh->setActiveBit(0);
        meshViewer->addObject(qmesh);
    }
}
void JQuadMesherDialog :: openQuadEditingDialog()
{
    if( quadVerdictDialog == nullptr)
        quadVerdictDialog.reset( new JQuadVerdictDialog(this));

    quadVerdictDialog->setViewManager( viewManager);
    quadVerdictDialog->setMesh( mesh);
    this->hide();
    quadVerdictDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMesherDialog :: openSketchQuadsDialog()
{
    if( sketchQuadsDialog == nullptr)
        sketchQuadsDialog.reset( new JSketchQuadsDialog(this));

    sketchQuadsDialog->setViewManager( viewManager);
    sketchQuadsDialog->setMesh( mesh);
    sketchQuadsDialog->show();
    this->hide();
}
void JQuadMesherDialog :: openMixedIntegerQuadsDialog()
{
    if( mixedIntegerQuadsDialog == nullptr)
        mixedIntegerQuadsDialog.reset( new JMixedIntegerQuadsDialog(this));

    mixedIntegerQuadsDialog->setViewManager( viewManager);
    mixedIntegerQuadsDialog->setMesh( mesh);
    this->hide();
    mixedIntegerQuadsDialog->show();
}
*/
