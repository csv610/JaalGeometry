#include "MeshToolsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshToolsDialog :: JMeshToolsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshToolsDialog :: ~JMeshToolsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: showEvent( QShowEvent *e)
{
    if( meshViewer == nullptr ) return;

    mesh = meshViewer->getCurrentMesh();
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );

    QDialog::showEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: setMesh()
{
    if( mesh == nullptr) return;


    bool val;
    val = primalMeshCheckBox->isChecked();
    mesh->setActiveBit(val);

    JMeshPtr dmesh;
    mesh->getAttribute("DualGraph", dmesh);
    if( dmesh ) {
        val = dualGraphCheckBox->isChecked();
        dmesh->setActiveBit(val);
        if( val ) meshViewer->addObject(dmesh);
    }

    JMeshPtr texmesh;
    mesh->getAttribute("TextureMesh", texmesh);
    if( texmesh ) {
        val = textureMeshCheckBox->isChecked();
        texmesh->setActiveBit(val);
        if( val ) meshViewer->addObject(texmesh);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openPickEntityDialog()
{

    if( pickEntityDialog.get() == nullptr )
        pickEntityDialog.reset(new JMeshEntityPickerDialog(this));

    pickEntityDialog->setViewManager( viewManager );
    pickEntityDialog->setMesh( mesh );
    pickEntityDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: warnMessage()
{
    QMessageBox msg;
    msg.setIcon(QMessageBox::Warning);
    msg.setText("No mesh object selected: Select one from the list");
    msg.setStandardButtons( QMessageBox::Ok);
    int ret = msg.exec();
    if( ret == QMessageBox::Ok ) return;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->getViewManager()->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openTriMesherDialog()
{
    if( triMesherDialog == nullptr )
        triMesherDialog.reset( new JTriMesherDialog(this) );

    meshViewer->setCurrentMesh(mesh);
    triMesherDialog->setViewManager( viewManager );
    triMesherDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openQuadMesherDialog()
{
    if( quadMesherDialog.get() == nullptr )
        quadMesherDialog.reset(new JQuadMesherDialog(this));

    meshViewer->setCurrentMesh(mesh);
    quadMesherDialog->setViewManager( viewManager );
    quadMesherDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openTetMesherDialog()
{
    if( tetMesherDialog.get() == nullptr )
        tetMesherDialog.reset(new JTetMesherDialog(this));

    tetMesherDialog->setViewManager( viewManager );
    tetMesherDialog->setMesh(mesh);
    tetMesherDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openHexMesherDialog()
{
    if( hexMesherDialog == nullptr )
        hexMesherDialog.reset(new JHexMesherDialog(this));

    hexMesherDialog->setViewManager( viewManager );
    hexMesherDialog->setMesh(mesh);
    hexMesherDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openMeshRenderDialog()
{
    if( meshRenderDialog == nullptr )
        meshRenderDialog.reset(new JMeshRenderDialog(this));

    meshRenderDialog->setViewManager( viewManager );
    meshRenderDialog->setMesh(mesh);
    meshRenderDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openAffineDialog()
{
    if( meshViewer == nullptr) return;

    if( affineDialog.get() == nullptr )
        affineDialog.reset(new JMeshAffineTransformsDialog(this));

    affineDialog->setViewManager( viewManager );
    affineDialog->setMesh(mesh);
    affineDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openRelationsDialog()
{
    if( meshViewer == nullptr) return;

    if( relationsDialog.get() == nullptr )
        relationsDialog.reset(new JMeshRelationsTableDialog(this));

    relationsDialog->setMesh( mesh );
    relationsDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openMeshGeomQualityDialog()
{
    if( meshGeomQualityDialog.get() == nullptr )
        meshGeomQualityDialog.reset(new JMeshGeometricQualityDialog(this));

    meshGeomQualityDialog->setViewManager( viewManager );
    meshGeomQualityDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openGenSimpleShapeDialog()
{
    if( genSimpleShapeDialog.get() == nullptr)
        genSimpleShapeDialog.reset(new JGenSimpleShapeDialog(this));

    genSimpleShapeDialog->setViewManager( viewManager );
    genSimpleShapeDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openTopoQueryDialog()
{
    if( meshViewer == nullptr) return;

    if( topoQueryDialog.get() == nullptr )
        topoQueryDialog.reset(new JMeshTopologyQueryDialog(this));

    topoQueryDialog->setViewManager( viewManager );
    topoQueryDialog->setMesh(mesh);
    topoQueryDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshInterpolationDialog()
{
    if( meshInterpolationDialog == nullptr )
        meshInterpolationDialog.reset(new JMeshInterpolationDialog(this));

    meshInterpolationDialog->setViewManager( viewManager );
    meshInterpolationDialog->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openSurfVecFieldDialog()
{
    if( surfVecFieldDialog == nullptr )
        surfVecFieldDialog.reset(new JMeshSurfaceVectorFieldDialog(this));

    surfVecFieldDialog->setViewManager( viewManager );
    surfVecFieldDialog->setMesh( mesh);
    surfVecFieldDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshDeformationDialog()
{
    if( meshDeformationDialog.get() == nullptr ) {
        meshDeformationDialog.reset(new JMeshDeformationDialog(this));
    }

    meshDeformationDialog->setViewManager( viewManager );
    meshDeformationDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshSegmentationDialog()
{
    if( segmentationDialog.get() == nullptr )
        segmentationDialog.reset(new JMeshSegmentationDialog(this));

    segmentationDialog->setViewManager( viewManager );
    segmentationDialog->setMesh(mesh);
    segmentationDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshBooleanDialog()
{
    if( meshBooleanDialog == nullptr )
        meshBooleanDialog.reset(new JMeshBooleanDialog(this));

    meshBooleanDialog->setViewManager( viewManager );
    meshBooleanDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openSurfaceParameterizationDialog()
{
    if( surfaceParameterizationDialog == nullptr )
        surfaceParameterizationDialog.reset(new JSurfaceParameterizationDialog(this));

    surfaceParameterizationDialog->setViewManager( viewManager );
    surfaceParameterizationDialog->setMesh( mesh );
    surfaceParameterizationDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openMeshFeaturesDialog()
{
    if( meshFeaturesDialog.get() == nullptr )
        meshFeaturesDialog.reset(new JMeshFeaturesDialog(this));

    meshFeaturesDialog->setViewManager( viewManager );
    meshFeaturesDialog->setMesh(mesh);
    meshFeaturesDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshGeodesicsDialog()
{
    if( meshGeodesicsDialog == nullptr )
        meshGeodesicsDialog.reset(new JMeshGeodesicsDialog(this));

    meshGeodesicsDialog->setViewManager( viewManager );
    meshGeodesicsDialog->setMesh(mesh);
    meshGeodesicsDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openMeshMeanCurvatureFlowDialog()
{
    if( meanCurvatureFlowDialog == nullptr )
        meanCurvatureFlowDialog.reset(new JMeshMeanCurvatureFlowDialog(this));

    meanCurvatureFlowDialog->setViewManager( viewManager );
    meanCurvatureFlowDialog->setMesh(mesh);
    meanCurvatureFlowDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: closeDialog()
{
    this->close();
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: selectMesh( QModelIndex index)
{
    int  row = index.row();
    mesh = meshViewer->getMesh(row);
    JBoundingBox box;
    mesh->getAttribute("AxisBoundingBox", box);
    viewManager->setCenter(box.getCenter());
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openSuggestiveContoursDialog()
{
    if(suggestiveContoursDialog == nullptr)
        suggestiveContoursDialog.reset( new JSuggestiveContoursDialog(this));

    suggestiveContoursDialog->setViewManager( viewManager );
    suggestiveContoursDialog->setMesh( mesh );
    suggestiveContoursDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openMagnifyingLensDialog()
{
    if(magnifyingLensDialog == nullptr)
        magnifyingLensDialog.reset( new JMagnifyingLensDialog(this));

    magnifyingLensDialog->setViewManager( viewManager );
    magnifyingLensDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openMeshDualGrapherDialog()
{
    if( dualMeshDialog == nullptr)
        dualMeshDialog.reset( new JMeshDualGrapherDialog(this));

    dualMeshDialog->setViewManager( viewManager );
    dualMeshDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openMeshSlicerDialog()
{
    if( meshSlicerDialog == nullptr)
        meshSlicerDialog.reset(new JMeshSlicerDialog(this));

    meshSlicerDialog->setViewManager( viewManager );
    meshSlicerDialog->setMesh(mesh);

    meshSlicerDialog->show();
    this->hide();

}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: resetCamera()
{
    if( mesh == nullptr) return;
    meshViewer->lookAt(mesh);

    /*
        viewManager->resetView( JViewDirection:: FRONT_VIEW);

        JBoundingBox box = mesh->getGeometry()->getBoundingBox();
        Point3D  pC   = box.getCenter();
        double   len  = box.getMaxLength();

        qglviewer::Vec  vec;
        vec[0] = pC[0];
        vec[1] = pC[1];
        vec[2] = pC[2] + 2.0*len;
        viewManager->camera()->setPosition(vec);

        vec[0] = pC[0];
        vec[1] = pC[1];
        vec[2] = pC[2];
        viewManager->camera()->lookAt(vec);
        viewManager->camera()->setSceneCenter(vec);

        vec[0] = 0.0;
        vec[1] = 1;
        vec[2] = 0;
        viewManager->camera()->setUpVector(vec);
        meshViewer->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: ambientOcclusion()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();
  
    cout << " Triagulated first " << endl;
    exit(0);

    JMeshEigenMatrix mat;
    mat.setMesh(mesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    MatrixXd N;
    igl::per_vertex_normals(V,F,N);

    Eigen::VectorXd AO;
    // Compute ambient occlusion factor using embree
    igl::embree::ambient_occlusion(V,F,V,N,500,AO);
    AO = 1.0 - AO.array();

    const RowVector3d color(0.9,0.85,0.9);
    MatrixXd C = color.replicate(V.rows(),1);
    for (unsigned i=0; i<C.rows(); ++i) {
        C.row(i) *= AO(i); //std::min<double>(AO(i)+0.2,1);
    }

    JNodeRenderPtr vAttrib;
    size_t numnodes = mesh->getSize(0);
    JColor rgb;
    size_t index = 0;
    for( size_t i = 0; i <  numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", vAttrib);
            rgb[0] = C.coeff(index,0);
            rgb[1] = C.coeff(index,1);
            rgb[2] = C.coeff(index,2);
            vAttrib->color = rgb;
            index++;
        }
    }

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->setSurfaceShade(JRender::SMOOTH_SHADE);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshSpectrumDialog()
{
    if( meshSpectrumDialog == nullptr )
        meshSpectrumDialog.reset(new JMeshSpectrumDialog(this));

    meshSpectrumDialog->setViewManager( viewManager );
    meshSpectrumDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshSkeletonDialog()
{
    if( meshSkeletonDialog == nullptr )
        meshSkeletonDialog.reset(new JMeshSkeletonDialog(this));

    meshSkeletonDialog->setViewManager(viewManager);
    meshSkeletonDialog->setMesh(mesh);
    meshSkeletonDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshOptDialog()
{
    if( meshOptDialog == nullptr )
        meshOptDialog.reset(new JMeshOptDialog(this));

    meshOptDialog->setViewManager( viewManager );
    meshOptDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshDualDialog()
{
    if( meshdualsDialog.get() == nullptr )
        meshdualsDialog.reset(new JMeshDualsDialog(this));

    meshdualsDialog->setViewManager( viewManager );
    meshdualsDialog->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshTopoQualityDialog()
{
    if( meshTopoQualityDialog == nullptr )
        meshTopoQualityDialog.reset(new JMeshTopologyQualityDialog(this));

    meshTopoQualityDialog->setViewManager( viewManager );
    meshTopoQualityDialog->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshRefine2DDialog()
{
    if( refine2dDialog == nullptr )
        refine2dDialog.reset(new JMeshRefine2DDialog(this));

    refine2dDialog->setViewManager( viewManager );
    refine2dDialog->setMesh(mesh);
    refine2dDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshRefine3DDialog()
{
    if( refine3dDialog == nullptr )
        refine3dDialog.reset(new JMeshRefine3DDialog(this));

    refine3dDialog->setViewManager( viewManager );
    refine3dDialog->setMesh(mesh);
    refine3dDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openMeshPartitionDialog()
{
    if( meshpartitionDialog == nullptr )
        meshpartitionDialog.reset(new JMeshPartitionDialog(this));

    meshpartitionDialog->setViewManager( viewManager );
    meshpartitionDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: openContourEditingDialog()
{
    if( contourEditingDialog == nullptr )
        contourEditingDialog.reset(new JContourEditingDialog(this));

    contourEditingDialog->setViewManager( viewManager );
    contourEditingDialog->setMesh(mesh);
    contourEditingDialog->show();
    this->hide();
}


/*
///////////////////////////////////////////////////////////////////////////////
void JaalMainWindow :: openPatchRemeshDialog()
{
        if(patchQuadmeshingDialog.get() == nullptr )
            patchQuadmeshingDialog.reset(new JPatchQuadmeshingDialog(this));

        patchQuadmeshingDialog->setViewManager( viewer );
        patchQuadmeshingDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openTetMesherDialog()
{
    if( tetmesherDialog.get() == nullptr )
        tetmesherDialog.reset(new JTetMesherDialog(this));

    tetmesherDialog->setViewManager( viewer );
    tetmesherDialog->show();
}
///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openTriMesherDialog()
{
    if( trimesherDialog.get() == nullptr )
        trimesherDialog.reset( new JTriMesherDialog(this) );

    trimesherDialog->setViewManager( viewer );
    trimesherDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openSaveAnimationDialog()
{
    if( saveAnimationDialog.get() == nullptr )
        saveAnimationDialog.reset(new JSaveAnimationDialog());

    saveAnimationDialog->setViewManager( viewer );
    saveAnimationDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
void JaalMainWindow :: openMorseAnalysisDialog()
{
    if( morseAnalysisDialog.get() == nullptr )
        morseAnalysisDialog.reset(new JMorseAnalysisDialog());

    morseAnalysisDialog->setViewManager( viewer);
    morseAnalysisDialog->show();
}
///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openTrimeshCleanupDialog()
{
    if( trimeshCleanupDialog.get() == nullptr )
        trimeshCleanupDialog.reset(new JTrimeshCleanupDialog());

    trimeshCleanupDialog->setViewManager( viewer);
    trimeshCleanupDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JaalMainWindow :: openTriSimplifyMeshDialog()
{
    if( triSimplifyDialog.get() == nullptr )
        triSimplifyDialog.reset(new JTriSimplificationDialog());

    triSimplifyDialog->setViewManager( viewer);
    triSimplifyDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JaalMainWindow :: openHeatConductionDialog()
{
    if( heatConductionDialog == nullptr )
        heatConductionDialog.reset(new JHeatConductionDialog(this));

    heatConductionDialog->setViewManager( viewer);
    heatConductionDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openMeshTangleDialog()
{
    if( meshTangleDialog.get() == nullptr )
        meshTangleDialog.reset(new JMeshTangleDialog(this));

    meshTangleDialog->setViewManager( viewer);
    meshTangleDialog->show();

}

///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openObjectsListDialog()
{
    if( objectsListDialog.get() == nullptr )
        objectsListDialog.reset(new JObjectsListDialog(this));

    objectsListDialog->setViewManager( viewer);
    objectsListDialog->setType(0);
    objectsListDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openQuadMesherDialog()
{
    if( quadmesherDialog.get() == nullptr )
        quadmesherDialog.reset(new JQuadMesherDialog(this));

    quadmesherDialog->setViewManager( viewer);
    quadmesherDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JaalMainWindow :: openTangleFEMTestsDialog()
{
    if( tanglefemtestsDialog == nullptr )
        tanglefemtestsDialog.reset(new JTangleFEMTestsDialog(this));

    tanglefemtestsDialog->setViewManager( viewer);
    tanglefemtestsDialog->show();
}
*/
///////////////////////////////////////////////////////////////////////////////
/*
void JaalMainWindow :: openMeshContoursDialog()
{
    if( meshContoursDialog == nullptr )
        meshContoursDialog.reset(new JMeshContoursDialog(this));

    meshContoursDialog->setViewManager( viewer);
    meshContoursDialog->show();
}
*/
///////////////////////////////////////////////////////////////////////////////

/*

void JaalMainWindow :: openMeshSpectrumDialog()
{
    if( meshSpectrumDialog == nullptr )
        meshSpectrumDialog.reset(new JMeshSpectrumDialog(this));

    meshSpectrumDialog->setViewManager( viewer);
    meshSpectrumDialog->show();
}

*/
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: makeConnections()
{
    PushButton( genTriMeshPushButton,  [=] {openTriMesherDialog();});
    PushButton( genQuadMeshPushButton, [=] {openQuadMesherDialog();});
    PushButton( genTetMeshPushButton,  [=] {openTetMesherDialog();});
    PushButton( genHexMeshPushButton,  [=] {openHexMesherDialog();});

    PushButton( contourEditingPushButton,  [=] {openContourEditingDialog();});
    PushButton( pickPushButton,       [=] {openPickEntityDialog();});
    PushButton( relationsPushButton,  [=] {openRelationsDialog();});
    PushButton( dualMeshPushButton,   [=] {openMeshDualGrapherDialog();});
    PushButton( topologyQueryPushButton, [=] {openTopoQueryDialog();});
    PushButton( magnifyingLensPushButton, [=] {openMagnifyingLensDialog();});
    PushButton( affineTransformPushButton, [=] {openAffineDialog();});
    PushButton( ambientOcclusionPushButton, [=] {ambientOcclusion();});
    PushButton( meshSlicerPushButton, [=] {openMeshSlicerDialog();});
    PushButton( suggestiveContoursPushButton, [=] {openSuggestiveContoursDialog();});
    PushButton( meshSkeletonPushButton, [=] {openMeshSkeletonDialog();});
    PushButton( meshRenderPushButton, [=] {openMeshRenderDialog();});
    PushButton( meshPartitionPushButton, [=] {openMeshPartitionDialog();});
    PushButton( meshSegmentationPushButton, [=] {openMeshSegmentationDialog();});
    PushButton( surfParameterizationPushButton, [=] {openSurfaceParameterizationDialog();});
    PushButton( meshInterpolationPushButton, [=] {openMeshInterpolationDialog();});
    PushButton( meshBooleanPushButton, [=] {openMeshBooleanDialog();});
    PushButton( meshGeodesicsPushButton, [=] {openMeshGeodesicsDialog();});
    PushButton( meshOptPushButton, [=] {openMeshOptDialog();});
    PushButton( surfVecFieldPushButton, [=] {openSurfVecFieldDialog();});
    PushButton( meanCurvatureFlowPushButton, [=] {openMeshMeanCurvatureFlowDialog();});
    PushButton( meshFeaturesPushButton, [=] {openMeshFeaturesDialog();});
    PushButton( meshRefine2DPushButton, [=] {openMeshRefine2DDialog();});
    PushButton( meshGeomQualityPushButton, [=] {openMeshGeomQualityDialog();});

//  PushButton( fitBoxPushButton,     [=] {fitBoundBox();});
//  PushButton( resetCameraPushButton,[=] {resetCamera();});
//  PushButton( boxColorPushButton,   [=] {openBoxColorDialog();});
//  PushButton( fitSpherePushButton, [=] {fitBoundSphere();});
//  ComboBox( boundingComboBox, [=] { setEnclosure(); });
//  CheckBox( boundingCheckBox,    [=] {setEnclosure();});

    CheckBox( primalMeshCheckBox,  [=] {setMesh();});
    CheckBox( dualGraphCheckBox,   [=] {setMesh();});
    CheckBox( textureMeshCheckBox, [=] {setMesh();});

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
/*
void JMeshToolsDialog :: setListView()
{
    if( meshViewer == nullptr ) return;

    QStandardItemModel *model = new QStandardItemModel();
    int numMesh = meshViewer->getSize();
    for( int i = 0; i < numMesh; i++) {
        Mesh *m = meshViewer->getMesh(i);
        string name = m->getName();
        QStandardItem *item = new QStandardItem( QString(name.c_str()));
        item->setEditable(false);
        model->appendRow(item);
    }
    listView->setModel(model);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: removeUnrefNodes()
{
        string name = StdString(objectNameLineEdit->text());
        mesh = meshViewer->getMesh(name);

        if( mesh == nullptr) {
            warnMessage();
            return;
        }
        mesh->getTopology()->remove_unattached_nodes();
        meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: removeUnrefEntities()
{
        string name = StdString(objectNameLineEdit->text());
        mesh = meshViewer->getMesh(name);

        if( mesh == nullptr) {
            warnMessage();
            return;
        }

        mesh->getTopology()->remove_unattached_simplices();
        meshViewer->refreshDisplay();
}
void JMeshToolsDialog :: setToolsMode()
{
    if( meshViewer == nullptr ) return;

    QString qstr;
    qstr = renderModeComboBox->currentText();
    string str = StdString(qstr);

    if( str == "PointCloud")
        meshViewer->setToolsMode(JToolsMode::POINTCLOUD);

    if( str == "Wireframe")
        meshViewer->setToolsMode(JToolsMode::WIREFRAME);

    if( str == "FlatSurface")
        meshViewer->setToolsMode(JToolsMode::FLAT_SHADE);

    if( str == "SmootSurface")
        meshViewer->setToolsMode(JToolsMode::SMOOTH_SHADE);

    if( str == "Hiddenlines")
        meshViewer->setToolsMode(JToolsMode::HIDDENLINES);

    meshViewer->refreshDisplay();
}

void JMeshToolsDialog :: openMeshlistDialog()
{
    if( meshlistDialog.get() == nullptr )
        meshlistDialog.reset(new JObjectsListDialog(this));

    meshlistDialog->setViewManager( viewManager );
    meshlistDialog->setType(1);
    meshlistDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: fitBoundSphere()
{
    if( mesh == nullptr) return;

    JSphere sph = mesh->getGeometry()->getMinimumSphere();
    Point3D center = sph.getCenter();

    qglviewer::Vec c;

    c[0] = center[0];
    c[1] = center[1];
    c[2] = center[2];
    viewManager->camera()->fitSphere(c,sph.getRadius());
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: fitBoundBox()
{
    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    Point3D pmin = box.getCorner(0);
    Point3D pmax = box.getCorner(6);

    qglviewer::Vec vmin, vmax;
    vmin[0] =  pmin[0];
    vmin[1] =  pmin[1];
    vmin[2] =  pmin[2];

    vmax[0] =  pmax[0];
    vmax[1] =  pmax[1];
    vmax[2] =  pmax[2];
    viewManager->camera()->fitBoundingBox(vmin,vmax);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshToolsDialog :: openBoxColorDialog()
{
    QColor color = QColorDialog::getColor();
    float rgb[3];
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;

        JEdgeToolsPtr eAttrib;
        for( size_t i = 0; i < edges.size(); i++) {
            edges[i]->getAttribute("Tools", eAttrib);
            eAttrib->color[0] = rgb[0];
            eAttrib->color[1] = rgb[1];
            eAttrib->color[2] = rgb[2];
        }
        meshViewer->updateBuffers( mesh );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshToolsDialog :: setEnclosure()
{
    if( meshViewer == nullptr ) return;

    QString qstr;
    qstr = boundingComboBox->currentText();
    string str = qstr.toUtf8().constData();

    bool val = boundingCheckBox->isChecked();

    if( str == "AxisAlignedBox")
        meshViewer->displayEnclosure(val, JMeshViewer::AXIS_ALIGNED_BOUNDING_BOX);

    if( str == "MinBox") {
        JHexahedronPtr minBox =  mesh->getGeometry()->getMinimumBox();
        mesh->setAttribute("MinBoundingBox", minBox);
        meshViewer->displayEnclosure(val, JMeshViewer::MINIMUM_BOUNDING_BOX);
    }

    if( str == "MinSphere")
        meshViewer->displayEnclosure(val, JMeshViewer::MINIMUM_SPHERE);

    meshViewer->refreshDisplay();
}
*/

