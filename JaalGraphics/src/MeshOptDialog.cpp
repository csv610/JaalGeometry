#include "MeshOptDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshOptDialog :: JMeshOptDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    meshViewer = nullptr;
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshOptDialog :: ~JMeshOptDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshOptDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    JNodeRenderPtr attrib;
    JColor redColor;
    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 1.0;

    JColor greenColor;
    greenColor[0] = 0.0;
    greenColor[1] = 1.0;
    greenColor[2] = 0.0;
    greenColor[3] = 1.0;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", attrib);
            if( vtx->hasAttribute("Constraint")) {
                attrib->color = redColor;
                attrib->display = 1;
                attrib->scale   = 1.5;
            } else {
                attrib->color = greenColor;
                attrib->display = 1;
                attrib->scale   = 1.0;
            }
        }
    }
    mesh->getGeometry()->getCoordsArray(orgCoords,l2g);
    meshViewer->refreshDisplay();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = -1;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshOptDialog :: showEvent( QShowEvent *)
{
    if( meshViewer == nullptr) return;
    setMesh( meshViewer->getCurrentMesh() );
}
///////////////////////////////////////////////////////////////////////////////
void JMeshOptDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: setConstraints()
{
    if( mesh == nullptr ) return;

    if( constraintsDialog == nullptr )
        constraintsDialog.reset(new JMeshConstraintsDialog(this));

    constraintsDialog->setViewManager( viewManager );
    viewManager->attach(constraintsDialog.get());

    constraintsDialog->setMesh(mesh);

    constraintsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshOptDialog :: untangle()
{
    if( mesh == nullptr ) return;

    if( untangleDialog == nullptr )
        untangleDialog.reset(new JMeshUntangleDialog(this));

    untangleDialog->setViewManager( viewManager );
    untangleDialog->setMesh(mesh);

    untangleDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: optimize()
{
    if( mesh == nullptr ) return;

    size_t numFixed = 0;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            if( vtx->hasAttribute("Constraint")) numFixed++;
        }
    }

    QMessageBox msg;
    if( numFixed == 0) {
        msg.setIcon(QMessageBox::Warning);
        msg.setText(tr("There are no constraints in the mesh. Need at least one constraints"));
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok) return;
    }

    vector<double> qData;

    bool with_simplices = withSimplicesCheckBox->isChecked();

    /*
        if( !with_simplices) {
            int topDim = mesh->getTopology()->getDimension();
            if( topDim == 2 ) {
                int etype = mesh->getTopology()->getElementsType(2);
                if( etype == JFace::QUADRILATERAL ) {
                    meshViewer->getViewManager()->refreshDisplay();
                    if( mesh->getGeometry()->count_concave_faces() ) {
                        QMessageBox msg;
                        msg.setIcon(QMessageBox::Warning);
                        msg.setText(tr("Quadmesh contains concave faces (or possibly doublets), optimization will be performed "
                                       "on a triangle mesh which is temporarly constructred. Since this might have large space "
                                       "and time complexity, it is advice to use less number of iterations. Doublets are always"
                                       "problemetc, so remove them before the optimization"));
                        msg.setStandardButtons( QMessageBox::Cancel | QMessageBox::Ok);
                        int ret = msg.exec();
                        if( ret == QMessageBox::Cancel ) return;
                    }
                }
            }
        }
    */

    int norm = JMeshNonlinearOptimization::LP_NORM;
    int normVal = 2;
    if( LinfNormRadioButton->isChecked() )
        norm = JMeshNonlinearOptimization::LINF_NORM;
    else {
        normVal = pnormSpinBox->value();
    }

    int niter = numNonLinearIterationsSpinBox->value();
    bool global_patch = globalPatchRadioButton->isChecked();

    nonlinearOpt.setNumIterations(niter);
    nonlinearOpt.setNorm(norm);
    nonlinearOpt.setNormVal(normVal);
    nonlinearOpt.setPatchType(global_patch);

    JWaitCursor waitCursor;
    waitCursor.start();

    cout << "Info: Start non-linear optimization " << endl;
    nonlinearOpt.setMesh(mesh);
    nonlinearOpt.improveQuality(with_simplices);

    meshViewer->updateGeometryBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshOptDialog :: nonlinear()
{
    QString qstr;
    string str;

    qstr = nonlinearAlgoComboBox->currentText();
    str  = qstr.toUtf8().constData();

    if( str == "Quasi Newton")
        nonlinearOpt.setAlgorithm(JMeshNonlinearOptimization::QUASI_NEWTON);

    if( str == "Conjugate Gradient")
        nonlinearOpt.setAlgorithm(JMeshNonlinearOptimization::CONJUGATE_GRADIENT);

    if( str == "Feasible Newton")
        nonlinearOpt.setAlgorithm(JMeshNonlinearOptimization::FEASIBLE_NEWTON);

    if( str == "Steepest Descent")
        nonlinearOpt.setAlgorithm(JMeshNonlinearOptimization::STEEPEST_DESCENT);

    if( str == "Trust Region")
        nonlinearOpt.setAlgorithm(JMeshNonlinearOptimization::TRUST_REGION);

    if( str == "Smart Laplace")
        nonlinearOpt.setAlgorithm(JMeshNonlinearOptimization::SMART_LAPLACIAN);

    qstr = nonlinearQualityComboBox->currentText();
    str  = qstr.toUtf8().constData();

    nonlinearOpt.setQualityMetric(-1);
    if( str == "Inverse Mean Ratio" )  {
        nonlinearOpt.setQualityMetric(JMeshNonlinearOptimization::INVERSE_MEAN_RATIO);
    }

    if( str == "Mean Ratio" )  {
        nonlinearOpt.setQualityMetric(JMeshNonlinearOptimization::MEAN_RATIO);
    }

    if( str == "Condition Number" ) {
        nonlinearOpt.setQualityMetric(JMeshNonlinearOptimization::CONDITION_NUMBER);
    }

    if( str == "Edge Length" )  {
        nonlinearOpt.setQualityMetric(JMeshNonlinearOptimization::EDGE_LENGTH);
    }

    optimize();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: getOriginal()
{
    if( meshViewer == nullptr ) return;
    if( mesh ) mesh->getGeometry()->setCoordsArray(orgCoords, l2g);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: openLaplaceDialog()
{
    if( laplaceDialog.get() == nullptr )
        laplaceDialog.reset(new JMeshLaplaceSmoothingDialog(this));

    laplaceDialog->setViewManager( viewManager );
    laplaceDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: openLloydDialog()
{
    if( lloydRelaxationDialog == nullptr )
        lloydRelaxationDialog.reset(new JLloydRelaxationDialog(this));

    lloydRelaxationDialog->setViewManager( viewManager );
    lloydRelaxationDialog->setMesh( mesh );
    lloydRelaxationDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: openMeanCurvatureFlowDialog()
{
    if( meanCurvatureFlowDialog == nullptr ) {
        meanCurvatureFlowDialog.reset(new JMeshMeanCurvatureFlowDialog(this));
    }

    meanCurvatureFlowDialog->setViewManager( viewManager );
    meanCurvatureFlowDialog->setMesh( mesh );
    meanCurvatureFlowDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshOptDialog :: reparamCurves()
{
    if( mesh == nullptr) return;

    vector<JEdgeSequence> boundCurves;

    mesh->getTopology()->getBoundary(boundCurves);
    int ncurves = boundCurves.size();

    for( int i = 0; i < ncurves; i++)
        JEdgeGeometry::makeUniform( boundCurves[i] );

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: smoothCurves()
{
    if( mesh == nullptr) return;

    int niter = smoothCurvesSpinBox->value();

    vector<JEdgeSequence> boundCurves;

    mesh->getTopology()->getBoundary(boundCurves);
    int ncurves = boundCurves.size();

    for( int i = 0; i < ncurves; i++)
        JEdgeGeometry::smooth( boundCurves[i], niter);

    meshViewer->updateBuffers(mesh);

}
///////////////////////////////////////////////////////////////////////////////
void JMeshOptDialog :: openShapeOpDialog()
{
    if( shapeOpDialog == nullptr )
        shapeOpDialog.reset(new JShapeOpDialog(this));

    shapeOpDialog->setViewManager( viewManager );
    shapeOpDialog->setMesh( mesh );
    shapeOpDialog->show();
    this->hide();

}
///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: openUntangleDialog()
{
    if( untangleDialog == nullptr )
        untangleDialog.reset(new JMeshUntangleDialog(this));

    untangleDialog->setViewManager( viewManager );
    untangleDialog->setMesh( mesh );
    untangleDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: closeDialog()
{
/*
    if( meshViewer ) {
        meshViewer->displayAll(0,1);
        meshViewer->refreshDisplay();
    }
*/
    parentWidget()->show();
    close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshOptDialog :: makeConnections()
{
    PushButton( laplacianPushButton, [=] {openLaplaceDialog();});
    PushButton( voronoiOptPushButton,[=] {openLloydDialog();});
    PushButton( orgmeshPushButton,   [=] {getOriginal();});
    PushButton( untanglePushButton,  [=] {openUntangleDialog();});
    PushButton( shapeOptPushButton,  [=] {openShapeOpDialog();});
    PushButton( reparamCurvePushButton,   [=] {reparamCurves();});
    PushButton( applyNonlinearPushButton, [=] {nonlinear();});
    PushButton( setConstraintsPushButton, [=] {setConstraints();});
    PushButton( smoothBoundaryCurvesPushButton, [=] {smoothCurves();});
    PushButton( meanCurvatureFlowPushButton,   [=] {openMeanCurvatureFlowDialog();});

    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
