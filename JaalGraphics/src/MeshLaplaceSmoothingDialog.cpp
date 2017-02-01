#include "MeshLaplaceSmoothingDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshLaplaceSmoothingDialog :: JMeshLaplaceSmoothingDialog( QWidget *parent) : QDialog(parent)
{
    mesh = nullptr;
    meshViewer = nullptr;

    setupUi(this);
    makeConnections();

    lambdaLineEdit->setText( QString::number(0.6307) );
    muLineEdit->setText( QString::number(-0.6734) );

    maxResidueLineEdit->setText( QString::number(1.0E-06) );

    edges_weight_assigned = 0;
}

///////////////////////////////////////////////////////////////////////////////

JMeshLaplaceSmoothingDialog :: ~JMeshLaplaceSmoothingDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: init()
{
    if( laplace == nullptr )
        laplace.reset( new JLaplaceMeshSmoother() );

    assert( viewManager ) ;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");

    meshViewer = nullptr;
    if( c ) meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    mesh = meshViewer->getCurrentMesh();

    setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    initialized = 0;

    if(mesh ) {
        string name = mesh->getName();
        objectNameLineEdit->setText( QString(name.c_str() ) );
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: initMesh()
{
    if( mesh == nullptr || initialized == 1) return;

    laplace->setMesh(mesh);
    mesh->getGeometry()->getCoordsArray(orgCoords,l2g);
    setConstraints();

    displayNegativeElements();
    initialized = 1;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshLaplaceSmoothingDialog :: showEvent( QShowEvent *)
{
    setConstraints();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: assignColors()
{
    if( mesh == nullptr ) return;

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
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            int err = vtx->getAttribute("Render", attrib);
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
    assert( meshViewer );
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshLaplaceSmoothingDialog :: setEdgesWeight()
{
    if( meshViewer == nullptr || mesh == nullptr) return;

    QString qstr = edgeWeightComboBox->currentText();
    string  str  = qstr.toUtf8().constData();

    int edgeweight = -1;
    if( str == "Combinatorial") edgeweight = JLaplaceMeshSmoother::COMBINATORIAL;
    if( str == "Length")        edgeweight = JLaplaceMeshSmoother::EDGE_LENGTH;
    if( str == "Tutte")         edgeweight = JLaplaceMeshSmoother::TUTTE;
    if( str == "Cotangent")     edgeweight = JLaplaceMeshSmoother::COTANGENT;
    if( str == "Normalized")    edgeweight = JLaplaceMeshSmoother::NORMALIZED;
    if( str == "FloaterMeanValue") edgeweight = JLaplaceMeshSmoother::FLOATER_MEAN_VALUE;

    JWaitCursor waitCursor;
    waitCursor.start();

    bool update = updateWeightCheckBox->isChecked();
    laplace->setEdgesWeight( edgeweight, update);

    edges_weight_assigned = 1;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: smooth()
{
    if( meshViewer == nullptr || mesh == nullptr) return;

    size_t nCount = mesh->getNumAttributes("Constraint", 0);
    if( nCount == 0) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Input Mesh does not have any constraint ");
        msg.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
        int ret = msg.exec();
        if( ret == QMessageBox::Cancel) return;
    }

    initMesh();

    edges_weight_assigned = 0;
    if( !edges_weight_assigned ) setEdgesWeight();

    int inversion = checkInversionCheckBox->isChecked();
    laplace->checkInversion(inversion);

    int niter = numIterationsSpinBox->value();
    laplace->setNumIterations( niter);

    QString qstr;

    qstr = maxResidueLineEdit->text();
    double tol = qstr.toDouble();
    laplace->setTolerance(tol);

    bool val = applyTaubinStepsCheckBox->isChecked();
    laplace->setTaubinSteps(val);
    qstr = lambdaLineEdit->text() ;
    laplace->setLambda(qstr.toDouble() );

    if( val ) {
        qstr = muLineEdit->text() ;
        laplace->setMu(qstr.toDouble() );
    }

    if( explicitRadioButton->isChecked() ) {
        laplace->setNumericalMethod( JLaplaceMeshSmoother::EXPLICIT_METHOD );
    }

    if( implicitRadioButton->isChecked() ) {
        laplace->setNumericalMethod( JLaplaceMeshSmoother::IMPLICIT_METHOD );
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    laplace->smoothAll();
    meshViewer->updateBuffers(mesh);

    double maxdiff = laplace->getMaxResidual();
    currResidueLineEdit->setText( QString::number(maxdiff) );

    displayNegativeElements();
}

////////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: freeConstraints()
{
/*
    if( meshViewer == nullptr || mesh == nullptr ) return;
    mesh->deleteNodeAttribute("Constraint");
    meshViewer->displayAll(0,1);
*/
}
///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: getOriginal()
{
    if( mesh == nullptr ) return;
    mesh->getGeometry()->setCoordsArray(orgCoords, l2g);
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshLaplaceSmoothingDialog :: openConstraintsDialog()
{
    if( meshConstraintsDialog == nullptr )
        meshConstraintsDialog.reset(new JMeshConstraintsDialog(this));
    viewManager->attach(meshConstraintsDialog.get());

    meshConstraintsDialog->setViewManager( viewManager );
    meshConstraintsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: closeDialog()
{
    shrink_to_zero( orgCoords);
    shrink_to_zero( l2g );

    if( laplace ) laplace->clearAll();

    edges_weight_assigned = 0;

    this->hide();
    this->parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: setConstraints()
{
    if( mesh == nullptr) return;

    size_t nCount = mesh->getNumAttributes("Constraint", 0);
    numConstraintsLineEdit->setText( QString::number(nCount) );
    assignColors();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: makeConnections()
{
    ComboBox( edgeWeightComboBox, [=] {setEdgesWeight();});
    PushButton( applyPushButton,  [=] {smooth();});
    PushButton( orgmeshPushButton, [=] {getOriginal();});
    PushButton( constraintsPushButton,  [=] {openConstraintsDialog();});
    PushButton( closePushButton,  [=] {closeDialog();});
//  connect( preserveBoundaryCheckBox,  SIGNAL( toggled( bool ) ) ,  this, SLOT( setConstraints() ));
}
///////////////////////////////////////////////////////////////////////////////

void JMeshLaplaceSmoothingDialog :: displayNegativeElements()
{
#ifdef CSV
    if( mesh == nullptr) return;

    int topDim = mesh->getTopology()->getDimension();
    size_t nCount = 0;
    if( topDim == 2 ) {
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                if( JFaceGeometry::getSignedArea(face) < 0.0) nCount++;
            }
        }
    }
    numNegativeLineEdit->setText( QString::number(nCount) );

    if( displayNegativeCheckBox->isChecked() ) {
//      glDisable( GL_CULL_FACE);
        size_t numfaces = mesh->getSize(2);
        JFaceRenderPtr fAttrib;
        JColor color;
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        color[3] = 0.8;
        nCount = 0;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Render", fAttrib);
                fAttrib->display = 0;
                if( JFaceGeometry::isInverted(face) ) {
                    fAttrib->color = color;
                    fAttrib->display = 1;
                    nCount++;
                }
            }
        }
        //  glEnable( GL_CULL_FACE);
    }
    meshViewer->updateBuffers(mesh);
#endif

}
