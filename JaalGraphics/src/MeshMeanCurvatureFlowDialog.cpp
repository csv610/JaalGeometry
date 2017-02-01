#include "MeshMeanCurvatureFlowDialog.hpp"
#include "MeshAffineTransforms.hpp"

JMeshMeanCurvatureFlowDialog :: JMeshMeanCurvatureFlowDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    timeStepLineEdit->setText( QString::number(0.000001) );
    progressBar->setTextVisible(0);
    rescale     = 0;
}

///////////////////////////////////////////////////////////////////////////////

JMeshMeanCurvatureFlowDialog :: ~JMeshMeanCurvatureFlowDialog()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshMeanCurvatureFlowDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
};

///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlowDialog :: initMesh()
{
    if( orgMesh == nullptr || initialized == 1) return;

    meanFlow.setMesh(orgMesh);
    orgMesh->setActiveBit(0);
    mesh = meanFlow.getSurfaceMesh();

    meshViewer->addObject(mesh);
    meshViewer->setCurrentMesh(mesh);

    if( cotanRadioButton->isChecked() )
        meanFlow.setAlgorithm( JMeshMeanCurvatureFlow::IGL_METHOD);

    if( keenanRadioButton->isChecked() )
        meanFlow.setAlgorithm( JMeshMeanCurvatureFlow::KEENAN_METHOD);

    /*
        bool val = intrinsicDelaunayCheckBox->isChecked();
        meanFlow.setIntrinsicDelaunayMesh(val);
    */

    currStep = 0;
    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    xlen0 = box.getLength(0);
    area0 = mesh->getGeometry()->getSurfaceArea();
    initAreaLineEdit->setText(QString::number(area0));

    int euler = mesh->getTopology()->getEulerCharacteristic();
    eulerCharLineEdit->setText(QString::number(euler));

    initialized = 1;
    viewManager->resetView(JViewDirection::FRONT_VIEW);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshMeanCurvatureFlowDialog :: setMesh( const JMeshPtr &m)
{
    orgMesh = m;
    initialized = 0;
    if( orgMesh  ) {
        string name = orgMesh->getName();
        objectNameLineEdit->setText(QString(name.c_str()));
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshMeanCurvatureFlowDialog :: displayDelaunay()
{
    if( meshViewer == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JDelaunayMesh2D delmesh;
    delmesh.setMesh(mesh);
    JEdgeColorPtr edgeColor(new JDelaunayEdgeColor);
//    JEdgeColor::assign(mesh, edgeColor);
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlowDialog :: startAllOver()
{
    currStep = 0;
    meanFlow.restart();
    mesh->getGeometry()->setFacesNormal();
    meshViewer->updateBuffers(mesh);
    initAreaLineEdit->setText(QString::number(area0));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlowDialog :: startFlow()
{
    if( !initialized ) initMesh();

    QString qstr = timeStepLineEdit->text();
    double  dt   = qstr.toDouble();
    meanFlow.setTimeStep(dt);

    JMeshAffineTransform affine;
    affine.setMesh(mesh);

    int numIter = numIterationsSpinBox->value();
    progressBar->setRange(0, numIter-1);

    JFaceRenderPtr  fAttrib;
    Vec3F normal;
    JColor color;
    for( int j = 0; j < numIter; j++) {
        progressBar->setValue(j);
        meanFlow.nextStep();
        /*
                if( preserveBoundingBoxCheckBox->isChecked() ) {
                    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
                    double xlen1  = box.getLength(0);
                    double sfactor = xlen0/xlen1;
                    affine.scale(sfactor, sfactor, sfactor);
                    affine.toCenter();
                }
        */
        area1 = mesh->getGeometry()->getSurfaceArea();
        if( fabs(area0-area1)/area0 < 0.01) {
            cout << "Warning: Surface area below threashold " << endl;
        }

        currAreaLineEdit->setText(QString::number(area1));
        mesh->getGeometry()->setFacesNormal();
        size_t numfaces = mesh->getSize(2);
//#pragma omp parallel for
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Render", fAttrib);
                face->getAttribute("Normal", normal);
                color[0] = 0.5*(normal[0]+1.0);
                color[1] = 0.5*(normal[1]+1.0);
                color[2] = 0.5*(normal[2]+1.0);
                fAttrib->color = color;
            }
        }
        meshViewer->updateBuffers(mesh);
        glFlush();
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshMeanCurvatureFlowDialog :: projectOnSphere()
{
    if( !initialized ) initMesh();

    meanFlow.projectOnSphere();
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshMeanCurvatureFlowDialog :: makeConnections()
{
//  PushButton( displayNonDelaunayPushButton, [=] {displayDelaunay();});
    PushButton( sphProjectionPushButton, [=] {projectOnSphere();});
    PushButton( resetPushButton, [=] {startAllOver();});
    PushButton( applyPushButton, [=] {startFlow();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
