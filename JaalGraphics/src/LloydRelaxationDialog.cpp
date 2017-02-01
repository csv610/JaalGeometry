#include "LloydRelaxationDialog.hpp"

JLloydRelaxationDialog :: JLloydRelaxationDialog( QWidget *parent) : QDialog(parent)
{
    viewManager = nullptr;
    meshViewer  = nullptr;

    setupUi(this);
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

JLloydRelaxationDialog :: ~JLloydRelaxationDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JLloydRelaxationDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JLloydRelaxationDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    initialized = 0;
    if( mesh ) {
        string name = mesh->getName();
        objectNameLineEdit->setText(QString(name.c_str()));
    }

    dualGraph.reset();
}
///////////////////////////////////////////////////////////////////////////////

void JLloydRelaxationDialog :: initMesh()
{
    if( mesh == nullptr || initialized == 1) return;

    mesh->getGeometry()->getCoordsArray(orgCoords,l2g);
    mopt.setMesh(mesh);

    JMeshRenderPtr  mrender;
    if( dualGraph == nullptr) {
        string name = mesh->getName();
        dualGraph = mopt.getDualGraph();
        dualGraph->setName("dual_" + name);
        meshViewer->addObject(dualGraph);
        dualGraph->getAttribute("Render", mrender);
        mrender->displayEntity[1] = displayDualCheckBox->isChecked();
    }

    initialized = 1;
}

///////////////////////////////////////////////////////////////////////////////

void JLloydRelaxationDialog :: smoothAll()
{
    JWaitCursor waitCursor;
    waitCursor.start();

    initMesh();

    int numiters = numIterSpinBox->value();
    mopt.setNumIterations(numiters);

    QString qstr;
    qstr = toleranceLineEdit->text();
    double tol = qstr.toDouble();
    mopt.setTolerance(tol);

    bool p = preserveBoundaryCheckBox->isChecked();
    mopt.setPreserveBoundary(p);

    mopt.smoothAll();

    numiters = mopt.getNumIterations();
    numIterSpinBox->setValue( numiters);
    double val = mopt.getMaxResidue();
    maxResidueLineEdit->setText( QString::number(val) );

    meshViewer->updateBuffers(mesh);
    meshViewer->updateBuffers(dualGraph);
}
///////////////////////////////////////////////////////////////////////////////

void JLloydRelaxationDialog :: checkDisplay()
{
    if( meshViewer == nullptr) return;
    initMesh();

    JMeshRenderPtr  mrender;
    if( mesh ) {
        mesh->getAttribute("Render", mrender);
        mrender->displayEntity[1] = displayPrimalCheckBox->isChecked();
    }

    if( dualGraph ) {
        dualGraph->getAttribute("Render", mrender);
        mrender->displayEntity[1] = displayDualCheckBox->isChecked();
    }

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLloydRelaxationDialog :: getOriginal()
{
    if( mesh == nullptr) return;

    mesh->getGeometry()->setCoordsArray(orgCoords, l2g);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JLloydRelaxationDialog :: closeDialog()
{
    if( dualGraph) {
        dualGraph->deleteAll();
        if( meshViewer ) meshViewer->removeObject(dualGraph);
        dualGraph.reset();
        meshViewer->refreshDisplay();
    }

    shrink_to_zero( orgCoords);
    shrink_to_zero( l2g );
    parentWidget()->show();
    close();
}
///////////////////////////////////////////////////////////////////////////////

void JLloydRelaxationDialog :: makeConnections()
{
    connect( applyPushButton, SIGNAL( clicked() ), this, SLOT( smoothAll() ));

    connect( displayPrimalCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( checkDisplay() ));
    connect( displayDualCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( checkDisplay() ));
    connect( orgmeshPushButton, SIGNAL( clicked() ), this, SLOT( getOriginal() ));

    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
