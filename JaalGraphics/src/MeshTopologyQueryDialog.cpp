#include "MeshTopologyQueryDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JMeshTopologyQueryDialog :: JMeshTopologyQueryDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    fontScaleLineEdit->setText( QString::number(1) );
}

///////////////////////////////////////////////////////////////////////////////
JMeshTopologyQueryDialog :: ~JMeshTopologyQueryDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: getInfo()
{
    JWaitCursor waitCursor;
    waitCursor.start();

    JMeshTopologyPtr topology = mesh->getTopology();
    int topdim = topology->getDimension();
    topdimLineEdit->setText( QString::number(topdim) );

    int elemType = topology->getElementsType(topdim);

    if( elemType )
        homogeneousCheckBox->setChecked(1);
    else
        homogeneousCheckBox->setChecked(0);

    if( topdim == 2 ) {
        switch( elemType ) {
        case JFace::TRIANGLE:
            elemTypeLineEdit->setText( QString("Triangle"));
            break;
        case JFace::QUADRILATERAL:
            elemTypeLineEdit->setText( QString("Quadrilateral"));
            break;
        }
    }

    if( topdim == 3 ) {
        switch( elemType ) {
        case JCell::TETRAHEDRON:
            elemTypeLineEdit->setText( QString("Tetrahedra"));
            break;
        case JCell::HEXAHEDRON:
            elemTypeLineEdit->setText( QString("Hexahedra"));
            break;
        }
    }

    bool val;
    val = mesh->getTopology()->isConsistent();
    consistentCheckBox->setChecked(val);

    val = mesh->getTopology()->isClosed();
    closedCheckBox->setChecked(val);

    val = mesh->getTopology()->isSimple();
    manifoldCheckBox->setChecked(val);


    int eval = mesh->getTopology()->getEulerCharacteristic();
    eulerCharLineEdit->setText( QString::number(eval) );

    JNodeSequence nodedoublets;
    mesh->getTopology()->getDoublets(nodedoublets);

    nodeDoubletsLineEdit->setText( QString::number(nodedoublets.size()) );

    JEdgeSequence edgedoublets;
    mesh->getTopology()->getDoublets(edgedoublets);
    edgeDoubletsLineEdit->setText( QString::number(edgedoublets.size()) );
}


///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: setMesh( const JMeshPtr &m)
{
    initialized = 0;
    mesh = m;
    if( mesh == nullptr ) return;
    objectNameLineEdit->setText( QString(mesh->getName().c_str() ));
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: showEvent( QShowEvent *e)
{
    if( e == nullptr || mesh == nullptr) return;

    if( !initialized ) {
        getInfo();
        initialized = 1;
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: getConsistent()
{
    if( mesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    mesh->getTopology()->getConsistent();
    bool val;
    val = mesh->getTopology()->isConsistent();
    consistentCheckBox->setChecked(val);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: displayMatrix()
{
    if( meshViewer == nullptr ) return;

    if( displayMatrixCheckBox->isChecked() ) {
        if(matrixViewer == nullptr) {
            JViewComponentPtr c = JMatrixViewer::registerComponent(viewManager);
            matrixViewer = dynamic_pointer_cast<JMatrixViewer>(c);
//          matrixViewer->setViewManager(viewManager);
        }
        viewManager->attach(matrixViewer );
        int pointSize = pointSizeSpinBox->value();
        matrixViewer->setPointSize( pointSize );
        viewManager->deactivateComponents();
        matrixViewer->setActive(1);
    } else {
        if( matrixViewer ) {
            matrixViewer->setActive(0);
            viewManager->detach(matrixViewer );
            meshViewer->getViewManager()->activateComponents();
         }
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: displayPrimalMatrix()
{
    if( matrixViewer == nullptr) return;
    matrixViewer->setMesh( mesh );
    matrixViewer->setGraphType(JMeshTopology::PRIMAL_GRAPH);
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: displayDualMatrix()
{
    if( matrixViewer == nullptr) return;
    matrixViewer->setMesh( mesh );
    matrixViewer->setGraphType( JMeshTopology::DUAL_GRAPH);
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: setLapFonts()
{
    if( matrixViewer == nullptr ) return;

    QString str = fontScaleLineEdit->text();
    matrixViewer->setFontsScale( str.toDouble() );
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: setLapGrid()
{
    if( matrixViewer == nullptr ) return;
    bool val = gridCheckBox->isChecked();
    matrixViewer->setMatrixGrid(val);
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: setStepLabel()
{
    if( matrixViewer == nullptr ) return;

    QString str = stepNodeLabelLineEdit->text();
    matrixViewer->setStepLabels( str.toInt() );
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: setPointSize()
{
    if( matrixViewer == nullptr ) return;

    displayMatrix();

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: getBettiNumber()
{
    int numloops = mesh->getTopology()->getBettiNumber(BettiMethod::SMITH_NORMAL_FORM);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: openMeshComponentsDialog()
{
    if( meshComponentsDialog == nullptr)
        meshComponentsDialog.reset( new JMeshComponentsDialog(this) );

    meshComponentsDialog->setViewManager( viewManager );
    meshComponentsDialog->setMesh( mesh );
    meshComponentsDialog->show();
    this->hide();

}
///////////////////////////////////////////////////////////////////////////////
void JMeshTopologyQueryDialog :: getOrphaned()
{
    if( mesh == nullptr) return;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: closeDialog()
{
    if( meshViewer == nullptr ) return;
    meshViewer->getViewManager()->detach(matrixViewer );
    this->close();
    meshViewer->getViewManager()->activateComponents();
    parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshTopologyQueryDialog :: makeConnections()
{
    PushButton( getOrphanedPushButton,   [=] {getOrphaned();});
    PushButton( getConsistentPushButton, [=] {getConsistent();});
    PushButton( primalGraphPushButton,   [=] {displayPrimalMatrix();});
    PushButton( dualGraphPushButton,     [=] {displayDualMatrix();});
    PushButton( getBettiNumberPushButton, [=] {getBettiNumber();});
    PushButton( editComponentsPushButton, [=] {openMeshComponentsDialog();});

    CheckBox( displayMatrixCheckBox, [=] {displayMatrix();});
    CheckBox( gridCheckBox, [=] {setLapGrid();});

    SpinBoxi( pointSizeSpinBox,  [=] { displayMatrix();});
    LineEdit( fontScaleLineEdit,  [=] { setLapFonts();});
    LineEdit( stepNodeLabelLineEdit, [=] { setStepLabel(); });

    PushButton( closePushButton, [=] {closeDialog(); });
}
///////////////////////////////////////////////////////////////////////////////

