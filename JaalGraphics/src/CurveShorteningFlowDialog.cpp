#include "CurveShorteningFlowDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JCurveShorteningFlowDialog :: JCurveShorteningFlowDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    csflow.reset( new JCurveShorteningFlow);
}

///////////////////////////////////////////////////////////////////////////////

JCurveShorteningFlowDialog :: ~JCurveShorteningFlowDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JCurveShorteningFlowDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
}

///////////////////////////////////////////////////////////////////////////////

void JCurveShorteningFlowDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    csflow->setMesh(mesh);

    mesh->getGeometry()->getCoordsArray(orgCoords,l2g);

    int dim;
    dim = mesh->getTopology()->getDimension();
    if( dim == 1) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[2] = 0;

    JEdgeRenderPtr eAttrib;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", eAttrib);
            if( edge->isBoundary() )
                eAttrib->display = 1;
            else
                eAttrib->display = 0;
        }
    }
    currStep = 0;
}

///////////////////////////////////////////////////////////////////////////////

void JCurveShorteningFlowDialog :: startflow()
{
    bool animate = animationCheckBox->isChecked();
    ostringstream oss;
    viewManager->setSnapshotQuality(100);
    viewManager->setSnapshotFormat("PNG");

    JWaitCursor  wcursor;
    wcursor.start();
    int niter =  numIterationsSpinBox->value();
    for( int i = 0; i < niter; ++i) {
        qApp->processEvents();
        csflow->performOneStep();
        meshViewer->updateBuffers(mesh);
        double len = csflow->getCurveLength();
        curveLengthLineEdit->setText(QString::number(len) );
        if( animate ) {
            oss << "./Animation/csf" << std::setw(6) << std::setfill('0') << currStep << ".png";
            QString qstr = QString::fromStdString(oss.str());
            viewManager->saveSnapshot( qstr );
            currStep++;
            oss.str("");
            oss.clear();
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void JCurveShorteningFlowDialog :: getOriginal()
{
    if( mesh == nullptr) return;

    mesh->getGeometry()->setCoordsArray(orgCoords, l2g);

    meshViewer->updateBuffers(mesh);
    currStep = 0;
}

///////////////////////////////////////////////////////////////////////////////

void JCurveShorteningFlowDialog :: closeDialog()
{
    this->close();
}

///////////////////////////////////////////////////////////////////////////////

void JCurveShorteningFlowDialog :: makeConnections()
{
    PushButton( applyPushButton,  [=] {startflow();});
    PushButton( undoPushButton,   [=] {getOriginal(); });
    PushButton( closePushButton,  [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
