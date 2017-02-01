#include "SphereHexMesherDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JSphereHexMesherDialog :: JSphereHexMesherDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    selectMesh[0] = 0;
    selectMesh[1] = 0;
    sphHexMesher.reset(new JSphereHexMesher);
}

///////////////////////////////////////////////////////////////////////////////

JSphereHexMesherDialog :: ~JSphereHexMesherDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JSphereHexMesherDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    /*
        if( meshViewer == nullptr ) return;
        setMesh( meshViewer->getCurrentMesh() );
    */
}

///////////////////////////////////////////////////////////////////////////////

void JSphereHexMesherDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    topDim = mesh->getTopology()->getDimension();
    eulerChar = mesh->getTopology()->getEulerCharacteristic();
    if( eulerChar != 2 ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning : The input mesh is not topological sphere");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    sphHexMesher->setMesh(mesh);

    JMeshPtr cm;
    int err = mesh->getAttribute("CurvatureFlowMesh", cm);
    if( !err) {
        name = mesh->getName();
        objectNameLineEdit->setText(QString(name.c_str()));
        sphHexMesher->setSphere(mesh);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JSphereHexMesherDialog :: loadInputMesh()
{
    selectMesh[0] = 1;
    selectMesh[1] = 0;
    if( objectsListDialog == nullptr )
        objectsListDialog.reset(new JObjectsListDialog(this));

    objectsListDialog->setViewManager(viewManager );
    objectsListDialog->setType(0);
    objectsListDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JSphereHexMesherDialog :: loadSphereMesh()
{
    selectMesh[0] = 0;
    selectMesh[1] = 1;

    if( objectsListDialog == nullptr )
        objectsListDialog.reset(new JObjectsListDialog(this));

    objectsListDialog->setViewManager(viewManager);
    objectsListDialog->setType(0);
    objectsListDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JSphereHexMesherDialog :: showEvent( QShowEvent *e)
{
    if( e == nullptr) return;
    string name;
    if( selectMesh[0] ) {
        mesh = meshViewer->getCurrentMesh();
        name = mesh->getName();
        objectNameLineEdit->setText( QString(name.c_str()));
        sphHexMesher->setMesh(mesh);
    }

    if( selectMesh[1] ) {
        sphMesh = meshViewer->getCurrentMesh();
        name = sphMesh->getName();
        sphNameLineEdit->setText( QString(name.c_str()));
        sphHexMesher->setSphere(sphMesh);
    }

    selectMesh[0] = 0;
    selectMesh[1] = 0;
    selectMesh[2] = 0;
}

///////////////////////////////////////////////////////////////////////////////

void JSphereHexMesherDialog :: openCurvatureFlowDialog()
{
    if( meanCurvatureFlowDialog == nullptr)
        meanCurvatureFlowDialog.reset(new JMeshMeanCurvatureFlowDialog());

    meanCurvatureFlowDialog->setViewManager( viewManager );
    meanCurvatureFlowDialog->setMesh( mesh );
    meanCurvatureFlowDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JSphereHexMesherDialog :: openLaplaceDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JSphereHexMesherDialog :: openHarmonicDialog()
{
}
///////////////////////////////////////////////////////////////////////////////
void JSphereHexMesherDialog :: makeConnections()
{
    PushButton( curvatureFlowPushButton, [=] {openCurvatureFlowDialog();});
    PushButton( sphMeshPushButton, [=] {loadSphereMesh();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
