#include "TetMesherDialog.hpp"
#include "MeshImporter.hpp"
#include "MeshExporter.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
JTetMesherDialog :: JTetMesherDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    exactinit();
    tetgenOptionsDialog.reset( new JTetGenOptionsDialog() );
}
///////////////////////////////////////////////////////////////////////////////

void JTetMesherDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////

JTetMesherDialog :: ~JTetMesherDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh(meshViewer->getCurrentMesh() );;
}

///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JTetMesherDialog :: genNewMesh()
{
    newTetMesh = nullptr;

    static int ncount = 0;
    if( meshViewer == nullptr || mesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    AllTetMeshGenerator tetmesher;

    if( convexHullRadioButton->isChecked() ) {
        newTetMesh = tetmesher.getConvexHull(mesh);
        newTetMesh->setName("ConvexHull");
        meshViewer->addObject(newTetMesh);
        mesh->setActiveBit(0);
    }

    if( cdtMeshRadioButton->isChecked() ) {
        newTetMesh = tetmesher.getConstrainedMesh(mesh);
        newTetMesh->setName("ConstrainedDT");
        meshViewer->addObject(newTetMesh);
        mesh->setActiveBit(0);
    }

    if( qualityMeshRadioButton->isChecked() ) {
        string cmd = tetgenOptionsDialog->getOptions();
        tetmesher.setOptions(cmd);
        newTetMesh = tetmesher.getQualityMesh(mesh);
        newTetMesh->setName("QualityDT");
        meshViewer->addObject(newTetMesh);
        mesh->setActiveBit(0);
    }
}
///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: fromHexMesh()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    int ncount = mesh->getTopology()->countElementType(JCell::HEXAHEDRON);
    if( ncount ) {
        AllTetMeshGenerator alltets;
        JMeshPtr tetmesh = alltets.fromHexMesh(mesh);
        if( tetmesh == nullptr);
        meshViewer->addObject(tetmesh);
        mesh->setActiveBit(0);
    }

}
///////////////////////////////////////////////////////////////////////////////

void JTetMesherDialog :: closeDialog()
{
    parentWidget()->show();
    close();
}
///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: openTetGenOptionsDialog()
{
    if( tetgenOptionsDialog == nullptr)
        tetgenOptionsDialog.reset( new JTetGenOptionsDialog(this));

    tetgenOptionsDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: getBoundedDistortion()
{
    JWaitCursor waitCursor;
    waitCursor.start();

    QString qstr = maxDistortionLineEdit->text();
    double K = qstr.toDouble();
    if( K > 1.0) {
        JBoundedDistortionOptimizer opt;
        if( newTetMesh)
            opt.optimize(newTetMesh, K);
        else
            opt.optimize(mesh, K);
        meshViewer->updateBuffers(newTetMesh);
    }
}
///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: getStellarOpt()
{
    JMeshPtr currMesh = mesh;
    if( newTetMesh ) currMesh = newTetMesh;

    JMeshTRIExporter mexp;
    mexp.writeFile(currMesh, "stellaropt");
    string cmd = "Stellar stellaropt";

    JMeshTRIImporter mimp;
    JMeshPtr optMesh = mimp.readFile("stellaropt.1.ele");
    optMesh->setName("OptMesh");
    currMesh->setActiveBit(0);
    meshViewer->addObject(optMesh);
}
///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: getMesquiteOpt()
{
    if( meshOptDialog == nullptr)
        meshOptDialog.reset( new JMeshOptDialog( this ));

    meshOptDialog->setViewManager( viewManager );

    JMeshPtr currMesh = mesh;
    if( newTetMesh ) currMesh = newTetMesh;

    meshViewer->setCurrentMesh(currMesh);
    tetgenOptionsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog :: openTetViewerDialog()
{
    if( meshCutterDialog == nullptr)
        meshCutterDialog.reset( new JImplicitMeshCutterDialog( this ));

    meshCutterDialog->setViewManager( viewManager );

    meshCutterDialog->setMesh(mesh);
    meshCutterDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JTetMesherDialog :: makeConnections()
{
    PushButton( generatePushButton,  [=] {genNewMesh();});
    PushButton( hex2tetsPushButton,  [=] {fromHexMesh();});
    PushButton( tetgenOptionsPushButton,  [=] {openTetGenOptionsDialog();});
    PushButton( boundedDistortionPushButton, [=] {getBoundedDistortion();});
    PushButton( stellarOptPushButton,  [=] {getStellarOpt();});
    PushButton( mesquiteOptPushButton,  [=] {getMesquiteOpt();});
    PushButton( tetviewerPushButton,  [=] {openTetViewerDialog(); });
    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
void JTetMesherDialog::MyThread :: run()
{
    /*
        int err = system (cmd.c_str() );

        if( err == -1 ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Tetmesh command failed ");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) return;
        }

        JMeshVTKImporter meshimp;
        JMeshPtr volmesh = meshimp.readFile("tmp.1.vtk");
        if( volmesh == nullptr ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Tetmesh  file could not be read ");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) return;
        }

        JMeshQuality mq(volmesh);
        vector<double> vol;
        mq.getCellsQuality(JMeshQuality::VOLUME, vol);
        double maxvol = *max_element(vol.begin(), vol.end() );
        dialog->maxvolLineEdit->setText( QString::number(maxvol) );
        dialog->opt_r_checkBox->setChecked(true);
        dialog->meshViewer->addObject( volmesh );
    */
}

