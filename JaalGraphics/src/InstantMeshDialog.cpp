#include "InstantMeshDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JInstantMeshDialog :: JInstantMeshDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JInstantMeshDialog :: ~JInstantMeshDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JInstantMeshDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JInstantMeshDialog :: showEvent( QShowEvent *event)
{
    if( meshViewer ) setMesh( meshViewer->getCurrentMesh() );
    QDialog::showEvent(event);

    if( meshtype == 3) {
        rot6RadioButton->setChecked(true);
        pos6RadioButton->setChecked(true);
    } else {
        rot4RadioButton->setChecked(true);
        pos4RadioButton->setChecked(true);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JInstantMeshDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    int numNodes = mesh->getSize(0);
    inNodesLineEdit->setText( QString::number(numNodes));
    expectedNodesLineEdit->setText( QString::number(numNodes));

    int numFaces = mesh->getSize(2);
    inFacesLineEdit->setText( QString::number(numFaces));
    expectedFacesLineEdit->setText( QString::number(numFaces));

    JMeshIO mio;
    mio.saveAs( mesh, "tmp1.obj");

    inEulerCharacteristic = mesh->getTopology()->getEulerCharacteristic();
}

///////////////////////////////////////////////////////////////////////////////
void JInstantMeshDialog :: genMesh()
{
    if( newMesh  ) meshViewer->removeObject( newMesh );

    QString qstr;
    ostringstream oss;
    oss << "imesh -o tmp2.obj ";

    if( rot2RadioButton->isChecked() )
        oss << " -r 2 " ;

    if( rot4RadioButton->isChecked() )
        oss << " -r 4 " ;

    if( rot6RadioButton->isChecked() )
        oss << " -r 6 " ;

    if( pos4RadioButton->isChecked() )
        oss << " -p 4 " ;

    if( pos6RadioButton->isChecked() )
        oss << " -p 6 " ;

    if( expectedNodesRadioButton->isChecked() )  {
        qstr = expectedNodesLineEdit->text();
        if( triquadDominantCheckBox->isChecked() || meshtype == 3)
            oss << " -v " << qstr.toInt();
        else
            oss << " -v " << ceil(qstr.toInt()/4.0);
    }

    if( expectedFacesRadioButton->isChecked() )  {
        qstr = expectedFacesLineEdit->text();
        if( triquadDominantCheckBox->isChecked() || meshtype == 3)
            oss << " -f " << qstr.toInt();
        else
            oss << " -f " << ceil(qstr.toInt()/4.0);
    }

    oss << " -S " << numSmoothSpinBox->value();
    oss << " -t " << numThreadsSpinBox->value();

    if( intrinsicCheckBox->isChecked() )
        oss << " -i ";

    if( triquadDominantCheckBox->isChecked() )
        oss << " -D ";

    if( alignBoundaryCheckBox->isChecked() )
        oss << " -b ";

    if( creaseCheckBox->isChecked() )
        oss << " -c " << creaseAngleSpinBox->value();

    oss << " tmp1.obj";

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render",mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;

    JWaitCursor waitCursor;
    waitCursor.start();
    string cmd = oss.str();
    system( cmd.c_str() );

    JMeshIO mio;
    newMesh = mio.readFile("tmp2.obj");
    meshViewer->addObject( newMesh );
    newMesh->setName("IQmesh");
    objectNameLineEdit->setText(QString("IQmesh"));

    expectedNodesLineEdit->setText( QString::number(newMesh->getSize(0)));
    expectedFacesLineEdit->setText( QString::number(newMesh->getSize(2)));


    waitCursor.stop();
    if( newMesh->getTopology()->getEulerCharacteristic() != inEulerCharacteristic)
    {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Euler Characteristic changed: Refine more ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JInstantMeshDialog :: rejectMesh()
{
    if( newMesh ) meshViewer->removeObject( newMesh );
    newMesh.reset();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render",mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 1;
    mrender->displayEntity[2] = 2;
    mrender->displayEntity[3] = 0;

    meshViewer->refreshDisplay();

    expectedNodesLineEdit->setText( QString::number(mesh->getSize(0)));
    expectedFacesLineEdit->setText( QString::number(mesh->getSize(2)));

    string name =mesh->getName();
    objectNameLineEdit->setText(QString( name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////

void JInstantMeshDialog :: closeDialog()
{
    if( newMesh )  {
        meshViewer->setCurrentMesh(newMesh);

        if( deleteInMeshCheckBox->isChecked() ) {
            meshViewer->removeObject(mesh);
            mesh->deleteAll();
        } else {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render",mrender);
        mrender->displayEntity[0] = 0;
        mrender->displayEntity[1] = 0;
        mrender->displayEntity[2] = 0;
        mrender->displayEntity[3] = 0;
        }
    }
   
    this->close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////
void JInstantMeshDialog :: openRenderMeshDialog()
{
    if( newMesh == nullptr) return;

    if( meshRenderDialog == nullptr )
        meshRenderDialog.reset(new JMeshRenderDialog(this));

    meshRenderDialog->setViewManager( viewManager );
    meshRenderDialog->setMesh(newMesh);
    meshRenderDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JInstantMeshDialog :: showOrgBoundary()
{
}

void JInstantMeshDialog :: makeConnections()
{
    PushButton( applyPushButton,   [=] { genMesh(); });
    PushButton( rejectPushButton,  [=] { rejectMesh(); });
    PushButton( renderMeshPushButton,   [=] { openRenderMeshDialog(); });
    CheckBox( showOrgBoundaryCheckBox,  [=] { showOrgBoundary(); });
    PushButton( closePushButton,   [=] { closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
