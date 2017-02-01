#include "MeshBooleanDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshBooleanDialog :: JMeshBooleanDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    meshA = nullptr;
    meshB = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    selectMesh[0] = 0;
    selectMesh[1] = 0;
}

///////////////////////////////////////////////////////////////////////////////

JMeshBooleanDialog :: ~JMeshBooleanDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshBooleanDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshBooleanDialog :: setMesh( const JMeshPtr &mA, const JMeshPtr &mB)
{
    meshA = mA;
    if( meshA == nullptr ) return ;

    string name = meshA->getName();
    meshALineEdit->setText(QString(name.c_str()));

    meshB = mB;
    if( meshB == nullptr ) return ;

    name = meshB->getName();
    meshBLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JMeshBooleanDialog :: loadMeshA()
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
void JMeshBooleanDialog :: loadMeshB()
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
void JMeshBooleanDialog :: applyOp()
{
    QString qstr;
    qstr = boolOpComboBox->currentText();

    if( meshA == nullptr || meshB == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    if( meshBoolean == nullptr)
        meshBoolean.reset( new JMeshBoolean);

    string opname = StdString(qstr);
    int    opid  =  JMeshBoolean::getOp( opname);

    meshBoolean->setMesh(meshA, meshB);
    meshC = meshBoolean->getMesh(opid);
    if( meshC ) {
        meshC->setName("MeshBoolean");
        meshViewer->addObject(meshC);
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshBooleanDialog :: showEvent( QShowEvent *e)
{
    if( e == nullptr) return;
    string name;
    if( selectMesh[0] ) {
        meshA = meshViewer->getCurrentMesh();
        name = meshA->getName();
        meshALineEdit->setText( QString(name.c_str()));
    }

    if( selectMesh[1] ) {
        meshB = meshViewer->getCurrentMesh();
        name = meshB->getName();
        meshBLineEdit->setText( QString(name.c_str()));
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshBooleanDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshBooleanDialog :: makeConnections()
{
    PushButton( meshAPushButton, [=] {loadMeshA(); });
    PushButton( meshBPushButton, [=] {loadMeshB(); });
    PushButton( applyPushButton, [=] {applyOp(); });
    PushButton( closePushButton, [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
