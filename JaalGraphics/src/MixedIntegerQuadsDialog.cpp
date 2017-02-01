#include "MixedIntegerQuadsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMixedIntegerQuadsDialog :: JMixedIntegerQuadsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    mixedIntegerQuads.reset(new JMixedIntegerQuads);
}

///////////////////////////////////////////////////////////////////////////////

JMixedIntegerQuadsDialog :: ~JMixedIntegerQuadsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMixedIntegerQuadsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMixedIntegerQuadsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JMixedIntegerQuadsDialog :: checkDisplay()
{
    bool val;
    val = displayCrossFieldCheckBox->isChecked();
}
///////////////////////////////////////////////////////////////////////////////
void JMixedIntegerQuadsDialog :: generateQuads()
{
    if( mesh == nullptr) return;
    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));
    mixedIntegerQuads->setMesh(mesh);

    uvMesh  = mixedIntegerQuads->getUVMesh();  // First to be called ...
    meshViewer->addObject(uvMesh);
    mesh->setActiveBit(0);

    /*
      mesh->setAttribute("UVMesh", uvMesh);
      crossField          = mixedIntegerQuads->getCrossField();
      bisectorField       = mixedIntegerQuads->getBisectorField();
      bisectCombinedField = mixedIntegerQuads->getCombinedBisectorField();
    */
    checkDisplay();

    QApplication::restoreOverrideCursor();
}
///////////////////////////////////////////////////////////////////////////////
void JMixedIntegerQuadsDialog :: makeConnections()
{
    CheckBox( displayCrossFieldCheckBox,    [=] {checkDisplay();});
    CheckBox( displayBisectorFieldCheckBox, [=] {checkDisplay();});
    CheckBox( displayGlobalParameterCheckBox, [=] {checkDisplay();});
    CheckBox( displaySingularitiesCheckBox,  [=] {checkDisplay();});
    CheckBox( displayBisectorFieldCombinedCheckBox,  [=] {checkDisplay();});

    PushButton( generatePushButton, [=] {generateQuads();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
