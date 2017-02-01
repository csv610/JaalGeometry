#include "ContourEditingDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JContourEditingDialog :: JContourEditingDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JContourEditingDialog :: ~JContourEditingDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JContourEditingDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JContourEditingDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
void JContourEditingDialog :: reparameterize()
{
   extractContours();
   cout << "Reparam " << endl;
}
///////////////////////////////////////////////////////////////////////////////
void JContourEditingDialog :: closeDialog()
{
    this->close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////
void JContourEditingDialog :: smoothing()
{
   extractContours();
   cout << "Smooth" << endl;
}

///////////////////////////////////////////////////////////////////////////////
void JContourEditingDialog :: curveShortening()
{
   extractContours();
   cout << "Curve " << endl;
}

///////////////////////////////////////////////////////////////////////////////
void JContourEditingDialog :: openPolySimplifyDialog()
{
}

void JContourEditingDialog :: extractContours()
{
    if( nullptr == mesh) return;
    if( contours.empty() ) mesh->getTopology()->getBoundary( contours );
}
///////////////////////////////////////////////////////////////////////////////

void JContourEditingDialog :: makeConnections()
{
    PushButton( reparamPushButton,  [=] {reparameterize(); });
    PushButton( polygonSimplifyPushButton,  [=] {openPolySimplifyDialog(); });
    PushButton( smoothPushButton,   [=] {smoothing(); });
    PushButton( curveShorteningPushButton,   [=] {curveShortening(); });
    PushButton( closePushButton,    [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
