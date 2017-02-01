#include "MeshUtilsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshUtilsDialog :: JMeshUtilsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
//    paramLinesColor.reset(new QuadParamLinesColor);
}

///////////////////////////////////////////////////////////////////////////////

JMeshUtilsDialog :: ~JMeshUtilsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshUtilsDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;
//   mesh = meshViewer->getMesh();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshUtilsDialog :: getQuadParamLines()
{
    if( mesh == nullptr ) return;
    QuadParametricLines qlines;
    qlines.setLines(mesh);

//  meshViewer->getEdgeDraw()->setColorMethod( paramLinesColor );
//  meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshUtilsDialog :: closeDialog()
{
//   meshViewer->getDrawEdge()->setColorMethod( nullptr );
    meshViewer->refreshDisplay();
    this->close();
}

void JMeshUtilsDialog :: makeConnections()
{
    PushButton( getParamlinesPushButton,  [=] {getQuadParamLines();});
    PushButton( closePushButton, [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
