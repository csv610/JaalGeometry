#include "MorseAnalysisDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMorseAnalysisDialog :: JMorseAnalysisDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();
    mesh = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMorseAnalysisDialog :: ~JMorseAnalysisDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMorseAnalysisDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    mesh = meshViewer->getCurrentMesh();
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );

}
///////////////////////////////////////////////////////////////////////////////
void JMorseAnalysisDialog :: genHeight()
{
    if( meshViewer == nullptr ) return;

    string name = StdString(objectNameLineEdit->text());
    mesh = meshViewer->getMesh(name);
    if( mesh == nullptr) return;

    QString qstr = heightComboBox->currentText();
    int dir = 0;
    if( qstr == "X-Coord") dir = 0;
    if( qstr == "Y-Coord") dir = 1;
    if( qstr == "Z-Coord") dir = 2;

    JMorseAnalysis morse;
    morse.setMesh(mesh);
    morse.setHeightDirection(dir);

    JMeshRender::setNodeScalarFieldColor(mesh, "MorseScalar");

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    mrender->setSurfaceShade(JRender::SMOOTH_SHADE);
    meshViewer->updateBuffers(mesh);
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMorseAnalysisDialog :: genContours()
{
    if( meshViewer == nullptr) return;

    if( mesh == nullptr ) return ;

    JWaitCursor wCursor;
    wCursor.start();

    JMarchingTriangles march;
    march.setMesh(mesh, "MorseScalar");

    march.getContours(10, contours);

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMorseAnalysisDialog :: makeConnections()
{
    PushButton( applyPushButton,   [=] {genHeight();});
    PushButton( contourPushButton, [=] {genContours();});
    PushButton( closePushButton,   [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
