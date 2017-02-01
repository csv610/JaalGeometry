#include "MeshDualGrapherDialog.hpp"

JMeshDualGrapherDialog :: JMeshDualGrapherDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshDualGrapherDialog :: ~JMeshDualGrapherDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualGrapherDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshDualGrapherDialog :: newGraph()
{
    if( meshViewer == nullptr ) return;

    JMeshPtr mesh = meshViewer->getCurrentMesh();
    if( mesh == nullptr ) return;

    JMeshDualGraph dgrapher;
    dgrapher.setMesh(mesh);

    int pos;

    if( centerRadioButton->isChecked() )
        pos = JMeshDualGraph::DUAL_AT_CENTER;

    if( circumcenterRadioButton->isChecked() )
        pos = JMeshDualGraph::DUAL_AT_CIRCUMCENTER;

    if( center_circumcenterRadioButton->isChecked() )
        pos = JMeshDualGraph::DUAL_AT_CIRCUMCENTER_OR_CENTER;

    if( orthocenterRadioButton->isChecked() )
        pos = JMeshDualGraph::DUAL_AT_ORTHOCENTER;

    if( orthocenterRadioButton->isChecked() )
        pos = JMeshDualGraph::DUAL_AT_ORTHOCENTER;

    /*
        if( hodgeCenterRadioButton->isChecked() )
            pos = JMeshDualGraph::DUAL_AT_HODGE_CENTER;
    */

    dgrapher.setDualPosition(pos);

    bool hanger = 0;
    if( boundaryHangersCheckBox->isChecked() )
        hanger = 1;
    dgrapher.setBoundaryNodes(hanger);

    /*
            bool midnodes = 0;
            if( midnodesCheckBox->isChecked() )
                midnodes = 1;
            dualgrapher.setMidNodes(midnodes);
    */
    dGraph = dgrapher.getGraph();

    if( dGraph) meshViewer->addObject(dGraph);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualGrapherDialog :: openMeshRenderDialog()
{
    if( meshRenderDialog == nullptr )
        meshRenderDialog.reset(new JMeshRenderDialog( this ));

    meshRenderDialog->setViewManager( viewManager );
    meshRenderDialog->setMesh(dGraph);
    meshRenderDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshDualGrapherDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}

void JMeshDualGrapherDialog :: makeConnections()
{
    PushButton( renderPushButton, [=] { openMeshRenderDialog(); });
    PushButton( generatePushButton, [=] { newGraph(); });
    PushButton( closePushButton,    [=] { closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
