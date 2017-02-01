#include "HeatConductionDialog.hpp"

JHeatConductionDialog :: JHeatConductionDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
    picker = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JHeatConductionDialog :: ~JHeatConductionDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JHeatConductionDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    picker = meshViewer->getEntityPicker();
    if( picker ) {
        picker->setMode(1);
    }

/*
    if( mesh == nullptr ) return ;
    mesh->buildRelations(0,0);
    heatSolver.setMesh(mesh);
    mesh->getTopology()->searchBoundary();
*/
}

///////////////////////////////////////////////////////////////////////////////

void JHeatConductionDialog ::solve()
{
    heatSolver.solve();
}

///////////////////////////////////////////////////////////////////////////////

void JHeatConductionDialog :: applyBoundConditions()
{
    if( picker == nullptr ) return;

    JNodeSequence nodeSeq = picker->getPickedNodes();
    if( nodeSeq.empty() ) return;

    QString str = dirichletBoundValueLineEdit->text() ;
    double val  = str.toDouble();

    for( size_t i = 0; i < nodeSeq.size(); i++)  {
        heatSolver.setBoundaryCondition( nodeSeq[i], 0, val);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JHeatConductionDialog :: newGroup()
{
    if( picker == nullptr ) return;
    picker->clearAll();
}

///////////////////////////////////////////////////////////////////////////////

void JHeatConductionDialog :: setLinearSolver()
{
    int niter = numIterationsSpinBox->value() ;
    heatSolver.setMaxIterations(niter);
}
///////////////////////////////////////////////////////////////////////////////

void JHeatConductionDialog :: makeConnections()
{
    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( close() ));
    connect( newGroupPushButton,  SIGNAL( clicked() ), this, SLOT( newGroup() ));
    connect( applyBoundConditionsPushButton,  SIGNAL( clicked() ), this, SLOT( applyBoundConditions() ));
    connect( solvePushButton,  SIGNAL( clicked() ), this, SLOT( solve() ));
}

///////////////////////////////////////////////////////////////////////////////
