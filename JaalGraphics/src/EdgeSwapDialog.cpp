#include "EdgeSwapDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JEdgeSwapDialog :: JEdgeSwapDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    currPos  = 0;
    curredge = nullptr;
    meshViewer  = nullptr;
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////
JEdgeSwapDialog :: ~JEdgeSwapDialog()
{ }

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

    if( mesh == nullptr ) return ;

    mesh->buildRelations(0,2);

    qedgeswapper.reset(new JSwapQuadEdge());
    qedgeswapper->setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: startswap()
{
    if( meshViewer == nullptr ) return;

    if( mesh == nullptr ) return ;

    mesh->buildRelations(0,2);

    JQuadCleanUp qClean;
    qClean.setMesh(mesh);
    qClean.vertex_degree_reduction();

    meshViewer->getViewManager()->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: selectiveEdge()
{
    /*
         int  val  = selectiveEdgeCheckBox->isChecked();
         if( meshViewer == nullptr ) return;
         entityPicker->setPickableEntity(1) ;
         entityPicker->setMode(val); // Select only one edge..
         if( mesh == nullptr ) return ;
         mesh->buildRelations(1,2);
    */
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: swapSelectiveEdge()
{
    if( meshViewer == nullptr ) return;

    JEdgeSequence edges = entityPicker->getPickedEdges();
    if( edges.empty() ) return;

    if( mesh == nullptr ) return ;

    JEdgePtr swapedge = edges[0];

    int err = qedgeswapper->applyAt(swapedge);
    if( !err ) meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: interactiveEdge()
{
    if( qedgeswapper == nullptr ) return;

    size_t nSize = mesh->getSize(1);
    curredge = nullptr;

    for( size_t i = 0; i < nSize; i++) {
        JEdgePtr e = mesh->getEdgeAt(i);
        if( qedgeswapper->isSwappable(e) ) {
            curredge = e;
            currPos  = i+1;
            return;
        }
    }

//     meshViewer->showSpecificEntity( curredge );
    if( curredge == nullptr ) {
        cout << "Info: No swappable edge found in the mesh " << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: swapCurrEdge()
{
    if( qedgeswapper == nullptr ) return;

    qedgeswapper->applyAt( curredge );

    curredge = nullptr;

    size_t nSize = mesh->getSize(1);
    for( size_t i = currPos; i < nSize; i++) {
        JEdgePtr e = mesh->getEdgeAt(i);
        if( qedgeswapper->isSwappable(e) ) {
            curredge = e;
            currPos  = i+1;
            return;
        }
    }

//     meshViewer->showSpecificEntity( curredge );

    if( curredge == nullptr ) {
        cout << "Info: No swappable edge found in the mesh " << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: skipCurrEdge()
{
    /*
        if( qedgeswapper == nullptr ) return;

        curredge = nullptr;

        size_t nSize = mesh->getSize(1);
        for( size_t i = currPos; i < nSize; i++) {
            JEdgePtr e = mesh->getEdgeAt(i);
            if( qedgeswapper->isSwappable(e) ) {
                curredge = e;
                currPos  = i+1;
                return;
            }
        }

    //   meshViewer->showSpecificEntity( curredge );

        if( curredge == nullptr ) {
            cout << "Info: No swappable edge found in the mesh " << endl;
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JEdgeSwapDialog :: accept()
{
    if( meshViewer == nullptr ) return;

    /*
         entityPicker->setPickableEntity(); // Select only one edge..
    //   meshViewer->getDrawNode()->setColorMethod( nullptr );
         this->hide();
         meshViewer->getViewManager()->refreshDisplay();
         selectiveEdgeCheckBox->setChecked(0);
         interactiveCheckBox->setChecked(0);
    */
}
///////////////////////////////////////////////////////////////////////////////


void JEdgeSwapDialog :: makeConnections()
{
    connect( startSwapPushButton,  SIGNAL( clicked() ), this, SLOT( startswap() ));
    connect( rotateEdgePushButton, SIGNAL( clicked() ), this, SLOT( swapSelectiveEdge() ));
    connect( selectiveEdgeCheckBox,  SIGNAL( toggled( bool ) ) , this, SLOT( selectiveEdge() ));
    connect( interactiveCheckBox,  SIGNAL( toggled( bool ) ) , this, SLOT( interactiveEdge() ));
}

///////////////////////////////////////////////////////////////////////////////
