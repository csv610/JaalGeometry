#include "MeshSingularityGraphDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshSingularityGraphDialog :: JMeshSingularityGraphDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshSingularityGraphDialog :: ~JMeshSingularityGraphDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraphDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraphDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;
    display_mesh_edges = 1;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    size_t numedges = mesh->getSize(1);
    if( numedges < 1) return;

    JEdgeRenderPtr eAttrib;
    edgeColor.resize(numedges);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", eAttrib);
            edgeColor[i] = eAttrib->color;
        }
    }
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 0;

    nodePicker = meshViewer->getEntityPicker();
    nodePicker->setMesh(mesh);
    viewManager->attach( this );

    singularGraph.reset( new JMeshSingularityGraph );
    singularGraph->setMesh(mesh);
    singularGraph->clear();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraphDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( nodePicker == nullptr || mesh == nullptr ) return;

    JNodeSequence nodeSeq = nodePicker->getPickedNodes();
    if( nodeSeq.empty() ) return;

    if( nodeSeq.size() > 1 ) {
        cout << "Only one nodes must be picked at a time " << endl;
        return;
    }

    singularGraph->setQuadPath( nodeSeq[0] );
    displayGraph();

    nodePicker->clearAll();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraphDialog :: openNodeAttributesDialog()
{
    if( mesh == nullptr || singularGraph == nullptr) return;

    JNodeSequence nodes = singularGraph->getNodes();
    if( nodes.empty() ) return;
     
    if( nodeAttribsDialog == nullptr )
        nodeAttribsDialog.reset(new JNodeAttributesDialog( this ));

    nodeAttribsDialog->setViewManager(viewManager);
    nodeAttribsDialog->setMesh(mesh);
    nodeAttribsDialog->setNodes(nodes);

    nodeAttribsDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraphDialog :: openEdgeAttributesDialog()
{
    if( mesh == nullptr || singularGraph == nullptr) return;

    JEdgeSequence edges = singularGraph->getEdges();
    if( edges.empty() ) return;

    if( edgeAttribsDialog == nullptr )
        edgeAttribsDialog.reset(new JEdgeAttributesDialog( this ));

    edgeAttribsDialog->setViewManager(viewManager);
    edgeAttribsDialog->setMesh(mesh);

    edgeAttribsDialog->setEdges(edges);

    edgeAttribsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraphDialog :: getGraph()
{
    if( mesh == nullptr) return;

    JWaitCursor wCursor;
    wCursor.start();

    int val = stopAtJunctionCheckBox->isChecked();

    singularGraph->setStop(val);
    singularGraph->setAllPaths();

    displayGraph();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSingularityGraphDialog :: closeDialog()
{
    if( mesh ) mesh->deleteEdgeAttribute("Interface");
    display_mesh_edges = 1;
    displayGraph();

    close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraphDialog :: displayMesh()
{
    display_mesh_edges = displayMeshCheckBox->isChecked();
    displayGraph();
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraphDialog :: makeConnections()
{
    CheckBox( displayMeshCheckBox,  [=] { displayMesh(); });
    PushButton( edgeAttributesPushButton, [=] { openEdgeAttributesDialog(); });
    PushButton( nodeAttributesPushButton, [=] { openNodeAttributesDialog(); });
    PushButton( getGraphPushButton, [=] { getGraph(); });
    PushButton( closePushButton,    [=] { closeDialog(); });
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSingularityGraphDialog :: displayGraph()
{
    JColor highlightColor = JEntityColor::getColor("Red");
    JEdgeRenderPtr eAttrib;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", eAttrib);
            eAttrib->display = 1;
            if( edge->hasAttribute("Interface") ) {
                eAttrib->color  =  highlightColor;
                eAttrib->lineWidth  =  3;
            } else {
                eAttrib->color  =  edgeColor[i];
                eAttrib->lineWidth  =  1;
                eAttrib->display = display_mesh_edges;
            }
        }

    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
