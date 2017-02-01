#include "MedialAxis2DDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMedialAxis2DDialog :: JMedialAxis2DDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMedialAxis2DDialog :: ~JMedialAxis2DDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMedialAxis2DDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JMedialAxis2DDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    if( mrender ) {
        mrender->displayEntity[2] = 0;
    }

    JEdgeSequence internalEdges;
    mesh->getTopology()->getInternal(internalEdges);
    JEdgeRenderPtr eAttrib;
    for( const JEdgePtr &edge : internalEdges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->display = 0;
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMedialAxis2DDialog :: getAxis()
{
    if( mesh == nullptr) return;

    JDelaunayMesh2D  del;
    del.setMesh(mesh);
    medialAxis = del.getMedialAxis();

    if( medialAxis == nullptr) return;

    meshViewer->addObject(medialAxis);

    JColor  r = JEntityColor::getColor("Red");
    JColor  g = JEntityColor::getColor("Green");
    JColor  b = JEntityColor::getColor("Blue");

    medialAxis->buildRelations(0,0);
    size_t numnodes = medialAxis->getSize(0);

    JNodeRenderPtr vAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = medialAxis->getNodeAt(i);
        vtx->getAttribute("Render", vAttrib);
        int degree = vtx->getNumRelations(0);
        switch( degree )
        {
        case 1:
            vAttrib->color = r;
            vAttrib->glyph = 1;
            break;
        case 2:
            vAttrib->color = g;
            break;
        default:
            vAttrib->color = b;
            vAttrib->glyph = 1;
            break;
        }
    }
    meshViewer->updateBuffers(medialAxis);
}

///////////////////////////////////////////////////////////////////////////////
void JMedialAxis2DDialog :: openNodeAttribsDialog()
{
    if( mesh == nullptr) return;

    if( nodeAttribDialog.get() == nullptr )
        nodeAttribDialog.reset(new JNodeAttributesDialog( this ));

    nodeAttribDialog->setViewManager(viewManager);

    if( medialAxis ) {
        nodeAttribDialog->setMesh(medialAxis);

        JNodeSequence seq = medialAxis->getNodes();
        nodeAttribDialog->setNodes(seq);

        nodeAttribDialog->show();
        this->hide();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMedialAxis2DDialog :: openEdgeAttribsDialog()
{
    if( mesh == nullptr) return;

    if( edgeAttribDialog.get() == nullptr )
        edgeAttribDialog.reset(new JEdgeAttributesDialog( this ));

    edgeAttribDialog->setViewManager(viewManager);
    if( medialAxis ) {
        edgeAttribDialog->setMesh(medialAxis);

        JEdgeSequence seq = medialAxis->getEdges();
        edgeAttribDialog->setEdges(seq);

        edgeAttribDialog->show();
        this->hide();
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMedialAxis2DDialog :: displayDelaunay()
{
    if( mesh == nullptr) return;

    bool val = delaunayMeshCheckBox->isChecked();

    if( delmesh ) {
        meshViewer->removeObject(delmesh);
        delmesh.reset();
        if( val == 0) return;
        return;
    }

    JMeshPtr tmpcopy = mesh->deepCopy();
    JDelaunayMesh2D   delaunay;
    delaunay.setMesh(tmpcopy);
    delmesh = delaunay.getSimpleMesh();
    meshViewer->addObject(delmesh);

    JMeshRenderPtr mrender;
    delmesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 1;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;

    meshViewer->updateBuffers(delmesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMedialAxis2DDialog :: closeDialog()
{
    if( delmesh ) {
        meshViewer->removeObject(delmesh);
        delmesh.reset();
    }
    parentWidget()->show();
    close();
}

///////////////////////////////////////////////////////////////////////////////

void JMedialAxis2DDialog :: makeConnections()
{
    connect( nodeAttribsPushButton,  SIGNAL( clicked() ), this, SLOT( openNodeAttribsDialog() ));
    connect( edgeAttribsPushButton,  SIGNAL( clicked() ), this, SLOT( openEdgeAttribsDialog() ));
    connect( delaunayMeshCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( displayDelaunay() ));
    connect( generatePushButton,  SIGNAL( clicked() ), this, SLOT( getAxis() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
