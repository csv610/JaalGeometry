#include "DoubletDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JDoubletDialog :: JDoubletDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
//  doubletColor.reset(new DoubletColor);
}
///////////////////////////////////////////////////////////////////////////////

JDoubletDialog :: ~JDoubletDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JDoubletDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JDoubletDialog :: assignColor()
{
    if( mesh == nullptr) return;

    size_t numnodes = mesh->getSize(0);

    JColor highlightColor, defaultColor;

    highlightColor[0] = 1.0;
    highlightColor[1] = 0.0;
    highlightColor[2] = 0.0;
    highlightColor[3] = 1.0;

    defaultColor[0] = 0.0;
    defaultColor[1] = 1.0;
    defaultColor[2] = 0.0;
    defaultColor[3] = 1.0;

    JNodeRenderPtr nAttrib;

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Render", nAttrib);
        if( !vtx->isBoundary()  && vtx->getNumRelations(2)  == 2) {
            nAttrib->color =  highlightColor;
            nAttrib->glyph = 1;
        } else {
            nAttrib->color =  defaultColor;
            nAttrib->glyph = 0;
        }
    }
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JDoubletDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    mesh->buildRelations(0,2);
//  meshViewer->getNodeDraw()->setColorMethod( doubletColor );

    doublet.setMesh(mesh);

    int nSize  = doublet.getDoublets().size();
    numDoubletsLineEdit->setText( QString::number(nSize) );

    assignColor();
}

///////////////////////////////////////////////////////////////////////////////

void JDoubletDialog :: setColor()
{
    if( meshViewer == nullptr ) return;
    if( mesh == nullptr ) return ;

    /*
    JNodeColor::assign(mesh, doubletColor.get() );
         mesh->buildRelations(0,2);
         meshViewer->getDrawNode()->setColorMethod( doubletColor.get() );
    */
}

///////////////////////////////////////////////////////////////////////////////

void JDoubletDialog :: removeAll()
{
    if( meshViewer == nullptr ) return;

    doublet.removeAll();

    int nSize = doublet.getDoublets().size();
    numDoubletsLineEdit->setText( QString::number(nSize) );

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JDoubletDialog :: makeConnections()
{
    connect( colorPushButton,  SIGNAL( clicked() ), this, SLOT( setColor() ));
    connect( removeAllPushButton, SIGNAL( clicked() ), this, SLOT( removeAll() ));
}

///////////////////////////////////////////////////////////////////////////////
