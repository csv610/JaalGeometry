#include "TriMeshCleanupDialog.hpp"
#include "DelaunayMesherDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JTrimeshCleanupDialog :: JTrimeshCleanupDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    creaseAngle =  30.0;
    manifoldEdgeColor = nullptr;
    creaseAngleLineEdit->setText( QString::number(0.0));
}

///////////////////////////////////////////////////////////////////////////////

JTrimeshCleanupDialog :: ~JTrimeshCleanupDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr) return;
//  setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

}
///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: makeDelaunay()
{
    if( meshViewer == nullptr || mesh == nullptr) return;

    QString qstr;
    qstr = creaseAngleLineEdit->text();
    creaseAngle = qstr.toDouble();

    JDelaunayMesh2D delmesh;
    delmesh.setMesh(mesh);
    delmesh.setCreaseAngle( creaseAngle );

    JWaitCursor waitCursor;
    waitCursor.start();

    int numflips = delmesh.getRemeshed();
    waitCursor.stop();

    edgeflipsLineEdit->setText( QString::number(numflips));

//  meshViewer->displayAll(1);
    displayDelaunay();
}
///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: displayDelaunay()
{
    if( meshViewer == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    int ncount = 0;

    JDelaunayMesh2D delmesh;
    delmesh.setMesh(mesh);
    if( nonDelaunayEdgesCheckBox->isChecked() ) {
        ncount = delmesh.countNonDelaunayEdges();
        nonDelaunayEdgesLineEdit->setText( QString::number(ncount));
    }

    JEdgeColorPtr edgeColor(new JDelaunayEdgeColor);
    edgeColor->setMesh(mesh);
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: flipEdges()
{
    if( meshViewer == nullptr ) return;

    if( mesh == nullptr ) return;

    entityPicker = meshViewer->getEntityPicker();

    if( entityPicker == nullptr ) return;

    JSwapTriEdge swapedges;
    swapedges.setMesh(mesh, JEdgeSwap::NO_RULE_FLIP);
    swapedges.setCreaseAngle(180);
    JEdgeSequence edges = entityPicker->getPickedEdges();

    JEdgeSequence newEdges;
    JFaceSequence edgefaces, newFaces;
    for( size_t i = 0; i < edges.size(); i++)
        swapedges.applyAt( edges[i] );

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: subdivideEdges()
{
}
///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: removeEdges()
{
    if( meshViewer == nullptr ) return;
    if( mesh == nullptr ) return;

    if( mesh->getAdjTable(1,2) == 0)
        mesh->buildRelations(1,2);

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    entityPicker = meshViewer->getEntityPicker();
    if( entityPicker == nullptr ) return;
    JEdgeSequence edges = entityPicker->getPickedEdges();

    JTriDecimator decimator;
    decimator.setMesh(mesh);

    JNodeSequence newNodes;
    JEdgeSequence newEdges;
    JFaceSequence newFaces;

    cout << "Debug Exit" << __LINE__ << endl;

    /*
        for( size_t i = 0; i < edges.size(); i++) {
            decimator.collapse( edges[i] );
            newNodes = decimator.getNewNodes();
            newEdges = decimator.getNewEdges();
            newFaces = decimator.getNewFaces();
            meshViewer->addObjects( newNodes );
            meshViewer->addObjects( newEdges );
            meshViewer->addObjects( newFaces );
        }
    */

    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: displayNonManifoldEdges()
{
    if( meshViewer == nullptr ) return;
    JEdgeSequence edges = mesh->getTopology()->getNonManifoldEdges();
}

///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: below5DegreeNodes()
{
    if( meshViewer == nullptr || mesh == nullptr ) return;

    JTriDecimator decimator;
    decimator.setMesh(mesh);
    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: above7DegreeNodes()
{
    if( meshViewer == nullptr || mesh == nullptr ) return;

    JTriDecimator decimator;
    decimator.setMesh(mesh);
    decimator.removeHighDegreeNodes();
    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: displayIrregularNodes()
{
    if( meshViewer == nullptr || mesh == nullptr ) return;

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    size_t numnodes = mesh->getSize(0);

    JNodeRenderPtr attrib;

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        if( vertex->isActive() ) {
            vertex->getAttribute("Render", attrib);
            attrib->display  = 1;
            attrib->color[0] = 0.0;
            attrib->color[1] = 1.0;
            attrib->color[2] = 0.0;
            attrib->color[3] = 1.0;
            attrib->scale    = 1.0;
            attrib->glyph    = JNodeRender::NODE_AS_POINT;
        }
    }

    if( displayDegree5CheckBox->isChecked() ) {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vertex = mesh->getNodeAt(i);
            if( vertex->isActive() ) {
                int nd  = vertex->getNumRelations(2);
                int err = vertex->getAttribute("Render", attrib);
                if( err == 0 && nd < 5) {
                    attrib->display  = 1;
                    attrib->color[0] = 1.0;
                    attrib->color[1] = 0.0;
                    attrib->color[2] = 0.0;
                    attrib->color[3] = 1.0;
                    attrib->scale    = 1.0;
                    attrib->glyph    = JNodeRender::NODE_AS_SPHERE;
                }
            }
        }
    }

    if( displayDegree7CheckBox->isChecked() ) {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vertex = mesh->getNodeAt(i);
            if( vertex->isActive() ) {
                int nd = vertex->getNumRelations(2);
                int err = vertex->getAttribute("Render", attrib);
                if( err == 0 && nd > 7) {
                    attrib->display  = 1;
                    attrib->color[0] = 0.0;
                    attrib->color[1] = 0.0;
                    attrib->color[2] = 1.0;
                    attrib->color[3] = 1.0;
                    attrib->scale    = 1.0;
                    attrib->glyph    = JNodeRender::NODE_AS_SPHERE;
                }
            }
        }
    }

    meshViewer->updateBuffers( mesh );
}

///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: displayCollapsables()
{
//  mesh = meshViewer->getMesh();

    if( mesh == nullptr ) return;

    if( mesh->getAdjTable(1,2) == 0)
        mesh->buildRelations(1,2);

    if( mesh->getAdjTable(0,2) == 0)
        mesh->buildRelations(0,2);

    JTriDecimator decimator;
    decimator.setMesh(mesh);
    JColor redColor;
    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 1.0;

    JEdgeRenderPtr eAttrib;

    bool check567 = 1;

    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            if( decimator.isCollapsable(edge, check567) ) {
                edge->getAttribute("Render", eAttrib);
                eAttrib->display = 1;
                eAttrib->color   = redColor;
                eAttrib->scale   = 1.5;
            }
        }
    }

    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: displayFlippables()
{
    if( mesh == nullptr ) return;

    JSwapTriEdge swapper;
    swapper.setMesh(mesh);

    JColor redColor;
    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 1.0;

    JEdgeRenderPtr eAttrib;

    bool check567 = 1;

    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            if( swapper.isSwappable(edge, check567) ) {
                edge->getAttribute("Render", eAttrib);
                eAttrib->display = 1;
                eAttrib->color   = redColor;
                eAttrib->scale   = 1.5;
            }
        }
    }
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: displayCircumCircles()
{
    if( meshViewer == nullptr ) return;
//   mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return;

    /*
        if( delaunayCirclesCheckBox->isChecked() ) {
            if( delaunayViewer == nullptr )  delaunayViewer.reset(new JTriDelaunayViewer());
            delaunayViewer->setMesh(mesh);
            meshViewer->getViewManager()->attach(delaunayViewer);
        } else
            meshViewer->getViewManager()->detach(delaunayViewer);
    */

    meshViewer->getViewManager()->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: openAdvancingfront()
{
    if( advfrontDialog.get() == nullptr )
        advfrontDialog.reset(new JTriAdvancingfrontCleanupDialog(this));

    advfrontDialog->setViewManager( viewManager );
    advfrontDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: openLloydDialog()
{
    if( lloydRelaxationDialog == nullptr )
        lloydRelaxationDialog.reset(new JLloydRelaxationDialog(this));

    lloydRelaxationDialog->setViewManager( viewManager );
    lloydRelaxationDialog->setMesh( mesh );
    lloydRelaxationDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JTrimeshCleanupDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JTrimeshCleanupDialog :: makeConnections()
{
    PushButton( makeDelaunayPushButton,  [=] {makeDelaunay();});
    PushButton( displayIrregularNodesPushButton,  [=] {displayIrregularNodes();});
    PushButton( lloydRelaxationPushButton,  [=] {openLloydDialog(); });
    PushButton( below5DegreePushButton,  [=] {below5DegreeNodes();});
    PushButton( above7DegreePushButton,  [=] {above7DegreeNodes();});

    CheckBox(  nonDelaunayEdgesCheckBox, [=] {displayDelaunay(); });

    PushButton( closePushButton,  [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
