#include "QuadMeshCleanupDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JQuadMeshCleanupDialog :: JQuadMeshCleanupDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JQuadMeshCleanupDialog :: ~JQuadMeshCleanupDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMeshCleanupDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshCleanupDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

}
///////////////////////////////////////////////////////////////////////////////

void JQuadMeshCleanupDialog :: searchSinglets()
{
    int n = 0;
    numSingletsLineEdit->setText( QString::number(n));

    if( mesh == nullptr) return;

    if( qSinglet == nullptr) qSinglet.reset( new JSinglet);
    qSinglet->setMesh(mesh);

    JNodeSequence singlets  = qSinglet->getSinglets();
    numSingletsLineEdit->setText( QString::number(singlets.size()));
    mesh->enumerate(0);

    size_t numnodes = mesh->getSize(0);
    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", nAttrib);
            nAttrib->display = 1;
            nAttrib->glyph   = 0;
        }
    }

    JColor highlightColor = JEntityColor::getColor("Red");

    for( const JNodePtr &vtx : singlets) {
        vtx->getAttribute("Render", nAttrib);
        nAttrib->color   =  highlightColor;
        nAttrib->glyph   = 1;
        nAttrib->display = 1;
    }

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMeshCleanupDialog :: removeSinglets()
{
    if( qSinglet == nullptr) return;
    qSinglet->setMesh(mesh);
//  qSinglet->removeAll();
    qSinglet->mergeAll();

    searchSinglets();
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMeshCleanupDialog :: searchDoublets()
{
    int n = 0;
    if( mesh ) {
        if( qDoublet == nullptr) qDoublet.reset( new JDoublet);
        qDoublet->setMesh(mesh);
        n = qDoublet->getSize();
    }
    numDoubletsLineEdit->setText( QString::number(n));

    if( mesh == nullptr) return;

    size_t numnodes = mesh->getSize(0);

    JColor highlightColor = JEntityColor::getColor("Red");

    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", nAttrib);
            nAttrib->display = 0;
            if( !vtx->isBoundary()  && vtx->getNumRelations(2) == 2) {
                nAttrib->color =  highlightColor;
                nAttrib->glyph   = 1;
                nAttrib->display = 1;
            }
        }
    }
    meshViewer->updateBuffers(mesh);

}
///////////////////////////////////////////////////////////////////////////////
void JQuadMeshCleanupDialog :: removeDoublets()
{
    if( qDoublet == nullptr) return;
    qDoublet->removeAll();

    int n = qDoublet->getSize();
    numDoubletsLineEdit->setText( QString::number(n));

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMeshCleanupDialog :: removeDiamonds()
{

}
///////////////////////////////////////////////////////////////////////////////

void JQuadMeshCleanupDialog :: searchDiamonds()
{
    int type = 33;

    numDiamondsLineEdit->setText( QString::number(0));

    if( mesh ) {
        if( qDiamond == nullptr) qDiamond.reset( new JDiamond);
        qDiamond->setMesh(mesh);
    }

    if( mesh == nullptr) return;

    size_t numfaces = mesh->getSize(2);

    JColor highlightColor = JEntityColor::getColor("Red");

    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", fAttrib);
            fAttrib->display = 0;
        }
    }

    JFaceSequence diamonds = qDiamond->getDiamonds(type);
    int nd = diamonds.size();
    numDiamondsLineEdit->setText( QString::number(nd));

    if( nd ) {
        for( const JFacePtr &face : diamonds)
            face->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
        fAttrib->color = highlightColor;
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JQuadMeshCleanupDialog :: swapEdges()
{

    mesh->buildRelations(0,2);

    size_t numnodes = mesh->getSize(0);
    JNodeSequence highdegree;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->getNumRelations(2) > 5) highdegree.push_back(vtx);
    }
    if (highdegree.empty() ) return;

    JSwapQuadEdge qswap;
    qswap.setMesh(mesh);
    for( const JNodePtr &vtx : highdegree ) {
        qswap.applyAt(vtx);
        meshViewer->updateBuffers(mesh);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMeshCleanupDialog :: openDualDialog()
{
    if( dualDialog == nullptr)
        dualDialog.reset( new JQuadMeshDualsDialog(this));

    dualDialog->setViewManager(viewManager);
    dualDialog->setMesh(mesh);
    dualDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMeshCleanupDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JQuadMeshCleanupDialog :: makeConnections()
{
    PushButton( dualOpsPushButton,         [=] {openDualDialog();});
    PushButton( swapEdgesPushButton,       [=] {swapEdges();});
    PushButton( searchSingletsPushButton,  [=] {searchSinglets();});
    PushButton( searchDoubletsPushButton,  [=] {searchDoublets();});
    PushButton( searchDiamondsPushButton,  [=] {searchDiamonds();});
    PushButton( removeSingletsPushButton,  [=] {removeSinglets();});
    PushButton( removeDoubletsPushButton,  [=] {removeDoublets();});
    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////

