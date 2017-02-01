#include "QuadEditingDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JQuadEditingViewer :: JQuadEditingViewer()
{
    colors[0][0] = 1.0;
    colors[0][1] = 0.0;
    colors[0][2] = 0.0;
    colors[0][3] = 0.0;

    colors[1][0] = 0.0;
    colors[1][1] = 1.0;
    colors[1][2] = 0.0;
    colors[1][3] = 0.0;

    colors[2][0] = 0.0;
    colors[2][1] = 0.0;
    colors[2][2] = 1.0;
    colors[2][3] = 0.0;

    colors[3][0] = 1.0;
    colors[3][1] = 0.0;
    colors[3][2] = 1.0;
    colors[3][3] = 0.0;

    colors[4][0] = 0.0;
    colors[4][1] = 1.0;
    colors[4][2] = 1.0;
    colors[4][3] = 0.0;
}
///////////////////////////////////////////////////////////////////////////////

void JQuadEditingViewer :: draw()
{
    if( singularNodes.empty() ) return;

    JNodeSequence vneighs;
    qglviewer::Vec from, to;
    Point3D xyz;

    for( const JNodePtr &v0 : singularNodes) {
        JNode::getRelations(v0, vneighs);
        xyz   = v0->getXYZCoords();
        from[0] = xyz[0];
        from[1] = xyz[1];
        from[2] = xyz[2];

        int index = 0;
        for( const JNodePtr &v1 : vneighs) {
            glColor3fv( &colors[index%5][0] );
            xyz   = v1->getXYZCoords();
            to[0] = xyz[0];
            to[1] = xyz[1];
            to[2] = xyz[2];
            viewManager->drawArrow(from, to);
            index++;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

JQuadEditingDialog :: JQuadEditingDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    stop_picking = 0;
}

///////////////////////////////////////////////////////////////////////////////

JQuadEditingDialog :: ~JQuadEditingDialog()
{
}
///////////////////////////////////////////////////////////////////////////////
int JQuadEditingDialog :: isQuadMesh()
{
    if( mesh == nullptr ) return 0;

    int dim = mesh->getTopology()->getDimension();
    if( dim == 2 ) {
        int etype = mesh->getTopology()->getElementsType(2);
        if( etype != 4 ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("At present singlet operations are only for all quad mesh ");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) {
                return 0;
            }
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

void JQuadEditingDialog :: init()
{
    if( viewManager == nullptr ) return;

    quadEditViewer.reset(new JQuadEditingViewer());
    quadEditViewer->setName("QuadEditViewer");
    quadEditViewer->setViewManager( viewManager);

    viewManager->attach( quadEditViewer );
    viewManager->attach( this );

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );

};
///////////////////////////////////////////////////////////////////////////////

void JQuadEditingDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    picker = meshViewer->getEntityPicker();
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 0;

    if( !isQuadMesh() ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("At present singlet operations are only for all quad mesh ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
    mstMesher.setMesh(mesh);
    mesh->buildRelations(0,0);

    faceColor.reset(new JFacePartitionColor);
    edgeColor.reset(new JEdgePartitionColor);

//    motorGraph.setMesh(mesh);

    displaySingularNodes();
}

///////////////////////////////////////////////////////////////////////////////
void JQuadEditingDialog :: assignPartitionColors()
{
    if( mesh == nullptr ) return;

    JFaceColorPtr faceColor(new JFacePartitionColor);
    JEdgeColorPtr edgeColor(new JEdgePartitionColor);

//    JFaceColor::assign(mesh, faceColor);
//    JEdgeColor::assign(mesh, edgeColor);

    meshViewer->updateBuffers( mesh );
}
///////////////////////////////////////////////////////////////////////////////

void JQuadEditingDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( stop_picking ) return;

    if( picker == nullptr || mesh == nullptr ) return;

    JNodeSequence nodeSeq = picker->getPickedNodes();
    if( nodeSeq.empty() ) return;

    const JNodePtr &vtx = nodeSeq[0];

    bool found = 0;
    for( const JNodePtr &v : pickedNodes)
        if( v == vtx) found = 1;
    if( !found) pickedNodes.push_back(vtx);

    JWaitCursor waitCursor;
    waitCursor.start();

/*
    if( moveRadioButton->isChecked() && pickedNodes.size() == 2 )  {
        motorGraph.getPartitions( pickedNodes);
        assignPartitionColors();
        stop_picking = 1;
    }

    if(collapseRadioButton->isChecked() && pickedNodes.size() == 3) {
        motorGraph.getPartitions( pickedNodes);
        assignPartitionColors();
        stop_picking = 1;
    }
*/

    if( quadEditViewer)
        quadEditViewer->setSingularNodes(pickedNodes);
}
///////////////////////////////////////////////////////////////////////////////
void JQuadEditingDialog :: collapse()
{

}

///////////////////////////////////////////////////////////////////////////////

void JQuadEditingDialog :: displaySingularNodes()
{
    if( mesh == nullptr) return;

    JNodeRenderPtr nAttrib;
    size_t numnodes = mesh->getSize(0);

    JColor  defaultColor, degree3Color, degree5Color;
    defaultColor[0] = 0.0;
    defaultColor[1] = 1.0;
    defaultColor[2] = 0.0;
    defaultColor[3] = 0.0;

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        int err = vtx->getAttribute("Render", nAttrib);
        assert(!err);
        nAttrib->display =  1;
        nAttrib->glyph   =  0;
        nAttrib->color   =  defaultColor;
    }

    JNodeSequence irrnodes = mstMesher.getSingularNodes();

    degree3Color[0] = 1.0;
    degree3Color[1] = 0.0;
    degree3Color[2] = 0.0;
    degree3Color[3] = 0.0;

    degree5Color[0] = 0.0;
    degree5Color[1] = 0.0;
    degree5Color[2] = 1.0;
    degree5Color[3] = 0.0;

    for( const JNodePtr &vtx : irrnodes) {
        int err = vtx->getAttribute("Render", nAttrib);
        assert( !err);
        nAttrib->glyph   =  1;
        nAttrib->display =  1;

        if( vtx->getNumRelations(2) < 4)
            nAttrib->color   =  degree3Color;

        if( vtx->getNumRelations(2) > 4)
            nAttrib->color   =  degree5Color;
    }

    meshViewer->updateBuffers(mesh);
}
//////////////////////////////////////////////////////////////////////////////////////////
void JQuadEditingDialog :: closeDialog()
{
    if( viewManager == nullptr) return;
    if( quadEditViewer)
        viewManager->detach( quadEditViewer );
    viewManager->detach( this );
    this->close();
}

//////////////////////////////////////////////////////////////////////////////////////////
void JQuadEditingDialog :: clearAll()
{
    if( mesh == nullptr) return;

    JColor color;
    color[0] = 1.0;
    color[1] = 1.0;
    color[2] = 1.0;
    color[3] = 1.0;

    size_t numfaces = mesh->getSize(2);
    JFaceRenderPtr fattrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        int err = face->getAttribute("Render", fattrib);
        if( !err) fattrib->color = color;
    }

    color[0] = 0.1;
    color[1] = 0.1;
    color[2] = 0.1;
    color[3] = 1.0;

    size_t numedges = mesh->getSize(1);
    JEdgeRenderPtr eattrib;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        int err = edge->getAttribute("Render", eattrib);
        if( !err) eattrib->color = color;
    }

    pickedNodes.clear();
    picker->clearAll();
    stop_picking = 0;
    quadEditViewer->setSingularNodes(pickedNodes);

    meshViewer->updateBuffers(mesh);
}

//////////////////////////////////////////////////////////////////////////////////////////

void JQuadEditingDialog :: makeConnections()
{
    PushButton( clearPushButton, [=] { clearAll();});
    PushButton( closePushButton, [=] { closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
