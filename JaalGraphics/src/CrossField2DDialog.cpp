#include "CrossField2DDialog.hpp"
#include <igl/jet.h>

///////////////////////////////////////////////////////////////////////////////

JCrossField2DDialog :: JCrossField2DDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    crossField.reset(new JCrossField2D);
    vecScaleLineEdit->setText( QString::number(1.0));
}

///////////////////////////////////////////////////////////////////////////////

JCrossField2DDialog :: ~JCrossField2DDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JCrossField2DDialog :: setLineWidth()
{
    int lineWidth = lineWidthSpinBox->value();
    JEdgeRenderPtr eAttrib;
    if( nodeVecField ) {
        size_t numedges = nodeVecField->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = nodeVecField->getEdgeAt(i);
            edge->getAttribute("Render", eAttrib);
            eAttrib->lineWidth = lineWidth;
        }
        meshViewer->updateBuffers( nodeVecField );
    }

    if( edgeVecField ) {
        size_t numedges = edgeVecField->getSize(1);
        for( size_t i = 0; i < numedges; i++) {
            const JEdgePtr &edge = edgeVecField->getEdgeAt(i);
            edge->getAttribute("Render", eAttrib);
            eAttrib->lineWidth = lineWidth;
        }
        meshViewer->updateBuffers( edgeVecField );
    }
}
///////////////////////////////////////////////////////////////////////////////

void JCrossField2DDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JCrossField2DDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    crossField->setMesh(mesh);

    JMeshQuality mq;
    mq.setMesh(mesh);
    veclen = 0.50*mq.getEdgesLength().getMedian();

}

///////////////////////////////////////////////////////////////////////////////
void JCrossField2DDialog :: setVecLength()
{
    QString qstr = vecScaleLineEdit->text();
    double scale = qstr.toDouble();

    double angle;
    Point3D  head, tail, arrow;
    head[0] = 0.0;
    head[1] = 0.0;
    head[2] = 0.0;

    arrow[0] = 0.0;
    arrow[1] = 0.0;
    arrow[2] = 0.0;

    JNodeSequence nodes;
    JEdgeSequence edges;
    JNodePtr v0, v1;

    if( nodeVecField ) {
        nodes = nodeVecField->getNodes();
        edges = nodeVecField->getEdges();
        size_t numnodes = mesh->getSize(0);
        size_t edgeid = 0;
        size_t nodeid = 0;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            int err = vtx->getAttribute("CrossFieldAngle", angle);
            if( !err) {
                tail = vtx->getXYZCoords();
                v0 = nodes[nodeid++];
                v0->setXYZCoords(tail);

                for( int j = 0; j < 4; j++) {
                    v1 = nodes[nodeid++];
                    double len = JNodeGeometry::getLength(v0,v1);

                    head[0] = tail[0] + scale*len*cos(angle + j*M_PI/2.0);
                    head[1] = tail[1] + scale*len*sin(angle + j*M_PI/2.0);
                    v1->setXYZCoords(head);

                    v1 = nodes[nodeid++];
                    arrow[0] = head[0] + 0.25*scale*len*cos( angle + j*M_PI/2.0 + 165*M_PI/180.0);
                    arrow[1] = head[1] + 0.25*scale*len*sin( angle + j*M_PI/2.0 + 165*M_PI/180.0);
                    v1->setXYZCoords(arrow);

                    v1 = nodes[nodeid++];
                    arrow[0] = head[0] + 0.25*scale*len*cos( angle + j*M_PI/2.0 + 195*M_PI/180.0);
                    arrow[1] = head[1] + 0.25*scale*len*sin( angle + j*M_PI/2.0 + 195*M_PI/180.0);
                    v1->setXYZCoords(arrow);
                }
            }
        }
        meshViewer->updateBuffers(nodeVecField);
    }



    /*
        if( edgeVecField ) {
            size_t numedges = edgeVecField->getSize(1);
            for( size_t i = 0; i < numedges; i++) {
                const JEdgePtr &edge = edgeVecField->getEdgeAt(i);
                const JNodePtr &v0   = edge->getNodeAt(0);
                const JNodePtr &v1   = edge->getNodeAt(1);
                const Point3D  &p0   = v0->getXYZCoords();
                const Point3D  &p1   = v1->getXYZCoords();
                double dx = p1[0] - p0[0];
                double dy = p1[1] - p0[1];
                pn[0] = p0[0] + scale*dx;
                pn[1] = p0[1] + scale*dy;
                v1->setXYZCoords(pn);
            }
            meshViewer->updateBuffers(edgeVecField);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JCrossField2DDialog :: drawVectors()
{
    bool val;
    val = edgeVectorsCheckBox->isChecked();
    if( edgeVecField ) {
        edgeVecField->setActiveBit(val);
    }

    val = nodeVectorsCheckBox->isChecked();
    if( nodeVecField ) {
        nodeVecField->setActiveBit(val);
    }

    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JCrossField2DDialog :: generateField()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    crossField->genVecField();

    double angle;
    Point3D  head, tail, arrow;
    head[0] = 0.0;
    head[1] = 0.0;
    head[2] = 0.0;

    arrow[0] = 0.0;
    arrow[1] = 0.0;
    arrow[2] = 0.0;

    JNodePtr v0, v1, v2;
    JEdgePtr vedge;
    vector<JColor> colors(4);
    colors[0] = JEntityColor::getColor("Red");
    colors[1] = JEntityColor::getColor("Green");
    colors[2] = JEntityColor::getColor("Blue");
    colors[3] = JEntityColor::getColor("Purple");

    JNodeSequence nodes;
    JEdgeSequence edges;

    edgeVecField = JMesh::newObject();
    size_t numedges = mesh->getSize(1);

    nodes = JNode::newObjects(13*numedges);
    edges = JEdge::newObjects(12*numedges);

    size_t nodeid = 0;
    size_t edgeid = 0;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        int err = edge->getAttribute("CrossFieldAngle", angle);
        if( !err) {
            tail = JEdgeGeometry::getMidPoint(edge);
            v0 = nodes[nodeid++];
            v0->setXYZCoords(tail);

            for( int j = 0; j < 4; j++) {
                head[0] = tail[0] + veclen*cos(angle + j*M_PI/2.0);
                head[1] = tail[1] + veclen*sin(angle + j*M_PI/2.0);
                v1 = nodes[nodeid++];
                v1->setXYZCoords(head);
                vedge = edges[edgeid++];
                vedge->setNodes(v0,v1);
                vedge->setAttribute("Color", colors[j] );

                arrow[0] = head[0] + 0.25*veclen*cos( angle + j*M_PI/2.0 + 165*M_PI/180.0);
                arrow[1] = head[1] + 0.25*veclen*sin( angle + j*M_PI/2.0 + 165*M_PI/180.0);
                v2 = nodes[nodeid++];
                v2->setXYZCoords(arrow);
                vedge = edges[edgeid++];
                vedge->setNodes(v1,v2);
                vedge->setAttribute("Color", colors[j] );

                arrow[0] = head[0] + 0.25*veclen*cos( angle + j*M_PI/2.0 + 195*M_PI/180.0);
                arrow[1] = head[1] + 0.25*veclen*sin( angle + j*M_PI/2.0 + 195*M_PI/180.0);
                v2 = nodes[nodeid++];
                v2->setXYZCoords(arrow);
                vedge = edges[edgeid++];
                vedge->setNodes(v1,v2);
                vedge->setAttribute("Color", colors[j] );
            }
        }
    }
    edgeVecField->addObjects(nodes);
    edgeVecField->addObjects(edges);
    meshViewer->addObject(edgeVecField);

    nodeVecField = JMesh::newObject();
    size_t numnodes = mesh->getSize(0);

    nodes = JNode::newObjects(13*numnodes);
    edges = JEdge::newObjects(12*numnodes);
    edgeid = 0;
    nodeid = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        int err = vtx->getAttribute("CrossFieldAngle", angle);
        if( !err) {
            tail = vtx->getXYZCoords();
            v0 = nodes[nodeid++];
            v0->setXYZCoords(tail);

            for( int j = 0; j < 4; j++) {
                head[0] = tail[0] + veclen*cos(angle + j*M_PI/2.0);
                head[1] = tail[1] + veclen*sin(angle + j*M_PI/2.0);
                v1 = nodes[nodeid++];
                v1->setXYZCoords(head);
                vedge = edges[edgeid++];
                vedge->setNodes(v0,v1);
                vedge->setAttribute("Color", colors[j] );

                arrow[0] = head[0] + 0.25*veclen*cos( angle + j*M_PI/2.0 + 165*M_PI/180.0);
                arrow[1] = head[1] + 0.25*veclen*sin( angle + j*M_PI/2.0 + 165*M_PI/180.0);
                v2 = nodes[nodeid++];
                v2->setXYZCoords(arrow);
                vedge = edges[edgeid++];
                vedge->setNodes(v1,v2);
                vedge->setAttribute("Color", colors[j] );

                arrow[0] = head[0] + 0.25*veclen*cos( angle + j*M_PI/2.0 + 195*M_PI/180.0);
                arrow[1] = head[1] + 0.25*veclen*sin( angle + j*M_PI/2.0 + 195*M_PI/180.0);
                v2 = nodes[nodeid++];
                v2->setXYZCoords(arrow);
                vedge = edges[edgeid++];
                vedge->setNodes(v1,v2);
                vedge->setAttribute("Color", colors[j] );
            }
        }
    }

    nodeVecField->addObjects(nodes);
    nodeVecField->addObjects(edges);
    meshViewer->addObject(nodeVecField);

    JColor eColor;
    JEdgeRenderPtr eAttrib;

    numedges = edgeVecField->getSize(1);
    for(size_t i  = 0; i < numedges; i++) {
        const JEdgePtr &edge = edgeVecField->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        edge->getAttribute("Color",  eColor);
        eAttrib->color = eColor;
        eAttrib->lineWidth = 2;
    }
    meshViewer->updateBuffers( edgeVecField );

    numedges = nodeVecField->getSize(1);
    for(size_t i  = 0; i < numedges; i++) {
        const JEdgePtr &edge = nodeVecField->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        edge->getAttribute("Color",  eColor);
        eAttrib->color = eColor;
        eAttrib->lineWidth = 2;
    }
    meshViewer->updateBuffers( nodeVecField );

    drawVectors();
}
///////////////////////////////////////////////////////////////////////////////

void JCrossField2DDialog :: alignEdges()
{
}

///////////////////////////////////////////////////////////////////////////////

void JCrossField2DDialog :: showScalarField()
{
    QString qstr = scalarFieldComboBox->currentText();
    string str   = StdString(qstr);

    vector<double> field;
    if( str == "Cos(4t)")
        field = crossField->getCos4tField();

    if( str == "Sin(4t)")
        field = crossField->getSin4tField();

    if( field.empty() ) return;
    vector<JColor> colors;
    JColorMap::jet(field,colors);
    JNodeRenderPtr vAttrib;
    JColor  color;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v = mesh->getNodeAt(i);
        v->getAttribute("Render", vAttrib);
        vAttrib->color = colors[i];
    }

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
//  mrender->setSurfaceShade(JRender::FLAT_SHADE);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JCrossField2DDialog :: closeDialog()
{
    if( nodeVecField ) {
        meshViewer->removeObject(nodeVecField);
        nodeVecField.reset();
    }

    if( edgeVecField ) {
        meshViewer->removeObject(edgeVecField);
        edgeVecField.reset();
    }
    meshViewer->refreshDisplay();
    close();
}
///////////////////////////////////////////////////////////////////////////////
void JCrossField2DDialog :: makeConnections()
{
    connect( nodeVectorsCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( drawVectors() ));
    connect( edgeVectorsCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( drawVectors() ));
    connect( faceVectorsCheckBox, SIGNAL( stateChanged(int) ), this, SLOT( drawVectors() ));

    connect( scalarFieldComboBox, SIGNAL( activated(int) ), this, SLOT( showScalarField() ));
    connect( vecScaleLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setVecLength() ));
    connect( lineWidthSpinBox, SIGNAL( valueChanged(int) ), this, SLOT( setLineWidth() ));

    connect( generatePushButton,  SIGNAL( clicked() ), this, SLOT( generateField() ));
    connect( alignEdgesPushButton, SIGNAL( clicked() ), this, SLOT( alignEdges() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
