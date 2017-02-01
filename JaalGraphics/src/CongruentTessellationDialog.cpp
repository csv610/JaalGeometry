#include "CongruentTessellationDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
JCongruentViewer  :: JCongruentViewer()
{
    mesh = nullptr;
    firstNode   = nullptr;
}
///////////////////////////////////////////////////////////////////////////////

void JCongruentViewer  :: draw()
{
    int nSize = paramCorners.size();

    glColor3f( 1.0, 0.8, 0.8);
    glLineWidth(2.0);
    Point3D p0,p1;


    glBegin(GL_LINES);
    for( int i = 0; i < nSize; i++) {
        p0 = paramCorners[i%nSize]->getXYZCoords();
        p1 = paramCorners[(i+1)%nSize]->getXYZCoords();
        glVertex3f( p0[0], p0[1],  0.0001);
        glVertex3f( p1[0], p1[1],0.0001);
    }
    glEnd();

    nSize = srcdstNodes.size()/2;

    if( firstNode ) {
        glBegin(GL_LINES);
        p0 = firstNode->getXYZCoords();
        p1 = paramCorners[0]->getXYZCoords();
        glVertex3f( p0[0], p0[1],  0.0001);
        glVertex3f( p1[0], p1[1],0.0001);
        glEnd();
    }

}

///////////////////////////////////////////////////////////////////////////////

JCongruentTessellationDialog :: JCongruentTessellationDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JCongruentTessellationDialog :: ~JCongruentTessellationDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JCongruentTessellationDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

//  mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;

    congruentViewer.reset( new JCongruentViewer );
    congruentViewer->setName("CongruentViewer");
    congruentViewer->setViewManager(viewManager);
    viewManager->attach( congruentViewer );
    viewManager->setMouseTracking(1);

    if( trianglesRadioButton->isChecked() )
        setBoundingTriangle();

    if( quadsRadioButton->isChecked() )
        setBoundingSquare();

    /*
        picker = meshViewer->getEntityPicker();
        if( picker ) {
            picker->setPickableEntity(0);
            picker->setMode(1);
        }
    */
    viewManager->attach( this );

    mesh->getTopology()->getBoundary(boundedges);
    JEdgeTopology::getChain( boundedges);
    assert( JEdgeTopology::isClosed( boundedges));
    assert( JEdgeTopology::isTopologicalSimple( boundedges));
}
///////////////////////////////////////////////////////////////////////////////

void JCongruentTessellationDialog :: setBoundingSquare()
{
    if( mesh == nullptr ) return;
    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
}

///////////////////////////////////////////////////////////////////////////////

void JCongruentTessellationDialog :: setBoundingTriangle()
{
    /*
        if( mesh == nullptr ) return;
        BoundingBox box = mesh->getGeometry()->getBoundingBox();

        double  base = box.getLength(1)*tan(M_PI/6.0) + 0.5*box.getLength(0);
        const Point3D &pC = box.getCenter();

        double xmin  = pC[0] - base;
        double xmax  = pC[0] + base;
        double ymin  = pC[1] - 0.5*box.getLength(1);
        double ymax  = ymin + fabs(xmin)*tan(M_PI/3.0);

        paramCorners.clear();

        paramCorners.resize(3);

        Point3D xyz;
        xyz[0] = xmin;
        xyz[1] = ymin;
        xyz[2] = pC[2];
        paramCorners[0] = JNode::newObject();
        paramCorners[0]->setXYZCoords(xyz);

        xyz[0] = xmax;
        xyz[1] = ymin;
        xyz[2] = pC[2];
        paramCorners[1] = JNode::newObject();
        paramCorners[1]->setXYZCoords(xyz);

        xyz[0] = pC[0];
        xyz[1] = ymax;
        xyz[2] = pC[2];
        paramCorners[2] = JNode::newObject();
        paramCorners[2]->setXYZCoords(xyz);
        congruentViewer->setParamCorners( paramCorners );

        int numSegments = boundedges.size();
        int numPoints   = numSegments/3;
        JNodeSequence nodes1, nodes2, nodes3;
        JEdge::generate_linear_nodes( paramCorners[0], paramCorners[1], numPoints, nodes1);
        JEdge::generate_linear_nodes( paramCorners[1], paramCorners[2], numPoints, nodes2);
        JEdge::generate_linear_nodes( paramCorners[2], paramCorners[3], numSegments-2*numPoints, nodes3);
    */
}

///////////////////////////////////////////////////////////////////////////////
void JCongruentTessellationDialog :: generate()
{

}
///////////////////////////////////////////////////////////////////////////////
void JCongruentTessellationDialog :: setParamShape()
{
}

///////////////////////////////////////////////////////////////////////////////

void JCongruentTessellationDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( !this->isVisible() || viewManager == nullptr ) return;

    if( picker ) {
        JNodeSequence nodeSeq = picker->getPickedNodes();
        if( nodeSeq.empty() ) return;
        if( !nodeSeq[0]->hasAttribute("TargetPos") ) {
            cout << "Warning: No Constraint target selected for movement " << endl;
            picker->clearAll();
            return;
        }

        if(!paramCorners.empty() ) {
            congruentViewer->setFirstNode( nodeSeq[0] );
        }
        meshViewer->refreshDisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JCongruentTessellationDialog :: closeDialog()
{
    viewManager->detach( congruentViewer );
    congruentViewer.reset();
    viewManager->detach( this );
    this->close();

}
///////////////////////////////////////////////////////////////////////////////

void JCongruentTessellationDialog :: makeConnections()
{
    connect( generatePushButton, SIGNAL( clicked() ), this, SLOT( generate() ));
    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
