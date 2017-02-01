#include "MeshConstraintsDialog.hpp"
#include "basic_math.hpp"

using namespace std;

int JMeshConstraintsDialog :: groupID = 1;
///////////////////////////////////////////////////////////////////////////////

JMeshConstraintsDialog :: JMeshConstraintsDialog( QWidget *parent) : QDialog(parent)
{
    viewManager = nullptr;
    picker      = nullptr;
    setupUi(this);
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

JMeshConstraintsDialog :: ~JMeshConstraintsDialog()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: init()
{
    assert( viewManager );

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    assert(c);

    // Get the mesh view, it containts entity picker..
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    // Only the nodes are pickable and you can pick many nodes ...
    picker = meshViewer->getEntityPicker();
    picker->setMode(2);

    /*
        picker->setPickableEntity(0);
        if( constraintsViewer == nullptr ) {
            constraintsViewer.reset( new JMeshConstraintsViewer() );
            constraintsViewer->setName("MeshConstraintsViewer");
        }
        viewManager->attach( constraintsViewer.get() );
    */

    viewManager->setMouseTracking(1);
    setMesh( meshViewer->getCurrentMesh() );

    // By default all the boundary nodes are constrainted.
    addBoundaryConstraints();
    countConstraints();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshConstraintsDialog :: showEvent( QShowEvent *event)
{
    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = 0;
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    picker->setMesh(mesh);
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}
///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: getPolygon()
{
    closedLassoCheckBox->setChecked(1);
    lassoPoints = polygonDialog->getPoints();

    if( lassoPoints.empty() ) return;
    lassoPoints.push_back( lassoPoints.front() );
    addWithinLasso();
    countConstraints();
    constraintsViewer->setLasso( lassoPoints );

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: add2Group(const JNodePtr &vtx)
{
    if( vtx == nullptr ) return;
    vtx->setAttribute("Constraint", groupID);
    if( constraintsViewer ) constraintsViewer->setConstraint(vtx, 1);
}

////////////////////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: parameterizeLasso()
{
    vector<Point3D> lassoCtrlPoints;

    size_t nSize = lassoPoints.size();
    if( nSize < 2) return;

    // Calculate the length of the lasso
    double len, lassolen = 0.0;
    for( size_t i = 0; i < nSize; i++) {
        const Point3D &p0 = lassoPoints[i];
        const Point3D &p1 = lassoPoints[(i+1)%nSize];
        len = JMath::length(p0, p1);
        lassolen += len;
    }

    // Lasso is always a closed curve. Check the distance between the last and  last
    // point, if the distance is not zero, insert one more in the last and make it
    // coincide with the first point.
    Point3D pfirst = lassoPoints.front();
    Point3D plast  = lassoPoints.back();
    len = JMath::length(pfirst, plast);
    if( len > 1.0E-05*lassolen) lassoPoints.push_back( pfirst );

    nSize = lassoPoints.size();

    // Calculate the parametric value ([0,1])  of each point on the lasso.
    vector<double> param(nSize);
    len = 0.0;
    param[0] = 0.0;
    for( size_t i = 0; i < nSize-1; i++) {
        const Point3D &p0 = lassoPoints[i];
        const Point3D &p1 = lassoPoints[i+1];
        len += JMath::length(p0, p1);
        param[i+1] = len/lassolen;
    }
    param[nSize-1] = 1.0;

    // Sample some points (uniformly ) on the lasso.
    int nparam = std::min(100,(int)(nSize-1));
    double dl = 1.0/(double)nparam;

    lassoCtrlPoints.clear();

    lassoCtrlPoints.push_back( lassoPoints[0] );
    Point3D xyz;
    for( int j = 1; j < nparam-1; j++) {
        double p = j*dl;
        for( size_t i = 0; i < nSize-1; i++) {
            if( param[i] >= p) {
                double t0 = param[i];
                double t1 = param[i+1];
                assert(t1 > t0);
                double t = (p - t0)/(t1 - t0);
                const Point3D &p0 = lassoPoints[i];
                const Point3D &p1 = lassoPoints[i+1];
                xyz[0] = (1-t)*p0[0] + t*p1[0];
                xyz[1] = (1-t)*p0[1] + t*p1[1];
                xyz[2] = (1-t)*p0[2] + t*p1[2];
                lassoCtrlPoints.push_back(xyz);
                break;
            }
        }
    }

    lassoCtrlPoints.push_back( lassoPoints[0] );
    lassoPoints = lassoCtrlPoints;

    constraintsViewer->setLasso(lassoCtrlPoints);
}

///////////////////////////////////////////////////////////////////////////////

int JMeshConstraintsDialog :: getSign( const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    double x1 = 0.0;
    double y1 = 0.0;

    double x2 = p1[0]-p0[0];
    double y2 = p1[1]-p0[1];

    double x3 = p2[0]-p0[0];
    double y3 = p2[1]-p0[1];

    double a  = 0.5*(-x2*y1+x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3);
    if( a < 0) return -1;
    if( a > 0) return  1;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

size_t JMeshConstraintsDialog :: countConstraints()
{
    set<int> aset;

    size_t nSize = mesh->getSize(0);
    int gid, nCount = 0;
    for( size_t i = 0; i < nSize; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( !vtx->isActive() ) continue;
        if( vtx->hasAttribute("Constraint") ) {
            vtx->getAttribute("Constraint", gid);
            aset.insert(gid);
            nCount++;
        }
    }

    groupID = *boost::max_element( aset ) + 1;

    numConstraintNodesLineEdit->setText( QString::number(nCount) );
    numSetsLineEdit->setText( QString::number(aset.size() ));
    return nCount;
}

///////////////////////////////////////////////////////////////////////////////

bool JMeshConstraintsDialog :: isWithinLasso( const Point3D &p0)
{
    static int nCount = 0;
    size_t nSize = lassoPoints.size();


    double angle = 0.0;
    for( int i = 0; i < nSize-1; i++) {
        double x1 = lassoPoints[i][0] - p0[0];
        double y1 = lassoPoints[i][1] - p0[1];
        double x2 = lassoPoints[i+1][0] - p0[0];
        double y2 = lassoPoints[i+1][1] - p0[1];
        angle    += JMath::Angle2D(x1,y1,x2,y2);
    }

    if(fabs(angle) < M_PI) return 0;
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: addWithinLasso()
{
    if( !closedLassoCheckBox->isChecked() )  return;

    size_t nSize = lassoPoints.size();
    if( nSize < 3) return;

    double xmin =  0.90*std::numeric_limits<double>::max();
    double ymin =  0.90*std::numeric_limits<double>::max();
    double xmax = -0.90*std::numeric_limits<double>::max();
    double ymax = -0.90*std::numeric_limits<double>::max();

    for( size_t i = 0; i < nSize; i++) {
        const Point3D &xyz = lassoPoints[i];
        xmin = std::min(xmin, xyz[0] );
        ymin = std::min(ymin, xyz[1] );
        xmax = std::max(xmax, xyz[0] );
        ymax = std::max(ymax, xyz[1] );
    }

    JNodeSequence fixednodes;
    nSize = mesh->getSize(0);
    for( size_t i = 0; i < nSize; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( !vtx->isActive() ) continue;
        const Point3D &xyz = vtx->getXYZCoords();
        if( xyz[0] < xmin || xyz[1] < ymin || xyz[0] > xmax || xyz[1] > ymax) continue;
        if( isWithinLasso(xyz)  ) {
            add2Group(vtx);
            fixednodes.push_back(vtx);
        }
    }
    groupID++;
    recordConstraints.push_back(fixednodes);

    lassoPoints.clear();
    constraintsViewer->clearLasso();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshConstraintsDialog :: addOnLasso()
{
    double xmin =  0.90*std::numeric_limits<double>::max();
    double ymin =  0.90*std::numeric_limits<double>::max();
    double xmax = -0.90*std::numeric_limits<double>::max();
    double ymax = -0.90*std::numeric_limits<double>::max();

    size_t nSize =  lassoPoints.size();
    for( size_t i = 0; i < nSize; i++) {
        const Point3D &xyz = lassoPoints[i];
        xmin = std::min(xmin, xyz[0] );
        ymin = std::min(ymin, xyz[1] );
        xmax = std::max(xmax, xyz[0] );
        ymax = std::max(ymax, xyz[1] );
    }

    JEdgeSequence feasible;
    size_t numEdges = mesh->getSize(1);
    for( size_t i = 0; i < numEdges; i++) {
        JEdgePtr edge = mesh->getEdgeAt(i);
        if( !edge->isActive() )  continue;
        const Point3D &p0 = edge->getNodeAt(0)->getXYZCoords();
        const Point3D &p1 = edge->getNodeAt(1)->getXYZCoords();
        if( min(p0[0], p1[0]) > xmax ) continue;
        if( min(p0[1], p1[1]) > ymax ) continue;
        if( max(p0[0], p1[0]) < xmin ) continue;
        if( max(p0[1], p1[1]) > ymin ) continue;
        feasible.push_back(edge);
    }
}
///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: eraseRegionPoints( const JNodePtr &vi)
{
    if( mesh == nullptr ) return;

    if( !vi->hasAttribute("Constraint") )  return;

    int id1 = -1, id2 = -1;
    vi->getAttribute("Constraint", id1);
    assert( id1 >= 0);

    size_t numnodes = mesh->getSize(0);
    for( size_t j = 0; j < numnodes; j++) {
        JNodePtr vj =  mesh->getNodeAt(j);
        if( vj->hasAttribute("Constraint") ) {
            vj->getAttribute("Constraint", id2);
            if( id1 == id2 ) {
                vj->deleteAttribute("Constraint");
                constraintsViewer->setConstraint(vj, 0);
            }
        }
    }
    constraintsViewer->deleteGroup(id1);

    countConstraints();
}

////////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: mousePressEvent(QMouseEvent *e)
{
    if( !this->isVisible() ) return;

    if( e->button() == Qt::LeftButton && ( e->modifiers() & Qt::ShiftModifier) )
        left_button_pressed = 1;

    assert( meshViewer );

    if( e->button() == Qt::RightButton) lassoPoints.clear();

    QDialog::mousePressEvent(e);
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    left_button_pressed = 0;

    if( !this->isVisible() || picker == nullptr ) return;

    JNodeSequence nodeSeq = picker->getPickedNodes();
    if( e->button() == Qt::LeftButton ) {
        if( addNodeRadioButton->isChecked() ) {
            if( lassoRadioButton->isChecked() ) {
                if( closedLassoCheckBox->isChecked() ) {
                    parameterizeLasso();
                    int nSize = lassoPoints.size();
                    int isign1 = getSign( lassoPoints[0], lassoPoints[nSize/4], lassoPoints[nSize/2] );
                    if( isign1 < 0)
                        std::reverse( lassoPoints.begin(), lassoPoints.end() ) ;
                    addWithinLasso();
                } else {
                    Point3D xyz;
                    int err = viewManager->getMouseCurrentXYZPosition(xyz);
                    if( !err) {
                        lassoPoints.push_back(xyz);
                        constraintsViewer->setLasso(lassoPoints);
                    }
                }
            }

            if( meshNodeRadioButton->isChecked() ) {
                picker->clearAll();
                for(size_t i = 0; i < nodeSeq.size(); i++) {
                    JNodePtr vtx = nodeSeq[i];
                    if( !vtx->hasAttribute("Constraint") ) {
                        add2Group(vtx);
                    }
                    groupID++;
                }
                if( !nodeSeq.empty() )
                    recordConstraints.push_back(nodeSeq);
            }
        }

        if( deleteNodeRadioButton->isChecked() ) {
            for(size_t i = 0; i < nodeSeq.size(); i++) {
                JNodePtr vtx = nodeSeq[i];
                if( vtx->hasAttribute("Constraint") ) {
                    eraseRegionPoints(vtx);
                }
            }
            picker->clearAll();
        }
    }
    countConstraints();

    if( e->button() == Qt::RightButton ) {
        JMeshEntityPickerPtr picker = meshViewer->getEntityPicker();
        if( picker ) picker->clearAll();
    }

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: mouseMoveEvent(QMouseEvent *e)
{
    if( !this->isVisible() || left_button_pressed == 0 ) return;

    Point3D xyz;
    if( e->type() == QEvent::MouseMove && left_button_pressed ) {
        if( lassoRadioButton->isChecked()) {
            int err = viewManager->getMouseCurrentXYZPosition(xyz);
            if( !err) {
                xyz[2] = 0.0;
                if( lassoPoints.empty() )
                    lassoPoints.push_back(xyz);
                else {
                    const Point3D &plast = lassoPoints.back();
                    double len = JMath::length(plast, xyz);
                    if( len > 0.0) lassoPoints.push_back(xyz);
                }
                constraintsViewer->setLasso(lassoPoints);
            }
        }
    }
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////


void JMeshConstraintsDialog :: clearAll()
{
    if( meshViewer == nullptr || mesh == nullptr ) return;

    JMeshEntityPickerPtr picker = meshViewer->getEntityPicker();
    if( picker ) picker->clearAll();

    lassoPoints.clear();

    mesh->deleteNodeAttribute("Constraint");
    assignColors();

    /*
        size_t nSize = mesh->getSize(0);
        for( size_t i = 0; i < nSize; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            constraintsViewer->setConstraint( vtx, 0);
        }
        constraintsViewer->clearAll();
    */
    countConstraints();

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: addBoundaryConstraints()
{
    if( mesh == nullptr ) return;

    mesh->getTopology()->searchBoundary();
    size_t numnodes = mesh->getSize(0);

    int gid;
    if( boundaryNodesCheckBox->isChecked() ) {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            int  err =  vtx->getAttribute("Constraint", gid);
        }
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isActive() && vtx->isBoundary() ) {
                vtx->setAttribute("Constraint", groupID);
            }
        }
        groupID++;
    } else {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isBoundary() ) {
                vtx->deleteAttribute("Constraint");
            }
        }
    }

    assignColors();
    if( meshViewer ) meshViewer->refreshDisplay();
}

//////////////////////////////////////////////////////////////////////////////////////
/*
void JMeshConstraintsDialog :: setGlyph()
{
    if( mesh == nullptr ) return;
    size_t numnodes = mesh->getSize(0);

    bool addConstraint = addConstraintRadioButton->isChecked();

    double  radius = meshViewer->getDrawNode()->getBallRadius();

    RenderNodeAttribute *attrib = nullptr;
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", attrib);
            attrib->glyph = DrawNode::NODE_AS_CIRCLE;
            attrib->ballRadius = radius;
        }
    }
}
*/
//////////////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: assignColors()
{
    if( mesh == nullptr ) return;

    JColor rColor = JEntityColor::getColor("Red");
    JColor gColor = JEntityColor::getColor("Green");

    JNodeRenderPtr attrib;
    numConstraintNodes = 0;

    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++)  {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", attrib);
            attrib->color = gColor;
            attrib->display = 1;
            if( vtx->hasAttribute("Constraint") ) {
                attrib->color = rColor;
                attrib->glyph = JNodeRender::NODE_AS_SPHERE;
                numConstraintNodes++;
            } else {
                attrib->glyph = JNodeRender::NODE_AS_POINT;
                attrib->display = 0;
            }
        }
    }
    meshViewer->updateBuffers(mesh);

    numConstraintNodesLineEdit->setText( QString::number(numConstraintNodes) );
    if( meshViewer ) meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: setSelectionMode()
{
    /*
        JMeshEntityPickerPtr picker = meshViewer->getEntityPicker();
        if( picker == nullptr ) return;

        picker->clearAll();

        if( meshNodeRadioButton->isChecked() ) {
            picker->setPickableEntity(0);
            picker->setMode( 2 );
        } else {
            picker->setPickableEntity(-1);
        }

        // If you created Lasso earlier, remove it ...
        constraintsViewer->clearLasso();
        viewManager->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void JMeshConstraintsDialog :: setSelectRegion()
{
}

////////////////////////////////////////////////////////////////////////////////
void JMeshConstraintsDialog :: openPolygonDialog()
{
    if(!polygonRadioButton->isChecked() ) {
        return;
    }

    if( polygonDialog == nullptr )
        polygonDialog.reset(new JPolygonDialog(this));

    viewManager->attach(polygonDialog.get());

    QObject::connect( polygonDialog.get(), SIGNAL(setPolygon()), this, SLOT(getPolygon()));

    polygonDialog->setViewManager( viewManager );
    polygonDialog->show();
    this->hide();
}
////////////////////////////////////////////////////////////////////////////////
void JMeshConstraintsDialog :: saveConstraints()
{
    if( mesh == nullptr ) return;

    QString qstr  = QFileDialog::getSaveFileName(this,
                    *new QString("Select Constraints File "),
                    lastSelectedDirectory,
                    *new QString( "File Format (*.dat)"));

    string filename = qstr.toUtf8().constData();

    ofstream ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) return;

    size_t numNodes = mesh->getSize(0);
    size_t nCount   = 0;
    for( size_t i = 0; i < numNodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() )
            if( vtx->hasAttribute("Constraint") ) nCount++;
    }

    Point3D xyz;
    int  gid;

    ofile << numNodes << " " << nCount << endl;
    for( size_t i = 0; i < numNodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->hasAttribute("Constraint") ) {
            vtx->getAttribute("Constraint",   gid);
            ofile << vtx->getID() << "  " << gid << endl;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: readConstraints()
{
    if( mesh == nullptr ) return;

    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Constraints File "),
                   lastSelectedDirectory,
                   *new QString( "File Format (*.dat)"));

    string filename = qstr.toUtf8().constData();

    ifstream ifile(filename.c_str(), ios::in);
    if( ifile.fail() ) return;

    Point3D xyz;
    int  vid, gid, numNodes, nCount;

    ifile >> numNodes >> nCount;
    for( int i = 0; i < nCount; i++) {
        ifile >> vid >> gid >> xyz[0] >> xyz[1] >> xyz[2];
        const JNodePtr &vtx = mesh->getNodeAt(vid);
        if( vtx ) {
            vtx->setAttribute("Constraint", gid);
            constraintsViewer->setConstraint(vtx, 1);
        }
    }
    viewManager->refreshDisplay();

}
//////////////////////////////////////////////////////////////////////////////////////
void JMeshConstraintsDialog :: lastDelete()
{
    if( recordConstraints.empty() ) return;

    JNodeSequence nodes = recordConstraints.back();
    for( size_t i = 0; i < nodes.size(); i++) {
        constraintsViewer->setConstraint(nodes[i], 0);
    }
    recordConstraints.pop_back();
    countConstraints();


    viewManager->refreshDisplay();
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: invertConstraints()
{
    if( mesh == nullptr) return;

    size_t numnodes = mesh->getSize(0);
    int val = 0;
    for( size_t i = 0; i < numnodes; i++)  {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->hasAttribute("Constraint" ) )
            vtx->deleteAttribute("Constraint");
        else
            vtx->setAttribute("Constraint", val);
    }
    assignColors();
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsDialog :: closeDialog()
{

#ifdef CSV
    //
    // When you close this dialog and retunr to the parent widget, then you
    // should be able to move the constaints and solve some equation. In the
    // parent widget you are allowed to move only one vertex set at a time.
    // The following will change the picking mode..

    /*
        if( meshViewer ) {
            JMeshEntityPickerPtr picker = meshViewer->getEntityPicker();
            if( picker  ) picker->setMode( 1 );
        }
    */


    // No lasso present, when the widget is closed...
    vector<Point3D> emptyLasso;
    constraintsViewer->setLasso(emptyLasso);

    // Detech the constraintViewer.
    viewManager->detach(constraintsViewer.get() );

    Color redColor;
    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 1.0;
    JNodeRenderPtr attrib;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", attrib);
            if( vtx->hasAttribute("Constraint") ) {
                attrib->display = 1;
                attrib->glyph   = JNodeRender::NODE_AS_SPLAT;
                attrib->color   = redColor;
            } else {
                attrib->display = 0;
                attrib->glyph   = JNodeRender::NODE_AS_POINT;
            }
        }
    }
    viewManager->refreshDisplay();
#endif

    emit setNewConstraints();
    parentWidget()->show();
    close();
}

/////////////////////////////////////////////////////////////////////////////////
void JMeshConstraintsDialog :: makeConnections()
{
    RadioButton( lassoRadioButton,  [=] {setSelectionMode(); });
    RadioButton( polygonRadioButton,  [=] {setSelectionMode(); });
    RadioButton( addNodeRadioButton,  [=] { setSelectionMode(); });
    RadioButton( deleteNodeRadioButton,  [=] {setSelectionMode(); });

    PushButton( invertPushButton,  [=] {invertConstraints(); });
    PushButton( polygonPushButton,   [=] { openPolygonDialog();});
    PushButton( clearAllPushButton,  [=] { clearAll();});
    PushButton( deletelastPushButton, [=] {lastDelete();});
    PushButton( saveConstraintsPushButton, [=] {saveConstraints(); });
    PushButton( readConstraintsPushButton, [=] {readConstraints(); });

    CheckBox( boundaryNodesCheckBox, [=] {addBoundaryConstraints(); });

    PushButton( closePushButton,  [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////

JMeshConstraintsViewer :: JMeshConstraintsViewer()
{
    mesh = nullptr;
    canvas = 0;

    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 0.5;

    greenColor[0] = 0.0;
    greenColor[1] = 1.0;
    greenColor[2] = 0.0;
    greenColor[3] = 0.5;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshConstraintsViewer :: draw()
{
    /*
        // From free hand drawing, the Z Componnent can be errorneous, therefore
        // use a canvas just below the mesh..
        if( canvas ) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glColor4f( 1.00, 1.00, 1.00, 0.0 );
            double l = 100;
            glBegin(GL_QUADS);
            glVertex3f( -l, -l , -0.0002);
            glVertex3f(  l, -l , -0.0002);
            glVertex3f(  l,  l , -0.0002);
            glVertex3f( -l,  l , -0.0002);
            glEnd();
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
        }
        glDisable(GL_LIGHTING);

        // Draw the lasso (including the Polygon).
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
        size_t numNodes = lassoPoints.size();
        if( numNodes > 1) {
            glLineWidth(2.0);
            glColor3f( 1.0, 0.0, 0.0);
            glBegin(GL_LINES);
            for( size_t i = 0; i < numNodes-1; i++) {
                const Point3D &p0  = lassoPoints[i];
                const Point3D &p1  = lassoPoints[i+1];
                glVertex3f( p0[0], p0[1], p0[2] + 0.001 );
                glVertex3f( p1[0], p1[1], p1[2] + 0.001 );
            }
            glEnd();
            glLineWidth(1.0);

            glPointSize(5.0);
            glColor3f( 0.0, 0.0, 1.0);
            glBegin(GL_POINTS);
            for( size_t i = 0; i < numNodes; i++) {
                const Point3D &p0  = lassoPoints[i];
                glVertex3f( p0[0], p0[1], p0[2] + 0.0001 );
            }
            glEnd();
            glPointSize(1.0);
        }
        glEnable( GL_LIGHTING );
    */
}

///////////////////////////////////////////////////////////////////////////////

int JMeshConstraintsViewer :: setConstraint( const JNodePtr &vtx, bool val)
{
    /*
        if( !vtx->isActive() ) return 1;

        JNodeRenderPtr attrib;
        int err = vtx->getAttribute("Render", attrib);
        if( err ) return 1;

        attrib->display = 1;
        attrib->color = greenColor;

        if( vtx->hasAttribute("Constraint") )
            attrib->color = redColor;
        else
            attrib->color = greenColor;
    */

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshConstraintsViewer :: setConstraints( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return 0;

    size_t numnodes = mesh->getSize(0);
    int nCount = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        int err = setConstraint(vtx);
        if( !err )  nCount++;
    }
    return nCount;
}

///////////////////////////////////////////////////////////////////////////////
