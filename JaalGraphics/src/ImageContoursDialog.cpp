#include "ImageContoursDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JImageContoursViewer :: JImageContoursViewer()
{
    edgeWidth  = 2.0;
    nodeRadius = 0.02;

    nodeColor[0] = 1.0;
    nodeColor[1] = 0.0;
    nodeColor[2] = 0.0;
    nodeColor[3] = 0.5;

    edgeColor[0] = 0.0;
    edgeColor[1] = 1.0;
    edgeColor[2] = 0.0;
    edgeColor[3] = 1.0;
    drawIDs = 0;
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursViewer :: draw( JContour *contour)
{
    if( contour == nullptr ) return;
    size_t numNodes = contour->points.size();
    if( numNodes < 1) return;

    Point3D p0, p1;

    int   numSlices = 32;
    if( numNodes > 0) {
        glColor4fv( &nodeColor[0] );
        int numSlices = 32;
        const Point3D &xyz  = contour->points.back();
        glPushMatrix();
        glTranslatef( xyz[0], xyz[1], xyz[2]+0.0001 );
        gluDisk( JNodeDraw::diskObj, 0.0, nodeRadius, numSlices, 1);
        glPopMatrix();
    }

    if( numNodes > 1 && edgeWidth > 0) {
        glLineWidth(edgeWidth);
        glColor4fv( &edgeColor[0] );
        glBegin(GL_LINES);
        for( size_t i = 0; i < numNodes-1; i++) {
            p0  = contour->points[i];
            p1  = contour->points[i+1];
            glVertex3f( p0[0], p0[1], p0[2] + 0.005 );
            glVertex3f( p1[0], p1[1], p1[2] + 0.005 );
        }
        if( contour->closed) {
            p0  = contour->points.front();
            p1  = contour->points.back();
            glVertex3f( p0[0], p0[1], p0[2] + 0.005 );
            glVertex3f( p1[0], p1[1], p1[2] + 0.005 );
        }
        glEnd();
        glLineWidth(1.0);
    }

    if( drawPoints ) {
        glColor4fv( &nodeColor[0] );
        int numSlices = 32;
        int numStacks = 32;
        for( size_t i = 0; i < numNodes-1; i++) {
            const Point3D &xyz  = contour->points[i];
            glPushMatrix();
            glTranslatef( xyz[0], xyz[1], xyz[2]);
            gluSphere( JNodeDraw::sphereObj, 0.25*nodeRadius, numSlices, numStacks);
            glPopMatrix();
        }
    }

    char number[128];
    if( drawIDs) {
        float fontScale = 0.001*FontsManager::Instance().getFontScale();
        FTFont *font     = FontsManager::Instance().getFont(nullptr);
        JColor fontColor = FontsManager::Instance().getColor();
        glColor4fv( &fontColor[0] );
        for( size_t i = 0; i < numNodes; i++) {
            const Point3D &xyz = contour->points[i];
            sprintf(number, "%ld", i);
            glPushMatrix();
            glTranslatef(xyz[0], xyz[1], xyz[2] );
            glScalef(fontScale, fontScale, fontScale);
            font->Render(number);
            glPopMatrix();
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////

void JImageContoursViewer :: draw()
{
    glDisable( GL_LIGHTING );
    glEnable( GL_DEPTH_TEST);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for( size_t i = 0; i < contours.size(); i++)
        draw(contours[i] );

    glEnable( GL_LIGHTING );
}

///////////////////////////////////////////////////////////////////////////////

JImageContoursDialog :: JImageContoursDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    contourViewer.reset( new JImageContoursViewer() );
    contourViewer->setName("ImageContourViewer");
    currContour = nullptr;

    edgeWidthLineEdit->setText( QString::number(2) );
    nodeRadiusLineEdit->setText( QString::number(0.005) );
}

///////////////////////////////////////////////////////////////////////////////

JImageContoursDialog :: ~JImageContoursDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: init()
{
    currContour = nullptr;

    if( viewManager == nullptr ) return;

    imageViewer = nullptr;
    JViewComponentPtr c = viewManager->getComponent("ImageViewer");
    if(c) {
        imageViewer = dynamic_pointer_cast<JImageViewer>(c);
    }

    // You need to capture mouse for defining the constraints, so register this
    // Object to the main class, so that events can be known to this class...
    viewManager->attach( this );
    viewManager->setMouseTracking(1);
    viewManager->setRotationAxis(JaalViewer::NO_ROTATION);

    contourViewer->setViewManager(viewManager);

    viewManager->attach( contourViewer );
    viewManager->camera()->setType(Camera::ORTHOGRAPHIC);
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: keyPressEvent( QKeyEvent *e)
{
    if( currContour == nullptr ) return;
    if( e->key() == Qt::Key_Return ) return;

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: mousePressEvent(QMouseEvent *e)
{
    leftButton = 0;
    if( e->button() == Qt::LeftButton)  leftButton = 1;
}

///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: mouseMoveEvent(QMouseEvent *e)
{
    if( !freehandCheckBox->isChecked() ) return;
    if( leftButton == 0 || currContour == nullptr ) return;
    if( !this->isVisible() || viewManager == nullptr ) return;

    Point3D xyz;
    int   err  = viewManager->getMouseCurrentXYZPosition( xyz);
    if( !err) {
        currContour->points.push_back(xyz);
        viewManager->refreshDisplay();
    }
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( leftButton == 0 || currContour == nullptr ) return;
    if( !this->isVisible() || viewManager == nullptr ) return;

    bool val = showAllCheckBox->isChecked();
    contourViewer->setDrawPoints(val);

    Point3D xyz;
    int err  = viewManager->getMouseCurrentXYZPosition( xyz);
    if( err) return;

    currContour->points.push_back(xyz);
    viewManager->refreshDisplay();
    leftButton = 0;
    numNodesLineEdit->setText( QString::number(currContour->points.size()));
}

///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: deleteAll()
{
    for( size_t i = 0; i < contours.size(); i++)
        contours[i]->points.clear();
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: deleteSegment()
{
    if( currContour == nullptr ) return;
    if( currContour->points.empty() ) return;

    currContour->points.pop_back();
    numNodesLineEdit->setText( QString::number(currContour->points.size()));
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: newContour()
{
    currContour = new JContour();
    contours.push_back( currContour);
    contourViewer->append(currContour);
    numNodesLineEdit->setText( QString::number(0));
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: closeContour()
{
    if( currContour == nullptr ) return;
    currContour->closed = 1;
    currContour = nullptr;
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: readFile()
{
    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Contour File "));
    string filename = qstr.toUtf8().constData();
    ifstream infile( filename.c_str(), ios::in);
    if( infile.fail() ) return;

    Point3D xyz;
    string str;
    int numContours, numPoints;

    infile >> str;
    if( str != "<NumContours>") return;
    infile >> numContours;
    infile >> str;
    if(str != "</NumContours>") return;

    for( int i = 0; i < numContours; i++) {
        infile >> str;
        if( str == "<Contour>") {
            newContour();
            infile >> str;
            assert( str == "<Closed>");
            infile >> currContour->closed;
            infile >> str;
            assert( str == "</Closed>");
            infile >> str;
            assert( str == "<NumPoints>");
            infile >>numPoints;
            infile >> str;
            assert( str == "</NumPoints>");
            for( int j = 0; j < numPoints; j++) {
                infile >> xyz[0] >> xyz[1] >> xyz[2];
                currContour->points.push_back(xyz);
            }
            infile >> str;
            assert( str == "</Contour>");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: saveOFF( string &filename)
{
    ofstream  ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) return;

    int numnodes = 0;
    int numedges = 0;
    for( size_t i = 0; i < contours.size(); i++) {
        numnodes += contours[i]->points.size();
        numedges += contours[i]->points.size() - 1;
        if( contours[i]->closed ) numedges++;
    }
    ofile << "OFF" << endl;
    ofile << numnodes << "  0 " <<  numedges <<  endl;
    double z = 0.0;
    for( size_t i = 0; i < contours.size(); i++) {
        numnodes = contours[i]->points.size();
        for( int j = 0; j < numnodes; j++) {
            const Point3D &xyz = contours[i]->points[j];
            ofile << xyz[0] << " " << xyz[1] <<  "  0.0 " << endl;
        }
    }

    int numsegments = 0;
    for( size_t i = 0; i < contours.size(); i++) {
        numnodes = contours[i]->points.size();
        if( numnodes ) {
            if( contours[i]->closed )
                numsegments += numnodes;
            else
                numsegments += numnodes-1;
        }
    }

    int index = 0;
    int start_index = 0;
    int segid = 0;
    for( size_t i = 0; i < contours.size(); i++) {
        numnodes = contours[i]->points.size();
        start_index = index;
        for( int j = 0; j < numnodes-1; j++) {
            ofile << "1 "  << " " << index << " " << index+1 << endl;
            index = index + 1;
        }
        if( contours[i]->closed ) {
            ofile << "1  " << index << " " << start_index  << endl;
            index++;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: saveXML( string &filename)
{
    int nCount = 0;
    for( size_t i = 0; i < contours.size(); i++)
        if( !contours[i]->points.empty()  ) nCount++;

    ofstream  ofile(filename.c_str(), ios::out);
    if( ofile.fail() ) return;

    ofile << "<NumContours>  " << nCount << " </NumContours>" << endl;
    double z = 0.0;
    for( size_t i = 0; i < contours.size(); i++) {
        if( !contours[i]->points.empty()  ) {
            ofile << "<Contour> " << endl;
            ofile << "<Closed> " << contours[i]->closed << " </Closed>" << endl;
            ofile << "<NumPoints> " << contours[i]->points.size() << " </NumPoints>" << endl;
            for( size_t j = 0; j < contours[i]->points.size(); j++) {
                const Point3D &xyz = contours[i]->points[j];
                ofile << xyz[0] << " " << xyz[1] << " " << z << endl;
            }
            ofile << "</Contour> " << endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: saveAs()
{
    static QString lastSelectedDirectory;
    QString qstr  = QFileDialog::getSaveFileName(this,
                    *new QString("Select File "),
                    lastSelectedDirectory,
                    *new QString( "Mesh Format (*.off *.xml)"));

    string filename = qstr.toUtf8().constData();

    JWaitCursor waitCursor;
    waitCursor.start();

    if (!filename.empty()) {
        if (filename.rfind(".xml")  != string::npos) saveXML(filename);
        if (filename.rfind(".off") != string::npos)  saveOFF(filename);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: mergeClosePoints()
{
    if( currContour == nullptr ) return;
    size_t numPoints = currContour->points.size();
    if( numPoints < 2) return;

    double sumlength = 0.0;
    for( size_t i = 0; i < numPoints-1; i++) {
        const Point3D &p0  = currContour->points[i];
        const Point3D &p1  = currContour->points[i+1];
        sumlength += JMath::length(p0, p1);
    }
    double totallength = sumlength;

    vector<Point3D> newpoints;
    int index = 0, jndex;
    while(1) {
        if( index >=  numPoints) break;
        const Point3D &p0  = currContour->points[index];
        newpoints.push_back(p0);
        jndex = index+1;
        for( size_t i = index+1; i < numPoints; i++) {
            const Point3D &p1  = currContour->points[i];
            if( JMath::length(p0, p1) > 1.0E-06*totallength) {
                jndex = i;
                break;
            }
        }
        index = jndex;
    }
    currContour->points = newpoints;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef CVD
void JImageContoursDialog :: uniformContour()
{
    if( currContour == nullptr ) return;
    size_t numPoints = currContour->points.size();
    if( numPoints < 3) return;

    mergeClosePoints();

    double sumlength = 0.0;
    for( size_t i = 0; i < numPoints-1; i++) {
        const Point3D &p0  = currContour->points[i];
        const Point3D &p1  = currContour->points[i+1];
        sumlength += JMath::length(p0, p1);
    }
    double totallength = sumlength;

    vector<double> param(numPoints);
    sumlength = 0.0;
    param[0] = 0.0;
    for( size_t i = 0; i < numPoints-1; i++) {
        const Point3D &p0  = currContour->points[i];
        const Point3D &p1  = currContour->points[i+1];
        sumlength += JMath::length(p0, p1);
        double p   = sumlength/totallength;
        if( p < 0.0) p = 0.0;
        if( p > 1.0) p = 1.0;
        param[i+1] =  p;
    }

    double dt = 1.0/(double)(numPoints-1);

    Point3D xyz;
    for( size_t i = 1; i < numPoints-1; i++) {
        double t = i*dt;
        for( int j = i-1; j < numPoints-1; j++) {
            if( t >= param[j] && t <= param[j+1] ) {
                double t0 = param[j];
                double t1 = param[j+1];
                if(t1-t0 > 0.0) {
                    double a  = (t-t0)/(t1-t0);
                    const Point3D &p0 = currContour->points[j];
                    const Point3D &p1 = currContour->points[j+1];
                    xyz[0] = (1-a)*p0[0] + a*p1[0];
                    xyz[1] = (1-a)*p0[1] + a*p1[1];
                    xyz[2] = (1-a)*p0[2] + a*p1[2];
                    newpoints[i] = xyz;
                }
                break;
            }
        }
    }
    for( size_t i = 0; i < numPoints; i++)
        currContour->points[i] = newpoints[i];

    sumlength = 0.0;
    for( size_t i = 0; i < numPoints-1; i++) {
        const Point3D &p0  = currContour->points[i];
        const Point3D &p1  = currContour->points[i+1];
        sumlength += JMath::length(p0, p1);
    }
    totallength = sumlength;

    sumlength = 0.0;
    for( size_t i = 1; i < numPoints-1; i++) {
        const Point3D &p0  = currContour->points[i];
        const Point3D &p1  = currContour->points[i+1];
        sumlength += JMath::length(p0, p1);
        double p   = sumlength/totallength;
        if( fabs(p - i*dt) > 1.0E-06)  {
            cout << "Param " << p << " Expected " << i*dt << endl;
            /*
                        QMessageBox msg;
                        msg.setIcon(QMessageBox::Warning);
                        msg.setText(tr("Some numerical error in unformization of the contour"));
                        msg.setStandardButtons( QMessageBox::Ok);
                        msg.exec();
            */
        }
    }

    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: smoothContour()
{
    if( currContour == nullptr ) return;
    size_t numPoints = currContour->points.size();
    if( numPoints < 3) return;

    vector<Point3D> newpoints(numPoints);

    for( size_t i = 0; i < numPoints; i++)
        newpoints[i] = currContour->points[i];

    Point3D xyz;
    for( size_t i = 1; i < numPoints-1; i++) {
        const Point3D &p0  = currContour->points[i];
        const Point3D &p1  = currContour->points[(i+2)%numPoints];
        xyz[0] = 0.5*(p0[0] + p1[0] );
        xyz[1] = 0.5*(p0[1] + p1[1] );
        xyz[2] = 0.5*(p0[2] + p1[2] );
        newpoints[(i+1)%numPoints] = xyz;
    }

    for( size_t i = 0; i < numPoints; i++)
        currContour->points[i] = newpoints[i];
    viewManager->refreshDisplay();
}
#endif

///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: setNodeRadius()
{
    QString str = nodeRadiusLineEdit->text() ;
    double radius  = str.toDouble();
    contourViewer->setNodeRadius(radius);
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: setEdgeWidth()
{
    QString str = edgeWidthLineEdit->text() ;
    double width  = str.toDouble();
    contourViewer->setEdgeWidth(width);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: setNodeColor()
{
    QColor color = QColorDialog::getColor();
    JColor rgb;
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    contourViewer->setNodeColor(rgb);
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: setEdgeColor()
{
    QColor color = QColorDialog::getColor();
    JColor rgb;
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    contourViewer->setEdgeColor(rgb);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: displayIDs()
{
    bool val = enumCheckBox->isChecked();
    contourViewer->setDrawID(val);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JImageContoursDialog :: genEdgeMesh()
{
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    JMeshViewerPtr meshViewer;
    if(c) meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

    JMeshPtr edgemesh = JMesh::newObject();
    int numContours = contours.size();

    int numPoints, numEdges;
    for( int ic = 0; ic < numContours; ic++) {
        numPoints = contours[ic]->points.size();
        if( numPoints < 1) continue;
        JNodeSequence nodes = JNode::newObjects( numPoints);
        for( size_t j = 0; j < numPoints; j++) {
            const Point3D &xyz  = contours[ic]->points[j];
            nodes[j]->setXYZCoords(xyz);
        }
        edgemesh->addObjects(nodes);
        if( contours[ic]->closed)
            numEdges = numPoints;
        else
            numEdges = numPoints-1;

        JEdgeSequence edges = JEdge::newObjects(numEdges);
        for( int j = 0; j < numPoints-1; j++)
            edges[j]->setNodes( nodes[j], nodes[j+1] );

        if( contours[ic]->closed)
            edges[numEdges-1]->setNodes( nodes.front(),  nodes.back() );
        edgemesh->addObjects(edges);
    }
    meshViewer->addObject(edgemesh);
    contourViewer->clearAll();
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: closeDialog()
{
    if( viewManager ) {
        viewManager->setMouseTracking(0);
        viewManager->setRotationAxis(JaalViewer::FREE_ROTATION);
    }
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JImageContoursDialog :: makeConnections()
{
    connect( nodeColorPushButton, &QPushButton::clicked, this, &JImageContoursDialog::setNodeColor);
    connect( edgeColorPushButton, &QPushButton::clicked, this, &JImageContoursDialog::setEdgeColor);
    connect( deleteAllPushButton, &QPushButton::clicked, this, &JImageContoursDialog::deleteAll );

    connect( deleteLastPushButton, SIGNAL( clicked() ), this, SLOT( deleteSegment() ));
    connect( newContourPushButton, SIGNAL( clicked() ), this, SLOT( newContour() ));
//  connect( uniformPushButton, SIGNAL( clicked() ), this, SLOT( uniformContour() ));
    connect( closeContourPushButton, SIGNAL( clicked() ), this, SLOT( closeContour() ));
    connect( savePushButton, SIGNAL( clicked() ), this, SLOT( saveAs() ));
    connect( nodeRadiusLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setNodeRadius() ));
    connect( edgeWidthLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setEdgeWidth() ));
    connect( enumCheckBox,  SIGNAL( toggled( bool ) ) ,  this, SLOT( displayIDs() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
