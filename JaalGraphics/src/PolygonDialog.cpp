#include "PolygonDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JPolygonViewer :: JPolygonViewer()
{
    color[0]  = 0.2;
    color[1]  = 0.2;
    color[2]  = 0.2;
    color[3]  = 1.0;
    lineWidth = 1.0;
    numSides  = 4;
    canvas    = 0;
}

///////////////////////////////////////////////////////////////////////////////
void JPolygonViewer :: genPoints()
{
    double dx = endPos[0] - startPos[0];
    double dy = endPos[1] - startPos[1];
    double r  = sqrt(dx*dx + dy*dy );

    double dtheta = 2.0*M_PI/(double)numSides;
    double theta0 = atan2(dy, dx);

    polyPoints.resize(numSides);

    Point3D p0 = endPos;
    Point3D p1 = endPos;
    p0[2] = 0.0;
    p1[2] = 0.0;

    double theta = theta0;
    for( int i = 0; i < numSides; i++) {
        p1[0] = startPos[0] + r*cos(theta);
        p1[1] = startPos[1] + r*sin(theta);
        p0 = p1;
        theta += dtheta;
        polyPoints[i] = p1;
    }
}
///////////////////////////////////////////////////////////////////////////////

void JPolygonViewer :: drawPolygon()
{
    genPoints();

    glDisable(GL_LIGHTING );
    glDisable(GL_BLEND );

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINES);
    glLineWidth(2);
    glColor3f( 1.0, 0.1, 0.1);

    glBegin(GL_LINES);
    for( int i = 0; i < numSides; i++) {
        int j = (i+1)%numSides;
        glVertex3f( polyPoints[i][0], polyPoints[i][1], 0.002);
        glVertex3f( polyPoints[j][0], polyPoints[j][1], 0.002);
    }
    glEnd();

    glPointSize(5);
    glColor3f( 0.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    glVertex3f( startPos[0], startPos[1], 0.002);
    glEnd();

    glColor3f( 0.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    glVertex3f( endPos[0], endPos[1], 0.002);
    glEnd();
}

///////////////////////////////////////////////////////////////////////////////

void JPolygonViewer :: draw()
{
    if( canvas ) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3f( 0.8, 0.8, 0.8);
        double l = 10;
        glBegin(GL_QUADS);
        glVertex3f( -l, -l, -0.001 );
        glVertex3f(  l, -l, -0.001 );
        glVertex3f(  l,  l, -0.001 );
        glVertex3f( -l,  l, -0.001 );
        glEnd();
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }

    if( started ) drawPolygon();
}

///////////////////////////////////////////////////////////////////////////////

JPolygonDialog :: JPolygonDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();
    viewManager = nullptr;
    left_button_pressed = 0;

    numSidesLineEdit->setText( QString::number(4) );
    lengthLineEdit->setText( QString::number(1.0) );
}

///////////////////////////////////////////////////////////////////////////////

JPolygonDialog :: ~JPolygonDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JPolygonDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( viewManager ) viewManager->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JPolygonDialog :: setViewManager( JaalViewer *v)
{
    viewManager = v;
    polyViewer.reset( new JPolygonViewer );
    polyViewer->setName("PolygonViewer");
    viewManager->attach(polyViewer);

    viewManager->setMouseTracking(1);
    viewManager->setRotationAxis(JaalViewer::NO_ROTATION);
}
///////////////////////////////////////////////////////////////////////////////

void JPolygonDialog :: genShape()
{
    QString qstr;
    if( interactiveRadioButton->isChecked() ) {
        polyPoints =  polyViewer->getPoints();
    }

    else {
        qstr = numSidesLineEdit->text();
        int n = qstr.toInt();

        qstr = lengthLineEdit->text();
        double l = qstr.toDouble();

        JFacePtr newface;
        switch(n) {
        case 3:
            newface = JTriangle::getCanonical(l);
            break;
        case 4:
            newface = JQuadrilateral::getCanonical(l);
            break;
        default:
            newface = JPolygon::getCanonical(n,l);
            break;
        }
        JNodeSequence nodes = newface->getNodes();
        polyPoints.resize(n);
        for( int i = 0; i < n; i++) {
            JNodePtr v = newface->getNodeAt(i);
            polyPoints[i] = v->getXYZCoords();
        }
    }

    emit setPolygon();
}

///////////////////////////////////////////////////////////////////////////////

void JPolygonDialog :: mousePressEvent(QMouseEvent *e)
{
    left_button_pressed = 0;
    if( !this->isVisible() || polyViewer == nullptr ) return;
    if( canonicalRadioButton->isChecked() ) return;

    if( e->button() == Qt::LeftButton && ( e->modifiers() & Qt::ShiftModifier) )
        left_button_pressed = 1;

    /*
        if( left_button_pressed ) {
            const Point3D &vpos  = viewManager->getMouseXYZPosition();
            polyViewer->setStartPos(vpos);
            viewManager->refreshDisplay();
        }
    */

    QDialog::mousePressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JPolygonDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( !this->isVisible() || polyViewer == nullptr ) return;
    if( canonicalRadioButton->isChecked() ) return;

    /*
        if( left_button_pressed ) {
            const Point3D &vpos  = viewManager->getMouseXYZPosition();
            polyViewer->setEndPos(vpos);
            viewManager->refreshDisplay();
        }
    */

    QDialog::mouseReleaseEvent(e);
    left_button_pressed = 0;
}

///////////////////////////////////////////////////////////////////////////////
void JPolygonDialog :: mouseMoveEvent(QMouseEvent *e)
{
    if( !this->isVisible() || polyViewer == nullptr || left_button_pressed == 0) return;
    if( canonicalRadioButton->isChecked() ) return;

    /*
        const Point3D &vpos  = viewManager->getMouseXYZPosition();
        polyViewer->setEndPos(vpos);
    */

    viewManager->refreshDisplay();
}

////////////////////////////////////////////////////////////////////////////////////

void JPolygonDialog :: setNumSides()
{
    QString qstr = numSidesLineEdit->text();
    int n = qstr.toInt();
    if( polyViewer ) polyViewer->setNumSides(n);
}
////////////////////////////////////////////////////////////////////////////////////
void JPolygonDialog :: closeDialog()
{
    viewManager->detach(this);
    viewManager->detach(polyViewer);
    viewManager->setRotationAxis(JaalViewer::FREE_ROTATION);
    viewManager->setMouseTracking(1);
    this->parentWidget()->show();
    this->close();
}

////////////////////////////////////////////////////////////////////////////////////
void JPolygonDialog :: setGenMode()
{
    if( interactiveRadioButton->isChecked() ) {
        polyViewer->setCanvas(1);
    } else {
        viewManager->setMouseTracking(0);
        polyViewer->setCanvas(0);
    }
    viewManager->setMouseTracking(1);
}

////////////////////////////////////////////////////////////////////////////////////

void JPolygonDialog :: makeConnections()
{
    RadioButton( interactiveRadioButton, [=] { setGenMode();});
    RadioButton( canonicalRadioButton,   [=] { setGenMode();});

    LineEdit( numSidesLineEdit,  [=] { setNumSides();});

    PushButton( applyPushButton, [=] { genShape();});
    PushButton( closePushButton, [=] { closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
