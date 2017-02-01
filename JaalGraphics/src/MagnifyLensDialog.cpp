#include "MagnifyLensDialog.hpp"
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
JMagnifyingLens:: JMagnifyingLens()
{
    width  = 0;
    height = 0;
    center[0] = 0;
    center[1] = 0;
    active    = 1;
    connection[0] = 0;
    connection[1] = 1;
    connection[2] = 2;
    connection[3] = 3;
    texID     = 0;
    shape     = 0;     // 0 is Circle, 1 is rectangle ...
    lensRadius  = 50;
    projRadius  = 200;
    borderWidth = 2;
    borderColor[0] = 1.0;
    borderColor[1] = 0.0;
    borderColor[2] = 0.0;
    borderColor[3] = 1.0;
    filled         = 0;
}

///////////////////////////////////////////////////////////////////////////////
void JMagnifyingLens :: draw()
{
    if( !active ) return;
    viewManager->camera()->setType(Camera::ORTHOGRAPHIC);

    if( texID ) {
        glClear(GL_DEPTH_BUFFER_BIT);
        glDisable(GL_BLEND);
        glBindTexture(GL_TEXTURE_2D, texID);
        double U = 1.0;
        double V = 1.0;
        int winHeight = viewManager->height();
        int x = projCenter[0] - 0.5*projWidth;
        int y = projCenter[1] - 0.5*projHeight;
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        glEnable(GL_TEXTURE_2D);
        viewManager->startScreenCoordinatesSystem(true);
        glBegin(GL_QUADS);
        glTexCoord2f(0.0, 0.0);
        glVertex2i(x, y);
        glTexCoord2f(U,   0.0);
        glVertex2i(x + projWidth, y);
        glTexCoord2f(U,   V  );
        glVertex2i(x+projWidth, y+projHeight);
        glTexCoord2f(0.0, V  );
        glVertex2i(x, y + projHeight);
        glEnd();
        viewManager->stopScreenCoordinatesSystem();
        glDisable(GL_TEXTURE_2D);
        glEnable(GL_BLEND);
    }

    drawROI();
    viewManager->camera()->setType(Camera::PERSPECTIVE);
}

///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLens :: drawROI()
{
    if( shape == 1)
        rectangleROI();
    else
        circleROI();
}
///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLens :: rectangleROI()
{
    int x, y;

    glPushMatrix();
    int winHeight = viewManager->height();
    int winWidth  = viewManager->width();

    viewManager->startScreenCoordinatesSystem();
    glDisable( GL_LIGHTING);


    //  Draw the Lens ....
    vector<Point2I> quad1(4);
    x = center[0] - 0.5*width;
    y = winHeight- center[1] - 0.5*height;

    quad1[0][0] = x;
    quad1[0][1] = y;

    quad1[1][0] = x;
    quad1[1][1] = y + height;

    quad1[2][0] = x + width;
    quad1[2][1] = y + height;

    quad1[3][0] = x + width;
    quad1[3][1] = y;

    glDisable(GL_BLEND);
    glLineWidth(borderWidth);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3fv( &borderColor[0] );
    glBegin(GL_QUADS);
    glVertex2iv(&quad1[0][0]);
    glVertex2iv(&quad1[1][0]);
    glVertex2iv(&quad1[2][0]);
    glVertex2iv(&quad1[3][0]);
    glEnd();

    if( filled) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_BLEND);
        glColor4f( 0.0, 0.0, 1.0, 0.2);
        glBegin(GL_QUADS);
        glVertex2iv(&quad1[0][0]);
        glVertex2iv(&quad1[1][0]);
        glVertex2iv(&quad1[2][0]);
        glVertex2iv(&quad1[3][0]);
        glEnd();
        glDisable(GL_BLEND);
    }

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    x = projCenter[0] - 0.5*projWidth;
    y = winHeight- projCenter[1] - 0.5*projHeight;

    vector<Point2I> quad2(4);
    quad2[0][0] = x;
    quad2[0][1] = y;

    quad2[1][0] = x;
    quad2[1][1] = y + projHeight;

    quad2[2][0] = x + projWidth;
    quad2[2][1] = y + projHeight;

    quad2[3][0] = x + projWidth;
    quad2[3][1] = y;

    glLineWidth(borderWidth);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f( 1.0, 0.2, 0.2);
    glBegin(GL_QUADS);
    glVertex2iv(&quad2[0][0]);
    glVertex2iv(&quad2[1][0]);
    glVertex2iv(&quad2[2][0]);
    glVertex2iv(&quad2[3][0]);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f( 1.0, 1.0, 1.0);
    glBegin(GL_QUADS);
    glVertex2iv(&quad2[0][0]);
    glVertex2iv(&quad2[1][0]);
    glVertex2iv(&quad2[2][0]);
    glVertex2iv(&quad2[3][0]);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    for( int ipass = 0; ipass < 2; ipass++) {
        for (x = 5; x < 10; x++) {
            for (y = 5; y < 10; y++) {
                glColor4f (0.8f,0.8f,0.8f, 1.0/(ipass+1));
                glBegin(GL_QUADS);
                glVertex2i(quad2[0][0]+x, quad2[0][1]+y);
                glVertex2i(quad2[1][0]+x, quad2[1][1]+y);
                glVertex2i(quad2[2][0]+x, quad2[2][1]+y);
                glVertex2i(quad2[3][0]+x, quad2[3][1]+y);
                glEnd();
            }
        }
    }
    glDisable(GL_BLEND);

    glColor3f( 1.0, 0.2, 0.2);
    glBegin(GL_LINES);
    for( int i = 0; i < 4; i++) {
        if(connection[i] >= 0) {
            glVertex2iv( &quad1[i][0] );
            int j = connection[i];
            glVertex2iv( &quad2[j][0] );
        }
    }
    glEnd();
    glEnable( GL_LIGHTING);
    viewManager->stopScreenCoordinatesSystem();
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
void JMagnifyingLens :: drawCircle( int xs, int ys)
{
    int N = 100;
    double x, y, dtheta = 2.0*M_PI/N;
    int winHeight  = viewManager->height();
    glColor3f( 1.0, 1.0, 1.0);

    glBegin(GL_TRIANGLE_FAN);
    for( int i = 0; i < N; i++) {
        x = xs + projCenter[0] + projRadius*cos(i*dtheta);
        y = ys + winHeight - projCenter[1];
        glVertex2d(x,y);
        x = xs + projCenter[0] + projRadius*cos(i*dtheta);
        y = ys + winHeight - (projCenter[1] + projRadius*sin(i*dtheta));
        glVertex2d(x,y);
        x = xs + projCenter[0] + projRadius*cos((i+1)*dtheta);
        y = ys + winHeight - (projCenter[1] + projRadius*sin((i+1)*dtheta));
        glVertex2d(x,y);
    }
    glEnd();
}

///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLens :: circleROI()
{
    int winHeight  = viewManager->height();

    glPushMatrix();
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    viewManager->startScreenCoordinatesSystem();
    glDisable( GL_LIGHTING);

    int N = 100;
    double x, y, dtheta = 2.0*M_PI/N;
    glLineWidth(borderWidth);

    if( filled ) {
        glColor4f( 0.0, 0.0, 0.5, 0.5);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glBegin(GL_TRIANGLE_FAN);
        for( int i = 0; i < N; i++) {
            glVertex2d(center[0], winHeight-center[1]);
            x = center[0] + lensRadius*cos(i*dtheta);
            y = winHeight - (center[1] + lensRadius*sin(i*dtheta));
            glVertex2d(x, y);
            x = center[0] + lensRadius*cos((i+1)*dtheta);
            y = winHeight - (center[1] + lensRadius*sin((i+1)*dtheta));
            glVertex2d(x, y);
        }
        glEnd();
    }

    glColor4fv( &borderColor[0] );
    // Draw the Lens ...
    glBegin(GL_LINE_LOOP);
    for( int i = 0; i < N; i++) {
        x = center[0] + lensRadius*cos(i*dtheta);
        y = winHeight - (center[1] + lensRadius*sin(i*dtheta));
        glVertex2d(x,y);
    }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glBegin(GL_LINE_LOOP);
    for( int i = 0; i < N; i++) {
        double  x = projCenter[0] + projRadius*cos(i*dtheta);
        double y = winHeight - (projCenter[1] + projRadius*sin(i*dtheta));
        glVertex2d(x,y);
    }
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f( 1.0, 1.0, 1.0);
    drawCircle(0,0);

    glEnable( GL_LIGHTING);
    viewManager->stopScreenCoordinatesSystem();
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLens :: captureRegion()
{
    qglviewer::Vec  prevCameraPos = viewManager->camera()->position();

    int winHeight = viewManager->height();
    int xc = center[0] - 0.5*width;
    int yc = winHeight - (center[1] + 0.5*height);

    int viewport[4];
    int scissor[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    glGetIntegerv(GL_SCISSOR_BOX,scissor);

    if( shape == 0) {
        projWidth  = 2.0*projRadius;
        projHeight = 2.0*projRadius;
    }

    int xp = projCenter[0] - 0.5*projWidth;
    int yp = projCenter[1] - 0.5*projHeight;
    glViewport(xp,yp, projWidth, projHeight);
    glScissor(xp,yp, projWidth, projHeight);

    xc = center[0] - 0.5*width;
    yc = winHeight - center[1] - 0.5*height;
    QRect qr( xc, yc, width, height);
    viewManager->camera()->fitScreenRegion(qr);
    glPushMatrix();
    viewManager->drawScene();
    glPopMatrix();

    viewManager->refreshDisplay();

    QImage qimg = viewManager->grabFrameBuffer();
    int xt = projCenter[0] - 0.5*projWidth;
    int yt = winHeight - (projCenter[1] + 0.5*projHeight);
    QImage qsmall = qimg.copy(xt, yt, projWidth, projHeight);
    if( shape == 0) {
        QRgb value = qRgb(255,255,255);
        double r2 = 0.25*projWidth*projWidth;
        for( int j = 0; j < projHeight; j++) {
            for( int i = 0; i < projWidth; i++) {
                double dx = i - 0.5*projWidth;
                double dy = j - 0.5*projHeight;
                if( dx*dx + dy*dy > r2) {
                    qsmall.setPixel(i,j, value);
                }
            }
        }
    }

    texID = viewManager->bindTexture(qsmall);
    qsmall.save( QString("frame.png") );
    JImage img;
    img.readFrom("frame.png");
    texID = img.genTexture();

    glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
    viewManager->camera()->setPosition( prevCameraPos);

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

JMagnifyingLensDialog :: JMagnifyingLensDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;

    xlensSizeSpinBox->setValue(50);
    ylensSizeSpinBox->setValue(50);

    xprojSizeSpinBox->setValue(200);
    yprojSizeSpinBox->setValue(200);
}

///////////////////////////////////////////////////////////////////////////////

JMagnifyingLensDialog :: ~JMagnifyingLensDialog()
{
    viewManager->freezeView(0);
}

///////////////////////////////////////////////////////////////////////////////
void JMagnifyingLensDialog :: mouseMoveEvent(QMouseEvent *e)
{
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLensDialog :: init()
{
    if( viewManager == nullptr ) return;
    viewManager->freezeView(1);

    int w = viewManager->width();
    int h = viewManager->height();

    xlensCenterSpinBox->setMaximum(w);
    ylensCenterSpinBox->setMaximum(h);

    xlensCenterSpinBox->setValue(w/2.0);
    ylensCenterSpinBox->setValue(h/2.0);
    xlensSizeSpinBox->setValue(50);
    ylensSizeSpinBox->setValue(50);
    lensRadiusSpinBox->setValue(25);

    xprojCenterSpinBox->setMaximum(w);
    yprojCenterSpinBox->setMaximum(h);
    xprojCenterSpinBox->setValue(110);
    yprojCenterSpinBox->setValue(110);
    projRadiusSpinBox->setValue(100);

    xprojSizeSpinBox->setMaximum(w);
    yprojSizeSpinBox->setMaximum(h);
    xprojSizeSpinBox->setValue(200);
    yprojSizeSpinBox->setValue(200);
}

///////////////////////////////////////////////////////////////////////////////
void JMagnifyingLensDialog :: addNewLens()
{
    currLens.reset(new JMagnifyingLens());
    currLens->viewManager   = viewManager;

    int shape = 0;
    if( squareRadioButton->isChecked() ) shape = 1;
    currLens->shape = shape;

    int x = xprojSizeSpinBox->value();
    int y = yprojSizeSpinBox->value();

    currLens->projWidth   = x;
    currLens->projHeight  = y;

    currLens->projCenter[0] = 0.5*x;
    currLens->projCenter[1] = 0.5*y;
    currLens->lensRadius    = min(0.*x, 0.5*y);
    xprojCenterSpinBox->setValue(0.5*x);
    yprojCenterSpinBox->setValue(0.5*y);

    int winWidth  = viewManager->width();
    int winHeight = viewManager->height();
    xlensCenterSpinBox->setValue(0.5*winWidth);
    ylensCenterSpinBox->setValue(0.5*winHeight);
    currLens->width     = 0.25*currLens->projWidth;
    currLens->height    = 0.25*currLens->projHeight;
    currLens->center[0]  = 0.5*winWidth;
    currLens->center[1]  = 0.5*winHeight;

    if( currLens->shape == 0) {
        currLens->projWidth  = 2.0*currLens->projRadius;
        currLens->projHeight = 2.0*currLens->projRadius;
        currLens->width      = 2.0*currLens->lensRadius;
        currLens->height     = 2.0*currLens->lensRadius;
    }

    lensList.push_back(currLens);
    currentLensSpinBox->setMaximum( lensList.size()-1);
    viewManager->addObject(currLens);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLensDialog :: captureRegion()
{
    if( currLens == nullptr) return;

    if( currLens->texID ) glDeleteTextures(1, &currLens->texID);

    for( size_t i = 0; i < lensList.size(); i++)
        lensList[i]->active = 0;

    viewManager->freezeView(0);
    currLens->captureRegion();
    viewManager->freezeView(1);

    for( size_t i = 0; i < lensList.size(); i++)
        lensList[i]->active = 1;

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMagnifyingLensDialog :: actionMouseEvent(int id)
{
}

///////////////////////////////////////////////////////////////////////////////
void JMagnifyingLensDialog :: updateCamera()
{
    if( locked ) return;
    if( currLens == nullptr) return;

    currLens->filled      = filledCheckBox->isChecked();
    currLens->borderWidth = borderWidthSpinBox->value();
    currLens->width       = xlensSizeSpinBox->value();
    currLens->height      = ylensSizeSpinBox->value();
    currLens->center[0]   = xlensCenterSpinBox->value();
    currLens->center[1]   = ylensCenterSpinBox->value();
    currLens->lensRadius  = lensRadiusSpinBox->value();

    currLens->projWidth     = xprojSizeSpinBox->value();
    currLens->projHeight    = yprojSizeSpinBox->value();
    currLens->projCenter[0] = xprojCenterSpinBox->value();
    currLens->projCenter[1] = yprojCenterSpinBox->value();

    currLens->projRadius    = projRadiusSpinBox->value();
    currLens->connection[0] = connect0SpinBox->value();
    currLens->connection[1] = connect1SpinBox->value();
    currLens->connection[2] = connect2SpinBox->value();
    currLens->connection[3] = connect3SpinBox->value();

    if( currLens->shape == 0) {
        currLens->projWidth  = 2.0*currLens->projRadius;
        currLens->projHeight = 2.0*currLens->projRadius;
        currLens->width      = 2.0*currLens->lensRadius;
        currLens->height     = 2.0*currLens->lensRadius;
    }

    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLensDialog :: loadCamera()
{
    if( currLens ) {
        currLens->borderWidth = 2;
    }

    int lid = currentLensSpinBox->value();
    if( lid < 0 || lid >= lensList.size() ) return;

    currLens = lensList[lid];

    if( currLens == nullptr) return;

    currLens->borderWidth = 3;
    locked = 1;
    xlensSizeSpinBox->setValue( currLens->width);
    ylensSizeSpinBox->setValue( currLens->height);
    xlensCenterSpinBox->setValue( currLens->center[0]);
    ylensCenterSpinBox->setValue( currLens->center[1]);
    lensRadiusSpinBox->setValue( currLens->lensRadius);

    xprojSizeSpinBox->setValue(currLens->projWidth);
    yprojSizeSpinBox->setValue(currLens->projHeight);
    xprojCenterSpinBox->setValue(currLens->projCenter[0]);
    yprojCenterSpinBox->setValue(currLens->projCenter[1]);
    projRadiusSpinBox->setValue( currLens->projRadius);

    currLens->connection[0] = connect0SpinBox->value();
    currLens->connection[1] = connect1SpinBox->value();
    currLens->connection[2] = connect2SpinBox->value();
    currLens->connection[3] = connect3SpinBox->value();
    locked = 0;

    viewManager->refreshDisplay();
}


///////////////////////////////////////////////////////////////////////////////
void JMagnifyingLensDialog :: deleteLens()
{
    int lid = currentLensSpinBox->value();
    if( lid >= lensList.size() ) return;

    currLens = lensList[lid];
    currLens->deleted = 1;
    boost::remove_erase(lensList, currLens);
    if( lensList.empty() )
        currLens = nullptr;
    else
        currLens = lensList.back();

    viewManager->refreshDisplay();
    currentLensSpinBox->setMaximum( lensList.size()-1);
}
///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLensDialog :: closeDialog()
{
    viewManager->removeAllLens();
    lensList.clear();
    currLens.reset();
    viewManager->freezeView(0);
    this->close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////

void JMagnifyingLensDialog :: makeConnections()
{
    SpinBoxi( xlensCenterSpinBox,  [=] { updateCamera();} );
    SpinBoxi( ylensCenterSpinBox,  [=] { updateCamera();} );
    SpinBoxi( xlensSizeSpinBox,    [=] { updateCamera();} );
    SpinBoxi( ylensSizeSpinBox,    [=] { updateCamera();} );
    SpinBoxi( lensRadiusSpinBox,   [=] { updateCamera();} );

    SpinBoxi( xprojCenterSpinBox,  [=] {updateCamera();});
    SpinBoxi( yprojCenterSpinBox,  [=] {updateCamera();});

    SpinBoxi( xprojSizeSpinBox,    [=] {updateCamera();});
    SpinBoxi( yprojSizeSpinBox,    [=] {updateCamera();});
    SpinBoxi( projRadiusSpinBox,   [=] {updateCamera();});

    SpinBoxi( connect0SpinBox,  [=] {updateCamera();});
    SpinBoxi( connect1SpinBox,  [=] {updateCamera();});
    SpinBoxi( connect2SpinBox,  [=] {updateCamera();});
    SpinBoxi( connect3SpinBox,  [=] {updateCamera();});

    SpinBoxi( currentLensSpinBox, [=] {loadCamera();});
    SpinBoxi( borderWidthSpinBox,  [=] {updateCamera();});
    RadioButton( squareRadioButton, [=] {updateCamera(); });
    RadioButton( circleRadioButton, [=] {updateCamera();});

    PushButton( newLensPushButton,  [=] {addNewLens();});
    PushButton( deleteLensPushButton, [=] {deleteLens(); });
    PushButton( captureRegionPushButton, [=] {captureRegion(); });
    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
