#include "CameraDialog.hpp"
#include "MeshViewer.hpp"

///////////////////////////////////////////////////////////////////////////////

JCameraDialog :: JCameraDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JCameraDialog :: ~JCameraDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JCameraDialog :: setViewManager( JaalViewer *v) {
    viewManager = v;
    viewManager->attach(this);
    setValues();
}
///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( viewManager ) viewManager->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: setProjection()
{
    int proj = 0;
    if( perspectiveRadioButton->isChecked() )  proj = 1;
    viewManager->setProjection( proj );
}

///////////////////////////////////////////////////////////////////////////////
void JCameraDialog ::  wheelEvent( QWheelEvent *e)
{
   setValues();
}
///////////////////////////////////////////////////////////////////////////////
void JCameraDialog :: setValues()
{
    if( viewManager == nullptr) return;

    float val;
    qglviewer::Vec vec;

    vec = viewManager->camera()->position();
    xcameraLineEdit->setText( QString::number(vec.x) );
    ycameraLineEdit->setText( QString::number(vec.y) );
    zcameraLineEdit->setText( QString::number(vec.z) );
    double dist =  sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
    distanceLineEdit->setText( QString::number(dist) );

    val = viewManager->camera()->sceneRadius();
    sceneRadiusLineEdit->setText( QString::number(val) );

    val = viewManager->camera()->fieldOfView();
    fovLineEdit->setText( QString::number(val) );

    val = viewManager->camera()->aspectRatio();
    aspectRatioLineEdit->setText( QString::number(val) );

    vec = viewManager->camera()->sceneCenter();
    xcenterLineEdit->setText( QString::number(vec.x) );
    ycenterLineEdit->setText( QString::number(vec.y) );
    zcenterLineEdit->setText( QString::number(vec.z) );

    val = viewManager->camera()->zNear();
    znearLineEdit->setText( QString::number(val) );
    val = viewManager->camera()->zFar();
    zfarLineEdit->setText( QString::number(val) );

    val = viewManager->camera()->screenWidth();
    screenWidthLineEdit->setText( QString::number(val) );

    val = viewManager->camera()->screenHeight();
    screenHeightLineEdit->setText( QString::number(val) );
}

///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: entireScene()
{
    setViewSide();
    setValues();
}

///////////////////////////////////////////////////////////////////////////////
void JCameraDialog :: setPosition()
{
    qglviewer::Vec   vec;
    QString qstr;
    qstr = xcameraLineEdit->text();
    vec[0] = qstr.toDouble();

    qstr = ycameraLineEdit->text();
    vec[1] = qstr.toDouble();

    qstr = zcameraLineEdit->text();
    vec[2] = qstr.toDouble();

    viewManager->camera()->setPosition(vec);
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: setRotationConstraint()
{
    QString qstr;
    qstr = rotationConstraintComboBox->currentText();
    string str = qstr.toUtf8().constData();

    if( str == "Free") viewManager->setRotationAxis(JaalViewer::FREE_ROTATION);
    if( str == "X-Axis") viewManager->setRotationAxis(JaalViewer::XAXIS_ROTATION);
    if( str == "Y-Axis") viewManager->setRotationAxis(JaalViewer::YAXIS_ROTATION);
    if( str == "Z-Axis") viewManager->setRotationAxis(JaalViewer::ZAXIS_ROTATION);
    if( str == "Forbidden") viewManager->setRotationAxis(JaalViewer::NO_ROTATION);
}

///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: setTranslationConstraint()
{
    QString qstr;
    qstr = translationConstraintComboBox->currentText();
    string str = qstr.toUtf8().constData();

    if( str == "Free")   viewManager->setTranslationAxis(JaalViewer::FREE_TRANSLATION);
    if( str == "X-Axis") viewManager->setTranslationAxis(JaalViewer::XAXIS_TRANSLATION);
    if( str == "Y-Axis") viewManager->setTranslationAxis(JaalViewer::YAXIS_TRANSLATION);
    if( str == "Z-Axis") viewManager->setTranslationAxis(JaalViewer::ZAXIS_TRANSLATION);
    if( str == "Forbidden") viewManager->setTranslationAxis(JaalViewer::NO_TRANSLATION);
}
///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: setSceneRadius()
{
    QString qstr = sceneRadiusLineEdit->text();
    viewManager->setSceneRadius( qstr.toDouble() );
    viewManager->showEntireScene();
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: setViewSide()
{
    if( viewManager == nullptr) return;

    QString qstr = viewSideComboBox->currentText();
    string str   = StdString(qstr);

    if( str == "Right")
        viewManager->resetView( JViewDirection:: RIGHT_VIEW);

    if( str == "Left")
        viewManager->resetView( JViewDirection:: LEFT_VIEW);

    if( str == "Top")
        viewManager->resetView( JViewDirection:: TOP_VIEW);

    if( str == "Bottom")
        viewManager->resetView( JViewDirection:: BOTTOM_VIEW);

    if( str == "Front")
        viewManager->resetView( JViewDirection:: FRONT_VIEW);

    if( str == "Back")
        viewManager->resetView( JViewDirection:: BACK_VIEW);
}
///////////////////////////////////////////////////////////////////////////////
void JCameraDialog :: setViewPoint()
{
    if( viewManager == nullptr) return;

    QString qstr = viewPointComboBox->currentText();
    string str   = StdString(qstr);

    if( str == "000")
        viewManager->resetView( JViewDirection:: CUBE_POINT_0_VIEW);

    if( str == "100")
        viewManager->resetView( JViewDirection:: CUBE_POINT_1_VIEW);

    if( str == "110")
        viewManager->resetView( JViewDirection:: CUBE_POINT_2_VIEW);

    if( str == "010")
        viewManager->resetView( JViewDirection:: CUBE_POINT_3_VIEW);

    if( str == "001")
        viewManager->resetView( JViewDirection:: CUBE_POINT_4_VIEW);

    if( str == "101")
        viewManager->resetView( JViewDirection:: CUBE_POINT_5_VIEW);

    if( str == "111")
        viewManager->resetView( JViewDirection:: CUBE_POINT_6_VIEW);

    if( str == "011")
        viewManager->resetView( JViewDirection:: CUBE_POINT_7_VIEW);
}

///////////////////////////////////////////////////////////////////////////////
void JCameraDialog :: setCamera()
{
    QString qstr;
    qstr = znearLineEdit->text();
    double znear = qstr.toDouble();

    double d  = viewManager->camera()->distanceToSceneCenter();
    double zcoeff = (d-znear)/viewManager->sceneRadius();
    viewManager->camera()->setZClippingCoefficient(zcoeff);
    viewManager->refreshDisplay();
    setValues();
}
///////////////////////////////////////////////////////////////////////////////
void JCameraDialog :: setWindowSize()
{
    QString qstr;
    qstr = screenWidthLineEdit->text();
    int w = qstr.toInt();

    qstr = screenHeightLineEdit->text();
    int h = qstr.toInt();

// viewManager->camera()->setScreenWidthAndAHeight(w,h);
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: mouseMoveEvent(QMouseEvent *e)
{
    if( e == nullptr) return;
    setValues();
}
///////////////////////////////////////////////////////////////////////////////
void JCameraDialog :: closeDialog()
{
    viewManager->detach(this);
    this->close();
}
///////////////////////////////////////////////////////////////////////////////

void JCameraDialog :: makeConnections()
{
    connect( rotationConstraintComboBox, SIGNAL( activated(int) ), this, SLOT( setRotationConstraint() ));
    connect( translationConstraintComboBox, SIGNAL( activated(int) ), this, SLOT( setTranslationConstraint() ));
    connect( orthoRadioButton,  SIGNAL( toggled( bool ) ) , this, SLOT( setProjection() ));
    connect( perspectiveRadioButton,  SIGNAL( toggled( bool ) ) , this, SLOT( setProjection() ));
    connect( xcameraLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setPosition() ));
    connect( ycameraLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setPosition() ));
    connect( zcameraLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setPosition() ));
    connect( distanceLineEdit, SIGNAL( editingFinished() ) , this, SLOT( setPosition() ));
    connect( screenWidthLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setWindowSize() ));
    connect( screenHeightLineEdit, SIGNAL( editingFinished() ) , this, SLOT( setWindowSize() ));

    connect( znearLineEdit,    SIGNAL( editingFinished() ) , this, SLOT( setCamera() ));

    connect( sceneRadiusLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setSceneRadius() ));

    connect( viewSideComboBox, SIGNAL( activated(int) ), this, SLOT( setViewSide() ));
    connect( viewPointComboBox, SIGNAL( activated(int) ), this, SLOT( setViewPoint() ));
    connect( entireScenePushButton,  SIGNAL( clicked() ), this, SLOT( entireScene() ));
//  connect( refreshPushButton,  SIGNAL( clicked() ), this, SLOT( setValues() ));
    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}

///////////////////////////////////////////////////////////////////////////////
