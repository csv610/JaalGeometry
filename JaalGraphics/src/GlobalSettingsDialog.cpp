#include "GlobalSettingsDialog.hpp"
#include "gl2ps.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JGlobalSettingsDialog :: JGlobalSettingsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    viewManager = nullptr;
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

JGlobalSettingsDialog :: ~JGlobalSettingsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: init()
{
    if( viewManager == nullptr) return;
    int w = viewManager->width();
    int h = viewManager->height();

    /*
       winWidthLineEdit->setText( QString::number(w) );
       winHeightLineEdit->setText( QString::number(h) );
    */
}

///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( viewManager == nullptr ) return;

    if( e->key() == Qt::Key_Return ) {
        viewManager->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////
void JGlobalSettingsDialog :: openCameraDialog()
{
    if( cameraDialog == nullptr) {
        cameraDialog.reset( new JCameraDialog(this));
    }
    cameraDialog->setViewManager( viewManager );

    cameraDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JGlobalSettingsDialog :: openFontsDialog()
{
    if( fontsDialog == nullptr) {
        fontsDialog.reset( new JFontsDialog(this));
        fontsDialog->setViewManager( viewManager );
    }

    fontsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: openLightsDialog()
{
    if( lightsDialog == nullptr) {
        lightsDialog.reset( new JLightsDialog(this));
    }
    lightsDialog->setViewManager( viewManager );
    lightsDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JGlobalSettingsDialog :: openFloorDialog()
{
    if( floorDialog == nullptr)
        floorDialog.reset( new JSceneFloorDialog(this));

    floorDialog->setViewManager( viewManager );
    floorDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: closeDialog()
{
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: setBackgroundColor()
{
    QColor color = QColorDialog::getColor();
    float rgb[3];
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    viewManager->setBackgroundColor( rgb );
}
///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: checkAxis()
{
    bool axis = axisCheckBox->isChecked();
    viewManager->setAxis(axis);
}
///////////////////////////////////////////////////////////////////////////////
void JGlobalSettingsDialog :: checkBoundingBox()
{
    bool axis = boundingBoxCheckBox->isChecked();
    viewManager->setBoundingBox(axis);
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: openScreenShotDialog()
{
    if( screenShotDialog == nullptr)
        screenShotDialog.reset( new JScreenShotDialog(this));

    screenShotDialog->setViewManager( viewManager );
    screenShotDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JGlobalSettingsDialog :: makeConnections()
{
    connect( screenShotPushButton, &QPushButton::clicked, [=] {openScreenShotDialog();} );
    connect( lightsPushButton,     &QPushButton::clicked, [=] {openLightsDialog();});
    connect( fontsPushButton,      &QPushButton::clicked, [=] {openFontsDialog();});
    connect( floorPushButton,      &QPushButton::clicked, [=] {openFloorDialog();});
    connect( cameraPushButton,     &QPushButton::clicked, [=] {openCameraDialog();});
    connect( backgroundColorPushButton, &QPushButton::clicked, [=] {setBackgroundColor();});

    connect( boundingBoxCheckBox,    SIGNAL( toggled( bool ) ) , this, SLOT( checkBoundingBox() ));
    connect( axisCheckBox,    SIGNAL( toggled( bool ) ) , this, SLOT( checkAxis() ));

    connect( closePushButton, &QPushButton::clicked,  [=] {closeDialog();});
}
///////////////////////////////////////////////////////////////////////////////
