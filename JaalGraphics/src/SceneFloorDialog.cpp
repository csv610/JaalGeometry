#include "SceneFloorDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JSceneFloorDialog :: JSceneFloorDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;

    distanceLineEdit->setText( QString::number(0.0) );
    lengthLineEdit->setText( QString::number(5.0) );
    pattern = 1;
}

///////////////////////////////////////////////////////////////////////////////

JSceneFloorDialog :: ~JSceneFloorDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JSceneFloorDialog :: setViewManager( JaalViewer *v)
{
    viewManager = v;
    if( v == nullptr) return;
    sceneFloor = viewManager->getSceneFloor();
    checkDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JSceneFloorDialog :: setColor()
{
    JColor rgb;
    QColor color = QColorDialog::getColor();
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    rgb[3] = 1.0;
    if( sceneFloor) {
        sceneFloor->setColor(rgb);
        viewManager->refreshDisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////
void JSceneFloorDialog :: setPattern()
{
    if( viewManager == nullptr ) return;
}

///////////////////////////////////////////////////////////////////////////////
void JSceneFloorDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        viewManager->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JSceneFloorDialog :: checkDisplay()
{
    if( sceneFloor == nullptr) return;

    bool val = displayCheckBox->isChecked();
    sceneFloor->setStatus(val);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JSceneFloorDialog :: setLength()
{
    if( sceneFloor == nullptr ) return;

    QString qstr = lengthLineEdit->text();
    double d = qstr.toDouble();

    sceneFloor->setLength(d);
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JSceneFloorDialog :: setLines()
{
    if( sceneFloor == nullptr ) return;

    int d = numGridsSpinBox->value();

    sceneFloor->setNumLines(d);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JSceneFloorDialog :: setDistance()
{
    if( sceneFloor == nullptr ) return;

    QString qstr = distanceLineEdit->text();
    double d = qstr.toDouble();

    sceneFloor->setDistance(d);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JSceneFloorDialog :: setDirection()
{
    if( sceneFloor == nullptr ) return;
    if( xfloorRadioButton->isChecked()) sceneFloor->setStdPlane(0);
    if( yfloorRadioButton->isChecked()) sceneFloor->setStdPlane(1);
    if( zfloorRadioButton->isChecked()) sceneFloor->setStdPlane(2);

    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JSceneFloorDialog :: makeConnections()
{
    RadioButton( xfloorRadioButton, [=] { setDirection();});
    RadioButton( yfloorRadioButton, [=] { setDirection();});
    RadioButton( zfloorRadioButton, [=] { setDirection();});

    LineEdit( distanceLineEdit, [=] {setDistance(); });
    LineEdit( lengthLineEdit,   [=] {setLength();});

    CheckBox( displayCheckBox,   [=] {checkDisplay();});

    PushButton( colorPushButton, [=] {setColor();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
