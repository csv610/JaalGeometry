#include "LightsDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JLightsDialog :: JLightsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    xposLineEdit->setText( QString::number(0.0) );
    yposLineEdit->setText( QString::number(0.0) );
    zposLineEdit->setText( QString::number(1.0) );

    xspotLineEdit->setText( QString::number(0.0) );
    yspotLineEdit->setText( QString::number(0.0) );
    zspotLineEdit->setText( QString::number(-1.0) );

    spotExponentLineEdit->setText( QString::number(0.0) );

    bulbsizeLineEdit->setText( QString::number(0.1) );
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: init()
{
    if( viewManager == nullptr) return;
    lights = viewManager->getLights();
}

///////////////////////////////////////////////////////////////////////////////

JLightsDialog :: ~JLightsDialog()
{
}
///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        viewManager->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: ambientColor()
{
    JColor rgba;
    QColor color = QColorDialog::getColor();
    rgba[0] = color.red()/255.0;
    rgba[1] = color.green()/255.0;
    rgba[2] = color.blue()/255.0;
    rgba[3] = 1.0;

    if( lights) {
        int lnum = lightSpinBox->value();
        lights->setAmbientColor( lnum, rgba );
        viewManager->refreshDisplay();
    }
}
///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: ambientModel()
{
    JColor rgba;

    QColor color = QColorDialog::getColor();
    rgba[0] = color.red()/255.0;
    rgba[1] = color.green()/255.0;
    rgba[2] = color.blue()/255.0;
    rgba[3] = 1.0;

    if( lights) {
        lights->setAmbientModel( rgba );
        viewManager->refreshDisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: diffuseColor()
{
    JColor rgba;

    QColor color = QColorDialog::getColor();
    rgba[0] = color.red()/255.0;
    rgba[1] = color.green()/255.0;
    rgba[2] = color.blue()/255.0;
    rgba[3] = 1.0;

    if( lights) {
        int lnum = lightSpinBox->value();
        lights->setDiffuseColor( lnum, rgba );
        viewManager->refreshDisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: specularColor()
{
    JColor rgba;

    QColor color = QColorDialog::getColor();
    rgba[0] = color.red()/255.0;
    rgba[1] = color.green()/255.0;
    rgba[2] = color.blue()/255.0;
    rgba[3] = 1.0;

    if( lights) {
        int lnum = lightSpinBox->value();
        lights->setSpecularColor( lnum, rgba );
        viewManager->refreshDisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: switchLight()
{
    bool val;
    val = lightSwitchCheckBox->isChecked();
    lights->Switch(val);

    val = light0CheckBox->isChecked();
    lights->Switch( 0, val );

    val = light1CheckBox->isChecked();
    lights->Switch( 1, val );

    val = light2CheckBox->isChecked();
    lights->Switch( 2, val );

    val = light3CheckBox->isChecked();
    lights->Switch( 3, val );

    val = light4CheckBox->isChecked();
    lights->Switch( 4, val );

    val = light5CheckBox->isChecked();
    lights->Switch( 5, val );

    val = light6CheckBox->isChecked();
    lights->Switch( 6, val );

    val = light7CheckBox->isChecked();
    lights->Switch( 7, val );

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: setPosition()
{
    Point4F light_position;
    QString qstr;

    qstr = xposLineEdit->text();
    light_position[0] = qstr.toDouble();

    qstr = yposLineEdit->text();
    light_position[1] = qstr.toDouble();

    qstr = zposLineEdit->text();
    light_position[2] = qstr.toDouble();

    bool val = directionalRadioButton->isChecked();

    if( val  )
        light_position[3] = 0.0;
    else
        light_position[3] = 1.0;

    if( lights) {
        int lnum = lightSpinBox->value();
        lights->setPosition( lnum, light_position );
        viewManager->refreshDisplay();
    }

}
///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: setLightType()
{
    int lnum = lightSpinBox->value();
    Point4F pos = lights->getPosition( lnum );

    bool val = directionalRadioButton->isChecked();
    if( val )
        pos[3] = 0.0;
    else
        pos[3] = 1.0;

    lights->setPosition( lnum, pos );
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: lightModel()
{
    bool val;

    val = twosidedCheckBox->isChecked();
    lights->setTwoSided(val);

    val = localViewerCheckBox->isChecked();
    lights->setLocalViewer(val);
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
/*
void JLightsDialog :: constAttenuate()
{
     if( lights == nullptr) return;

     QString qstr;
     qstr = constantAttenuateLineEdit->text();

     int lid = lightSpinBox->value();
     GLenum light = lightID[lid];
     JLights.getInstance()setConstantAttenuation(light, qstr.toDouble() );
     viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JLightsDialog :: linearAttenuate()
{
     if( lights == nullptr) return;

     QString qstr;
     qstr = linearAttenuateLineEdit->text();
     int lid = lightSpinBox->value();
     GLenum light = lightID[lid];
     JLights.getInstance()setLinearAttenuation(light, qstr.toDouble() );
     viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: quadraticAttenuate()
{
     if( lights == nullptr) return;

     QString qstr;
     qstr = quadraticAttenuateLineEdit->text();
     int lid = lightSpinBox->value();
     GLenum light = lightID[lid];
     JLights.getInstance()setQuadraticAttenuation(light, qstr.toDouble() );
     viewManager->refreshDisplay();

}
*/
///////////////////////////////////////////////////////////////////////////////
void JLightsDialog :: spotDirection()
{
    Point3F dir;

    QString qstr;

    qstr = xspotLineEdit->text();
    dir[0] = qstr.toDouble();

    qstr = yposLineEdit->text();
    dir[1] = qstr.toDouble();

    qstr = zposLineEdit->text();
    dir[2] = qstr.toDouble();

    int lnum = lightSpinBox->value();
//   JLights::getInstance().setSpotDirection( lnum, dir );
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: spotExponent()
{
    QString qstr;
    qstr = spotExponentLineEdit->text();
    int lnum = lightSpinBox->value();
//   JLights::getInstance().setSpotExponent( lnum, qstr.toDouble() );
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: spotCutoff()
{
    double angle = spotCutoffSpinBox->value();

    int lnum = lightSpinBox->value();
//   JLights::getInstance().setSpotCutoff( lnum, angle );
    viewManager->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JLightsDialog :: getLightNum()
{
    /*
         int lnum = lightSpinBox->value();
         Point4F pos = JLights::getInstance().getPosition( lnum );
         xposLineEdit->setText( QString::number(pos[0]) );
         yposLineEdit->setText( QString::number(pos[1]) );
         zposLineEdit->setText( QString::number(pos[2]) );
    */
}

void JLightsDialog :: lightBulbs()
{
    /*
         bool val;
         val = bulbsCheckBox->isChecked();
         JLights::getInstance().setLightBulbs(val);

         QString qstr;
         qstr = bulbsizeLineEdit->text();
         JLights::getInstance().setBulbRadius( qstr.toDouble() );
         viewManager->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////
void JLightsDialog :: makeConnections()
{

    PushButton( ambientPushButton, [=] {ambientColor();});
    PushButton( diffusePushButton, [=] {diffuseColor();});
    PushButton( specularPushButton, [=] {specularColor();});
    PushButton( ambientModelPushButton, [=] {ambientModel(); });

    CheckBox( light0CheckBox, [=] {switchLight();});
    CheckBox( light1CheckBox, [=] {switchLight();});
    CheckBox( light2CheckBox, [=] {switchLight();});
    CheckBox( light3CheckBox, [=] {switchLight();});
    CheckBox( light4CheckBox, [=] {switchLight();});
    CheckBox( light5CheckBox, [=] {switchLight();});
    CheckBox( light6CheckBox, [=] {switchLight();});
    CheckBox( light7CheckBox, [=] {switchLight();});
    CheckBox( bulbsCheckBox,  [=] {lightBulbs();});
    CheckBox( twosidedCheckBox,  [=] {lightModel();});
    CheckBox( localViewerCheckBox, [=] {lightModel();});
    CheckBox( lightSwitchCheckBox, [=] {switchLight();});

    SpinBoxi( lightSpinBox, [=] {getLightNum();});
    SpinBoxi( spotCutoffSpinBox, [=] {spotCutoff();});

    LineEdit( xposLineEdit,  [=] {setPosition();});
    LineEdit( yposLineEdit,  [=] {setPosition();});
    LineEdit( zposLineEdit,  [=] {setPosition();});
    LineEdit( bulbsizeLineEdit,  [=] {lightBulbs();});
    LineEdit( xspotLineEdit,  [=] {spotDirection();});
    LineEdit( yspotLineEdit,  [=] {spotDirection();});
    LineEdit( zspotLineEdit,  [=] {spotDirection();});
    LineEdit( spotExponentLineEdit,  [=] {spotExponent();});

    PushButton( closePushButton, [=] {close();});
}
///////////////////////////////////////////////////////////////////////////////
