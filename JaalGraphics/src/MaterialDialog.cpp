#include "MaterialDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMaterialDialog :: JMaterialDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    meshViewer = nullptr;
    material   = nullptr;

    nameLineEdit->setText(QString("unnamed"));
}

///////////////////////////////////////////////////////////////////////////////

JMaterialDialog :: ~JMaterialDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMaterialDialog :: setColor()
{
    QColor color = QColorDialog::getColor();

    if( meshViewer == nullptr ) return;

    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
}
///////////////////////////////////////////////////////////////////////////////

void JMaterialDialog :: setAmbient()
{
    setColor();

    if( material) material->setAmbientColor(rgb);
    if( meshViewer) meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMaterialDialog :: setDiffuse()
{
    setColor();

    if( material) material->setDiffuseColor(rgb);
    if( meshViewer) meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMaterialDialog :: setSpecular()
{
    setColor();
    if( material) material->setSpecularColor(rgb);
    if( meshViewer) meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMaterialDialog :: setEmission()
{
    setColor();
    if( material) material->setEmissionColor(rgb);
    if( meshViewer) meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMaterialDialog :: setStdMaterial()
{
    if( predefinedCheckBox->isChecked() ) {
        QString qs = predefinedComboBox->currentText();
        string str = qs.toUtf8().constData();
        if( material) material->setMaterial(str);
    }
    if( meshViewer) meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMaterialDialog :: makeConnections()
{
    connect( predefinedCheckBox, SIGNAL( toggled(bool) ) , this, SLOT( setStdMaterial() ));
    connect( predefinedComboBox, SIGNAL( activated(int) ), this, SLOT( setStdMaterial() ));
    connect( ambientPushButton,  SIGNAL( clicked() ), this, SLOT( setAmbient() ));
    connect( diffusePushButton,  SIGNAL( clicked() ), this, SLOT( setDiffuse() ));
    connect( emissionPushButton, SIGNAL( clicked() ), this, SLOT( setEmission() ));
    connect( specularPushButton, SIGNAL( clicked() ), this, SLOT( setSpecular() ));
    connect( closePushButton,    SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
