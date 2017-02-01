#include "FontsDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

JFontsDialog :: JFontsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();
    viewManager = nullptr;

    fontSizeLineEdit->setText( QString::number(1) );
}

///////////////////////////////////////////////////////////////////////////////

void JFontsDialog :: setFontSize()
{
    if( viewManager == nullptr ) return;

    QString str = fontSizeLineEdit->text();
    FontsManager::Instance().setFontScale( str.toDouble() );

    viewManager->updateGL();
}

///////////////////////////////////////////////////////////////////////////////
void JFontsDialog :: xRotate()
{
    if( viewManager == nullptr ) return;
    double angle = (double) xangleSpinBox->value();

    Point3F rot;
    rot[0] = angle;
    rot[1] = 0.0;
    rot[2] = 0.0;
    FontsManager::Instance().setRotateAngles( rot );

    viewManager->updateGL();
}
///////////////////////////////////////////////////////////////////////////////

void JFontsDialog :: yRotate()
{
    if( viewManager == nullptr ) return;

    double angle = (double) yangleSpinBox->value();

    Point3F rot;
    rot[0] = 0.0;
    rot[1] = angle;
    rot[2] = 0.0;
    FontsManager::Instance().setRotateAngles( rot );

    viewManager->updateGL();
}
///////////////////////////////////////////////////////////////////////////////

void JFontsDialog :: zRotate()
{
    if( viewManager == nullptr ) return;

    double angle = (double) zangleSpinBox->value();

    Point3F rot;
    rot[0] = 0.0;
    rot[1] = 0.0;
    rot[2] = angle;
    FontsManager::Instance().setRotateAngles( rot );

    viewManager->updateGL();
}
///////////////////////////////////////////////////////////////////////////////


void JFontsDialog :: setColor()
{
    if( viewManager == nullptr ) return;

    QColor color = QColorDialog::getColor();
    JColor rgba;
    rgba[0] = color.red()/255.0;
    rgba[1] = color.green()/255.0;
    rgba[2] = color.blue()/255.0;
    rgba[3] = 1.0;
    FontsManager::Instance().setColor( rgba );
    viewManager->updateGL();
}
///////////////////////////////////////////////////////////////////////////////

void JFontsDialog :: makeConnections()
{
    connect( xangleSpinBox, SIGNAL( valueChanged(int) ), this, SLOT( xRotate() ));
    connect( yangleSpinBox, SIGNAL( valueChanged(int) ), this, SLOT( yRotate() ));
    connect( zangleSpinBox, SIGNAL( valueChanged(int) ), this, SLOT( zRotate() ));
    connect( colorPushButton,  SIGNAL( clicked() ), this, SLOT( setColor() ));
    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( close() ));
    connect( fontSizeLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setFontSize() ));
}

///////////////////////////////////////////////////////////////////////////////
