#include "FaceAttributesDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JFaceAttributesDialog :: JFaceAttributesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

void JFaceAttributesDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JFaceAttributesDialog :: setMesh(const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}
///////////////////////////////////////////////////////////////////////////////
void JFaceAttributesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    update_default_values = 0;
    numFacesLineEdit->setText( QString::number(faces.size()) );
}

///////////////////////////////////////////////////////////////////////////////
void JFaceAttributesDialog :: setFaces(JFaceSequence &seq)
{
    faces = seq;

    int nCount = 0;
    JFaceRenderPtr attrib;
    for( size_t i = 0; i < faces.size(); i++) {
        int err = faces[i]->getAttribute("Render", attrib);
        if( !err) {
            if( attrib->display ) nCount++;
        }
    }

    bool val = 0;
    if( nCount) val = 1;
    displayCheckBox->setChecked( val );

    if( mesh ) meshViewer->updateBuffers( mesh );
    numFacesLineEdit->setText( QString::number(faces.size()) );
    numVisibleLineEdit->setText( QString::number(nCount) );
}

///////////////////////////////////////////////////////////////////////////////

void JFaceAttributesDialog :: setColor()
{
    QColor color = QColorDialog::getColor();
    if( meshViewer == nullptr ) return;

    float rgb[3];
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;

    JFaceRenderPtr fAttrib;
    if( backRadioButton->isChecked() ) {
        for( size_t i = 0; i < faces.size(); i++) {
            faces[i]->getAttribute("Render", fAttrib);
            if( fAttrib) {
                fAttrib->backColor[0] = rgb[0];
                fAttrib->backColor[1] = rgb[1];
                fAttrib->backColor[2] = rgb[2];
                fAttrib->backColor[3] = 1.0;
                fAttrib->faceSide     = 2;
            }
        }
        return;
    }

    if( frontRadioButton->isChecked() ) {
        for( size_t i = 0; i < faces.size(); i++) {
            faces[i]->getAttribute("Render", fAttrib);
            if( fAttrib) {
                fAttrib->color[0] = rgb[0];
                fAttrib->color[1] = rgb[1];
                fAttrib->color[2] = rgb[2];
                fAttrib->color[3] = 1.0;
                fAttrib->faceSide   = 1;
            }
        }
    }

    if( frontbackRadioButton->isChecked() ) {
        for( size_t i = 0; i < faces.size(); i++) {
            faces[i]->getAttribute("Render", fAttrib);
            if( fAttrib) {
                fAttrib->color[0] = rgb[0];
                fAttrib->color[1] = rgb[1];
                fAttrib->color[2] = rgb[2];
                fAttrib->color[3] = 1.0;
                fAttrib->faceSide = 0;
            }
        }
    }
    if( mesh ) meshViewer->updateBuffers( mesh );
}

//////////////////////////////////////////////////////////////////////////////

void JFaceAttributesDialog :: setMaterial()
{
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JFaceAttributesDialog :: checkDisplay()
{
    JFaceRenderPtr attrib;
    bool val = displayCheckBox->isChecked();
    size_t nCount = 0;
    for( size_t i = 0; i < faces.size(); i++) {
        int err = faces[i]->getAttribute("Render", attrib);
        if(!err) {
            attrib->display = val;
            if( val ) nCount++;
        }
    }

    if( mesh ) meshViewer->updateBuffers(mesh);
    numVisibleLineEdit->setText( QString::number(nCount) );
}
///////////////////////////////////////////////////////////////////////////////


void JFaceAttributesDialog :: makeConnections()
{
    CheckBox( displayCheckBox,  [=] { checkDisplay();});
    PushButton( colorPushButton, [=] { setColor();});
    PushButton( materialPushButton, [=] { setMaterial();});
    PushButton( closePushButton,   [=] { close();});
}

///////////////////////////////////////////////////////////////////////////////
