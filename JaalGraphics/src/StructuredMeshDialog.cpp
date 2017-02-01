#include "StructuredMeshDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JStructuredMeshDialog :: JStructuredMeshDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    oldmesh = nullptr;
    newmesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    xorgLineEdit->setText( QString::number(0.0) );
    yorgLineEdit->setText( QString::number(0.0) );
    zorgLineEdit->setText( QString::number(0.0) );

    xlengthLineEdit->setText( QString::number(1.0) );
    ylengthLineEdit->setText( QString::number(1.0) );
    zlengthLineEdit->setText( QString::number(1.0) );

    xdimLineEdit->setText( QString::number(2) );
    ydimLineEdit->setText( QString::number(2) );
    zdimLineEdit->setText( QString::number(1) );

    fixedSpacingLineEdit->setText( QString::number(1.0));
}

///////////////////////////////////////////////////////////////////////////////

JStructuredMeshDialog :: ~JStructuredMeshDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JStructuredMeshDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr) {
        c  = JMeshViewer::registerComponent(viewManager);
        meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    }
    assert( meshViewer );
}

///////////////////////////////////////////////////////////////////////////////
void JStructuredMeshDialog :: setOrigin( double x, double y, double z )
{
    xorgLineEdit->setText( QString::number(x) );
    yorgLineEdit->setText( QString::number(y) );
    zorgLineEdit->setText( QString::number(z) );
}
///////////////////////////////////////////////////////////////////////////////

void JStructuredMeshDialog :: setLength( double x, double y, double z )
{
    xlengthLineEdit->setText( QString::number(x) );
    ylengthLineEdit->setText( QString::number(y) );
    zlengthLineEdit->setText( QString::number(z) );
}
///////////////////////////////////////////////////////////////////////////////

void JStructuredMeshDialog :: setDimension( int x, int y, int z )
{
    xdimLineEdit->setText( QString::number(x) );
    ydimLineEdit->setText( QString::number(y) );
    zdimLineEdit->setText( QString::number(z) );
}
///////////////////////////////////////////////////////////////////////////////
void JStructuredMeshDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return )
    {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JStructuredMeshDialog :: genMesh()
{
    int dim[3];
    double length[3];
    double origin[3];
    QString qstr;

    // Collect length ...
    qstr = xlengthLineEdit->text();
    length[0] = qstr.toDouble();

    qstr = ylengthLineEdit->text();
    length[1] = qstr.toDouble();

    qstr = zlengthLineEdit->text();
    length[2] = qstr.toDouble();

    // Collect Origin ...
    qstr = xorgLineEdit->text();
    origin[0] = qstr.toDouble();

    qstr = yorgLineEdit->text();
    origin[1] = qstr.toDouble();

    qstr = zorgLineEdit->text();
    origin[2] = qstr.toDouble();

    if( fixedSpacingCheckBox->isChecked() )
    {
        qstr = fixedSpacingLineEdit->text();
        double dx = qstr.toDouble();
        dim[0]  = max(1,(int)(length[0]/dx));
        dim[1]  = max(1,(int)(length[1]/dx));
        dim[2]  = max(1,(int)(length[2]/dx));
        xdimLineEdit->setText( QString::number(dim[0]) );
        ydimLineEdit->setText( QString::number(dim[1]) );
        zdimLineEdit->setText( QString::number(dim[2]) );
    }
    else
    {
        // Collect dimension ...
        qstr = xdimLineEdit->text();
        dim[0] = qstr.toInt();

        qstr = ydimLineEdit->text();
        dim[1] = qstr.toInt();

        qstr = zdimLineEdit->text();
        dim[2] = qstr.toInt();
    }

    JWaitCursor wcursor;
    wcursor.start();

    if( dim[2] == 1 ) {
        newmesh = AllQuadMeshGenerator::getStructuredMesh(dim, length, origin, texCoords);
    } else {
        newmesh = AllHexMeshGenerator::getStructuredMesh(dim, length, origin);
    }

    if(meshViewer) meshViewer->addObject( newmesh );
}
///////////////////////////////////////////////////////////////////////////////
void JStructuredMeshDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JStructuredMeshDialog :: makeConnections()
{
    PushButton( applyPushButton,  [=] {genMesh(); });
    PushButton( closePushButton,  [=] {closeDialog(); });
}
///////////////////////////////////////////////////////////////////////////////
