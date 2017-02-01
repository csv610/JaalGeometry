#include "MeshAffineTransformsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshAffineTransformsDialog :: JMeshAffineTransformsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;

    xtranslateLineEdit->setText( QString::number(0.0) );
    ytranslateLineEdit->setText( QString::number(0.0) );
    ztranslateLineEdit->setText( QString::number(0.0) );

    angleLineEdit->setText( QString::number(0.0) );

    xscaleLineEdit->setText( QString::number(1.0) );
    yscaleLineEdit->setText( QString::number(1.0) );
    zscaleLineEdit->setText( QString::number(1.0) );

    xcenterLineEdit->setText( QString::number(0.0) );
    ycenterLineEdit->setText( QString::number(0.0) );
    zcenterLineEdit->setText( QString::number(0.0) );

    xlengthLineEdit->setText( QString::number(0.0) );
    ylengthLineEdit->setText( QString::number(0.0) );
    zlengthLineEdit->setText( QString::number(0.0) );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransformsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    affine.reset( new JMeshAffineTransform);
}

////////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransformsDialog :: geomInfo()
{
    if( mesh == nullptr) return;

    Point3D pC = mesh->getGeometry()->getCenter();
    xcenterLineEdit->setText( QString::number( pC[0] ) );
    ycenterLineEdit->setText( QString::number( pC[1] ) );
    zcenterLineEdit->setText( QString::number( pC[2] ) );

    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    mesh->setAttribute("AxisBoundingBox", box);
    xlengthLineEdit->setText( QString::number( box.getLength(0)) );
    ylengthLineEdit->setText( QString::number( box.getLength(1)) );
    zlengthLineEdit->setText( QString::number( box.getLength(2)) );
}
////////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransformsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    objectNameLineEdit->setText( QString(mesh->getName().c_str() ) );
    geomInfo();
    affine->setMesh(mesh);
}
////////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransformsDialog :: closeDialog()
{
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransformsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransformsDialog :: applyTransform()
{
    if( mesh == nullptr ) return;

    if( permanentChangesCheckBox->isChecked() == 0) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    QString qstr;
    double xt, yt, zt, angle;

    affine->setMesh(mesh);

    if( translateRadioButton->isChecked() ) {
        qstr = xtranslateLineEdit->text();
        xt = qstr.toDouble();
        qstr = ytranslateLineEdit->text();
        yt = qstr.toDouble();
        qstr = ztranslateLineEdit->text();
        zt = qstr.toDouble();
        affine->translate(xt, yt, zt);
    }

    if( rotateRadioButton->isChecked() ) {
        qstr  = angleLineEdit->text();
        angle = M_PI*qstr.toDouble()/180.0;
        if( xRotateRadioButton->isChecked() )
            affine->rotate(angle, 0);

        if( yRotateRadioButton->isChecked() )
            affine->rotate(angle, 1);

        if( zRotateRadioButton->isChecked() )
            affine->rotate(angle, 2);
    }

    if( scaleRadioButton->isChecked() ) {
        qstr = xscaleLineEdit->text();
        xt = qstr.toDouble();

        qstr = yscaleLineEdit->text();
        yt = qstr.toDouble();

        qstr = zscaleLineEdit->text();
        zt = qstr.toDouble();
        affine->scale(xt, yt, zt);
    }

    if( fixedLengthRadioButton->isChecked() ) {
        qstr = setXlenLineEdit->text();
        xt = qstr.toDouble();

        qstr = setYlenLineEdit->text();
        yt = qstr.toDouble();

        qstr = setZlenLineEdit->text();
        zt = qstr.toDouble();
        affine->scale(xt, yt, zt);
    }
    geomInfo();
    meshViewer->updateGeometryBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransformsDialog :: normalize()
{
    if( mesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    affine->normalize();
    meshViewer->updateGeometryBuffers( mesh);

    geomInfo();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransformsDialog :: toCenter()
{
    if( mesh == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    Point3D pC = mesh->getGeometry()->getCenter();

    affine->translate(-pC[0], -pC[1], -pC[2] );

    meshViewer->updateGeometryBuffers( mesh);

    geomInfo();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshAffineTransformsDialog :: setNewNodeID()
{
    if( mesh == nullptr ) return;

    QString qstr;
    qstr = nodeIDLineEdit->text();

    int id  = qstr.toInt();
    if( id >= mesh->getSize(0) ) return;
    JNodePtr v = mesh->getNodeAt(id);

    Point3D p3d = v->getXYZCoords();

    updateXCoordLineEdit->setText( QString::number(p3d[0]) );
    updateYCoordLineEdit->setText( QString::number(p3d[1]) );
    updateZCoordLineEdit->setText( QString::number(p3d[2]) );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransformsDialog :: changeNodeCoords()
{
    if( mesh == nullptr ) return;

    QString qstr;
    qstr = nodeIDLineEdit->text();
    size_t id  = qstr.toInt();
    if( id >= mesh->getSize(0) ) return;
    JNodePtr v = mesh->getNodeAt(id);
    if( !v->isActive() ) return;

    Point3D p3d;
    qstr = updateXCoordLineEdit->text();
    p3d[0]   = qstr.toDouble();

    qstr = updateYCoordLineEdit->text();
    p3d[1]   = qstr.toDouble();

    qstr = updateZCoordLineEdit->text();
    p3d[2]   = qstr.toDouble();
    v->setXYZCoords(p3d);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshAffineTransformsDialog :: makeConnections()
{
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
    connect( applyPushButton, SIGNAL( clicked() ), this, SLOT( applyTransform() ));
    connect( normalizePushButton, SIGNAL( clicked() ), this, SLOT( normalize() ));
    connect( toCenterPushButton, SIGNAL( clicked() ), this, SLOT( toCenter() ));
    connect( updateNodeCoordsPushButton, SIGNAL( clicked() ), this, SLOT( changeNodeCoords() ));
    connect( nodeIDLineEdit,  SIGNAL( editingFinished() ) , this, SLOT( setNewNodeID() ));
}

///////////////////////////////////////////////////////////////////////////////
