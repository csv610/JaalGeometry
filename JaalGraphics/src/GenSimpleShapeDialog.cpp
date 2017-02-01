#include "GenSimpleShapeDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JGenSimpleShapeDialog :: JGenSimpleShapeDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    orgmesh = nullptr;
    newmesh = nullptr;
    meshViewer = nullptr;
    setDefault();
}

///////////////////////////////////////////////////////////////////////////////

JGenSimpleShapeDialog :: ~JGenSimpleShapeDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: setDefault()
{
    helixNLineEdit->setText( QString::number(200) );
    helixMLineEdit->setText( QString::number(50) );
    helixTLineEdit->setText( QString::number(4) );
    helixRLineEdit->setText( QString::number(0.5) );

    torusNLineEdit->setText( QString::number(50) );
    torusMLineEdit->setText( QString::number(20) );
    torusRLineEdit->setText( QString::number(0.2) );

    kleinBottleNLineEdit->setText( QString::number(100) );
    kleinBottleMLineEdit->setText( QString::number(20) );

    knotNLineEdit->setText( QString::number(200) );
    knotMLineEdit->setText( QString::number(20) );
    knotRLineEdit->setText( QString::number(0.2) );

    sphereNLineEdit->setText( QString::number(4));
    sphereMLineEdit->setText( QString::number(4));
}
///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: init()
{
    if( viewManager == nullptr) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c == nullptr) c = JMeshViewer::registerComponent(viewManager);
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////
int JGenSimpleShapeDialog ::genSphere()
{
    static int nCount = 0;
    QString nstr = sphereNLineEdit->text();
    QString mstr = sphereMLineEdit->text();

    int n = nstr.toInt();
    int m = mstr.toInt();
    if( n < 1 || m < 1 ) return 1;

    ostringstream oss;
    if( subdivSphereCheckBox->isChecked() ) {
        if( m > 8) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Do you really want high refinement(2nd Box): May be too slow");
            msg.setStandardButtons( QMessageBox::Cancel | QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Cancel ) return 1;
        }

        if( n ==4 || n == 6 || n == 8 || n == 12 || n == 20 ) {
            oss << "mesh_make ssphere " << n << " " << m << " " << "model.off";
        } else {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("The valid values in the first box are 4 6 8 12 and 20 ");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) return 1;
        }
    }

    if( !subdivSphereCheckBox->isChecked() )
        oss << "mesh_make sphere " << n << " " << m << " " << "model.off";

    int err = system( oss.str().c_str() );
    if( err < 0) return 1;

    newmesh = JMeshIO::readFile( "model.off" );
    ostringstream oss1;
    oss1 << "Sphere" << nCount++;
    newmesh->setName( oss1.str() );

    if( meshViewer ) {
        meshViewer->addObject( newmesh );
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog :: genHelix()
{
    static int nCount = 0;

    QString mstr = helixMLineEdit->text();
    QString nstr = helixNLineEdit->text();
    QString tstr = helixTLineEdit->text();
    QString rstr = helixRLineEdit->text();

    int m = mstr.toInt();
    int n = nstr.toInt();
    int t = tstr.toInt();
    double r = rstr.toDouble();

    if( n < 1 || m < 1 ) return 1;

    ostringstream oss;
    oss << "mesh_make helix " << n << " " << m << " " << t << " " << r  << " model.off";

    int err = system( oss.str().c_str() );
    if( err < 0) return 1;

    newmesh = JMeshIO::readFile( "model.off" );

    ostringstream oss1;
    oss1 << "Helix" << nCount++;
    newmesh->setName( oss1.str() );

    if( meshViewer ) {
        meshViewer->addObject( newmesh );
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JGenSimpleShapeDialog :: genTorus()
{
    static int nCount = 0;

    QString mstr = torusMLineEdit->text();
    QString nstr = torusNLineEdit->text();
    QString rstr = torusRLineEdit->text();

    int m = mstr.toInt();
    int n = nstr.toInt();
    double r = rstr.toDouble();

    if( n < 1 || m < 1 ) return 1;

    ostringstream oss;
    oss << "mesh_make torus " << n << " " << m << " " <<  r  << " model.off";

    int err = system( oss.str().c_str() );
    if( err < 0) return 1;

    newmesh = JMeshIO::readFile( "model.off" );

    ostringstream oss1;
    oss1 << "Torus" << nCount++;
    newmesh->setName( oss1.str() );

    if( meshViewer ) {
        meshViewer->addObject( newmesh );
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog ::genKleinBottle()
{
    static int nCount = 0;
    QString mstr = kleinBottleMLineEdit->text();
    QString nstr = kleinBottleNLineEdit->text();

    int m = mstr.toInt();
    int n = nstr.toInt();

    if( n < 1 || m < 1 ) return 1;

    ostringstream oss;
    oss << "mesh_make klein " << n << " " << m << " model.off";

    int err = system( oss.str().c_str() );
    if( err < 0) return 1;

    newmesh = JMeshIO::readFile( "model.off" );
    ostringstream oss1;
    oss1 << "KleinBottle" << nCount++;
    newmesh->setName( oss1.str() );

    if( meshViewer ) {
        meshViewer->addObject(newmesh);
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog ::genKnot()
{
    static int nCount = 0;

    QString nstr = knotNLineEdit->text();
    QString mstr = knotMLineEdit->text();
    QString rstr = knotRLineEdit->text();

    int n = nstr.toInt();
    int m = mstr.toInt();
    double r = rstr.toDouble();

    ostringstream oss;

    oss << "mesh_make knot " << n <<  " " << m << " " << r <<  " model.off";
    int err = system( oss.str().c_str() );
    if( err < 0) return 1;

    newmesh =  JMeshIO::readFile("model.off");

    ostringstream oss1;
    oss1 << "Knot" << nCount++;
    newmesh->setName( oss1.str() );

    if( meshViewer ) {
        meshViewer->addObject( newmesh );
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: genShape()
{
    if( meshViewer == nullptr ) return;

    if( helixRadioButton->isChecked() ) {
        genHelix();
        return;
    }

    if( torusRadioButton->isChecked() ) {
        genTorus();
        return;
    }

    if( kleinBottleRadioButton->isChecked() ) {
        genKleinBottle();
        return;
    }

    if( knotRadioButton->isChecked() ) {
        genKnot();
        return;
    }

    if( sphereRadioButton->isChecked() ) {
        genSphere();
        return;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: acceptShape()
{
    if( newmesh == nullptr ) {
        this->close();
        return;
    }

    /*
        // When you accept the new shape, you delete the original mesh. The new
        // shape is already with the meshViewer.
        if( orgmesh ) {
            orgmesh->deleteAll();
    //      delete orgmesh;
            orgmesh = nullptr;
        }
    */
    this->close();
}
///////////////////////////////////////////////////////////////////////////////
void JGenSimpleShapeDialog :: cancelShape()
{
    // When you wish you reject the newmesh, the meshViewer must be rolled back
    // to the original mesh.

    if( newmesh == nullptr ) {
        this->close();
        return;
    }

    if( orgmesh == newmesh ) {
        this->close();
        return;
    }

    newmesh->deleteAll();
//  delete newmesh;
//   newmesh = nullptr;

//  meshViewer->setNewMesh( orgmesh );
    this->close();
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: openLSystemDialog()
{
    if( lsystemDialog == nullptr )
        lsystemDialog.reset(new LSystemDialog(this));
    lsystemDialog->setViewManager( viewManager );
    lsystemDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: getPolygon()
{
    const vector<Point3D> &polyPoints = polygonDialog->getPoints();
    int nSize = polyPoints.size();

    if( nSize < 3) return;

    JNodeSequence nodes(nSize);
    for( int i = 0; i < nSize; i++) {
        JNodePtr v = JNode::newObject();
        v->setXYZCoords( polyPoints[i] );
        nodes[i] = v;
    }

    JFacePtr newface;
    switch(nSize) {
    case 3:
        newface = JTriangle::newObject( nodes);
        break;
    case 4:
        newface = JQuadrilateral::newObject( nodes);
        break;
    default:
        newface = JPolygon::newObject( nodes);
        break;
    }

    if( newface ) {
        JMeshPtr newmesh = JMesh::newObject();
        newmesh->addObjects(nodes);
        newmesh->addObject(newface);
        meshViewer->addObject(newmesh);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: openPolygonDialog()
{
    if( polygonDialog == nullptr ) {
        polygonDialog.reset(new JPolygonDialog(this));
        viewManager->attach( polygonDialog.get() );
    }

    QObject::connect( polygonDialog.get(), SIGNAL(setPolygon()), this, SLOT(getPolygon()));

    polygonDialog->setViewManager( viewManager );
    polygonDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: openPolyhedraDialog()
{
    if( polyhedraDialog == nullptr )
        polyhedraDialog.reset(new JPolyhedraDialog(this));
    polyhedraDialog->setViewManager( viewManager );
    polyhedraDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: openCurveDialog()
{
    if( curveDialog == nullptr )
        curveDialog.reset(new JCurveGenDialog(this));
    curveDialog->setViewManager( viewManager );
    curveDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: openKnotsDialog()
{
    if( knotsDialog == nullptr )
        knotsDialog.reset(new JKnotsDialog(this));
    knotsDialog->setViewManager( viewManager );
    knotsDialog->show();
    this->hide();
}


///////////////////////////////////////////////////////////////////////////////
void JGenSimpleShapeDialog :: openImageEdgesDialog()
{
    if( imageEdgesDialog == nullptr )
        imageEdgesDialog.reset( new JImageContoursDialog( this));

    imageEdgesDialog->setViewManager( viewManager);
    imageEdgesDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JGenSimpleShapeDialog :: makeConnections()
{
    PushButton( knotPushButton,       [=] {openKnotsDialog();});
    PushButton( curvePushButton,      [=] {openCurveDialog();});
    PushButton( polygonPushButton,    [=] {openPolygonDialog();});
    PushButton( polyhedraPushButton,  [=] {openPolyhedraDialog();});
    PushButton( lsystemPushButton,    [=] {openLSystemDialog();});
    PushButton( imageEdgesPushButton, [=] {openImageEdgesDialog();});

    PushButton( generatePushButton,   [=] {genShape();});
    PushButton( acceptPushButton,     [=] {acceptShape();});
    PushButton( cancelPushButton,     [=] {cancelShape();});
}

///////////////////////////////////////////////////////////////////////////////

/*
int JGenSimpleShapeDialog :: genPlane()
{
     QString mstr =   planeXLineEdit->text();
     QString nstr =   planeYLineEdit->text();

     int m = mstr.toInt();
     int n = nstr.toInt();

     if( m < 2 || n < 2 ) return 1;

     ostringstream oss;
     oss << "mesh_make plane " << m << " " << n << " " << "model.off";
     cout << oss.str() << endl;
     int err = system( oss.str().c_str() );

     if( err < 0) return 1;

     if( newmesh == nullptr ) newmesh = new Mesh;
     if( newmesh)  newmesh->deleteAll();

     newmesh->readFromFile( "model.off" );

     if( meshViewer ) {
          meshViewer->setNewMesh( newmesh );
          newmesh = nullptr;
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog :: genDisc()
{
     QString nstr = discNLineEdit->text();
     QString mstr = discNLineEdit->text();

     int n = nstr.toInt();
     int m = mstr.toInt();
     if( n < 1 || m < 1 ) return 1;

     ostringstream oss;
     oss << "mesh_make disc " << n << " " << m << " " << "model.off";
     cout << oss.str() << endl;

     int err = system( oss.str().c_str() );
     if( err < 0) return 1;

     if( newmesh == nullptr ) newmesh = new Mesh;
     if( newmesh)  newmesh->deleteAll();

     newmesh->readFromFile( "model.off" );

     if( meshViewer ) {
          meshViewer->setNewMesh( newmesh );
          newmesh = nullptr;
     }

     return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog :: genCube()
{
     QString qstr;
     int dim[3];

     qstr = xcubeLineEdit->text();
     dim[0] = qstr.toInt();

     qstr = ycubeLineEdit->text();
     dim[1] = qstr.toInt();

     qstr = zcubeLineEdit->text();
     dim[2] = qstr.toInt();

     if( newmesh == nullptr ) newmesh = new Mesh;
     if( newmesh)  newmesh->deleteAll();

     newmesh = AllHexMeshGenerator::getStructuredMesh( dim );

     if( meshViewer ) {
          meshViewer->setNewMesh( newmesh );
          newmesh = nullptr;
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog :: genCyl()
{
     QString nstr = cylinderNLineEdit->text();
     QString mstr = cylinderMLineEdit->text();
     QString rstr = cylinderRLineEdit->text();

     int n = nstr.toInt();
     int m = mstr.toInt();
     double r = rstr.toDouble();

     ostringstream oss;

     if( capCylinderCheckBox->isChecked() )
          oss << "mesh_make ccyl " << n <<  " " << m << " " << r <<  " model.off";
     else
          oss << "mesh_make cyl " << n <<  " " << m << " " << r <<  " model.off";

     int err = system( oss.str().c_str() );
     if( err < 0) return 1;

     if( newmesh == nullptr ) newmesh = new Mesh;
     if( newmesh)  newmesh->deleteAll();

     newmesh->readFromFile( "model.off" );

     if( meshViewer ) {
          meshViewer->setNewMesh( newmesh );
          newmesh = nullptr;
     }

     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog :: genCone()
{
     QString nstr = coneNLineEdit->text();
     QString mstr = coneMLineEdit->text();
     QString rstr = coneRLineEdit->text();

     int n = nstr.toInt();
     int m = mstr.toInt();
     double r = rstr.toDouble();

     ostringstream oss;

     if( capConeCheckBox->isChecked() )
          oss << "mesh_make ccone " << n <<  " " << m << " " << r <<  " model.off";
     else
          oss << "mesh_make cone " << n <<  " " << m << " " << r <<  " model.off";

     int err = system( oss.str().c_str() );
     if( err < 0) return 1;

     if( newmesh == nullptr ) newmesh = new Mesh;
     if( newmesh)  newmesh->deleteAll();

     newmesh->readFromFile( "model.off" );

     if( meshViewer ) {
          meshViewer->setNewMesh( newmesh );
          newmesh = nullptr;
     }
     return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JGenSimpleShapeDialog ::genPlatonic()
{
    static int nCount = 0;
    QString nstr = platonicLineEdit->text();

    int n = nstr.toInt();

    ostringstream oss;

    oss << "mesh_make platonic " << n <<  " model.off";
    int err = system( oss.str().c_str() );
    if( err < 0) return 1;

    newmesh = JMeshIO::readFile( "model.off" );
    ostringstream oss1;
    oss1 << "Platonic" << nCount++;
    newmesh->setName( oss1.str() );

    if( meshViewer ) {
        meshViewer->addNewMesh( newmesh );
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

*/
