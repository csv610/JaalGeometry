#include "PolyhedraDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JPolyhedraDialog :: JPolyhedraDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    mesh    = nullptr;
    newmesh = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JPolyhedraDialog :: ~JPolyhedraDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JPolyhedraDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}
///////////////////////////////////////////////////////////////////////////////

void JPolyhedraDialog :: platonicSolid()
{
    QString qs = platonicComboBox->currentText();
    string str = qs.toUtf8().constData();

    int n = 0;
    if( str == "Cube")         n = 6;
    if( str == "Tetrahedron")  n = 4;
    if( str == "Octahedron")   n = 8;
    if( str == "Dodecahedron") n = 12;
    if( str == "Icosahedron")  n = 20;

    if( n  == 0 ) return;

    if( newmesh ) newmesh->setActiveBit(0);
    newmesh = JMesh::newObject();

    if( n ==  6) {
        JHexahedronPtr hex = JHexahedron::getCanonical();
        newmesh->addObjects( hex->getNodes() );
        for( int i = 0; i < 6; i++)
            newmesh->addObject(hex->getFaceAt(i));
    } else {
        ostringstream oss;
        oss << "mesh_make platonic " << n << " model.off";

        int err = system( oss.str().c_str() );
        if( err < 0) return;

        newmesh = JMeshIO::readFile( "model.off" );
    }
    if( meshViewer ) meshViewer->addObject( newmesh );
}

////////////////////////////////////////////////////////////////////////////////
void JPolyhedraDialog :: archimedSolid()
{
    QString qs = archimedComboBox->currentText();
    string str = qs.toUtf8().constData();
    cout << str << endl;
}
////////////////////////////////////////////////////////////////////////////////

void JPolyhedraDialog :: genRegularPolytopes()
{
    if( platonicRadioButton->isChecked() ) platonicSolid();
    if( archimedRadioButton->isChecked() ) archimedSolid();
}

///////////////////////////////////////////////////////////////////////////////
void JPolyhedraDialog :: genBoySurface()
{
    if( meshViewer == nullptr ) return;

    JBoySurface boy;

    if( min1RadioButton->isChecked() ) {
        newmesh = boy.getMin1();
    }

    if( min2RadioButton->isChecked() ) {
        cout << "Min-I Boy Surface " << endl;
    }

    if( min2RadioButton->isChecked() ) {
        cout << "Min-I Boy Surface " << endl;
    }

    if( aperyRadioButton->isChecked() ) {
        newmesh = boy.getApery(20,20);
    }

    if( bryantRadioButton->isChecked() ) {
        newmesh = boy.getBryant(10,10);
    }

//   meshViewer->addObject( newmesh );
}


///////////////////////////////////////////////////////////////////////////////
void JPolyhedraDialog :: makeConnections()
{
    PushButton( applyRegularPushButton, [=] {genRegularPolytopes();});
    PushButton( applyBoyPushButton,     [=] {genBoySurface();});
    PushButton( closePushButton,        [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
