#include "CongruentSetDialog.hpp"
#include "AllTriMeshGenerator.hpp"

JCongruentSetDialog :: JCongruentSetDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    angleTolLineEdit->setText( QString::number(0) );
}

///////////////////////////////////////////////////////////////////////////////

JCongruentSetDialog :: ~JCongruentSetDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JCongruentSetDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;
//  mesh = meshViewer->getMesh();
    cset.setMesh(mesh);

    mesh->deleteNodeAttribute("Constraint");
}

///////////////////////////////////////////////////////////////////////////////

void JCongruentSetDialog :: generateSet()
{
    if( mesh == nullptr ) return;

    QString qstr = angleTolLineEdit->text();

    cset.setTolerance( qstr.toDouble() );
    cset.apply();

    int nsize = cset.getSetSize();
    cout << "Info: Creating Congruent file " << endl;
    ofstream ofile( "congru.dat", ios::out);
    vector<double> angles;
    for( int i = 0; i < nsize; i++) {
        JFaceSequence faces = cset.getSet(i);
        ofile << "Group " << i <<  " Size " << faces.size() << endl;
        for( int j = 0; j < faces.size(); j++) {
            JFaceGeometry::getAngles( faces[j], angles, ANGLE_IN_DEGREES);
            sort( angles.begin(), angles.end() );
            ofile << std::fixed << std::setprecision(10) << angles[0] << "  " << angles[1] << " " << angles[2] << endl;
        }
    }

    numGroupsLineEdit->setText( QString::number(nsize) );

    displayGroupSpinBox->setMinimum(0);
    displayGroupSpinBox->setMaximum(nsize-1);

    JNodeRenderPtr nAttrib;
    JColor greenColor, redColor;
    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 1.0;

    greenColor[0] = 0.0;
    greenColor[1] = 1.0;
    greenColor[2] = 0.0;
    greenColor[3] = 1.0;

    size_t nSize = mesh->getSize(0);
    for( size_t i = 0; i < nSize; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Render", nAttrib);

        if( vtx->hasAttribute("Constraint") ) {
            nAttrib->display = 1;
            nAttrib->color   = redColor;
            nAttrib->glyph   = JNodeRender::NODE_AS_SPLAT;
        } else {
            nAttrib->display = 0;
            nAttrib->color   = greenColor;
            nAttrib->glyph   = JNodeRender::NODE_AS_POINT;
        }
    }

    displayGroup();
}

///////////////////////////////////////////////////////////////////////////////

void JCongruentSetDialog :: displayGroup()
{
    if( mesh == nullptr ) return;

//    meshViewer->displayAll(1,1);

    JFaceRenderPtr fAttrib;
    JColor clr;
    clr[0] = 0.2;
    clr[1] = 1.0;
    clr[2] = 0.3;
    clr[3] = 1.0;

    size_t nSize = mesh->getSize(2);
    for( size_t i = 0; i < nSize; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
        fAttrib->color   = clr;
    }
    int gid  = displayGroupSpinBox->value();
    JFaceSequence faces = cset.getSet(gid);
    numElementsLineEdit->setText( QString::number(faces.size()) );

    clr[0] = 1.0;
    clr[1] = 0.2;
    clr[2] = 0.2;
    clr[3] = 1.0;
    nSize = faces.size();
    vector<double> angles;
    for( size_t i = 0; i < faces.size(); i++) {
        faces[i]->getAttribute("Render", fAttrib);
//      FaceGeometry::getAngles( faces[i], angles );
        fAttrib->display = 1;
        fAttrib->color   = clr;
    }

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JCongruentSetDialog :: getCongruent()
{
    AllTriMeshGenerator::getCongruentMesh(mesh);
    generateSet();
}

///////////////////////////////////////////////////////////////////////////////
void JCongruentSetDialog :: makeConnections()
{
    connect( applyPushButton,  SIGNAL( clicked() ), this, SLOT( generateSet() ));
    connect( displayGroupSpinBox, SIGNAL( valueChanged(int) ), this, SLOT( displayGroup() ));
    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( close() ));
    connect( getCongruentPushButton,  SIGNAL( clicked() ), this, SLOT( getCongruent() ));
}

///////////////////////////////////////////////////////////////////////////////
