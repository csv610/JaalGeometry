#include "MeshSubdivisionDialog.hpp"
///////////////////////////////////////////////////////////////////////////////

JMeshSubdivisionDialog :: JMeshSubdivisionDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshSubdivisionDialog :: ~JMeshSubdivisionDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSubdivisionDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSubdivisionDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );
}
///////////////////////////////////////////////////////////////////////////////

void JMeshSubdivisionDialog :: getRefinedMesh()
{
    if( mesh == nullptr) return;
    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2) {
        cout << "Warning: Subdivion schenmes can ba applied to surface mesh only" << endl;
        return;
    }

    JMeshSubdivision msubdiv;

    QString qstr = algorithmComboBox->currentText();
    string  str  = StdString(qstr);

    if( str == "CatmullClark")
        msubdiv.setAlgorithm(JMeshSubdivision::CATMULL_CLARK);

    if( str == "DooSabin")
        msubdiv.setAlgorithm(JMeshSubdivision::DOOSABIN);

    if( str == "Sqrt3")
        msubdiv.setAlgorithm(JMeshSubdivision::SQRT3);

    if( str == "Loop") {
        int entity = mesh->getTopology()->getElementsType(2);
        if( entity != JFace::TRIANGLE) {
            cout << "Warning; : :oop subdivion is applicable to triangle mesh only " << endl;
            return;
        }
        msubdiv.setAlgorithm(JMeshSubdivision::LOOP);
    }

    int iter = numIterationsSpinBox->value();
    msubdiv.setIterations(iter);

    JWaitCursor wcursor;
    wcursor.start();

    msubdiv.setMesh(mesh);
    JMeshPtr newmesh = msubdiv.getSubdivided();
    if( newmesh ) {
        stringstream ss;
        ss << mesh->getName() << iter;
        newmesh->setName( ss.str() ) ;

        // Remove old mesh ..
        mesh->deleteAll();
        meshViewer->removeObject(mesh);

        // Add new mesh ..
        meshViewer->addObject(newmesh);
        setMesh(newmesh);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSubdivisionDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSubdivisionDialog :: makeConnections()
{
    PushButton( applyPushButton,  [=] {getRefinedMesh(); });
    PushButton( closePushButton,  [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
