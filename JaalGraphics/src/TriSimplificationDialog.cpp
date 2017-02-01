#include "TriSimplificationDialog.hpp"

JTriSimplificationDialog :: JTriSimplificationDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JTriSimplificationDialog :: ~JTriSimplificationDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JTriSimplificationDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

//   mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;
}

///////////////////////////////////////////////////////////////////////////////

void JTriSimplificationDialog :: write_smf()
{
    if( mesh == nullptr ) return;

    ofstream ofile("infile.smf", ios::out);
    if( ofile.fail() ) return;

    size_t numnodes  = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            Point3D xyz = vtx->getXYZCoords();
            ofile << "v " << xyz[0] << "  " << xyz[1] << " " << xyz[2] << endl;
        }
    }

    size_t numfaces  = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            ofile << "f " <<  face->getNodeAt(0)->getID() + 1 << " "
                  <<  face->getNodeAt(1)->getID() + 1 << " "
                  <<  face->getNodeAt(2)->getID() + 1 << endl;
        }
    }
//  smfmodel.reset(new MxStdModel(100,100));
}

void JTriSimplificationDialog :: startup_and_init()
{
//   smfReader.reset(new MxSMFReader);
}


void JTriSimplificationDialog :: accept()
{
//   meshViewer->getDrawNode()->setColorMethod( nullptr );
//   this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JTriSimplificationDialog :: makeConnections()
{
    PushButton(closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
