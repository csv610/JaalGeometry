#include "HexMesherDialog.hpp"
#include "MeshImporter.hpp"
#include "MeshExporter.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
JHexMesherDialog :: JHexMesherDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    exactinit();
}
///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
JHexMesherDialog :: ~JHexMesherDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JHexMesherDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh(meshViewer->getCurrentMesh() );;
}

///////////////////////////////////////////////////////////////////////////////
void JHexMesherDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: openMeshStackDialog()
{
    if( meshStackDialog.get() == nullptr )
        meshStackDialog.reset(new JMeshStackDialog(this));

    meshStackDialog->setViewManager( viewManager );
    meshStackDialog->show();
}

void JHexMesherDialog :: openLegoMeshDialog()
{
    if( legoMeshDialog.get()== nullptr )
        legoMeshDialog.reset(new JLegoBuilderDialog(this) );

    legoMeshDialog->setViewManager( viewManager );
    legoMeshDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: openStructuredMeshDialog()
{
    if( structMeshDialog.get() == nullptr )
        structMeshDialog.reset(new JStructuredMeshDialog(this));

    structMeshDialog->setViewManager( viewManager );
    structMeshDialog->show();
}
///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: openBernHexOpsDialog()
{
    if( bernHexOpsDialog.get() == nullptr )
        bernHexOpsDialog.reset(new JBernHexOpsDialog(this));

    bernHexOpsDialog->setViewManager( viewManager );
    bernHexOpsDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->close();
}

///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: allTet2HexConversion()
{
    /*
        if( viewer == 0 ) return;

        Mesh *mesh = meshViewer->getMesh();
        if( mesh == 0 ) return;

        int topDim = mesh->getTopology()->getDimension();

        int elemType;
        if( topDim == 2 ) {
            elemType = mesh->getTopology()->getElementsType(2);
            if( elemType != Face::TRIANGLE ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("All Hex conversion require triangular surface mesh. Should I use default converter ?");
                msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel);
                int ret = msg.exec();
                if( ret == QMessageBox::Cancel ) return;
                if( elemType == Face::QUADRILATERAL ) {
                    Mesh *trimesh = AllTriMeshGenerator::getFromQuadMesh( mesh, 2);
                    meshViewer->setNewMesh(trimesh);
                    mesh->deleteAll();
                    delete mesh;
                    mesh = trimesh;
                }
            }
            Mesh *tetmesh = AllTetMeshGenerator::getTetMesh(mesh);
            mesh->deleteAll();
            delete mesh;
            meshViewer->setNewMesh(tetmesh);
        }

        mesh = meshViewer->getMesh();
        elemType = mesh->getTopology()->getElementsType(3);

        if( elemType == Cell::HEXAHEDRON ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("The input mesh is already a hex mesh, you may try cleanup to improve the quality.");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) return;
        }

        if( elemType == Cell::TETRAHEDRON ) {
            Mesh *hexmesh = AllHexMeshGenerator::Tetrahedral2Hexahedral(mesh);
            mesh = hexmesh;
            meshViewer->setNewMesh(mesh);
        }
    */
}

void JHexMesherDialog :: getGeodeTemplate()
{
    if( meshViewer == nullptr ) return;

    JMeshPtr m = AllHexMeshGenerator::getGeodeTemplate();
    meshViewer->addObject(m);
}
///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: getSchneiderTemplate()
{
    if( meshViewer == nullptr ) return;
    JMeshPtr m = AllQuadMeshGenerator::SchneiderPyramid();
    meshViewer->addObject(m);
}
///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: openPolyCubesDialog()
{
    if( polycubesDialog == nullptr )
        polycubesDialog.reset(new JPolyCubesDialog(this));
    polycubesDialog->setViewManager( viewManager );
    polycubesDialog->show();
}
///////////////////////////////////////////////////////////////////////////////
void JHexMesherDialog :: getOctreeMesh()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();
    QString qstr = numNodesPerCellLineEdit->text();
    int nnodes  = qstr.toInt();

    JMeshOctree octree;
    octree.setNumOfPointsPerCell(nnodes);
    octree.setMesh(mesh);
    JMeshPtr  octmesh = octree.getVoxels();
    meshViewer->addObject(octmesh);
}

///////////////////////////////////////////////////////////////////////////////

void JHexMesherDialog :: openSphHexMesherDialog()
{
    if( sphHexMesherDialog == nullptr )
        sphHexMesherDialog.reset(new JSphereHexMesherDialog(this));
    sphHexMesherDialog->setViewManager( viewManager );
    sphHexMesherDialog->setMesh(mesh);
    sphHexMesherDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JHexMesherDialog :: makeConnections()
{
    connect( genOctreePushButton, SIGNAL( clicked() ), this, SLOT( getOctreeMesh() ));
    connect( sphereHexMeshPushButton, SIGNAL( clicked() ), this, SLOT( openSphHexMesherDialog() ));

    connect( structuredPushButton, SIGNAL( clicked() ), this, SLOT( openStructuredMeshDialog() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( closeDialog() ));
}
///////////////////////////////////////////////////////////////////////////////
