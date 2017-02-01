#include "LegoBuilderDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JLegoBuilderDialog :: JLegoBuilderDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    oldmesh = nullptr;
    newmesh = nullptr;

    viewManager = nullptr;
    meshViewer  = nullptr;

    nxLineEdit->setText( QString::number(1) );
    nyLineEdit->setText( QString::number(1) );
    nzLineEdit->setText( QString::number(1) );
}

///////////////////////////////////////////////////////////////////////////////

JLegoBuilderDialog :: ~JLegoBuilderDialog()
{
}

///////////////////////////////////////////////////////////////////////////////


int JLegoBuilderDialog :: getCellID( int i, int j, int k)
{
    return k*boxDim[0]*boxDim[1] + j*boxDim[1] + i;
}

///////////////////////////////////////////////////////////////////////////////
void JLegoBuilderDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

//  oldmesh = meshViewer->getMesh();
    QString qstr;
    qstr = nxLineEdit->text();
    boxDim[0] = qstr.toInt();

    qstr = nyLineEdit->text();
    boxDim[1] = qstr.toInt();

    qstr = nzLineEdit->text();
    boxDim[2] = qstr.toInt();

    /*
        boxDim[0] = 5;
        boxDim[1] = 5;
        boxDim[2] = 5;

        int dim[]  = {boxDim[0]+1, boxDim[1]+1, boxDim[2]+1};
        double length[] = {1,1,1};
        double origin[] = {0.0, 0.0, 0.0};

        if( newmesh ) {
            newmesh->deleteAll();
            delete newmesh;
        }
    */
    Point3D p3d;
    JMeshPtr legomesh = JMesh::newObject();
    JEdgePtr edge;

    JNodeSequence nodes;

    /*
       p3d[0] = 0.0;
       p3d[1] = 0.0;
       p3d[2] = 0.0;
       Vertex *vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       p3d[0] = 1.0;
       p3d[1] = 0.0;
       p3d[2] = 0.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       edge = Edge::newObject( legomesh->getNodeAt(0), legomesh->getNodeAt(1));
       legomesh->addObject(edge);

       p3d[0] =  1.0;
       p3d[1] =  0.0;
       p3d[2] = -1.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       edge = Edge::newObject( legomesh->getNodeAt(1), legomesh->getNodeAt(2));
       legomesh->addObject(edge);

       p3d[0] =  0.0;
       p3d[1] =  0.0;
       p3d[2] = -1.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       edge = Edge::newObject( legomesh->getNodeAt(2), legomesh->getNodeAt(3));
       legomesh->addObject(edge);

       edge = Edge::newObject( legomesh->getNodeAt(3), legomesh->getNodeAt(0));
       legomesh->addObject(edge);

       p3d[0] =  0.0;
       p3d[1] =  0.0;
       p3d[2] =  1.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       edge = Edge::newObject( legomesh->getNodeAt(0), legomesh->getNodeAt(4));
       legomesh->addObject(edge);

       p3d[0] =   0.0;
       p3d[1] =  -1.0;
       p3d[2] =   1.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       edge = Edge::newObject( legomesh->getNodeAt(4), legomesh->getNodeAt(5));
       legomesh->addObject(edge);

       p3d[0] =   0.0;
       p3d[1] =  -1.0;
       p3d[2] =   0.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       edge = Edge::newObject( legomesh->getNodeAt(5), legomesh->getNodeAt(6));
       legomesh->addObject(edge);

       edge = Edge::newObject( legomesh->getNodeAt(6), legomesh->getNodeAt(0));
       legomesh->addObject(edge);

       p3d[0] =   0.0;
       p3d[1] =   1.0;
       p3d[2] =   0.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);

       edge = Edge::newObject( legomesh->getNodeAt(0), legomesh->getNodeAt(7));
       legomesh->addObject(edge);

       p3d[0] =  -1.0;
       p3d[1] =   1.0;
       p3d[2] =   0.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);
       edge = Edge::newObject( legomesh->getNodeAt(7), legomesh->getNodeAt(8));
       legomesh->addObject(edge);

       p3d[0] =  -1.0;
       p3d[1] =   0.0;
       p3d[2] =   0.0;
       vtx = Vertex::newObject();
       vtx->setXYZCoords(p3d);
       legomesh->addObject(vtx);
       edge = Edge::newObject( legomesh->getNodeAt(8), legomesh->getNodeAt(9));
       legomesh->addObject(edge);

       edge = Edge::newObject( legomesh->getNodeAt(9), legomesh->getNodeAt(0));
       legomesh->addObject(edge);

       meshViewer->setNewMesh( legomesh );
       return;
    */

    JFaceSequence faces;

    JCellPtr cell;
    cell = JHexahedron::getCanonical(1.0);
    nodes = cell->getNodes();
    faces = cell->getFaces();
    legomesh->addObjects(nodes);
    legomesh->addObjects(faces);
    JCellPtr centercell = cell;
    legomesh->addObject(centercell); // Cell 1

    JFacePtr face;
    JMeshAffineTransform af;
    Point3D avg;

    face = JQuadrilateral::getCanonical(1);
    face->getAvgXYZ(avg);
    nodes = face->getNodes();
    af.scale(nodes, sqrt(2), 1, 1);
    af.rotate(nodes, 45, 1);
    af.translate(nodes, 2.0, 0.0, 0.0);
    face->getAvgXYZ(avg);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);
    legomesh->addObject(cell); // Cell 1

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(1);
    nodes[1] = legomesh->getNodeAt(5);
    nodes[2] = legomesh->getNodeAt(6);
    nodes[3] = legomesh->getNodeAt(2);
    nodes[4] = legomesh->getNodeAt(8);
    nodes[5] = legomesh->getNodeAt(9);
    nodes[6] = legomesh->getNodeAt(10);
    nodes[7] = legomesh->getNodeAt(11);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell); // Cell 2

    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, sqrt(2), 1, 1);
    af.rotate(nodes, -45, 1);
    af.translate(nodes, 2.0, 0.0, -2.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(8);
    nodes[1] = legomesh->getNodeAt(9);
    nodes[2] = legomesh->getNodeAt(10);
    nodes[3] = legomesh->getNodeAt(11);
    nodes[4] = legomesh->getNodeAt(12);
    nodes[5] = legomesh->getNodeAt(13);
    nodes[6] = legomesh->getNodeAt(14);
    nodes[7] = legomesh->getNodeAt(15);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell); // Cell 3

    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, sqrt(2), 1, 1);
    af.rotate(nodes, 45, 1);
    af.translate(nodes, 0.0, 0.0, -2.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(12);
    nodes[1] = legomesh->getNodeAt(13);
    nodes[2] = legomesh->getNodeAt(14);
    nodes[3] = legomesh->getNodeAt(15);
    nodes[4] = legomesh->getNodeAt(17);
    nodes[5] = legomesh->getNodeAt(16);
    nodes[6] = legomesh->getNodeAt(19);
    nodes[7] = legomesh->getNodeAt(18);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell); // Cell 4

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(0);
    nodes[1] = legomesh->getNodeAt(1);
    nodes[2] = legomesh->getNodeAt(2);
    nodes[3] = legomesh->getNodeAt(3);
    nodes[4] = legomesh->getNodeAt(16);
    nodes[5] = legomesh->getNodeAt(17);
    nodes[6] = legomesh->getNodeAt(18);
    nodes[7] = legomesh->getNodeAt(19);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(1);
    nodes[1] = legomesh->getNodeAt(8);
    nodes[2] = legomesh->getNodeAt(11);
    nodes[3] = legomesh->getNodeAt(2);
    nodes[4] = legomesh->getNodeAt(17);
    nodes[5] = legomesh->getNodeAt(12);
    nodes[6] = legomesh->getNodeAt(15);
    nodes[7] = legomesh->getNodeAt(18);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
// legomesh->addObject(cell);
    legomesh->addObject(centercell); // Cell 1

    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, 1, sqrt(2), 1);
    af.rotate(nodes, -45, 0);
    af.translate(nodes, 0.0, 0.0,  2.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(4);
    nodes[1] = legomesh->getNodeAt(5);
    nodes[2] = legomesh->getNodeAt(6);
    nodes[3] = legomesh->getNodeAt(7);
    nodes[4] = legomesh->getNodeAt(20);
    nodes[5] = legomesh->getNodeAt(21);
    nodes[6] = legomesh->getNodeAt(22);
    nodes[7] = legomesh->getNodeAt(23);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, 1,sqrt(2), 1);
    af.rotate(nodes, 45, 0);
    af.translate(nodes, 0.0, -2.0,  2.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(20);
    nodes[1] = legomesh->getNodeAt(21);
    nodes[2] = legomesh->getNodeAt(22);
    nodes[3] = legomesh->getNodeAt(23);
    nodes[4] = legomesh->getNodeAt(27);
    nodes[5] = legomesh->getNodeAt(26);
    nodes[6] = legomesh->getNodeAt(25);
    nodes[7] = legomesh->getNodeAt(24);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, 1, sqrt(2), 1);
    af.rotate(nodes, -45, 0);
    af.translate(nodes, 0.0, -2.0,  0.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(24);
    nodes[1] = legomesh->getNodeAt(25);
    nodes[2] = legomesh->getNodeAt(26);
    nodes[3] = legomesh->getNodeAt(27);
    nodes[4] = legomesh->getNodeAt(28);
    nodes[5] = legomesh->getNodeAt(29);
    nodes[6] = legomesh->getNodeAt(30);
    nodes[7] = legomesh->getNodeAt(31);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(31);
    nodes[1] = legomesh->getNodeAt(30);
    nodes[2] = legomesh->getNodeAt(29);
    nodes[3] = legomesh->getNodeAt(28);
    nodes[4] = legomesh->getNodeAt(4);
    nodes[5] = legomesh->getNodeAt(5);
    nodes[6] = legomesh->getNodeAt(1);
    nodes[7] = legomesh->getNodeAt(0);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(27);
    nodes[1] = legomesh->getNodeAt(31);
    nodes[2] = legomesh->getNodeAt(4);
    nodes[3] = legomesh->getNodeAt(20);
    nodes[4] = legomesh->getNodeAt(26);
    nodes[5] = legomesh->getNodeAt(30);
    nodes[6] = legomesh->getNodeAt(5);
    nodes[7] = legomesh->getNodeAt(21);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
// legomesh->addObject(cell);
    legomesh->addObject(centercell); // Cell 1


    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, 1, sqrt(2), 1);
    af.rotate(nodes, 90, 1);
    af.rotate(nodes, 45, 2);
    af.translate(nodes, 0.0, 2.0,  0.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(2);
    nodes[1] = legomesh->getNodeAt(6);
    nodes[2] = legomesh->getNodeAt(7);
    nodes[3] = legomesh->getNodeAt(3);
    nodes[4] = legomesh->getNodeAt(35);
    nodes[5] = legomesh->getNodeAt(34);
    nodes[6] = legomesh->getNodeAt(33);
    nodes[7] = legomesh->getNodeAt(32);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, 1, sqrt(2), 1);
    af.rotate(nodes, 90, 1);
    af.rotate(nodes, -45, 2);
    af.translate(nodes, -2.0, 2.0,  0.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(32);
    nodes[1] = legomesh->getNodeAt(33);
    nodes[2] = legomesh->getNodeAt(34);
    nodes[3] = legomesh->getNodeAt(35);
    nodes[4] = legomesh->getNodeAt(36);
    nodes[5] = legomesh->getNodeAt(37);
    nodes[6] = legomesh->getNodeAt(38);
    nodes[7] = legomesh->getNodeAt(39);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    face = JQuadrilateral::getCanonical(1.0);
    nodes = face->getNodes();
    af.scale(nodes, 1, sqrt(2), 1);
    af.rotate(nodes, 90, 1);
    af.rotate(nodes, 45, 2);
    af.translate(nodes, -2.0, 0.0,  0.0);
    legomesh->addObjects(nodes);
    legomesh->addObject(face);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(36);
    nodes[1] = legomesh->getNodeAt(37);
    nodes[2] = legomesh->getNodeAt(38);
    nodes[3] = legomesh->getNodeAt(39);
    nodes[4] = legomesh->getNodeAt(43);
    nodes[5] = legomesh->getNodeAt(42);
    nodes[6] = legomesh->getNodeAt(41);
    nodes[7] = legomesh->getNodeAt(40);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(0);
    nodes[1] = legomesh->getNodeAt(4);
    nodes[2] = legomesh->getNodeAt(7);
    nodes[3] = legomesh->getNodeAt(3);
    nodes[4] = legomesh->getNodeAt(40);
    nodes[5] = legomesh->getNodeAt(41);
    nodes[6] = legomesh->getNodeAt(42);
    nodes[7] = legomesh->getNodeAt(43);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
    legomesh->addObject(cell);

    cell = JHexahedron::newObject();
    nodes.resize(8);
    nodes[0] = legomesh->getNodeAt(3);
    nodes[1] = legomesh->getNodeAt(7);
    nodes[2] = legomesh->getNodeAt(33);
    nodes[3] = legomesh->getNodeAt(32);
    nodes[4] = legomesh->getNodeAt(43);
    nodes[5] = legomesh->getNodeAt(42);
    nodes[6] = legomesh->getNodeAt(37);
    nodes[7] = legomesh->getNodeAt(36);
    cell->setNodes(nodes);
    faces = cell->getFaces();
    legomesh->addObjects(faces);
//   legomesh->addObject(cell);

    meshViewer->addObject( legomesh );

}
///////////////////////////////////////////////////////////////////////////////
void JLegoBuilderDialog :: attachAtFace()
{
    if( meshViewer == nullptr ) return;
    JMeshEntityPickerPtr picker = meshViewer->getEntityPicker();
    if( picker == nullptr ) return;
    JFaceSequence facesPicked = picker->getPickedFaces();

    if( facesPicked.empty() ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("No faces picked: Enable picker and select");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

}
///////////////////////////////////////////////////////////////////////////////
void JLegoBuilderDialog :: newCubeNodes()
{
}
///////////////////////////////////////////////////////////////////////////////
void JLegoBuilderDialog :: newCubeIJK()
{
}
///////////////////////////////////////////////////////////////////////////////
void JLegoBuilderDialog :: makeConnections()
{
    connect( initPushButton,  SIGNAL( clicked() ), this, SLOT( init() ));
    connect( closePushButton,  SIGNAL( clicked() ), this, SLOT( close() ));
    connect( newCubeIJKPushButton,  SIGNAL( clicked() ), this, SLOT( newCubeIJK() ));
    connect( newCubeNodesPushButton,  SIGNAL( clicked() ), this, SLOT( newCubeNodes() ));
    connect( pickedFacePushButton,  SIGNAL( clicked() ), this, SLOT( attachAtFace() ));
}
///////////////////////////////////////////////////////////////////////////////
