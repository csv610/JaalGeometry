#include "InteractiveMeshingDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JInteractiveMeshingDialog :: JInteractiveMeshingDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();

    oldmesh   = nullptr;
    mesh      = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JInteractiveMeshingDialog :: ~JInteractiveMeshingDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JInteractiveMeshingDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}
///////////////////////////////////////////////////////////////////////////////

void JInteractiveMeshingDialog :: keyPressEvent( QKeyEvent *e)
{
    if( meshViewer == nullptr ) return;

    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JInteractiveMeshingDialog :: setNode()
{
    if( oldmesh != nullptr && appendMeshCheckBox->isChecked() )
        mesh = oldmesh;

    if( mesh == nullptr ) {
        mesh = JMesh::newObject();
        meshViewer->addObject(mesh);
    }

    bool embedded = embedCheckBox->isChecked();

    QString qstr;
    Point3D p3d;
    qstr = xCoordLineEdit->text();
    p3d[0]   = qstr.toDouble();

    qstr = yCoordLineEdit->text();
    p3d[1]   = qstr.toDouble();

    qstr = zCoordLineEdit->text();
    p3d[2]   = qstr.toDouble();

    JNodePtr newnode = JNode::newObject();
    newnode->setXYZCoords(p3d);
    bool val = 1;
    newnode->setAttribute("Display", val);
    mesh->addObject(newnode);
}
///////////////////////////////////////////////////////////////////////////////
void JInteractiveMeshingDialog :: genNewEdge()
{
    if( oldmesh != nullptr && appendMeshCheckBox->isChecked() )
        mesh = oldmesh;

    if( mesh == nullptr ) {
        mesh = JMesh::newObject();
        meshViewer->addObject(mesh);
    }

    if( nodelist.size() != 2 ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("New edge require atleast two nodes");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    JNodePtr v0 = mesh->getNodeAt( nodelist[0] );
    JNodePtr v1 = mesh->getNodeAt( nodelist[1] );
    JEdgePtr edge = JEdge::newObject(v0,v1);
    mesh->addObject(edge);
    bool val = 1;
    edge->setAttribute("Display", val);
}
///////////////////////////////////////////////////////////////////////////////

void JInteractiveMeshingDialog :: genNewFace()
{
    if( oldmesh != nullptr && appendMeshCheckBox->isChecked() )
        mesh = oldmesh;


    if( mesh == nullptr ) {
        mesh = JMesh::newObject();
        meshViewer->addObject(mesh);
    }

    if( nodelist.size() < 3 ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("New edge require atleast three nodes");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
    nodeseq.clear();
    int nsize = nodelist.size();
    for( int i = 0; i < nsize; i++) {
        JNodePtr vtx = mesh->getNodeAt(nodelist[i] );
        if(vtx ) nodeseq.push_back(vtx);
    }

    if( (int)nodeseq.size() == nsize ) {
        JFacePtr face = JFace::newObject(nodeseq);
        mesh->addObject(face);
        bool val = 1;
        face->setAttribute("Display", val);
    }
}
///////////////////////////////////////////////////////////////////////////////

void JInteractiveMeshingDialog :: genNewCell()
{
    if( oldmesh != nullptr && appendMeshCheckBox->isChecked() )
        mesh = oldmesh;

    if( mesh == nullptr ) {
        mesh = JMesh::newObject();
        meshViewer->addObject(mesh);
    }

    if( nodelist.size() < 4 ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("New edge require atleast four nodes");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if(ret == QMessageBox::Ok ) return;
    }
    nodeseq.clear();

    int nsize = nodelist.size();
    for( int i = 0; i < nsize; i++) {
        JNodePtr vtx = mesh->getNodeAt(nodelist[i] );
        if(vtx ) nodeseq.push_back(vtx);
    }

    if( (int)nodeseq.size() == nsize ) {
        JCellPtr cell = JCell::newObject(nodeseq);
        mesh->addObject(cell);
        bool val = 1;
        cell->setAttribute("Display", val);
    }
}
///////////////////////////////////////////////////////////////////////////////

void JInteractiveMeshingDialog :: setEntity()
{
    QString qstr;
    qstr = nodeseqLineEdit->text();
    string str = qstr.toUtf8().constData();
    nodelist.clear();

    StringTokenizer  token(str);
    StringTokenizer::iterator iter;
    for(iter = token.begin(); iter != token.end(); ++iter) {
        string s = *iter;
        nodelist.push_back(atoi(s.c_str()));
    }
    if( edgeRadioButton->isChecked() ) genNewEdge();
    if( faceRadioButton->isChecked() ) genNewFace();
    if( cellRadioButton->isChecked() ) genNewCell();
}

///////////////////////////////////////////////////////////////////////////////

void JInteractiveMeshingDialog :: makeConnections()
{
    connect( applyNodePushButton,  SIGNAL( clicked() ), this, SLOT( setNode() ));
    connect( applyEntityPushButton,  SIGNAL( clicked() ), this, SLOT( setEntity() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
