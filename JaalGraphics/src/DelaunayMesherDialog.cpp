#include "DelaunayMesherDialog.hpp"

int JDelaunayEdgeColor :: assign( const JEdgePtr &edge)
{
    if( edge == nullptr ) return 1;
    if( !edge->isActive() )  return 2;

    JEdgeRenderPtr eAttrib;
    edge->getAttribute("Render", eAttrib);

    if( JDelaunayMesh2D::isDelaunay(edge) ) {
        eAttrib->color = defaultColor;
        eAttrib->scale = 1.0;
    } else {
        eAttrib->color = highlightColor;
        eAttrib->scale = 1.5;
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int  JDelaunayFaceColor :: assign(const JFacePtr &face)
{
    highlightColor[0]  = 0.0;
    highlightColor[1]  = 1.0;
    highlightColor[2]  = 0.0;
    highlightColor[3]  = alpha;

    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);
    fAttrib->color = highlightColor;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JDelaunayMesherDialog :: JDelaunayMesherDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JDelaunayMesherDialog :: ~JDelaunayMesherDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JDelaunayMesherDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JDelaunayMesherDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
}
///////////////////////////////////////////////////////////////////////////////

void JDelaunayMesherDialog :: openTopoQualityDialog()
{
    if( topoQualityDialog.get() == nullptr )
        topoQualityDialog.reset(new JMeshTopologyQualityDialog(this));

    topoQualityDialog->setViewManager( viewManager );
    topoQualityDialog->setMesh(mesh);
    topoQualityDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JDelaunayMesherDialog :: openGeomQualityDialog()
{
    if( geomQualityDialog.get() == nullptr )
        geomQualityDialog.reset(new JMeshGeometricQualityDialog(this));

    geomQualityDialog->setViewManager( viewManager );
    geomQualityDialog->show();
    geomQualityDialog->setMesh( newMesh );
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JDelaunayMesherDialog :: generate()
{
    if( mesh == nullptr ) return;

    int topDim = mesh->getTopology()->getDimension();

    if( topDim > 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning Triangle mesh generation is for 2-manifold only");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    mesher.reset( new JDelaunayMesh2D);
    mesher->setMesh(mesh);

    QString qstr;
    if( opt_a_checkBox->isChecked() ) {
        qstr = maxAreaLineEdit->text();
        mesher->setMaxArea(qstr.toDouble());
    }

    if( newMesh ) {
        meshViewer->removeObject(newMesh);
        newMesh->deleteSurfaceMesh();
    }

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;

    newMesh = mesher->getQualityMesh();
    meshViewer->addObject( newMesh );
}
////////////////////////////////////////////////////////////////////////////////////
void JDelaunayMesherDialog :: closeDialog()
{
   this->close();
   parentWidget()->show();
}
////////////////////////////////////////////////////////////////////////////////////

void JDelaunayMesherDialog :: makeConnections()
{
    PushButton( geomQualityPushButton, [=] {openGeomQualityDialog();});
    PushButton( topoQualityPushButton, [=] {openTopoQualityDialog();});
    PushButton( applyPushButton, [=] {generate();});
    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
