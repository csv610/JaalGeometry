#include "MeshOptBoundaryLayerDialog.hpp"

void JContourNormalsViewer :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    Vec3D normal;
    Point3D pmid;

    eNormalsTail.clear();
    eNormalsHead.clear();

    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i  < numedges; i++) {
        const JEdgePtr &edge  = mesh->getEdgeAt(i);
        int err = edge->getAttribute("ContourNormal", normal);
        if( !err ) {
            pmid = JEdgeGeometry::getMidPoint(edge);
            eNormalsTail.push_back(pmid[0]);
            eNormalsTail.push_back(pmid[1]);
            eNormalsTail.push_back(pmid[2]);
            eNormalsHead.push_back(pmid[0] + scale*normal[0]);
            eNormalsHead.push_back(pmid[1] + scale*normal[1]);
            eNormalsHead.push_back(pmid[2] + scale*normal[2]);
        }
    }

    vNormalsTail.clear();
    vNormalsHead.clear();
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i  < numnodes; i++) {
        const JNodePtr &vtx  = mesh->getNodeAt(i);
        int err = vtx->getAttribute("ContourNormal", normal);
        if( !err ) {
            pmid = vtx->getXYZCoords();
            vNormalsTail.push_back(pmid[0]);
            vNormalsTail.push_back(pmid[1]);
            vNormalsTail.push_back(pmid[2]);
            vNormalsHead.push_back(pmid[0] + scale*normal[0]);
            vNormalsHead.push_back(pmid[1] + scale*normal[1]);
            vNormalsHead.push_back(pmid[2] + scale*normal[2]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////

void JContourNormalsViewer :: draw()
{
    glLineWidth(1.0);
    glColor3f( 1.0, 0.2, 0.2);

    if( displayEdgeNormals) {
        size_t numedges = eNormalsTail.size()/3;
        glBegin(GL_LINES);
        for( size_t i = 0; i < numedges; i++) {
            glVertex3fv( &eNormalsTail[3*i] ) ;
            glVertex3fv( &eNormalsHead[3*i] ) ;
        }
        glEnd();

        glPointSize(5);
        glColor3f( 0.0, 1.0, 0.0);

        glBegin(GL_POINTS);
        for( size_t i = 0; i < numedges; i++) {
            glVertex3fv( &eNormalsTail[3*i] ) ;
        }
        glEnd();
        glColor3f( 0.0, 0.0, 1.0);
        glBegin(GL_POINTS);
        for( size_t i = 0; i < numedges; i++) {
            glVertex3fv( &eNormalsHead[3*i] ) ;
        }
        glEnd();
        glPointSize(1);
    }

    if( displayNodeNormals) {
        size_t numnodes = vNormalsTail.size()/3;
        glBegin(GL_LINES);
        for( size_t i = 0; i < numnodes; i++) {
            glVertex3fv( &vNormalsTail[3*i] ) ;
            glVertex3fv( &vNormalsHead[3*i] ) ;
        }
        glEnd();
    }
}
///////////////////////////////////////////////////////////////////////////////

JMeshOptBoundaryLayerDialog :: JMeshOptBoundaryLayerDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshOptBoundaryLayerDialog :: ~JMeshOptBoundaryLayerDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshOptBoundaryLayerDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshOptBoundaryLayerDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    double elen = mesh->getGeometry()->getMeanEdgeLength();
    boundaryOffsetLineEdit->setText( QString::number(0.90*elen));


}

///////////////////////////////////////////////////////////////////////////////////
void JMeshOptBoundaryLayerDialog :: displayNormals()
{
    if( contourNormalsViewer == nullptr) {
        contourNormalsViewer.reset ( new JContourNormalsViewer );
        contourNormalsViewer->setName("ContourNormalsViewer");
    }

    contourNormalsViewer->setMesh(mesh);

    int  val1 = showEdgeNormalsCheckBox->isChecked();
    int  val2 = showNodeNormalsCheckBox->isChecked();

    if( val1 + val2 > 0)
        viewManager->attach( contourNormalsViewer);
    else
        viewManager->detach( contourNormalsViewer);

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////////

void JMeshOptBoundaryLayerDialog :: optLayer()
{
    JMeshOptBoundaryLayer  mopt;
    mopt.setMesh(mesh);

    QString str = boundaryOffsetLineEdit->text();
    double  d   = str.toDouble();
    mopt.setOffset(d);
    JWaitCursor waitCursor;
    waitCursor.start();

    mopt.setNormals();

    displayNormals();

    mopt.update();
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshOptBoundaryLayerDialog :: makeConnections()
{
    CheckBox( showEdgeNormalsCheckBox, [=] {displayNormals();});
    CheckBox( showNodeNormalsCheckBox, [=] { displayNormals();});

    PushButton( optLayerPushButton, [=] {optLayer();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////
