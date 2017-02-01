#include "MeshRenderDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshRenderDialog :: JMeshRenderDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);

    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshRenderDialog :: ~JMeshRenderDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: showEvent( QShowEvent *e)
{
    if( meshViewer == nullptr ) return;

/*
    mesh = meshViewer->getCurrentMesh();
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );
*/
    QDialog::showEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: warnMessage()
{
    QMessageBox msg;
    msg.setIcon(QMessageBox::Warning);
    msg.setText("No mesh object selected: Select one from the list");
    msg.setStandardButtons( QMessageBox::Ok);
    int ret = msg.exec();
    if( ret == QMessageBox::Ok ) return;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->getViewManager()->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: openMeshlistDialog()
{
    if( meshlistDialog.get() == nullptr )
        meshlistDialog.reset(new JObjectsListDialog(this));

    meshlistDialog->setViewManager( viewManager );
    meshlistDialog->setType(1);
    meshlistDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: openNodesDialog()
{
    if( meshViewer == nullptr ) return;

    string name = StdString(objectNameLineEdit->text());
    mesh = meshViewer->getMesh(name);
    if( mesh == nullptr) {
        warnMessage();
        return;
    }

    if( meshNodesDialog.get() == nullptr )
        meshNodesDialog.reset(new JMeshNodesDialog(this));

    meshNodesDialog->setViewManager( viewManager );
    meshNodesDialog->setMesh( mesh );

    meshNodesDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: openEdgesDialog()
{
    if( meshViewer == nullptr) return;

    string name = StdString(objectNameLineEdit->text());
    mesh = meshViewer->getMesh(name);

    if( mesh == nullptr) {
        warnMessage();
        return;
    }

    if( meshEdgesDialog.get() == nullptr )
        meshEdgesDialog.reset(new JMeshEdgesDialog(this));

    meshEdgesDialog->setViewManager( viewManager );
    meshEdgesDialog->setMesh( mesh );

    meshEdgesDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: openFacesDialog()
{
    if( meshViewer == nullptr) return;

    string name = StdString(objectNameLineEdit->text());
    mesh = meshViewer->getMesh(name);

    if( mesh == nullptr) {
        warnMessage();
        return;
    }

    if( meshFacesDialog.get() == nullptr)
        meshFacesDialog.reset(new JMeshFacesDialog(this));

    meshFacesDialog->setViewManager( viewManager );
    meshFacesDialog->setMesh( mesh );

    meshFacesDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: openCellsDialog()
{
    if( meshViewer == nullptr) return;

    string name = StdString(objectNameLineEdit->text());
    mesh = meshViewer->getMesh(name);

    if( mesh == nullptr) {
        warnMessage();
        return;
    }

    if( meshCellsDialog.get() == nullptr )
        meshCellsDialog.reset(new JMeshCellsDialog(this));

    meshCellsDialog->setViewManager( viewManager );
    meshCellsDialog->setMesh( mesh );

    meshCellsDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////


void JMeshRenderDialog :: closeDialog()
{
    this->close();
    parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshRenderDialog :: openMaterialDialog()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );
}
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshRenderDialog :: selectMesh( QModelIndex index)
{
    int  row = index.row();
    mesh = meshViewer->getMesh(row);
    JBoundingBox box;
    mesh->getAttribute("AxisBoundingBox", box);
    viewManager->setCenter(box.getCenter());
}
*/

///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: resetCamera()
{
    if( mesh == nullptr) return;
    meshViewer->lookAt(mesh);

/*
    viewManager->resetView( JViewDirection:: FRONT_VIEW);

    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    Point3D  pC   = box.getCenter();
    double   len  = box.getMaxLength();

    qglviewer::Vec  vec;
    vec[0] = pC[0];
    vec[1] = pC[1];
    vec[2] = pC[2] + 2.0*len;
    viewManager->camera()->setPosition(vec);

    vec[0] = pC[0];
    vec[1] = pC[1];
    vec[2] = pC[2];
    viewManager->camera()->lookAt(vec);
    viewManager->camera()->setSceneCenter(vec);

    vec[0] = 0.0;
    vec[1] = 1;
    vec[2] = 0;
    viewManager->camera()->setUpVector(vec);
    meshViewer->refreshDisplay();
*/
}

///////////////////////////////////////////////////////////////////////////////

    /*
void JMeshRenderDialog :: ambientOcclusion()
{
    if( mesh == nullptr) return;


        JWaitCursor waitCursor;
        waitCursor.start();

        JMeshEigenMatrix mat;
        mat.setMesh(mesh);
        Eigen::MatrixXd V = mat.getNodeMatrix();
        Eigen::MatrixXi F = mat.getFaceMatrix();

        MatrixXd N;
        igl::per_vertex_normals(V,F,N);

        Eigen::VectorXd AO;
        // Compute ambient occlusion factor using embree
        igl::embree::ambient_occlusion(V,F,V,N,500,AO);
        AO = 1.0 - AO.array();

        const RowVector3d color(0.9,0.85,0.9);
        MatrixXd C = color.replicate(V.rows(),1);
        for (unsigned i=0; i<C.rows(); ++i) {
            C.row(i) *= AO(i); //std::min<double>(AO(i)+0.2,1);
        }

        JNodeRenderPtr vAttrib;
        size_t numnodes = mesh->getSize(0);
        JColor rgb;
        size_t index = 0;
        for( size_t i = 0; i <  numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->isActive() ) {
                vtx->getAttribute("Render", vAttrib);
                rgb[0] = C.coeff(index,0);
                rgb[1] = C.coeff(index,1);
                rgb[2] = C.coeff(index,2);
                vAttrib->color = rgb;
                index++;
            }
        }
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->setSurfaceShade(JRender::SMOOTH_SHADE);
        meshViewer->updateBuffers(mesh);
}
    */

///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: fitBoundSphere()
{
    if( mesh == nullptr) return;

    JSphere sph = mesh->getGeometry()->getMinSphere();
    Point3D center = sph.getCenter();

    qglviewer::Vec c;

    c[0] = center[0];
    c[1] = center[1];
    c[2] = center[2];
    viewManager->camera()->fitSphere(c,sph.getRadius());
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: fitBoundBox()
{
    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    Point3D pmin = box.getCorner(0);
    Point3D pmax = box.getCorner(6);

    qglviewer::Vec vmin, vmax;
    vmin[0] =  pmin[0];
    vmin[1] =  pmin[1];
    vmin[2] =  pmin[2];

    vmax[0] =  pmax[0];
    vmax[1] =  pmax[1];
    vmax[2] =  pmax[2];
    viewManager->camera()->fitBoundingBox(vmin,vmax);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: openBoxColorDialog()
{
    QColor color = QColorDialog::getColor();
    float rgb[3];
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;

    /*
        JEdgeRenderPtr eAttrib;
        for( size_t i = 0; i < edges.size(); i++) {
            edges[i]->getAttribute("Render", eAttrib);
            eAttrib->color[0] = rgb[0];
            eAttrib->color[1] = rgb[1];
            eAttrib->color[2] = rgb[2];
        }
        meshViewer->updateBuffers( mesh );
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: setEnclosure()
{
}

void JMeshRenderDialog :: setPOVRayScene()
{
   JMeshPOVExporter mexp;

   qglviewer::Vec pos;
   Point3D xyz;

   pos = viewManager->camera()->position();
   xyz[0] = pos[0];
   xyz[1] = pos[1];
   xyz[2] = pos[2];
   mexp.setCamera(xyz);

   pos = viewManager->camera()->viewDirection();
   xyz[0] = pos[0];
   xyz[1] = pos[1];
   xyz[2] = pos[2];
   mexp.setLookAt(xyz);

   Point4F lpos;
   glGetLightfv( GL_LIGHT0, GL_POSITION, &lpos[0] );
   pos[0] = lpos[0];
   pos[1] = lpos[1];
   pos[2] = lpos[2];
   mexp.addLight(xyz);

   mexp.writeFile(mesh, "model.pov");
   if( povRayImageCheckBox->isChecked() ) {
  
       ostringstream oss;
       oss << "povray ";
       oss << "+H" << viewManager->height() << " ";
       oss << "+W" << (int)(4.0*viewManager->height()/3.0)  << " ";
       oss << " model.pov";
       system( oss.str().c_str() );
       system( "display model.png &");
   }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshRenderDialog :: makeConnections()
{
    PushButton( resetCameraPushButton,[=] {resetCamera();});
    PushButton( boxColorPushButton,   [=] {openBoxColorDialog();});
    PushButton( meshNodesPushButton,  [=] {openNodesDialog();});
    PushButton( meshEdgesPushButton,  [=] {openEdgesDialog();});
    PushButton( meshFacesPushButton,  [=] {openFacesDialog();});
    PushButton( meshCellsPushButton,  [=] {openCellsDialog();});
    PushButton( povrayPushButton,     [=] {setPOVRayScene();});
//  PushButton( fitSpherePushButton, [=] {fitBoundSphere();});

    ComboBox( boundingComboBox, [=] { setEnclosure(); });
    CheckBox( boundingCheckBox,    [=] {setEnclosure();});

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////


