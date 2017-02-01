#include "MeshRefine2DDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshRefine2DDialog :: JMeshRefine2DDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr) return;

    entityPicker = meshViewer->getEntityPicker();
    setMesh( meshViewer->getCurrentMesh() );

    viewManager->attach( this );
    entityPicker->setMode(2);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 2;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );
    entityPicker->setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: mouseReleaseEvent(QMouseEvent *)
{
    if( entityPicker == nullptr || mesh == nullptr ) return;
    JFaceSequence faceSeq = entityPicker->getPickedFaces();
    if( faceSeq.empty() ) return;
    for( const JFacePtr &f : faceSeq )
        selectedFaces.insert(f);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: refinetrimesh()
{
    if( mesh == nullptr ) return;

    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The topological dimension of the mesh must be two");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JFaceSequence faces2refine;
    faces2refine = entityPicker->getPickedFaces();

    JWaitCursor waitCursor;
    waitCursor.start();

    JTriRefiner triRefiner;
    triRefiner.setMesh(mesh);
    if( faces2refine.empty() )
        triRefiner.refineAll(refine_tri_type);
    else
        triRefiner.refine(faces2refine, refine_tri_type);

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: refinetri13()
{
    refine_tri_type = 13;
    refinetrimesh();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: refinetri14()
{
    refine_tri_type = 14;
    refinetrimesh();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: refinetri16()
{
    refine_tri_type = 16;
    refinetrimesh();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: tri2quads()
{

}
///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: refineQuad14()
{
    if( meshViewer == nullptr ) return;

    string name = StdString(objectNameLineEdit->text());

    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The topological dimension of the mesh must be two");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    JQuadRefiner quadRefiner;
    quadRefiner.setMesh(mesh);

    if( !selectedFaces.empty() ) {
        JFaceSequence faces;
        boost::copy(selectedFaces, back_inserter(faces));
        quadRefiner.refineAll(faces, 14);
        selectedFaces.clear();
        entityPicker->clearAll();
    } else
        quadRefiner.refineAll(14);

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();
    mesh->enumerate(0);
    mesh->enumerate(1);
    mesh->enumerate(2);

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: refineQuad15()
{
    if( meshViewer == nullptr ) return;

    string name = StdString(objectNameLineEdit->text());

    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The topological dimension of the mesh must be two");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    JQuadRefiner quadRefiner;
    quadRefiner.setMesh(mesh);

    if( !selectedFaces.empty() ) {
        JFaceSequence faces;
        boost::copy(selectedFaces, back_inserter(faces));
        quadRefiner.refineAll(faces, 15);
        selectedFaces.clear();
        entityPicker->clearAll();
    } else
        quadRefiner.refineAll(15);

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: quadtri2()
{
    if( meshViewer == nullptr ) return;

    string name = StdString(objectNameLineEdit->text());

    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The topological dimension of the mesh must be two");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    bool randomize = randomQuadDiagonalCheckBox->isChecked();

    JMeshPtr trimesh = AllTriMeshGenerator::getFromQuadMesh( mesh, 2, randomize);
    meshViewer->removeObject(mesh);
    meshViewer->addObject(trimesh);
    mesh = trimesh;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRefine2DDialog :: quadtri4()
{
    if( meshViewer == nullptr ) return;

    string name = StdString(objectNameLineEdit->text());

    int topDim = mesh->getTopology()->getDimension();
    if( topDim != 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("The topological dimension of the mesh must be two");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    JMeshPtr trimesh = AllTriMeshGenerator::getFromQuadMesh( mesh, 4);
    trimesh->setName( mesh->getName() );
    meshViewer->removeObject(mesh);
    meshViewer->addObject(trimesh);
    mesh = trimesh;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRefine2DDialog :: insertPillows()
{
    if( meshViewer == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    string name = StdString(objectNameLineEdit->text());

    JQuadRefiner qrefiner;
    qrefiner.setMesh(mesh);
    qrefiner.insert_boundary_pillows();

    JLloydMeshOptimizer mopt;
    mopt.setMesh(mesh);
    mopt.setNumIterations(10);
    mopt.smoothAll();

    mesh->getGeometry()->setFacesNormal();
    mesh->getGeometry()->setNodesNormal();

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRefine2DDialog :: quadBlocks()
{
    if( meshViewer == nullptr ) return;

    string name = StdString(objectNameLineEdit->text());

    int dim[2];
    QString qstr;
    qstr = quadNLineEdit->text();
    dim[0] = qstr.toInt();

    qstr = quadMLineEdit->text();
    dim[1] = qstr.toInt();

    JQuadRefiner quadRefiner;
    quadRefiner.setMesh(mesh);
    quadRefiner.refineAll(dim);

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: refineEdge()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: flipEdge()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: setRefineSet()
{
    if( meshViewer == nullptr ) return;
    /*
        JMeshEntityPickerPtr entityPicker = meshViewer->getEntityPicker();
        if( entityPicker == nullptr ) return;

        if( apply2AllCheckBox->isChecked() ) {
            entityPicker->setMode(-1);
            entityPicker->clearAll();
        } else {
            entityPicker->setMode(2);
            entityPicker->setPickableEntity(2);
        }
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: openTri2QuadsDialog()
{
    if( quadmesherDialog == nullptr)
        quadmesherDialog.reset( new JQuadMesherDialog(this));

    meshViewer->setCurrentMesh(mesh);
    quadmesherDialog->setViewManager(viewManager);
    quadmesherDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: openMeshSubdivisionDialog()
{
    if( meshSubdivisionDialog == nullptr)
        meshSubdivisionDialog.reset( new JMeshSubdivisionDialog(this));

    meshSubdivisionDialog->setViewManager(viewManager);
    meshSubdivisionDialog->setMesh(mesh);
    meshSubdivisionDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshRefine2DDialog :: closeDialog()
{
    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = -1;
    }
    entityPicker->clearAll();

    viewManager->detach( this );

    parentWidget()->show();
    close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshRefine2DDialog :: makeConnections()
{
    PushButton( tri13PushButton,  [=] {refinetri13();} );
    PushButton( tri14PushButton,  [=] {refinetri14();});
    PushButton( tri16PushButton,  [=] {refinetri16();});

    PushButton( quad14PushButton, [=] {refineQuad14();});
    PushButton( quad15PushButton, [=] {refineQuad15();});

    PushButton( quadTri2PushButton, [=] {quadtri2();});
    PushButton( quadTri4PushButton, [=] {quadtri4();});
    PushButton( quadBlocksPushButton, [=] {quadBlocks(); });

    PushButton( pillowPushButton,   [=] {insertPillows(); });

    PushButton( subdivisionPushButton, [=] {openMeshSubdivisionDialog(); });
    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////

void JMeshRefine2DDialog :: projectOnCircle()
{
    /*
        if( meshViewer == nullptr ) return;
        JMeshPtr mesh = meshViewer->getMesh();

        if( mesh == nullptr ) return;

        mesh->getTopology()->search_boundary();

        Point3D pCenter = mesh->getGeometry()->getCenter();
        size_t nSize = mesh->getSize(0);

        Vec3D vec;
        Point3D xyz;
        double rad = 0.0, absval;

        for( size_t i = 0; i < nSize; i++) {
            JNodePtr vtx = mesh->getNodeAt(i);
            if( vtx->isBoundary() ) {
                xyz = vtx->getXYZCoords();
                vec[0] =  xyz[0] - pCenter[0];
                vec[1] =  xyz[1] - pCenter[1];
                vec[2] =  xyz[2] - pCenter[2];
                rad    =  max(rad, vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
            }
        }

        rad = sqrt( rad );

        if( rad <= 0.0) return;

        for( size_t i = 0; i < nSize; i++) {
            JNodePtr vtx = mesh->getNodeAt(i);
            if( vtx->isBoundary() ) {
                xyz = vtx->getXYZCoords();
                vec[0]  =  xyz[0] - pCenter[0];
                vec[1]  =  xyz[1] - pCenter[1];
                vec[2]  =  xyz[2] - pCenter[2];
                absval  =  1.0/sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
                vec[0]  *= absval;
                vec[1]  *= absval;
                vec[2]  *= absval;

                xyz[0]  =  pCenter[0] + rad*vec[0];
                xyz[1]  =  pCenter[1] + rad*vec[1];
                xyz[2]  =  pCenter[2] + rad*vec[2];
                vtx->setXYZCoords(xyz);
            }
        }
    */
}
