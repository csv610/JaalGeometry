#include "MeshSlicerDialog.hpp"

JMeshSlicerDialog :: JMeshSlicerDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    distanceLineEdit->setText( QString::number(0.0));
}

///////////////////////////////////////////////////////////////////////////////

JMeshSlicerDialog :: ~JMeshSlicerDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSlicerDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}


///////////////////////////////////////////////////////////////////////////////
void JMeshSlicerDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    meshViewer->addObject( plane );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSlicerDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    slicer.setMesh(mesh);
    if( plane == nullptr) {
        plane = JMesh::newObject();
        meshViewer->addObject(plane);
        JMeshRenderPtr mrender;
        plane->getAttribute("Render", mrender);
        mrender->setTransparency(1);
        mrender->setAlpha(0.8);
    }

    // Create  a mesh with only one quadrilateral.
    JBoundingBox box = mesh->getGeometry()->getBoundingBox();
    double len = box.getMaxLength();
    meshCenter = box.getCenter();
    JFacePtr face = JQuadrilateral::getCanonical(meshCenter, len);

    JCellPtr minBox =  mesh->getGeometry()->getMinBox();
    boundBox = JMesh::newObject();
    boundBox->addObjects( minBox->getNodes() );
    boundBox->addObjects( minBox->getEdges() );
    meshViewer->addObject( boundBox );

    plane->clearAll();
    plane->addObjects( face->getNodes() );
    plane->addObject( face );
    Vec3F  n;
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 1.0;
    face->setAttribute("Normal", n);
    meshViewer->updateBuffers(plane);

    JColor red = JEntityColor::getColor("Red");
    JEdgeRenderPtr eAttrib;
    for( int i = 0; i < 4; i++) {
        const JEdgePtr &e = face->getEdgeAt(i);
        e->getAttribute("Render", eAttrib);
        eAttrib->lineWidth = 3.0;
        eAttrib->color = red;
    }
    setPlane();
    contourMesh = JMesh::newObject();
    meshViewer->addObject(contourMesh);
    meshViewer->setCurrentMesh(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSlicerDialog :: setDisplay()
{
    JMeshRenderPtr mrender;
    bool val;

    val = displayPlaneCheckBox->isChecked();
    if( plane ) {
        plane->getAttribute("Render", mrender);
        mrender->displayEntity[0] = val;
        mrender->displayEntity[1] = val;
        mrender->displayEntity[2] = val;
    }

    val = displayBoxCheckBox->isChecked();
    if( boundBox ) {
        boundBox->getAttribute("Render", mrender);
        mrender->displayEntity[0] = val;
        mrender->displayEntity[1] = val;
        mrender->displayEntity[2] = val;
    }

    val = displayMeshCheckBox->isChecked();
    if( mesh ) {
        mesh->getAttribute("Render", mrender);
        mrender->displayEntity[0] = val;
        mrender->displayEntity[1] = val;
        mrender->displayEntity[2] = val;
    }

    val = displayContoursCheckBox->isChecked();
    if( contourMesh ) {
        contourMesh->getAttribute("Render", mrender);
        mrender->displayEntity[0] = val;
        mrender->displayEntity[1] = val;
        mrender->displayEntity[2] = val;
    }
    meshViewer->refreshDisplay();

}
///////////////////////////////////////////////////////////////////////////////
void JMeshSlicerDialog :: setPlane()
{
    if( plane == nullptr) return;

    JMeshAffineTransform affine;
    affine.setMesh(plane);
    affine.toCenter();

    double nx = nxSpinBox->value();
    double ny = nySpinBox->value();
    double nz = nzSpinBox->value();

    QString qstr = distanceLineEdit->text();
    double d = qstr.toDouble();

    Vec3F srcVec, dstVec;
    JFacePtr thisface = plane->getFaceAt(0);
    int err = thisface->getAttribute("Normal", srcVec);

    double absval = sqrt( nx*nx + ny*ny + nz*nz);
    if( absval < 1.0E-06) return;

    dstVec[0] = nx/absval;
    dstVec[1] = ny/absval;
    dstVec[2] = nz/absval;
    meshViewer->alignAlong(plane, srcVec, dstVec);

    Point3D center;
    center[0] = meshCenter[0] + d*dstVec[0];
    center[1] = meshCenter[1] + d*dstVec[1];
    center[2] = meshCenter[2] + d*dstVec[2];
    affine.translate( center[0], center[1], center[2]);
    plane->getGeometry()->setFacesNormal();

    Point3D centroid;
    JFaceGeometry::getCentroid( plane->getFaceAt(0), centroid);
    meshViewer->getViewManager()->setCenter(centroid);

    meshViewer->updateBuffers(plane);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSlicerDialog :: sliceMesh()
{
    if( plane == nullptr) return;

    if( !appendCheckBox->isChecked() )
        contourMesh->deleteAll();

    JWaitCursor waitCursor;
    waitCursor.start();

    JFacePtr pface = plane->getFaceAt(0);
    Vec3F n;
    pface->getAttribute("Normal", n);

    Point3D centroid;
    JFaceGeometry::getCentroid(pface, centroid);

    Vec3D normal;
    normal[0] = n[0];
    normal[1] = n[1];
    normal[2] = n[2];
    vector<JEdgeSequence> cedges = slicer.getContours(normal, centroid);

    int numloops = cedges.size();
    JNodeSequence cnodes;
    for( int i = 0; i < numloops; i++) {
        JMeshTopology::getEntitySet( cedges[i], cnodes);
        contourMesh->addObjects(cnodes);
        contourMesh->addObjects(cedges[i]);
    }
    meshViewer->updateBuffers(contourMesh);

//    JColor color = JEntityColor::getColor("Black");
    JEdgeRenderPtr eAttrib;
    for( int i = 0; i < numloops; i++) {
        for( const JEdgePtr &edge : cedges[i] )  {
            edge->getAttribute("Render", eAttrib);
            eAttrib->lineWidth = 3.0;
//          eAttrib->color     = color;
        }
    }
    meshViewer->updateBuffers(contourMesh);

}
///////////////////////////////////////////////////////////////////////////////
void JMeshSlicerDialog :: clearAll()
{
    if(contourMesh) {
        contourMesh->deleteAll();
        meshViewer->updateBuffers(contourMesh);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSlicerDialog :: openEdgeAttribDialog()
{
    if( contourMesh == nullptr) return;

    if( edgeAttribDialog.get() == nullptr )
        edgeAttribDialog.reset(new JEdgeAttributesDialog( this ));

    edgeAttribDialog->setViewManager(viewManager);
    edgeAttribDialog->setMesh(contourMesh);

    JEdgeSequence eseq;
    size_t numedges = contourMesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &e = contourMesh->getEdgeAt(i);
        if( e->isActive() && !e->isBoundary()  ) eseq.push_back(e);
    }
    edgeAttribDialog->setEdges(eseq);

    edgeAttribDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshSlicerDialog :: closeDialog()
{
    if( plane ) {
        meshViewer->removeObject(plane);
        plane.reset();
    }
    this->close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSlicerDialog :: makeConnections()
{
    SpinBoxd( nxSpinBox, [=] {setPlane();});
    SpinBoxd( nySpinBox, [=] {setPlane();});
    SpinBoxd( nzSpinBox, [=] {setPlane();});

    LineEdit( distanceLineEdit,  [=] {setPlane();});

    CheckBox(displayPlaneCheckBox, [=] {setDisplay();});
    CheckBox(displayBoxCheckBox,   [=] {setDisplay();});
    CheckBox(displayMeshCheckBox,  [=] {setDisplay();});
    CheckBox(displayContoursCheckBox,  [=] {setDisplay();});

    PushButton(contourEdgesPushButton, [=] {openEdgeAttribDialog();});
    PushButton(applyPushButton, [=] {sliceMesh();});
    PushButton(clearAllPushButton, [=] {clearAll();});
    PushButton(closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
