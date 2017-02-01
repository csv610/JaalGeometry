#include "TriMesherDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
GLUquadricObj* JTriDelaunayViewer :: diskObj   = gluNewQuadric();

void JTriDelaunayViewer :: setCircles()
{
    circumcircles.clear();
    if( mesh == nullptr ) return;

    size_t numfaces = mesh->getSize(2);

    JCircle circle;
    Point3D center, param;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            Point3D pa =  face->getNodeAt(0)->getXYZCoords();
            Point3D pb =  face->getNodeAt(1)->getXYZCoords();
            Point3D pc =  face->getNodeAt(2)->getXYZCoords();
            TriCircumCenter2D( &pa[0], &pb[0], &pc[0], &center[0], &param[0]);
            double radius = JMath::length(pa, center);
            double r2     = JMath::length(pb, center);
            double r3     = JMath::length(pc, center);
            assert( fabs(r2-radius) < 1.0E-06);
            assert( fabs(r3-radius) < 1.0E-06);
            circle.setCenter(center);
            circle.setRadius(radius);
            circumcircles.push_back(circle);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JTriDelaunayViewer :: draw()
{
    glColor4f( 1.0, 0.0, 0.0, 1.0);
    glDisable(GL_LIGHTING);

    glEnable( GL_DEPTH );
    glEnable( GL_BLEND);
    glEnable( GL_DEPTH_TEST);
    glEnable( GL_CULL_FACE);
    glEnable( GL_POLYGON_OFFSET_FILL);
    glPolygonOffset( 1.0, 1.0);

    if( displayCircumCircles ) {
        gluQuadricDrawStyle( diskObj, GLU_FILL);
        for( size_t i = 0; i < circumcircles.size(); i++) {
            Point3D xyz = circumcircles[i].getCenter();
            double  rad = circumcircles[i].getRadius();
            glPushMatrix();
            glTranslatef( xyz[0], xyz[1], xyz[2]);
            gluDisk( diskObj, 0.0, rad, 32, 1);
            glPopMatrix();
        }

        glLineWidth(1.0);
        glColor4f( 1.0, 1.0, 1.0, 1.0);
        gluQuadricDrawStyle( diskObj, GLU_SILHOUETTE);
        for( size_t i = 0; i < circumcircles.size(); i++) {
            Point3D xyz = circumcircles[i].getCenter();
            double  rad = circumcircles[i].getRadius();
            glPushMatrix();
            glTranslatef( xyz[0], xyz[1], xyz[2]);
            gluDisk( diskObj, 0.0, rad, 32, 1);
            glPopMatrix();
        }
    }

    if( displayMedialAxis ) {
        glLineWidth(2.0);
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
        size_t numfaces = kdt->getSize(2);
        glBegin( GL_TRIANGLES);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = kdt->getFaceAt(i);
            const Point3D &pa = f->getNodeAt(0)->getXYZCoords();
            const Point3D &pb = f->getNodeAt(1)->getXYZCoords();
            const Point3D &pc = f->getNodeAt(2)->getXYZCoords();
            glVertex2f( pa[0], pa[1] );
            glVertex2f( pb[0], pb[1] );
            glVertex2f( pc[0], pc[1] );
        }
        glEnd();
        glLineWidth(1.0);

        glPointSize(2.0);
        glColor4f( 1.0, 1.0, 1.0, 1.0);
        glBegin( GL_POINTS);
        for( size_t i = 0; i < medialPoints.size(); i++) {
            const Point3D &xyz =  medialPoints[i];
            glVertex3f( xyz[0], xyz[1], 0.001);
        }
        glEnd();
    }

    if( displayConvexHull) {
        size_t numedges = convexHull.size();

        glPointSize(10.0);
        glColor4f( 1.0, 0.0, 0.0, 1.0);
        glBegin( GL_POINTS);
        for( size_t i = 0; i < numedges; i++) {
            const Point3D &p0 = convexHull[i]->getNodeAt(0)->getXYZCoords();
            const Point3D &p1 = convexHull[i]->getNodeAt(1)->getXYZCoords();
            glVertex2f( p0[0], p0[1] );
            glVertex2f( p1[0], p1[1] );
        }
        glEnd();

        glColor4f( 1.0, 1.0, 0.0, 1.0);
        glLineWidth(2.0);
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);
        glBegin( GL_LINES);
        for( size_t i = 0; i < numedges; i++) {
            const Point3D &p0 = convexHull[i]->getNodeAt(0)->getXYZCoords();
            const Point3D &p1 = convexHull[i]->getNodeAt(1)->getXYZCoords();
            glVertex2f( p0[0], p0[1] );
            glVertex2f( p1[0], p1[1] );
        }
        glEnd();
    }

    return;
}

////////////////////////////////////////////////////////////////////////////////

JTriMesherDialog :: JTriMesherDialog( QWidget *parent) : QDialog(parent)
{
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    setupUi(this);
//  minAngleLineEdit->setText( QString::number(30 ));
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

JTriMesherDialog :: ~JTriMesherDialog()
{ }

///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c ) meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    creaseAngle = 0.0;
}
///////////////////////////////////////////////////////////////////////////////

void JTriMesherDialog :: showEvent( QShowEvent *event)
{
    if( meshViewer ) setMesh( meshViewer->getCurrentMesh() );
    QDialog::showEvent(event);
}
///////////////////////////////////////////////////////////////////////////////

void JTriMesherDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;
    objectNameLineEdit->setText( QString(mesh->getName().c_str() ));

}
///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: genNewMesh()
{
    if( mesh == nullptr ) return;
    string str = StdString(newStructComboBox->currentText());

    trimesher.setMesh(mesh);

    if( str == "ConstraintDelaunay") {
        if( cdt ) meshViewer->removeObject(cdt);
        cdt = trimesher.getSimpleMesh();
        if( cdt ) {
            JMeshRenderPtr mrender;
            mesh->getAttribute("Render", mrender);
            mrender->displayEntity[0] = 0;
            mrender->displayEntity[1] = 0;
            mrender->displayEntity[2] = 0;
            mrender->displayEntity[3] = 0;
            meshViewer->addObject(cdt);
        }
    }

    if( str == "ConvexHull") {
        if( chull ) meshViewer->removeObject(chull);
        chull = trimesher.getConvexHull();
        meshViewer->addObject(cdt);

    }

    if( str == "MedialAxis") {
        if( medialAxis ) meshViewer->removeObject(medialAxis);
        medialAxis = trimesher.getMedialAxis();
        if( medialAxis ) {
            JMeshRenderPtr mrender;
            mesh->getAttribute("Render", mrender);
            mrender->displayEntity[0] = 0;
            mrender->displayEntity[1] = 1;
            mrender->displayEntity[2] = 0;
            mrender->displayEntity[3] = 0;
            JEdgeRenderPtr eAttrib;
            size_t numedges = mesh->getSize(1);
            for( size_t i = 0; i < numedges; i++) {
                const JEdgePtr &edge = mesh->getEdgeAt(i);
                edge->getAttribute("Render", eAttrib);
                eAttrib->display = edge->isBoundary();
            }
            meshViewer->updateBuffers(mesh);
        }
        meshViewer->addObject(medialAxis);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JTriMesherDialog :: rejectNewMesh()
{
    if( cdt ) meshViewer->removeObject( cdt );
    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->displayEntity[0] = 0;
        mrender->displayEntity[1] = 0;
        mrender->displayEntity[2] = 0;
        mrender->displayEntity[3] = 0;
        meshViewer->refreshDisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////

void JTriMesherDialog :: openCleanupDialog()
{
    if( cleanupDialog == nullptr )
        cleanupDialog.reset(new JTrimeshCleanupDialog(this));

    cleanupDialog->setViewManager( viewManager );
    cleanupDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JTriMesherDialog :: openQualityDelaunayDialog()
{
    if( delaunayDialog == nullptr )
        delaunayDialog.reset(new JDelaunayMesherDialog(this));

    delaunayDialog->setViewManager( viewManager );
    delaunayDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

bool JTriMesherDialog :: event(QEvent *e)
{
    if( meshViewer == nullptr ) return 1;

    JMeshEntityPickerPtr picker = meshViewer->getEntityPicker();
    if( picker == nullptr ) return 1;

    /*
        JNodeSequence nodeSeq = picker->getPickedNodes();
        if( e->type() == QEvent::MouseMove && !nodeSeq.empty() ) {
            Point3D xyz = meshViewer->getViewManager()->getMouseCurrentWorldPosition();
            JNodePtr moveVertex = nodeSeq[0];
            moveVertex->setXYZCoords(xyz);
        }
    */
    QDialog::event(e);
    meshViewer->refreshDisplay();
    return 0;
}
///////////////////////////////////////////////////////////////////////////////

void JTriMesherDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JTriMesherDialog :: closeDialog()
{
    if( meshViewer == nullptr ) return;

    meshViewer->getViewManager()->detach(this);
    meshViewer->getViewManager()->setMouseTracking(0);
    this->close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: advancingFront()
{
    if( adfront) meshViewer->removeObject(adfront);
    JMeshGmshExporter  mexp;
    mexp.writeFile(mesh, "model.geo");

    system("gmsh model.geo -2 -algo front2d");

    // First write the geometry into "geo" format"
    adfront = JMeshIO::readFile( "model.msh");
    if( adfront == nullptr) return;
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 0;
    mrender->displayEntity[3] = 0;
    meshViewer->addObject(adfront);
}
///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: openIsotropicDialog()
{
    if( isotropicDialog == nullptr )
        isotropicDialog.reset(new JIsotropicTriMeshDialog(this));

    isotropicDialog->setViewManager( viewManager );
    isotropicDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: openInstantMeshDialog()
{
    if( instantMeshDialog == nullptr )
        instantMeshDialog.reset(new JInstantMeshDialog(this));

    meshViewer->setCurrentMesh(mesh);
    instantMeshDialog->setViewManager( viewManager );
    instantMeshDialog->setMeshType(3);
    instantMeshDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: getIntrinsicDelaunayMesh()
{
}

///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: openSurfReconstructionDialog()
{
    if( surfReconDialog == nullptr )
        surfReconDialog.reset(new JSurfaceReconstructionDialog(this));

    surfReconDialog->setViewManager( viewManager );
    surfReconDialog->show();
    surfReconDialog->setMesh(mesh);
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: rejectMesh()
{
}
///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: openContourEditingDialog()
{
    if( contourEditingDialog == nullptr )
        contourEditingDialog.reset(new JContourEditingDialog(this));

    contourEditingDialog->setViewManager( viewManager );
    contourEditingDialog->show();
    contourEditingDialog->setMesh(mesh);
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JTriMesherDialog :: openHolesFillDialog()
{
    if( holesFillDialog == nullptr )
        holesFillDialog.reset(new JMeshHolesFillDialog(this));

    holesFillDialog->setViewManager( viewManager );
    holesFillDialog->show();
    holesFillDialog->setMesh(mesh);
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////


void JTriMesherDialog :: makeConnections()
{
    PushButton( qualityDelaunayMesherPushButton,  [=] {openQualityDelaunayDialog();} );
    PushButton( newStructPushButton,  [=] {genNewMesh();} );
    PushButton( meshCleanupPushButton,[=] {openCleanupDialog();} );
    PushButton( advancingFrontPushButton, [=] {advancingFront();} );
    PushButton( isotropicMeshPushButton, [=] {openIsotropicDialog();} );
    PushButton( instantMeshPushButton, [=] {openInstantMeshDialog();} );
    PushButton( poissonPushButton, [=] {openSurfReconstructionDialog();} );
//  PushButton( intrinsicDelaunayMeshPushButton, [=] {getIntrinsicDelaunayMesh();} );
    PushButton( contourEditingPushButton, [=] {openContourEditingDialog();} );
    PushButton( holesFillingPushButton, [=] {openHolesFillDialog();} );
    PushButton( rejectPushButton, [=] {rejectMesh();} );

    PushButton( closePushButton, [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
void JTriMesherDialog :: genNewMesh()
{
    if( mesh == nullptr) return;

    int topDim = mesh->getTopology()->getDimension();

    if( topDim > 2) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Triangle mesh generation is for 2-manifold only");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    int elemType = mesh->getTopology()->getElementsType(2);

    ostringstream oss;
    QString qstr;
    oss << "triangle -";

    // Quality mesh generation ...
    if( opt_q_checkBox->isChecked() ) {
        oss << "q";

        //Set the min angle for Quality triangulation ....
        qstr = minAngleLineEdit->text();
        float minAngle = max(5.0, qstr.toDouble() );
        oss <<minAngle;
    }

    // Conforming Delaunay triangulation:  use this switch if you want to
    // ensure that all the triangles in the mesh are Delaunay, and not
    // merely constrained Delaunay; or if you want to ensure that all the
    // Voronoi vertices lie within the triangulation.  (Some finite volume
    // methods have this requirement.)  This switch invokes Ruppert's
    // original algorithm, which splits every subsegment whose diametral
    // circle is encroached.  It usually increases the number of vertices
    // and triangles.

    if( opt_D_checkBox->isChecked() )
        oss << "D";

    // Numbers all items starting from zero (rather than one)
    if( opt_z_checkBox->isChecked() )
        oss << "z";

    // No new vertices on the boundary.
    if( opt_Y_checkBox->isChecked() )
        oss << "Y";

    //Check the consistency of the final mesh.
    if( opt_C_checkBox->isChecked() )
        oss << "C";

    //Quiet: Suppresses all explanation of what Triangle is doing
    if( opt_Q_checkBox->isChecked() )
        oss << "Q";

    // Refines a previously generated mesh
    if( opt_r_checkBox->isChecked() )
        oss << "r";

    // Generates second-order subparametric elements with six nodes each
    if( opt_o2_checkBox->isChecked() )
        oss << "o2";

    // Use an incremental rather than a divide-and-conquer algorithm
    if( opt_i_checkBox->isChecked() )
        oss << "i";

    // Imposes a maximum triangle area.
    if( opt_a_checkBox->isChecked() ) {
        oss << "a";
        qstr = maxAreaLineEdit->text();
        double maxArea = max(0.0000001, qstr.toDouble() );
        oss <<maxArea;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    JEdgeSequence edges;
    JNodeSequence nodes;
    JMeshTRIExporter exporter;
    exporter.setDimension(2);

    if( topDim == 1 ) {
        oss.str("");
        oss << "triangle -p";
        mesh->getNodes( nodes );
        mesh->getEdges( edges );
        write_node_file(nodes);
        write_poly_file(edges);
        oss << " tmp";
    }

    if( topDim == 2 ) {
        if( interiorPointsCheckBox->isChecked() )  {
            exporter.writeNodes(mesh, "tmp.node");
            exporter.writeFaces(mesh, "tmp.ele");
            oss << " tmp";
        } else {
            oss.str("");
            oss << "triangle -p";
            mesh->getTopology()->getBoundary( nodes );
            mesh->getTopology()->getBoundary( edges );
            write_node_file(nodes);
            write_poly_file(edges);
            oss << " tmp";
        }
    }

    string cmd = oss.str();
    int stat = system( cmd.c_str() );
    waitCursor.stop();

    if( stat < 0) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Triangle Software probably not installed or not in the path");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    newMesh = Mesh::newObject();
    int err = newmesh->readFromFile( "tmp.1.ele");
    if( !err ) {
        newMesh->getTopology()->search_boundary();
        meshViewer->setNewMesh( newmesh );
        JMeshQuality mq(newmesh);
        vector<double> area;
        mq.getFacesQuality(JMeshQuality::AREA, area);
        double maxarea = *max_element(area.begin(), area.end() );
        maxAreaLineEdit->setText( QString::number(maxarea) );
    } else {
        newmesh->deleteAll();
        newmesh = nullptr;
    }
}
#endif

/*
void JTriMesherDialog :: baseTriangulated()
{
   ostringstream oss;
   oss << "triangle -pYQ  tmp";

  string cmd = oss.str();
    int stat = system( cmd.c_str() );

   int numfaces = mesh->getSize(2);
   for( size_t i  = 0; i < numfaces; i++) {
        Face *face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int nnodes = face->getSize(0);
            if( nnode > 3) {
                face->getNodes( nodes );
                face->getEdges( edges );
                write_node_file(nodes);
                write_poly_file(edges);
                oss << " tmp";
}
*/

/*
void JTriMesherDialog :: getMedialAxis()
{
    JEdgeSequence edges;
    JNodeSequence nodes;

    mesh->getTopology()->getBoundary( nodes );
    mesh->getTopology()->getBoundary( edges );
    write_node_file(nodes);
    write_poly_file(edges);

    ostringstream oss;
    QString qstr;
    oss << "triangle -p";
    oss << " tmp";

    string cmd = oss.str();
    int stat = system( cmd.c_str() );

    kdt = JMeshIO::readFile( "tmp.1.ele");
    size_t numfaces = kdt->getSize(2);
    vector<Point3D> medialPoints;
    if( numfaces ) medialPoints.resize(numfaces);
    double center[3];
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr f = kdt->getFaceAt(i);
        const Point3D &pa = f->getNodeAt(0)->getXYZCoords();
        const Point3D &pb = f->getNodeAt(1)->getXYZCoords();
        const Point3D &pc = f->getNodeAt(2)->getXYZCoords();
        TriCircumCenter2D( &pa[0], &pb[0], &pc[0], center);
        medialPoints[i][0] = center[0];
        medialPoints[i][1] = center[1];
        medialPoints[i][2] = 0.0;
    }
    meshViewer->displayAll(0);
    delaunayViewer->drawMedialAxis(1);
    delaunayViewer->setConstraintMesh( kdt );
    delaunayViewer->setMedial( medialPoints);
    meshViewer->refreshDisplay();
}
*/
///////////////////////////////////////////////////////////////////////////////

/*
void JTriMesherDialog :: getConvexHull()
{
    JNodeSequence nodes;

    if( mesh == nullptr) return;

    mesh->getTopology()->getBoundary( nodes );
    write_node_file(nodes);

    ostringstream oss;
    QString qstr;
    oss << "triangle ";
    oss << " tmp.node";

    string cmd = oss.str();
    int stat = system( cmd.c_str() );

    JMeshPtr newmesh = JMeshIO::readFile( "tmp.1.ele");
    JEdgeSequence boundedges;
    newmesh->getTopology()->getBoundary( boundedges );

    JMeshPtr ch = JMesh::newObject();

    string name = mesh->getName() + "_ConvexHull";
    ch->setName(name);

    meshViewer->addObject(ch);
}

void JTriMesherDialog :: paveMesh()
{
        if( meshViewer == nullptr) return;
        mesh = meshViewer->getMesh();
        if( mesh == nullptr ) return;

        cout << "Before Pave " << mesh->getSize(0) << " " << mesh->getSize(2) << endl;
        cout << mesh->getTopology()->getBoundarySize(0) << endl;

        int nlayer = numlayersSpinBox->value();

        JMeshPaver paver;
        paver.setMesh(mesh);
        paver.setNumBoundaryLayers(nlayer);
        paver.execute();
        meshViewer->setNewMesh(mesh);

        cout << "After Pave " << mesh->getSize(0) << " " << mesh->getSize(2) << endl;
        cout << mesh->getTopology()->getBoundarySize(0) << endl;
        meshViewer->refreshDisplay();
}
*/
/*
void JTriMesherDialog :: checkMaxArea()
{
    if( !opt_a_checkBox->isChecked() ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Maximum area checkbox is still disabled: Area constraint will not be applied");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
}
*/
///////////////////////////////////////////////////////////////////////////////
