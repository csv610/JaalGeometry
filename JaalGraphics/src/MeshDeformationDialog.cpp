#include "MeshDeformationDialog.hpp"
#include "MeshQuality.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshDeformViewer :: JMeshDeformViewer()
{
    arrows = 0;
    viewManager   = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDeformViewer :: draw()
{
    if( mesh == nullptr ) return;

    glDisable( GL_LIGHTING );
    glDisable( GL_BLEND );

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    size_t numNodes = mesh->getSize(0);
    arrows  = 1;

    Point3D p0, p1;
    if( arrows ) {
        // Draw Arrows from the source vertices to the destination position.
        glColor3f( 1.0, 0.0, 1.0);
        glLineWidth(3.0);
        int nCount =0;
        for( size_t i = 0; i < numNodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            if( vtx->hasAttribute("TargetPos") ) {
                p0 = vtx->getXYZCoords();
                vtx->getAttribute("TargetPos", p1);
                glBegin( GL_LINES);
                glVertex3f( p0[0], p0[1], p0[2] + 0.001 );
                glVertex3f( p1[0], p1[1], p1[2] + 0.001 );
                nCount++;
                glEnd();
            }
        }
    }

    glColor3f( 1.0, 1.0, 1.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_QUADS);
    glVertex3f( -500, -500.0, -0.00001);
    glVertex3f(  500, -500.0, -0.00001);
    glVertex3f(  500,  500.0, -0.00001);
    glVertex3f( -500,  500.0, -0.00001);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    /*
        JNodeRenderPtr attrib;
        GLUquadricObj* diskObj = JNodeDraw::diskObj;
        for( size_t i = 0; i < numNodes; i++) {
            JNodePtr vtx = mesh->getNodeAt(i);
            vtx->getAttribute("Render", attrib);
            if( vtx->hasAttribute("TargetPos") ) {
                attrib->display = 1;
    //          attrib->glyph   = JNodeDraw::NODE_AS_CIRCLE;
                attrib->color   = srcColor;
            }
        }

        glDisable(GL_LIGHTING);
        numNodes = lassoPoints.size();
        if( numNodes > 1) {
            glLineWidth(2.0);
            glColor3f( 1.0, 0.0, 0.0);
            glBegin(GL_LINES);
            for( size_t i = 0; i < numNodes-1; i++) {
                const Point3D &p0  = lassoPoints[i];
                const Point3D &p1  = lassoPoints[i+1];
                glVertex3f( p0[0], p0[1], p0[2] + 0.001 );
                glVertex3f( p1[0], p1[1], p1[2] + 0.001 );
            }
            glEnd();
            glLineWidth(1.0);

            glPointSize(5.0);
            glColor3f( 0.0, 0.0, 1.0);
            glBegin(GL_POINTS);
            for( size_t i = 0; i < numNodes; i++) {
                const Point3D &p0  = lassoPoints[i];
                glVertex3f( p0[0], p0[1], p0[2] + 0.0001 );
            }
            glEnd();
            glPointSize(1.0);
        }


        glEnable( GL_LIGHTING );
    */
}

///////////////////////////////////////////////////////////////////////////////

JMeshDeformationDialog :: JMeshDeformationDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    meshViewer  = nullptr;
    imageViewer = nullptr;
    picker      = nullptr;
//    thread.reset( new MyThread );
//    limDeformer.reset( new JLocallyInjectiveMap );
    initmesh = 0;
    maxDistLineEdit->setText( QString::number(0.0) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshDeformationDialog :: ~JMeshDeformationDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: getNewConstraints()
{
    countNodes();
    assignColors();

//  if( deformViewer) deformViewer->setConstraints(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: showEvent( QShowEvent *)
{
    if( meshViewer == nullptr) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: setHandleAsConstraint()
{
    if( moveVertex == nullptr) return;

    int id = 0;
    JNodeRenderPtr attrib;
    Point3D pos = moveVertex->getXYZCoords();
    moveVertex->setAttribute("Constraint", id);
    moveVertex->getAttribute("Render", attrib);
    attrib->color[0]  = 1.0;
    attrib->color[1]  = 0.0;
    attrib->color[2]  = 0.0;
    attrib->display   = 1;
    meshViewer->updateBuffers(mesh);

    fixedNodes.insert(moveVertex);
    vector<int>  vnodes;
    for( const JNodePtr &vtx : fixedNodes)
        vnodes.push_back(vtx->getID() );

    limDeformer->setConstraints(vnodes) ;

    moveVertex = nullptr;
    picker->clearAll();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: keyPressEvent( QKeyEvent *e)
{
    // After selecting the vertex, if the user presses "P" Key then the
    // vertex becomes constraints vertex and changes color to red..
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }

    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: init()
{
    if( viewManager == nullptr) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");

    meshViewer = nullptr;
    if( c ) meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

    c = viewManager->getComponent("ImageViewer");
    if(c) {
        imageViewer = dynamic_pointer_cast<JImageViewer>(c);
    }

    // You need to capture mouse for defining the constraints, so register this
    // Object to the main class, so that events can be known to this class...
    viewManager->attach( this );

    // We need "DeformViewer" visual component, so that wr can modify teh
    // nodes, edges and face attributes...
    if( deformViewer == nullptr ) {
        deformViewer.reset( new JMeshDeformViewer() );
        deformViewer->setName("MeshDeformViewer");
        deformViewer->setViewManager(viewManager);

    }
    viewManager->attach(deformViewer);

    limDeformer.reset( new JLocallyInjectiveMap );

    setMesh( meshViewer->getCurrentMesh() );

    displaySrcRadioButton->setChecked(true);

    picker = meshViewer->getEntityPicker();
    picker->setMode(1);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    initmesh = 0;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    sourceMeshLineEdit->setText(QString(name.c_str()));

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 0;

    assignColors();
    countNodes();
    getInvertedElements();

    initMesh();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog ::  initMesh()
{
    entityDim = mesh->getTopology()->getDimension();

    if( entityDim == 2 ) {
        if( mesh->getTopology()->getElementsType(2) == JFace::QUADRILATERAL) {
            JMeshPtr tmpmesh = mesh->deepCopy();
            AllTriMeshGenerator alltri;
            simplicialMesh = alltri.getFromQuadMesh(tmpmesh, 4);
        } else
            simplicialMesh = mesh;
        assert( simplicialMesh->getTopology()->getElementsType(2) == JFace::TRIANGLE);
        limDeformer->setMesh(simplicialMesh);
//        limDeformer->initSolve();
    }

    if( entityDim == 3 ) {
        if( mesh->getTopology()->getElementsType(3) == JCell::HEXAHEDRON )  {
            JMeshPtr tmpmesh = mesh->deepCopy();
            AllTetMeshGenerator alltet;
            simplicialMesh = alltet.fromHexMesh(tmpmesh);
        } else
            simplicialMesh = mesh;
        assert( simplicialMesh->getTopology()->getElementsType(3) == JCell::TETRAHEDRON);
        limDeformer->setMesh(simplicialMesh);
//        limDeformer->initSolve();
    }

    mesh->getGeometry()->getCoordsArray(orgCoords,l2g);

    int val = 0;
    mesh->setAttribute("PickableEntity", val);

    if( entityDim == 2 ) {
        size_t numfaces = mesh->getSize(2);
        JFaceRenderPtr attrib;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Render", attrib);
                attrib->display = 0;
            }
        }
    }

    if( deformViewer ) deformViewer->setMesh( mesh );

    initmesh = 1;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: mousePressEvent(QMouseEvent *e)
{
    left_button_pressed = 0;
    if( e->button() == Qt::LeftButton) left_button_pressed = 1;
    if( meshViewer) meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: mouseMoveEvent(QMouseEvent *e)
{
    if( moveVertex == nullptr ) return;

    Point3D  xyz;
    if( e->modifiers() == Qt::ShiftModifier) {
        int err = viewManager->getMouseCurrentXYZPosition(xyz);
        if( !err) {
            xyz[2] = 0.0;   // we are dealing with 2D deformation ...
            moveVertex->setAttribute("TargetPos", xyz);
        }

        if( continuousDeformCheckBox->isChecked() )  {
            vector<int> vnodes;
            Eigen::Matrix<double,1,3>  pos;
            pos.coeffRef(0,0) = xyz[0];
            pos.coeffRef(0,1) = xyz[1];
            pos.coeffRef(0,2) = 0.0;
            vnodes.push_back( moveVertex->getID() );
            limDeformer->setConstraintsPosition( vnodes, pos);
            limDeformer->setMaxIterations(2);
            runSolver();
            glFinish();
        }
        viewManager->refreshDisplay();
    }
}
///////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    JNodeSequence nodeSeq = picker->getPickedNodes();
    if( nodeSeq.empty() ) return;

    JNodeRenderPtr attrib;
    if( left_button_pressed && (e->modifiers() == Qt::ShiftModifier) ) {
        if( moveVertex ) {
            moveVertex->deleteAttribute("TargetPos");
            moveVertex->getAttribute("Render", attrib);
            attrib->color[0]  = 0.0;
            attrib->color[1]  = 1.0;
            attrib->color[2]  = 0.0;
            attrib->glyph     = JNodeRender::NODE_AS_POINT;
            attrib->display   = 1;
            attrib->scale     = 1;
        }
        moveVertex = nodeSeq[0];
        const Point3D &pos = moveVertex->getXYZCoords();
        moveVertex->setAttribute("TargetPos", pos);
        moveVertex->getAttribute("Render", attrib);
        attrib->color[0] = 0.0;
        attrib->color[1] = 0.0;
        attrib->color[2] = 1.0;
        attrib->display =  1;
        attrib->glyph   = JNodeRender::NODE_AS_SPHERE;
        meshViewer->updateBuffers(mesh);

    }
    left_button_pressed = 0;
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: countNodes()
{
    if( mesh == nullptr);

    size_t numnodes = mesh->getSize(0);
    int ncount1 = 0;
    int ncount2 = 0;
    for(  size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->hasAttribute("Constraint")) ncount1++;
        if( vtx->hasAttribute("TargetPos")) ncount2++;
    }
    fixedNodesLineEdit->setText(QString::number(ncount1));
    targetNodesLineEdit->setText(QString::number(ncount2));
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: genUVCoords()
{
}

///////////////////////////////////////////////////////////////////////////////

double JMeshDeformationDialog :: getMaxDistance(const Point3D &xyz) const
{
    assert(nodeNeighs.size() );

    // Calculate the maximum distance the "Desired location" and the "Present location".

    double maxdist =  0.0;
    for( size_t i= 0; i < nodeNeighs.size(); i++) {
        double dist =  JMath::length2( nodeNeighs[i]->getXYZCoords(), xyz);
        maxdist = std::max(maxdist, dist);
    }

    return sqrt(maxdist);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: setMaxTargetDistance()
{
    maxDistLineEdit->setText( QString::number(0.0) );
    if( mesh == nullptr) return;

    size_t numNodes = mesh->getSize(0);
    double maxDist = 0.0;
    Point3D p1;
    for( size_t i = 0; i < numNodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            if( vtx->hasAttribute("TargetPos") ) {
                vtx->getAttribute("TargetPos", p1);
                const Point3D &p0 = vtx->getXYZCoords();
                maxDist = std::max(maxDist, JMath::length2(p0, p1));

            }
        }
    }
    maxDist = sqrt(maxDist);
    maxDistLineEdit->setText( QString::number(maxDist) );
}

///////////////////////////////////////////////////////////////////////////////

int JMeshDeformationDialog :: getInvertedElements()
{
    if( mesh == nullptr ) return 1;

    JColor redColor, greenColor;
    redColor[0] = 1.0;
    redColor[1] = 0.0;
    redColor[2] = 0.0;
    redColor[3] = 1.0;

    greenColor[0] = 0.0;
    greenColor[1] = 1.0;
    greenColor[2] = 0.0;
    greenColor[3] = 1.0;

    JFaceRenderPtr fattrib;
    JCellRenderPtr cattrib;

    int numInverted = 0;

    int dim = mesh->getTopology()->getDimension();

//  #pragma omp parallel for private(attrib) reduction(+:nCount)
    JMeshQuality quality;
    double value;

    if( dim == 2) {
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                value = quality.getJacobian(face);
                if( value < 0 ) numInverted++;

                int err = face->getAttribute("Render", fattrib);
                if( !err) {
                    if( value > 0)
                        fattrib->color = greenColor;
                    else
                        fattrib->color = redColor;
                }
            }
        }
    }

    if( dim == 3) {
        size_t numCells = mesh->getSize(3);
        for( size_t i = 0; i < numCells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                value =  quality.getJacobian(cell);
                if( value < 0 ) numInverted++;
                int err = cell->getAttribute("Render", cattrib);
                if( !err) {
                    if( value > 0)
                        cattrib->color = greenColor;
                    else
                        cattrib->color = redColor;
                }
            }
        }
    }

    return numInverted;
}

///////////////////////////////////////////////////////////////////////////////////a

void JMeshDeformationDialog :: runSolver()
{
    if( !initmesh ) initMesh();

    if( initmesh == 0) {
        cout << "Warning: Preprocessing not done: Deformation not starting " << endl;
        return;
    }

    size_t numnodes = mesh->getSize(0);
    Point3D xyz;
    #pragma omp parallel for private(xyz)
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = mesh->getNodeAt(i);
        const JNodePtr &v1 = simplicialMesh->getNodeAt(i);
        int err = v0->getAttribute("TargetPos", xyz);
        if( err)
            v1->deleteAttribute("TargetPos");
        else
            v1->setAttribute("TargetPos", xyz);
    }

    JWaitCursor waitCursor;
    waitCursor.start();
    injectiveSolver();
//  mesquiteSolver();
    meshViewer->updateBuffers(mesh);
    setMaxTargetDistance();
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: getWorstQuality()
{
    if( meshGeomQualityDialog == nullptr ) return;
    if( meshGeomQualityDialog->isVisible() ) {
        meshGeomQualityDialog->setHistogramBins(100);
        meshGeomQualityDialog->setValues();
    }
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: openMeshGeomQualityDialog()
{
    if( meshGeomQualityDialog == nullptr )
        meshGeomQualityDialog.reset( new JMeshGeometricQualityDialog( this));

    meshGeomQualityDialog->setViewManager( viewManager );
    meshGeomQualityDialog->show();
}
//////////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: openCurveShorteningDialog()
{
    if( curveShorteningDialog == nullptr )
        curveShorteningDialog.reset( new JCurveShorteningFlowDialog( this));

    curveShorteningDialog->setViewManager( viewManager );
    curveShorteningDialog->setMesh(mesh);

    curveShorteningDialog->show();
    this->hide();
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: mesquiteSolver()
{
    if( mesh == nullptr ) return;

    Point3D tpos, spos, xyz;
    size_t numnodes = mesh->getSize(0);
    vector<Point3D>  srcPoints, dstPoints;

    JNodeSequence  movenodes;

    int grp = 1;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        int err = vtx->getAttribute("TargetPos", tpos);
        if( !err)  {
            spos = vtx->getXYZCoords();
            movenodes.push_back(vtx);
            srcPoints.push_back(spos);
            dstPoints.push_back(tpos);
            vtx->setAttribute("Constraint", grp);
        }
    }
    meshOpt.setMesh(mesh);

    int n = 100;
    double dt = 1/(double)n;

    JLaplaceMeshSmoother mlap;
    mlap.setMesh(mesh);
    mlap.setNumIterations(5);

    int nCount;
    for( int i = 1; i < n; i++) {
        double t = i*dt;

        // Assign new position to the deformed nodes ...
        size_t vid = 0;
        for( const JNodePtr &vtx : movenodes) {
            spos   = vtx->getXYZCoords();
            xyz[0] = (1-t)*spos[0] + t*dstPoints[vid][0];
            xyz[1] = (1-t)*spos[1] + t*dstPoints[vid][1];
            xyz[2] = (1-t)*spos[2] + t*dstPoints[vid][2];
            vtx->setXYZCoords(xyz);
            vid++;
        }
        // Smooth all nodes except on the source ..
        mlap.smoothAll();

        // Check for inverted elements, if there are no inverted elements,
        // we can improve the shapes of the elements...
        nCount = getInvertedElements();
        if( nCount == 0) meshOpt.improveQuality();

        // if there are tangled elements, try to untangle them..
        if( nCount ) meshOpt.untangle();
        // Prepare for rendering ..
        mesh->getGeometry()->setFacesNormal();
        meshViewer->updateBuffers(mesh);
        // Check again for the inverted elments. if the mesh is tangle
        // free, we can go for the next iteation.

        nCount = getInvertedElements();
        if( nCount ) {
            mesh->deleteNodeAttribute("Constraint");
            meshOpt.untangle();
            nCount = getInvertedElements();
            if( nCount) {
                cout << "Info: Mesh not recovered from inverted elements " << endl;
                return;
            }
            for( const JNodePtr &vtx : movenodes)
                vtx->setAttribute("Constraint", grp);

        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: injectiveSolver()
{
    int stat;
    JWaitCursor waitCursor;
    waitCursor.start();

    stat = limDeformer->solve();
//   stat = limDeformer->initSolve();
//   stat = limDeformer->stepSolve();

    size_t numnodes = mesh->getSize(0);
    #pragma omp parallel for
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &v0 = simplicialMesh->getNodeAt(i);
        const Point3D  &p0 = v0->getXYZCoords();
        const JNodePtr &v1 = mesh->getNodeAt(i);
        v1->setXYZCoords(p0);
    }

    if( stat != -2) meshViewer->updateBuffers(mesh);
    return;

#ifdef USE_THREAD
    if( limDeformer == nullptr )
        limDeformer.reset( new JLocallyInjectiveMap );
    assert( limDeformer );
    limDeformer->setMesh(mesh);

    if( limParamsDialog == nullptr) {
        limParamsDialog.reset( new JLocallyInjectiveMapParamsDialog);
        limParamsDialog->setDeformer(limDeformer.get() );
    }
    limParamsDialog->setParams();

    if( thread ) {
        if( thread->isRunning() ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Previous thread is still running: Do you wish canel previous and start new thread" );
            msg.setStandardButtons( QMessageBox::Cancel | QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Cancel) return;
            thread->terminate();
        }
    }

    thread->deformer    = limDeformer;
    thread->viewManager = viewManager;
    thread->start();
#else
    limDeformer->stepSolve();
#endif

    if( optShapeCheckBox->isChecked() )
        mesquiteSolver();

    viewManager->refreshDisplay();
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: checkDisplay()
{
    if( limDeformer == nullptr ) return;
    bool val;

    /*
        val = displayMeshCheckBox->isChecked();
        if( meshViewer ) meshViewer->setActive(val);


        val = displayTextureCheckBox->isChecked();
        if( imageViewer ) imageViewer->setActive(val);

    */
    if( displayDeformedRadioButton->isChecked() ) {
        limDeformer->setDeformedMesh();
    }

    if( displaySrcRadioButton->isChecked() ) {
        if( mesh )
            mesh->getGeometry()->setCoordsArray(orgCoords, l2g);
        meshViewer->updateBuffers(mesh);
    }

    val = displayArrowsCheckBox->isChecked();
    if( deformViewer ) deformViewer->setArrows(val);

    getInvertedElements();

    meshViewer->refreshDisplay();
}
//////////////////////////////////////////////////////////////////////////////////////
/*
void JMeshDeformationDialog :: openConstraintsDialog()
{
    if( meshViewer == nullptr ) return;

    if( meshConstraintsDialog == nullptr )
        meshConstraintsDialog.reset(new JMeshConstraintsDialog(this));
    viewManager->attach(meshConstraintsDialog.get());

    QObject::connect( meshConstraintsDialog.get(), SIGNAL(setNewConstraints()), this, SLOT(getNewConstraints()));

    meshConstraintsDialog->setViewManager( viewManager );
    meshConstraintsDialog->setMesh(mesh);
    meshConstraintsDialog->show();
    this->hide();
}
*/

//////////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: resetData()
{
    if( limDeformer == nullptr || meshViewer == nullptr || mesh == nullptr) return;

    mesh->getGeometry()->setCoordsArray(orgCoords,l2g);
    meshViewer->updateBuffers(mesh);
    getInvertedElements();
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: getImageMesh()
{
    this->init();

#ifdef CSV

    if( imageViewer == nullptr )  return;
    JImagePtr img = imageViewer->getImage();
    if( img == nullptr ) return;

    if( meshViewer == nullptr )  return;
    if( mesh == nullptr ) return;
//  img->setMesh(mesh);
#endif
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: genImageMesh()
{
#ifdef CSV
    if( imageViewer == nullptr )  return;
    JImagePtr img = imageViewer->getImage();
    if( img == nullptr ) return;

    double len[2];
    double origin[2];

    if( structmeshDialog == nullptr ) {
        structmeshDialog.reset( new JStructuredMeshDialog( this ));
        QObject::connect( structmeshDialog.get(), SIGNAL(meshCreated()), this, SLOT(getImageMesh()));
    }

    structmeshDialog->setViewManager( viewManager );
    int width  = img->getWidth();
    int height = img->getHeight();
    if( width > height ) {
        len[0] = width/(double)height;
        len[1]  = 1.0;
    } else {
        len[0] = 1.0;
        len[1] = height/(double)width;
    }
    origin[0] = -0.5*len[0];
    origin[1] = -0.5*len[1];
    structmeshDialog->setOrigin(origin[0], origin[1], 0.0);
    structmeshDialog->setLength(len[0], len[1], 0.0);
    structmeshDialog->setUVCoords(1);
    structmeshDialog->show();
#endif
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: openLIMDialog()
{
    if( limParamsDialog == nullptr )
        limParamsDialog.reset( new JLocallyInjectiveMapParamsDialog);

    limParamsDialog->setDeformer(limDeformer.get() );
    limParamsDialog->show();
}

////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: dst2src()
{
    if( mesh == nullptr ) return;
    size_t numNodes = mesh->getSize(0);
    Point3D xyz;

    size_t index = 0;
    for( size_t i = 0; i < numNodes; i++) {
        JNodePtr vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            if( vtx->hasAttribute("TargetPos") ) {
                xyz[0] = orgCoords[3*index+0];
                xyz[1] = orgCoords[3*index+1];
                xyz[2] = orgCoords[3*index+2];
                vtx->setAttribute("TargePos", xyz);
            }
            index += 3;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: assignColors()
{

    if( mesh == nullptr) return;

    JNodeRenderPtr attrib;
    size_t numnodes = mesh->getSize(0);
    int val = 0;
    for( size_t i = 0; i < numnodes; i++)  {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Render", attrib);

        attrib->glyph   = JNodeRender::NODE_AS_POINT;
        attrib->color[0] = 0.0;
        attrib->color[1] = 1.0;
        attrib->color[2] = 0.0;
        attrib->display  = 1;

        if( vtx->hasAttribute("Constraint" ) ) {
            attrib->color[0] = 1.0;
            attrib->color[1] = 0.0;
            attrib->color[2] = 0.0;
            attrib->display =  1.0;
//            attrib->glyph   = JNodeDraw::NODE_AS_SPLAT;
        }

        if( vtx->hasAttribute("TargetPos" ) ) {
            attrib->color[0] = 0.0;
            attrib->color[1] = 0.0;
            attrib->color[2] = 1.0;
            attrib->display =  1;
//          attrib->glyph   = JNodeDraw::NODE_AS_SPLAT;
        }
    }
    meshViewer->updateBuffers(mesh);
}
////////////////////////////////////////////////////////////////////////////////
void JMeshDeformationDialog :: setFixedBoundary()
{
    bool val = fixedBoundaryCheckBox->isChecked();

    int id = 1;
    JNodeRenderPtr attrib;
    size_t numnodes = mesh->getSize(0);
    for( int i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Render", attrib);
        if( vtx->isBoundary() ) {
            if( val == 0) {
                attrib->glyph     = JNodeRender::NODE_AS_POINT;
                vtx->deleteAttribute("Constraint");
            }
            else {
                vtx->setAttribute("Constraint", id);
                attrib->color[0]  = 1.0;
                attrib->color[1]  = 0.0;
                attrib->color[2]  = 0.0;
                attrib->display   = 1;
                attrib->glyph     = JNodeRender::NODE_AS_SPHERE;
                fixedNodes.insert(vtx);
            }
        }
    }

    vector<int>  vnodes;
    for( const JNodePtr &vtx : fixedNodes)
        vnodes.push_back(vtx->getID() );
    limDeformer->setConstraints(vnodes) ;
    meshViewer->updateBuffers(mesh);
}
////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: closeDialog()
{
    if( viewManager) {
        viewManager->detach( this );
        viewManager->detach( deformViewer );
        viewManager->setMouseTracking(0);
        viewManager->setRotationAxis(JaalViewer::FREE_ROTATION);
        viewManager->camera()->setType(Camera::PERSPECTIVE);
    }

    moveVertex = nullptr;
    shrink_to_zero( orgCoords);
    shrink_to_zero( l2g );
    this->close();
}

////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: loadSource()
{
    static QString lastSelectedDirectory;
    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Mesh File "),
                   lastSelectedDirectory,
                   *new QString( "Mesh Format (*.xml *.off *obj)"));

    string fileName = qstr.toUtf8().constData();
    if (fileName.empty()) return;

    JWaitCursor waitCursor;
    waitCursor.start();
    mesh = JMeshIO::readFile(fileName);

    if( mesh == nullptr) return;

    meshViewer->addObject(mesh);

    sourceMeshLineEdit->setText(QString(mesh->getName().c_str()));
    displayTargetRadioButton->setChecked(false);
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: loadTarget()
{
    if( mesh == nullptr) {
        cout << "First load source mesh, then target mesh " << endl;
        return;
    }

    static QString lastSelectedDirectory;
    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Mesh File "),
                   lastSelectedDirectory,
                   *new QString( "Mesh Format (*.xml *.off *obj)"));

    string fileName = qstr.toUtf8().constData();
    if (fileName.empty()) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JMeshPtr targetMesh = JMeshIO::readFile(fileName);

    if( targetMesh == nullptr ) return;

    if( !targetMesh->getTopology()->isSameAs(mesh) ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning: Source and destination mesh have different toplogy");
        msg.setStandardButtons(QMessageBox::Ok);
        msg.exec();
        return;
    }

    JMeshPtr bdmesh = targetMesh;
    if( mesh->getTopology()->getDimension() == 3)
        bdmesh = mesh->getTopology()->getSurfaceMesh();

    if( bdmesh->getSize(2) != mesh->getSize(2) ) {
        cout << "Warning: #Boundary faces on two mesh are different " << endl;
        return;
    }

    if( bdmesh->getSize(0) != mesh->getSize(0) ) {
        cout << "Warning: #Boundary faces on two mesh are different " << endl;
        return;
    }
    meshViewer->addObject( targetMesh);

    /*
        int numnodes = bdmesh->getSize(0);
        for( int i = 0; i < numnodes; i++) {
            const JNodePtr &srcnode = mesh->getNodeAt(i);
            const JNodePtr &dstnode = bdmesh->getNodeAt(i);
            srcnode->setAttribute("TargetPos", dstnode->getXYZCoords() );
        }
        countNodes();
        setMaxTargetDistance();
        targetMeshLineEdit->setText(QString(targetMesh->getName().c_str()));
        displayTargetRadioButton->setChecked(true);
    */
}

////////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: makeConnections()
{
    PushButton( resetPushButton,       [=] {resetData();});
    PushButton( limParamsPushButton,   [=] {openLIMDialog();});
    PushButton( meshQualityPushButton, [=] {openMeshGeomQualityDialog();});
    PushButton( sourceMeshPushButton,  [=] {loadSource();});
    PushButton( runSolverPushButton,   [=] {runSolver();});
    PushButton( curveShorteningPushButton,  [=] {openCurveShorteningDialog();});
    PushButton( setHandleConstraintPushButton, [=] {setHandleAsConstraint();});

    CheckBox( displayMeshCheckBox, [=] {checkDisplay();});
    CheckBox( displayArrowsCheckBox, [=] {checkDisplay();});
    CheckBox( displayTextureCheckBox, [=] {checkDisplay();});

    RadioButton( displaySrcRadioButton, [=] {checkDisplay();});
    RadioButton( displayDeformedRadioButton, [=] {checkDisplay();});

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDeformationDialog :: src2dst()
{
    if( mesh == nullptr ) return;

    const JBoundingBox &box = mesh->getGeometry()->getBoundingBox();

    JEdgeSequence boundedges;
    mesh->getTopology()->getBoundary(boundedges);
    JEdgeTopology::getChain( boundedges);

    // Since the orientation of the first edge is the orientation of the triangle, so
    // the chain will be anti-clockwise, if the triangle mesh is consistent and all the
    // triangles have anti-clockwise orientation.
    // So, we need not check the orientation of the boundary nodes...

    JNodeSequence boundnodes;
    JEdgeTopology::getChainNodes( boundedges, boundnodes);

    /*
        QString qs = projectComboBox->currentText();
        string str = qs.toUtf8().constData();


        size_t nSize = boundnodes.size();
        Point3D xyz;

        vector<Point3D> points;
        if( str == "Circle") {
            // We select the largest circle enclosing the given shape...
            double  len1   =  JMath::length(box.getCorner(0), box.getCorner(2));
            double  len2   =  JMath::length(box.getCorner(1), box.getCorner(3));
            double  radius =  0.5*std::max(len1, len2);
            Point3D center =  box.getCenter();
            double dtheta = 2.0*M_PI/(double)nSize;
            points.resize(nSize);
            for( size_t i = 0; i < nSize; i++) {
                points[i][0] = center[0] + radius*cos( i*dtheta);
                points[i][1] = center[1] + radius*sin( i*dtheta);
            }
            deformViewer->setLasso(points);

            // User Can select the first vertex which will be mapped to the first
            // vertex of the parametric boundary. This point should be judiciously
            // picked, and the quality of the mapping will depend on this value..
            for( size_t i = 0; i < nSize; i++) {
                if( boundnodes[i] == moveVertex ) {
                    for( int j = 0; j <  nSize; j++)
                        boundnodes[(i+j)%nSize]->setAttribute("TargetPos", points[j]);
                    break;
                }
            }
        }
    */
}



