#include "MeshInterpolationDialog.hpp"
#include "MeshQuality.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshInterpolationDialog :: JMeshInterpolationDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    meshViewer  = nullptr;
    imageViewer = nullptr;
    picker      = nullptr;
    thread.reset( new MyThread );
    limDeformer.reset( new JLocallyInjectiveMap );
    maxDistLineEdit->setText( QString::number(0.0) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshInterpolationDialog :: ~JMeshInterpolationDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: init()
{
    if( viewManager == nullptr) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");

    if( c == nullptr) c = JMeshViewer::registerComponent(viewManager);
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    c = viewManager->getComponent("ImageViewer");
    if(c) imageViewer = dynamic_pointer_cast<JImageViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: setSource( const JMeshPtr &m)
{
    sourceMesh = m;
    if( sourceMesh == nullptr ) return;

    string name = sourceMesh->getName();
    sourceMeshLineEdit->setText(QString(name.c_str()));

    /*
        vector<double> dist;
        if(targetMesh) {
           dist = JMeshGeometry::getEuclideanDistance(srcMesh, dstMesh);
        }

        int nCount = getInvertedElements( sourceMesh );
        if( nCount ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Warning: Source mesh contains inverted elements" );
            msg.setStandardButtons( QMessageBox::Cancel | QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok) {
                sourceMesh.reset();
                return;
            }
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: setTarget( const JMeshPtr &m)
{
    targetMesh = m;
    if( targetMesh == nullptr ) return;

    string name = targetMesh->getName();
    targetMeshLineEdit->setText(QString(name.c_str()));

    /*
        vector<double> dist;
        if( sourceMesh ) {
           dist = JMeshGeometry::getEuclideanDistance(srcMesh, dstMesh);
        }
        int nCount = getInvertedElements( targetMesh );
        if( nCount ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Warning: Target mesh contains inverted elements" );
            msg.setStandardButtons( QMessageBox::Cancel | QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok) {
                sourceMesh.reset();
                return;
            }
        }
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: mousePressEvent(QMouseEvent *e)
{
    if( !this->isVisible() ) return;
    if( picker == nullptr ) return;
    if( e->button() == Qt::LeftButton) left_button_pressed = 1;
    if( meshViewer) meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: genUVCoords()
{
}

///////////////////////////////////////////////////////////////////////////////

double JMeshInterpolationDialog :: getMaxDistance(const Point3D &xyz) const
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

void JMeshInterpolationDialog :: setMaxTargetDistance()
{
    maxDistLineEdit->setText( QString::number(0.0) );

    if( interpolatedMesh == nullptr) return;

    size_t numNodes = interpolatedMesh->getSize(0);
    double maxDist = 0.0;
    Point3D p0, p1;
    #pragma omp parallel for private(p0,p1) reduction(max:maxDist)
    for( size_t i = 0; i < numNodes; i++) {
        JNodePtr vtx = interpolatedMesh->getNodeAt(i);
        if( vtx->isActive() ) {
            if( vtx->hasAttribute("TargetPos") ) {
                vtx->getAttribute("TargetPos", p1);
                p0 = vtx->getXYZCoords();
                maxDist = std::max(maxDist, JMath::length2(p0, p1));

            }
        }
    }
    maxDist = sqrt(maxDist);
    maxDistLineEdit->setText( QString::number(maxDist) );
}

///////////////////////////////////////////////////////////////////////////////

int JMeshInterpolationDialog :: getInvertedElements( const JMeshPtr &mesh)
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

    JMeshQuality quality;
    double value;

    if( dim == 2) {
        size_t numfaces = mesh->getSize(2);
//        #pragma omp parallel for private(value, fattrib) reduction(+:numInverted)
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
//        #pragma omp parallel for private(value, cattrib) reduction(+:numInverted)
        for( size_t i = 0; i < numCells; i++) {
            const JCellPtr &cell = mesh->getCellAt(i);
            if( cell->isActive() ) {
                value  = quality.getJacobian(cell);
                if( value < 0.0 ) numInverted++;
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
    numInvertedLineEdit->setText( QString::number(numInverted) );
    return numInverted;
}

///////////////////////////////////////////////////////////////////////////////////a

void JMeshInterpolationDialog :: deformSolver()
{
    double t = paramScrollBar->value()/100.0;

    if( sourceMesh == nullptr || targetMesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    if( elasticDeformRadioButton->isChecked() ) {
        if( elasticDeform == nullptr) {
            elasticDeform.reset( new DDG::JMeshElasticDeformation);
            elasticDeform->setSource(sourceMesh);
            elasticDeform->setTarget(targetMesh);
            interpolatedMesh = elasticDeform->getInterpolatedMesh(0.0);
            meshViewer->addObject(interpolatedMesh);
        }
        interpolatedMesh = elasticDeform->getInterpolatedMesh(t);
        meshViewer->updateBuffers(interpolatedMesh);
    }

    if( harmonicMapRadioButton->isChecked() ) {
        if( harmonicMap == nullptr)
            harmonicMap.reset( new JHarmonicMap());

        cout << "Harmonic " << endl;
        harmonicMap->setSource(sourceMesh);
        harmonicMap->setTarget(targetMesh);

        interpolatedMesh = harmonicMap->getDeformedMesh();
        meshViewer->addObject(interpolatedMesh);
    }

    /*
       injectiveSolver();
        mesquiteSolver();
        meshViewer->updateBuffers(mesh);
        setMaxTargetDistance();
    */

    sourceMesh->setActiveBit(0);
    targetMesh->setActiveBit(0);
    meshViewer->refreshDisplay();
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: getWorstQuality()
{
    if( meshGeomQualityDialog == nullptr ) return;
    if( meshGeomQualityDialog->isVisible() ) {
        meshGeomQualityDialog->setHistogramBins(100);
        meshGeomQualityDialog->setValues();
    }
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: getHarmonicMap()
{

}
//////////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: openMeshGeomQualityDialog()
{
    if( meshGeomQualityDialog == nullptr )
        meshGeomQualityDialog.reset( new JMeshGeometricQualityDialog( this));

    meshGeomQualityDialog->setViewManager( viewManager );
    meshGeomQualityDialog->show();
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: injectiveSolver()
{
    int stat;
    stat = limDeformer->solve();

    if( stat != -2) meshViewer->updateBuffers(interpolatedMesh);
    return;

    /*
        stat = limDeformer->stepSolve();
            int nCount;
            nCount = invertedElements();

            if( nCount == 0) {
                if( str == "Jaal")
                    mesquiteSolver();
                else
            }
            swatch.stop();

            invertedElements();
            displayDeformedRadioButton->setChecked(1);
            meshViewer->refreshDisplay();
            getWorstQuality();
        */


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

    /*
        if( optMeshCheckBox->isChecked() )
            mesquiteSolver();
    */

    viewManager->refreshDisplay();
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: checkDisplay()
{
    if( limDeformer == nullptr ) return;
    bool val;

    val = displayMeshCheckBox->isChecked();
    if( meshViewer ) meshViewer->setActive(val);


    val = displayTextureCheckBox->isChecked();
    if( imageViewer ) imageViewer->setActive(val);

    /*
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
    */
}
//////////////////////////////////////////////////////////////////////////////////////
/*
void JMeshInterpolationDialog :: resetData()
{
    if( limDeformer == nullptr || meshViewer == nullptr ) return;
    if( sourceMesh ) {
        sourceMesh->getGeometry()->setCoordsArray(orgCoords,l2g);
        meshViewer->updateBuffers(sourceMesh);
    }
}
*/


//////////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: getImageMesh()
{
    this->init();

    if( imageViewer == nullptr )  return;
    JImagePtr img = imageViewer->getImage();
    if( img == nullptr ) return;

    if( meshViewer == nullptr )  return;

//  if( mesh == nullptr ) return;
//  img->setMesh(mesh);
}

//////////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: openLIMDialog()
{
    if( limParamsDialog == nullptr )
        limParamsDialog.reset( new JLocallyInjectiveMapParamsDialog);

    limParamsDialog->setDeformer(limDeformer.get() );
    limParamsDialog->show();
}

////////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: closeDialog()
{
    if( viewManager) {
        viewManager->detach( this );
        viewManager->setMouseTracking(0);
        viewManager->setRotationAxis(JaalViewer::FREE_ROTATION);
        viewManager->camera()->setType(Camera::PERSPECTIVE);
    }

    moveVertex = nullptr;
    this->close();
}

////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: showEvent( QShowEvent *e)
{
    if( e == nullptr) return;

    JMeshPtr currmsh = meshViewer->getCurrentMesh();
    if( selectMesh[0] ) setSource(currmsh);
    if( selectMesh[1] ) setTarget(currmsh);
}
////////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: loadSource()
{
    selectMesh[0] = 1;
    selectMesh[1] = 0;
    if( meshListDialog == nullptr )
        meshListDialog.reset(new JObjectsListDialog(this));

    meshListDialog->setViewManager(viewManager );
    meshListDialog->setType(0);
    meshListDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: loadTarget()
{
    selectMesh[0] = 0;
    selectMesh[1] = 1;

    if( meshListDialog == nullptr )
        meshListDialog.reset(new JObjectsListDialog(this));

    meshListDialog->setViewManager(viewManager);
    meshListDialog->setType(0);
    meshListDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: getHausdorff()
{
    /*
       if( sourceMesh == nullptr || targetMesh == nullptr) return;

       int elemType;
       JMeshPtr asurf = sourceMesh->getTopology()->getSurfaceMesh();

       elemType = asurf->getTopology()->getElementsType(2);
       JMeshPtr  trimesh1 = asurf;
       if( elemType == JFace::QUADRILATERAL) {
           AllTriMeshGenerator alltri;
           trimesh1 = alltri.getFromQuadMesh(asurf, 2);
       }

       JMeshPtr bsurf = targetMesh->getTopology()->getSurfaceMesh();
       elemType = bsurf->getTopology()->getElementsType(2);
       JMeshPtr  trimesh2 = bsurf;
       if( elemType == JFace::QUADRILATERAL) {
           AllTriMeshGenerator alltri;
           trimesh2 = alltri.fromQuadMesh(bsurf, 2);
       }
    */

}
///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog :: makeConnections()
{
    PushButton( hausdorffPushButton,   [=] {getHausdorff();});
    PushButton( limParamsPushButton,   [=] {openLIMDialog();});
    PushButton( meshQualityPushButton, [=] {openMeshGeomQualityDialog();});

    CheckBox( displayMeshCheckBox,     [=] {checkDisplay();});
    CheckBox( displayTextureCheckBox,  [=] {checkDisplay();});
    CheckBox( displayArrowsCheckBox,   [=] {checkDisplay();});

//    RadioButton( paramScrollBar, [=]{deformSolver();});
    RadioButton( displaySrcRadioButton, [=] {checkDisplay();});
    RadioButton( displayDeformedRadioButton,  [=] {checkDisplay();});

    PushButton( sourceMeshPushButton,  [=] {loadSource();});
    PushButton( targetMeshPushButton,  [=] {loadTarget();});
    PushButton( closePushButton,       [=] {closeDialog();});
    PushButton( applyPushButton,       [=] {deformSolver();});
}

///////////////////////////////////////////////////////////////////////////////

void JMeshInterpolationDialog ::  initMesh()
{
    /*
        int dim = mesh->getTopology()->getDimension();

        if( dim  < 2) {
            cout << "Warning: Mesh deformation is only for surface and volume mesh " << endl;
            return;
        }

        if( dim == 2 ) {
            if( mesh->getTopology()->getElementsType(2) != JFace::TRIANGLE )  {
                cout << "Warning: Mesh deformation only for the triangle mesh " << endl;
                return;
            }
            if( !mesh->getTopology()->isDisc() ) {
                if( tetmesherDialog == nullptr)
                    tetmesherDialog.reset( new JTetMesherDialog(this) );
                tetmesherDialog->setViewManager( viewManager );
                tetmesherDialog->setMesh(mesh);
                tetmesherDialog->show();
                this->hide();
                return;
            }
        }

        if( dim == 3 ) {
            if( mesh->getTopology()->getElementsType(3) != JCell::TETRAHEDRON )  {
                cout << "Warning: Mesh deformation only for the tetrahedral mesh " << endl;
                return;
            }
        }

        int val = 0;
        mesh->setAttribute("PickableEntity", val);

        mesh->getGeometry()->getCoordsArray(orgCoords,l2g);
        size_t numfaces = mesh->getSize(2);
        JFaceRenderPtr attrib;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Render", attrib);
                attrib->display = 0;
            }
        }
        picker = meshViewer->getEntityPicker();

            limDeformer->setMesh(mesh);
        initmesh = 1;
    */
}
///////////////////////////////////////////////////////////////////////////////

/*
void JMeshInterpolationDialog :: moveMesh()
{
if( mesh == nullptr ) return;
size_t numNodes = mesh->getSize(0);

    Point3D p0, p1;
    JNodeSequence boundnodes;
    vector<Point3D>  initPos, finalPos;


    numNodes = boundnodes.size();
    if( numNodes < 1) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("There are no constraints in the mesh to deform");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }
*/

/*
    int numSteps = numStepsSpinBox->value();
    double dt = 1.0/(double)numSteps;

    Point3D xyz;
    for( int istep = 0; istep < numSteps; istep++) {
        double t = (istep+1)*dt;
        for( int i = 0; i <  numNodes; i++) {
            xyz[0] = (1-t)*initPos[i][0] + t*finalPos[i][0];
            xyz[1] = (1-t)*initPos[i][1] + t*finalPos[i][1];
            xyz[2] = (1-t)*initPos[i][2] + t*finalPos[i][2];
            boundnodes[i]->setAttribute("TargetPos", xyz);
        }
        int err = runSolver();
        if( err ) return;
        meshViewer->refreshDisplay();
    }
}
*/
///////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: mouseMoveEvent(QMouseEvent *e)
{
    /*
        if( moveVertex == nullptr ) return;
        if( !this->isVisible() ) return;

        if( !moveHandleRadioButton->isChecked() ) return;

        Point3D  xyz;
        if( left_button_pressed ) {
            int err = viewManager->getMouseCurrentXYZPosition(xyz);
            if( !err) {
                xyz[2] = 0.0;   // we are dealing with 2D deformation ...
                moveVertex->setAttribute("TargetPos", xyz);
                runSolver();
            }
       }
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshInterpolationDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    /*
        if( picker == nullptr) return;

        JNodeSequence nodeSeq = picker->getPickedNodes();
        if( nodeSeq.empty() ) return;

        JNodeRenderPtr attrib;
        if( left_button_pressed && (e->modifiers() == Qt::ShiftModifier) ) {
            if( newHandleRadioButton->isChecked() ) {
                if( moveVertex ) {
                    moveVertex->deleteAttribute("TargetPos");
                    moveVertex->getAttribute("Render", attrib);
                    attrib->color[0]  = 0.0;
                    attrib->color[1]  = 1.0;
                    attrib->color[2]  = 0.0;
                    attrib->glyph     = JNodeDraw::NODE_AS_POINT;
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
                attrib->glyph   = JNodeDraw::NODE_AS_SPLAT;
                meshViewer->updateBuffers(mesh);
                updateConstraints();
            }
            if( moveHandleRadioButton->isChecked() ) {
                Point3D xyz;
                int err = viewManager->getMouseCurrentXYZPosition(xyz);
                if( !err) {
                    xyz[2] = 0.0;   // we are dealing with 2D deformation ...
                    moveVertex->setAttribute("TargetPos", xyz);
                    runSolver();
                }
            }
        }
        viewManager->refreshDisplay();
    */
}


void JMeshInterpolationDialog :: genImageMesh()
{
    /*
        if( imageViewer == nullptr )  return;
        JImage *img = imageViewer->getImage();
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
    */
}

