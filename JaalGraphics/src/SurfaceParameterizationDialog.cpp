#include "SurfaceParameterizationDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JSurfaceParameterizationDialog :: JSurfaceParameterizationDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterizationDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterizationDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if(  nullptr == mesh ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////
int JSurfaceParameterizationDialog :: initMesh()
{
    if( nullptr == mesh ) return 1;

    int euler = mesh->getTopology()->getEulerCharacteristic();

    if( euler == 1) {
        numChartsLineEdit->setText( QString::number(1) );
        return 0;
    }

    meshPart.reset( new JMeshPartitioner);
    meshPart->setMesh(mesh);
    int numParts = meshPart->getNumPartitions();
    numChartsLineEdit->setText( QString::number(numParts) );
    if( numParts ) {
        for( int i = 0; i < numParts; ++i) {
            JMeshPtr apart = meshPart->getSubMesh(i);
            if( !apart->getTopology()->isDisk() ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("Submesh is not topological disk");
                msg.setStandardButtons( QMessageBox::Ok);
                int ret = msg.exec();
                if( ret == QMessageBox::Ok ) return 1;
            }
            return 0;
        }
    }

    vector<JEdgeSequence> edges;
    mesh->getTopology()->getBoundary(edges);
    if( edges.size() == 1) {
        numChartsLineEdit->setText( QString::number(1) );
        return 0;
    }

    openMeshPartitionDialog();
    return 1;
}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterizationDialog :: smoothUVBoundary()
{
}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterizationDialog :: untangleUVMesh()
{
    if( surfParam == nullptr) return;
    vector<JMeshPtr>  invertedUVMesh = surfParam->getInvertedUVMeshes();


}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterizationDialog :: setParameters()
{
    QString qstr;
    string  str;

    qstr = weightTypeComboBox->currentText();
    str  = StdString(qstr);

    int w = 0;
    if( str == "ShapePreserving") w  = 0;
    if( str == "Tutte") w  = 1;
    if( str == "HarmonicMap") w  = 2;
    if( str == "IntrinsicMap") w  = 3;
    if( str == "MeanValue") w  = 4;
    if( str == "LeastSquareConformal") w  = 5;
    if( str == "AsRigidAsPossible") w  = 6;
    surfParam->setWeight(w);

    qstr = boundaryMapComboBox->currentText();
    str  = StdString(qstr);

    int b;
    if( str == "Square") b = 0;
    if( str == "Circle") b = 1;
    if( str == "NaturalBoundary") b = 2;
    surfParam->setBoundary(b);

    bool val = areaRescaleCheckBox->isChecked();
    surfParam->setRescale(val);
}

///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterizationDialog :: displayCharts()
{
    if( uvMeshes.empty() ) return;

    JBoundingBox box;
    double maxlen = 0.0;

    for( const JMeshPtr &uvmesh : uvMeshes ) {
        box = uvmesh->getGeometry()->getBoundingBox();
        maxlen = max( maxlen, box.getLength(0));
        maxlen = max( maxlen, box.getLength(1));
    }
    int nx = xCellsSpinBox->value();
    int ny = yCellsSpinBox->value();
    maxlen *= 1.01;

    JMeshAffineTransform affine;
    int index = 0;
    for( int j = 0; j < ny; j++) {
        for( int i = 0; i < nx; i++) {
            if( index >= uvMeshes.size() ) break;
            JMeshPtr uvmesh = uvMeshes[index++];
            double   xc = (i+0.5)*maxlen;
            double   yc = (j+0.5)*maxlen;
            affine.setMesh(uvmesh);
            affine.setCenter(xc, yc, 0.0);
        }
    }

    meshViewer->lookAt( uvMeshes[0] );

    JFaceRenderPtr fAttrib;
    JColor color;

    JColor white = JEntityColor::getColor("White");

    color = white;
    for( const JMeshPtr &uvm : uvMeshes) {
        if( coloredChartsCheckBox->isChecked() ) {
            uvm->getFaceAt(0)->getAttribute("Render", fAttrib);
            color = fAttrib->color;
        }
        for( size_t j = 0; j < uvm->getSize(2); j++) {
            const JFacePtr &f  = uvm->getFaceAt(j);
            f->getAttribute("Render", fAttrib);
            fAttrib->color = color;
        }
        meshViewer->updateBuffers(uvm);
    }

    checkDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterizationDialog :: checkDisplay()
{
    bool val;
    JMeshRenderPtr mrender;

    val = displayChartsRadioButton->isChecked();
    for( const JMeshPtr &uvmesh : uvMeshes ) {
        uvmesh->getAttribute("Render", mrender);
        mrender->displayEntity[0] = 0;
        mrender->displayEntity[1] = val;
        mrender->displayEntity[2] = val;
        mrender->displayEntity[3] = 0;
    }

    val = displayMeshRadioButton->isChecked();
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = val;
    mrender->displayEntity[2] = val;
    mrender->displayEntity[3] = 0;

    if( displayChartsRadioButton->isChecked() )
        meshViewer->lookAt( uvMeshes[0] );
    else
        meshViewer->lookAt( mesh );

    viewManager->showEntireScene();
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceParameterizationDialog :: openMeshPartitionDialog()
{
    if( meshPartDialog == nullptr )
        meshPartDialog.reset( new JMeshPartitionDialog(this));

    meshPartDialog->setViewManager( viewManager );
    meshPartDialog->setMesh(mesh);
    meshPartDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterizationDialog :: getUVCharts()
{
    int err = initMesh();
    if( err ) return;

    JWaitCursor wCursor;
    wCursor.start();

    // Remove any prevously generated mesh ...
    for( const JMeshPtr &uvmesh : uvMeshes )
        meshViewer->removeObject(uvmesh);

    uvMeshes.clear();

    if( surfParam == nullptr)
        surfParam.reset( new JSurfaceParameterization);

    surfParam->setMesh(mesh);
    setParameters();
    surfParam->genCharts();

    JFaceColorPtr faceColor(new JFacePartitionColor);
    faceColor->setMesh(mesh);

    vector<JMeshPtr> submeshes = surfParam->getSubMeshes();

    JMeshPtr uvm;
    JFaceRenderPtr fAttrib;
    JColor color;

    for( const JMeshPtr &m : submeshes) {
        int err = m->getAttribute("UVMesh", uvm);
        if( !err ) {
            uvMeshes.push_back(uvm);
            meshViewer->addObject(uvm);
            m->getFaceAt(0)->getAttribute("Render", fAttrib);
            color = fAttrib->color;
            for( size_t j = 0; j < uvm->getSize(2); j++) {
                const JFacePtr &f  = uvm->getFaceAt(j);
                f->getAttribute("Render", fAttrib);
                fAttrib->color = color;
            }
            meshViewer->updateBuffers(uvm);
        }
    }

    int nc = submeshes.size();
    numChartsLineEdit->setText( QString::number( nc ));
    int nx = sqrt( nc );
    int ny = nc/nx + 1;
    xCellsSpinBox->setValue( nx );
    yCellsSpinBox->setValue( ny );

    displayCharts();
}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterizationDialog :: closeDialog()
{
    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &f = mesh->getFaceAt(i);
        f->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
    }

    // Remove any prevously generated mesh ...
    for( const JMeshPtr &uvmesh : uvMeshes )
        meshViewer->removeObject(uvmesh);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 1;
    mrender->displayEntity[1] = 1;
    mrender->displayEntity[2] = 1;
    mrender->displayEntity[3] = 0;

    meshViewer->lookAt( mesh );
    meshViewer->updateBuffers(mesh);

    this->close();
}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceParameterizationDialog :: makeConnections()
{
    RadioButton(displayChartsRadioButton,  [=] {checkDisplay(); });
    RadioButton(displayMeshRadioButton,  [=] {checkDisplay(); });
    PushButton( applyPushButton,  [=] {getUVCharts(); });
    SpinBoxi( xCellsSpinBox, [=] {displayCharts();});
    SpinBoxi( yCellsSpinBox, [=] {displayCharts();});

    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////

JSurfaceMap2DViewer  :: JSurfaceMap2DViewer()
{
    projCenter[0] = 300;
    projCenter[1] = 300;
    projWidth     = 600;
    projHeight    = 600;
    texID         = 0;
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceMap2DViewer  :: setPatch( const JMeshPtr &m)
{
    patch = m;
    if( patch == nullptr) return;

    JBoundingBox box = patch->getGeometry()->getBoundingBox();
    Point3D center = box.getCenter();
    glDeleteTextures(1, &texID);
    texID = 0;
    showChart = 0;
    int id = 0;
    patch->getAttribute("PatchID", id);
    mapPatches[id] = patch;
}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceMap2DViewer  :: getImage()
{
//  viewManager->refreshDisplay();
    int winHeight = viewManager->height();

    QImage qimg = viewManager->grabFrameBuffer();
    int xt = projCenter[0] - 0.5*projWidth;
    int yt = winHeight - (projCenter[1] + 0.5*projHeight);
    QImage qsmall = qimg.copy(xt, yt, projWidth, projHeight);
    qsmall.save( QString("frame.png") );
    /*
        QImage glImg = QGLWidget::convertToGLFormat(qsmall);
        glImg.save( QString("frame.png") );
    */
#ifdef CSV
    JImage img;
    img.readFrom("frame.png");
    texID = img.genTexture();
#endif
}
///////////////////////////////////////////////////////////////////////////////

void JSurfaceMap2DViewer  :: drawPatch()
{

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    qglviewer::Vec  tmppos;
    tmppos[0] =  0.0;
    tmppos[1] =  0.0;
    tmppos[2] =  1.0;
    viewManager->camera()->setPosition(tmppos);

    tmppos[0] = 0.0;
    tmppos[1] = 0.0;
    tmppos[2] = -1.0;
    viewManager->camera()->setViewDirection(tmppos);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(0.95, 0.95, 0.95);
    glBegin(GL_QUADS);
    glVertex2f(-1.0, -1.0);
    glVertex2f( 1.0, -1.0);
    glVertex2f( 1.0,  1.0);
    glVertex2f(-1.0,  1.0);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(1.0, 0.0, 0.0);
    glLineWidth(3.0);
    glBegin(GL_QUADS);
    glVertex3f(-0.99, -0.99, 0.0001);
    glVertex3f( 0.99, -0.99, 0.0001);
    glVertex3f( 0.99,  0.99, 0.0001);
    glVertex3f(-0.99,  0.99, 0.0001);
    glEnd();

    glLineWidth(1.0);
    glDisable(GL_CULL_FACE);

    JColor color;
    if( patch ) {
        patch->getAttribute("Color", color);
        glScalef( 1.75, 1.75, 1.75);
        glLineWidth(1.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glColor3fv(&color[0]);
        size_t numfaces = patch->getSize(2);
        glBegin(GL_TRIANGLES);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = patch->getFaceAt(i);
            for( int j = 0; j < 3; j++) {
                const Point3D &xyz = face->getNodeAt(j)->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], 0.0001 );
            }
        }
        glEnd();

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(0.0, 0.0, 0.0);
        glBegin(GL_TRIANGLES);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = patch->getFaceAt(i);
            for( int j = 0; j < 3; j++) {
                const Point3D &xyz = face->getNodeAt(j)->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], 0.0002 );
            }
        }
        glEnd();

        glColor3f(1.0, 0.0, 1.0);
        glPointSize(10);
        size_t numnodes = patch->getSize(0);
        glBegin(GL_POINTS);
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx= patch->getNodeAt(i);
            if( vtx->hasAttribute("PartitionCorner") ) {
                const Point3D &xyz = vtx->getXYZCoords();
                glVertex3f(xyz[0], xyz[1], 0.0005 );
            }
        }
        glEnd();
        glPointSize(1);
    }
    glEnable(GL_CULL_FACE);
    glScalef( 1.0, 1.0, 1.0);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}
///////////////////////////////////////////////////////////////////////////////////////
void JSurfaceMap2DViewer  :: drawImage()
{
    if( texID == 0) return;

    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_BLEND);
    glBindTexture(GL_TEXTURE_2D, texID);
    double U = 1.0;
    double V = 1.0;
    int winHeight = viewManager->height();
    int x = projCenter[0] - 0.5*projWidth;
    int y = projCenter[1] - 0.5*projHeight;

    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
    viewManager->startScreenCoordinatesSystem(true);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0);
    glVertex2f(x, y);

    glTexCoord2f(U,   0.0);
    glVertex2f(x + projWidth, y);

    glTexCoord2f(U,   V  );
    glVertex2f(x+projWidth, y+projHeight);
    glTexCoord2f(0.0, V  );
    glVertex2f(x, y + projHeight);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(1.0, 1.0, 0.0);
    glLineWidth(3.0);
    glBegin(GL_QUADS);
    glVertex2f(x,y);
    glVertex2f(x + projWidth, y);
    glVertex2f(x+projWidth, y+projHeight);
    glVertex2f(x, y + projHeight);
    glEnd();

    viewManager->stopScreenCoordinatesSystem();
}
///////////////////////////////////////////////////////////////////////////////////////
void JSurfaceMap2DViewer  :: init()
{
    if( patch == nullptr) return;

    //Backup ...
    qglviewer::Vec  prevCameraPos = viewManager->camera()->position();
    qglviewer::Vec  prevCameraDir = viewManager->camera()->viewDirection();

    // New Setting ....
    int winHeight = viewManager->height();

    int viewport[4];
    int scissor[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    glGetIntegerv(GL_SCISSOR_BOX,scissor);

    int xp = projCenter[0] - 0.5*projWidth;
    int yp = projCenter[1] - 0.5*projHeight;
    glViewport(xp,yp, projWidth, projHeight);
    glScissor(xp,yp, projWidth, projHeight);
    glClear(GL_DEPTH_BUFFER_BIT);

    glDisable(GL_LIGHTING);
    // Draw the stuff into new viewport ...
    viewManager->camera()->setType(Camera::ORTHOGRAPHIC);
    drawPatch();
    getImage();
    viewManager->camera()->setType(Camera::PERSPECTIVE);
    // Revert back ...
    glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);

    viewManager->camera()->setPosition( prevCameraPos);
    viewManager->camera()->setViewDirection(prevCameraDir);
    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JSurfaceMap2DViewer  :: drawChart()
{
    int numPatches = mapPatches.size();

    int nx = sqrt(numPatches) + 1;
    int ny = nx;

    // New Setting ....
    int winHeight = viewManager->height();

    int viewport[4];
    int scissor[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    glGetIntegerv(GL_SCISSOR_BOX,scissor);

    int w = 0.5*viewManager->width();
    int h = 0.5*viewManager->height();
    int d = max(w,h)+5;

    glViewport(10, 10, d, d);
    glScissor( 10, 10, d, d);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_LIGHTING);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    double xmin, ymin, xmax, ymax;
    xmin = -0.01;
    ymin = -0.01;
    xmax = nx+0.01;
    ymax = ny+0.01;

    glOrtho(xmin, xmax, ymin, ymax, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    qglviewer::Vec  tmppos;
    tmppos[0] =  0.0;
    tmppos[1] =  0.0;
    tmppos[2] =  1.0;
    viewManager->camera()->setPosition(tmppos);

    tmppos[0] = 0.0;
    tmppos[1] = 0.0;
    tmppos[2] = -1.0;
    viewManager->camera()->setViewDirection(tmppos);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(0.95, 0.95, 0.95);
    glBegin(GL_QUADS);
    glVertex2f( xmin, ymin);
    glVertex2f( xmax, ymin);
    glVertex2f( xmax, ymax);
    glVertex2f( xmin, ymax);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(1.0, 0.0, 0.0);
    glLineWidth(3.0);
    glBegin(GL_QUADS);
    glVertex2f( xmin, ymin);
    glVertex2f( xmax, ymin);
    glVertex2f( xmax, ymax);
    glVertex2f( xmin, ymax);
    glEnd();

    glLineWidth(1.0);
    glDisable(GL_CULL_FACE);
    JColor color;

    int ncount = 0;
    for( int iy = 0; iy < ny; iy++) {
        for( int ix = 0; ix < nx; ix++) {
            int id = iy*nx + ix;
            if( mapPatches.find(id) == mapPatches.end() ) continue;
            patch = mapPatches[id];
            if( patch )  {
                patch->getAttribute("Color", color);
                glPushMatrix();
                glTranslatef(1.0*ix+0.5, 1.0*iy+0.5, 0.0);

                glScalef( 0.95, 0.95, 0.95);
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glColor3fv(&color[0]);
                size_t numfaces = patch->getSize(2);
                glBegin(GL_TRIANGLES);
                for( size_t i = 0; i < numfaces; i++) {
                    const JFacePtr &face = patch->getFaceAt(i);
                    for( int j = 0; j < 3; j++) {
                        const Point3D &xyz = face->getNodeAt(j)->getXYZCoords();
                        glVertex3f(xyz[0], xyz[1], 0.0001 );
                    }
                }
                glEnd();
                glLineWidth(1.0);
                glColor3f(1.0, 0.0, 1.0);
                glPointSize(2);
                size_t numnodes = patch->getSize(0);
                glBegin(GL_POINTS);
                for( size_t i = 0; i < numnodes; i++) {
                    const JNodePtr &vtx= patch->getNodeAt(i);
                    if( vtx->hasAttribute("PartitionCorner") ) {
                        const Point3D &xyz = vtx->getXYZCoords();
                        glVertex3f(xyz[0], xyz[1], 0.0002 );
                    }
                }
                glEnd();
                glPointSize(1);
                glScalef( 1.0, 1.0, 1.0);
                glPopMatrix();

                /*
                                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                                glColor3f(0.0, 0.0, 0.0);
                                glBegin(GL_TRIANGLES);
                                for( size_t i = 0; i < numfaces; i++) {
                                    const JFacePtr &face = patch->getFaceAt(i);
                                    for( int j = 0; j < 3; j++) {
                                        const Point3D &xyz = face->getNodeAt(j)->getXYZCoords();
                                        glVertex3f(xyz[0], xyz[1], 0.0002 );
                                    }
                                }
                                glEnd();
                */
                glLineWidth(5.0);
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                glColor3f( 0.0, 0.0, 0.0);
                glBegin(GL_QUADS);
                glVertex3f( 1.0*ix,     1.0*iy,   0.0003);
                glVertex3f( 1.0*ix +1 , 1.0*iy,   0.0003);
                glVertex3f( 1.0*ix +1 , 1.0*iy+1, 0.0003);
                glVertex3f( 1.0*ix ,    1.0*iy+1, 0.0003);
                glEnd();
                glLineWidth(1.0);
            }
        }
    }
    glEnable(GL_CULL_FACE);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    // Revert back ...
    glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
}
///////////////////////////////////////////////////////////////////////////////
void JSurfaceMap2DViewer  :: draw()
{
    if( viewManager == nullptr) return;

    if( showChart ) {
        drawChart();
        return;
    }

    if( texID == 0) init();

    // New Setting ....
    int winHeight = viewManager->height();

    int viewport[4];
    int scissor[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    glGetIntegerv(GL_SCISSOR_BOX,scissor);

    int xp = projCenter[0] - 0.5*projWidth;
    int yp = projCenter[1] - 0.5*projHeight;
    glViewport(xp,yp, projWidth, projHeight);
    glScissor(xp,yp, projWidth, projHeight);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_LIGHTING);

    drawImage();

    // Revert back ...
    glScissor(scissor[0],scissor[1],scissor[2],scissor[3]);
    glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
}

