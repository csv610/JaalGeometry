#include "AlphaMSTQuadMeshDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JAlphaMSTQuadMeshDialog :: JAlphaMSTQuadMeshDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    alphaMST.reset( new JAlphaMSTQuadMesh);
    ballRadiusLineEdit->setText(QString::number(0.0001));
}

///////////////////////////////////////////////////////////////////////////////

JAlphaMSTQuadMeshDialog :: ~JAlphaMSTQuadMeshDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMeshDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    entityPicker = meshViewer->getEntityPicker();
    viewManager->attach( this );
    entityPicker->setMode(1);
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: showEvent( QShowEvent *e)
{
    if( mesh == nullptr) return;
    entityPicker->setMesh(mesh);
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 0;
    entityPicker->setMode(1);
    resetColor();
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: enablePicking()
{
    if( mesh ) {
        entityPicker->setMesh(mesh);
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = 0;
        entityPicker->setMode(1);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: setMesh( const JMeshPtr &m)
{
    mesh.reset();;

    if( m == nullptr ) return ;

    mesh = m;
    int dim = mesh->getTopology()->getDimension();
    int etype = mesh->getTopology()->getElementsType(2);
    if( dim != 2  || etype != JFace::QUADRILATERAL) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("AlphaMSTQuad mesh requires Quadrilateral mesh");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
    entityPicker->setMesh(mesh);

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    alphaMST->setMesh(mesh);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 0;

    int nCount = mesh->getTopology()->getNumIrregularNodes();
    numMeshSingularitiesLineEdit->setText(QString::number(nCount));

    size_t numfaces = mesh->getActiveSize(2);
    numFacesLineEdit->setText( QString::number(numfaces));

    mesh->pruneFaces();
    faceColors.resize(numfaces);
    JFaceRenderPtr fAttrib;
    numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", fAttrib);
            fAttrib->display = 0;
            faceColors[i] = fAttrib->color;
        }
    }
    double elen = mesh->getGeometry()->getMeanEdgeLength();
    boundaryEdgeLengthLineEdit->setText( QString::number(elen));

    expectedEdgeLengthLineEdit->setText( QString::number(elen));

    JMeshQuality mq;
    mq.setMesh(mesh);
    vector<double> q =  mq.getEdgesQuality(JMeshQuality::EDGE_LENGTH,
                                           JMeshEntity::INTERNAL_ENTITY);
    double sum  = boost::accumulate(q, 0.0);
    double avg  = sum/q.size();
    internalEdgeLengthLineEdit->setText( QString::number(avg));

    mesh->buildRelations(0,2);
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: getNewPatch( const Point3D &center)
{
    JCircle circle;
    if( fixedRadiusCheckBox->isChecked() ) {
        double r = circleRadiusSpinBox->value();
        circle.setCenter(center);
        circle.setRadius(r);
        alphaMST->setCircle(circle);
    } else
        alphaMST->setCenter(center);

    alphaMST->buildPatch();

    circle = alphaMST->getCircle();
    circleRadiusSpinBox->setValue(circle.getRadius());

    int nsingular = alphaMST->getNumSingularities();
    numPatchSingularitiesLineEdit->setText( QString::number(nsingular));

    patchNodes = alphaMST->getBoundaryNodes();
    boundNodesLineEdit->setText( QString::number(patchNodes.size()));

    mesh->enumerate(0);
    JColor grey = JEntityColor::getColor("LigthGray");
    grey[0] = 0.9;
    grey[1] = 0.9;
    grey[2] = 0.9;
    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", fAttrib);
            fAttrib->color   = grey;
            fAttrib->display = 0;
        }
    }

    JColor black = JEntityColor::getColor("Black");
    JEdgeRenderPtr eAttrib;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", eAttrib);
            eAttrib->color = black;
            eAttrib->lineWidth = 1.0;
        }
    }

    JColor pink = JEntityColor::getColor("Pink");
    JFaceSequence patchfaces = alphaMST->getFaces();
    for( const JFacePtr &face : patchfaces) {
        face->getAttribute("Render", fAttrib);
        fAttrib->color = pink;
        fAttrib->display = 1;
    }
    numPatchFacesLineEdit->setText(QString::number(patchfaces.size()));

    JColor blue = JEntityColor::getColor("Blue");
    JEdgeSequence patchedges = alphaMST->getBoundaryEdges();
    for( const JEdgePtr &edge : patchedges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->color = blue;
        eAttrib->lineWidth = 3.0;
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: mouseReleaseEvent(QMouseEvent *)
{
    if( entityPicker == nullptr || mesh == nullptr ) return;

    JNodeSequence nodeSeq = entityPicker->getPickedNodes();
    if( nodeSeq.empty() ) return;

    const Point3D &center = nodeSeq[0]->getXYZCoords();
    getNewPatch(center);
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: openMedialAxisDialog()
{
    if( medialAxisDialog == nullptr)
        medialAxisDialog.reset( new JMedialAxis2DDialog(this));

    medialAxisDialog->setViewManager( viewManager );
    medialAxisDialog->setMesh(mesh);
    medialAxisDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: getNewCircle()
{
    if( medialCircles.empty() ) {
        JCircleCover cc;
        cc.setMesh(mesh);
        medialCircles = cc.getCircles();
        size_t numnodes = medialCircles.size();
        currCircleSpinBox->setMaximum(numnodes-1);
    }

    int id = currCircleSpinBox->value();
    const Point3D center = medialCircles[id].getCenter();
    getNewPatch( center );
    currCircleSpinBox->setValue(id+1);
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: resetColor()
{
    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", fAttrib);
            fAttrib->color = faceColors[i];
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMeshDialog :: remesh()
{
    if( alphaMST->isEmpty() ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    int id = 0;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->setAttribute("Partition", id);
        }
    }

//  double alpha = adaptFactorSpinBox->value();
//  alphaMST->setAdaptationFactor(alpha);

    QString qstr = expectedEdgeLengthLineEdit->text();
    double  elen = qstr.toDouble();
    alphaMST->setExpectedEdgeLength(elen);

    alphaMST->remeshPatch();
    meshViewer->updateBuffers(mesh);  // Important to call just after remesh ....

    mesh->enumerate(0);
    mesh->enumerate(1);
    mesh->enumerate(2);
    JColor colors[6];
    colors[0] =  JEntityColor::getColor("White");
    colors[1] =  JEntityColor::getColor("Red");
    colors[2] =  JEntityColor::getColor("Green");
    colors[3] =  JEntityColor::getColor("Blue");
    colors[4] =  JEntityColor::getColor("Magenta");
    colors[5] =  JEntityColor::getColor("Yellow");

    mesh->pruneFaces();
    JFaceRenderPtr fAttrib;
    numfaces = mesh->getSize(2);
    faceColors.resize(numfaces);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Partition", id);
            face->getAttribute("Render", fAttrib);
            fAttrib->color = colors[id];
            faceColors[i]  = fAttrib->color;
        }
    }
    meshViewer->updateBuffers(mesh);
    numFacesLineEdit->setText( QString::number(numfaces));

    const JFaceSequence &newFaces = alphaMST->getNewFaces();
    numPatchFacesLineEdit->setText(QString::number(newFaces.size()));

    int nsingular = alphaMST->getNumSingularities();
    numPatchSingularitiesLineEdit->setText( QString::number(nsingular));

    int nCount = mesh->getTopology()->getNumIrregularNodes();
    numMeshSingularitiesLineEdit->setText(QString::number(nCount));
    showSingularities();

    JMeshQuality mq;
    mq.setMesh(mesh);
    vector<double> q =  mq.getEdgesQuality(JMeshQuality::EDGE_LENGTH,
                                           JMeshEntity::INTERNAL_ENTITY);
    double sum  = boost::accumulate(q, 0.0);
    double avg  = sum/q.size();
    internalEdgeLengthLineEdit->setText( QString::number(avg));
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: showAllCircles()
{
    bool val = showAllCirclesCheckBox->isChecked();
    if( !val ) {
        if( coverMesh ) {
            meshViewer->removeObject(coverMesh);
            coverMesh.reset();
        }
        meshViewer->refreshDisplay();
        return;
    }

    JCircleCover cc;
    cc.setMesh(mesh);
    coverMesh = cc.getCoverMesh();

    size_t numnodes = coverMesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = coverMesh->getNodeAt(i);
        Point3D  xyz = vtx->getXYZCoords();
        xyz[2] = 0.001;
        vtx->setXYZCoords(xyz);
    }
    meshViewer->addObject(coverMesh);

    JColor color = JEntityColor::getColor("LightBlue");
    color[0] = 0.5;
    color[1] = 0.5;
    color[2] = 1.0;
    color[3] = 0.2;

    JFaceRenderPtr fAttrib;
    size_t numfaces = coverMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = coverMesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->color = color;
    }

    JMeshRenderPtr mrender;
    coverMesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 1;
    mrender->transparent      = 1;
    mrender->transparencyMethod   = JRender::TRANSPARENT_BLEND;
    meshViewer->updateBuffers(coverMesh);
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog ::  showSingularities()
{
    if( mesh == nullptr) return;

    JNodeRenderPtr nAttrib;
    size_t numnodes = mesh->getSize(0);
    JColor defaultColor  = JEntityColor::getColor("Green");
    JColor degree3Color  = JEntityColor::getColor("Red");
    JColor degree5Color  = JEntityColor::getColor("Blue");

    bool val = showSingularitiesCheckBox->isChecked();
    if( !val ) {
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            int err = vtx->getAttribute("Render", nAttrib);
            nAttrib->display =  1;
            nAttrib->glyph   =  0;
            nAttrib->scale   =  1.0;
            nAttrib->color   =  defaultColor;
        }
        mstMesher.reset();
        meshViewer->updateBuffers(mesh);
        return;
    }

    if( mstMesher == nullptr)
        mstMesher.reset( new JMSTQuadMesher);

    mstMesher->setMesh(mesh);

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        int err = vtx->getAttribute("Render", nAttrib);
        assert(!err);
        nAttrib->display =  1;
        nAttrib->glyph   =  0;
        nAttrib->scale   =  1.0;
        nAttrib->color   =  defaultColor;
    }

    JNodeSequence irrnodes = mstMesher->getSingularNodes();

    QString qstr = ballRadiusLineEdit->text();
    double rad  = qstr.toDouble();
    cout << "Radius " << endl;
    for( const JNodePtr &vtx : irrnodes) {
        int err = vtx->getAttribute("Render", nAttrib);
        assert( !err);

        if( vtx->getNumRelations(2) < 4) {
            nAttrib->display =  1;
            nAttrib->glyph   =  1;
            nAttrib->color   =  degree3Color;
            nAttrib->ballRadius =  rad;
        }

        if( vtx->getNumRelations(2) > 4) {
            nAttrib->display =  1;
            nAttrib->glyph   =  1;
            nAttrib->color   =  degree5Color;
            nAttrib->ballRadius =  rad;
        }
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: openMeshOptDialog()
{
    if( meshOptDialog == nullptr)
        meshOptDialog.reset( new JMeshOptDialog(this) );

    mesh->pruneFaces();
    size_t numfaces = mesh->getSize(2);
    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        faceColors[i] = fAttrib->color;
    }
    clearAll();
    meshViewer->setCurrentMesh(mesh);
    meshOptDialog->setViewManager( viewManager );
    meshOptDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMeshDialog :: sweepAll()
{
    if( medialCircles.empty() ) {
        JCircleCover cc;
        cc.setMesh(mesh);
        medialCircles = cc.getCircles();
        size_t numnodes = medialCircles.size();
        currCircleSpinBox->setMaximum(numnodes-1);
    }

    for( int i = 0; i < medialCircles.size(); i++) {
        qApp->processEvents();
        getNewCircle();
        remesh();
    }
}
///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: clearAll()
{
    if( mesh == nullptr) return;

    JColor grey = JEntityColor::getColor("LigthGray");
    grey[0] = 0.9;
    grey[1] = 0.9;
    grey[2] = 0.9;
    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", fAttrib);
            fAttrib->color = grey;
            fAttrib->display = 1;
        }
    }

    JColor black = JEntityColor::getColor("Black");
    JEdgeRenderPtr eAttrib;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", eAttrib);
            eAttrib->color     = black;
            eAttrib->lineWidth = 1.0;
            eAttrib->display   = 1;
        }
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMeshDialog :: openPolyMesherDialog()
{
    if( patchNodes.empty() ) return;
/*

    if( polyMesherDialog == nullptr)
        polyMesherDialog.reset( new JPolygonQuadMesherDialog(this));

    polyMesherDialog->setViewManager( viewManager );
    polyMesherDialog->setMesh(mesh);
//  polyMesherDialog->setPatchCenters(patchNodes);
    polyMesherDialog->show();
    this->hide();
*/
}
///////////////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMeshDialog :: hidePatch()
{
    if(alphaMST == nullptr) return;

    bool val = hidePatchElementsCheckBox->isChecked();

    JFaceRenderPtr fAttrib;
    JFaceSequence patchfaces = alphaMST->getFaces();
    for( const JFacePtr &face : patchfaces) {
        face->getAttribute("Render", fAttrib);
        fAttrib->display = val;
    }

    JEdgeSequence iedges;
    JNodeSequence inodes;
    mesh->getTopology()->getSubmeshInternal(patchfaces, iedges, inodes);

    JEdgeRenderPtr eAttrib;
    for( const JEdgePtr &edge : iedges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->display = val;
    }

    JNodeRenderPtr vAttrib;
    for( const JNodePtr &vtx : inodes) {
        vtx->getAttribute("Render", vAttrib);
        vAttrib->display = val;
    }

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMeshDialog :: showCorners()
{
}
///////////////////////////////////////////////////////////////////////////////
void JAlphaMSTQuadMeshDialog :: closeDialog()
{
    clearAll();

    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = -1;
    }
    alphaMST->clear();

    viewManager->detach( this );

    parentWidget()->show();
    close();
}

///////////////////////////////////////////////////////////////////////////////

void JAlphaMSTQuadMeshDialog :: makeConnections()
{
    LineEdit( ballRadiusLineEdit, [=] {showSingularities();});
    CheckBox( showAllCirclesCheckBox, [=]{showAllCircles();});
    CheckBox( hidePatchElementsCheckBox, [=]{hidePatch();});
    CheckBox( showPatchCornersCheckBox, [=]{showCorners();} );
    CheckBox( showSingularitiesCheckBox, [=]{showSingularities();});

    PushButton( pickPushButton,  [=]{enablePicking();});
    PushButton( medialAxisPushButton,  [=]{openMedialAxisDialog(); });
    PushButton( nextCirclePushButton,  [=]{getNewCircle(); });
    PushButton( remeshPushButton,  [=]{remesh(); });
    PushButton( meshOptPushButton,  [=]{openMeshOptDialog(); });
    PushButton( clearPushButton,  [=]{clearAll(); });
    PushButton( removeBoundSingularitiesPushButton, [=]{openPolyMesherDialog(); });
    PushButton( patchColorPushButton,  [=]{resetColor(); });
    PushButton( sweepAllPushButton, [=]{sweepAll();});
    PushButton( closePushButton, [=]{closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
