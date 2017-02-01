#include "MSTQuadMesherDialog.hpp"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace Jaal;

///////////////////////////////////////////////////////////////////////////////

JMSTQuadMesherDialog :: JMSTQuadMesherDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
    numSmoothSpinBox->setValue(1);
    numQuadsWantedLineEdit->setText( QString::number(1) );
    minSingularNodesLineEdit->setText( QString::number(3) );
}

///////////////////////////////////////////////////////////////////////////////
JMSTQuadMesherDialog :: ~JMSTQuadMesherDialog()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: makeConnections()
{
    ComboBox( showSingularitiesComboBox, [=] {displaySingularNodes();});
    SpinBoxd( radiusSingularitySpinBox, [=] {displaySingularNodes();});

    PushButton( applyCanonicalPushButton, [=] {unittest(); });
    PushButton( smoothMeshPushButton,  [=] {smoothMesh(); });
    PushButton( untanglePushButton,  [=] {untangleMesh();});
    PushButton( repeatPushButton,  [=] {repeatSearch(); });
    PushButton( clearPushButton, [=] {clearAll();});
    PushButton( checkPointPushButton, [=] {checkPoint();});
    PushButton( patchRemeshPushButton, [=] {remeshPatch(); });
    PushButton( randomPushButton, [=] {randomSegments();});
    PushButton( savePatternPushButton, [=] {savePattern(); });
    PushButton( recoveryPushButton, [=] {checkPointRecovery();});
    PushButton( closePushButton, [=] {closeDialog(); });
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c == nullptr)
        c = JMeshViewer::registerComponent(viewManager);
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    nodePicker = meshViewer->getEntityPicker();
    viewManager->attach( this );

    viewManager->getLights()->Switch(0);
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    int topDim = mesh->getTopology()->getDimension();

    if( topDim != 2 ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Quadmesher in only for 2D mesh");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        mesh = nullptr;
        return;
    }

    int elemType = mesh->getTopology()->getElementsType(2);
    if( elemType != JFace::QUADRILATERAL ) {
        cout << "ElemType " << elemType << endl;
        cout << "Debug: Current mesh " << meshViewer->getCurrentMesh()->getName() << endl;
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Quadmesher works only for pure quadmesh");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        mesh = nullptr;
        return;
    }
    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    mesh->buildRelations(0,2);

    mesh->getTopology()->searchBoundary();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 0;
    nodePicker->setMesh( mesh );

    mstMesher.setMesh(mesh);

    JNodeSequence irrnodes = mstMesher.getSingularNodes();
    numSearchSpinBox->setValue(irrnodes.size());

    displaySingularNodes();
    displayMesh();
/*

    double elen = mesh->getGeometry()->getMeanBoundaryEdgeLength();
    double area = mesh->getGeometry()->getSurfaceArea();
    int    ntarget = area/(elen*elen);
    numTargetFacesLineEdit->setText( QString::number(ntarget) );
    numCurrentFacesLineEdit->setText( QString::number(numfaces) );

*/
}

///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: displayMesh()
{
    JFaceRenderPtr fAttrib;
    JColor whitish;
    whitish[0] = 0.99;
    whitish[1] = 0.99;
    whitish[2] = 0.99;
    whitish[3] = 1.0;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
        fAttrib->color   = whitish;
    }
    meshViewer->updateBuffers(mesh);

}

void JMSTQuadMesherDialog :: setPatchCenters( const JNodeSequence &v)
{
    specifiedPatchCenters = v;
    int niter = v.size();
    numSearchSpinBox->setMaximum(niter);
    numSearchSpinBox->setValue(0);
}
///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: displaySingularNodes()
{
    if( mesh == nullptr) return;

    double radius = radiusSingularitySpinBox->value();
    if( radius < 0.0) return;

    JNodeRenderPtr nAttrib;
    size_t numnodes = mesh->getSize(0);

    JColor  defaultColor, degree3Color, degree5Color;
    defaultColor = JEntityColor::getColor("Green");

    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        int err = vtx->getAttribute("Render", nAttrib);
        assert(!err);
        nAttrib->display =  1;
        nAttrib->glyph   =  0;
        nAttrib->scale   =  1.0;
        nAttrib->color   =  defaultColor;
    }

    JNodeSequence irrnodes = mstMesher.getSingularNodes();
    totalSingularNodesLineEdit->setText( QString::number(irrnodes.size()) );

    degree3Color  = JEntityColor::getColor("Red");
    degree5Color  = JEntityColor::getColor("Blue");

    QString qstr = showSingularitiesComboBox->currentText();
    std::string str = StdString(qstr);

    bool showLow  = 1;
    bool showHigh = 1;

    int lowDegree = 3;
    if( str == "Lowest Degree") {
        for( const JNodePtr &vtx : irrnodes)
            lowDegree = min(lowDegree, vtx->getNumRelations(2));
        showHigh  = 0;
    }

    int highDegree = 5;
    if( str == "Highest Degree") {
        for( const JNodePtr &vtx : irrnodes)
            highDegree = max(highDegree, vtx->getNumRelations(2));
        showLow  = 0;
    }

    for( const JNodePtr &vtx : irrnodes) {
        int err = vtx->getAttribute("Render", nAttrib);
        assert( !err);

        if( showLow) {
            if( vtx->getNumRelations(2) <= lowDegree) {
                nAttrib->display =  1;
                nAttrib->glyph   =  1;
                nAttrib->ballRadius  =  radius;
                nAttrib->color   =  degree3Color;
            }
        }

        if( showHigh) {
            if( vtx->getNumRelations(2) >= highDegree) {
                nAttrib->display =  1;
                nAttrib->glyph   =  1;
                nAttrib->ballRadius  =  radius;
                nAttrib->color   =  degree5Color;
            }
        }
    }

    if(defectivePatch) {
        JColor highlightColor;
        highlightColor[0] = 0.0;
        highlightColor[1] = 1.0;
        highlightColor[2] = 0.0;
        highlightColor[3] = 1.0;

        JNodeSequence corners = defectivePatch->getCorners();
        for( const JNodePtr &node : corners) {
            node->getAttribute("Render", nAttrib);
            nAttrib->color      =  highlightColor;
            nAttrib->ballRadius =  radius;
            nAttrib->glyph      =  1;
            nAttrib->display    =  1;
            nAttrib->scale      =  1.2;
        }
    }
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: openSingularityGraphDialog()
{
    if( meshSingularityGraphDialog == nullptr)
        meshSingularityGraphDialog.reset( new JMeshSingularityGraphDialog(this));
    meshSingularityGraphDialog->setViewManager(viewManager);
    meshSingularityGraphDialog->setMesh(mesh);
    meshSingularityGraphDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: getDefectAt()
{
    getDefectivePatch();
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: clearAll()
{
    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display = 1;
    }

    JColor  highlightColor;
    highlightColor[0] = 0.3;
    highlightColor[1] = 0.3;
    highlightColor[2] = 0.3;
    highlightColor[3] = 1.0;

    size_t numedges = mesh->getSize(1);
    JEdgeRenderPtr eAttrib;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        eAttrib->color     =  highlightColor;
        eAttrib->lineWidth =  1;
        eAttrib->display   =  1;
    }

    size_t numnodes = mesh->getSize(0);
    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &node = mesh->getNodeAt(i);
        node->getAttribute("Render", nAttrib);
        nAttrib->display = 0;
        nAttrib->glyph   = 0;
    }
    defectivePatch.reset();

    displaySingularNodes();
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: displayPatch()
{
    if( defectivePatch == nullptr) return;

    int n = defectivePatch->getNumSingularNodes(0);
    numSingularNodesInPatchLineEdit->setText( QString::number(n) );

    vector<int> segments = defectivePatch->getSegments();
    int nsum = 0;
    spinBox1->setValue(0);
    spinBox2->setValue(0);
    spinBox3->setValue(0);
    spinBox4->setValue(0);
    spinBox5->setValue(0);
    spinBox6->setValue(0);

    for( int i = 0; i < segments.size(); i++) {
        int nval = segments[i];
        if( i == 0) spinBox1->setValue(nval);
        if( i == 1) spinBox2->setValue(nval);
        if( i == 2) spinBox3->setValue(nval);
        if( i == 3) spinBox4->setValue(nval);
        if( i == 4) spinBox5->setValue(nval);
        nsum += segments[i];
    }
    totalBoundNodesLineEdit->setText( QString::number(nsum) );

    JNodeSequence irrnodes = mstMesher.getSingularNodes();
    totalSingularNodesLineEdit->setText( QString::number(irrnodes.size()) );

    JColor  highlightColor;

    JFaceRenderPtr fAttrib;
    highlightColor[0] = 0.8;
    highlightColor[1] = 0.8;
    highlightColor[2] = 0.8;
    highlightColor[3] = 1.0;
    JFaceSequence faces = defectivePatch->getFaces();
    for( const JFacePtr &face : faces) {
        face->getAttribute("Render", fAttrib);
        fAttrib->display =  1;
        fAttrib->color   = highlightColor;
    }

    numQuadsWantedLineEdit->setText( QString::number(faces.size()) );

    highlightColor[0] = 1.0;
    highlightColor[1] = 0.0;
    highlightColor[2] = 0.0;
    highlightColor[3] = 1.0;

    JEdgeRenderPtr eAttrib;
    JEdgeSequence boundedges = defectivePatch->getBoundEdges();
    for( const JEdgePtr &edge : boundedges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->color     =  highlightColor;
        eAttrib->lineWidth  =  4;
    }

    if( changeSceneCenterCheckBox->isChecked() ) {
        JMeshAffineTransform affineTrans;
        JNodePtr seed = defectivePatch->getSeedNode();
        if( seed ) {
            affineTrans.setMesh(mesh);
            const Point3D &p3d = seed->getXYZCoords();
            affineTrans.translate( -p3d[0], -p3d[1], -p3d[2] );
        }
    }

    displaySingularNodes();
}

/////////////////////////////////////////////////////////////////////////////////

int JMSTQuadMesherDialog :: getDefectivePatch()
{
    clearAll();

    JWaitCursor waitCursor;
    waitCursor.start();

    if( !specifiedPatchCenters.empty() ) {
        cout << "Get Specified Defective Patch " << endl;
        JNodePtr seed = specifiedPatchCenters.back();
        specifiedPatchCenters.pop_back();
        defectivePatch = mstMesher.getDefectivePatch(seed);
    } else
        defectivePatch = mstMesher.getAnyDefectivePatch();

    waitCursor.stop();

    if( defectivePatch == nullptr) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("No valid interior patch found ");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return 1;
    }

    displayPatch();

    return 0;
}

//////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: showEvent(QShowEvent *)
{
    mylogfile.open("singular.dat", ios::out);
    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = 0;
    }
}
//////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( nodePicker == nullptr || mesh == nullptr ) return;

    JNodeSequence nodeSeq = nodePicker->getPickedNodes();
    if( nodeSeq.empty() ) return;

    if( nodeSeq.size() > 1 ) {
        cout << "Only one nodes must be picked at a time " << endl;
        return;
    }

    /*
        bool val = geodesicPathCheckBox->isChecked();
        mesher.setGeodesicPath(val);
    */

    QString str = minSingularNodesLineEdit->text() ;
    int minSingular  = str.toInt();

    mstMesher.setMinSingularNodes(minSingular);

    clearAll();
    const JNodePtr &vtx = nodeSeq[0];

    defectivePatch = mstMesher.getDefectivePatch(vtx);

    displayPatch();

    nodePicker->clearAll();
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: smoothBoundary()
{
    if( mesh == nullptr) return;

    JMeshContour  jc;
    jc.setMesh( mesh );
    jc.smooth(1);
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: errorMessage()
{
    QMessageBox msg;
    msg.setIcon(QMessageBox::Warning);
    msg.setText(" Patch not quadmeshable");
    msg.setStandardButtons( QMessageBox::Ok);
    int ret = msg.exec();
    if( ret == QMessageBox::Ok ) return;
}

///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: savePattern()
{

    if( pattern.empty() ) return;

    pattern = pattern + ".png";
    viewManager->setSnapshotFormat("PNG");

    QString qstr = QString::fromStdString(pattern);
    //viewManager->setSnapshotFileName( qstr );
    viewManager->saveSnapshot( qstr );

}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: remeshPatch()
{
    static size_t remeshCounter = 0;
    if( defectivePatch == nullptr) return;

    if( !defectivePatch->isValid() ) return;

    bool inverted;
    inverted = inversionCheckBox->isChecked();
    mstMesher.setHandleInvertedElements(inverted);

    JWaitCursor waitCursor;
    waitCursor.start();

    int err = mstMesher.remesh(defectivePatch);
    if( !err) {
        mesh->getGeometry()->setFacesNormal();
        mesh->getGeometry()->setNodesNormal();
        JMeshPtr mtemplate = defectivePatch->getTemplateMesh();
        if( mtemplate ) {
            JFaceRenderPtr fAttrib;
            JColor highlightColor;
            highlightColor[0] = 0.8;
            highlightColor[1] = 0.8;
            highlightColor[2] = 0.8;
            highlightColor[3] = 1.0;
            size_t nfaces = mtemplate->getSize(2);
            for( size_t i = 0; i < nfaces; i++) {
                const JFacePtr &face = mtemplate->getFaceAt(i);
                fAttrib.reset( new JFaceRender );
                fAttrib->display =  1;
                fAttrib->color   = highlightColor;
                face->setAttribute("Render", fAttrib);
            }
            mesh->enumerate(0);
        }
        meshViewer->updateBuffers(mesh);
        JNodeSequence irrnodes = mstMesher.getSingularNodes();
        totalSingularNodesLineEdit->setText( QString::number(irrnodes.size()) );

        displaySingularNodes();
        defectivePatch->clear();
        pickedFaces.clear();

        remeshCounter++;
        mylogfile << remeshCounter << " " << irrnodes.size() << endl;
    }

    size_t numfaces = mesh->getActiveSize(2);
    numCurrentFacesLineEdit->setText( QString::number(numfaces) );
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: smoothMesh()
{
    if( mesh == nullptr) return;

    int niter  = numSmoothSpinBox->value();

    JWaitCursor waitCursor;
    waitCursor.start();

    JLloydMeshOptimizer mopt;
    mopt.setMesh(mesh);
    mopt.setNumIterations(niter);
    mopt.smoothAll();

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: untangleMesh()
{
    if( mesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    meshViewer->getFaceDraw()->setCulling(1);
    meshViewer->getFaceDraw()->setFrontFace(GL_CCW);
    meshViewer->getFaceDraw()->setLights(0);

    JMeshUntangle untangle;
    untangle.setMesh(mesh);

    JMeshPtr inflated = untangle.getInflatedMesh();
    if( inflated ) {
        inflated->setName("InflatedMesh");
        meshViewer->addObject(inflated);
        inflated->setActiveBit(0);
        untangle.startBackProjection();
        meshViewer->updateBuffers(mesh);
    }
    meshViewer->getFaceDraw()->setCulling(0);
    meshViewer->getFaceDraw()->setFrontFace(GL_CW);
    meshViewer->getFaceDraw()->setLights(1);
}


///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: reorderSingularities()
{
    QString qstr = removalPolicyComboBox->currentText();
    string str   = StdString(qstr);
    if( str == "Random")
        mstMesher.setSingularityRemovalPolicy(0);

    if( str == "BoundaryFirst")
        mstMesher.setSingularityRemovalPolicy(1);

    if( str == "InteriorFirst")
        mstMesher.setSingularityRemovalPolicy(2);
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: repeatSearch()
{
    int niter  = numSearchSpinBox->value();

    JWaitCursor waitCursor;
    waitCursor.start();
    reorderSingularities();

    for( int i = 0; i < niter; i++) {
        numSearchSpinBox->setValue(niter -1 - i );
        qApp->processEvents();
        int err = getDefectivePatch();
        if( err ) return;
        remeshPatch();
    }
    displayMesh();

    JNodeSequence irrnodes = mstMesher.getSingularNodes();
    numSearchSpinBox->setValue(irrnodes.size());
}

///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: getMedialAxis()
{
    if( mesh == nullptr ) return;

    JDelaunayMesh2D delmesh;
    delmesh.setMesh(mesh);
    JMeshPtr medial = delmesh.getMedialAxis();
    meshViewer->addObject(medial);
}
///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: randomSegments()
{
    vector<int> segments(6);
    int maxval = 100;
    int total;
    segments[0] = 0;
    segments[1] = 0;
    segments[2] = 0;
    segments[3] = 0;
    segments[4] = 0;
    segments[5] = 0;

    QString qs  = templateShapeComboBox->currentText();
    string shape = qs.toUtf8().constData();

    if( boundNodesCheckBox->isChecked() ) {
        qs     = totalBoundNodesLineEdit->text();
        maxval = qs.toInt();
        total  = maxval;
        if( shape == "Triangle" ) {
            maxval -= 2;
            segments[0] = JMath::random_value(1, maxval);
            maxval -= segments[0];
            segments[1] = JMath::random_value(1, maxval);
            segments[2] = total - segments[0] - segments[1];
        }

        if( shape == "Quad" ) {
            maxval -= 3;
            segments[0] = JMath::random_value(1, maxval);
            maxval -= segments[0];
            segments[1] = JMath::random_value(1, maxval);
            maxval -= segments[1];
            segments[2] = JMath::random_value(1, maxval);
            segments[3] = total - segments[0] - segments[1] - segments[2];
        }

        if( shape == "Pentagon" ) {
            maxval -= 3;
            segments[0] = JMath::random_value(1, maxval);
            maxval -= segments[0];
            segments[1] = JMath::random_value(1, maxval);
            maxval -= segments[1];
            segments[2] = JMath::random_value(1, maxval);
            maxval -= segments[2];
            segments[3] = JMath::random_value(1, maxval);
            segments[4] = total - segments[0] - segments[1] - segments[2] -segments[3];
        }

        if( shape == "Hexagon" ) {
            maxval -= 3;
            segments[0] = JMath::random_value(1, maxval);
            maxval -= segments[0];
            segments[1] = JMath::random_value(1, maxval);
            maxval -= segments[1];
            segments[2] = JMath::random_value(1, maxval);
            maxval -= segments[2];
            segments[3] = JMath::random_value(1, maxval);
            maxval -= segments[3];
            segments[4] = JMath::random_value(1, maxval);
            segments[5] = total - segments[0] - segments[1] - segments[2] -segments[3] -segments[4];
        }
    } else {
        segments[0] = JMath::random_value(1, maxval);
        segments[1] = JMath::random_value(1, maxval);
        segments[2] = JMath::random_value(1, maxval);
        segments[3] = JMath::random_value(1, maxval);
        segments[4] = JMath::random_value(1, maxval);
        segments[5] = JMath::random_value(1, maxval);
    }

    bool atleastone = 1;
    if( shape == "Triangle" ) {
        segments.resize(3);
        segments[3] = 0.0;
        segments[4] = 0.0;
        segments[5] = 0.0;
        for(int i = 0; i < 3; i++)
            if( segments[i] < 1) atleastone = 0;;
    }

    if( shape == "Quad" ) {
        segments.resize(4);
        segments[4] = 0.0;
        segments[5] = 0.0;
        for(int i = 0; i < 4; i++) if( segments[i] < 1) atleastone = 0;;
    }

    if( shape == "Pentagon" ) {
        segments.resize(5);
        segments[5] = 0.0;
        for(int i = 0; i < 5; i++) if( segments[i] < 1) atleastone = 0;;
    }

    if( shape == "Hexagon" ) {
        segments.resize(6);
        for(int i = 0; i < 6; i++) if( segments[i] < 1) atleastone = 0;;
    }

    int sum = 0;
    for( int i = 0; i < segments.size(); i++)
        sum += segments[i];
    if( sum%2 ) segments[0]++;

    spinBox1->setValue( segments[0] );
    spinBox2->setValue( segments[1] );
    spinBox3->setValue( segments[2] );
    spinBox4->setValue( segments[3] );
    spinBox5->setValue( segments[4] );
    spinBox6->setValue( segments[5] );

    if( atleastone) unittest();
}
///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: checkPoint()
{
    if( mesh == nullptr) return;
    JMeshIO::saveAs(mesh, "checkpoint.off");
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: checkPointRecovery()
{
    JMeshPtr m = JMeshIO::readFile("checkpoint.off");
    if( m == nullptr) return;
    meshViewer->removeObject(mesh);
    meshViewer->addObject(m);
    mesh = m;
}
///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: closeDialog()
{
    if( testmesh ) {
        testmesh->deleteAll();
        meshViewer->removeObject(testmesh);
    }

    if( singularGraph) {
        meshViewer->removeObject(singularGraph);
        singularGraph->deleteAll();
        singularGraph.reset();
    }

    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = -1;
    }

    clearAll();
    viewManager->detach( this );
    parentWidget()->show();
    close();
    mylogfile.close();
}

///////////////////////////////////////////////////////////////////////////////

void JMSTQuadMesherDialog :: unittest()
{
    if( meshViewer == nullptr ) return;

    if( testmesh ) {
        testmesh->deleteAll();
        testmesh->setActiveBit(0);
        meshViewer->removeObject(testmesh);
        assert( testmesh->getSize(0) == 0);
        assert( testmesh->getSize(1) == 0);
        assert( testmesh->getSize(2) == 0);
    }

    vector<int> segments(6);
    segments[0] = spinBox1->value();
    segments[1] = spinBox2->value();
    segments[2] = spinBox3->value();
    segments[3] = spinBox4->value();
    segments[4] = spinBox5->value();
    segments[5] = spinBox6->value();

    JWaitCursor waitCursor;
    waitCursor.start();

    QString qs  = templateShapeComboBox->currentText();
    string shape = qs.toUtf8().constData();

    ostringstream oss;

    if( shape == "Triangle") {
        segments.resize(3);
        testmesh = mstMesher.getTemplate( segments);
        oss << "t";
        for( int i = 0; i < 3; i++) oss << "e" << segments[i];
    }

    if( shape == "Quad") {
        segments.resize(4);
        testmesh = mstMesher.getTemplate( segments);
        oss << "q";
        for( int i = 0; i < 4; i++) oss << "e" << segments[i];
    }

    if( shape == "Pentagon") {
        segments.resize(5);
        testmesh = mstMesher.getTemplate( segments );
        oss << "p";
        for( int i = 0; i < 5; i++) oss << "e" << segments[i];
    }

    if( shape == "Hexagon") {
        segments.resize(6);
        testmesh = mstMesher.getTemplate(segments);
        oss << "h";
        for( int i = 0; i < 6; i++) oss << "e" << segments[i];
    }

    waitCursor.stop();

    if( testmesh == nullptr ) {
        errorMessage();
        return;
    }

    if( testmesh->getSize(2) == 0) {
        errorMessage();
        return;
    }

    JFaceColorPtr faceColor(new JFacePartitionColor);
    JEdgeColorPtr edgeColor(new JEdgePartitionColor);

    int nCount = testmesh->getTopology()->getNumIrregularNodes();
    numSingularNodesLineEdit->setText(QString::number(nCount));

    waitCursor.start();

    meshViewer->addObject( testmesh );

    size_t numQuads = testmesh->getSize(2);
    numQuadsGeneratedLineEdit->setText( QString::number(numQuads) );

    JLloydMeshOptimizer mopt;
    mopt.setMesh(testmesh);
    mopt.setNumIterations(10);
    mopt.smoothAll();

    testmesh->buildRelations(0,2);
    testmesh->getTopology()->searchBoundary();
    testmesh->enumerate(0);

    JColor degree3Color, degree5Color, boundColor, defaultColor;

    degree3Color[0] = 1.0;
    degree3Color[1] = 0.0;
    degree3Color[2] = 0.0;
    degree3Color[3] = 1.0;

    degree5Color[0] = 0.0;
    degree5Color[1] = 0.0;
    degree5Color[2] = 1.0;
    degree5Color[3] = 1.0;

    boundColor[0] = 0.0;
    boundColor[1] = 1.0;
    boundColor[2] = 0.0;
    boundColor[3] = 1.0;

    defaultColor[0] = 0.1;
    defaultColor[1] = 0.1;
    defaultColor[2] = 0.1;
    defaultColor[3] = 1.0;

    size_t numnodes = testmesh->getSize(0);
    JNodeRenderPtr nAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = testmesh->getNodeAt(i);
        if( vtx->isActive() ) {
            int err = vtx->getAttribute("Render", nAttrib);
            assert( !err);
            nAttrib->glyph   =  1;
            nAttrib->display =  1;
            nAttrib->ballRadius =  0.02;
            nAttrib->pointSize  =  1;
            nAttrib->color   =  defaultColor;
            if( vtx->getNumRelations(2) < 4) {
                nAttrib->color   =  degree3Color;
                nAttrib->ballRadius =  0.04;
                nAttrib->glyph   =  1;
                nAttrib->display =  1;
            }

            if( vtx->getNumRelations(2) > 4) {
                nAttrib->color   =  degree5Color;
                nAttrib->ballRadius =  0.04;
                nAttrib->glyph   =  1;
                nAttrib->display =  1;
            }
        }
    }

    defaultColor[0] = 1.0;
    defaultColor[1] = 1.0;
    defaultColor[2] = 1.0;
    defaultColor[3] = 1.0;
    size_t numfaces = testmesh->getSize(2);
    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = testmesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->color = defaultColor;
    }

    meshViewer->updateBuffers(testmesh);
    pattern = oss.str();
}


/*
void JMSTQuadMesherDialog :: getPatch()
{
    if( mesh == nullptr) return;

    size_t numedges = mesh->getSize(1);
    JEdgeRenderPtr attrib;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            int err = edge->getAttribute("Render", attrib);
            if( !err )  {
                attrib->color[0] = 0;
                attrib->color[1] = 0;
                attrib->color[2] = 0;
                attrib->lineWidth = 1.0;
            }
        }
    }

    int id = patchSpinBox->value();
    JMeshPartitioner mp;
    mp.setMesh(mesh);
    patch = mp.getRegion(id);

    mp.getRegionBoundary( id, patchboundary);

    numedges = patchboundary.size();
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = patchboundary[i];
        if( edge->isActive() ) {
            int err = edge->getAttribute("Render", attrib);
            if( !err )  {
                attrib->color[0] = 1;
                attrib->color[1] = 1;
                attrib->color[2] = 1;
                attrib->lineWidth = 2.0;
            }
        }
    }
    if( meshViewer) meshViewer->updateBuffers(mesh);
}
*/

void JMSTQuadMesherDialog :: detectCorners()
{
    if( mesh == nullptr) return;

    /*
        JMeshContour jc;
        jc.setBoundaryOf( mesh );
        jc.setCornerAngle(15);
        JNodeSequence nodes = jc.getCorners();
        Color highlightColor;
        highlightColor[0] = 1.0;
        highlightColor[1] = 0.0;
        highlightColor[2] = 0.0;
        highlightColor[3] = 0.0;

        JNodeRenderPtr nAttrib;

        for( const JNodePtr &vtx : nodes) {
            vtx->getAttribute("Render", nAttrib);
            nAttrib->color   =  highlightColor;
            nAttrib->glyph   =  1;
            nAttrib->display =  1;
        }
        meshViewer->updateBuffers(mesh);
    */
}


/*
void JMSTQuadMesherDialog :: getPath()
{
    djkPath.setMesh(mesh);
    if( pickedNodes.size() < 2 ) return;

    int nsize = pickedNodes.size();
    JNodePtr vsrc = pickedNodes[nsize-2];
    JNodePtr vdst = pickedNodes[nsize-1];
    JNodeSequence  nodes = djkPath.getPath( vsrc, vdst);
}
*/
//////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: addPatch()
{
    /*
        if( pickedNodes.size() < 3 ) return;

        for( const JNodePtr &vtx : pickedNodes) {
            if( singularNodes.find(vtx) == singularNodes.end() ) {
                singularNodes.insert(vtx);
                singularGraph->addObject(vtx);
            }
        }

        JFacePtr f = JFace::newObject(pickedNodes);
        mesh->addObject(f);
        singularGraph->addObject(f);
        Vec3F normal;
        normal[0] = 0.0;
        normal[1] = 0.0;
        normal[2] = 1.0;
        f->setAttribute("Normal", normal);

        mesh->enumerate(1);
        mesh->enumerate(2);
        meshViewer->updateBuffers(mesh);
    */

    pickedNodes.clear();
    nodePicker->clearAll();
}
///////////////////////////////////////////////////////////////////////////////
void JMSTQuadMesherDialog :: simplifyBoundary()
{
    if( mesh == nullptr) return;

    cout << "Exit Now " << endl;
    exit(0);

    /*
        JMeshContour jc;
        jc.setMesh( mesh );
        JNodeSequence nodes = jc.getSimplified();

        JColor highlightColor;
        highlightColor[0] = 1.0;
        highlightColor[1] = 0.0;
        highlightColor[2] = 0.0;
        highlightColor[3] = 0.0;

        JNodeRenderPtr nAttrib;

        size_t numnodes = mesh->getSize(0);
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            vtx->getAttribute("Render", nAttrib);
            nAttrib->glyph   =  0;
            nAttrib->display =  0;
        }

        for( const JNodePtr &vtx : nodes) {
            vtx->getAttribute("Render", nAttrib);
            nAttrib->color   =  highlightColor;
            nAttrib->glyph   =  1;
            nAttrib->display =  1;
        }
        meshViewer->updateBuffers(mesh);
    */
}

/*
void JMSTQuadMesherDialog :: remesh()
{
    if( mesh == nullptr ) return;

    if( quadmesh ) {
        quadmesh->deleteAll();
        meshViewer->removeObject(quadmesh);
    }

    int nside = 0;
    if( triRadioButton->isChecked()  ) nside = 3;
    if( quadRadioButton->isChecked() ) nside = 4;
    if( pentaRadioButton->isChecked()) nside = 5;
    if( hexaRadioButton->isChecked() ) nside = 6;
}
*/


/*
void JMSTQuadMesherDialog :: openPartitionDialog()
{
    if( partitionDialog == nullptr)
        partitionDialog.reset( new JMeshPartitionDialog(this) );
    partitionDialog->setViewManager(viewManager);
    partitionDialog->setMesh(mesh);
    partitionDialog->show();
    this->hide();
}
*/

void JMSTQuadMesherDialog :: setSeedSelect()
{
    if( mesh == nullptr ) return;


    /*
        if( seedNodeRadioButton->isChecked() )
            pick_entity = 0;
        if( seedFaceRadioButton->isChecked() )
            pick_entity = 2;
    */

}
