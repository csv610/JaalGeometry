#include "PolyCubesDialog.hpp"
#include "NormalSmoothing.hpp"


JPolyCubesDialog :: JPolyCubesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JPolyCubesDialog :: ~JPolyCubesDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c == nullptr)
        c = JMeshViewer::registerComponent(viewManager);
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    assert( meshViewer );

    entityPicker = meshViewer->getEntityPicker();
    /*
        entityPicker->setMode(2);
        Color  white;
        white[0] = 1.0;
        white[1] = 1.0;
        white[2] = 1.0;
        white[3] = 1.0;
        entityPicker->setHighlightColor(white,2);
    */
    setMesh( meshViewer->getCurrentMesh() );
}
///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: setMesh( const JMeshPtr &m)
{
    modelMesh = m;
    if( modelMesh == nullptr ) return;

    int top = modelMesh->getTopology()->getDimension();
    if( top == 3) {
        cells = modelMesh->getCells();
        modelMesh->deleteCells();
    }

    polyCubes.setModelMesh(modelMesh);

    string name = modelMesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );
    objectNameLineEdit->setText( QString(modelMesh->getName().c_str() ) );
    assignColors( modelMesh );
}

///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: checkSide()
{
    if( modelMesh == nullptr ) return;

    bool val[6];
    val[0] = side1CheckBox->isChecked();
    val[1] = side2CheckBox->isChecked();
    val[2] = side3CheckBox->isChecked();
    val[3] = side4CheckBox->isChecked();
    val[4] = side5CheckBox->isChecked();
    val[5] = side6CheckBox->isChecked();

    int side;
    JFaceRenderPtr attrib;
    size_t numfaces = modelMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        JFacePtr face = modelMesh->getFaceAt(i);
        if( face->isActive() ) {
            int err = face->getAttribute("CubeSide", side);
            if( !err) {
                face->getAttribute("Render", attrib);
                attrib->display = val[side];
            }
        }
    }
    meshViewer->updateBuffers(modelMesh);
}

///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: checkTopology()
{
    int err = polyCubes.buildCubicalTopology();

    if( err ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("Warning Cublical topology is not correct");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }
}
///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: assignColor( const JFacePtr &face)
{
    if( !face->isActive() ) return;

    int side;
    int err = face->getAttribute("CubeSide", side);
    if( err ) return;

    JFaceRenderPtr attrib;
    face->getAttribute("Render", attrib);

    JColor  color;
    color[0] = 0.3;
    color[1] = 0.3;
    color[2] = 0.3;
    color[3] = 1.0;

    countfaces[side]++;
    switch(side)
    {
    case JHexahedron::LEFT_SIDE:
        color[0] = 1.0;
        color[1] = 1.0;
        color[2] = 0.0;
        color[3] = 1.0;
        break;
    case JHexahedron::RIGHT_SIDE:
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 0.0;
        color[3] = 1.0;
        break;
    case JHexahedron::BOTTOM_SIDE:
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 1.0;
        color[3] = 1.0;
        break;
    case JHexahedron::TOP_SIDE:
        color[0] = 0.0;
        color[1] = 1.0;
        color[2] = 0.0;
        color[3] = 1.0;
        break;
    case JHexahedron::BACK_SIDE:
        color[0] = 1.0;
        color[1] = 0.0;
        color[2] = 1.0;
        color[3] = 1.0;
        break;
    case JHexahedron::FRONT_SIDE:
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 1.0;
        color[3] = 1.0;
        break;
    }
    attrib->color = color;
}

///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: assignColor( const JEdgePtr &edge)
{
    if( !edge->isActive() ) return;

    JEdgeRenderPtr attrib;
    int err = edge->getAttribute("Render", attrib);
    if( err ) return;

    JColor color;
    if( edge->hasAttribute("Interface"))
    {
        color[0] = 1.0;
        color[1] = 1.0;
        color[2] = 1.0;
        color[3] = 1.0;
        attrib->scale = 3.0;
    } else {
        color[0] = 0.3;
        color[1] = 0.3;
        color[2] = 0.3;
        color[3] = 1.0;
        attrib->scale = 1.0;
    }
    attrib->color = color;
}
///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: assignColor( const JNodePtr &vtx)
{
    if( !vtx->isActive() ) return;

    JNodeRenderPtr attrib;
    int err = vtx->getAttribute("Render", attrib);
    if( err ) return;

    JColor color;
    if( vtx->hasAttribute("Interface"))
    {
        color[0] = 1.0;
        color[1] = 1.0;
        color[2] = 1.0;
        color[3] = 1.0;
        attrib->scale = 1.0;
        attrib->glyph = 1;
    } else {
        color[0] = 0.3;
        color[1] = 0.3;
        color[2] = 0.3;
        color[3] = 1.0;
        attrib->scale = 1.0;
        attrib->glyph = 0;
    }

    if( vtx->hasAttribute("PartitionCorner"))
    {
        color[0] = 1.0;
        color[1] = 1.0;
        color[2] = 1.0;
        color[3] = 1.0;
        attrib->scale = 1.0;
        attrib->glyph = 1;
    }
    attrib->color = color;
}

///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: straighten( const JEdgeSequence &eseq)
{

}
///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: flatten( const JMeshPtr &patch)
{

}

///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: assignColors( const JMeshPtr &amesh)
{
    if( amesh == nullptr) return;

    size_t numfaces = amesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++)
        assignColor( amesh->getFaceAt(i) );

    size_t numedges = amesh->getSize(1);
    for( size_t i = 0; i < numedges; i++)
        assignColor( amesh->getEdgeAt(i) );

    size_t numnodes = amesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++)
        assignColor( amesh->getNodeAt(i) );

    meshViewer->updateBuffers(amesh);
}

/////////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: initialSegments()
{
    string name = StdString(objectNameLineEdit->text());
    setMesh(meshViewer->getMesh(name));
    if( modelMesh == nullptr ) return ;

    JWaitCursor waitCursor;
    waitCursor.start();

    modelMesh->getGeometry()->setFacesNormal();

    polyCubes.initialSegmentation();

    /*
           int numParts = mp.getNumPartitions();
                displayPatchSpinBox->setMinimum(0);
                displayPatchSpinBox->setMaximum(numParts-1);

                int numInterfaces = mp.getNumInterfaces();
                displayInterfaceSpinBox->setMinimum(0);
                displayInterfaceSpinBox->setMaximum(numInterfaces-1);
        */

    assignColors( modelMesh );
}
///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: regionGrow()
{
    string name = StdString(objectNameLineEdit->text());
    modelMesh = meshViewer->getMesh(name);

    if( modelMesh == nullptr ) return ;

    JWaitCursor waitCursor;
    waitCursor.start();

    polyCubes.regionGrow();

    assignColors( modelMesh );
}
///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: optSurfmesh()
{
    JMeshNormalSmoothing smooth;
    string name = StdString(objectNameLineEdit->text());
    modelMesh = meshViewer->getMesh(name);
    smooth.setMesh(modelMesh);
    smooth.execute();

//  polyCubes.setMesh(mesh);

    assignColors( modelMesh );
    meshViewer->updateBuffers(modelMesh);

    QApplication::restoreOverrideCursor();

}
///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: openTetmesherDialog()
{
    if( tetmesherDialog == nullptr) tetmesherDialog.reset(new JTetMesherDialog(this));
    tetmesherDialog->setViewManager(viewManager);
    tetmesherDialog->setMesh(modelMesh);
    tetmesherDialog->show();
}
///////////////////////////////////////////////////////////////////////////////


void JPolyCubesDialog :: optVolmesh()

{
    /*
        if( laplaceDialog == nullptr) laplaceDialog.reset(new JMeshLaplaceSmoothingDialog(this));
        laplaceDialog->setViewManager(viewManager);
        laplaceDialog->setMesh(mesh);
        laplaceDialog->show();
    */
}
///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: openPaintingDialog()

{
    if( meshpaintingDialog == nullptr) meshpaintingDialog.reset(new JMeshPaintingDialog(this));
    meshpaintingDialog->setViewManager(viewManager);
    meshpaintingDialog->setMesh(modelMesh);
    meshpaintingDialog->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: realign()
{
    /*
        if( entityPicker == nullptr || mesh == nullptr) return;
        JFaceSequence faces = entityPicker->getPickedFaces();
        string str = StdString( alignComboBox->currentText());

        int side = -1;
        if( str == "+X")  side = Hexahedron::RIGHT_SIDE;
        if( str == "-X")  side = Hexahedron::LEFT_SIDE;
        if( str == "+Y")  side = Hexahedron::TOP_SIDE;
        if( str == "-Y")  side = Hexahedron::BOTTOM_SIDE;
        if( str == "+Z")  side = Hexahedron::FRONT_SIDE;
        if( str == "-Z")  side = Hexahedron::BACK_SIDE;
        assert( side >= 0);

        size_t numfaces = faces.size();
        for( size_t i = 0; i < numfaces; i++) {
            faces[i]->setAttribute("Partition", side);
            assignColor( faces[i] );
        }
        entityPicker->clearAll();
    */

    meshViewer->updateBuffers(modelMesh);
}
///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: saveAs()
{
    if( modelMesh == nullptr ) return;

    static QString lastSelectedDirectory;

    QString qstr  = QFileDialog::getSaveFileName(this,
                    *new QString("Select Mesh File "),
                    lastSelectedDirectory,
                    *new QString( "Mesh Format (*.xml)"));

    size_t numfaces = modelMesh->getSize(2);

    /*
        JFaceRenderPtr attrib;
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Render", attrib);
                face->setAttribute("Color",  attrib->color);
            }
        }
    */
    JFace::registerAttribute("CubeSide",  "int");
    JFace::registerAttribute("Partition", "int");

    JWaitCursor waitCursor;
    waitCursor.start();

    string meshFileName = StdString(qstr);
    if (!meshFileName.empty()) {
        JMeshXMLExporter mexp;
        mexp.addFaceAttribute("CubeSide");
        mexp.addFaceAttribute("Partition");
        mexp.writeFile(modelMesh, meshFileName);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: smoothInterface()
{
    /*
        if( mesh == nullptr) return;

        int pid, interfaceID = displayInterfaceSpinBox->value();

        JEdgeSequence edges;
        size_t numedges = mesh->getSize(1);
        JEdgeRenderPtr eAttrib;
        for( size_t i  = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                edge->getAttribute("Render", eAttrib);
                int err = edge->getAttribute("Interface", pid );
                if( !err &&  pid == interfaceID ) edges.push_back(edge);
            }
        }
        JEdgeGeometry::smooth( edges);
        meshViewer->updateGeometryBuffers(mesh);
    */
}
///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: alignAlongXYZ()
{
    if( modelMesh == nullptr) return;

    polyCubes.alignAlongXYZPlanes();

    modelMesh->enumerate(0);
    modelMesh->enumerate(1);
    modelMesh->enumerate(2);

    meshViewer->updateBuffers(modelMesh);
}
///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: generate()
{

}
///////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: showPatch()
{
    if( modelMesh == nullptr) return;
    int id = showPatchSpinBox->value();

    JFaceRenderPtr attrib;

    size_t numfaces = modelMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face  = modelMesh->getFaceAt(i);
        if( face->isActive() ) {
            int err = face->getAttribute("Render", attrib);
            attrib->display = 0;
        }
    }

    JMeshPartitioner mp;
    mp.setMesh(modelMesh);

    JMeshPtr submesh = mp.getSubMesh(id);
    numfaces = submesh->getSize(2);
    if( numfaces  == 0) return;

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face  = submesh->getFaceAt(i);
        if( face->isActive() ) {
            int err = face->getAttribute("Render", attrib);
            attrib->display = 1;
        }
    }

    meshViewer->updateBuffers(modelMesh);
}

///////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: setCubeSide()
{
    if( modelMesh == nullptr) return;
    int id = showPatchSpinBox->value();

    JFaceRenderPtr attrib;

    size_t numfaces = modelMesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face  = modelMesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", attrib);
            attrib->display = 0;
        }
    }

    QString qstr = cubeSideComboBox->currentText();
    string str = qstr.toUtf8().constData();

    int side;
    if( str == "Right")  side  = JHexahedron::RIGHT_SIDE;
    if( str == "Left")   side  = JHexahedron::LEFT_SIDE;
    if( str == "Top")    side  = JHexahedron::TOP_SIDE;
    if( str == "Buttom") side  = JHexahedron::BOTTOM_SIDE;
    if( str == "Front")  side  = JHexahedron::FRONT_SIDE;
    if( str == "Back")   side  = JHexahedron::BACK_SIDE;

    JMeshPartitioner mp;
    mp.setMesh(modelMesh);

    JMeshPtr submesh = mp.getSubMesh(id);
    numfaces = submesh->getSize(2);
    if( numfaces  == 0) return;

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face  = submesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", attrib);
            face->setAttribute("CubeSide", side);
            attrib->display = 1;
        }
    }
    assignColors( modelMesh );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: loadPolyCubes()
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

    JMeshPtr m = JMeshIO::readFile(fileName);

    if( m ) {
        if( polyCubes.setPolycubes(m) ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Warning Polycube mesh has different topology: polycubes discarded");
            msg.setStandardButtons( QMessageBox::Ok);
            msg.exec();
        } else {
            polyMesh1 = m;
            meshViewer->addObject(polyMesh1);
            assignColors(m);
            assignColors(modelMesh);
            displayPolycubes1CheckBox->setChecked(true);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: showMesh()
{
    bool val;

    val = displayModelCheckBox->isChecked();
    if( modelMesh ) {
        modelMesh->setActiveBit(val);
        meshViewer->updateBuffers(modelMesh);
    }

    val = displayPolycubes1CheckBox->isChecked();
    if( polyMesh1) {
        polyMesh1->setActiveBit(val);
        meshViewer->updateBuffers(polyMesh1);
    }

    val = displayPolycubes2CheckBox->isChecked();
    if( polyMesh2) {
        polyMesh2->setActiveBit(val);
        meshViewer->updateBuffers(polyMesh2);
    }

    val = displayCubicalMeshCheckBox->isChecked();
    if( cubicalMesh ) {
        cubicalMesh->setActiveBit(val);
        meshViewer->updateBuffers(cubicalMesh);
    }

    val = displayHexMeshCheckBox->isChecked();
    if( hexMesh ) {
        hexMesh->setActiveBit(val);
        meshViewer->updateBuffers(hexMesh);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void JPolyCubesDialog :: integerSnap()
{
    /*
        polyMesh2 = polyCubes.integerSnap();
        if( polyMesh2) {
        meshViewer->updateBuffers(polyMesh2);
        }
    */
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void JPolyCubesDialog :: makeConnections()
{
    PushButton( alignPushButton, [=] {alignAlongXYZ(); });
    PushButton( regionGrowPushButton,  [=] {regionGrow(); });
    PushButton( initialSegmentsPushButton, [=] {initialSegments(); });
    PushButton( checkTopologyPushButton, [=] {checkTopology(); });

    CheckBox( side1CheckBox, [=] {checkSide();});
    CheckBox( side2CheckBox, [=] {checkSide();});
    CheckBox( side3CheckBox, [=] {checkSide();});
    CheckBox( side4CheckBox, [=] {checkSide();});
    CheckBox( side5CheckBox, [=] {checkSide();});
    CheckBox( side6CheckBox, [=] {checkSide();});

    CheckBox( displayModelCheckBox,  [=] {showMesh();});
    CheckBox( displayHexMeshCheckBox, [=] {showMesh();});
    CheckBox( displayPolycubes1CheckBox, [=] {showMesh();});
    CheckBox( displayPolycubes2CheckBox, [=] {showMesh();});
    CheckBox( displayCubicalMeshCheckBox, [=] {showMesh();});

    PushButton( loadPolyCubesPushButton, [=] {loadPolyCubes();});
    PushButton( integerSnapPushButton, [=] {integerSnap();});
    PushButton( showPatchPushButton,  [=] {showPatch();});
    PushButton( assignSidePushButton, [=] {setCubeSide();});
    PushButton( savePushButton,  [=] {saveAs();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////

/*
void JPolyCubesDialog :: displayInterface()
{
    size_t numedges = mesh->getSize(1);
    JEdgeRenderPtr eAttrib;

    if( allInterfaceRadioButton->isChecked() ) {
        for( size_t i  = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                assignColor(edge);
                edge->getAttribute("Render", eAttrib);
                eAttrib->display = 1;
            }
        }
        meshViewer->updateBuffers(mesh);
    }

    if( noneInterfaceRadioButton->isChecked() ) {
        meshViewer->displayAll(mesh, 1, 0);
    }

    if( selectiveInterfaceRadioButton->isChecked() ) {
        size_t numnodes = mesh->getSize(0);
        JNodeRenderPtr nAttrib;
        for( size_t i  = 0; i < numnodes; i++) {
            const JNodePtr &vertex = mesh->getNodeAt(i);
            if( vertex->isActive() ) {
                vertex->getAttribute("Render", nAttrib);
                nAttrib->display = 1;
                nAttrib->glyph   = JNodeDraw::NODE_AS_POINT;
            }
        }

        int pid;
        int interfaceID = displayInterfaceSpinBox->value();

        for( size_t i  = 0; i < numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() ) {
                edge->getAttribute("Render", eAttrib);
                eAttrib->display = 0;
                int err = edge->getAttribute("Interface", pid );
                if( !err )
                    if( pid == interfaceID ) eAttrib->display = 1;
            }
        }
        meshViewer->updateBuffers(mesh);
    }
    meshViewer->refreshDisplay();
}
*/
/////////////////////////////////////////////////////////////////////////////////

/*
void JPolyCubesDialog :: displayPatch()
{
    if( allPatchRadioButton->isChecked() ) {
        meshViewer->displayAll(mesh, 2, 1);
    }

    if( nonePatchRadioButton->isChecked() ) {
        meshViewer->displayAll(mesh, 2, 0);
    }
    if( selectivePatchRadioButton->isChecked() ) {

        int partID = displayPatchSpinBox->value();
        size_t numnodes = mesh->getSize(0);
        JNodeRenderPtr nAttrib;
        for( size_t i  = 0; i < numnodes; i++) {
            const JNodePtr &vertex = mesh->getNodeAt(i);
            if( vertex->isActive() ) {
                vertex->getAttribute("Render", nAttrib);
                nAttrib->display = 1;
                nAttrib->glyph   = JNodeDraw::NODE_AS_POINT;
            }
        }
        size_t numfaces = mesh->getSize(2);
        int pid = 0;
        JFaceRenderPtr fAttrib;
        for( size_t i  = 0; i < numfaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            if( face->isActive() ) {
                face->getAttribute("Partition", pid );
                face->getAttribute("Render", fAttrib);
                if( pid == partID )
                    fAttrib->display = 1;
                else
                    fAttrib->display = 0;
            }
        }
        meshViewer->updateBuffers(mesh);
    }
meshViewer->refreshDisplay();
}

void JPolyCubesDialog :: setMesh( const string &name)
{
    JMeshPtr orgmesh = meshViewer->getMesh(name);
    if( orgmesh == nullptr) return;
    meshViewer->displayAll(orgmesh, 0);
    mesh = orgmesh->deepCopy();
    string name2 = orgmesh->getName() + "_polymesh";
    objectNameLineEdit->setText( QString(name2.c_str() ) );
    meshViewer->addMesh(mesh);
}
*/

