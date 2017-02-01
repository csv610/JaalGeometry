#include "QuadMeshDualsDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
void JChordNodeColor :: generateColors( int n )
{
    for( int i = 0; i < n; i++) {
        if( mapColor.find(i) != mapColor.end() ) return;

        JColor clr;

        clr[3] = alpha;
        switch(i) {
        case 0:
            clr[0] = Red[0];
            clr[1] = Red[1];
            clr[2] = Red[2];
            break;
        case 1:
            clr[0] = Green[0];
            clr[1] = Green[1];
            clr[2] = Green[2];
            break;
        case 2:
            clr[0] = Blue[0];
            clr[1] = Blue[1];
            clr[2] = Blue[2];
            break;
        default:
            clr  = JEntityColor::getRandomColor();
            break;
        }
        mapColor[i] = clr;
    }
}

///////////////////////////////////////////////////////////////////////////////

int JChordNodeColor :: assign(const JNodePtr &vertex)
{
    if( vertex == nullptr ) return 1;

    JNodeRenderPtr nAttrib;
    int cid = 0;
    int err = vertex->getAttribute("ChordID", cid);
    if( !err) {
        vertex->getAttribute("Render", nAttrib);
        nAttrib->color = mapColor[cid];
        return 0;
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

void JChordEdgeColor :: generateColors(int n)
{
    for( int i = 0; i < n; i++) {
        if( mapColor.find(i) != mapColor.end() ) return;

        JColor clr;

        clr[3] = alpha;
        switch(i) {
        case 0:
            clr[0] = Red[0];
            clr[1] = Red[1];
            clr[2] = Red[2];
            break;
        case 1:
            clr[0] = Green[0];
            clr[1] = Green[1];
            clr[2] = Green[2];
            break;
        case 2:
            clr[0] = Blue[0];
            clr[1] = Blue[1];
            clr[2] = Blue[2];
            break;
        default:
            clr  = JEntityColor::getRandomColor();
            break;
        }
        mapColor[i] = clr;
    }
}

///////////////////////////////////////////////////////////////////////////////
int JChordEdgeColor :: assign( const JEdgePtr &edge)
{
    if( edge == nullptr ) return 1;
    JEdgeRenderPtr eAttrib;

    int cid = 0;
    int err = edge->getAttribute("ChordID", cid);
    if( !err) {
        edge->getAttribute("Render", eAttrib);
        if( eAttrib) eAttrib->color = mapColor[cid];
        return 0;
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

void JChordFaceColor :: generateColors( int n )
{
    for( int i = 0; i < n; i++) {
        if( mapColor.find(i) != mapColor.end() ) return;

        JColor clr;

        clr[3] = alpha;
        switch(i) {
        case 0:
            clr[0] = Red[0];
            clr[1] = Red[1];
            clr[2] = Red[2];
            break;
        case 1:
            clr[0] = Green[0];
            clr[1] = Green[1];
            clr[2] = Green[2];
            break;
        case 2:
            clr[0] = Blue[0];
            clr[1] = Blue[1];
            clr[2] = Blue[2];
            break;
        default:
            clr  = JEntityColor::getRandomColor();
            break;
        }
        mapColor[i] = clr;
    }
}


///////////////////////////////////////////////////////////////////////////////
int JChordFaceColor :: assign( const JFacePtr &face)
{
    if( face == nullptr ) return 1;

    JFaceRenderPtr fAttrib;
    int cid = 0;
    int err = face->getAttribute("ChordID", cid);
    if( !err) {
        face->getAttribute("Render", fAttrib);
        if( fAttrib) fAttrib->color = mapColor[cid];
        return 0;
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////

JQuadMeshDualsDialog :: JQuadMeshDualsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    qdual = nullptr;

    nextEdgeID = 0;
    nextFaceID = 0;

    dualname  = "Chord";
    chordFaceColor.reset( new JChordFaceColor);
    chordEdgeColor.reset( new JChordEdgeColor);
}

///////////////////////////////////////////////////////////////////////////////

JQuadMeshDualsDialog :: ~JQuadMeshDualsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c ) meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr) return;
    entityPicker = meshViewer->getEntityPicker();

    entityPicker->setMode(1);
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMeshDualsDialog :: showEvent( QShowEvent *event)
{
    if( mesh == nullptr) return;

    entityPicker->setMesh(mesh);
    enablePicking();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = nullptr;
    if( m == nullptr )  return;

    int nelem = m->getTopology()->getElementsType(2);
    if( nelem !=   JFace::QUADRILATERAL) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("All Quadrilateral mesh required ");
        msg.setStandardButtons( QMessageBox::Ok );
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
        return;
    }

    mesh = m;
    objectNameLineEdit->setText( QString(mesh->getName().c_str() ));
    qdual.reset();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: initMesh()
{
    if( mesh == nullptr) return;

    int topDim = mesh->getTopology()->getDimension();

    int chordDim  = mesh->getTopology()->getDimension();
    qdual.reset(new JQuadDual(mesh));

    clearChords();
    enablePicking();

    /*
            getMaxEdgeLength();
            double area = mesh->getGeometry()->getSurfaceArea();
            int    ntarget = area/(elen*elen);
            numTargetFacesLineEdit->setText( QString::number(ntarget) );

            size_t numfaces = mesh->getActiveSize(2);
            numCurrentFacesLineEdit->setText( QString::number(numfaces) );
        */
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: getStyle()
{
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) return;
    QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: clearChords()
{
    if( meshViewer == nullptr) return;

    JColor white = JEntityColor::getColor("White");

    JFaceRenderPtr fAttrib;
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->color   = white;
        fAttrib->display = 1;
    }

    JColor black = JEntityColor::getColor("Black");

    JEdgeRenderPtr eAttrib;
    size_t numedges = mesh->getSize(1);
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        edge->getAttribute("Render", eAttrib);
        eAttrib->color = black;
        eAttrib->scale  = 1.0;
        eAttrib->display = 1;
    }

    mesh->deleteEdgeAttribute("ChordID");
    mesh->deleteFaceAttribute("ChordID");

    qChords.clear();
    entityPicker->clearAll();
    meshViewer->updateBuffers(mesh);
}

////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: setChordID( const JQuadChordPtr &chord, int id)
{
    JFaceSequence  chordFaces = chord->getFaces();
    for( const JFacePtr &face : chordFaces)
        face->setAttribute("ChordID", id);

    JEdgeSequence  chordEdges = chord->getEdges();
    for( const JEdgePtr &edge : chordEdges)
        edge->setAttribute("ChordID", id);
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMeshDualsDialog :: displayChords()
{
    int numChords = qChords.size();
    chordFaceColor->generateColors(numChords);

    for( int i =0; i < numChords; i++)
        setChordID( qChords[i], i);

    chordFaceColor->setMesh(mesh);
// chordEdgeColor->setMesh(mesh);
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: getEdgeChords(JEdgeSequence & eseq)
{
    if( qdual == nullptr ) return;

    if( eseq.empty() ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("At least one edge must be picked for chord calculation");
        msg.setStandardButtons( QMessageBox::Ok);
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    qChords = qdual->getChords(eseq);
    displayChords();

    JColor black = JEntityColor::getColor("Black");
    JEdgeRenderPtr eAttrib;
    for( const JEdgePtr &edge: eseq) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->color = black;
        eAttrib->scale = 1.0;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: enablePicking()
{
    int val = -1;
    if( entityPicker == nullptr ) {
        cout << "Warning: Entity Picker not present " << endl;
        return;
    }

    if( dualname == "Chord")  val= 1;
    if( dualname == "Column") val= 2;

    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = val;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: modifyDual()
{
    if( mesh == nullptr || qChords.empty() ) return ;

    JDiceQuadChord dice;
    bool  diceType = completeDiceRadioButton->isChecked();
    dice.setCompleteDice(diceType);

    JFaceSequence newfaces;
    if( diceChordRadioButton->isChecked() ) {
        dice.setMesh(mesh);
        dice.setChord( qChords[0] );
        meshViewer->updateBuffers(mesh);
        newfaces = dice.getNewFaces();
        if( newfaces.empty() ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Information);
            msg.setText("Dicing over: no new faces");
            msg.setStandardButtons(QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok)  return;
        }
        JColor green = JEntityColor::getColor("Green");
        JFaceRenderPtr fAttrib;
        for( const JFacePtr &face : newfaces) {
            face->getAttribute("Render", fAttrib);
            fAttrib->color   = green;
            fAttrib->display = 1;
        }
    }

#ifdef CSV
    if( shrinkChordRadioButton->isChecked() ) {
        dChords[0]->shrink_parallel_edges();
    }
    if( removeChordRadioButton->isChecked() ) {
        JNodeSequence nodes;
        dChords[0]->get_parallel_nodes(nodes);
        bool boundary = 0;
        for( size_t i = 0; i < nodes.size(); i++) {
            if( nodes[i]->isBoundary() )  boundary = 1;
        }

        if( boundary) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("Chord is touching the boundary:");
            msg.setStandardButtons( QMessageBox::Cancel | QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Cancel)  return;
        }

        dChords[0]->delete_parallel_edges();
        meshViewer->displayAll(1);
    }

    if( dualViewer ) dualViewer->refreshDisplay();
#endif

    getMaxEdgeLength();

    entityPicker->clearAll();

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JQuadMeshDualsDialog :: mouseReleaseEvent(QMouseEvent *)
{
    if( entityPicker == nullptr || mesh == nullptr ) return;
    getNewDual();
}
///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: getNewDual()
{
    if( mesh == nullptr)  return;

    if( dualname == "Chord") {
        JEdgeSequence eseq = entityPicker->getPickedEdges();
        if( !eseq.empty() ) getEdgeChords(eseq);
    }

    if( dualname == "Column") {
        JFaceSequence fseq = entityPicker->getPickedFaces();
        if( !fseq.empty() ) qChords = qdual->getColumns(fseq);
    }
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMeshDualsDialog :: diceAllChords()
{
    if( mesh == nullptr) return;

    double maxlen = 0.0;
    JEdgeSequence seedEdges(1);
    size_t   numedges;

    JDiceQuadChord dice;
    bool  diceType = completeDiceRadioButton->isChecked();
    dice.setCompleteDice(diceType);

    for( int j = 0; j < 1; j++) {
        double maxlen = 0.0;
        numedges = mesh->getSize(1);
        dice.setMesh(mesh);
        for( size_t i = 0; i <  numedges; i++) {
            const JEdgePtr &edge = mesh->getEdgeAt(i);
            if( edge->isActive() && !edge->isBoundary() ) {
                double len = JEdgeGeometry::getLength(edge);
                if( len > maxlen ) {
                    maxlen = len;
                    seedEdges[0] = edge;
                }
            }
        }
        getEdgeChords(seedEdges);
        dice.setChord( qChords[0] );

        /*
                nonlinearOpt.setMesh(mesh);
                nonlinearOpt.setAlgorithm(JMeshNonlinearOptimization::SMART_LAPLACIAN);
                nonlinearOpt.improveQuality();
        */
        meshViewer->updateGeometryBuffers(mesh);
        meshViewer->updateBuffers(mesh);
        getMaxEdgeLength();
    }

    size_t numfaces = mesh->getActiveSize(2);
    numCurrentFacesLineEdit->setText( QString::number(numfaces) );

}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: getAllChords()
{
    if( mesh == nullptr ) return;
    if( qdual == nullptr) initMesh();

    JWaitCursor waitCursor;
    waitCursor.start();

    clearChords();

    qChords = qdual->getAllChords();
    displayChords();

    numChordsLineEdit->setText( QString::number( qChords.size() ));

}
///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: nextSeeds()
{
    if( mesh == nullptr )  return;

    clearChords();

    static int ncnt = 0;

    if( dualname == "Chord" ) {
        size_t numEdges = mesh->getSize(1);
        if( nextEdgeID >= numEdges ) nextEdgeID = 0;
        nextSeedLineEdit->setText( QString::number(nextEdgeID) );
        JEdgePtr edge = mesh->getEdgeAt(nextEdgeID);
        JEdgeSequence eseq(1);
        eseq[0] = edge;
        getEdgeChords( eseq );
        nextEdgeID++;
    }
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMeshDualsDialog :: selectDual()
{
    clearChords();
    QString qs = dualComboBox->currentText();
    dualname = qs.toUtf8().constData();
    enablePicking();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog ::  getMaxEdgeLength()
{
    if( mesh == nullptr) return;

    JMeshQuality mq;
    vector<double> elen;
    mq.setMesh(mesh);

    elen = mq.getEdgesQuality(JMeshQuality::EDGE_LENGTH, 0, 1);
    maxEdgeLengthLineEdit->setText( QString::number(elen.back()));

    double el = mesh->getGeometry()->getMeanEdgeLength();
    meanEdgeLengthLineEdit->setText( QString::number(el));

}
///////////////////////////////////////////////////////////////////////////////
void JQuadMeshDualsDialog :: displayCyclicChords()
{
    if( qdual == nullptr) initMesh();
    clearChords();

    qChords = qdual->getAllCyclicChords();
    displayChords();

    numCyclesLineEdit->setText( QString::number(qChords.size()));
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: openBoundaryLayerDialog()
{
    if( boundaryLayerDialog == nullptr)
        boundaryLayerDialog.reset( new JMeshOptBoundaryLayerDialog(this));
    boundaryLayerDialog->setViewManager( viewManager );
    boundaryLayerDialog->setMesh( mesh );
    boundaryLayerDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JQuadMeshDualsDialog :: closeDialog()
{
    if( mesh ) {
        JMeshRenderPtr mrender;
        mesh->getAttribute("Render", mrender);
        mrender->pickableEntity = -1;
    }
    clearChords();

    viewManager->detach( this );
    parentWidget()->show();
    close();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadMeshDualsDialog :: makeConnections()
{
    PushButton( pickPushButton, [=] { enablePicking();});
    PushButton( displayCyclicChordsPushButton, [=] {displayCyclicChords();});
    PushButton( getAllChordsPushButton, [=] {getAllChords();} );
    PushButton( clearAllPushButton, [=] {clearChords();});
    PushButton( getPickedChordsPushButton, [=] {getNewDual();});
    PushButton( diceAllChordsPushButton, [=] {diceAllChords();});
    PushButton( modifyDualPushButton, [=] {modifyDual();});
    PushButton( optBoundaryLayerPushButton, [=] {openBoundaryLayerDialog();});
    PushButton( maxEdgeLengthPushButton, [=] {getMaxEdgeLength();});

//  PushButton( elementsCheckBox,  [=] {checkDisplay();});
//  PushButton( parEdgesCheckBox,  [=] {checkDisplay();});

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////


void JQuadMeshDualsDialog :: checkDisplay()
{
    /*
        bool val;
        val = beadsCheckBox->isChecked();
        dualViewer->displayBeads(val);

        val = elementsCheckBox->isChecked();
        dualViewer->displayElements(val);

        val = parEdgesCheckBox->isChecked();
        dualViewer->displayParallelEdges(val);

        val = intersectingFacesCheckBox->isChecked();
        dualViewer->displayIntersectingElements(2, val);

        val = intersectingCellsCheckBox->isChecked();
        dualViewer->displayIntersectingElements(3, val);

        val = touchingNodesCheckBox->isChecked();
        dualViewer->displayTouchingElements(0, val);

        val = touchingEdgesCheckBox->isChecked();
        dualViewer->displayTouchingElements(1, val);

        val = touchingFacesCheckBox->isChecked();
        dualViewer->displayTouchingElements(2, val);
        dualViewer->refreshDisplay();
    */
}
