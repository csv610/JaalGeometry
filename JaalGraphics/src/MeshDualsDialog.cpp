#include "MeshDualsDialog.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
void JChordCellColor :: setChordID( int id )
{
    chordID = id;

    if( mapColor.find(id) != mapColor.end() ) return;

    JColor clr;

    clr[3] = alpha;
    switch(id) {
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
    mapColor[id] = clr;
}

///////////////////////////////////////////////////////////////////////////////
int JChordCellColor :: assign( const JCellPtr &cell)
{
    color = mapColor[chordID];
    cell->setAttribute("Color", color);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JSheetColor :: setSheetID( int id )
{
    chordID = id;

    if( mapColor.find(id) != mapColor.end() ) return;

    JColor clr;

    switch(id) {
    case 0:
        clr = JEntityColor::getColor("Red");
        break;
    case 1:
        clr = JEntityColor::getColor("Green");
        break;
    case 2:
        clr = JEntityColor::getColor("Blue");
        break;
    case 3:
        clr = JEntityColor::getColor("Yellow");
        break;
    case 4:
        clr = JEntityColor::getColor("Olive");
        break;
    case 5:
        clr = JEntityColor::getColor("Navy");
        break;
    case 6:
        clr = JEntityColor::getColor("Teal");
        break;
    case 7:
        clr = JEntityColor::getColor("Purple");
        break;
    case 8:
        clr = JEntityColor::getColor("Gray");
        break;
    default:
        clr = JEntityColor::getRandomColor();
        break;
    }
    clr[3] = alpha;
    mapColor[id] = clr;
}
///////////////////////////////////////////////////////////////////////////////
int JSheetColor :: assign( const JFacePtr &face)
{

    if( face == nullptr ) return 1;

    int cid = 0;
    int err = face->getAttribute("SheetID", cid);
    if( !err) {
        color = mapColor[chordID];
        face->setAttribute("Color", color);
        return 0;
    }
    return 1;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshDualsDialog :: pickSeeds()
{
    int val = -1;
    if( entityPicker == nullptr ) {
        cout << "Warning: Entity Picker not present " << endl;
        return;
    }

    if( chordDim == 2) {
        if( dualname == "Chord")  val= 1;
        if( dualname == "Column") val= 2;
    }

    if( chordDim == 3) {
        if( dualname == "Sheet")  val = 1;
        if( dualname == "Sheaf")  val = 1;
        if( dualname == "Chord")  val = 2;
        if( dualname == "Column") val = 3;
    }

    mesh->setAttribute("PickableEntity", val);
}


JMeshDualsDialog :: JMeshDualsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;

    qdual = nullptr;
    hdual = nullptr;

    nextEdgeID = 0;
    nextFaceID = 0;
    nextCellID = 0;

    dualname  = "Chord";
}

///////////////////////////////////////////////////////////////////////////////

JMeshDualsDialog :: ~JMeshDualsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: init()
{
    if( viewManager == nullptr ) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");

    if( c == nullptr) return;

    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr || meshViewer == nullptr)  return;

    int topDim = mesh->getTopology()->getDimension();

    if( topDim == 2 ) {
        int etype = mesh->getTopology()->getElementsType(2);
        if( etype != JFace::QUADRILATERAL ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText(" A 2 Dimensional mesh should have all quad elements");
            msg.setStandardButtons( QMessageBox::Ok );
            msg.exec();
            mesh = nullptr;
            return;
        }
//       meshViewer->displayAll(2,1);
    }

    if( topDim == 3 ) {
        int etype = mesh->getTopology()->getElementsType(3);
        if( etype != JCell::HEXAHEDRON ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText(" A 3 Dimensional mesh should have all hex  elements");
            msg.setStandardButtons( QMessageBox::Ok );
            msg.exec();
            mesh = nullptr;
            return;
        }
//       meshViewer->displayAll(3,1);
    }

    entityPicker = meshViewer->getEntityPicker();
    chordDim  = mesh->getTopology()->getDimension();

    qdual.reset();
    hdual.reset();
    switch( chordDim  ) {
    case 2:
        qdual.reset(new JQuadDual(mesh));
        break;
    case 3:
        //        hdual = new HexDual(mesh);
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: getStyle()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: keyPressEvent( QKeyEvent *e)
{
   if( e->key() == Qt::Key_Return ) return;
   QDialog::keyPressEvent(e);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: getAllQuadChords()
{
    if( mesh == nullptr) return;

    int nelem = mesh->getTopology()->getElementsType(2);
    if( nelem != 4  ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("All Quadrilateral mesh required ");
        msg.setStandardButtons( QMessageBox::Ok );
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    qChords = qdual->getAllChords();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualsDialog :: getQuadChords(JEdgeSequence &eseq)
{
    if( qdual == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    qChords = qdual->getChords(eseq);
}

///////////////////////////////////////////////////////////////////////////////

#ifdef CSV
void JMeshDualsDialog :: getAllHexChords()
{
    if( mesh == nullptr || hdual == nullptr || dualViewer == nullptr) return;

    int nelem = mesh->getTopology()->getElementsType(3);
    if( nelem != 8  ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("All Hexmesh required ");
        msg.setStandardButtons( QMessageBox::Ok );
        msg.exec();
        return;
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    hdual->getAllChords(dChords);
    dualViewer->setChords( dChords );
    dualViewer->refreshDisplay();
}
#endif

///////////////////////////////////////////////////////////////////////////////
void JMeshDualsDialog :: getHexChords(JFaceSequence &)
{
#ifdef CSV
    if( hdual == nullptr ) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    meshViewer->displayAll(1,0);
    meshViewer->displayAll(2,0);
    meshViewer->displayAll(3,0);

    hdual->getChords(fseq, dChords);
    dualViewer->setChords( dChords );

#endif
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: getHexSheets(JEdgeSequence &eseq)
{
#ifdef CSV
    if( eseq.empty() || hdual == nullptr ) return;

    meshViewer->displayAll(1,0);
    meshViewer->displayAll(2,0);
    meshViewer->displayAll(3,0);

    JWaitCursor waitCursor;
    waitCursor.start();

    hdual->getSheets(eseq, dSheets);
    dualViewer->setSheets( dSheets );
#endif
}

///////////////////////////////////////////////////////////////////////////////


#ifdef CSV
void JMeshDualsDialog :: getHexChords(JCellSequence &cseq)
{
    JWaitCursor waitCursor;
    waitCursor.start();

    FaceSet fset;
    for( size_t i = 0; i < cseq.size(); i++) {
        fset.insert( cseq[i]->getFaceAt(0) );
        fset.insert( cseq[i]->getFaceAt(2) );
        fset.insert( cseq[i]->getFaceAt(4) );
    }

    JFaceSequence fseq;
    FaceSet::iterator fit;
    for( fit = fset.begin(); fit != fset.end(); ++fit)
        fseq.push_back(*fit);
    getHexChords(fseq);
}
#endif
///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: getAllHexSheets()
{
#ifdef CSV
    if( mesh == nullptr || hdual == nullptr ) return;

    int nelem = mesh->getTopology()->getElementsType(3);
    if( nelem != 8  ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("All Hexmesh required ");
        msg.setStandardButtons( QMessageBox::Ok );
        msg.exec();
        return;
    }
    JWaitCursor waitCursor;
    waitCursor.start();

    hdual->getAllSheets( dSheets );

    if( dualViewer ) dualViewer->setSheets(dSheets);
#endif
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: modifyDual()
{
    if( mesh == nullptr ) return ;

#ifdef CSV

    if( dChords.size() != 1) return;

    dChords[0]->setMesh(mesh);

    if( flipChordRadioButton->isChecked() ) {
        dChords[0]->flip();
    }

    if( diceChordRadioButton->isChecked() ) {
        dChords[0]->diceQuads();
        dChords[0]->deleteSegments();
        meshViewer->displayAll(1);
    }


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
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: clearChords()
{
    if( meshViewer == nullptr) return;
    /*

    //    bool nval = meshViewer->isEnabled(0);
    //    bool eval = meshViewer->isEnabled(1);
        bool fval = 0;

        JNodeRenderPtr nAttrib;
        JEdgeRenderPtr eAttrib;
        JFaceRenderPtr fAttrib;

        JFaceSequence fseq;
        for( size_t i = 0; i < qChords.size(); i++) {
            fseq = qChords[i]->getFaces();
            for( size_t j = 0; j < fseq.size(); j++) {
                const JFacePtr &f = fseq[j];
                f->getAttribute("Render", fAttrib);
                fAttrib->display = fval;
                for( int k = 0; k < f->getSize(1); k++) {
                    const JEdgePtr &e = f->getEdgeAt(k);
                    e->getAttribute("Render", eAttrib);
                    eAttrib->display = eval;
                }

                for( int k = 0; k < f->getSize(0); k++) {
                    const JNodePtr &v = f->getNodeAt(k);
                    v->getAttribute("Render", nAttrib);
                    nAttrib->display = nval;
                }
            }
    //      qChords[i]->deleteSegments();
        }
        qChords.clear();
    */
}

/////////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: resetPrimalMesh()
{
//  if( meshViewer == nullptr || dualViewer == nullptr ) return;
//  meshViewer->displayAll(1);
    clearChords();
}

///////////////////////////////////////////////////////////////////////////////

/*
int JMeshDualsDialog :: getSheaves( const JEdgePtr &edge)
{
    JFaceSequence fseq, fneighs;

    edge->getRelations( fneighs );

    for( size_t  i  = 0; i < fneighs.size(); i++)
        if( fneighs[i]->isBoundary() ) fseq.push_back( fneighs[i] );

    if( fseq.size() == 2 ) {
        getHexChords( fseq );
        return 0;
    }
    return 1;
}
*/
///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: getSomeDual()
{
    if( mesh == nullptr)  return;

    clearChords();

    JEdgeSequence eseq;
    JFaceSequence fseq, fneighs;
    JCellSequence cseq;

    if( chordDim  == 2 ) {
        if( dualname == "Chord") {
            eseq = entityPicker->getPickedEdges();
            if( eseq.empty() ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("At least one edge must be picked for chord calculation");
                msg.setStandardButtons( QMessageBox::Ok);
                msg.exec();
                return;
            }
            qChords = qdual->getChords(eseq);

        }

        if( dualname == "Sheaf") {
            eseq = entityPicker->getPickedEdges();
            if( eseq.size() < 2 ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("At least two edges must be picked for sheaf calculation");
                msg.setStandardButtons( QMessageBox::Ok);
                msg.exec();
                return;
            }
            qChords = qdual->getChords(eseq);
        }

        if( dualname == "Column") {
            fseq = entityPicker->getPickedFaces();
            if( fseq.empty() ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("At least one face must be picked for column calculation");
                msg.setStandardButtons( QMessageBox::Ok);
                msg.exec();
                return;
            }
            qChords = qdual->getColumns(fseq);
        }
    }

#ifdef CSV
    if( chordDim == 3 ) {
        if( dualname == "Sheaf") {
            if( mesh->getAdjTable(1,2) == 0) mesh->buildRelations(1,2);
            eseq = entityPicker->getPickedEdges();
            if( eseq.size() < 2  ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("At least two adjacenct faces must be picked for sheaf calculation");
                msg.setStandardButtons( QMessageBox::Ok);
                msg.exec();
                return;
            }
            getSheaves( eseq[0] );
        }

        if( dualname == "Chord" ) {
            fseq = entityPicker->getPickedFaces();
            if( fseq.empty() ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("At least one face must be picked for chord calculation");
                msg.setStandardButtons( QMessageBox::Ok);
                msg.exec();
                return;
            }
            getHexChords( fseq );
        }

        if( dualname == "Column" ) {
            cseq = entityPicker->getPickedCells();
            if( cseq.empty() ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("At least one cell must be picked for  column calculation");
                msg.setStandardButtons( QMessageBox::Ok);
                msg.exec();
                return;
            }
            getHexChords( cseq );
        }

        if( dualname == "Sheet" ) {
            eseq = entityPicker->getPickedEdges();
            if( eseq.empty() ) {
                QMessageBox msg;
                msg.setIcon(QMessageBox::Warning);
                msg.setText("At least one edge must be picked for  Sheet calculation");
                msg.setStandardButtons( QMessageBox::Ok);
                msg.exec();
                return;
            }
            getHexSheets( eseq );
        }
    }

    if( entityPicker ) entityPicker->clearAll();

    dualViewer->refreshDisplay();
#endif
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: getAllChords()
{
    if( chordDim == 2 ) getAllQuadChords();
//    if( chordDim == 3 ) getAllHexChords();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshDualsDialog :: checkDisplay()
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
///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: nextSeeds()
{
    if( mesh == nullptr )  return;

    clearChords();

    static int ncnt = 0;

    if( dualname == "Chord" ) {
        if( chordDim == 2 ) {
            size_t numEdges = mesh->getSize(1);
            if( nextEdgeID >= numEdges ) nextEdgeID = 0;
            nextSeedLineEdit->setText( QString::number(nextEdgeID) );
            JEdgePtr edge = mesh->getEdgeAt(nextEdgeID);
            JEdgeSequence eseq(1);
            eseq[0] = edge;
            getQuadChords( eseq );
            nextEdgeID++;
        }

        if( chordDim == 3 ) {
            size_t numFaces = mesh->getSize(2);
            if( nextFaceID >= numFaces ) nextFaceID = 0;
            nextSeedLineEdit->setText( QString::number(nextFaceID) );
            JFacePtr face = mesh->getFaceAt( nextFaceID);
            JFaceSequence fseq(1);
            fseq[0] = face;
//          getHexChords( fseq );
            nextFaceID++;
        }
    }

#ifdef CSV
    if( dualname == "Column" ) {
        size_t numCells = mesh->getSize(3);
        if( nextCellID >= numCells ) nextCellID = 0;
        nextSeedLineEdit->setText( QString::number(nextCellID) );
        JFaceSequence fseq(3);
        Cell *cell = mesh->getCellAt( nextCellID);
        fseq[0] = cell->getFaceAt(0);
        fseq[1] = cell->getFaceAt(2);
        fseq[2] = cell->getFaceAt(4);
        getHexChords( fseq );
        nextCellID++;
    }

    if( dualname == "Sheaf" ) {
        if( mesh->getAdjTable(1,2) == 0) mesh->buildRelations(1,2);
        size_t numEdges = mesh->getSize(2);
        if( nextEdgeID >= numEdges ) nextEdgeID = 0;
        for( size_t i = nextEdgeID; i < numEdges; i++) {
            Edge *edge = mesh->getEdgeAt(i);
            if(getSheaves( edge ) == 0) {
                nextSeedLineEdit->setText( QString::number(i) );
                nextEdgeID = i+1;
                return;
            }
        }
        nextEdgeID = numEdges;
    }

    if( dualname == "Sheet" ) {
        if( mesh->getAdjTable(1,3) == 0) mesh->buildRelations(1,3);
        size_t numEdges = mesh->getSize(1);
        if( nextEdgeID >= numEdges ) nextEdgeID = 0;
        Edge *edge = mesh->getEdgeAt(nextEdgeID);
        JEdgeSequence eseq(1);
        eseq[0] = edge;
        getHexSheets( eseq );
        nextEdgeID++;
        nextSeedLineEdit->setText( QString::number(nextEdgeID) );
    }

    dualViewer->refreshDisplay();
    meshViewer->refreshDisplay();
#endif
}

///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: selectDual()
{
    clearChords();
    QString qs = dualComboBox->currentText();
    dualname = qs.toUtf8().constData();
    pickSeeds();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog ::  displayDuals()
{
    /*
        if( dualViewer == nullptr ) return;
        bool v = displayDualCheckBox->isChecked();
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshDualsDialog :: reject()
{
    clearChords();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshDualsDialog :: makeConnections()
{
    CheckBox( elementsCheckBox,  [=] {checkDisplay();});
    CheckBox( parEdgesCheckBox,  [=] {checkDisplay();});
    CheckBox( pickSeedsCheckBox,  [=] {pickSeeds();});
    CheckBox( displayDualCheckBox,  [=] {displayDuals();});

    PushButton( getAllChordsPushButton,  [=] {getAllChords();} );
    PushButton( clearAllPushButton, [=] {resetPrimalMesh();});
    PushButton( getPickedChordsPushButton, [=] {getSomeDual();});
}

///////////////////////////////////////////////////////////////////////////////
