#include "MeshWavefrontsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshWavefrontsDialog :: JMeshWavefrontsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshWavefrontsDialog :: ~JMeshWavefrontsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshWavefrontsDialog :: keyPressEvent( QKeyEvent *e)
{
    if( e->key() == Qt::Key_Return ) {
        if( meshViewer ) meshViewer->refreshDisplay();
        return;
    }
    QDialog::keyPressEvent(e);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshWavefrontsDialog :: setWave( int wid )
{
    bool val = 0;
    if( boundaryPeelRadioButton->isChecked() ) {
        if( wid > 1) {
            val = 0;
            for( size_t i = 0; i < cellwaves[wid-1].size(); i++) {
                JCellPtr cell = cellwaves[wid-1][i];
                for( int j = 0; j < 6; j++) {
                    JFacePtr face = cell->getFaceAt(j);
                    face->setAttribute("Display", val);
                }
                for( int j = 0; j < 12; j++) {
                    JEdgePtr edge = cell->getEdgeAt(j);
                    edge->setAttribute("Display", val);
                }
                cell->setAttribute("Display", val);
            }

        }
    }

    val = 1;
    for( size_t i = 0; i < cellwaves[wid].size(); i++) {
        JCellPtr cell = cellwaves[wid][i];
        cell->setAttribute("Display", val);
        for( int j = 0; j < 6; j++) {
            JFacePtr face = cell->getFaceAt(j);
            face->setAttribute("Display", val);
        }
        for( int j = 0; j < 12; j++) {
            JEdgePtr edge = cell->getEdgeAt(j);
            edge->setAttribute("Display", val);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JMeshWavefrontsDialog :: init()
{

    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

//   mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;

    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor));

    int layerid;

    if( bfsRadioButton->isChecked() )
        mesh->getTopology()->set_cells_wavefront(0);
    else
        mesh->getTopology()->set_cells_wavefront(1);

    cellwaves.clear();
    bool val = 0;

    size_t numCells = mesh->getSize(3);
    for( int i = 0; i < numCells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        cell->setAttribute("Display", val);
        for( int j = 0; j < 6; j++) {
            JFacePtr face = cell->getFaceAt(j);
            face->setAttribute("Display", val);
        }
        for( int j = 0; j < 12; j++) {
            JEdgePtr edge = cell->getEdgeAt(j);
            edge->setAttribute("Display", val);
        }
    }

    for( size_t i = 0; i < numCells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            cell->getAttribute("Layer", layerid);
            cellwaves[layerid].push_back(cell);
            cell->setAttribute("Display", val);
        }
    }

    currWaveID = 1;
    maxWaveID  = cellwaves.size();
    maxLevelsLineEdit->setText( QString::number(maxWaveID) );

    if( insideOutRadioButton->isChecked() ) {
        currWaveID = maxWaveID;
        waveSteps  = -1;
    } else {
        currWaveID = 1;
        waveSteps  = 1;
    }
    waveIDLineEdit->setText( QString::number(currWaveID) );

    val = 1;
    for( int i = 0; i < cellwaves[currWaveID].size(); i++) {
        JCellPtr cell = cellwaves[currWaveID][i];
        cell->setAttribute("Display", val);
        for( int j = 0; j < 6; j++) {
            JFacePtr face = cell->getFaceAt(j);
            face->setAttribute("Display", val);
        }
        for( int j = 0; j < 12; j++) {
            JEdgePtr edge = cell->getEdgeAt(j);
            edge->setAttribute("Display", val);
        }
    }

    QApplication::restoreOverrideCursor();
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshWavefrontsDialog :: setNewWave()
{
    QString qstr = waveIDLineEdit->text();
    currWaveID  = qstr.toInt();
    setWave( currWaveID );
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshWavefrontsDialog :: nextWave()
{
    if( currWaveID > maxWaveID || currWaveID < 1) return;


    currWaveID = currWaveID + waveSteps;
    waveIDLineEdit->setText( QString::number(currWaveID) );
    setWave( currWaveID );
    meshViewer->refreshDisplay();
}
///////////////////////////////////////////////////////////////////////////////
void JMeshWavefrontsDialog :: waveAnimation()
{
    if( animateCheckBox->isChecked() ) {
        bool val = 0;
        size_t numCells = mesh->getSize(3);
        for( size_t i = 0; i < numCells; i++) {
            JCellPtr cell = mesh->getCellAt(i);
            if( cell->isActive() )
                cell->setAttribute("Display", val);
        }
        for( int i = 1; i <= maxWaveID; i++) {
            currWaveID = i;
            nextWave();
        }
    }

}
///////////////////////////////////////////////////////////////////////////////
void JMeshWavefrontsDialog :: closeDialog()
{
    bool val = 1;
    size_t numCells = mesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        JCellPtr cell = mesh->getCellAt(i);
        if( cell->isActive() )
            cell->setAttribute("Display", val);
    }
    this->close();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshWavefrontsDialog :: makeConnections()
{
    PushButton( initPushButton,     [=] {init();});
    PushButton( nextWavePushButton, [=] {nextWave();});

    CheckBox( animateCheckBox,  [=] {waveAnimation();});
    LineEdit( waveIDLineEdit,   [=] {setNewWave();});

    PushButton( closePushButton,    [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
