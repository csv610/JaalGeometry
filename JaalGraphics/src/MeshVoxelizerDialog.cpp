#include "MeshVoxelizerDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshVoxelizerDialog :: JMeshVoxelizerDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager = nullptr;
    meshViewer  = nullptr;

    voxDimLineEdit->setText( QString::number(100) );
}

///////////////////////////////////////////////////////////////////////////////

JMeshVoxelizerDialog :: ~JMeshVoxelizerDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizerDialog :: showEvent( QShowEvent *e)
{
    if( meshViewer == nullptr ) return;

    mesh = meshViewer->getCurrentMesh();
    if( mesh == nullptr) return;
    string name = mesh->getName();
    objectNameLineEdit->setText( QString(name.c_str() ) );

    QDialog::showEvent(e);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizerDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
}

///////////////////////////////////////////////////////////////////////////////
void JMeshVoxelizerDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;

    string name = mesh->getName();
    if( !name.empty() )
        objectNameLineEdit->setText(QString(name.c_str()));
}

///////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizerDialog :: genMesh()
{
    JWaitCursor waitCursor;
    waitCursor.start();

    if( hexmesh)
        meshViewer->removeObject(hexmesh);

    JMeshVoxelizer voxeler;

    QString qstr;
    qstr = voxDimLineEdit->text() ;
    int nsize  = qstr.toInt();

    if( pointCloudCheckBox->isChecked() ) {
        JNodeSequence nodes = mesh->getNodes();
        voxmesh = voxeler.genVoxels(nodes, nsize);
    } else
        voxmesh = voxeler.genVoxels( mesh, nsize);

    if( voxmesh) {
        if( volumeVoxelsCheckBox->isChecked() )
            voxmesh->makeSolid();

        hexmesh = voxmesh->getModelMesh();
        assert( hexmesh );
        numVoxelsLineEdit->setText( QString::number(hexmesh->getSize(3)) );
        meshViewer->addObject(hexmesh);

        cellDim = voxmesh->getCellDimensions();
        xstartSpinBox->setRange( 0, cellDim[0] );
        ystartSpinBox->setRange( 0, cellDim[1] );
        zstartSpinBox->setRange( 0, cellDim[2] );

        xendSpinBox->setRange( 0, cellDim[0] );
        yendSpinBox->setRange( 0, cellDim[1] );
        zendSpinBox->setRange( 0, cellDim[2] );
        int nval = voxmesh->getNumSingularPoints();
        numSingularLineEdit->setText( QString::number(nval) );

    }
    QApplication::restoreOverrideCursor();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizerDialog :: openMeshlistDialog()
{
    if( meshlistDialog.get() == nullptr )
        meshlistDialog.reset(new JObjectsListDialog(this));

    /*
        meshlistDialog->setViewManager( viewManager );
        meshlistDialog->setType(1);
        meshlistDialog->show();
        this->hide();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizerDialog :: viewCell( const JCellPtr &cell, bool v)
{
    if( cell == nullptr) return;

    JFaceRenderPtr attrib;
    for(int i = 0; i < 6; i++) {
        const JFacePtr &face = cell->getFaceAt(i);
        int err = face->getAttribute("Render", attrib);
        if( !err ) attrib->display = v;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizerDialog :: slicer()
{
    if( voxmesh == nullptr) return;

    JMeshPtr bgmesh = voxmesh->getBackgroundMesh();

    size_t numCells = bgmesh->getSize(3);
    #pragma omp for
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = bgmesh->getCellAt(i);
        viewCell(cell,0);
    }

    int i0 = 0;
    int i1 = cellDim[0];

    int j0 = 0;
    int j1 = cellDim[1];

    int k0 = 0;
    int k1 = cellDim[2];

    if( xplanesCheckBox->isChecked() ) {
        i0 = xstartSpinBox->value();
        i1 = xendSpinBox->value();
    }

    if( yplanesCheckBox->isChecked() ) {
        j0 = ystartSpinBox->value();
        j1 = yendSpinBox->value();
    }

    if( zplanesCheckBox->isChecked() ) {
        k0 = zstartSpinBox->value();
        k1 = zendSpinBox->value();
    }

    #pragma omp for
    for( int k = k0; k < k1; k++) {
        for( int j = j0; j < j1; j++) {
            for( int i = i0; i < i1; i++) {
                size_t cid = k*cellDim[0]*cellDim[1] + j*cellDim[0] +i;
                if( cid < numCells) {
                    const JCellPtr &cell = bgmesh->getCellAt(cid);
                    if( cell->isActive() ) viewCell(cell,1);
                }
            }
        }
    }
    meshViewer->updateBuffers(hexmesh);
}

/////////////////////////////////////////////////////////////////////////////////////
void JMeshVoxelizerDialog :: saveAs()
{
    static QString lastSelectedDirectory;

    QString qstr  = QFileDialog::getSaveFileName(this,
                    *new QString("Select Mesh File "),
                    lastSelectedDirectory,
                    *new QString( "Mesh Format (*.xml)"));

    string fileName = StdString(qstr);
    if (!fileName.empty()) {
        if(voxmesh) voxmesh->saveAs(fileName);
    }
}
/////////////////////////////////////////////////////////////////////////////////////
void JMeshVoxelizerDialog :: readFile()
{
    static QString lastSelectedDirectory;
    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Mesh File "),
                   lastSelectedDirectory,
                   *new QString( "Mesh Format (*.xml)"));

    string fileName = qstr.toUtf8().constData();
    if (fileName.empty()) return;

    voxmesh = boost::shared_ptr<JVoxelMesh>( new JVoxelMesh);

    voxmesh->readFrom( fileName);
    hexmesh = voxmesh->getModelMesh();
    meshViewer->addObject(hexmesh);
}

/////////////////////////////////////////////////////////////////////////////////////
void JMeshVoxelizerDialog :: monotoneVoxels()
{
    if( voxmesh == nullptr) return;

    JWaitCursor waitCursor;
    waitCursor.start();

    JMonotoneVoxelizer mvoxelizer;
    mvoxelizer.setVoxelMesh( voxmesh );
    mvoxelizer.optimize();
    hexmesh = voxmesh->getModelMesh();

    int nval = voxmesh->getNumSingularPoints();
    numSingularLineEdit->setText( QString::number(nval) );
    hexmesh->getGeometry()->setFacesNormal();

    meshViewer->updateBuffers(hexmesh);
}
/////////////////////////////////////////////////////////////////////////////////////

void JMeshVoxelizerDialog :: makeConnections()
{
    SpinBoxi( xstartSpinBox,  [=] { slicer();});
    SpinBoxi( xendSpinBox,    [=] { slicer();});
    SpinBoxi( ystartSpinBox,  [=] { slicer();});
    SpinBoxi( yendSpinBox,    [=] { slicer();});
    SpinBoxi( zstartSpinBox,  [=] { slicer();});
    SpinBoxi( zendSpinBox,    [=] { slicer();});

    PushButton( saveAsPushButton,   [=] {saveAs();});
    PushButton( readFilePushButton, [=] {readFile();});

    PushButton( monotoneVoxelsPushButton, [=] {monotoneVoxels();});
    PushButton( meshlistPushButton, [=] {openMeshlistDialog();});
    PushButton( applyPushButton, [=] {genMesh();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////

