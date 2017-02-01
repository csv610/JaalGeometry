#include "MeshCellsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshCellsDialog :: JMeshCellsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    viewManager   = nullptr;
    meshViewer    = nullptr;

    numCellsLineEdit->setText( QString::number(0) );
//  explodeValLineEdit->setText( QString::number(0.8) );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshCellsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    /*
         cellIDSlider->setMinimum(0);
         cellIDSlider->setMaximum(numCells-1);
         bool val = meshViewer->isEnabled(3);
         displayCellsCheckBox->setChecked( 1 );
         meshViewer->displayAll(2,val);
    //   lowerFacesCheckBox->setChecked( 1 );
    */

}
///////////////////////////////////////////////////////////////////////////////

void JMeshCellsDialog :: setMesh(const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    size_t numCells = mesh->getActiveSize(3);
    numCellsLineEdit->setText( QString::number(numCells) );
    setNumVisible();


    size_t numtets = mesh->getTopology()->countElementType( JCell::TETRAHEDRON);
    numTetsLineEdit->setText( QString::number(numtets) );

    size_t numHexs = mesh->getTopology()->countElementType( JCell::HEXAHEDRON);
    numHexsLineEdit->setText( QString::number(numHexs) );

    size_t numPolys = mesh->getTopology()->countElementType( JCell::POLYHEDRON);
    numPolyLineEdit->setText( QString::number(numPolys) );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: deleteAll()
{

}

///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: reverseAll()
{
    if( mesh == nullptr) return;
    mesh->getTopology()->reverseAll();
    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: checkDisplay()
{
    if( mesh == nullptr) return;
    bool val = displayCellsCheckBox->isChecked();
    JCellRenderPtr attrib;

    size_t nCount = 0;
    size_t numCells = mesh->getSize(3);
    for( size_t i = 0; i < numCells; i++) {
        const JCellPtr &cell = mesh->getCellAt(i);
        if( cell->isActive() ) {
            cell->getAttribute("Render", attrib);
            attrib->display = val;
            if( val) nCount++;
        }
    }
    numVisibleLineEdit->setText( QString::number(nCount) );
    if( meshViewer) meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshCellsDialog :: openAttribListDialog()
{
    if( attribListDialog == nullptr )
        attribListDialog.reset(new JMeshEntityAttribListDialog( this ));

    attribListDialog->setMesh( mesh, 3);
    attribListDialog->show();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: setNumVisible()
{
    if( meshViewer == nullptr ) return;
    size_t nCount = meshViewer->getNumVisible(3);
    numVisibleLineEdit->setText( QString::number(nCount) );

    bool val = 1;
    if( nCount)
        displayCellsCheckBox->setChecked(val);
    else
        displayCellsCheckBox->setChecked(!val);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: saveSurfmesh()
{
    static QString lastSelectedDirectory;

    if( mesh == nullptr) return;

    JMeshPtr surfmesh = mesh->getTopology()->getSurfaceMesh();
    if( surfmesh == nullptr) return;

    QString qstr  = QFileDialog::getSaveFileName(this,
                    *new QString("Select Mesh File "),
                    lastSelectedDirectory,
                    *new QString( "Mesh Format (*.off *.obj *.xml *.vtk *.ele)"));

    string meshFileName = StdString(qstr);

    if (!meshFileName.empty()) {
        JMeshIO::saveAs(surfmesh, meshFileName);
    }

}
///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: closeDialog()
{
    /*
         if( exploded_mesh) {
              exploded_mesh->deleteAll();
              delete exploded_mesh;
              exploded_mesh = nullptr;
    //        meshViewer->setNewMesh( mesh );
         }
    */

    this->parentWidget()->show();
    this->close();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: openMeshSlicerDialog()
{
    if( meshSlicerDialog == nullptr )
        meshSlicerDialog.reset(new JMeshSlicerDialog(this));
    meshSlicerDialog->setViewManager( viewManager );
    meshSlicerDialog->show();
    this->hide();
}

//////////////////////////////////////////////////////////////////////////////

void JMeshCellsDialog :: openWavefrontDialog()
{
    if( wavefrontDialog == nullptr )
        wavefrontDialog.reset(new JMeshWavefrontsDialog());
    wavefrontDialog->setViewManager( viewManager );
    wavefrontDialog->show();
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshCellsDialog :: makeConnections()
{
    PushButton( saveSurfmeshPushButton, [=] {saveSurfmesh();});
    PushButton( meshSlicerPushButton, [=] {openMeshSlicerDialog(); });
    PushButton( reverseAllPushButton, [=] {reverseAll();});
    PushButton( attribListPushButton, [=] {openAttribListDialog();});

    CheckBox( enumCheckBox, [=] {checkState();});
    CheckBox( displayCellsCheckBox, [=] {checkDisplay(); });

    PushButton( closePushButton, [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshCellsDialog :: activate( Cell *cell)
{
     if( cell == nullptr ) return;
     if( !cell->isActive() )  return;

     bool val = 1;
     cell->setAttribute("Display", val);
     for( int j = 0; j < cell->getSize(2); j++) {
          Face *face = cell->getFaceAt(j);
          face->setAttribute("Display", val);
     }

     if( changeCenterCheckBox->isChecked() ) {
          const Point3D xyz = cell->getNodeAt(0)->getXYZCoords();
          AffineTransform af(mesh);
          af.translate(-xyz[0], -xyz[1], -xyz[2]);
     }

     if( realignCheckBox->isChecked() ) {
          meshViewer->alignAlong( cell->getEdgeAt(0), 0, 0);
     }
}
*/


/*
void JMeshCellsDialog :: explode()
{
     QString qstr = explodeValLineEdit->text();
     if( exploded_mesh ) exploded_mesh->deleteAll();

     exploded_mesh =  mesh->getGeometry()->explode(qstr.toDouble() );
//   meshViewer->setNewMesh(exploded_mesh);
}
*/
///////////////////////////////////////////////////////////////////////////////

void JMeshCellsDialog :: checkState()
{
    /*
         JExit();
         if( meshViewer == nullptr ) return;
         bool val = 1;

         size_t numCells = mesh->getSize(3);
         for( size_t i = 0; i < numCells; i++)  {
              Cell *cell = mesh->getCellAt(i);
              cell->setAttribute("Display", val);
         }

         if( val ) {
              val = lowerEdgesCheckBox->isChecked();
              drawCell->display_lower_entity(1, val);
         }
    */
}

///////////////////////////////////////////////////////////////////////////////
/*
void JMeshCellsDialog :: showOne()
{
     if( mesh == nullptr ) return;

     Cell *cell;

     if( showOneCheckBox->isChecked() ) {
          meshViewer->displayAll(3,0);
          meshViewer->displayAll(2,0);

          size_t currCellID = cellIDSlider->value();
          if( currCellID < mesh->getSize(3) ) {
               cell =  mesh->getCellAt( currCellID ) ;
               activate( cell );
          }
     } else {
          meshViewer->displayAll(3,1);
          meshViewer->displayAll(2,1);
     }

     setNumVisible();
     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////
/*
void JMeshCellsDialog :: getCellID()
{
     if( mesh == nullptr ) return;

     Cell *cell = mesh->getCellAt( currCellID ) ;
     deactivate(cell);

     currCellID  = cellIDSlider->value();

     cell = mesh->getCellAt( currCellID ) ;
     activate(cell);

     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////

/*
void JMeshCellsDialog :: nextCell()
{
     if( mesh == nullptr ) return;

     if( !showOneCheckBox->isChecked() ) return;

     Cell *cell = mesh->getCellAt( currCellID );
     deactivate(cell);

     size_t numCells = mesh->getSize(3);

     for( size_t i = currCellID+1; i < numCells; i++) {
          cell = mesh->getCellAt(i);
          if( cell->isActive() ) {
               currCellID = i;
               activate(cell);
               break;
          }
     }

     cellIDSlider->setValue(currCellID);

     meshViewer->refreshDisplay();
}
*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

