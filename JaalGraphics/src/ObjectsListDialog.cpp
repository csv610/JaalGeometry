#include "ObjectsListDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JObjectsListDialog :: JObjectsListDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
    tableWidget->setColumnCount(2);

    QStringList labels;
    labels << tr("ObjectName");
    labels << tr("Filename");
    tableWidget->setHorizontalHeaderLabels(labels);
    tableWidget->setShowGrid(true);
    tableWidget->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);

    QPalette *palette = new QPalette();
    palette->setColor(QPalette::Highlight, Qt::lightGray);
    tableWidget->setPalette( *palette);

    meshViewer  = nullptr;
    imageViewer = nullptr;
    listObjects[0] = 1;
    listObjects[1] = 1;
}

///////////////////////////////////////////////////////////////////////////////

JObjectsListDialog :: ~JObjectsListDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JObjectsListDialog :: setType(int t )
{
    listObjects[0] = 1;
    listObjects[1] = 1;

    if( t == MESH_OBJECTS  ) listObjects[1] = 0;
    if( t == IMAGE_OBJECTS ) listObjects[0] = 0;

    filltable();
}

///////////////////////////////////////////////////////////////////////////////

void JObjectsListDialog :: filltable()
{
    int numMesh = 0, numImage = 0;

    if( meshViewer && listObjects[0] ) {
        numMesh = meshViewer->getSize();
    }

    if( imageViewer && listObjects[1] )
        numImage = imageViewer->getSize();

    tableWidget->setRowCount(numMesh + numImage);

    int index = 0;
    if( listObjects[0] ) {
        for( int i = 0; i < numMesh; i++) {
            const JMeshPtr &mesh = meshViewer->getMesh(i);
            string name = mesh->getName();
            QTableWidgetItem *item1 = new QTableWidgetItem(name.c_str() );
            if(mesh->isActive() )
                item1->setCheckState(Qt::Checked);
            else
                item1->setCheckState(Qt::Unchecked);
            tableWidget->setItem(index, 0, item1);

            string fname = mesh->getFileName();
            if( !fname.empty() ) {
                QTableWidgetItem *item2 = new QTableWidgetItem(fname.c_str() );
                tableWidget->setItem(index, 1, item2);
            }
            index++;
        }
    }

    if( listObjects[1] ) {
        for( int i = 0; i < numImage; i++) {
            JImagePtr img = imageViewer->getImage(i);
            string name = img->getName();
            QTableWidgetItem *item = new QTableWidgetItem(name.c_str() );
            item->setCheckState(Qt::Checked);
            tableWidget->setItem(index, 0, item);
            index++;
        }
    }
    tableWidget->resizeColumnsToContents();
}
///////////////////////////////////////////////////////////////////////////////

void JObjectsListDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer  = dynamic_pointer_cast<JMeshViewer>(c);

    c = viewManager->getComponent("ImageViewer");
    imageViewer = dynamic_pointer_cast<JImageViewer>(c);

    filltable();
}

///////////////////////////////////////////////////////////////////////////////

void JObjectsListDialog :: getCell(int i, int j)
{
    currMesh  = nullptr;
    currImage = nullptr;
    QTableWidgetItem *item = tableWidget->item(i,j);

    bool bitVal = 0;
    switch( item->checkState() ) {
    case Qt::Checked:
        bitVal = 1;
        break;
    case Qt::Unchecked:
        bitVal = 0;
        break;
    }

    JMeshRenderPtr mrender;
    if( listObjects[0] )  {
        item = tableWidget->item(i,0);
        string meshname = item->text().toUtf8().constData();
        currMesh = meshViewer->getMesh(meshname);
        currMesh->getAttribute("Render", mrender);
        mrender->display = bitVal;
        meshViewer->setCurrentMesh(currMesh);
    }

    /*
        if( item->text() == QString("Image")) {
            item = tableWidget->item(i,1);
            string imagename = item->text().toUtf8().constData();
            currImage = imageViewer->getImage(imagename);
            currImage->setActiveBit(bitVal);
        }
    */

    viewManager->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JObjectsListDialog :: deleteObject()
{
    if( currMesh ) {
        meshViewer->removeObject(currMesh);
        currMesh = nullptr;
    }

    /*
        if( currImage ) {
            imageViewer->remove(currImage);
            currImage = nullptr;
        }
    */

    filltable();
}

///////////////////////////////////////////////////////////////////////////////
void
JObjectsListDialog :: saveMesh()
{
    if( currMesh == nullptr) return;

//  static QString lastSelectedDirectory;

    QString qstr;
    string name = currMesh->getFileName();

    if( name.empty() ) {
        name = "unnamed_mesh.off";
    }

    if( currMesh->getTopology()->getDimension() == 3) {
        qstr  = QFileDialog::getSaveFileName(this,
                                             *new QString("Select Mesh File "),
                                             QString::fromStdString( name ),
                                             *new QString( "Mesh Format (*.xml *.vtk *.ele)"));
    }

    if( currMesh->getTopology()->getDimension() == 2) {
        qstr  = QFileDialog::getSaveFileName(this,
                                             *new QString("Select Mesh File "),
                                             QString::fromStdString( name ),
                                             *new QString( "Mesh Format (*.off *.obj *.xml *.vtk *.ele)"));
    }

    if( currMesh->getTopology()->getDimension() == 1) {
        qstr  = QFileDialog::getSaveFileName(this,
                                             *new QString("Select Mesh File "),
                                             QString::fromStdString( name ),
                                             *new QString( "Mesh Format (*.off *.xml *.vtk)"));
    }

    JWaitCursor waitCursor;
    waitCursor.start();

    string meshFileName = StdString(qstr);
    if (!meshFileName.empty()) {
        JMeshIO::saveAs(currMesh, meshFileName);
        currMesh->setFileName(meshFileName);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JObjectsListDialog :: saveObject()
{
    if( currMesh ) {
        saveMesh();
    }

    /*
        if( currImage ) {
            saveImage();
        }
    */

}
///////////////////////////////////////////////////////////////////////////////
void JObjectsListDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->close();
}
///////////////////////////////////////////////////////////////////////////////
void JObjectsListDialog :: renameObject()
{
    QString qstr = renameLineEdit->text();
    if( currMesh) currMesh->setName( StdString(qstr)) ;
    filltable();
}

///////////////////////////////////////////////////////////////////////////////
void JObjectsListDialog :: makeConnections()
{
    connect( tableWidget,  SIGNAL( cellClicked(int,int) ), this, SLOT( getCell(int, int) ));

    PushButton( savePushButton,   [=] {saveObject();});
    PushButton( renamePushButton, [=] {renameObject();});
    PushButton( deletePushButton, [=] {deleteObject();});
    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
