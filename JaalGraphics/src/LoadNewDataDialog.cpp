#include "LoadNewDataDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JLoadNewDataDialog :: JLoadNewDataDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    viewManager = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JLoadNewDataDialog :: ~JLoadNewDataDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JLoadNewDataDialog :: init()
{
}

///////////////////////////////////////////////////////////////////////////////

void JLoadNewDataDialog :: loadNewMesh()
{
    if( viewManager == nullptr) return;

    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    if( c == nullptr)
        c = JMeshViewer::registerComponent(viewManager);

    JMeshViewerPtr meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Mesh File "),
                   lastSelectedDirectory,
                   *new QString( "Mesh Format (*.xml *.off *.obj *.node *.ele  *.vtk *cub  *.mesh *.msh *.smf *.poly)"));

    string meshFileName = qstr.toUtf8().constData();
    if (meshFileName.empty()) return;

    JWaitCursor cursor;
    cursor.start();

    ThreadWork1 work;
    work.name = meshFileName;
    work.meshViewer = meshViewer;
    work.run();

    /*
    //  readThread.reset( new std::thread(work));
        JMeshPtr msh = JMeshIO::readFile(meshFileName);
        if( msh) {
            objectnameLineEdit->setText( QString::fromStdString(msh->getName()));
            msh->setFileName( meshFileName);
            meshViewer->addObject(msh);
        }
    */
//    QApplication::restoreOverrideCursor();
}

///////////////////////////////////////////////////////////////////////////////

void JLoadNewDataDialog :: loadNewImage()
{
    static int imageCounter = 1;

    JViewComponentPtr c = viewManager->getComponent("ImageViewer");
    if( c == nullptr)
        c = JImageViewer::registerComponent(viewManager);

    JImageViewerPtr imageViewer = dynamic_pointer_cast<JImageViewer>(c);
    if( imageViewer == nullptr ) return;

    QString qstr = QFileDialog::getOpenFileName(this,
                   *new QString("Select Image File "),
                   lastSelectedDirectory,
                   *new QString( "Image Format (*.jpg *.JPG *.png *.gif *tif)"));

    string imageFileName = StdString(qstr);

    JImagePtr image = imageViewer->getImage(imageFileName );
    if( image ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("An image already exists. ");
        msg.setInformativeText("Do you want to replace it ? ");
        msg.setStandardButtons( QMessageBox::Ok | QMessageBox::Cancel);
        int ret = msg.exec();
        if( ret == QMessageBox::Cancel ) return;
    }

    if (!imageFileName.empty()) {
        imageViewer->addImage(imageFileName);
    }

    imageCounter++;
}

////////////////////////////////////////////////////////////////////////////////

void JLoadNewDataDialog :: loadNewData()
{
    if( meshRadioButton->isChecked() )  loadNewMesh();
    if( imageRadioButton->isChecked() ) loadNewImage();
    this->hide();
}

////////////////////////////////////////////////////////////////////////////////
void JLoadNewDataDialog :: closeDialog()
{
    objectnameLineEdit->clear();
    this->hide();
}
////////////////////////////////////////////////////////////////////////////////

void JLoadNewDataDialog :: makeConnections()
{
    PushButton(  loadPushButton,  [=] {loadNewData();});
    PushButton(  closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
