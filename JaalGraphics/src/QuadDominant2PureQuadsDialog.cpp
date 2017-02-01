#include "QuadDominant2PureQuadsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JQuadDominant2PureQuadsDialog :: JQuadDominant2PureQuadsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    pureQuadsMesher.reset( new JQuadDominant2PureQuadsMesher);
    srcSpinBox->setValue( -1 );
    dstSpinBox->setValue( -1 );
}

///////////////////////////////////////////////////////////////////////////////

JQuadDominant2PureQuadsDialog :: ~JQuadDominant2PureQuadsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsDialog :: setColors()
{
    if( mesh == nullptr) return;

    size_t numfaces = mesh->getSize(2);

    JColor white = JEntityColor::getColor("White");
    JColor red  = JEntityColor::getColor("Red");
    JColor blue = JEntityColor::getColor("Blue");

    size_t numTris   = 0;
    size_t numPolys  = 0;
    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            int err = face->getAttribute("Render", fAttrib);
            int nn = face->getSize(0);
            fAttrib->color = white;
            if( nn != 4) {
                fAttrib->display = 1;
                if( nn == 3) {
                    fAttrib->color   = red;
                    numTris++;
                } else {
                    fAttrib->color   = blue;
                    numPolys++;
                }
            }
        }
    }
    numTrianglesLineEdit->setText(QString::number(numTris));
    numPolygonsLineEdit->setText(QString::number(numPolys));
    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    viewManager->attach(this);
    entityPicker = meshViewer->getEntityPicker();
    entityPicker->setMode(1);
}
///////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    int elemType = mesh->getTopology()->getElementsType(2);
    if( elemType == JFace::QUADRILATERAL) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Information);
        msg.setText("All elements are quadrilateral");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    setColors();
    pureQuadsMesher->setMesh(mesh);

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 2;
    entityPicker->setMesh(mesh);
    size_t numfaces = mesh->getSize(2);
    srcSpinBox->setMaximum(numfaces-1);
    dstSpinBox->setMaximum(numfaces-1);
}

///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog :: mouseReleaseEvent(QMouseEvent *e)
{
    if( !manualCheckBox->isChecked() )  return;
    JFaceSequence faces = entityPicker->getPickedFaces();

    if( faces.empty() ) return;

    if( srcRadioButton->isChecked() ) {
        if( faces[0]->getSize(0) != 4)
            srcSpinBox->setValue( faces[0]->getID() );
    }

    if( dstRadioButton->isChecked() ) {
        if( faces[0]->getSize(0) != 4)
            dstSpinBox->setValue( faces[0]->getID() );
    }

}
///////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsDialog :: getAllStrips()
{
    JWaitCursor waitCursor;
    waitCursor.start();
    vector<JFaceSequence> allStrips = pureQuadsMesher->getAllStrips();

    setColors();

    JColor green = JEntityColor::getColor("Green");
    JFaceRenderPtr fAttrib;

    for( size_t j = 0; j < allStrips.size(); j++) {
        newStrip = allStrips[j];
        int numfaces = newStrip.size();
        if( numfaces >= 2) {
            for( int i = 1; i < numfaces-1; i++) {
                const JFacePtr &face = newStrip[i];
                int err = face->getAttribute("Render", fAttrib);
                fAttrib->color = green;
            }
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog :: getNewStrip()
{
    if( pureQuadsMesher == nullptr) return;

    int fid;
    JFacePtr fsrc, fdst;
    if( manualCheckBox->isChecked() ) {
        fid = srcSpinBox->value();
        if( fid >= 0) fsrc = mesh->getFaceAt(fid);
        fid = dstSpinBox->value();
        if( fid >= 0) fdst = mesh->getFaceAt(fid);
        if( (fsrc != nullptr)  && (fdst == nullptr)) {
            newStrip = pureQuadsMesher->getNewStrip(fsrc);
        }
        if( (fsrc != nullptr)  && (fdst != nullptr))
            newStrip = pureQuadsMesher->getNewStrip(fsrc, fdst);
    } else {
        newStrip = pureQuadsMesher->getNewStrip();
    }

    size_t numfaces = newStrip.size();
    if( numfaces < 2) return;
    setColors();

    JColor green = JEntityColor::getColor("Green");
    JFaceRenderPtr fAttrib;
    for( int i = 1; i < numfaces-1; i++) {
        const JFacePtr &face = newStrip[i];
        int err = face->getAttribute("Render", fAttrib);
        fAttrib->color = green;
    }
    meshViewer->alignAlong(mesh, newStrip[0] );
    meshViewer->updateBuffers(mesh);

}
///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog :: remeshStrip()
{
    /*
        if( pureQuadsMesher == nullptr) return;
        pureQuadsMesher->remeshCurrentStrip();
        faceSeq.clear();
        setColors();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDominant2PureQuadsDialog :: remeshAllStrips()
{
    /*
        if( pureQuadsMesher == nullptr) return;
        size_t numNonQuads = pureQuadsMesher->getNumOfNonQuads();

        for( size_t i = 0; i < 50; i++) {
            getMergeStrip();
            remeshCurrentStrip();
        }
    */
}
///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog ::  refineAll()
{
}
///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog ::  enumFaces()
{

}
///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog :: clearAll()
{

    size_t numfaces = mesh->getSize(2);

    JColor white = JEntityColor::getColor("White");
    JFaceRenderPtr fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( face->isActive() ) {
            face->getAttribute("Render", fAttrib);
            fAttrib->color   = white;
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog ::  closeDialog()
{
    this->close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////
void JQuadDominant2PureQuadsDialog :: makeConnections()
{
    PushButton( getNewStripPushButton, [=] {getNewStrip();});
    PushButton( remeshStripPushButton, [=] {remeshStrip();});
    PushButton( getAllStripsPushButton, [=] {getAllStrips();});
    PushButton( remeshAllStripsPushButton, [=] {remeshAllStrips();});
    PushButton( refineAllPushButton, [=] {refineAll();});
    PushButton( clearPushButton, [=] {clearAll();});

    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
