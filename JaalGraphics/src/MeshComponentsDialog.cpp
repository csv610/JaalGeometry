#include "MeshComponentsDialog.hpp"


///////////////////////////////////////////////////////////////////////////////
JMeshComponentsDialog :: JMeshComponentsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    jc.reset( new JMeshComponents() );
}

///////////////////////////////////////////////////////////////////////////////

JMeshComponentsDialog :: ~JMeshComponentsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshComponentsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshComponentsDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    jc->setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshComponentsDialog :: searchComponents()
{
    JWaitCursor waitCursor;
    waitCursor.start();

    jc->searchComponents();

    int numComponents = jc->getNumComponents();

    numComponentsLineEdit->setText( QString::number(numComponents) );
    displayComponentSpinBox->setMaximum( numComponents-1);
    aCompSpinBox->setMaximum( numComponents-1);
    bCompSpinBox->setMaximum( numComponents-1);

    if( numComponents < 2) return;

    if( numComponents > colors.size() ) {
        colors.resize(numComponents);
        for( int i = 0; i < numComponents; i++)
            colors[i] = JEntityColor::getRandomColor();
    }

    JFaceRenderPtr fAttrib;
    for( int i = 0; i < numComponents; i++) {
        JMeshPtr submesh = jc->getComponent(i);
        if( submesh) {
            size_t numfaces = submesh->getSize(2);
            JFaceRenderPtr  attrib;
            for( size_t j = 0; j < numfaces; j++) {
                const JFacePtr &face = submesh->getFaceAt(j);
                face->getAttribute("Render", fAttrib);
                fAttrib->color = colors[i];
            }
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshComponentsDialog :: mergeComponents()
{
}
///////////////////////////////////////////////////////////////////////////////
void JMeshComponentsDialog :: removeComponent()
{
    /*
        int id = displayComponentSpinBox->value();
        JMeshPtr submesh =  jc->getComponent(id);
        if( submesh == nullptr) return;
        mesh->getTopology()->removeSubmesh(submesh);
        meshViewer->updateBuffers(mesh);
    */
}
///////////////////////////////////////////////////////////////////////////////
void JMeshComponentsDialog :: mergeAll()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMeshComponentsDialog :: displayComponent()
{
    cout << "HELLO " << endl;
    if( mesh == nullptr) return;

    size_t numfaces = mesh->getSize(2);
    JFaceRenderPtr  fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display  = 0;
    }

    int id = displayComponentSpinBox->value();

    JMeshPtr submesh =  jc->getComponent(id);
    if( submesh == nullptr) return;

    numfaces = submesh->getSize(2);

    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = submesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display  = 1;
    }

    numNodesLineEdit->setText( QString::number(submesh->getSize(0)) );
    numEdgesLineEdit->setText( QString::number(submesh->getSize(1)) );
    numFacesLineEdit->setText( QString::number(submesh->getSize(2)) );
    numCellsLineEdit->setText( QString::number(submesh->getSize(3)) );

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshComponentsDialog :: select2Components()
{
    if( mesh == nullptr) return;

    size_t numfaces = mesh->getSize(2);
    JFaceRenderPtr  fAttrib;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display  = 0;
    }

    JMeshPtr submesh;
    int id1 = aCompSpinBox->value();
    submesh =  jc->getComponent(id1);
    if( submesh == nullptr) return;

    numfaces = submesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = submesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display  = 1;
    }

    int id2 = bCompSpinBox->value();
    submesh =  jc->getComponent(id2);
    if( submesh == nullptr) return;

    numfaces = submesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = submesh->getFaceAt(i);
        face->getAttribute("Render", fAttrib);
        fAttrib->display  = 1;
    }

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshComponentsDialog :: closeDialog()
{
    this->close();
    parentWidget()->show();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshComponentsDialog :: makeConnections()
{
    PushButton( searchPushButton,  [=] {searchComponents(); });
    PushButton( removePushButton,  [=] {removeComponent(); });
    SpinBoxi (  displayComponentSpinBox, [=] {displayComponent(); });
    SpinBoxi (  aCompSpinBox, [=] {select2Components(); });
    SpinBoxi (  bCompSpinBox, [=] {select2Components(); });
    PushButton( closePushButton, [=] { closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
