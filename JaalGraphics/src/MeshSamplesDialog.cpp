#include "MeshSamplesDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshSamplesDialog :: JMeshSamplesDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
    percentSpinBox->setValue(0.01);
    meshSamples.reset( new JMeshSamples);
}

///////////////////////////////////////////////////////////////////////////////

JMeshSamplesDialog :: ~JMeshSamplesDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSamplesDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSamplesDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    double val  = percentSpinBox->value();
    if( nodeRadioButton->isChecked() ) numSamples = val*mesh->getSize(0);
    if( edgeRadioButton->isChecked() ) numSamples = val*mesh->getSize(1);
    if( faceRadioButton->isChecked() ) numSamples = val*mesh->getSize(2);
    if( cellRadioButton->isChecked() ) numSamples = val*mesh->getSize(3);
    numSamplesLineEdit->setText( QString::number(numSamples));

    meshSamples->setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshSamplesDialog :: getSamples()
{
    QString qstr =  algoComboBox->currentText();
    string  str  =  StdString(qstr);
    if( str == "Random")     algo = JMeshSamples::RANDOM;
    if( str == "Geodesic")   algo = JMeshSamples::GEODESIC;
    if( str == "Biharmonic") algo = JMeshSamples::BIHARMONIC;
    if( str == "Octree")     algo = JMeshSamples::BIHARMONIC;
    if( str == "Regular")    algo = JMeshSamples::REGULAR;
    meshSamples->setAlgorithm(algo);

    qstr = numSamplesLineEdit->text();
    numSamples = qstr.toInt();

    if( nodeRadioButton->isChecked() ) getNodeSamples();
    if( edgeRadioButton->isChecked() ) getEdgeSamples();
    if( faceRadioButton->isChecked() ) getFaceSamples();
    if( cellRadioButton->isChecked() ) getCellSamples();

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSamplesDialog :: getNodeSamples()
{
    JWaitCursor waitCursor;
    waitCursor.start();
    JNodeSequence nodes = meshSamples->getNodeSamples(numSamples);
    numSamples = nodes.size();
    numSamplesLineEdit->setText( QString::number(numSamples) );

    JNodeRenderPtr vAttrib;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", vAttrib);
            vAttrib->display = 0;
        }
    }

    JColor highlight;
    highlight[0] = 1.0;
    highlight[1] = 0.0;
    highlight[2] = 0.0;
    highlight[2] = 1.0;
    for( const JNodePtr &vtx : nodes) {
        vtx->getAttribute("Render", vAttrib);
        vAttrib->display = 1;
        vAttrib->color   = highlight;
        vAttrib->glyph   = 1;
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshSamplesDialog :: getEdgeSamples()
{
}

void JMeshSamplesDialog :: getFaceSamples()
{
}

void JMeshSamplesDialog :: getCellSamples()
{
}

void JMeshSamplesDialog :: makeConnections()
{
    PushButton( applyPushButton, [=] {getSamples();});
    PushButton( closePushButton, [=] {close();});
}

///////////////////////////////////////////////////////////////////////////////

