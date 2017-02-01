#include "MeshNormalClustersDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
JMeshNormalClustersDialog :: JMeshNormalClustersDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
}

///////////////////////////////////////////////////////////////////////////////

JMeshNormalClustersDialog :: ~JMeshNormalClustersDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalClustersDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalClustersDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));
    normalClusters.reset( new JMeshNormalClusters);
    normalClusters->setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNormalClustersDialog :: getClusters()
{
    if( normalClusters == nullptr) return;

    int nlevel = binLevelSpinBox->value();
    normalClusters->setBinLevels(nlevel);
    normalClusters->createClusters();

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 1;
    mrender->displayEntity[3] = 0;

    JFacePartitionColor faceColor;
    faceColor.setMesh(mesh);

    JEdgePartitionColor edgeColor;
    edgeColor.setMesh(mesh);

    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshNormalClustersDialog :: openNormalsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JMeshNormalClustersDialog :: showSphere()
{
    meshViewer->removeObject( sphMesh);

    if( mesh == nullptr) return;

    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->displayEntity[0] = 0;
    mrender->displayEntity[1] = 0;
    mrender->displayEntity[2] = 1;
    mrender->displayEntity[3] = 0;

    meshViewer->refreshDisplay();

    if( !showSphereCheckBox->isChecked() ) return;

    int nlevel = binLevelSpinBox->value();
    normalClusters->setBinLevels(nlevel);
    sphMesh = normalClusters->getQuantizedSphere();
    meshViewer->addObject(sphMesh); 

    JFacePartitionColor faceColor;
    faceColor.setMesh(sphMesh);

    mrender->displayEntity[2] = 0;  // Mesh rendering switched off...
    meshViewer->updateBuffers(sphMesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNormalClustersDialog :: closeDialog()
{
    this->close();
    parentWidget()->show();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshNormalClustersDialog :: makeConnections()
{
    CheckBox( showSphereCheckBox,   [=] { showSphere(); });
    PushButton( applyPushButton,    [=] { getClusters(); });
    PushButton( getNormalsPushButton,  [=] {openNormalsDialog();});
    PushButton( closePushButton,     [=] {closeDialog(); });
}

///////////////////////////////////////////////////////////////////////////////
