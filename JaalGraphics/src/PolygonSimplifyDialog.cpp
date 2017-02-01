#include "PolygonSimplifyDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JPolygonSimplifyDialog :: JPolygonSimplifyDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
//  simplify.reset( new JPolygonSimplify);
    toleranceLineEdit->setText(QString::number(0.001));
}

///////////////////////////////////////////////////////////////////////////////

JPolygonSimplifyDialog :: ~JPolygonSimplifyDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JPolygonSimplifyDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
}

///////////////////////////////////////////////////////////////////////////////

void JPolygonSimplifyDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    vector<JEdgeSequence> boundedges;
    mesh->getTopology()->getBoundary(bounedges);
    int nsize = boundedges.size();
    
/*
    int numnodes = bndnodes.size();
    numNodes1LineEdit->setText(QString::number(numnodes));

//    simplify->setMesh(mesh);

    JColor black = JEntityColor::getColor("Black");
    JNodeRenderPtr vAttrib;
    for( const JNodePtr &vtx : bndnodes) {
        vtx->getAttribute("Render", vAttrib);
        vAttrib->glyph = 0;
        vAttrib->pointSize = 3;
        vAttrib->color = black;
    }
*/
}

///////////////////////////////////////////////////////////////////////////////
void JPolygonSimplifyDialog :: getSimplified()
{
/*
    QString  qstr = algoComboBox->currentText();
    string   str  = StdString(qstr);
    int      algo = JPolygonSimplify::getAlgorithmID(str);
    if( algo < 0) return;

    if( simplifiedMesh ) {
        simplifiedMesh->deleteAll();
        meshViewer->removeObject( simplifiedMesh );
    }

    qstr = toleranceLineEdit->text();
    double tol = qstr.toDouble();
    simplifiedMesh = simplify->getSimplified(algo, tol);
    meshViewer->addObject(simplifiedMesh);

    JNodeSequence bndnodes;
    simplifiedMesh->getTopology()->getBoundary(bndnodes);
    int numnodes = bndnodes.size();
    numNodes2LineEdit->setText(QString::number(numnodes));

    JColor red = JEntityColor::getColor("Red");
    JNodeRenderPtr vAttrib;
    for( const JNodePtr &vtx : bndnodes) {
        vtx->getAttribute("Render", vAttrib);
        vAttrib->glyph = 0;
        vAttrib->pointSize = 5;
        vAttrib->color = red;
    }

    JEdgeSequence bndedges;
    simplifiedMesh->getTopology()->getBoundary(bndedges);
    JColor blue = JEntityColor::getColor("Blue");

    JEdgeRenderPtr eAttrib;
    for( const JEdgePtr &edge : bndedges) {
        edge->getAttribute("Render", eAttrib);
        eAttrib->lineWidth = 2.0;
        eAttrib->color  = blue;
    }

    meshViewer->updateBuffers( simplifiedMesh );
    displayMesh();
*/
}

///////////////////////////////////////////////////////////////////////////////

void JPolygonSimplifyDialog :: displayMesh()
{
    bool val;
    val = displayOriginalCheckBox->isChecked();
    meshViewer->displayAll(mesh, val);

    if( simplifiedMesh ) {
        val = displaySimplifiedCheckBox->isChecked();
        meshViewer->displayAll(simplifiedMesh, val);
    }

    meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JPolygonSimplifyDialog :: acceptMesh()
{
}

///////////////////////////////////////////////////////////////////////////////
void JPolygonSimplifyDialog :: closeDialog()
{
/*
    if( simplifiedMesh ) {
        simplifiedMesh->deleteAll();
        meshViewer->removeObject( simplifiedMesh );
        meshViewer->refreshDisplay();
    }
    close();
*/
}
///////////////////////////////////////////////////////////////////////////////
void JPolygonSimplifyDialog :: makeConnections()
{
    CheckBox( displayOriginalCheckBox,   [=] {displayMesh();});
    CheckBox( displaySimplifiedCheckBox, [=] {displayMesh();});

    PushButton( applyPushButton,  [=] {getSimplified();});
    PushButton( acceptPushButton, [=] {acceptMesh();});
    PushButton( closePushButton,  [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
