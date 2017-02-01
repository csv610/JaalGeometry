#include "SingletDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
SingletColor ::  SingletColor()
{
    highlightColor[0] = 1.0;
    highlightColor[1] = 0.0;
    highlightColor[2] = 0.0;
    highlightColor[3] = alpha;

    defaultColor[0] = 0.0;
    defaultColor[1] = 1.0;
    defaultColor[2] = 0.0;
    defaultColor[3] = alpha;
}
///////////////////////////////////////////////////////////////////////////////

int SingletColor::assign(const JNodePtr &vertex)
{
    JNodeRenderPtr nAttrib;
    vertex->getAttribute("Render", nAttrib);

    if( JSinglet::isSinglet(vertex) )
        nAttrib->color =  highlightColor;
    else
        nAttrib->color =  nosingletColor;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JSingletDialog :: JSingletDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;

    singlets.reset( new JSinglet);
    singletColor.reset( new SingletColor );
}

///////////////////////////////////////////////////////////////////////////////

JSingletDialog :: ~JSingletDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JSingletDialog :: setColor()
{
    JColor rgb;
    QColor color = QColorDialog::getColor();
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    rgb[3] = 1.0;
    if( singletColor) singletColor->setHighlightColor(rgb);

    if( meshViewer == nullptr ) return;
    meshViewer->getViewManager()->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

int JSingletDialog :: isQuadMesh()
{
    if( mesh == nullptr ) return 0;

    int dim = mesh->getTopology()->getDimension();
    if( dim == 2 ) {
        int etype = mesh->getTopology()->getElementsType(2);
        if( etype != 4 ) {
            QMessageBox msg;
            msg.setIcon(QMessageBox::Warning);
            msg.setText("At present singlet operations are only for all quad mesh ");
            msg.setStandardButtons( QMessageBox::Ok);
            int ret = msg.exec();
            if( ret == QMessageBox::Ok ) {
                return 0;
            }
        }
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////
void JSingletDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;
    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////

void JSingletDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    if( !isQuadMesh() ) {
        QMessageBox msg;
        msg.setIcon(QMessageBox::Warning);
        msg.setText("At present singlet operations are only for all quad mesh ");
        msg.setStandardButtons( QMessageBox::Ok);
        int ret = msg.exec();
        if( ret == QMessageBox::Ok ) return;
    }

    mesh->buildRelations(0,2);
    singlets->setMesh(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JSingletDialog :: removeAll()
{
    if( singlets == nullptr) return;
    singlets->removeAll();

    searchSinglets();
}

///////////////////////////////////////////////////////////////////////////////
void JSingletDialog :: searchSinglets()
{
    JNodeSequence nodes = singlets->getSinglets();
    int nSize = nodes.size();
    numSingletsLineEdit->setText( QString::number(nSize) );

    size_t numnodes = mesh->getSize(0);

    JNodeRenderPtr vAttrib;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        if( vtx->isActive() ) {
            vtx->getAttribute("Render", vAttrib);
            vAttrib->display = 0;
            vAttrib->glyph   = 0;
        }
    }

    JColor highlight = JEntityColor::getColor("Red");
    for( const JNodePtr &vtx : nodes) {
        vtx->getAttribute("Render", vAttrib);
        vAttrib->display = 1;
        vAttrib->glyph   = 1;
        vAttrib->color   = highlight;
        cout << vtx->getID() << endl;
    }

    meshViewer->updateBuffers(mesh);

}
///////////////////////////////////////////////////////////////////////////////

void JSingletDialog :: makeConnections()
{
    PushButton( searchPushButton,    [=] {searchSinglets(); });
    PushButton( singletColorButton,  [=] {setColor();});
    PushButton( removeAllPushButton, [=] {removeAll();});
    PushButton( closePushButton,     [=] {close(); });
}

///////////////////////////////////////////////////////////////////////////////
