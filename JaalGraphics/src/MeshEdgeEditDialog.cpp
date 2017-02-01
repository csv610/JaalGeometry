#include "MeshEdgeEditDialog.hpp"

///////////////////////////////////////////////////////////////////////////////

JMeshEdgeEditDialog :: JMeshEdgeEditDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();

    mesh = nullptr;
    viewManager = nullptr;
    meshViewer  = nullptr;
}

///////////////////////////////////////////////////////////////////////////////

JMeshEdgeEditDialog :: ~JMeshEdgeEditDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgeEditDialog :: setColor()
{
    /*
        JColor rgb;
        QColor color = QColorDialog::getColor();
        rgb[0] = color.red()/255.0;
        rgb[1] = color.green()/255.0;
        rgb[2] = color.blue()/255.0;
        rgb[3] = 1.0;
        if( singletColor) singletColor->setHighlightColor(rgb);

        if( meshViewer == nullptr ) return;
        meshViewer->getViewManager()->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////

void JMeshEdgeEditDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    if( meshViewer == nullptr ) return;

    picker = meshViewer->getEntityPicker();
    if( picker ) picker->setMode(2);

    setMesh( meshViewer->getCurrentMesh() );
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgeEditDialog :: setInfo()
{
    size_t numEdges = mesh->getSize(1);

    vector<double> elen;
    elen.reserve(numEdges);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            double l = JEdgeGeometry::getLength(edge);
            elen.push_back(l);
        }
    }
    double minlen = *boost::min_element(elen);
    double maxlen = *boost::max_element(elen);


    minLengthLineEdit->setText(QString::number(minlen));
    maxLengthLineEdit->setText(QString::number(maxlen));
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEdgeEditDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr ) return ;

    string name = mesh->getName();
    objectNameLineEdit->setText(QString(name.c_str()));

    JMeshRenderPtr  mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 1;

    setInfo();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgeEditDialog :: getMinEdges()
{
    if( mesh == nullptr) return;
    QString qstr;

    JColor red  = JEntityColor::getColor("Red");
    JColor gray = JEntityColor::getColor("Gray");

    size_t numEdges = mesh->getSize(1);
    JEdgeRenderPtr   eAttrib;
    qstr = getMinLenLineEdit->text();
    double minlen  = qstr.toDouble();
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", eAttrib);
            eAttrib->color = gray;
            eAttrib->scale = 1.0;
            double l = JEdgeGeometry::getLength(edge);
            if( l <= minlen ) {
                eAttrib->color = red;
                eAttrib->scale = 2.0;
            }
        }
    }

    meshViewer->updateBuffers(mesh);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgeEditDialog :: getMaxEdges()
{
    if( mesh == nullptr) return;

    QString qstr;

    JColor red  = JEntityColor::getColor("Red");
    JColor gray = JEntityColor::getColor("Gray");

    size_t numEdges = mesh->getSize(1);
    JEdgeRenderPtr   eAttrib;
    qstr = getMaxLenLineEdit->text();
    double maxlen  = qstr.toDouble();

    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->isActive() ) {
            edge->getAttribute("Render", eAttrib);
            eAttrib->color = gray;
            eAttrib->scale = 1.0;
            double l = JEdgeGeometry::getLength(edge);
            if( l >= maxlen ) {
                eAttrib->color = red;
                eAttrib->scale = 2.0;
            }
        }
    }
    meshViewer->updateBuffers(mesh);
}
///////////////////////////////////////////////////////////////////////////////
void JMeshEdgeEditDialog :: refine()
{

}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgeEditDialog :: collapse()
{
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgeEditDialog :: deleteInternal()
{
    /*
        if( mesh == nullptr ) return;

        mesh->deleteEdges( JMeshEntity::ANY_ENTITY );
        meshViewer->displayAll(1, 0);

        size_t numedges = mesh->getSize(1);
        numEdgesLineEdit->setText( QString::number(numedges) );

        meshViewer->refreshDisplay();
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgeEditDialog :: deleteAll()
{
    /*
        if( mesh == nullptr ) return;

        mesh->deleteEdges( JMeshEntity::ANY_ENTITY );
        meshViewer->displayAll(1, 0);

        size_t numedges = mesh->getSize(1);
        numEdgesLineEdit->setText( QString::number(numedges) );

        meshViewer->refreshDisplay();
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshEdgeEditDialog :: closeDialog()
{
    this->parentWidget()->show();
    this->hide();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshEdgeEditDialog :: makeConnections()
{
    connect( refinePushButton, SIGNAL( clicked() ), this, SLOT( refine() ));
    connect( minLengthPushButton, SIGNAL( clicked() ), this, SLOT( getMinEdges() ));
    connect( maxLengthPushButton, SIGNAL( clicked() ), this, SLOT( getMaxEdges() ));
    connect( closePushButton, SIGNAL( clicked() ), this, SLOT( close() ));
}

///////////////////////////////////////////////////////////////////////////////
