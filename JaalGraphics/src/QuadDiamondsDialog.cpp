#include "QuadDiamondsDialog.hpp"

///////////////////////////////////////////////////////////////////////////////
QuadDiamondColor ::  QuadDiamondColor()
{
    highlightColor[0] = 1.0;
    highlightColor[1] = 0.0;
    highlightColor[2] = 0.0;
    highlightColor[3] = 1.0;

    defaultColor[0] = 0.0;
    defaultColor[1] = 1.0;
    defaultColor[2] = 0.0;
    defaultColor[3] = 1.0;
}
///////////////////////////////////////////////////////////////////////////////

int QuadDiamondColor :: assign( const JFacePtr &face)
{
    int pos, type = 33;

    JFaceRenderPtr fAttrib;
    face->getAttribute("Render", fAttrib);

    if(JDiamond::isDiamond(face, pos, type) ) {
        fAttrib->color = highlightColor;
    } else
        fAttrib->color = defaultColor;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

JQuadDiamondsDialog :: JQuadDiamondsDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    mesh = nullptr;

//    diamondColor.reset( new QuadDiamondColor);
//    irregularColor.reset(new IrregularNodeColor);
//    meshoptDialog.reset(new JMeshOptDialog(this));
}

///////////////////////////////////////////////////////////////////////////////

JQuadDiamondsDialog :: ~JQuadDiamondsDialog()
{
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDiamondsDialog :: init()
{
    if( viewManager == nullptr ) return;
    JViewComponentPtr c = viewManager->getComponent("MeshViewer");
    meshViewer = dynamic_pointer_cast<JMeshViewer>(c);

    if( meshViewer == nullptr ) return;

//  mesh = meshViewer->getMesh();
    if( mesh == nullptr ) return ;

//    meshViewer->displayAll(2,1);
//    meshViewer->displayAll(0,1);
//  meshViewer->getDrawFace()->setShade(DrawFace::FLAT_SHADE);

    searchDiamonds();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDiamondsDialog :: searchDiamonds()
{
    /*
        if( mesh == nullptr ) return;

        mesh->buildRelations(0,2);
        size_t nSize = 0;
        if( diamond33CheckBox->isChecked() )
            nSize += Diamond::getSize(mesh, 33);

        if( diamond35CheckBox->isChecked() )
            nSize += Diamond::getSize(mesh, 35);

        numDiamondsLineEdit->setText( QString::number(nSize) );

        JFaceColor::assign(mesh, diamondColor);
        JNodeColor::assign(mesh, irregularColor);
    */
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDiamondsDialog :: setColor()
{
    QColor color = QColorDialog::getColor();

    JColor rgb;
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    rgb[3] = 1.0;

    diamondColor->setHighlightColor( rgb );

    if( meshViewer ) meshViewer->refreshDisplay();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDiamondsDialog :: removeAll()
{
    if( mesh == nullptr ) return ;

    mesh->buildRelations(0,2);

    JQuadCleanUp qClean;
    qClean.setMesh(mesh);
    /*
    qClean.removeDiamonds();

    //  meshViewer->setNewMesh(mesh);

        size_t nSize = 0;
        if( diamond33CheckBox->isChecked() )
            nSize += Diamond::getSize(mesh, 33);

        if( diamond35CheckBox->isChecked() )
            nSize += Diamond::getSize(mesh, 35);
        numDiamondsLineEdit->setText( QString::number(nSize) );
    */
}

///////////////////////////////////////////////////////////////////////////////
void JQuadDiamondsDialog :: incrementalRemove()
{
}

///////////////////////////////////////////////////////////////////////////////
void JQuadDiamondsDialog :: accept()
{
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDiamondsDialog :: reject()
{
    /*
         meshViewer->getDrawFace()->setColorMethod( nullptr );
         meshViewer->getDrawNode()->setColorMethod( nullptr );
    */
    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDiamondsDialog :: meshOpt()
{
    meshoptDialog->setViewManager( viewManager );
    meshoptDialog->show();
}

///////////////////////////////////////////////////////////////////////////////

void JQuadDiamondsDialog :: makeConnections()
{
    PushButton( colorPushButton,  [=] {setColor();});
    PushButton( removeAllPushButton, [=] {removeAll(); });
    PushButton( meshOptPushButton, [=] {meshOpt();});
    PushButton( incrementalRemovePushButton, [=] {incrementalRemove();});
    PushButton( searchPushButton, [=] {searchDiamonds();});
}

///////////////////////////////////////////////////////////////////////////////
