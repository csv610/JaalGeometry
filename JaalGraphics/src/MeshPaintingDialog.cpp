#include "MeshPaintingDialog.hpp"

void JMeshPainter :: actionMouseEvent( int id )
{
    if( viewManager == nullptr) return;

    if(picking) {
        const Point2I &pCurr  = viewManager->getMouseCurrentPixelPosition();
        JFaceSequence picked = pickingTexture->getPickedFaces( pCurr);
        JFaceRenderPtr  fattrib;
        for( JFacePtr face  : picked) {
            face->getAttribute("Render", fattrib);
            fattrib->color = color;
            faceSet.insert(face);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshPainter :: addEntities( const vector<size_t> &entities)
{
    if( mesh == nullptr ) return;

    if( entity_type == 0)  {
        JNodeRenderPtr  nattrib;
        for( size_t i : entities) {
            JNodePtr vtx = mesh->getNodeAt(i);
            vtx->getAttribute("Render", nattrib);
            nattrib->color = color;
            nodeSet.insert(vtx);
        }
    }

    if( entity_type == 1)  {
        JEdgeRenderPtr  eattrib;
        for( size_t i : entities) {
            JEdgePtr edge = mesh->getEdgeAt(i);
            edge->getAttribute("Render", eattrib);
            eattrib->color = color;
            edgeSet.insert(edge);
        }
    }

    if( entity_type == 2)  {
        JFaceRenderPtr  fattrib;
        for( size_t i : entities) {
            JFacePtr face = mesh->getFaceAt(i);
            face->getAttribute("Render", fattrib);
            fattrib->color = color;
            faceSet.insert(face);
        }
    }
    meshViewer->updateBuffers(mesh);
}
//////////////////////////////////////////////////////////////////////////////////////

void JMeshPainter :: drawRectBrush()
{
    viewManager->startScreenCoordinatesSystem();
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);

    glBlendFunc(GL_ONE, GL_ONE);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    const Point2I &pCurr  = viewManager->getMouseCurrentPixelPosition();

    int x0 =  (int)(pCurr[0] - 0.5*brushsize);
    int y0 =  (int)(pCurr[1] - 0.5*brushsize);

    int x1 =  (int)(pCurr[0] + 0.5*brushsize);
    int y1 =  (int)(pCurr[1] + 0.5*brushsize);

    glColor4f(0.0, 0.0, 0.3f, 0.3f);
    glBegin(GL_QUADS);
    glVertex2i( x0, y0 );
    glVertex2i( x1, y0 );
    glVertex2i( x1, y1 );
    glVertex2i( x0, y1 );
    glEnd();

    glLineWidth(2.0);
    glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
    glBegin(GL_LINE_LOOP);

    glVertex2i( x0, y0 );
    glVertex2i( x1, y0 );
    glVertex2i( x1, y1 );
    glVertex2i( x0, y1 );
    glEnd();

    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
    viewManager->stopScreenCoordinatesSystem();
}

//////////////////////////////////////////////////////////////////////////////////////

void JMeshPainter :: drawCircleBrush()
{
    viewManager->startScreenCoordinatesSystem();
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);

    glBlendFunc(GL_ONE, GL_ONE);

    const Point2I &pCurr  = viewManager->getMouseCurrentPixelPosition();
    int x0 =  pCurr[0];
    int y0 =  pCurr[1];
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glColor4f(0.0, 0.0, 0.3f, 0.3f);

    const unsigned int triangles = 50; // number of triangles
    const float twoPi = 2.0*M_PI;

    glBegin(GL_TRIANGLE_FAN);
    glVertex2i( x0, y0);
    double delta = twoPi / (double)triangles;
    for(unsigned int i = 0; i <= triangles; i++) {
        int x1 =  x0 + brushsize*cos(i*delta);
        int y1 =  y0 + brushsize*sin(i*delta);
        glVertex2i( x1, y1 );
    }
    glEnd();

    glLineWidth(2.0);
    glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
    glBegin(GL_LINE_LOOP);
    for(unsigned int i = 0; i <= triangles; i++) {
        int x1 =  x0 + brushsize*cos(i*delta);
        int y1 =  y0 + brushsize*sin(i*delta);
        glVertex2i( x1, y1 );
    }

    glEnd();

    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
    viewManager->stopScreenCoordinatesSystem();
}
//////////////////////////////////////////////////////////////////////////////////////

void JMeshPainter :: draw()
{
    if( viewManager == nullptr) return;

    if( !isActive() ) return;

    if( picking ) {
        if( brushtype == CIRCLE_BRUSH)
            drawCircleBrush();
        else
            drawRectBrush();
    }

    if( meshViewer == nullptr ) return;

    glDisable(GL_BLEND);
    glEnable( GL_DEPTH_TEST);

    if( entity_type == 0)  {
        JNodeDraw *nodeDraw = meshViewer->getNodeDraw();
        for( JNodePtr vtx : nodeSet)  {
            nodeDraw->draw(vtx);
        }
    }

    if( entity_type == 1)  {
        JEdgeDraw *edgeDraw = meshViewer->getEdgeDraw();
        for( JEdgePtr edge : edgeSet)  {
            edgeDraw->draw(edge);
        }
    }

    if( entity_type == 2)  {
        JFaceDraw *faceDraw = meshViewer->getFaceDraw();
        faceDraw->preRender();
        for( JFacePtr face : faceSet)  {
            faceDraw->draw(face);
        }
    }
}
///////////////////////////////////////////////////////////////////////////////

JMeshPaintingDialog :: JMeshPaintingDialog( QWidget *parent) : QDialog(parent)
{
    setupUi(this);
    makeConnections();
    painter.reset( new JMeshPainter);
}

///////////////////////////////////////////////////////////////////////////////

JMeshPaintingDialog :: ~JMeshPaintingDialog()
{
}

///////////////////////////////////////////////////////////////////////////////
void JMeshPaintingDialog :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    if( mesh == nullptr) return;
    JMeshRenderPtr mrender;
    mesh->getAttribute("Render", mrender);
    mrender->pickableEntity = 2;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshPaintingDialog :: setBrushColor()
{
    QColor color = QColorDialog::getColor();

    JColor rgb;
    rgb[0] = color.red()/255.0;
    rgb[1] = color.green()/255.0;
    rgb[2] = color.blue()/255.0;
    rgb[3] = 1.0;
    painter->setBrushColor(rgb);
}

///////////////////////////////////////////////////////////////////////////////
void JMeshPaintingDialog :: setBrushType()
{
    if( squareBrushRadioButton->isChecked() ) painter->setBrushType( JMeshPainter::QUAD_BRUSH);
    if( circleBrushRadioButton->isChecked() ) painter->setBrushType( JMeshPainter::CIRCLE_BRUSH);

}

///////////////////////////////////////////////////////////////////////////////
void JMeshPaintingDialog :: setBrushSize()
{
    QString str = brushSizeLineEdit->text() ;
    double val  = str.toDouble();
    painter->setBrushSize(val);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPaintingDialog :: init()
{
    if( viewManager == nullptr ) return;
    painter->setViewManager(viewManager);
    viewManager->attach(painter );

    painter->setViewManager(viewManager) ;
    freezeModel();
}

///////////////////////////////////////////////////////////////////////////////
void JMeshPaintingDialog :: beginSelection()
{
    if( viewManager == nullptr ) return;
    viewManager->deactivateComponents();
    viewManager->refreshDisplay();
    painter->setActive(1);
    painter->genTexture(mesh);
    viewManager->activateComponents();
}
///////////////////////////////////////////////////////////////////////////////

void JMeshPaintingDialog :: freezeModel()
{
    if( viewManager == nullptr ) return;

    bool freeze = freezemodelCheckBox->isChecked();
    viewManager->freezeView( freeze );

    if( freeze )  {
        this->beginSelection();
    }  else
        painter->setActive(0);
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPaintingDialog :: closeDialog()
{
    if( viewManager == nullptr ) return;

    viewManager->detach(painter );
    viewManager->freezeView( 0 );

    this->hide();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshPaintingDialog :: makeConnections()
{
    RadioButton( squareBrushRadioButton,  [=] {setBrushType();});
    RadioButton( circleBrushRadioButton,  [=] {setBrushType();});

    LineEdit(  brushSizeLineEdit,  [=] { setBrushSize();});
    CheckBox(  freezemodelCheckBox,  [=] {freezeModel();});

    PushButton( brushColorPushButton,    [=] {setBrushColor();});
    PushButton( closePushButton, [=] {closeDialog();});
}

///////////////////////////////////////////////////////////////////////////////
