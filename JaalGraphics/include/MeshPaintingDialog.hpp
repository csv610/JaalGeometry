#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshPaintingDialog.hpp"
#include "MeshViewer.hpp"
#include "PickingTexture.hpp"

///////////////////////////////////////////////////////////////////////////////

class JMeshPainter : public JViewComponent
{
public:
    static const int QUAD_BRUSH   = 1;
    static const int CIRCLE_BRUSH = 0;

    JMeshPainter() {
        brushtype = CIRCLE_BRUSH;
        brushsize = 50;
        name      = "MeshPainter";
        picking    = 0;
        entity_type = 2;
        pickingTexture.reset( new JPickingTexture);
    }
    ~JMeshPainter() {}

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        pickingTexture->setViewManager(v);
        JViewComponentPtr c = viewManager->getComponent("PrimalMeshViewer");
        meshViewer = dynamic_pointer_cast<JMeshViewer>(c);
    }

    void setEntityType( int e ) {
        entity_type = e;
    }
    void setBrushColor( const JColor &c) {
        color = c;
    }
    void setBrushType( int t) {
        brushtype = t;
    }
    void setBrushSize( double l) {
        brushsize = l;
    }

    void genTexture( const JMeshPtr &m) {
        mesh = m;
        pickingTexture->genTexture(mesh);
    }

    void setPicking( bool v ) {
        picking = v;
    }
    void addEntities( const vector<size_t> &e);

    void actionMouseEvent(int id);

    void draw();
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    int     brushtype;
    double  brushsize;
    JColor  color;
    int     entity_type;
    int     picking;
    boost::scoped_ptr<JPickingTexture> pickingTexture;

    JNodeSet  nodeSet;
    JEdgeSet  edgeSet;
    JFaceSet  faceSet;

    void drawCircleBrush();
    void drawRectBrush();
};

///////////////////////////////////////////////////////////////////////////////

class JMeshPaintingDialog : public QDialog, public Ui::MeshPaintingDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    JMeshPaintingDialog( QWidget *parent = 0);
    ~JMeshPaintingDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

    void beginSelection();

private slots:

    void setBrushColor();
    void setBrushSize();
    void setBrushType();
    void freezeModel();

    void closeDialog();

private:
    JMeshPtr    mesh;
    boost::shared_ptr<JMeshPainter> painter;

    void init();
    void makeConnections();

};

