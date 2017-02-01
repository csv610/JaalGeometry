#pragma once

#include <QDialog>

#include "Ui_MeshEdgesDialog.hpp"
#include "MeshViewer.hpp"
#include "EdgeAttributesDialog.hpp"
#include "MeshEdgeEditDialog.hpp"

class MeshEntityAttributesDialog;

class JMeshEdgesDialog : public QDialog, public Ui::MeshEdgesDialog {
    Q_OBJECT

public:
    static const int  WIREFRAME         = 0;
    static const int  BACKFACE_CULLING  = 1;
    static const int  FRONTLINES        = 2;
    static const int  DEPTH_TESTING     = 3;

    JMeshEdgesDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

    void setMesh(const JMeshPtr &m);

private slots:
    void showEvent(QShowEvent *e);
    void keyPressEvent( QKeyEvent *e);

    void countBoundEdges();
    void getHiddenlines();

    void setInternal();
    void setBoundary();
    void getNonManifoldEdges();

    void closeDialog();
    void setNumVisible();
    void checkState();
    void checkDisplay();

    void lookAt();

    void saveBoundary();
     

    /*
         void checkLower();
         void checkEdges();
         void setTransparency();
         void showOne();
         void getEdgeID();
         void nextEdge();
         void showNeighs();
         void alignAlongAxis( bool refresh = 1 );
         void changeCenter( bool refresh = 1 );
    */

    void setDefaultColorMethod();
    void openAttribListDialog();
    void openEditDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    JEdgeDraw  *drawEdge;

    boost::scoped_ptr<JEdgeColor> edgeColor;
    boost::scoped_ptr<JEdgeAttributesDialog>  attribDialog;
    boost::scoped_ptr<JMeshEntityAttribListDialog> attribListDialog;
    boost::scoped_ptr<JMeshEdgeEditDialog> editDialog;

    int  renderStyle;
    JEdgePtr currEdge;

    void reset_neighs( bool val);
    void activate(JEdgePtr edge );

    void init();
    void delAttrib();
    void makeConnections();
};

