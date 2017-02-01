#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshSingularityGraphDialog.hpp"
#include "MeshViewer.hpp"
#include "NodeAttributesDialog.hpp"
#include "EdgeAttributesDialog.hpp"

class JMeshSingularityGraphDialog : public QDialog, public Ui::MeshSingularityGraphDialog {
    Q_OBJECT

public:
    JMeshSingularityGraphDialog( QWidget *parent = 0);
    ~JMeshSingularityGraphDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void openNodeAttributesDialog();
    void openEdgeAttributesDialog();

    void  mouseReleaseEvent( QMouseEvent *e);

    void getGraph();
    void displayMesh();

    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshEntityPickerPtr nodePicker;

    JMeshPtr mesh;
    JMeshSingularityGraph  mgraph;

    bool display_mesh_edges = 1;
    vector<JColor> edgeColor;
    boost::scoped_ptr<JMeshSingularityGraph>  singularGraph;
    boost::scoped_ptr<JNodeAttributesDialog> nodeAttribsDialog;
    boost::scoped_ptr<JEdgeAttributesDialog> edgeAttribsDialog;
    
    void init();
    void displayGraph();
    
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
