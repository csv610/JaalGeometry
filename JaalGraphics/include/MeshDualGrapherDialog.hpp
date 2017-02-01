#pragma once

#include <QDialog>

#include "Ui_MeshDualGrapherDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshDualGraph.hpp"
#include "MeshRenderDialog.hpp"

class JMeshDualGrapherDialog : public QDialog, public Ui::MeshDualGrapherDialog
{
    Q_OBJECT

public:
    JMeshDualGrapherDialog( QWidget *parent = 0);
    ~JMeshDualGrapherDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void  newGraph();
    void  openMeshRenderDialog();
    void  closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JMeshPtr dGraph;
    boost::scoped_ptr<JMeshRenderDialog> meshRenderDialog;

    void init();
    void makeConnections();
};

