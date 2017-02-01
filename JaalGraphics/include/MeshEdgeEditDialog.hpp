#pragma once

#include <QDialog>

#include "Ui_MeshEdgeEditDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"

class JMeshEdgeEditDialog : public QDialog, public Ui::MeshEdgeEditDialog {
    Q_OBJECT

public:
    JMeshEdgeEditDialog( QWidget *parent = 0);
    ~JMeshEdgeEditDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void getMinEdges();
    void getMaxEdges();
    void refine();
    void collapse();
    void deleteAll();
    void deleteInternal();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshEntityPickerPtr picker;

    JMeshPtr mesh;

    void setColor();

    void init();
    void setInfo();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
