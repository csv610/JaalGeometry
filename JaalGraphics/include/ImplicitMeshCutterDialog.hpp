#pragma once

#include <QDialog>
#include "Ui_ImplicitMeshCutterDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"

class JImplicitMeshCutterDialog : public QDialog, public Ui::ImplicitMeshCutterDialog {
    Q_OBJECT

public:
    JImplicitMeshCutterDialog( QWidget *parent = 0);
    ~JImplicitMeshCutterDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void setPlane();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
