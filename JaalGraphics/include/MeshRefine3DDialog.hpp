#pragma once

#include <QDialog>

#include "Ui_MeshRefine3DDialog.hpp"

#include "MeshViewer.hpp"
#include "MeshRefine.hpp"
#include "BernHexOps.hpp"

class JMeshRefine3DDialog : public QDialog, public Ui::MeshRefine3DDialog {
    Q_OBJECT

public:
    JMeshRefine3DDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void insertPillows();
    void refineHex17();
    void refineHex18();
    void genHexBlocks();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;

    JMeshPtr hexmesh, tetmesh;

    void init();
    void makeConnections();
};
