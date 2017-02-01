#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_LegoBuilderDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshAffineTransforms.hpp"

class JLegoBuilderDialog : public QDialog, public Ui::LegoBuilderDialog {
    Q_OBJECT

public:
    JLegoBuilderDialog( QWidget *parent = 0);
    ~JLegoBuilderDialog();

    void setViewManager( JaalViewer *v) {
        viewManager  = v;
        init();
    }

private slots:
    void newCubeIJK();
    void newCubeNodes();
    void attachAtFace();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    int boxDim[3];
    JMeshPtr oldmesh, newmesh;

    void init();
    void makeConnections();
    int  getCellID( int i, int j, int k);
};
