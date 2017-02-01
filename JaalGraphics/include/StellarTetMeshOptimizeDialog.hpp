#pragma once

#include <QDialog>

#include "Ui_StellarTetMeshOptDialog.hpp"
#include "MeshViewer.hpp"


class JStellarTetMeshOptimizeDialog : public QDialog, public Ui::StellarTetMeshOptDialog {
    Q_OBJECT

public:
    JStellarTetMeshOptimizeDialog( QWidget *parent = 0);
    ~JStellarTetMeshOptimizeDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void optmesh();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;

    void init();
    void genConfigFile();
    void makeConnections();

};
