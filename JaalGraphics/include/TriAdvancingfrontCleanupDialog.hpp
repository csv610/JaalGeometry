#pragma once

#include <QDialog>

#include "Ui_TriAdvancingfrontCleanupDialog.hpp"
#include "MeshViewer.hpp"

class JTriAdvancingfrontCleanupDialog : public QDialog, public Ui::TriAdvancingfrontCleanupDialog {
    Q_OBJECT

public:
    JTriAdvancingfrontCleanupDialog( QWidget *parent = 0);
    ~JTriAdvancingfrontCleanupDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;

    void init();
    void makeConnections();
};

