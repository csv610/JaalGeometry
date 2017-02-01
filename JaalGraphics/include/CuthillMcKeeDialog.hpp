#pragma once

#include <QDialog>

#include "Ui_CuthillMcKeeDialog.hpp"
#include "MeshViewer.hpp"

class JCuthillMcKeeDialog : public QDialog, public Ui::CuthillMcKeeDialog {
    Q_OBJECT

public:
    JCuthillMcKeeDialog( QWidget *parent = 0);
    ~JCuthillMcKeeDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void applyAlgorithm();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
