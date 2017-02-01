#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshSpacePartitionsDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"

class JMeshSpacePartitionsDialog : public QDialog, public Ui::MeshSpacePartitionsDialog {
    Q_OBJECT

public:
    JMeshSpacePartitionsDialog( QWidget *parent = 0);
    ~JMeshSpacePartitionsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void genPartitions();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
