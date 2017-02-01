#pragma once

#include <QDialog>

#include "Ui_EditPartitionRegionDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshPartitioner.hpp"

class JEditPartitionRegionDialog : public QDialog, public Ui::EditPartitionRegionDialog {
    Q_OBJECT

public:
    JEditPartitionRegionDialog( QWidget *parent = 0);
    ~JEditPartitionRegionDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:

private:
    JaalViewer *viewManager;
    JMeshViewer *meshViewer;

    QStandardItemModel *model;

    JMeshPtr mesh;

    JMeshPartitioner mp;

    void init();
    void makeConnections();
    void filltable();
};
