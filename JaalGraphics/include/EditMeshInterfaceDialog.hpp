#pragma once

#include <QDialog>

#include "Ui_EditMeshInterfaceDialog.hpp"
#include "MeshViewer.hpp"

class JEditMeshInterfaceDialog : public QDialog, public Ui::EditMeshInterfaceDialog {
    Q_OBJECT

public:
    JEditMeshInterfaceDialog( QWidget *parent = 0);
    ~JEditMeshInterfaceDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

private slots:

private:
    JaalViewer *viewManager;
    JMeshViewer *meshViewer;
    JMeshPtr mesh;

    QStandardItemModel *model;

    JMeshPartitioner mp;

    struct Item
    {
        bool operator < ( const Item &rhs) const {
            return numSegments < rhs.numSegments;
        }
        int          id;
        int          numSegments;
        double       length;
    };

    void init();
    void makeConnections();
    void filltable();

};
