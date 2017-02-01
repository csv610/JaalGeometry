#pragma once

#include <QDialog>

#include "Ui_ObjectsListDialog.hpp"
#include "MeshViewer.hpp"
#include "ImageViewer.hpp"

class JObjectsListDialog : public QDialog, public Ui::ObjectsListDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    static const int  ALL_OBJECTS   = 0;
    static const int  MESH_OBJECTS  = 1;
    static const int  IMAGE_OBJECTS = 2;

    JObjectsListDialog( QWidget *parent = 0);
    ~JObjectsListDialog();

    void setViewManager( JaalViewer *v)
    {
        viewManager = v;
        init();
    }

    void setType(int t = 0);

    JMeshPtr   getSelectedMesh() const {
        return currMesh;
    }
    JImagePtr  getSelectedImage() const {
        return currImage;
    }

private slots:
    void getCell(int i, int j);
    void deleteObject();
    void saveObject();
    void closeDialog();
    void renameObject();

private:
    JMeshViewerPtr  meshViewer;
    JImageViewerPtr imageViewer;

    bool   listObjects[2];

    JMeshPtr   currMesh;
    JImagePtr  currImage;

    void init();
    void makeConnections();
    void filltable();
    void saveMesh();
};
