#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshSkeletonShapesDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshSkeleton.hpp"

class JMeshSkeletonShapesDialog : public QDialog, public Ui::MeshSkeletonShapesDialog {
    Q_OBJECT

public:
    JMeshSkeletonShapesDialog( QWidget *parent = 0);
    ~JMeshSkeletonShapesDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);
    void setSkeleton( const JMeshSkeletonPtr &s) {
         meshSkel = s;
    }

private slots:
    void getContours();
    void closeDialog();
    
private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JMeshSkeletonPtr meshSkel;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
