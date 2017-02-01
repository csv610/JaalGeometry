#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_RingQuadsDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshSkeletonDialog.hpp"

class JRingQuadsDialog : public QDialog, public Ui::RingQuadsDialog {
    Q_OBJECT

public:
    JRingQuadsDialog( QWidget *parent = 0);
    ~JRingQuadsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void openMeshSkeletonDialog();
    void getDisks();
    void getQuadMesh();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::scoped_ptr<JMeshSkeletonDialog> skeletonDialog;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
