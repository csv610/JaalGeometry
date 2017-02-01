#pragma once

#include <QDialog>
#include "Ui_MeshSegmentationDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshSegmentation.hpp"
#include "MeshGeodesicsDialog.hpp"
#include "MeshSDF.hpp"
#include "MeshPartitionDialog.hpp"    // For coloring the faces..

class JMeshSegmentationDialog : public QDialog, public Ui::MeshSegmentationDialog {
    Q_OBJECT

public:
    JMeshSegmentationDialog( QWidget *parent = 0);
    ~JMeshSegmentationDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void openMeshGeodesicsDialog();
    void openMeshPartitionDialog();

    void segment();
    void keyPressEvent( QKeyEvent *e);

    void getSDFSegments();
    void closeDialog();

private:
    JMeshPtr mesh;
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    boost::scoped_ptr<JMeshGeodesicsDialog>  meshGeodesicsDialog;
    boost::shared_ptr<JEdgeColor>  edgeColor;
    boost::shared_ptr<JFaceColor>  faceColor;
    boost::shared_ptr<JMeshPartitionDialog>  meshPartDialog;

    void init();
    void makeConnections();
};
