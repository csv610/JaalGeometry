#pragma once

#include <QDialog>

#include "Ui_MeshPartitionDialog.hpp"

#include "MeshFacesPartitionsDialog.hpp"
#include "MeshEdgesPartitionsDialog.hpp"

#include "MeshViewer.hpp"
#include "GraphColor.hpp"
#include "DrawMeshEntity.hpp"
#include "MeshNormalClustersDialog.hpp"
#include "MeshPartitionColors.hpp"
#include "MeshFacesClustering.hpp"


class JMeshPartitionDialog : public QDialog, public Ui::MeshPartitionDialog {
    Q_OBJECT

public:
    JMeshPartitionDialog( QWidget *parent = 0);

    ~JMeshPartitionDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void optPartition();

    void applyAlgorithm();
    void openCellsPartitionsDialog();
    void openFacesPartitionsDialog();
    void openEdgesPartitionsDialog();
    void openNodesPartitionsDialog();
    void openNormalClustersDialog();
    void getRegionGrowingClusters();
    void removeZigZagInterfaces();

    void savePartitions();

    void clearAll();

    void getTopologicalDisks();
    void keyPressEvent( QKeyEvent *e);
    void closeDialog();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;

//  boost::scoped_ptr<JMeshNodesPartitionsDialog> nodesPartitionsDialog;
    boost::scoped_ptr<JMeshFacesPartitionsDialog> facesPartitionsDialog;
    boost::scoped_ptr<JMeshEdgesPartitionsDialog> edgesPartitionsDialog;
    boost::scoped_ptr<JMeshNormalClustersDialog> normalClustersDialog;

    JNodeColorPtr nodeColor;
    JEdgeColorPtr edgeColor;
    JFaceColorPtr faceColor;

    void init();
    void makeConnections();

    void metisPartition();
    void convexPartition();
    void assignColors();

};

