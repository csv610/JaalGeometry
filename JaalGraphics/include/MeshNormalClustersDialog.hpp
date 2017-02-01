#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshNormalClustersDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshNormalClusters.hpp"
#include "MeshPartitionColors.hpp"

////////////////////////////////////////////////////////////////////////////////

class JMeshNormalClustersDialog : public QDialog, public Ui::MeshNormalClustersDialog {
    Q_OBJECT

public:
    JMeshNormalClustersDialog( QWidget *parent = 0);
    ~JMeshNormalClustersDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void getClusters();
    void openNormalsDialog();
    void showSphere();
    void closeDialog();
    
private:
    JaalViewer *viewManager = nullptr;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    JMeshPtr sphMesh;

    boost::scoped_ptr<JMeshNormalClusters>  normalClusters;
    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
