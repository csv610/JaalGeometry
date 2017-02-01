#pragma once

#include <QDialog>

#include "Ui_MeshVoxelizerDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshVoxelizer.hpp"
#include "ObjectsListDialog.hpp"

class JMeshVoxelizerDialog : public QDialog, public Ui::MeshVoxelizerDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    JMeshVoxelizerDialog( QWidget *parent = 0);
    ~JMeshVoxelizerDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

private slots:
    void genMesh();
    void slicer();
    void showEvent(QShowEvent *e);
    void openMeshlistDialog();
    void saveAs();
    void readFile();
    void monotoneVoxels();

private:
    JMeshViewerPtr meshViewer;
    vector<int> cellDim;

    JMeshPtr mesh;          // The original model which we need to voxelized..
    JMeshPtr hexmesh;       // The voxelized cubical mesh ...
    JVoxelMeshPtr voxmesh;  // The voxel mesh. It contains both the model mesh and
    // the background mesh ...

    boost::scoped_ptr<JObjectsListDialog> meshlistDialog;

    void init();
    void makeConnections();
    void viewCell( const JCellPtr &cell, bool v);
};

