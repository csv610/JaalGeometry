#pragma once

#include <QDialog>

#include "Ui_MeshLaplaceSmoothingDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshLaplacian.hpp"
#include "MeshDualGrapherDialog.hpp"
#include "MeshConstraintsDialog.hpp"

class JMeshConstraintsDialog;

class JMeshLaplaceSmoothingDialog : public QDialog, public Ui::MeshLaplaceSmoothingDialog
{

    Q_OBJECT

public:

    JMeshLaplaceSmoothingDialog( QWidget *parent = 0);
    ~JMeshLaplaceSmoothingDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh(const JMeshPtr &m);

protected:
    void  showEvent(QShowEvent *e);

private slots:

    void smooth();
    void setEdgesWeight();
    void getOriginal();
    void closeDialog();
    void openConstraintsDialog();
    void setConstraints();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr  meshViewer;

    bool   initialized = 0;
    JMeshPtr mesh;
    vector<double> orgCoords;
    vector<size_t> l2g;
    bool edges_weight_assigned;
    boost::scoped_ptr<JLaplaceMeshSmoother> laplace;
    boost::shared_ptr<JMeshConstraintsDialog> meshConstraintsDialog;

    void init();
    void initMesh();
    void assignColors();
    void makeConnections();
    int  setPrimalConstraints();
    void freeConstraints();
    void setConstraints( JMeshPtr submesh, size_t &nFixed, size_t &nFree);
    void displayNegativeElements();
};

