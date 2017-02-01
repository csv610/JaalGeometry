#pragma once

#include <QDialog>
#include "Ui_SphereHexMesherDialog.hpp"

#include "MeshViewer.hpp"
#include "MeshMeanCurvatureFlowDialog.hpp"
#include "MeshLaplaceSmoothingDialog.hpp"
#include "SphereHexMesher.hpp"
#include "ObjectsListDialog.hpp"

class JSphereHexMesherDialog : public QDialog, public Ui::SphereHexMesherDialog {
    Q_OBJECT

public:
    JSphereHexMesherDialog( QWidget *parent = 0);
    ~JSphereHexMesherDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:
    void showEvent(QShowEvent *e);

private slots:
    void openCurvatureFlowDialog();
    void openLaplaceDialog();
    void openHarmonicDialog();
    void loadInputMesh();
    void loadSphereMesh();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    int selectMesh[3];
    int eulerChar;
    int topDim;

    JMeshPtr mesh, sphMesh;

    boost::scoped_ptr<JObjectsListDialog> objectsListDialog;
    boost::scoped_ptr<JMeshMeanCurvatureFlowDialog> meanCurvatureFlowDialog;
    boost::scoped_ptr<JMeshLaplaceSmoothingDialog>  meshLaplaceSmoothingDialog;
    boost::scoped_ptr<JSphereHexMesher>             sphHexMesher;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
