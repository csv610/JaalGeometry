#pragma once

#include <QDialog>

#include "Ui_MeshMeanCurvatureFlowDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshMeanCurvatureFlow.hpp"
#include "DelaunayMesherDialog.hpp"

class JMeshMeanCurvatureFlowDialog : public QDialog, public Ui::MeshMeanCurvatureFlowDialog {
    Q_OBJECT

    struct ThreadWork
    {
        ThreadWork() {
            finished = -1;
        }
        int     finished;
        std::string  name;
        JMeshViewer *meshViewer;

        void  run();

        void  operator() ()
        {
            run();
        }
    };

public:
    JMeshMeanCurvatureFlowDialog( QWidget *parent = 0);
    ~JMeshMeanCurvatureFlowDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void  displayDelaunay();
    void  startAllOver();
    void  startFlow();
    void  projectOnSphere();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr  meshViewer;

private:
    JMeshPtr orgMesh, mesh;
    JMeshMeanCurvatureFlow  meanFlow;

    double xlen0, area0, area1;
    bool   rescale;
    bool   initialized = 0;
    int    currStep = 0;

    void init();
    void initMesh();
    void makeConnections();
};
