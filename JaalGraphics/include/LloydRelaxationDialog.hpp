#pragma once

#include <QDialog>

#include "Ui_LloydRelaxationDialog.hpp"
#include "LloydRelaxationDialog.hpp"
#include "MeshViewer.hpp"

class JLloydRelaxationDialog : public QDialog, public Ui::LloydRelaxationDialog {
    Q_OBJECT

public:
    JLloydRelaxationDialog( QWidget *parent = 0);
    ~JLloydRelaxationDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void smoothAll();
    void getOriginal();
    void closeDialog();
    void checkDisplay();

private:
    JaalViewer *viewManager;

    JMeshViewerPtr meshViewer;
    boost::shared_ptr<JMeshViewer> dualViewer;

    vector<double> orgCoords;
    vector<size_t> l2g;

    JMeshPtr mesh, dualGraph;
    bool initialized = 0;
    JLloydMeshOptimizer mopt;

    void init();
    void initMesh();
    void makeConnections();
};

