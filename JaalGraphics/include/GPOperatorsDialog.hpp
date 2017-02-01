#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_GPOperatorsDialog.hpp"
#include "MeshViewer.hpp"
#include "QuadCleanUp.hpp"

class JGPOperatorsDialog : public QDialog, public Ui::GPOperatorsDialog {
    Q_OBJECT

public:
    JGPOperatorsDialog( QWidget *parent = 0);
    ~JGPOperatorsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:

    void init();

    void getSeeds();
    void shiftLeftOp();
    void shiftRightOp();
    void collapseOp();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    JEdgePtr currEdge;
    boost::shared_ptr<JQuadDual> qdual;

    void makeConnections();

};

