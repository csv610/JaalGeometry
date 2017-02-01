#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_BernHexOpsDialog.hpp"
#include "MeshViewer.hpp"
#include "BernHexOps.hpp"

using namespace Jaal;

class JBernHexOpsDialog : public QDialog, public Ui::BernHexOpsDialog {
    Q_OBJECT

public:
    JBernHexOpsDialog( QWidget *parent = 0);
    ~JBernHexOpsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void init();
    void patternSearch();
    void applyOp();
    void cancelOp();
    void demoElem1();
    void demoElem2();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh, dummymesh;
    JBernHexOps    bernOp;
    vector<JHexahedronPtr> localCells;

    void  demoOp(int t);
    void  demoOp1(int t);
    void  demoOp2(int t);
    void  demoOp3(int t);
    void  demoOp4(int t);
    void  demoOp5(int t);
    void  demoOp6(int t);

    void  searchOp1();
    void  searchOp2();
    void  searchOp3();
    void  searchOp4();
    void  searchOp5();
    void  searchOp6();

    void  applyOp1();
    void  applyOp2();
    void  applyOp3();
    void  applyOp4();
    void  applyOp5();
    void  applyOp6();

    bool isHexmesh();
    void makeConnections();
};
