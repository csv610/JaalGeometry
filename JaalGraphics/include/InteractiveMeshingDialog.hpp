#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_InteractiveMeshingDialog.hpp"
#include "MeshViewer.hpp"

class JInteractiveMeshingDialog : public QDialog, public Ui::InteractiveMeshingDialog {
    Q_OBJECT

public:
    JInteractiveMeshingDialog( QWidget *parent = 0);
    ~JInteractiveMeshingDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void setNode();
    void setEntity();
    void keyPressEvent( QKeyEvent *e);

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr oldmesh, mesh;
    JNodeSequence nodeseq;
    vector<int>  nodelist;

    void init();
    void genNewEdge();
    void genNewFace();
    void genNewCell();

    void makeConnections();
};
