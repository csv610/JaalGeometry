#pragma once

#include <QDialog>

#include "Ui_LSystemDialog.hpp"
#include "MeshViewer.hpp"
#include "AllTriMeshGenerator.hpp"
#include "AllQuadMeshGenerator.hpp"

using namespace Jaal;

class LSystemDialog : public QDialog, public Ui::LSystemDialog {

    Q_OBJECT

public:
    static const int SIERPINSKI_CARPET    = 0;
    static const int SIERPINSKI_PYRAMID   = 1;
    static const int SIERPINSKI_TRIANGLES = 2;
    static const int MENGER_SPONGE        = 3;

    LSystemDialog( QWidget *parent = 0);

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void keyPressEvent( QKeyEvent *e);

    void genShape();

    void applyDialog();
    void cancelDialog();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    void init();
    void makeConnections();
};

