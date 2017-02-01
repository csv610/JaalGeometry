#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_QuadMeshCleanupDialog.hpp"
#include "MeshViewer.hpp"
#include "Singlet.hpp"
#include "Doublet.hpp"
#include "Diamond.hpp"

#include "QuadMeshDualsDialog.hpp"

class JQuadMeshCleanupDialog : public QDialog, public Ui::QuadMeshCleanupDialog {
    Q_OBJECT

public:
    JQuadMeshCleanupDialog( QWidget *parent = 0);
    ~JQuadMeshCleanupDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void searchSinglets();
    void searchDoublets();
    void searchDiamonds();

    void removeSinglets();
    void removeDoublets();
    void removeDiamonds();
    void openDualDialog();
    void swapEdges();
    void closeDialog();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr  meshViewer;

private:
    JMeshPtr mesh;
    boost::scoped_ptr<JSinglet> qSinglet;
    boost::scoped_ptr<JDoublet> qDoublet;
    boost::scoped_ptr<JDiamond> qDiamond;

    boost::scoped_ptr<JQuadMeshDualsDialog> dualDialog;

    void init();
    double creaseAngle;
    void makeConnections();
};
