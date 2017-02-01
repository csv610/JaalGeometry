#pragma once

#include <QDialog>

#include "Ui_MeshTopologyQueryDialog.hpp"
#include "MeshViewer.hpp"
#include "MatrixViewer.hpp"
#include "MeshComponentsDialog.hpp"

class JMeshTopologyQueryDialog : public QDialog, public Ui::MeshTopologyQueryDialog {
    Q_OBJECT
public:
    JMeshTopologyQueryDialog( QWidget *parent = 0);

    ~JMeshTopologyQueryDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

protected:
    void showEvent( QShowEvent *e);

private slots:
    void getConsistent();
    void displayMatrix();
    void displayPrimalMatrix();
    void displayDualMatrix();
    void closeDialog();
    void setLapGrid();
    void setLapFonts();
    void setStepLabel();
    void setPointSize();
    void getBettiNumber();
    void openMeshComponentsDialog();
    void getOrphaned();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;
    JMatrixViewerPtr matrixViewer;

    bool  initialized;
    JMeshPtr mesh;
    boost::scoped_ptr<JNodeDegreeColor> nodeColor;
    boost::scoped_ptr<JMeshComponentsDialog>  meshComponentsDialog;

    void getInfo();
    void init();
    void makeConnections();
};
