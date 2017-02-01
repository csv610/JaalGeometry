#pragma once

#include <QDialog>

#include "Ui_SurfaceReconstructionDialog.hpp"
#include "MeshViewer.hpp"
#include "SurfaceReconstruction.hpp"

class JSurfaceReconstructionDialog : public QDialog, public Ui::SurfaceReconstructionDialog {
    Q_OBJECT

public:
    JSurfaceReconstructionDialog( QWidget *parent = 0);
    ~JSurfaceReconstructionDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void genMesh();
    void rejectMesh();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;
    JMeshPtr mesh;
    JMeshPtr newMesh;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
