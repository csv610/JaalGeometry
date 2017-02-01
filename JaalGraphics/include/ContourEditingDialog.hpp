#pragma once

#include <QDialog>

#include "Ui_ContourEditingDialog.hpp"
#include "MeshViewer.hpp"
#include "PolygonSimplifyDialog.hpp"

class JContourEditingDialog : public QDialog, public Ui::ContourEditingDialog {
    Q_OBJECT

public:
    JContourEditingDialog( QWidget *parent = 0);
    ~JContourEditingDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void reparameterize();
    void curveShortening();
    void smoothing();

    void openPolySimplifyDialog();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    vector<JEdgeSequence> contours;
    std::unique_ptr<JPolygonSimplifyDialog> polySimplifyDialog;

    void init();
    void makeConnections();
    void extractContours();
};
////////////////////////////////////////////////////////////////////////////////
