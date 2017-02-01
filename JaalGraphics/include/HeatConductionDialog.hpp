#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_HeatConductionDialog.hpp"
#include "MeshViewer.hpp"
#include "HeatConduction3D.hpp"


class JHeatConductionDialog : public QDialog, public Ui::HeatConductionDialog {
    Q_OBJECT
    JaalViewer *viewManager;
public:
    JHeatConductionDialog( QWidget *parent = 0);
    ~JHeatConductionDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void solve();
    void applyBoundConditions();
    void newGroup();
    void setLinearSolver();

private:
    JMeshPtr mesh;
    JMeshViewerPtr meshViewer;
    JMeshEntityPickerPtr picker;
    HeatConduction3D  heatSolver;
    void init();
    void makeConnections();
};
