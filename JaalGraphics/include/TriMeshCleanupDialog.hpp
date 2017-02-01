#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_TriMeshCleanupDialog.hpp"
#include "MeshViewer.hpp"
#include "DelaunayMesh.hpp"
#include "TriDecimator.hpp"
#include "SwapEdges.hpp"
#include "TriAdvancingfrontCleanupDialog.hpp"
#include "LloydRelaxationDialog.hpp"

class JTrimeshCleanupDialog : public QDialog, public Ui::TriMeshCleanupDialog {
    Q_OBJECT

public:
    JTrimeshCleanupDialog( QWidget *parent = 0);
    ~JTrimeshCleanupDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void makeDelaunay();
    void displayDelaunay();
    void displayNonManifoldEdges();
    void displayIrregularNodes();
    void displayCircumCircles();
    void displayCollapsables();
    void displayFlippables();
    void openLloydDialog();

    void flipEdges();
    void subdivideEdges();
    void removeEdges();

    void below5DegreeNodes();
    void above7DegreeNodes();
    void openAdvancingfront();

    void closeDialog();

private:
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    ManifoldEdgeColor *manifoldEdgeColor;
    JMeshEntityPickerPtr entityPicker;

    boost::scoped_ptr<JLloydRelaxationDialog> lloydRelaxationDialog;
    boost::scoped_ptr<JTriAdvancingfrontCleanupDialog>  advfrontDialog;

    void init();
    double creaseAngle;
    void makeConnections();
};
