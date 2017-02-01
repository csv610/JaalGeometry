#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_ShapeOpDialog.hpp"
#include "MeshViewer.hpp"
#include "ShapeOpt.hpp"

class JShapeOpDialog : public QDialog, public Ui::ShapeOpDialog {
    Q_OBJECT

public:
    JShapeOpDialog( QWidget *parent = 0);
    ~JShapeOpDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

    void setMesh( const JMeshPtr &m);

private slots:
    void applyAreaConstraints();
    void applyCocircularConstraints();
    void applyCoplanarConstraints();
    void applyCrossFieldConstraints();
    void applyEdgeStrainConstraints();
    void applyFixedLengthConstraints();
    void applyLaplaceConstraints();
    void applyParallelogramConstraints();
    void applyTriangleStrainConstraints();
    void applyRectangleConstraints();
    void applyClosenessConstraints();
    void applyFixedBoundaryConstraints();
    void applySimilarityConstraints();
    void applyRigidConstraints();

    void displayConstraints();
    void solve();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::scoped_ptr<JShapeOptimizer>  shapeOpt;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
