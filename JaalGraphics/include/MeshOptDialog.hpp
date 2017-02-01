#pragma once

#include "Ui_MeshOptDialog.hpp"

#include <QDialog>
#include <QColorDialog>
#include "MeshViewer.hpp"
#include "MeshGeometricQualityDialog.hpp"
#include "MeshLaplaceSmoothingDialog.hpp"
#include "MeshMeanCurvatureFlowDialog.hpp"
#include "LloydRelaxationDialog.hpp"
#include "MeshUntangleDialog.hpp"
#include "MeshConstraintsDialog.hpp"
#include "ShapeOpDialog.hpp"

class JMeshLaplaceSmoothingDialog;
class JMeshConstraintsDialog;

class JMeshOptDialog : public QDialog, public Ui::MeshOptimizationDialog {
    Q_OBJECT

public:
    JMeshOptDialog( QWidget *parent = 0);
    ~JMeshOptDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

protected:
    virtual void  keyPressEvent( QKeyEvent *e);
    virtual void  showEvent( QShowEvent *e);

private slots:

    void getOriginal();
    void nonlinear();
    void setConstraints();
    void closeDialog();
    void untangle();
    void reparamCurves();
    void smoothCurves();

    void openLaplaceDialog();
    void openLloydDialog();
    void openMeanCurvatureFlowDialog();
    void openUntangleDialog();
    void openShapeOpDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;
    boost::scoped_ptr<JMeshQuality> meshQual;
    boost::scoped_ptr<JMeshGeometricQualityDialog> geomQualityDialog;
    boost::scoped_ptr<JMeshLaplaceSmoothingDialog> laplaceDialog;
    boost::scoped_ptr<JMeshConstraintsDialog> constraintsDialog;
    boost::scoped_ptr<JLloydRelaxationDialog> lloydRelaxationDialog;
    boost::scoped_ptr<JMeshMeanCurvatureFlowDialog> meanCurvatureFlowDialog;
    boost::scoped_ptr<JMeshUntangleDialog> untangleDialog;
    boost::scoped_ptr<JShapeOpDialog> shapeOpDialog;

    vector<double> orgCoords;
    vector<size_t> l2g;

    JFaceQualityColor *faceColor;
    JCellQualityColor *cellColor;

    JMeshNonlinearOptimization nonlinearOpt;

    void init();
    void setMesh( const JMeshPtr &m);
    void makeConnections();
    void optimize();
};
