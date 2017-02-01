#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_TangleFEMTestsDialog.hpp"
#include "MeshViewer.hpp"
#include "Elasticity2D.hpp"
#include "AllQuadMeshGenerator.hpp"
#include "AllTriMeshGenerator.hpp"
#include "MeshRefine.hpp"

class JTangleFEMTestsDialog : public QDialog, public Ui::TangleFEMTestsDialog {
    Q_OBJECT

    JaalViewer *viewManager;
public:
    JTangleFEMTestsDialog( QWidget *parent = 0);
    ~JTangleFEMTestsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }

private slots:
    void setParams();
    void solve();
    void setMesh();
    void setShapeOrder();
    void displayQuadraticNodes();

private:
    JMeshViewerPtr meshViewer;

    int randomTangle, checkTangle, meshType;
    int elemOrder, numGuassPnts, triGaussPnts;
    int  numGauss1;            // (1,3,7)
    int  numGauss2;            // (1,3,7)
    bool absJacobian;
    int  fieldType;
    int  shapeFamily;
    int  checkConvex;
    int  problemID;
    double youngModulus;
    double poissonRatio;

    JMeshPtr mesh;
//     JElasticity2D  elastic;

    void patchTest();

    void genRandomMesh();
    void genSingleNodeMesh();
    void genSimpleTangleMesh();
    void genRotateQuadMesh();
    void genBendingBeamMesh();
    void genSquareCircle1Mesh();
    void genSquareCircle2Mesh();

    void solveSingleNode();
    void solveSimpleTangle();
    void solveRandomTangle();
    void solveRotateQuads();
    void solveLinearElasticity();

    void analyticalSol( const Point3D &p, double &u, double &v);

    void init();
    void makeConnections();
};

