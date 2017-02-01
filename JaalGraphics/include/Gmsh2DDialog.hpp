#pragma once

#include <QDialog>

#include "Ui_Gmsh2DDialog.hpp"
#include "MeshViewer.hpp"
#include "MSTQuadMesher.hpp"

class JGmsh2DDialog : public QDialog, public Ui::Gmsh2DDialog {
    Q_OBJECT

public:
    JGmsh2DDialog( QWidget *parent = 0);
    ~JGmsh2DDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void generate();
    void reparamBoundary();
    void smoothBoundary();
    void rejectMesh();
    void closeDialog();
protected:
    void showEvent ( QShowEvent * event );

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    int     algorithm = 2;
    int     smoothCrossField = 0;
    int     crossFieldSmoothSteps = 0;
    double  minLen = 0;
    double  maxLen = 1E+22;
    int     sizeFromPoints    = 1;
    int     sizeFromCurvature = 0;
    int     recombinationAlgorithm = 0;
    int     recombineAll = 0;
    int     smoothSteps = 1;
    int     lloydSteps  = 0;
    int     subDivision = 0;
    double  maxAnisotrophy = 1.0E+33;
    double  minElementSize = 0.0;
    double  maxElementSize = 1.0E+22;
    double  lengthFactor   = 1.0;
    bool    halfSamples    = 0;
    bool    modifyBoundary = 0;

    JMeshPtr mesh;
    JMeshPtr newQuadMesh;

    void init();
    void makeConnections();
    void getOptions();
    void recoverBoundary();
};
////////////////////////////////////////////////////////////////////////////////
