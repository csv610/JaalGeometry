#pragma once

#include <QDialog>
#include <QColorDialog>
#include <QSlider>

#include <igl/parula.h>
#include <igl/jet.h>

#include "Ui_MeshFeaturesDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshCurvature.hpp"
#include "EdgeAttributesDialog.hpp"


#include "MeshSpectrum.hpp"

class JMeshFeaturesDialog : public QDialog, public Ui::MeshFeaturesDialog {
    Q_OBJECT

public:
    JMeshFeaturesDialog( QWidget *parent = 0);

    ~JMeshFeaturesDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void angleDefect();
    void sharpEdges();
    void setGaussianCurvature();
    void setMeanCurvature();
    void getEigenVectors();
    void displayEigenVector();
    void getCurvatureDirections();
    void openMinKAttribDialog();
    void openMaxKAttribDialog();
    void setNodesColor();

private:
    JMeshPtr mesh;
    JaalViewer  *viewManager;
    JMeshViewerPtr meshViewer;
    boost::scoped_ptr<JMeshSpectrum> meshSpectrum;
    boost::scoped_ptr<JEdgeAttributesDialog>  edgeAttribsDialog;

    void init();
    void setNodesColor( int id);
    void makeConnections();
    void setEdgeAttributes( const JMeshPtr &m);
};
