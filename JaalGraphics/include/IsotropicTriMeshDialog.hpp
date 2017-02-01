#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_IsotropicTriMeshDialog.hpp"
#include "MeshViewer.hpp"
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>
#include "IsotropicMesh.hpp"


////////////////////////////////////////////////////////////////////////////////

class JIsotropicTriMeshDialog : public QDialog, public Ui::IsotropicTriMeshDialog {
    Q_OBJECT

public:
    JIsotropicTriMeshDialog( QWidget *parent = 0);
    ~JIsotropicTriMeshDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:

    void genMesh();
    void setSharpEdges();
    void getEdgesLength();
    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

    JMeshPtr mesh;

    boost::scoped_ptr<QwtPlotGrid>  grid;
    boost::scoped_ptr<QwtPlotCurve> plotcurve;
    boost::scoped_ptr<QwtPlotMagnifier> magnifier;
    boost::scoped_ptr<QwtPlotPanner>    panner;

    vector<double>  xData, yData, qData;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
