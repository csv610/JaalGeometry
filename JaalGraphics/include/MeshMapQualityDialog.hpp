#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshMapQualityDialog.hpp"
#include "MeshViewer.hpp"

#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>

class JMeshMapQualityDialog : public QDialog, public Ui::MeshMapQualityDialog {
    Q_OBJECT

public:
    JMeshMapQualityDialog( QWidget *parent = 0);
    ~JMeshMapQualityDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setSourceMesh( const JMeshPtr &m);
    void setTargetMesh( const JMeshPtr &m);

private slots:
    void checkQuality();

private:
    JaalViewer *viewManager;

private:
    boost::scoped_ptr<QwtPlotGrid>       grid;
    boost::scoped_ptr<QwtPlotCurve>      plotcurve;
    boost::scoped_ptr<QwtPlotMagnifier>  magnifier;
    boost::scoped_ptr<QwtPlotPanner>     panner;

    boost::scoped_ptr<QwtPlotCurve> curve1, curve2;

    JMeshPtr sourceMesh, targetMesh;
    vector<double> xData, aData, bData;

    void init();
    void setValues( const vector<double> &a, const vector<double> &b);
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
