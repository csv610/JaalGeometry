#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshSpectrumDialog.hpp"
#include "MeshViewer.hpp"
#include "MeshSpectrum.hpp"

#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>

#include <igl/jet.h>


////////////////////////////////////////////////////////////////////////////////

class JMeshSpectrumDialog : public QDialog, public Ui::MeshSpectrumDialog {
    Q_OBJECT

public:
    JMeshSpectrumDialog( QWidget *parent = 0);
    ~JMeshSpectrumDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void genVectors();
    void viewVector();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;
    boost::scoped_ptr<JMeshSpectrum>  meshSpectrum;

    boost::scoped_ptr<QwtPlotGrid>  grid;
    boost::scoped_ptr<QwtPlotCurve> plotcurve;
    boost::scoped_ptr<QwtPlotMagnifier>  magnifier;
    boost::scoped_ptr<QwtPlotPanner>     panner;

    void init();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
