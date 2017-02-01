#pragma once

#include <QDialog>
#include <QColorDialog>

#include "Ui_MeshFacesPartitionsDialog.hpp"
#include "MeshViewer.hpp"

#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>

class JMeshFacesPartitionsDialog : public QDialog, public Ui::MeshFacesPartitionsDialog {
    Q_OBJECT

public:
    JMeshFacesPartitionsDialog( QWidget *parent = 0);
    ~JMeshFacesPartitionsDialog();

    void setViewManager( JaalViewer *v) {
        viewManager = v;
        init();
    }
    void setMesh( const JMeshPtr &m);

private slots:
    void viewAll();
    void displayPatch();
    void checkDisk();
    void randomColors();

    void closeDialog();

private:
    JaalViewer *viewManager;
    JMeshViewerPtr meshViewer;

private:
    JMeshPtr mesh;

    boost::scoped_ptr<JMeshPartitioner> meshpart;

    boost::scoped_ptr<QwtPlotCurve> plotcurve;
    boost::scoped_ptr<QwtPlotGrid>  grid;
    boost::scoped_ptr<QwtPlotMagnifier> magnifier;
    boost::scoped_ptr<QwtPlotPanner>    panner;

    void init();
    void  initMesh();
    void makeConnections();
};
////////////////////////////////////////////////////////////////////////////////
